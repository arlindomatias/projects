library(tidyverse)
library(httr2)
library(jsonlite)


# ── 1. Leitura e parsing do arquivo WoS ──────────────────────
parse_wos <- function(filepath) {
  lines <- readLines(filepath, encoding = "UTF-8", warn = FALSE)
  
  records <- list()
  current <- list(AB = character())
  current_field <- NULL
  
  for (line in lines) {
    if (nchar(line) < 2) next
    
    tag <- substr(line, 1, 2)
    content <- if (nchar(line) > 3) trimws(substr(line, 4, nchar(line))) else ""
    
    if (tag == "ER") {
      # Fim do registro
      if (length(current$AB) > 0) {
        current$AB <- paste(current$AB, collapse = " ")
        records <- append(records, list(current))
      }
      current <- list(AB = character())
      current_field <- NULL
      next
    }
    
    if (tag == "AB") {
      current_field <- "AB"
      current$AB <- c(current$AB, content)
    } else if (tag == "TI") {
      current_field <- "TI"
      current$TI <- content
    } else if (tag == "PY") {
      current$PY <- content
      current_field <- NULL
    } else if (tag == "  " && !is.null(current_field) && current_field == "AB") {
      # Continuação do campo AB (linhas iniciadas com espaço)
      current$AB <- c(current$AB, content)
    } else {
      current_field <- NULL
    }
  }
  
  bind_rows(records)
}

raw <- "~/github/private/wos.txt"
corpus <- parse_wos(raw)
cat("Registros lidos:", nrow(corpus), "\n")


# ── Função de anotação via BERN2 ──────────────────────────────
annotate_bern2 <- function(text, record_id) {
  
  resp <- tryCatch({
    request("http://bern2.korea.ac.kr/plain") %>%
      req_body_json(list(text = text)) %>%
      req_timeout(60) %>%
      req_retry(
        max_tries    = 3,
        is_transient = ~ resp_status(.x) %in% c(429, 500, 503)
      ) %>%
      req_perform()
  }, error = function(e) {
    message("Erro no registro ", record_id, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(resp)) return(tibble())
  if (resp_status(resp) != 200) {
    message("Status ", resp_status(resp), " no registro ", record_id)
    return(tibble())
  }
  
  # ── Sanitizar JSON: substitui NaN por null antes de parsear ──
  json_text <- resp %>%
    resp_body_string() %>%
    str_replace_all("\\bNaN\\b", "null")  # NaN → null (válido em JSON)
  
  result <- tryCatch(
    fromJSON(json_text, simplifyVector = FALSE),
    error = function(e) {
      message("Erro ao parsear JSON no registro ", record_id, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(result)) return(tibble())
  
  annotations <- result$annotations
  if (is.null(annotations) || length(annotations) == 0) return(tibble())
  
  map_dfr(annotations, function(ann) {
    tibble(
      record_id = as.character(record_id),
      mention   = ann$mention %||% NA_character_,
      type      = ann$obj     %||% NA_character_,
      ncbi_id   = if (!is.null(ann$id) && length(ann$id) > 0)
        paste(ann$id, collapse = "; ")
      else NA_character_,
      start     = ann$span$begin %||% NA_integer_,
      end       = ann$span$end   %||% NA_integer_
    )
  })
}

# Operador %||% (null coalescing) — caso não esteja disponível
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Preparar corpus ───────────────────────────────────────────
corpus_clean <- corpus %>%
  filter(!is.na(AB), AB != "") %>%
  mutate(
    record_id = as.character(row_number()),
    TI        = if_else(is.na(TI), "", TI),
    # Concatena título + abstract para melhor contexto
    full_text = paste0(TI, ". ", AB)
  )

# ── Teste com um registro ─────────────────────────────────────
teste <- annotate_bern2(
  text      = corpus_clean$full_text[1],
  record_id = "1"
)
print(teste)

# ── Loop completo com checkpoint ──────────────────────────────
# Checkpoint: salva a cada 50 registros para não perder progresso
checkpoint_file <- "bern2_checkpoint.rds"

# Retoma de onde parou se o arquivo existir
if (file.exists(checkpoint_file)) {
  checkpoint    <- readRDS(checkpoint_file)
  all_annot     <- checkpoint$all_annot
  start_from    <- checkpoint$last_record + 1
  message("Retomando do registro ", start_from)
} else {
  all_annot  <- list()
  start_from <- 1
}

for (j in start_from:nrow(corpus_clean)) {
  
  if (j %% 10 == 0) message("Processando registro ", j, "/", nrow(corpus_clean))
  
  Sys.sleep(0.4)  # max ~2.5 req/s — conservador para não ser bloqueado
  
  result <- annotate_bern2(
    text      = corpus_clean$full_text[j],
    record_id = corpus_clean$record_id[j]
  )
  
  all_annot[[j]] <- result
  
  # Salva checkpoint a cada 50 registros
  if (j %% 50 == 0) {
    saveRDS(list(all_annot = all_annot, last_record = j), checkpoint_file)
    message("  Checkpoint salvo no registro ", j)
  }
}

# Salva checkpoint final
saveRDS(list(all_annot = all_annot, last_record = nrow(corpus_clean)), checkpoint_file)

# ── Consolidar resultados ─────────────────────────────────────
annotations_raw <- bind_rows(all_annot)

# ── Filtrar genes/proteínas ───────────────────────────────────
# Tipos retornados pelo BERN2: gene, disease, drug,
# species, cell_line, cell_type, DNA, RNA
gene_annotations <- annotations_raw %>%
  filter(str_detect(type, regex("gene|protein", ignore_case = TRUE)))

cat("Total de menções de genes/proteínas:", nrow(gene_annotations), "\n")

# ── Contagem por mention ──────────────────────────────────────
counts_by_mention <- gene_annotations %>%
  count(mention, sort = TRUE) %>%
  rename(n_mentions = n)

# ── Contagem por NCBI ID (agrupa sinônimos) ───────────────────
counts_by_id <- gene_annotations %>%
  group_by(ncbi_id) %>%
  summarise(
    n_mentions     = n(),
    n_articles     = n_distinct(record_id),
    mentions_found = paste(sort(unique(mention)), collapse = "; ")
  ) %>%
  arrange(desc(n_mentions))

print(counts_by_mention, n = 50)
print(counts_by_id, n = 50)

# ── Exportar ──────────────────────────────────────────────────
write_csv(counts_by_mention, "bern2_counts_mention.csv")
write_csv(counts_by_id,      "bern2_counts_by_id.csv")

# Tabela completa com título e ano do artigo
full_table <- gene_annotations %>%
  left_join(
    corpus_clean %>% select(record_id, TI, PY),
    by = "record_id"
  )

write_csv(full_table, "bern2_full_annotations.csv")
