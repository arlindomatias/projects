library(plumber)
library(jsonlite)
library(ggplot2)
library(dplyr)

setwd("C:/Users/User/Documents/github/projects/website")
#* @apiTitle Serviço de Gráficos Estatísticos
#* @apiDescription Recebe grupos numéricos e retorna um PNG com pontos individuais, média e erro padrão.


# =========================================
# Função auxiliar para calcular a média e SEM
# =========================================
calc_summary <- function(values) {
  tibble(
    mean = mean(values),
    sem  = sd(values) / sqrt(length(values))
  )
}

# =========================================
# Endpoint principal: retorna um PNG
# =========================================

#* Gera gráfico para grupos
#* @post /plots
#* @serializer contentType list(type="image/png")
function(req, res) {
  
  body <- req$postBody
  
  if (is.null(body) || nchar(body) == 0) {
    res$status <- 400
    return(list(error = "Body vazio."))
  }
  
  parsed <- tryCatch(jsonlite::fromJSON(body), error = function(e) NULL)
  
  if (is.null(parsed$groups)) {
    res$status <- 400
    return(list(error = "Formato esperado: { groups: [[...],[...],...] }"))
  }
  
  groups <- parsed$groups
  
  # transformar em data frame longo
  df <- do.call(rbind, lapply(seq_along(groups), function(i) {
    tibble(
      group = paste0("Grupo ", i),
      value = as.numeric(groups[[i]])
    )
  }))
  
  # resumo para média + SEM
  summary_df <- df %>%
    group_by(group) %>%
    summarise(
      mean = mean(value),
      sem  = sd(value) / sqrt(n())
    )
  
  # gerar a imagem
  img <- tempfile(fileext = ".png")
  
  # gráfico
  p <- ggplot(df, aes(x = group, y = value)) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
    geom_point(data = summary_df, aes(y = mean), color = "red", size = 3) +
    geom_errorbar(
      data = summary_df,
      aes(y = mean, ymin = mean - sem, ymax = mean + sem),
      width = 0.2,
      color = "red"
    ) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Pontos individuais, média e erro padrão",
      x = "Grupo",
      y = "Valor"
    )
  
  # salvar PNG
  ggsave(img, p, width = 6, height = 4, dpi = 150)
  
  # retornar o arquivo como binário (PNG)
  readBin(img, "raw", n = file.info(img)$size)
}

#pr <- plumb("api.R")
#pr$run(host="0.0.0.0", port=8000)
