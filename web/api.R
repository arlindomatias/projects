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

# -------------------------
# CORS
# -------------------------
#* @filter cors
function(req, res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  res$setHeader("Access-Control-Allow-Methods", "*")
  res$setHeader("Access-Control-Allow-Headers", "*")
  if (req$REQUEST_METHOD == "OPTIONS") {
    return(plumber::forward())
  }
  plumber::forward()
}:forward()


#* @post /plots
function(req, res){
  body <- jsonlite::fromJSON(req$postBody)
  groups <- body$groups
  
  df <- data.frame(
    value = unlist(groups),
    group = rep(seq_along(groups), times = sapply(groups, length))
  )
  
  p <- ggplot(df, aes(x = factor(group), y = value)) +
    geom_jitter(width = 0.1) +
    stat_summary(fun = mean, geom = "point", color = "red", size = 3) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)
  
  tmp <- tempfile(fileext = ".png")
  ggsave(tmp, p, width = 6, height = 4, dpi = 150)
  
  res$setHeader("Content-Type", "image/png")
  readBin(tmp, "raw", n = file.info(tmp)$size)
}

plumber::pr("api.R")$run(host="0.0.0.0", port=8000)

pr <- plumb("api.R")
pr$run(host="0.0.0.0", port=8000)
