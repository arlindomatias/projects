library(dplyr)


calculo_id <- function() {
  fam <- readline("Digite o tempo de interação com o objeto familiar:\n")
  fam <-  gsub(",", ".", fam) %>%
    as.numeric() 
  mod <- readline("Digite o tempo de interação com o objeto modificado:\n")
  mod <- gsub(",", ".", mod) %>%
    as.numeric() 
  soma <- fam + mod
  if(soma < 10) {
    ID = (mod - fam) / 10} else
    {ID = (mod - fam) / soma}
  print(paste("ID = ", ID))
inicio()
}

calculo_y <- function() {
  ent <- as.numeric(readline("Digite o número de entradas: \n"))
  alt <- as.numeric(readline("Digite o número de alternâncias: \n"))
  err <- as.numeric(readline("Digite o número de erros: \n"))
  pa = alt / (ent - 2)
  pe = err / ent
  print(paste("Porcentagem de alternancias =", pa))
  print(paste("Porcentagem de erros = ", pe))
inicio()
}

calcular_alternancias <- function() {
  input <- readline("Digite a sequência \n")
  sequencia <- unlist(strsplit(toupper(gsub("\\|", "", input)), split = ""))
  entradas <- length(sequencia)
  alternancias_possiveis <- (entradas - 2)
  alternancias_observadas <- list()
  repeticoes <- 0
  
  
  if (alternancias_possiveis > 0) {
    for (i in 1:(entradas - 2)) {
      tripla <- sequencia[i:(i + 2)]
      if (length(unique(tripla)) == 3) {
        alternancias_observadas <- append(alternancias_observadas, list(tripla))
      }
    }
  }
  
  taxa_alternancia <- if (alternancias_possiveis > 0) {
    (length(alternancias_observadas) / alternancias_possiveis) * 100
  } else {
    0
  }
  
  for (i in 1:(entradas - 1)) {
    if (sequencia[i] == sequencia[i + 1]) {
      repeticoes <- repeticoes + 1
    }
  }  
  
  print(paste(
    "Total de entradas:", entradas))
  print(paste(
    "Erros de alternância:", repeticoes/alternancias_possiveis))
  print(paste(
    "Alternâncias observadas (n):", length(alternancias_observadas)))
  print(paste(  
    "Taxa de alternância (%)", round(taxa_alternancia, 2)
  ))
inicio()
}

calculo_baco <- function() {
  peso <- readline("Digite o peso do animal em gramas: ") %>%
  as.numeric()
  baco <- readline("Digite o peso do baço em miligramas: ") %>%
  as.numeric()
  per_baco <- (baco / (peso * 1000)) * 100
  print(paste("Porcentagem do peso do baço é", per_baco,"%"))
inicio()
}

text_to_hex <- function(texto) {
  texto <- readline("Insira o texto:")
  hex <- paste0(charToRaw(texto), collapse = "")
  print(hex)
inicio()
}

hex_to_text <- function(hex) {
  hex <- readline("Insira o código:")
  bytes <- sapply(seq(1, nchar(hex), by = 2), function(i) {
    as.raw(strtoi(substr(hex, i, i+1), base = 16))
  })
  texto <- rawToChar(bytes)
  print(texto)
inicio()
}


inicio <- function() {
  cat("\n O que deseja calcular? \n",
      "\n",
    "1: Índice de discriminação;\n",
    "2: Porcentagem de alternâncias manual;\n",
    "3: Porcentagem de alternâncias automática;\n",
    "4: Porcentagem de peso do baço;\n",
    "5: Converter texto para código hexadecimal;\n",
    "6: Converter código hexadecimal para texto;\n",
    "7: Sair. \n",
    "\n"
  )
  
  type <- readline("Digite a opção: \n")

  if(type == "1") {
    calculo_id()
  } else if(type == "2") {
    calculo_y() 
  } else if(type == "3") {
    calcular_alternancias()
  } else if(type == "4") {
    calculo_baco() 
  } else if(type == "5") {
    text_to_hex()
  } else if(type == "6") {
    hex_to_text()
  } else if(type == "7") {
    return(NULL)
    break
  } else {
    cat("Opção não encontrada")
    return(NULL)
  }
}

inicio()