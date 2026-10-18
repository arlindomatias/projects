library(dplyr)
library(docstring)  # install.packages("docstring") - lets you read the descriptions with ?function_name


calculate_di <- function() {
  #' Discrimination index (DI)
  #'
  #' This function calculates the discrimination index for the novel object
  #' recognition test. It asks for the interaction time with the familiar and
  #' the modified object (comma or dot as decimal separator) and prints
  #' DI = (modified - familiar) / total. If the total exploration time is below
  #' 10 s, the denominator is fixed at 10. Then it returns to the main menu.

  fam <- readline("Enter the interaction time with the familiar object:\n")
  fam <- gsub(",", ".", fam) %>%
    as.numeric()
  mod <- readline("Enter the interaction time with the modified object:\n")
  mod <- gsub(",", ".", mod) %>%
    as.numeric()
  total <- fam + mod
  if (total < 10) {
    DI <- (mod - fam) / 10
  } else {
    DI <- (mod - fam) / total
  }
  print(paste("DI = ", DI))
  start()
}

calculate_y <- function() {
  #' Y-maze percentages (manual)
  #'
  #' This function calculates the Y-maze percentages from manually entered
  #' counts. It asks for the number of entries, alternations and errors, and
  #' prints the percentage of alternations (alternations / (entries - 2)) and
  #' the percentage of errors (errors / entries). Then it returns to the main
  #' menu.

  entries <- as.numeric(readline("Enter the number of entries: \n"))
  alternations <- as.numeric(readline("Enter the number of alternations: \n"))
  errors <- as.numeric(readline("Enter the number of errors: \n"))
  pct_alt <- alternations / (entries - 2)
  pct_err <- errors / entries
  print(paste("Percentage of alternations =", pct_alt))
  print(paste("Percentage of errors = ", pct_err))
  start()
}

calculate_alternations <- function() {
  #' Y-maze alternations (automatic)
  #'
  #' This function calculates the Y-maze alternations automatically from the
  #' arm sequence. It asks for the sequence of arm letters (the "|" character
  #' is ignored), counts the entries, slides a window of three consecutive
  #' entries and counts the triplets with three different arms (observed
  #' alternations). It also counts consecutive repetitions and prints the total
  #' entries, the alternation errors, the observed alternations and the
  #' alternation rate (%). Then it returns to the main menu.

  input <- readline("Enter the sequence \n")
  sequence <- unlist(strsplit(toupper(gsub("\\|", "", input)), split = ""))
  entries <- length(sequence)
  possible_alternations <- (entries - 2)
  observed_alternations <- list()
  repetitions <- 0
  if (possible_alternations > 0) {
    for (i in 1:(entries - 2)) {
      triplet <- sequence[i:(i + 2)]
      if (length(unique(triplet)) == 3) {
        observed_alternations <- append(observed_alternations, list(triplet))
      }
    }
  }
  alternation_rate <- if (possible_alternations > 0) {
    (length(observed_alternations) / possible_alternations) * 100
  } else {
    0
  }
  for (i in 1:(entries - 1)) {
    if (sequence[i] == sequence[i + 1]) {
      repetitions <- repetitions + 1
    }
  }
  print(paste("Total entries:", entries))
  print(paste("Alternation errors:", repetitions / possible_alternations))
  print(paste("Observed alternations (n):", length(observed_alternations)))
  print(paste("Alternation rate (%)", round(alternation_rate, 2)))
  start()
}

text_to_hex <- function(text) {
  #' Text to hexadecimal
  #'
  #' This function converts text to hexadecimal code. It asks for a text,
  #' turns each character into its byte value in hexadecimal and prints the
  #' result as a single string. Then it returns to the main menu.

  text <- readline("Enter the text:")
  hex <- paste0(charToRaw(text), collapse = "")
  print(hex)
  start()
}

hex_to_text <- function(hex) {
  #' Hexadecimal to text
  #'
  #' This function converts hexadecimal code back to text. It asks for a
  #' hexadecimal string, reads it two characters at a time as bytes, rebuilds
  #' the original text and prints it. Then it returns to the main menu.

  hex <- readline("Enter the code:")
  bytes <- sapply(seq(1, nchar(hex), by = 2), function(i) {
    as.raw(strtoi(substr(hex, i, i + 1), base = 16))
  })
  text <- rawToChar(bytes)
  print(text)
  start()
}

start <- function() {
  #' Main menu
  #'
  #' This function shows the main menu and runs the option typed by the user.
  #' Each option calls the matching function, which returns to this menu when
  #' it finishes. The exit option (or an invalid option) ends the program.

  cat("\n What would you like to calculate? \n",
      "\n",
      "1: Discrimination index;\n",
      "2: Manual alternation percentage;\n",
      "3: Automatic alternation percentage;\n",
      "4: Convert text to hexadecimal code;\n",
      "5: Convert hexadecimal code to text;\n",
      "6: Exit. \n",
      "\n"
  )
  type <- readline("Enter the option: \n")

  if (type == "1") {
    calculate_di()
  } else if (type == "2") {
    calculate_y()
  } else if (type == "3") {
    calculate_alternations()
  } else if (type == "4") {
    text_to_hex()
  } else if (type == "5") {
    hex_to_text()
  } else if (type == "6") {
    return(NULL)
  } else {
    cat("Option not found")
    return(NULL)
  }
}

start()