# --- Library related to graphics formatting ---
#
#
# Author: Paolo Rossi
# 
# ----------------------------------------------

#======== Scientific format ========

require(knitr)
require(kableExtra)

sci_to_string <- function(x, digits = 1) {
  s <- formatC(x, format = "e", digits = digits)
  
  # Rimuove 'e+0' o 'e+00' → 'e'
  s <- gsub("e\\+0*", "e", s)
  
  # Rimuove 'e-0' o 'e-00' → 'e-'
  s <- gsub("e-0*", "e-", s)
  
  # Converte in LaTeX
  s <- gsub("e(-?\\d+)", " \\\\times 10^{\\1}", s)
  
  # Aggiunge i delimitatori LaTeX
  s <- paste0("$", s, "$")
  return(s)
}

to_latex_sci <- function(x, digits = 2) {
  exp <- floor(log10(abs(x)))
  mantissa <- round(x / 10^exp, digits)
  paste0("$", mantissa, " \\times 10^{", exp, "}$")
}

# Stampa "10^ordine" per gli assi dei grafici: estrae solo l'esponente
notazione_scientifica <- function(x) {
  return(parse(text = paste0("10^", log10(x))))
}