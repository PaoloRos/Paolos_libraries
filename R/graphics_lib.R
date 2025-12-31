# --- Library related to graphics formatting ---
#
#
# Author: Paolo Rossi
#
# Last upload: 12/31/2025
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

notazione_scientifica <- function(x) {
  # It prints "10^order" for the plot axes: it extracts the exponent only
  return(parse(text = paste0("10^", log10(x))))
}

#======== Bode Plot ========

bode_plot <- function(tf, w, iu = 1)  
{
  # Parameters class:
  # transfer funciont: class(tf) = "tf", "ss", "zpk"
  # frequency vector (logharitmic): class(w) = "numeric"
  # inputs number: class(iu) = "integer"
  
  # Frequency response
  b <- bode(tf, w, iu)
  
  df <- tibble(
    f = b$w/(2*pi),
    mag = b$mag,
    phase = b$phase
  )
  
  # Rescaling the secondary axis
  phase_min <- min(df$phase)
  phase_max <- max(df$phase)
  mag_min   <- min(df$mag)
  mag_max   <- max(df$mag)
  
  a <- (mag_max - mag_min)/(phase_max - phase_min)
  c <- mag_min - a*phase_min
  
  df$phase_scaled <- a*df$phase + c
    
  # Output plot
  df %>% ggplot(aes(x = f)) +
    geom_line(aes(y = mag, color = "Magnitude")) +
    geom_line(aes(y = phase_scaled, color = "Phase"), linetype = "dashed") +
    scale_x_log10(labels = notazione_scientifica) +
    annotation_logticks(sides = "b") +
    scale_y_continuous(
      name = "Magnitude (dB)",
      sec.axis = sec_axis(~ (. - c)/a, name = "Phase (deg)")  # trasformazione inversa
    ) +
    labs(
      x = "Frequency (Hz)", color = "Quantity"
    )
}