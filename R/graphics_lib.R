# --- Library related to graphics formatting ---
#
#
# Author: Paolo Rossi
#
# Last upload: 24/01/2025
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

ggbodeplot_continous <- function(tf, fs, fmin=1, fmax=NULL, df=0.01) {
  
  # Bode's plot in continous-time domain
  # `tf`: transfer function
  # `fs`: sampling freq.
  
  if (is.null(fmax)) fmax <- fs/2
  
  # vector of points for each order of magnitude (OOM):
  pts <- 10^seq(0, 1, df) %>% tail(-1)
  # vector of OOMs:
  ooms <- 10^(floor(log10(fmin)):ceiling(log10(fmax)-1))
  # combine pts and ooms:
  freqs <- as.vector(pts %o% ooms)
  
  # warning: bode wants pulsation!
  bode(tf, freqs*2*pi) %>% {
    tibble(
      f=.$w/(2*pi), 
      `Magnitude (dB)`=.$mag, 
      `Phase (deg)`=.$phase)} %>%
    pivot_longer(-f) %>% 
    ggplot(aes(x=f, y=value)) +
    geom_line() +
    scale_x_log10(
      breaks = 10^seq(floor(log10(fmin)), ceiling(log10(fmax))),
      minor_breaks=scales::minor_breaks_n(10), 
      labels= ~ latex2exp::TeX(paste0("$10^{", log10(.), "}$"))
    ) +
    facet_wrap(~name, nrow=2, scales="free") +
    labs(x="Frequency (Hz)", y="")
}

ggbodeplot_digital <- function(tf, fs, fmin=1e-3, fmax=NULL, npts=2000, xaxis = c("frequency", "omega")) {
  
  # Bode's plot in discrete-time domain
  # `tf`: transfer function
  # `fs`: sampling freq.
  # `xaxis`: defines whether plot the physical freq. or the discrete angular one
  
  xaxis <- match.arg(xaxis)
  if (is.null(fmax)) fmax <- fs/2
  
  f <- 10^seq(log10(fmin), log10(fmax), length.out = npts)  # physical frequency
  w <- 2*pi*f/fs                                            # discrete angular frequency (rad/sample) [0,π]
  
  H <- freqz(tf$b, tf$a, w)                               # digital frequency response
  
  xtitle = ifelse(xaxis == "frequency", "Frequency (Hz)", latex2exp::TeX("Discrete \\omega (rad/sample)"))
  
  p <- tibble(
    f = H$f * ifelse(xaxis == "frequency", fs/(2*pi), 1),    # H$f returns the previous `w`
    `Magnitude (dB)` = 20*log10(Mod(H$h)),
    `Phase (deg)` = Arg(H$h)*180/pi
  ) %>%
    pivot_longer(-f) %>%
    ggplot(aes(x=f, y=value)) +
    geom_line() +
    scale_x_log10(
      breaks = 10^seq(floor(log10(fmin)), ceiling(log10(fmax)), by = 1),
      minor_breaks = scales::minor_breaks_n(10),
      labels = ~ latex2exp::TeX(paste0("$10^{", log10(.), "}$"))
    ) +
    facet_wrap(~name, nrow=2, scales="free") +
    labs(x=xtitle, y="") 
  
  if (xaxis=="omega") {
    p <- p + scale_x_log10(
      breaks = 10^seq(floor(log10(2*pi*fmin/fs)), ceiling(log10(2*pi*fmax/fs)), by = 1),
      minor_breaks = scales::minor_breaks_n(10),
      labels = ~ latex2exp::TeX(paste0("$10^{", log10(.), "}$"))
    ) + 
      geom_vline(xintercept = pi, linetype = 4, color = "blue") +
      annotate("text", x = pi, y = Inf, label = expression(pi), vjust = 1.2, hjust = 1.2)
  }
  
  # output
  p
}