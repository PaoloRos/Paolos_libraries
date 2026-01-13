# --- Library related to statistics [under development] ---
#
#
# Author: Paolo Rossi
#
# Last upload: 12/31/2025
# 
# ----------------------------------------------
# 

#======== Inclusione librerie ======== 

require(tidyverse)
require(roxygen2)

# 0 ======================== 0

example_url <- function(example) {
  url = paste0("https://paolobosetti.quarto.pub/data/", example)
  return(url)
}

#======== Criterio di Chauvenet ======== 

require(glue)

chauvenet <- function(x, threshold=0.5){
  abs.diff <- abs(x-mean(x))/sd(x) #vettore delle differenze
  s0 <- max(abs.diff)
  i0 <- which.max(abs.diff) #posizione della max. diff.
  
  #calcolo la CDF+
  freq <- length(x) * pnorm(s0, lower.tail = F) #pnorm = pdf funz. normale
  
  result <- list(
    s0 = s0,
    index = i0,
    value = x[i0],
    reject = freq < threshold #TRUE = RIUFIUTO IL PUNTO perché freq < threshold
  )
  samp <- deparse(substitute(x)) #si prende il nome dell'argomento passato
  print(glue("Chauvenet's criterion for sample {samp}"))
  print(glue("Suspect outlier: {i0}, value {x[i0]}")) #passo stringa -> metto espressione R tra {}
  print(glue("Expected frequency: {freq}, threshold: {threshold}"))
  print(glue("Decision: {d}", d=ifelse(result$reject, "reject it", "keep it")))
  invisible(result) #non vedo il risultato, anziché di return()
}

# 0 ======================== 0

#======== Normalità dei residui ======== 

resNormality <- function(df) { #, resid = "resid"
  p <- shapiro.test(df$resid)$p.value
  
  plot <- df %>%
    ggplot(aes(sample=resid)) +
      geom_qq_line(color="red") +
      geom_qq() +
      labs(x = "Quantili teorici", y = "Quantili campionari")
  
  plot + 
    annotate("text", x = 0, y = 0, 
             label = paste("Shapiro-Wilk p-value:", round(p, 4)),
             size = 5, hjust = 0)
}

# 0 ======================== 0

#======== Plot dei residui ======== 

#Residuals plot nei parametri tabella tibble `t` e colonna `n`, nome della colonna  da stampare "nome"
residPlot <- function(t, n, nome="inserire nome", unita=NULL) {
  if (!is.null(unita)) {
    nome_y <- paste0("Residui (", unita, ")")
  } else {
    nome_y <- "Residui"
  }
  t %>%
    ggplot(aes(x={{n}}, y=resid)) +
    geom_hline(yintercept = 0, color = "black", linetype = "dotted", linewidth = 0.5) +
    geom_point() +
    labs(x = nome, y = nome_y)
}

# {{}} -> permette di far interpretare `n` come il nome di una colonna, piuttosto che di una variabile

# 0 ======================== 0

#======== Plot dei quantili ======== 

#sistemare

quantilePlot <- function(t) {
  t %>%
    ggplot(aes(sample=resid)) +
    geom_qq() +
    geom_qq_line(color="red") +
    labs(x = "Quantili teorici", y = "Quantili campionari")
}

# 0 ======================== 0

#======== Uncertainty calculation ========

# Usa roxygen2 per la documentazine

# There are three different kind of uncertainty:
#   1. Standard uncertainty, defined as the standard deviation of the sample mean:
#     u_x = s_x/sqrt(n), where s_x = std. deviation and n = number of the sample elements

#   2. Relative uncertainty, defined as the u_x divided by the absolute value of the sample mean:
#     u_x,rel = u_x/n (adimensional)

#   3. Extended uncertainty, defined as the u_x multiplied by a coverage factor:
#     U_x = k_n * u_x, where k_n depends on n.
#     If n < 50 elements, then k_n is caculated as the quantile α/2 of the T.Student's distribution with, in general, n DoF.
#     Whereas, if n >= 50 the T.Student is undistinguishable from the Normal distribution, then is used that last one distribution on the quantile α/2.
#     The probability α = 1 - conf, where conf is the confidential level (between 0 and 1) of the uncertainty.
uncertainty <- function(vec, type = "standard", conf = 0.95){
  #eliminates every NA
  vec <- na.omit(vec)
  
  n <- length(vec)
  
  #Std uncertainty
  ux <- sd(vec)/sqrt(n)
  
  if(type == "standard")
    return(ux)
  
  if(type == "relative")
    return(ux/mean(vec))
  
  if(type == "extended"){
    #coverage factor
    alpha = 1 - conf
    if(n<50) {
      kn <- qt(alpha/2, df = n, lower.tail = F)
    }
    else{
      kn <- qnorm(alpha/2, lower.tail = F)
    }
    return(kn*ux)
  }
  
  return(NULL) #if the user digits a wrong  type
}

# 0 ======================== 0

#======== T test.2 ========

# This test checks automatically whether the two samples are omoscedastic

t.test2 <- function(
    x, y = NULL,
    alternative = c("two.sided", "less", "greater"),
    mu = 0,
    paired = FALSE,
    conf.level = 0.95
) {
  x <- na.omit(x)
  y <- na.omit(y)
  
  if(!is.null(y) && !paired)
    vt <- var.test(x,y)
  
  #cat("condizionale:", vt$p.value > 1 - conf.level)
  
  t.test(
    x, y,
    alternative,
    mu,
    paired,
    var.equal = (vt$p.value > 1 - conf.level),
    conf.level
  )
}

# 0 ======================== 0

#======== Uncertainty calculation ========

#Formatta in notazione scientifica i valori negli assi dei grafici

notazione_scientifica <- function(x) {
  return(parse(text = paste0("10^", log10(x))))
}

# 0 ======================== 0