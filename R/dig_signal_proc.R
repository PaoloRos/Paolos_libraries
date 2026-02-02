# --- Library related to digital signal processing ---
#
#
# Author: Paolo Rossi
#
# Last upload: 24/01/2025
# 
# ----------------------------------------------

# --- Signal operations ---

# Scalar product
`%ps%` <- function(a, b) mean(a*b)

# Norm of a vector
norm_ps <- function(a) sqrt(a %ps% a)

# Scalar normalized product between n decimated vectors: i.e. cosine similarity
psDnN <- function(a, b, n = 1) {
  stopifnot(n >= 1, n %% 1 == 0)
  
  idx <-  seq(1, min(length(a), length(b)), by = n)
  
  a_d <- a[idx]
  b_d <- b[idx]
  
  mean(a_d * b_d) / (norm_ps(a_d)*norm_ps(b_d))
}

extract_plateaus <- function(df, fs, win = 50, min_time = 1, thr_coef = 5) {
  
  # This functions extracts plateaus discarding transitions.
  # df = data.frame, needs to contain columns of time `t` and input signal `input`!
  # fs = sampling frequency
  # win = window of the rolling mean
  # min_time = minimum time of plateaus' duration
  # thr_coef = tolerance of the MAD threashold (lower thr_coef --> narrower tolerance)
  
  x  <- df$input
  tt <- df$t
  
  # Precompute rolling means (ChatGPT's algorithm)
  # compute the slope of the window:  slope = cov(t, x)/var(t)
  # This solution is faster than computing slope by calling lm()
  mx  <- rollmean(x,  win, fill = NA, align = "center")
  mt  <- rollmean(tt, win, fill = NA, align = "center")
  mxt <- rollmean(x * tt, win, fill = NA, align = "center")
  mt2 <- rollmean(tt^2, win, fill = NA, align = "center")
  
  # Closed-form slope
  slope <- (mxt - mx * mt) / (mt2 - mt^2)
  
  out <- df
  out$slope <- slope
  
  # MAD threshold
  mad_slope <- mad(out$slope, na.rm = TRUE)                 # mad() is implemented in C, so it's faster
  thr <- thr_coef * mad_slope
  
  out$is_plateau_raw <- abs(out$slope) < thr
  out$is_plateau_raw[is.na(out$is_plateau_raw)] <- FALSE    # fix NA values
  
  # Block detection
  change <- c(TRUE, diff(out$is_plateau_raw) != 0)          # Flag of transitions: slope[i+1] - slope[i] != 0?
  block  <- cumsum(change)                                  
  
  # Plateaus' duration
  min_len <- ceiling(min_time * fs)
  
  out <- out %>%
    mutate(block = block) %>%
    group_by(block) %>%
    mutate(
      n_block = n(),
      is_plateau = is_plateau_raw & n_block >= min_len      # final flag
    ) %>%
    ungroup() %>% 
    mutate(
      plateau_id = ifelse(is_plateau, cumsum(change & is_plateau), NA)
    )
  
  return(
    out %>% dplyr::filter(is_plateau) %>%
      select(-slope, -is_plateau_raw, -block, -n_block, -is_plateau)
  )
}

# --- Fourier Transform ---

# Harmonic components
cos_k <- function(t, f, k=1) cos(2*pi*f*k*t)
sin_k <- function(t, f, k=1) sin(2*pi*f*k*t)

compute_fft <- function(df, Ta, name_sig, monolateral = TRUE) {
  
  # Compute the FFT of a signal
  # df: data frame
  # Ta: acquisition time
  # name_sig: name of the signal contained in `df`
  
  N <- nrow(df)
  
  df <- df %>% mutate(
    n = row_number(),
    f = n/Ta,
    fft = fft({{ name_sig }}),
    mod = ifelse(monolateral, (2-(f==0)), 1) * Mod(fft) / length(n),
    phase = Arg(fft)/pi*180
  )
  
  if(monolateral)
    return(df %>% slice_head(n = floor(N/2)))
  
  df
}

# Upgraded version of `compute_fft`
spectrum_fft <- function(df, fs = NULL, name_sig, monolateral = TRUE, DC = TRUE) {
  
  # Compute the spectrum of a signal using FFT
  # df: signal data frame, it's required the time `t` vector!
  # fs: sampling frequency, it's required a uniform sampling!
  # name_sig: name of the signal contained in `df`. Must be a string!
  # monolateral: defines whether monolateral or bilateral spectrum
  # DC: defines whether consider or not the DC component
  
  # Note aboute the scale factor of the module:
  # To preserve the energy content, when the spectrum is monolateral
  # it's needed to double the module of the trasform, excluding
  # DC component and component at Nyquist frequency
  
  if (!is.character(name_sig) || length(name_sig) != 1 || nchar(name_sig) == 0)
    stop("`name_sig` must be a non-empty string of length 1.")
  
  if (!name_sig %in% names(df))
    stop(sprintf("Column '%s' not found in data frame.", name_sig))
  
  
  N <- nrow(df)
  
  if (is.null(fs)) {
    dt <- diff(df$t)
    dt_mean <- mean(dt)
    
    if (sd(dt) / dt_mean > 1e-6) {
      stop("Non-uniform sampling detected.")
    }
    
    fs <- 1 / dt_mean
  }
  
  out <- df %>% mutate(
    k = 0:(N-1),
    f = k*fs/N,
    fft = fft(df[[name_sig]]),
    mod = Mod(fft) / N * if (monolateral) {           # scale factor
      ifelse(k == 0 | (N %% 2 == 0 & k == N/2), 1, 2)
    } else {
      1
    },
    phase = Arg(fft)/pi*180
  )
  if(!DC)
    out <- out[-1,]
  
  if(monolateral)
    return(out %>% slice_head(n = floor(N/2) + 1))
  
  return(out)
}

dyn_compensation <- function(df, tf, name_sig = "y", padding = TRUE, alpha = 2) {
  
  # This function provides the dynamic compensation of an
  # LTI proper system
  
  # df: contains the signal(t, `name_sig`)
  # alpha: scale factor of padding
  # tf: transfer function, it must be a class `tf`
  
  if (!is.character(name_sig) || length(name_sig) != 1 || nchar(name_sig) == 0)
    stop("`name_sig` must be a non-empty string of length 1.")
  
  if (!name_sig %in% names(df))
    stop(sprintf("Column '%s' not found in data frame.", name_sig))
  
  # Sampling parameters
  dt <- diff(df$t)
  Ts <- mean(dt)            # sampling period
  
  if (sd(dt) / Ts > 1e-6) 
    stop("Non-uniform sampling detected.")
  
  fs <- 1 / Ts              # sampling frequency
  
  y <- df[[name_sig]]
  
  if(padding) {
    N <- nrow(df)
    M <- alpha*N            # padding length
    
    pad_left  <- floor((M - N)/2)
    pad_right <- ceiling((M - N)/2)
    t_start <- df$t[1] - pad_left * Ts
    
    out <- tibble(
      t = seq(from = t_start, by = Ts, length.out = M),
      y = c(
        rep(0, pad_left),
        y,
        rep(0, pad_right)
      )
    )
    
    stopifnot(isTRUE(all.equal(           # checks whether the signal is centered after padding
      out$y[(pad_left+1):(pad_left+N)],
      y,
      tolerance = 1e-10
      )))
  }
  
  freq <- tibble(
    k = 0:(M-1),
    f = ifelse(k <= M/2, k, k-M) * fs/M,   # FFT is symmetric respect to the frequency!
    w = 2*pi*f
  )
  
  Hiw <- freqresp(H, freq$w)
  
  out <- out %>% mutate(
    Y = fft(y),
    U = Y/Hiw,
    u = Re(fft(U, inverse = TRUE)) / M      # Re() to avoid numerical residuals in Im() component
                                            # fft() is not normalized, so it's divided by `M`
  )
  
  return(out[ (pad_left+1):(pad_left+N), ])
}

# --- Filters ---

gaussian_kernel <- function(
    x, sigma, radius = ceiling(3 * sigma), pad = c("none", "replicate", "reflect")
    ) {
  
  # Gaussian mask: law-pass filter, FIR and linear
  # sigma: number of samples to consider in the gauss. kernel
  
  stopifnot(is.numeric(x), length(x) > 0, is.numeric(sigma), sigma > 0)
  pad <- match.arg(pad)
  
  # Kernels calculation
  k <- seq(-radius, radius)
  h <- exp(-(k^2) / (2 * sigma^2))
  h <- h / sum(h)  
  
  if (pad == "none") {
    y <- stats::filter(x, h, sides = 2) # convolution
    return(as.numeric(y))
  }
  
  # Padding of signal tails
  if (pad == "replicate") {
    x_pad <- c(rep(x[1], radius), x, rep(x[length(x)], radius))
  } else { # reflect
    
    left  <- rev(x[2:(radius + 1)])
    right <- rev(x[(length(x) - radius):(length(x) - 1)])
    x_pad <- c(left, x, right)
  }
  
  # Convolution and truncation of the excess tails
  y_pad <- stats::filter(x_pad, h, sides = 2)
  y <- y_pad[(radius + 1):(radius + length(x))]
  
  as.numeric(y)
}
