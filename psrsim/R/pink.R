## pink noise generation algorithm,
## code is taken from {tuneR} package (fn not exported)

#' Simulate response data with pink noise profile
#'
#' @param N Number of trials.
#' @param alpha Value of alpha (default 1).
#' @returns A vector with response values.
#' @export
pink <- function(N, alpha = 1) {
  ## Fourier frequencies
  f <- seq(from = 0, to = pi, length.out = (N/2 + 1))[-c(1, (N/2 + 1))] 
  f_ <- 1/f^alpha # power law
  RW <- sqrt(0.5 * f_) * rnorm(N/2 - 1) # for the real part
  IW <- sqrt(0.5 * f_) * rnorm(N/2 - 1) # for the imaginary part
  fR <- complex(real = c(rnorm(1), RW, rnorm(1), RW[(N/2 - 1):1]),
                imaginary = c(0, IW, 0, -IW[(N/2 - 1):1]),
                length.out = N)
  ## Those complex numbers that are to be back transformed for
  ## Fourier Frequencies 0, 2pi/N, 2*2pi/N, ..., pi, ..., 2pi-1/N
  ## Choose in a way that frequencies are complex-conjugated and
  ## symmetric around pi 0 and pi do not need an imaginary part
  reihe <- fft(fR, inverse = TRUE) # back into the time domain
  return(Re(reihe)) # don't need imaginary part
}
