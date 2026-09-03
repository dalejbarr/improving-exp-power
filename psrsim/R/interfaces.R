## R interfaces to C++ functions so as to provide error checking

#' Simulate data from a one-factor design
#'
#' @param etastep Step number to determine eta^2 value (1-7).
#' @param nPart Number of participants (25, 40, or 60).
#' @param nLevels Number of levels of the independent variable (2, 4, 8).
#' @param nReps Number of repetitions within each level for each participant (6, 12, 24).
#' @param nSubblocks Number of subblocks.
#' @param mu Intercept.
#' @param etaSQ Eta-squared (measure of the effect size).
#' @param prop_rint Proportion of variance accounted for by random intercepts.
#' @param prop_rslp Proportion of variance accounted for by random slopes.
#' @param errstr Error structure ("noauto", "expDecay", "pinkNoise", "randomWalk", or "mixed")
#' @returns A data frame with simulated data.
#' @export
simdata_1f <- function(etastep, nPart, nLevels, nReps, nSubblocks,
                       mu = 0,
                       etaSQ = etaSQ_1f(nPart, nLevels, nReps)[etastep],
                       prop_rint = .35,
                       prop_rslp = .11,
                       errstr = c("noauto",
                                  "expDecay",
                                  "pinkNoise",
                                  "randomWalk",
                                  "mixed")) {

  if (!etastep %in% seq_len(formals(etaSQ_1f)$length.out)) {
    stop("etastep must be an integer between 1 and ",
         formals(etaSQ_1f)$length.out)
  }
  
  if (!(nSubblocks %in% numbers::divisors(nReps))) {
    stop("'nSubblocks' must be an integer divisor of 'nReps'")
  }

  this_errstr <- match.arg(errstr)
  
  estr_id <- which(c("noauto",
                     "expDecay",
                     "pinkNoise",
                     "randomWalk",
                     "mixed") == this_errstr)

  simdata_1f_int(nPart, nLevels, nReps, nSubblocks,
                 mu, etaSQ, prop_rint, prop_rslp,
                 estr_id)
}

#' Simulate data from a 2x2 factorial design
#'
#' @param etastep Step number to determine eta^2 value (1-7).
#' @param nPart Number of participants (28, 40, or 50).
#' @param nReps Number of repetitions within each level for each participant (10, 16).
#' @param nSubblocks Number of subblocks.
#' @param mu Intercept.
#' @param etaSQ_row Eta-squared (measure of the effect size) for the row factor ("A").
#' @param etaSQ_col Eta-squared for the column factor ("B").
#' @param etaSQ_rcx Eta-squared for the interaction.
#' @param prop_rint Proportion of variance accounted for by random intercepts.
#' @param prop_rslp_row Proportion of variance accounted for by random slopes for row.
#' @param prop_rslp_col Proportion of variance accounted for by random slopes for column.
#' @param prop_rslp_rcx Proportion of variance accounted for by random slopes for interaction.
#' @param errstr Error structure ("noauto", "expDecay", "pinkNoise", "randomWalk", or "mixed")
#' @returns A data frame with simulated data.
#' @export
simdata_2x2ww <- function(etastep, nPart, nReps, nSubblocks,
                          mu = 0,
                          etaSQ_row = etaSQ_2x2ww(nPart, nReps)[etastep],
                          etaSQ_col = etaSQ_2x2ww(nPart, nReps)[etastep],
                          etaSQ_rcx = etaSQ_2x2ww(nPart, nReps)[etastep],
                          prop_rint = .11,
                          prop_rslp_row = .11,
                          prop_rslp_col = .11,
                          prop_rslp_rcx = .11,
                          errstr = c("noauto",
                                     "expDecay",
                                     "pinkNoise",
                                     "randomWalk",
                                     "mixed")) {

  if (!etastep %in% seq_len(formals(etaSQ_2x2ww)$length.out)) {
    stop("etastep must be an integer between 1 and ",
         formals(etaSQ_2x2ww)$length.out)
  }
  
  if (!(nSubblocks %in% numbers::divisors(nReps))) {
    stop("'nSubblocks' must be an integer divisor of 'nReps'")
  }

  this_errstr <- match.arg(errstr)
  
  estr_id <- which(c("noauto",
                     "expDecay",
                     "pinkNoise",
                     "randomWalk",
                     "mixed") == this_errstr)

  simdata_2x2ww_int(nPart, nReps, nSubblocks,
                    mu,
                    etaSQ_row, etaSQ_col, etaSQ_rcx,
                    prop_rint,
                    prop_rslp_row, prop_rslp_col, prop_rslp_rcx,
                    estr_id)
}

etaSQ_1f <- function(nPart, # c(25L, 40L, 60L),
                     nLevels, # c(2L, 4L, 8L),
                     nReps, # c(6L, 12L, 24L),
                     length.out = 7L) {
  etamax <-
    c(0.23029872, 0.21840330, 0.20396078, 0.18432309, 0.15264338, 0.14230249,
      0.14756355, 0.12519984, 0.11236103,
      0.18841444, 0.16232683, 0.15264338, 0.14230249, 0.12519984, 0.11236103,
      0.11570436, 0.09785193, 0.08972179,
      0.15264338, 0.13114877, 0.13114877, 0.11895377, 0.09785193, 0.08972179,
      0.09785193, 0.08077747, 0.07335734)

  np_opts <- c(25L, 40L, 60L)
  lv_opts <- c(2L, 4L, 8L)
  nr_opts <- c(6L, 12L, 24L)
  
  if (!nPart %in% np_opts) {
    stop("'nPart' must be one of: ",
         paste(np_opts, collapse = ", "))
  }
  if (!nLevels %in% lv_opts) {
    stop("'nLevels' must be one of: ",
         paste(lv_opts, collapse = ", "))
  }
  if (!nReps %in% nr_opts) {
    stop("'nReps' must be one of: ",
         paste(nr_opts, collapse = ", "))
  }

  np <- which(np_opts == nPart)
  lv <- which(lv_opts == nLevels)
  nr <- which(nr_opts == nReps)

  offset_np <- (np - 1L) * length(lv_opts) * length(nr_opts) + 1L
  offset_lv <- (lv - 1L) * length(nr_opts)

  offset <- offset_np + offset_lv + nr - 1L
  
  this_max <- etamax[offset]

  ix <- 0:(length.out - 1L)
  (ix * this_max / 6) * (ix * this_max / 6)
}

etaSQ_2x2ww <- function(nPart, nReps, length.out = 7L) {
  eta_max <- c(0.2351495, 0.2259243, 0.1956559,
               0.1844661, 0.1725521, 0.1662758)
  
  etaSQ_mx <- sapply(eta_max,
                     \(.x) ((0:(length.out - 1L)) * .x / (length.out - 1L))^2)

  np_opts <- c(28L, 40L, 50L)
  nr_opts <- c(10L, 16L)
  
  if (!nPart %in% np_opts) {
    stop("'nPart' must be one of: ",
         paste(np_opts, collapse = ", "))
  }
  if (!nReps %in% nr_opts) {
    stop("'nReps' must be one of: ",
         paste(nr_opts, collapse = ", "))
  }

  np <- which(np_opts == nPart)
  nr <- which(nr_opts == nReps)

  offset_np <- (np - 1L) * length(nr_opts) + 1L
  offset <- offset_np + nr - 1L
  
  etaSQ_mx[, offset]
}
