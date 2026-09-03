options(tidyverse.quiet=TRUE)
library("tidyverse")
suppressPackageStartupMessages({
  library("mgcv")
  library("psrsim")
})

g_bases <- c("ps") ## p-splines

vers <- packageVersion("psrsim")
if (vers < "0.0.1.9000") { ## old version was 0.0.0.9000
  stop("needs version '0.0.1.900' of psrsim or higher")
} else {
  prefix <- "B"
}
message("using psrsim version '", vers, "'")

outpath <- if (interactive()) {
             "results-choose-num-basis-full"
           } else {
             if (is.na(commandArgs(TRUE)[1])) {
               "need to supply output directory as command-line argument"
             } else {
               commandArgs(TRUE)[1]
             }
           }

if (!dir.exists(outpath)) stop("output directory '", outpath, "' does not exist")

errstrs <- c("expDecay", "pinkNoise", "randomWalk", "mixed")

versions <- tibble(errstr = factor(errstrs, levels = errstrs)) |>
  crossing(bfn = g_bases)

todo <- crossing(tibble(np = 25),
                 tibble(nl = c(2L, 4L, 8L)),
                 tibble(nr = c(6L, 12L, 24L)),
                 versions) |>
  mutate(ns = 1L)

while (1) {
  for (.i in seq_len(nrow(todo))) {
    seed <- sample(1:.Machine$integer.max, 1L)
    set.seed(seed)
    
    outfile <- sprintf("%s%d_%02d_%d_%02d_%s_%s.rds",
                       prefix, seed,
                       todo[["np"]][.i],
                       todo[["nl"]][.i],
                       todo[["nr"]][.i],
                       todo[["bfn"]][.i],
                       todo[["errstr"]][.i])

    n_trials <- todo[["nl"]][.i] * todo[["nr"]][.i]
    n_starting <- (n_trials <= 12L) * 6L +
      as.integer((1/3) * n_trials) * (n_trials > 12L)
        
    dat <- simdata_1f(etastep = 1L,
                      nPart = todo[["np"]][.i],
                      nLevels = todo[["nl"]][.i],
                      nReps = todo[["nr"]][.i],
                      nSubblocks = 1L,
                      errstr = as.character(todo[["errstr"]][.i])) |>
      mutate(id = factor(id),
             cond = factor(cond))

    res <- lapply(n_starting, \(.b) {
      cat(sprintf("fitting %s %02d-%d-%02d with %d bases (%s)... ",
                  todo[["errstr"]][.i],
                  todo[["np"]][.i], todo[["nl"]][.i],
                  todo[["nr"]][.i], .b,
                  todo[["bfn"]][.i]))
      t1 <- Sys.time()
      mod_bam <- tryCatch(
        error = function(cnd) NULL, {
          bam(dv ~ cond +
                s(id, order, bs = "fs",
                  xt = list(bs = todo[["bfn"]][.i]),
                  m = 1, k = .b) +
                s(id, cond, bs = "re"),
              data = dat)
        })
      t1_diff <- Sys.time() - t1

      if (!is.null(mod_bam)) {
        edf <- tryCatch(
          error = function(cnd) NA_real_, {
            summary(mod_bam)[["s.table"]]["s(id,order)", "edf"] /
              todo[["np"]][.i]
          })
      } else {
        NA_real_
      }
      
      cat("edf = ", edf, " ", format(t1_diff), "\n", sep = "")
      tibble(basis = .b,
             edf = edf,
             time = t1_diff)
    })

    saveRDS(bind_rows(res), file.path(outpath, outfile))
    message("wrote results to '", file.path(outpath, outfile), "'")
  }
}
