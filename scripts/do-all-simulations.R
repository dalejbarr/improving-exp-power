options(tidyverse.quiet=TRUE)
library("tidyverse")
suppressPackageStartupMessages({
  library("mgcv")
  library("psrsim")
})

vers <- packageVersion("psrsim")
if (vers < "0.0.2.9000") { 
  stop("you are using an old version of psrsim, ", vers)
}
prefix <- "D"

g_basis_type <- "ps"

message("using psrsim version '", vers, "'")

seeds <- sample(1:.Machine$integer.max, 1L)
.ff <- as.integer(commandArgs(TRUE)[1])
.gg <- commandArgs(TRUE)[2]
if (is.na(.ff) || is.na(.gg)) {
    stop("need two command-line args:\n", "  no-monte-carlo-runs  output-path")
}
nmc <- .ff
outpath <- .gg

if (!dir.exists(outpath)) {
  dir.create(outpath)
}

errstrs <- c("expDecay", "pinkNoise", "randomWalk", "mixed")

basis_k <- read_rds(file.path("..", "results", "choose-basis-fn_targets-25.rds"))

## make sure we have the full data!
stopifnot(all(basis_k[["N"]] > 200L))

versions <- crossing(errstr = factor(errstrs, levels = errstrs)) |>
  crossing(etastep = 4L)

designs_pre <- crossing(tibble(np = 25L),
                        tibble(nl = c(2L, 4L, 8L)),
                        tibble(nr = c(6L, 12L, 24L)),
                        versions)  |>
  mutate(ns = 1L)

designs <- bind_rows(designs_pre,
                     designs_pre |>
                     mutate(ns = as.integer(nr / 2)),
                     designs_pre |>
                     mutate(ns = nr)) |>
  arrange(np, nl, nr, errstr, ns) |>
  select(np, nl, nr, ns, errstr, etastep)

# format string for progress reporting
fmt <- sprintf("(%%0%dd/%%0%dd)",
               nrow(designs) |> as.character() |> nchar(),
               nrow(designs) |> as.character() |> nchar())

todo <- designs |>
  inner_join(basis_k |> select(np, nl, nr, errstr, targ_k),
             join_by(np, nl, nr, errstr))

set.seed(seeds[1])
message(nmc, " Monte Carlo runs with ",
        "seed ", seeds[1], " and ",
        "np=", todo[["np"]][1], ", ",
        "nl=", todo[["nl"]][1], ", ",
        "nr=", todo[["nr"]][1], ", ",
        "ns=", todo[["ns"]][1])

outfile <- sprintf("%s%s_%010d_%05d_%02d_%d_%02d_%s.rds",
                   prefix,
                   Sys.info()["nodename"], seeds, nmc,
                   todo[["np"]][1],
                   todo[["nl"]][1],
                   todo[["nr"]][1],
                   Sys.time())

master_start_time <- Sys.time()
allres <- lapply(seq_len(nrow(todo)), \(.i) {
  r <- todo[.i, ]

  mstr <- sprintf("nl=%02d nr=%02d ns=%02d eta=%d errs=%s bs=%s k=%02d",
                  r$nl, r$nr, r$ns, r$etastep, r$errstr, g_basis_type,
                  r$targ_k)
  cat("  ", sprintf(fmt, .i, nrow(todo)), " ",
      mstr, "... ", sep = "")
  start_time <- Sys.time()
  
  res <- t(sapply(seq_len(nmc), \(.x) {
    dat <- simdata_1f(etastep = r$etastep, nPart = r$np,
                      nLevels = r$nl, nReps = r$nr,
                      nSubblocks = r$ns,
                      errstr = as.character(r$errstr)) |>
      mutate(id = factor(id),
             cond = factor(cond))

    pv_aov <- tryCatch(
      error = function(cnd) NA_real_,
      {
        dagg <- with(dat, aggregate(dv ~ cond + id, FUN = mean))
        .av <- summary(aov(dv ~ cond + Error(id), data = dagg))
        .av[["Error: Within"]][[1]][, "Pr(>F)"][1]
      })
    
    mod_lmm <- tryCatch(
      error = function(cnd) NULL,
      bam(dv ~ cond +
            s(id, bs = "re") +
            s(id, cond, bs = "re"),
          data = dat))
    
    pv_lmm <- if (!is.null(mod_lmm)) {
      tryCatch(
        error = function(cnd) NA_real_,
        anova(mod_lmm)[["pTerms.table"]][1, "p-value"]
      )
    } else {
      NA_real_
    }

    mod_gam <- 
      tryCatch(
        error = function(cnd) NULL,
        bam(dv ~ cond +
              s(id, cond, bs = "re") +
              s(order, id, bs = "fs", xt = list(bs = g_basis_type),
                m = 1, k = r$targ_k),
            data = dat)
      )
    
    pv_gam <- if (!is.null(mod_gam)) {
      tryCatch(
        error = function(cnd) NA_real_,
        anova(mod_gam)[["pTerms.table"]][1, "p-value"]
      )
    } else {
      NA_real_
    }
    
    c(p_aov = pv_aov, p_lmm = pv_lmm, p_gam = pv_gam)
  }))
  
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  cat(as.double(time_taken), attr(time_taken, "units"), "\n")
  res
})

message("finished! processing results")

result <- todo |>
  mutate(p_aov = map(allres, \(.x) .x[, "p_aov"]),
         p_lmm = map(allres, \(.x) .x[, "p_lmm"]),
         p_gam = map(allres, \(.x) .x[, "p_gam"]))

saveRDS(result, file = file.path(outpath, outfile))
message("results saved to '",
        file.path(outpath, outfile), "'")

master_diff <- Sys.time() - master_start_time
message("total time elapsed: ",
        as.double(master_diff), " ", attr(master_diff, "units"), "\n")

