options(tidyverse.quiet=TRUE)
library("tidyverse")
suppressMessages(requireNamespace("binom"))

## Agresti-Coull confidence intervals for the mean of a binomially distributed
## random variable
agresti_coull <- function(.dat, .x, .n, .c = .99) {
  .dat |>
    bind_cols(
      map2(.dat |> pull({{.x}}),
           .dat |> pull({{.n}}),
           function(.x1, .n1) {
             binom::binom.confint(.x1, .n1, .c,
                                  methods = "agresti-coull") |>
               select(lower, upper)
           }) |>
	bind_rows()
    )
}

## run a chi-square test for independence for a 2x2 matrix
chi_sq_test <- function(r1c1, r1c2,
                        r2c1, r2c2, alpha = .01) {
  result <- chisq.test(rbind(c(r1c1, r1c2),
                             c(r2c1, r2c2)))
  result$p.value < alpha
  }

chi_sq_test_p <- function(r1c1, r1c2,
                          r2c1, r2c2, alpha = .01) {
  result <- chisq.test(rbind(c(r1c1, r1c2),
                             c(r2c1, r2c2)))
  result$p.value
}

dat_pow <- read_rds(file.path("..", "results", "results-power.rds")) |>
  filter(analysis != "ANOVA") |>
  mutate(nmc = map_int(pvals, \(.x) length(.x)),
	 ncov = map_int(pvals, \(.x) sum(!is.na(.x))),
	 nsig = map_int(pvals, \(.x) sum(.x < .05, na.rm = TRUE)),
	 estr = fct_recode(estr,
			   "learning effect" = "exponential decay",
			   "random walk" = "gaussian random walk"),
	 case = case_when(analysis == "LMEM" & psr == "baseline" ~ "baseline",
			  analysis == "LMEM" & psr != "baseline" ~
                            paste0("PSR", substr(psr, 1, 1)),
			  analysis == "GAMM" & psr == "baseline" ~ "WRI",
			  analysis == "GAMM" & psr != "baseline" ~
                            paste0("WRI+PSR", substr(psr, 1, 1))),
	 power = nsig / nmc) |>
  select(-psr) |>
  select(nl, nr, estr, case, nmc, nsig, power) |>
  arrange(nl, nr, estr, case) |>
  agresti_coull(nsig, nmc)

dat_bl <- dat_pow |>
  filter(case == "baseline")

dat_sum <- dat_pow |>
  filter(case != "baseline") |>
  inner_join(dat_bl |>
	     select(nl, nr, estr, bl_power = power, bl_0 = lower, bl_1 = upper),
	     join_by(nl, nr, estr)) |>
  mutate(adv = 100 * (power - bl_power) / bl_power,
	 adv_0 = 100 * (lower - bl_1) / bl_1,
	 adv_1 = 100 * (upper - bl_0) / bl_0) |>
  select(-bl_power)

## stage 1:
## is the correction method better than baseline power?
dat_isEffective <- dat_pow |>
  mutate(notsig = nmc - nsig) |>
  filter(case != "baseline") |>
  inner_join(dat_bl |>
               mutate(bl_notsig = nmc - nsig) |>
               select(nl, nr, estr, bl_nsig = nsig, bl_notsig),
	     join_by(nl, nr, estr))  |>
  mutate(sigdiff_bl = pmap_lgl(list(nsig, notsig,
                                    bl_nsig, bl_notsig),
                               chi_sq_test),
         p_bl = pmap_dbl(list(nsig, notsig,
                              bl_nsig, bl_notsig),
                         chi_sq_test_p),
         dir = case_when(sigdiff_bl & (nsig > bl_nsig) ~ "better",
                         sigdiff_bl & (nsig < bl_nsig) ~ "worse"),
         adv = 100 * ((nsig / bl_nsig) - 1),
         diff = nsig - bl_nsig) |>
  select(-lower, -upper, -notsig, -bl_notsig)
                              
## stage 2:
## head-to-head comparison of PSR and WRI
dat_cmp <- dat_sum |>
  filter(case %in% c("WRI", "PSRm")) |>
  mutate(notsig = nmc - nsig) |>
  select(nl, nr, estr, case, nsig, notsig, power) |>
  pivot_wider(id_cols = c(nl, nr, estr),
              names_from = case,
              values_from = c(nsig, notsig, power)) |>
  mutate(sigdiff = pmap_lgl(list(notsig_PSRm, nsig_PSRm,
                                 notsig_WRI, nsig_WRI), chi_sq_test),
         pval = pmap_dbl(list(notsig_PSRm, nsig_PSRm,
                              notsig_WRI, nsig_WRI), chi_sq_test_p),
         winner = case_when(nsig_PSRm > nsig_WRI ~ "PSRm",
                            nsig_WRI > nsig_PSRm ~ "WRI",
                            .default = "tie"),
         win_by = case_when(nsig_WRI > nsig_PSRm ~ 100 * (nsig_WRI - nsig_PSRm) / nsig_PSRm,
                            nsig_PSRm > nsig_WRI ~ 100 * (nsig_PSRm - nsig_WRI) / nsig_WRI,
                            .default = 0)) |>
  arrange(winner, estr, desc(win_by))  |>
  select(-starts_with("nsig"), -starts_with("notsig"))
  
## stage 3:
## does combining both ever do better than the best of the two?
dat_combined_v_best <- dat_sum |>
  filter(case %in% c("WRI", "PSRm", "WRI+PSRm")) |>
  mutate(notsig = nmc - nsig) |>
  select(nl, nr, estr, case, nsig, notsig, power) |>
  pivot_wider(id_cols = c(nl, nr, estr),
              names_from = case,
              values_from = c(nsig, notsig, power)) |>
  mutate(best = case_when(nsig_PSRm > nsig_WRI ~ "PSRm",
                          nsig_WRI > nsig_PSRm ~ "WRI",
                          .default = "tie"),
         best_ns = case_when(best == "PSRm" ~ nsig_PSRm,
                             best == "WRI" ~ nsig_WRI,
                             .default = nsig_PSRm),
         best_nots = case_when(best == "PSRm" ~ notsig_PSRm,
                               best == "WRI" ~ notsig_WRI,
                               .default = notsig_PSRm),
         sigdiff_c = pmap_lgl(list(best_ns, best_nots,
                                   `nsig_WRI+PSRm`, `notsig_WRI+PSRm`),
                              chi_sq_test),
         winner2 = case_when(sigdiff_c & (`nsig_WRI+PSRm` > best_ns) ~ "combined",
                             sigdiff_c & (best_ns > `nsig_WRI+PSRm`) ~ best,
                             .default = NA_character_),
         win_by2 = 100 * (`nsig_WRI+PSRm` - best_ns) / best_ns)  |>
  select(nl, nr, estr, best, sigdiff_c, winner2, win_by2)

## dat_cmp |> arrange(desc(win_by)) |> print(n = +Inf)

## dat_combined_v_best |> arrange(winner2) |> print(n = +Inf)

## dat_isEffective |> arrange(desc(sigdiff_bl), estr, case, desc(adv), desc(diff)) |> print(n = +Inf)

.dat_psig <- dat_isEffective |>
  filter(case %in% c("PSRm", "WRI")) |>
  group_by(case) |>
  summarize(nsig = sum(sigdiff_bl),
            ntot = n(),
            perc_sig = 100 * nsig / ntot)

dat_casesum <- .dat_psig |>
  inner_join(dat_isEffective |>
             filter(sigdiff_bl, case %in% c("PSRm", "WRI")) |>
             group_by(case) |>
             summarize(min = min(adv),
                       median = median(adv),
                       max = max(adv),
                       .groups = "drop"),
             join_by(case))

.dat_psig_s <- dat_isEffective |>
  filter(case == "PSRs") |>
  summarize(nsig = sum(sigdiff_bl),
            ntot = n(),
            perc_sig = 100 * nsig / ntot)

dat_casesum_s <- .dat_psig_s |>
  bind_cols(dat_isEffective |>
            filter(sigdiff_bl, case == "PSRs") |>
            summarize(min = min(adv),
                      median = median(adv),
                      max = max(adv),
                      .groups = "drop"))

rm(.dat_psig, .dat_psig_s)

## number of cases (out of 36) where there was a significant difference
cmp_ndiff <- dat_cmp |>
  filter(sigdiff) |>
  nrow()

cmp_adv_wim <- dat_cmp |>
  filter(sigdiff, winner == "WRI") |>
  nrow()

cmp_adv_psr <- cmp_ndiff - cmp_adv_wim

nmc <- dat_pow |>
  pull(nmc) |>
  unique()
