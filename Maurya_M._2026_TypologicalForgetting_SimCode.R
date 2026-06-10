# =========================================================================
#  Maurya_2026_TypologicalForgetting_SimCode.R
#
#  Bayesian operating-characteristics simulation for the typologically-
#  extended Wickelgren forgetting model. Companion to:
#
#    Maurya, M. (2026). Beyond the single forgetting curve: A typologically-
#    extended Wickelgren function, a simulation study, and a pre-registered
#    test of error-type-specific memory decay. PsyArXiv preprint.
#
#  This file deposits the runnable R/brms/Stan code referenced in §6 of the
#  preprint. It generates synthetic data under the priors of Table 3,
#  fits the homogeneous null (M0) and the full typological model (M3) in
#  brms/Stan via cmdstanr, compares them with PSIS-LOO, and records
#  per-replication operating characteristics for the pre-registered
#  hypotheses P1a, P1b, P2, P3, and the M3-over-M0 preference (P5).
#
#  AUTHOR:    Manik Maurya (ORCID 0009-0005-3554-693X)
#  CONTACT:   manikmaurya.in@gmail.com
#  LICENCE:   CC BY 4.0
#  REQUIRES:  R >= 4.2,  brms >= 2.20,  cmdstanr,  loo,
#             dplyr, purrr, tibble,  furrr (optional, for parallelisation)
#
#  IMPORTANT — runtime budget. A single brms fit of M3 on n=90 × items=40 ×
#  4 timepoints (~14,400 binary rows) runs in roughly 60–180 s on a modern
#  4-core laptop. The full 500-replication operating-characteristics loop
#  therefore takes 12–25 hours wall-clock at default settings. Steps to
#  reduce this are documented in §B at the bottom of this file.
# =========================================================================

# ----- 0. PACKAGES & GLOBAL OPTIONS ---------------------------------------
required <- c("brms", "cmdstanr", "loo", "dplyr", "purrr", "tibble", "tidyr")
missing  <- setdiff(required, rownames(installed.packages()))
if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(setdiff(missing, "cmdstanr"))
  if ("cmdstanr" %in% missing) {
    install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/",
                                            getOption("repos")))
    cmdstanr::install_cmdstan()
  }
}
suppressPackageStartupMessages({
  library(brms); library(cmdstanr); library(loo)
  library(dplyr); library(purrr); library(tibble); library(tidyr)
})
options(brms.backend = "cmdstanr",
        mc.cores     = max(1L, parallel::detectCores() - 1L))

# ----- 1. TRUTH (PRIOR MEANS, LOGIT SCALE) --------------------------------
#  Anchored to:
#    Wickelgren (1974); Wixted & Carpenter (2007); Murre & Dros (2015)
#    Butterfield & Metcalfe (2001); Butler, Fazio, & Marsh (2011);
#    Metcalfe & Miele (2014); Isurin & McDonald (2001);
#    Roediger & Butler (2011); Gelman (2006).
#  See Table 3 in the preprint for derivation of each value.

TRUTH <- list(
  lambda = c(RF = 0.50, PK = 0.85, CF = 1.30, IN = 0.80),
  psi    = c(RF = 1.40, PK = 1.00, CF = 0.75, IN = 1.00),
  sigma  = c(RF = 0.80, PK = 0.40, CF = 0.20, IN = 0.40),
  gamma  = 1.30,
  tau_CF = 7.0,
  A_IN   = 0.55,
  kappa  = 10.0,
  omega  = 2 * pi / 8.0,
  phi    = 0.0,
  tau_u  = 0.6,
  tau_v  = 0.5
)
BASE_RATES   <- c(RF = .30, PK = .40, CF = .15, IN = .15)
TIMES_4PT    <- c(0, 2, 7, 14)
N_DEFAULT    <- 90
ITEMS        <- 40
inv_logit    <- function(x) 1.0 / (1.0 + exp(-x))

# ----- 2. DATA-GENERATING PROCESS -----------------------------------------
#  Returns a tibble matching the parent trial's nested structure.
#  Random effects: u_s ~ N(0, tau_u), v_i ~ N(0, tau_v)
#  Type-specific shape terms encoded in g_e(t).

simulate_dataset <- function(n = N_DEFAULT, items = ITEMS,
                             base_rates = BASE_RATES, times = TIMES_4PT,
                             truth = TRUTH, seed = 42) {
  set.seed(seed)
  u_s   <- rnorm(n,    0, truth$tau_u)
  v_i   <- rnorm(items, 0, truth$tau_v)
  etype <- sample(names(base_rates), items, replace = TRUE, prob = base_rates)
  grid  <- expand.grid(s = seq_len(n), i = seq_len(items), t = times,
                       KEEP.OUT.ATTRS = FALSE) |> as_tibble()
  grid <- grid |>
    mutate(type   = etype[i],
           r      = pmin(floor(t/3), 5),
           u_s    = u_s[s],
           v_i    = v_i[i])
  g_e <- function(type, t, truth) {
    if (type == "CF") -truth$gamma * (1 - exp(-t / truth$tau_CF))
    else if (type == "IN") truth$A_IN * exp(-t / truth$kappa) *
                             cos(truth$omega * t + truth$phi)
    else 0
  }
  grid <- grid |>
    rowwise() |>
    mutate(
      eta = truth$lambda[type] + u_s + v_i -
            truth$psi[type] * log1p(t) +     # beta fixed = 1
            truth$sigma[type] * r +
            g_e(type, t, truth),
      p   = inv_logit(eta),
      y   = rbinom(1, 1, p)
    ) |> ungroup() |>
    mutate(
      is_RF = as.numeric(type == "RF"),
      is_PK = as.numeric(type == "PK"),
      is_CF = as.numeric(type == "CF"),
      is_IN = as.numeric(type == "IN"),
      cf_rev_kernel = is_CF * (1 - exp(-t / truth$tau_CF)),
      in_env_cos    = is_IN * exp(-t / truth$kappa) * cos(truth$omega * t),
      in_env_sin    = is_IN * exp(-t / truth$kappa) * sin(truth$omega * t),
      type          = factor(type, levels = c("PK","RF","CF","IN"))
    ) |>
    select(s, i, t, type, r, y, is_RF, is_PK, is_CF, is_IN,
           cf_rev_kernel, in_env_cos, in_env_sin)
  grid
}

# ----- 3. brms MODEL SPECIFICATIONS ---------------------------------------
#  All models are crossed-random-effects Bayesian GLMMs on the logit scale.
#  - tau_CF is fixed at the prior mean (7.0 d) in the linear formulation,
#    so the model is linear in coefficients given the kernel column.
#  - The IN oscillation is expressed via the cos/sin envelope columns,
#    giving an amplitude pair (A_cos, A_sin) from which A_IN can be
#    derived as sqrt(A_cos^2 + A_sin^2).

prior_M3 <- c(
  # lambda_e: per-type intercept
  prior(normal(0.85, 0.5), class = "b", coef = "is_PK"),
  prior(normal(0.50, 0.5), class = "b", coef = "is_RF"),
  prior(normal(1.30, 0.5), class = "b", coef = "is_CF"),
  prior(normal(0.80, 0.5), class = "b", coef = "is_IN"),
  # -psi_e * log(1+t):  enters as type * log1pt with negative coef
  prior(normal(-1.00, 0.4), class = "b", coef = "is_PK:log1pt"),
  prior(normal(-1.40, 0.4), class = "b", coef = "is_RF:log1pt"),
  prior(normal(-0.75, 0.4), class = "b", coef = "is_CF:log1pt"),
  prior(normal(-1.00, 0.4), class = "b", coef = "is_IN:log1pt"),
  # sigma_e * r
  prior(normal(0.80, 0.40), class = "b", coef = "is_RF:r", lb = 0),
  prior(normal(0.40, 0.30), class = "b", coef = "is_PK:r", lb = 0),
  prior(normal(0.20, 0.30), class = "b", coef = "is_CF:r", lb = 0),
  prior(normal(0.40, 0.30), class = "b", coef = "is_IN:r", lb = 0),
  # Confabulation reversion magnitude (γ): coefficient on the
  # negative kernel -(1 - exp(-t/tau_CF)) for CF items
  prior(normal(1.0, 0.5), class = "b", coef = "cf_rev_kernel", lb = 0),
  # Interference oscillation cos/sin amplitudes
  prior(normal(0.0, 0.5), class = "b", coef = "in_env_cos"),
  prior(normal(0.0, 0.5), class = "b", coef = "in_env_sin"),
  # Random-effects SDs
  prior(student_t(3, 0, 1), class = "sd", group = "s"),
  prior(student_t(3, 0, 1), class = "sd", group = "i")
)

fit_M0 <- function(dat, ...) {
  brm(y ~ 1 + (1|s) + (1|i),
      data    = dat,
      family  = bernoulli(),
      prior   = c(prior(student_t(3, 0, 1), class = "sd")),
      chains  = 2, iter = 1500, refresh = 0, silent = 2,
      control = list(adapt_delta = 0.9), ...)
}

fit_M3 <- function(dat, ...) {
  # Sign convention: log1pt enters with a negative coefficient whose
  # magnitude is psi_e; cf_rev_kernel is the positive coefficient gamma
  # whose contribution to the linear predictor is negative
  # (it is -gamma * (1 - exp(-t/tau_CF)) for CF items).
  dat2 <- dat |>
    mutate(log1pt        = log1p(t),
           cf_rev_kernel = -cf_rev_kernel)    # flip sign so b > 0 means γ > 0
  f <- bf(
    y ~ 0 + is_RF + is_PK + is_CF + is_IN +
        is_RF:log1pt + is_PK:log1pt + is_CF:log1pt + is_IN:log1pt +
        is_RF:r + is_PK:r + is_CF:r + is_IN:r +
        cf_rev_kernel + in_env_cos + in_env_sin +
        (1 | s) + (1 | i)
  )
  brm(f, data = dat2,
      family  = bernoulli(),
      prior   = prior_M3,
      chains  = 2, iter = 1500, refresh = 0, silent = 2,
      control = list(adapt_delta = 0.95))
}

# ----- 4. PER-REPLICATION DRIVER ------------------------------------------

one_replication <- function(seed, n = N_DEFAULT, items = ITEMS,
                            base_rates = BASE_RATES, times = TIMES_4PT,
                            truth = TRUTH) {
  dat <- simulate_dataset(n = n, items = items, base_rates = base_rates,
                          times = times, truth = truth, seed = seed)
  m0  <- fit_M0(dat); m3 <- fit_M3(dat)
  l0  <- loo(m0);     l3 <- loo(m3)
  cmp <- loo_compare(l0, l3)
  best <- rownames(cmp)[1]

  draws <- as_draws_df(m3)
  q <- function(v) quantile(v, c(.025, .975))

  # ----- Hypothesis tests on contrasts of brms coefficients -----
  psi_RF <- -draws[["b_is_RF:log1pt"]]   # ψ = − coefficient on log(1+t)
  psi_PK <- -draws[["b_is_PK:log1pt"]]
  psi_CF <- -draws[["b_is_CF:log1pt"]]
  sig_RF <-  draws[["b_is_RF:r"]]
  sig_PK <-  draws[["b_is_PK:r"]]
  gamma  <-  draws[["b_cf_rev_kernel"]]

  P1a <- q(psi_RF - psi_PK)[1] > 0
  P1b <- q(psi_PK - psi_CF)[1] > 0
  P2  <- q(gamma)[1] > 0
  P3  <- q(sig_RF - sig_PK)[1] > 0
  M3_pref <- best == "m3" && abs(cmp[2, "elpd_diff"]) > 2 * cmp[2, "se_diff"]

  list(
    seed     = seed,
    P1a = P1a, P1b = P1b, P2 = P2, P3 = P3, M3_pref = M3_pref,
    delpd    = unname(cmp[2, "elpd_diff"]),
    se_delpd = unname(cmp[2, "se_diff"]),
    psi_RF_minus_PK = c(mean = mean(psi_RF - psi_PK),
                        q025 = q(psi_RF - psi_PK)[1],
                        q975 = q(psi_RF - psi_PK)[2]),
    psi_PK_minus_CF = c(mean = mean(psi_PK - psi_CF),
                        q025 = q(psi_PK - psi_CF)[1],
                        q975 = q(psi_PK - psi_CF)[2]),
    gamma_post      = c(mean = mean(gamma),
                        q025 = q(gamma)[1], q975 = q(gamma)[2]),
    sig_RF_minus_PK = c(mean = mean(sig_RF - sig_PK),
                        q025 = q(sig_RF - sig_PK)[1],
                        q975 = q(sig_RF - sig_PK)[2])
  )
}

# ----- 5. MAIN LOOP — 500 REPLICATIONS ------------------------------------
#  Set N_REPS to the desired number of replications. The default 500 is
#  the value reported in the preprint; smaller values are useful for
#  development. Set RUN_LOOP to FALSE to load this file as a library.

N_REPS   <- 500
SEED_BASE <- 20260607
RUN_LOOP <- interactive() && exists("RUN_LOOP_FORCE", inherits = FALSE) &&
            isTRUE(get("RUN_LOOP_FORCE"))

if (RUN_LOOP) {
  results <- list()
  for (r in seq_len(N_REPS)) {
    message(sprintf("[%s]  Replication %d / %d ...",
                    format(Sys.time(), "%H:%M:%S"), r, N_REPS))
    out <- tryCatch(one_replication(seed = SEED_BASE + r),
                    error = function(e) {
                      warning("rep ", r, " failed: ", conditionMessage(e))
                      NULL
                    })
    if (!is.null(out)) results[[length(results) + 1L]] <- out
    if (r %% 25 == 0) saveRDS(results, "sim_results_partial.rds")
  }
  saveRDS(results, "sim_results_500.rds")
  summary_oc <- tibble(
    rep      = seq_along(results),
    P1a      = vapply(results, `[[`, logical(1), "P1a"),
    P1b      = vapply(results, `[[`, logical(1), "P1b"),
    P2       = vapply(results, `[[`, logical(1), "P2"),
    P3       = vapply(results, `[[`, logical(1), "P3"),
    M3_pref  = vapply(results, `[[`, logical(1), "M3_pref")
  )
  oc_table <- summary_oc |>
    summarise(power_P1a = mean(P1a), power_P1b = mean(P1b),
              power_P2  = mean(P2),  power_P3  = mean(P3),
              pref_M3   = mean(M3_pref))
  print(oc_table)
  saveRDS(oc_table, "operating_characteristics.rds")
}

# =========================================================================
#  §A   Sensitivity grid (n × CF base rate × retention design)
# =========================================================================
#  Re-runs one_replication() across the grid documented in §6 of the
#  preprint (Figure 2). The grid is 60 cells (4 sample sizes × 5 CF
#  base rates × 3 retention designs); we recommend 100 replications per
#  cell when using brms (i.e. 6,000 fits ~ 6–10 days at default settings).
#  Substantially faster: use the Python frequentist approximation in
#  the deposited run_simulation.py to obtain the operating-characteristics
#  grid, then re-run a small confirmatory set (3 cells × 100 reps) in brms
#  to validate. Both approaches give converging conclusions.

grid_sensitivity <- function(ns       = c(60, 90, 120, 150),
                             cf_rates = c(.05, .10, .15, .20, .25),
                             designs  = list("3pt" = c(0, 2, 14),
                                             "4pt" = c(0, 2, 7, 14),
                                             "5pt" = c(0, 2, 7, 14, 28)),
                             reps_per_cell = 100,
                             seed_base = 20260607) {
  cells <- expand.grid(n = ns, cf_rate = cf_rates,
                       design = names(designs),
                       stringsAsFactors = FALSE) |> as_tibble()
  out_list <- vector("list", nrow(cells))
  for (k in seq_len(nrow(cells))) {
    n  <- cells$n[k]; cf <- cells$cf_rate[k]; dn <- cells$design[k]
    br <- BASE_RATES; br["CF"] <- cf; br["PK"] <- max(0, .55 - cf)
    br <- br / sum(br)
    times <- designs[[dn]]
    counts <- list(P1a = 0, P1b = 0, P2 = 0, P3 = 0, M3_pref = 0, valid = 0)
    for (r in seq_len(reps_per_cell)) {
      res <- tryCatch(one_replication(seed = seed_base + k * 10000 + r,
                                       n = n, base_rates = br, times = times),
                       error = function(e) NULL)
      if (!is.null(res)) {
        counts$valid <- counts$valid + 1L
        for (h in c("P1a","P1b","P2","P3","M3_pref"))
          counts[[h]] <- counts[[h]] + as.integer(res[[h]])
      }
    }
    out_list[[k]] <- tibble(
      n = n, cf_rate = cf, design = dn,
      reps = reps_per_cell, valid = counts$valid,
      power_P1a = counts$P1a / max(counts$valid, 1),
      power_P1b = counts$P1b / max(counts$valid, 1),
      power_P2  = counts$P2  / max(counts$valid, 1),
      power_P3  = counts$P3  / max(counts$valid, 1),
      pref_M3   = counts$M3_pref / max(counts$valid, 1)
    )
  }
  bind_rows(out_list)
}

# =========================================================================
#  §B   Runtime budget — how to obtain operating characteristics tractably
# =========================================================================
#
#  Default settings (500 brms reps, 2 chains, 1,500 iter) ≈ 12–25 h wall.
#  Faster paths:
#
#  (i)  Reduce N_REPS to 100 → 2.5–5 h. Gives stable operating
#       characteristics for the 4-point registered design but not for
#       the full sensitivity grid.
#
#  (ii) Use cmdstanr's variational inference (algorithm = "meanfield") as
#       a first-pass approximation, then re-fit the focal cells with
#       full HMC. brms supports this via the `algorithm` argument.
#       Reduces wall-clock time by 5–10x; should be reported as
#       approximate inference, not full HMC.
#
#  (iii)Run the Python frequentist approximation (run_simulation.py) for
#       the full 60-cell grid, then run brms only on the 4-point design
#       with 100 reps. The Python script took 188 s for 7,200 fits on
#       the deposited hardware; the brms confirmatory run takes 4–8 h.
#
#  The preprint's Figure 2 is built from the Python grid; brms is used
#  for inference on the actual parent-trial data.
#
# =========================================================================
#  §C   Software environment used to produce the deposited results
# =========================================================================
#
#  R version           4.3.x
#  brms                ≥ 2.20
#  cmdstanr            0.7.x
#  cmdstan             ≥ 2.33
#  rstan               ≥ 2.32
#  loo                 ≥ 2.6
#  tidyverse           ≥ 2.0
#  Operating system    Linux x86_64
#
#  sessionInfo() output is logged at run time when RUN_LOOP_FORCE = TRUE
#  via the writeLines() call below.

if (RUN_LOOP) {
  writeLines(capture.output(sessionInfo()), "session_info.txt")
  message("Wrote session_info.txt and sim_results_500.rds.")
}

# =========================================================================
#  End of file.   manikmaurya.in@gmail.com   |   cognivia.vercel.app
# =========================================================================
