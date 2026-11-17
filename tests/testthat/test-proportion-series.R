context("rptProportion series-level original scale")

# Regression test for the original-scale repeatability of binomial
# (proportion) data with n > 1 trials per observation.
#
# With the default scale_obs = "trial", R_org uses the single-trial residual
# variance Ep * (1 - Ep), i.e. it is the repeatability of single binary
# decisions. The repeatability of the observed series-level proportions
# (successes / trials) is much larger when each observation has many trials.
# scale_obs = "series" divides the binomial residual variance by the
# harmonic mean of the number of trials per observation.

suppressWarnings(RNGversion("3.5.0"))

sim_series <- function(n_ind, n_rep, n_trials, VA, beta0, seed) {
        set.seed(seed)
        ind <- factor(rep(seq_len(n_ind), each = n_rep))
        a <- stats::rnorm(n_ind, 0, sqrt(VA))
        p <- stats::plogis(beta0 + a[as.integer(ind)])
        succ <- stats::rbinom(length(p), n_trials, p)
        data.frame(ind = ind, succ = succ, fail = n_trials - succ)
}

# True repeatability of the observed proportions y / n (series level),
# obtained by numerical integration over the logit-normal individual effects:
# R = Var(pi_i) / (Var(pi_i) + E[pi_i (1 - pi_i)] / n)
true_series_R <- function(VA, beta0, n_trials) {
        f <- function(a, g) g(stats::plogis(beta0 + a)) * stats::dnorm(a, 0, sqrt(VA))
        m1 <- stats::integrate(f, -Inf, Inf, g = function(p) p)$value
        m2 <- stats::integrate(f, -Inf, Inf, g = function(p) p^2)$value
        var_between <- m2 - m1^2
        var_within <- (m1 - m2) / n_trials
        var_between / (var_between + var_within)
}

VA <- 1
beta0 <- 0
n_trials <- 10
R_true <- true_series_R(VA, beta0, n_trials)

dat <- sim_series(n_ind = 60, n_rep = 4, n_trials = n_trials, VA = VA,
                  beta0 = beta0, seed = 2026)

fit_trial <- suppressWarnings(rptProportion(cbind(succ, fail) ~ (1 | ind),
        grname = "ind", data = dat, nboot = 0, npermut = 0))
fit_series <- suppressWarnings(rptProportion(cbind(succ, fail) ~ (1 | ind),
        grname = "ind", data = dat, nboot = 0, npermut = 0, scale_obs = "series"))

R_trial_org <- fit_trial$R["R_org", "ind"]
R_series_org <- fit_series$R["R_org", "ind"]
# Spearman-Brown step-up of the single-trial R to a series of n trials
R_sb <- n_trials * R_trial_org / (1 + (n_trials - 1) * R_trial_org)

test_that("default (trial) original-scale R underestimates the series-level R", {
        expect_true(R_true > 0.6)
        expect_lt(R_trial_org, R_true / 2)
})

test_that("scale_obs = 'series' recovers the true series-level R", {
        expect_equal(R_series_org, R_true, tolerance = 0.1)
        # roughly the Spearman-Brown step-up of the trial-level R
        expect_equal(R_series_org, R_sb, tolerance = 0.1)
})

test_that("scale_obs = 'series' leaves the link-scale R unchanged", {
        expect_equal(fit_series$R["R_link", "ind"], fit_trial$R["R_link", "ind"])
})

test_that("default output is unchanged by the new argument", {
        fit_default_explicit <- suppressWarnings(rptProportion(cbind(succ, fail) ~ (1 | ind),
                grname = "ind", data = dat, nboot = 0, npermut = 0, scale_obs = "trial"))
        expect_equal(fit_default_explicit$R, fit_trial$R)
})

test_that("rpt() passes scale_obs through for datatype = 'Proportion'", {
        fit_rpt <- suppressWarnings(rpt(cbind(succ, fail) ~ (1 | ind), grname = "ind",
                data = dat, datatype = "Proportion", nboot = 0, npermut = 0,
                scale_obs = "series"))
        expect_equal(fit_rpt$R, fit_series$R)
})

test_that("bootstrap samples use the series-level formula", {
        set.seed(1)
        fit_boot <- suppressWarnings(rptProportion(cbind(succ, fail) ~ (1 | ind),
                grname = "ind", data = dat, nboot = 3, npermut = 0, scale_obs = "series"))
        expect_true(all(unlist(fit_boot$R_boot_org) > 2 * R_trial_org))
})
