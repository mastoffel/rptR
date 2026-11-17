context("rptProportion series-level R with unequal numbers of trials")

# Unbalanced designs: the number of trials n_i differs across observations.
# Estimators of the original-scale R compared here:
#   scale_obs = "trial"      single-trial residual Ep(1-Ep) (default)
#   scale_obs = "series"     Ep(1-Ep) / n_h, n_h = harmonic mean of n_i
#   scale_obs = "series_obs" mean_i[p_i(1-p_i) / n_i], fitted p_i
#   scale_obs = "series", n_ref = 10   Ep(1-Ep) / 10
# True values come from numerical integration over the logit-normal
# individual effects (pi = plogis(beta0 + a), a ~ N(0, VA)):
#   Vb = Var(pi), W = E[pi (1 - pi)]
#   trial level:           Vb / (Vb + W)
#   series, actual n mix:  Vb / (Vb + W * mean(1 / n_i))
#   series of n_ref trials: Vb / (Vb + W / n_ref)
# Tolerances are about 2.5 SD of each estimator over 500 simulated data sets
# of the same design (see the PR description).

suppressWarnings(RNGversion("3.5.0"))

VA <- 1
beta0 <- -0.5
n_ref <- 10

true_R <- function(n) {
        f <- function(a, k) stats::plogis(beta0 + a)^k * stats::dnorm(a, 0, sqrt(VA))
        m1 <- stats::integrate(f, -Inf, Inf, k = 1)$value
        m2 <- stats::integrate(f, -Inf, Inf, k = 2)$value
        Vb <- m2 - m1^2
        W <- m1 - m2
        c(trial = Vb / (Vb + W), mix = Vb / (Vb + W * mean(1 / n)),
          ref = Vb / (Vb + W / n_ref))
}

sim_unbal <- function(n, seed) {
        set.seed(seed)
        ind <- factor(rep(1:60, each = 4))
        a <- stats::rnorm(60, 0, sqrt(VA))
        s <- stats::rbinom(240, n, stats::plogis(beta0 + a[as.integer(ind)]))
        data.frame(ind = ind, s = s, f = n - s)
}

R_org <- function(dat, ...) {
        suppressWarnings(rptProportion(cbind(s, f) ~ (1 | ind), grname = "ind",
                data = dat, nboot = 0, npermut = 0, ...))$R["R_org", "ind"]
}

estimate_all <- function(dat) {
        c(trial = R_org(dat), series = R_org(dat, scale_obs = "series"),
          series_obs = R_org(dat, scale_obs = "series_obs"),
          ref = R_org(dat, scale_obs = "series", n_ref = n_ref))
}

# scenario A: n_i uniform on 2..30
set.seed(101)
n_A <- sample(2:30, 240, replace = TRUE)
dat_A <- sim_unbal(n_A, seed = 102)
tru_A <- true_R(n_A)
est_A <- estimate_all(dat_A)

# scenario B: strongly skewed, most n_i = 2, a few n_i = 50
set.seed(201)
n_B <- ifelse(stats::runif(240) < 0.9, 2, 50)
dat_B <- sim_unbal(n_B, seed = 202)
tru_B <- true_R(n_B)
est_B <- estimate_all(dat_B)

test_that("unequal trials (2..30): each estimator recovers its own target", {
        expect_lt(abs(est_A[["trial"]] - tru_A[["trial"]]), 0.1)
        expect_lt(abs(est_A[["series"]] - tru_A[["mix"]]), 0.1)
        expect_lt(abs(est_A[["series_obs"]] - tru_A[["mix"]]), 0.1)
        expect_lt(abs(est_A[["ref"]] - tru_A[["ref"]]), 0.1)
        # the trial-level R is far below the series-level R
        expect_lt(est_A[["trial"]], tru_A[["mix"]] / 2)
})

test_that("strongly skewed trials (mostly 2, a few 50): each estimator recovers its own target", {
        expect_lt(abs(est_B[["trial"]] - tru_B[["trial"]]), 0.1)
        expect_lt(abs(est_B[["series"]] - tru_B[["mix"]]), 0.15)
        expect_lt(abs(est_B[["series_obs"]] - tru_B[["mix"]]), 0.15)
        expect_lt(abs(est_B[["ref"]] - tru_B[["ref"]]), 0.15)
        # with mostly 2 trials, R for the actual mix is far below R for 10 trials
        expect_gt(est_B[["ref"]] - est_B[["series"]], 0.15)
})

test_that("without fixed effects series_obs equals series (all p_i = Ep)", {
        expect_equal(est_A[["series_obs"]], est_A[["series"]], tolerance = 1e-10)
        expect_equal(est_B[["series_obs"]], est_B[["series"]], tolerance = 1e-10)
})

test_that("n_ref equal to the harmonic mean reproduces the default series value", {
        n_h <- 1 / mean(1 / n_A)
        expect_equal(R_org(dat_A, scale_obs = "series", n_ref = n_h), est_A[["series"]],
                tolerance = 1e-10)
})

test_that("series_obs differs from series when fitted proportions vary with fixed effects", {
        set.seed(301)
        dat_C <- sim_unbal(n_A, seed = 302)
        dat_C$x <- rep(c(0, 1), 120)
        dat_C$s <- stats::rbinom(240, n_A, stats::plogis(-0.5 + 2 * dat_C$x +
                stats::rnorm(60)[as.integer(dat_C$ind)]))
        dat_C$f <- n_A - dat_C$s
        fit <- function(...) suppressWarnings(rptProportion(cbind(s, f) ~ x + (1 | ind),
                grname = "ind", data = dat_C, nboot = 0, npermut = 0, ...))$R
        R_series <- fit(scale_obs = "series")
        R_obs <- fit(scale_obs = "series_obs")
        expect_false(isTRUE(all.equal(R_obs["R_org", "ind"], R_series["R_org", "ind"])))
        expect_equal(R_obs["R_link", "ind"], R_series["R_link", "ind"])
        for (e in c("latent", "liability")) {
                R_e <- fit(scale_obs = "series_obs", expect = e)["R_org", "ind"]
                expect_true(is.finite(R_e) && R_e > 0 && R_e < 1)
        }
})

test_that("n_ref is validated", {
        for (bad in list(0, -1, NA_real_, Inf, c(5, 10), "10")) {
                expect_error(R_org(dat_A, scale_obs = "series", n_ref = bad), "n_ref")
        }
        expect_error(R_org(dat_A, n_ref = 10), "scale_obs")
        expect_error(R_org(dat_A, scale_obs = "series_obs", n_ref = 10), "scale_obs")
})

test_that("rpt() passes n_ref and series_obs through", {
        fit_rpt <- function(...) suppressWarnings(rpt(cbind(s, f) ~ (1 | ind), grname = "ind",
                data = dat_A, datatype = "Proportion", nboot = 0, npermut = 0, ...))$R["R_org", "ind"]
        expect_equal(fit_rpt(scale_obs = "series", n_ref = n_ref), est_A[["ref"]])
        expect_equal(fit_rpt(scale_obs = "series_obs"), est_A[["series_obs"]])
})

test_that("bootstrap refits use the n_ref formula", {
        set.seed(2)
        fit_boot <- suppressWarnings(rptProportion(cbind(s, f) ~ (1 | ind), grname = "ind",
                data = dat_B, nboot = 3, npermut = 0, scale_obs = "series", n_ref = n_ref))
        expect_equal(fit_boot$R["R_org", "ind"], est_B[["ref"]])
        expect_true(all(unlist(fit_boot$R_boot_org) > 2 * est_B[["trial"]]))
})
