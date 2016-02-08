
context("rptPoisson")

# Set a seed for reproducibility of the randomization 
set.seed(23)

nind = 40
nrep = 10
latmu = 0
latbv = 0.3
latgv = 0.3
latrv = 0.2
indid = factor(rep(1:nind, each=nrep))
groid = factor(rep(1:nrep, nind))
latim = rep(rnorm(nind, 0, sqrt(latbv)), each=nrep)
latgm = rep(rnorm(nrep, 0, sqrt(latgv)), nind)
latvals = latmu + latim + latgm + rnorm(nind*nrep, 0, sqrt(latrv))
expvals = exp(latvals)
obsvals = rpois(nind*nrep, expvals)
beta0 = latmu
beta0 = log(mean(obsvals))
md = data.frame(obsvals, indid, groid)

R_est <- rptPoisson(formula = obsvals ~ (1|indid), grname = c("indid"), 
        data = md, nboot = 0, link = "log", npermut = 0, parallel = FALSE)

test_that("repeatability point estimate works", {
        expect_that(is.numeric(unlist(R_est$R)), is_true()) 
        expect_equal(R_est$R["R_org", ], 0.147, tolerance = 0.01)
        expect_equal(R_est$R["R_link", ], 0.160, tolerance = 0.01)
})

R_est <- rptPoisson(formula = obsvals ~ (1|indid) + (1|groid), grname = c("indid", "groid"), 
        data = md, nboot = 0, link = "log", npermut = 0, parallel = FALSE)

test_that("repeatability point estimate works for more than one group", {
        expect_that(is.numeric(unlist(R_est$R)), is_true()) 
        expect_equal(R_est$R["R_org", ]$indid, 0.262, tolerance = 0.01)
        expect_equal(R_est$R["R_org", ]$groid, 0.343, tolerance = 0.01)
        expect_equal(R_est$R["R_link", ]$indid, 0.267, tolerance = 0.01)
        expect_equal(R_est$R["R_link", ]$groid, 0.333, tolerance = 0.01)
})

