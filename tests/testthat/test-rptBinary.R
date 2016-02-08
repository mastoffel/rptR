
context("rptBinary")

# Set a seed for reproducibility of the randomization 
set.seed(23)

nind = 80
nrep = 30 # a bit higher
latmu = 0
latbv = 0.3
latgv = 0.1
latrv = 0.2
indid = factor(rep(1:nind, each=nrep))
groid = factor(rep(1:nrep, nind))
latim = rep(rnorm(nind, 0, sqrt(latbv)), each=nrep)
latgm = rep(rnorm(nrep, 0, sqrt(latgv)), nind)
latvals = latmu + latim + latgm + rnorm(nind*nrep, 0, sqrt(latrv))
expvals = VGAM::logit(latvals, inverse = TRUE)
obsvals = rbinom(nind*nrep, 1, expvals)
beta0 = latmu
beta0 = VGAM::logit(mean(obsvals))
md = data.frame(obsvals, indid, groid)

R_est <- rptBinary(formula = obsvals ~ (1|indid), grname = c("indid"), 
        data = md, nboot = 0, link = "logit", npermut = 0, parallel = FALSE)

test_that("repeatability point estimate works", {
        expect_that(is.numeric(unlist(R_est$R)), is_true()) 
        expect_equal(R_est$R["R_org", ], 0.065, tolerance = 0.01)
        expect_equal(R_est$R["R_link", ], 0.078, tolerance = 0.01)
})

R_est <- rptPoisson(formula = obsvals ~ (1|indid) + (1|groid), grname = c("indid", "groid"), 
        data = md, nboot = 0, link = "log", npermut = 0, parallel = FALSE)

test_that("repeatability point estimate works for more than one group", {
        expect_that(is.numeric(unlist(R_est$R)), is_true()) 
        expect_equal(R_est$R["R_org", ]$indid, 0.0148, tolerance = 0.01)
        expect_equal(R_est$R["R_org", ]$groid, 0, tolerance = 0.01)
        expect_equal(R_est$R["R_link", ]$indid, 0.0263, tolerance = 0.01)
        expect_equal(R_est$R["R_link", ]$groid, 0, tolerance = 0.01)
})