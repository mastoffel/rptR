# rpt.remlLMM tests

# Set a seed for reproducibility of the randomization 
set.seed(23)
# simulation of a 3-level factor
alpha = 10    # Simulated (i.e. "true") causal effect of predictor on response
n = 99
n.per.group = n/3
pred = factor(rep(LETTERS[1:3],each=n.per.group))
group.dev = c(-0.3, -0.1, +0.4)
epsilon = rnorm(n, mean=0, sd=1)
resp = alpha + group.dev[as.numeric(pred)] + epsilon
md = data.frame(resp, pred)

test_that("NSE works", {
        # non-standard eval
        expect_equal(rpt.remlLMM(data = md, y = resp, groups = pred, nboot = 0, npermut = 0)$R, 0.0995, tolerance = 0.01)
})

test_that("SE works", {
        # non-standard eval
        expect_equal(rpt.remlLMM_(data = md, y = "resp", groups = "pred", nboot = 0, npermut = 0)$R, 0.0995, tolerance = 0.01)
})

test_that("bootstrapping works", {
        # non-standard eval
        expect_equal(length(rpt.remlLMM(data = md, y = resp, groups = pred, nboot = 10, npermut = 0)$R.boot), 10)
        expect_equal(is.numeric(rpt.remlLMM(data = md, y = resp, groups = pred, nboot = 10, npermut = 0)$R.boot), TRUE)
})
test_that("permutation works", {
        # non-standard eval
        expect_equal(length(rpt.remlLMM(data = md, y = resp, groups = pred, nboot = 0, npermut = 10)$R.permut), 10)
        expect_equal(is.numeric(rpt.remlLMM(data = md, y = resp, groups = pred, nboot = 0, npermut = 10)$R.permut), TRUE)
})

