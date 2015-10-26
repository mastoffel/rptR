# rpt.aov tests

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
        expect_equal(rpt.aov(data = md, y = resp, groups = pred,  npermut = 0)$R, 0.0995, tolerance = 0.01)
})

test_that("SE works", {
        expect_equal(rpt.aov_(data = md, y = "resp", groups = "pred",  npermut = 0)$R, 0.0995, tolerance = 0.01)
})

test_that("Permutation test works", {
        expect_equal(length(rpt.aov(data = md, y = resp, groups = pred,  npermut = 10)$R.permut), 10)
        expect_equal(is.numeric(rpt.aov(data = md, y = resp, groups = pred,  npermut = 10)$R.permut), TRUE)
})

# test_that("Parallelization works", {
#         expect_equal(length(rpt.aov(data = md, y = resp, groups = pred,  npermut = 10, parallel = TRUE)$R.permut), 10)
#         expect_equal(is.numeric(rpt.aov(data = md, y = resp, groups = pred,  npermut = 10)$R.permut), TRUE)
# })
