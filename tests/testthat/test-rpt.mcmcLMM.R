
# rpt.mcmcLMM tests
context("rpt.mcmcLMM")
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
        expect_equal(is.numeric(rpt.mcmcLMM(data = md, y = resp, groups = pred)$R), TRUE)

})

test_that("SE works", {
        expect_equal(is.numeric(rpt.mcmcLMM_(data = md, y = "resp", groups = "pred")$R), TRUE)
})

