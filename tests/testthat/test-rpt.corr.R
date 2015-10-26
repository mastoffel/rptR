context("rpt.corr")

# Set a seed for reproducibility of the randomization 
set.seed(23)
# simulation of a 3-level factor
alpha = 10    # Simulated (i.e. "true") causal effect of predictor on response
n = 98
n.per.group = 2
pred = factor(rep(1:(n/2),each=n.per.group))
group.dev = c(sample(seq(from = -0.4, to = 0.4, by = 0.1), size = n/2, replace = TRUE))
epsilon = rnorm(n, mean=0, sd=1)
resp = alpha + group.dev[as.numeric(pred)] + epsilon
md = data.frame(resp, pred)

test_that("NSE works", {
        expect_equal(rpt.corr(data = md, resp, pred)$R, 0.141, tolerance = 0.01) 
})

test_that("SE works", {
        expect_equal(rpt.corr_( data = md, "resp", "pred")$R, 0.141, tolerance = 0.01) 
})

test_that("Permutations work", {
        expect_equal(length(rpt.corr(data = md, y = resp, group = pred, npermut = 10)$R.permut), 10)
})

