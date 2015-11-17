context("rpt.binomGLMM.add")

# Set a seed for reproducibility of the randomization 
set.seed(23)
# simulation of a 3-level factor
alpha = 10    # Simulated (i.e. "true") causal effect of predictor on response
n = 100
n.per.group = 2
pred = factor(sample(rep(1:n/2, each=n.per.group)))
resp= rbinom(n, 1, 0.3)
# resp = alpha + group.dev[as.numeric(pred)] + epsilon
md = data.frame(resp, pred)


test_that("NSE works", {
        expect_equal(is.numeric(rpt.binomGLMM.add(data = md, resp, pred, CI = 0.95, prior = NULL,
                                                  verbose = FALSE)$R.link), TRUE)
        expect_equal(is.numeric(rpt.binomGLMM.add(data = md, resp, pred, CI = 0.95, prior = NULL,
                                                  verbose = FALSE)$R.org), TRUE)
        
})

test_that("SE works", {
        expect_equal(is.numeric(rpt.binomGLMM.add_(data = md, "resp", "pred", CI = 0.95, prior = NULL,
                                                  verbose = FALSE)$R.link), TRUE)
        expect_equal(is.numeric(rpt.binomGLMM.add_(data = md, "resp", "pred", CI = 0.95, prior = NULL,
                                                  verbose = FALSE)$R.org), TRUE)
        
})