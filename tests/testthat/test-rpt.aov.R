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


test_that("all input formats work and point estimates for repeatability correct", {
        # non-standard eval
        expect_equal(rpt.aov(y = resp, groups = pred, data = md, npermut = 0)$R, 0.0995, tolerance = 0.01)
        expect_equal(rpt.aov(resp, pred, md, npermut = 0)$R, 0.0995, tolerance = 0.01)
        # standard eval
        expect_equal(rpt.aov(y = "resp", groups = "pred", data = md, npermut = 0)$R, 0.0995, tolerance = 0.01)
        expect_equal(rpt.aov("resp", "pred", md, npermut = 0)$R, 0.0995, tolerance = 0.01)
        # no df
        attach(md)
        expect_equal(rpt.aov(resp, pred, npermut = 0)$R, 0.0995, tolerance = 0.01)
        detach(md)
})

