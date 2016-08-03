
context("rptBinary")

# load data
data(BeetlesMale)

# Set a seed for reproducibility of the randomization 
set.seed(23)


########## checks for one random effect

# run with one random effect, no boot, no permut
R_est_1 <- rptBinary(Colour ~ (1|Population),  grname=c("Population"), data=BeetlesMale, nboot=0, npermut=0)

test_that("rpt estimation works for one random effect, no boot, no permut, no parallelisation, logit link", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        expect_equal(R_est_1$R["R_org", ], 0.1858031, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", ], 0.2232935, tolerance = 0.001)
})

test_that("LRT works", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        expect_equal(R_est_1$P$LRT_P, 8.656602e-15, tolerance = 0.001)
})


# run with one random effect, boot, no permut
R_est_2 <- rptBinary(Colour ~ (1|Population),  grname=c("Population"), data=BeetlesMale, nboot=2, npermut=0)

test_that("rpt estimation works for one random effect, boot, no permut, no parallelisation, logit link", {
        
        expect_that(is.numeric(unlist(R_est_2$R)), is_true()) 
        expect_equal(R_est_2$R["R_org", ], 0.1858031, tolerance = 0.001)
        expect_equal(R_est_2$R["R_link", ], 0.2232935, tolerance = 0.001)
        
        # original scale
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["2.5%"]), 0.1793817, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["97.5%"]), 0.2019976, tolerance = 0.001)
        # link scale
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["2.5%"]), 0.2189531, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["97.5%"]), 0.2374219, tolerance = 0.001)
        
})



# run with one random effect, no boot, permut
R_est_3 <- rptBinary(Colour ~ (1|Population),  grname=c("Population"), data=BeetlesMale, nboot=0, npermut=5)

test_that("rpt estimation works for one random effect, no boot, permut, no parallelisation, logit link", {
        
        expect_that(is.numeric(unlist(R_est_3$R)), is_true()) 
        expect_equal(R_est_3$R["R_org", ], 0.1858031, tolerance = 0.001)
        expect_equal(R_est_3$R["R_link", ], 0.2232935, tolerance = 0.001)
        
        # original scale
        expect_equal(R_est_3$P$P_permut_org, 0.2, tolerance = 0.001)
        # link scale
        expect_equal(R_est_3$P$P_permut_link, 0.2, tolerance = 0.001)
        
})

# run with one random effect, 1 boot, no permut
R_est_4 <- rptBinary(Colour ~ (1|Population),  grname=c("Population"), data=BeetlesMale, nboot=0, npermut=1)

test_that("rpt estimation works for one random effect, no boot, permut, no parallelisation, logit link", {
        
        expect_that(is.numeric(unlist(R_est_3$R)), is_true()) 
        expect_equal(R_est_3$R["R_org", ], 0.1858031, tolerance = 0.001)
        expect_equal(R_est_3$R["R_link", ], 0.2232935, tolerance = 0.001)
        
        # original scale
        expect_equal(R_est_3$P$P_permut_org, 0.2, tolerance = 0.001)
        # link scale
        expect_equal(R_est_3$P$P_permut_link, 0.2, tolerance = 0.001)
        
})




########## checks for two random effects


