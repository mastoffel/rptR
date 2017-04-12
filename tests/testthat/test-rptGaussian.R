
context("rptGaussian")

# Set a seed for reproducibility of the randomization 
set.seed(23)

# load data
data(BeetlesBody)

########## checks for one random effect
# run with one random effect, no boot, no permut
R_est_1 <- rptGaussian(BodyL ~ (1|Population), grname="Population", data=BeetlesBody, nboot=0,
                       npermut=0)

test_that("rpt estimation works for one random effect, no boot, no permut, no parallelisation, logit link", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        expect_equal(as.numeric(R_est_1$R), 0.2985548, tolerance = 0.001)
 
})

test_that("LRT works", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        expect_equal(R_est_1$P$LRT_P,  9.59e-62, tolerance = 0.001)
})


# run with one random effect, boot, no permut
R_est_2 <-  rptGaussian(BodyL ~ (1|Population), grname="Population", data=BeetlesBody, nboot=2,
                        npermut=0)

test_that("rpt estimation works for one random effect, boot, no permut, no parallelisation, logit link", {
        expect_equal(as.numeric(R_est_2$CI_emp["2.5%"]), 0.1881897 , tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp["97.5%"]), 0.2463285, tolerance = 0.001)
})


# run with one random effect, no boot, permut
R_est_3 <-  rptGaussian(BodyL ~ (1|Population), grname="Population", data=BeetlesBody, nboot=0,
                        npermut=5)

test_that("rpt estimation works for one random effect, no boot, permut, no parallelisation, logit link", {
        # original scale
        expect_equal(R_est_3$P$P_permut, 0.2, tolerance = 0.001)
})



########## checks for two random effects

# run with two random effects, no boot, no permut
R_est_4 <- rptGaussian(BodyL ~ (1|Container) + (1|Population), grname=c("Container", "Population"), data=BeetlesBody, nboot=0,
                       npermut=0)

test_that("rpt estimation works for two random effect, no boot, no permut, no parallelisation, logit link", {
        expect_that(is.numeric(unlist(R_est_4$R)), is_true()) 
        # 1st random effect
        expect_equal(R_est_4$R[[1]], 0.478317, tolerance = 0.001)
        # 2nd random effect
        expect_equal(R_est_4$R[[2]],  0.2561723, tolerance = 0.001)
})

test_that("LRTs works", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        expect_equal(R_est_4$P[1, "LRT_P"], 2.292582e-138, tolerance = 0.001)
        expect_equal(R_est_4$P[2, "LRT_P"], 2.156259e-07, tolerance = 0.001)

})

# variance components sum up to one
R_est_5 <- rptGaussian(BodyL ~ (1|Container) + (1|Population), 
        grname=c("Container", "Population", "Overdispersion"), 
        data=BeetlesBody, nboot=0, npermut=0)

test_that("random effect components sum to up to one", {
        expect_equal(sum(R_est_5$R), 1)
})

# run with two random effects, boot, no permut
R_est_6 <- rptGaussian(BodyL ~ (1|Container) + (1|Population), grname=c("Container", "Population"),
                       data=BeetlesBody, nboot=2, npermut=0)

test_that("rpt estimation works for two random effect, boot, no permut, no parallelisation, logit link", {
        # 1st random effect
        expect_equal(as.numeric(R_est_6$CI_emp[1, "2.5%"]),   0.3837608, tolerance = 0.001)
        expect_equal(as.numeric(R_est_6$CI_emp[1, "97.5%"]), 0.5347381, tolerance = 0.001)
        # 2nd random effect
        expect_equal(as.numeric(R_est_6$CI_emp[2, "2.5%"]), 0.1587494, tolerance = 0.001)
        expect_equal(as.numeric(R_est_6$CI_emp[2, "97.5%"]), 0.3445943, tolerance = 0.001)
})



# run with two random effects, no boot, permut
R_est_7 <- rptGaussian(BodyL ~ (1|Container) + (1|Population), grname=c("Container", "Population"), data=BeetlesBody, nboot=0,
        npermut=5)

test_that("rpt estimation works for two random effect, no boot, permut, no parallelisation, logit link", {
        # original scale
        expect_equal(R_est_7$P$P_permut[1], 0.2, tolerance = 0.001)
        expect_equal(R_est_7$P$P_permut[2], 0.4, tolerance = 0.001)
        
})


# Run with Variance, Residual and Overdisp
R_est_8 <- rptGaussian(BodyL ~ (1|Container) + (1|Population), 
        grname=c("Container", "Population", "Residual", "Overdispersion"), data=BeetlesBody, nboot=3,
        npermut=3, ratio = FALSE)

test_that("Variance estimation works for two random effects with residual and overdispersion and boot and permut", {
        expect_that(is.numeric(unlist(R_est_8$R)), is_true()) 
        # 1st random effect
        expect_equal(R_est_8$R[[1]], 2.205629, tolerance = 0.001)
        # 2nd random effect
        expect_equal(R_est_8$R[[2]],  1.181269, tolerance = 0.001)
        # test bootstraps
        expect_equal(R_est_8$CI_emp["Container", "2.5%"], 1.772372, tolerance = 0.001)
        expect_equal(R_est_8$CI_emp["Population", "2.5%"],  1.296191, tolerance = 0.001)
        expect_equal(R_est_8$CI_emp["Residual", "2.5%"], 1.227603, tolerance = 0.001)
        expect_equal(R_est_8$CI_emp["Overdispersion", "2.5%"],  1.227603, tolerance = 0.001)
})



# test whether repeatabilities are equal for grouping factors independent of residual and overdispersion specification

R_est_9 <- rptGaussian(BodyL ~ (1|Container) + (1|Population), 
        grname=c("Container", "Population", "Residual", "Overdispersion"), data=BeetlesBody, nboot=0,
        npermut=0, ratio = TRUE)

R_est_10 <- rptGaussian(BodyL ~ (1|Container) + (1|Population), 
        grname=c("Container", "Population"), data=BeetlesBody, nboot=0,
        npermut=0, ratio = TRUE)
        
test_that("Repeatabilities are equal for grouping factors independent of residual and overdispersion specification", {
        expect_false(any(R_est_9$R[c(1,2)] == R_est_10$R) == FALSE)
})


# test wether estimates are independent of sequence

R_est_11 <- rptGaussian(BodyL ~ (1|Container) + (1|Population), 
        grname=c("Population", "Container"), data=BeetlesBody, nboot=0,
        npermut=0, ratio = TRUE)


test_that("Repeatabilities are equal when grname has a different grname sequence", {
        expect_equal(R_est_10$R$Container, R_est_11$R$Container)
        expect_equal(R_est_10$R$Population, R_est_11$R$Population)
})

test_that("LRTs are equal for different different sequence in grname argument", {
        expect_equal(R_est_10$P["Container", "LRT_P"], R_est_11$P["Container", "LRT_P"])
        expect_equal(R_est_10$P["Population", "LRT_P"], R_est_11$P["Population", "LRT_P"])
})



R_est_12 <- rptGaussian(BodyL ~ (1|Population) + (1|Container), 
        grname=c("Container", "Population"), data=BeetlesBody, nboot=0,
        npermut=0, ratio = TRUE)

test_that("Repeatabilities are equal when random effects in formula have a different order", {
        expect_equal(R_est_10$R$Container, R_est_12$R$Container)
        expect_equal(R_est_10$R$Population, R_est_12$R$Population)
})

test_that("LRTs are equal when random effects in formula have a different order", {
        expect_equal(R_est_10$P["Container", "LRT_P"], R_est_12$P["Container", "LRT_P"])
        expect_equal(R_est_10$P["Population", "LRT_P"], R_est_12$P["Population", "LRT_P"])
})

# random slopes
R_est_13 <- rptGaussian(BodyL ~ Treatment + Sex + (1 + Treatment|Population),  
            grname=c("Population"), data=BeetlesBody, nboot=0, npermut=0)

test_that("Random Slope repeatability point estimate is correct", {
        expect_equal(R_est_13$R$Population, 0.5015875, tolerance = 0.001)
})

test_that("Random slopes fitted without correlation (i.e. with grname occuring in multiple random effect terms) throw an error", {
        expect_error(rptGaussian(BodyL ~ Treatment + Sex + (1|Population) + (0 + Treatment|Population),  
                grname=c("Population"), data=BeetlesBody, nboot=0, npermut=0))
})


