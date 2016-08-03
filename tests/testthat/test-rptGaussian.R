
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
        
        # original scale
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

# # run with one random effect, boot, permut
# R_est_4 <- rptPoisson(Egg ~ Treatment + (1|Container), grname=c("Container"), data = BeetlesFemale,
#                       nboot=2, npermut=2)
# 
# test_that("rpt estimation works for one random effect, boot, permut, no parallelisation, logit link", {
# 
#         # original scale
#         expect_equal(R_est_4$P$P_permut_org, 0.5, tolerance = 0.001)
#         # link scale
#         expect_equal(R_est_4$P$P_permut_link, 0.5, tolerance = 0.001)
#         
#         
# 
# })




########## checks for two random effects

# run with one random effect, no boot, no permut
R_est_1 <- rptGaussian(BodyL ~ (1|Container) + (1|Population), grname=c("Container", "Population"), data=BeetlesBody, nboot=0,
                       npermut=0)

test_that("rpt estimation works for two random effect, no boot, no permut, no parallelisation, logit link", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        # 1st random effect
        expect_equal(R_est_1$R[[1]], 0.478317, tolerance = 0.001)
        # 2nd random effect
        expect_equal(R_est_1$R[[2]],  0.2561723, tolerance = 0.001)
})

test_that("LRTs works", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        expect_equal(R_est_1$P[1, "LRT_P"], 2.292582e-138, tolerance = 0.001)
        expect_equal(R_est_1$P[2, "LRT_P"], 2.156259e-07, tolerance = 0.001)
})


# run with one random effect, boot, no permut
R_est_2 <- rptGaussian(BodyL ~ (1|Container) + (1|Population), grname=c("Container", "Population"), data=BeetlesBody, nboot=2,
                       npermut=0)

test_that("rpt estimation works for two random effect, boot, no permut, no parallelisation, logit link", {
        
        # 1st random effect
        expect_equal(as.numeric(R_est_2$CI_emp[1, "2.5%"]),   0.3837608, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp[1, "97.5%"]), 0.5347381, tolerance = 0.001)
        # 2nd random effect
        expect_equal(as.numeric(R_est_2$CI_emp[2, "2.5%"]), 0.1587494, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp[2, "97.5%"]), 0.3445943, tolerance = 0.001)
        
})



# run with one random effect, no boot, permut
R_est_3 <- rptGaussian(BodyL ~ (1|Container) + (1|Population), grname=c("Container", "Population"), data=BeetlesBody, nboot=0,
        npermut=5)

test_that("rpt estimation works for two random effect, no boot, permut, no parallelisation, logit link", {
        # original scale
        expect_equal(R_est_3$P$P_permut[1], 0.2, tolerance = 0.001)
        expect_equal(R_est_3$P$P_permut[2], 0.4, tolerance = 0.001)
        
})

