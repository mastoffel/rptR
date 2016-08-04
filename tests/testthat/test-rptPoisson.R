
context("rptPoisson")

# Set a seed for reproducibility of the randomization 
set.seed(23)

# load data
data(BeetlesFemale)

########## checks for one random effect

# run with one random effect, no boot, no permut
R_est_1 <- rptPoisson(Egg ~ Treatment + (1|Container), grname=c("Container"), data = BeetlesFemale,
                      nboot=0, npermut=0)

test_that("rpt estimation works for one random effect, no boot, no permut, no parallelisation, logit link", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        expect_equal(R_est_1$R["R_org", ], 0.505, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", ], 0.482, tolerance = 0.001)
})

test_that("LRT works", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        expect_equal(R_est_1$P$LRT_P, 1.98e-46, tolerance = 0.001)
})


# run with one random effect, boot, no permut
R_est_2 <- rptPoisson(Egg ~ Treatment + (1|Container), grname=c("Container"), data = BeetlesFemale,
                      nboot=2, npermut=0)

test_that("rpt estimation works for one random effect, boot, no permut, no parallelisation, logit link", {
        
        expect_that(is.numeric(unlist(R_est_2$R)), is_true()) 
        expect_equal(R_est_2$R["R_org", ], 0.505, tolerance = 0.001)
        expect_equal(R_est_2$R["R_link", ], 0.482, tolerance = 0.001)
        
        # original scale
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["2.5%"]), 0.3648884, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["97.5%"]), 0.5191138, tolerance = 0.001)
        # link scale
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["2.5%"]), 0.3562329, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["97.5%"]), 0.496932, tolerance = 0.001)
        
})



# run with one random effect, no boot, permut
R_est_3 <- rptPoisson(Egg ~ Treatment + (1|Container), grname=c("Container"), data = BeetlesFemale,
                      nboot=0, npermut=5)

test_that("rpt estimation works for one random effect, no boot, permut, no parallelisation, logit link", {
        
        expect_that(is.numeric(unlist(R_est_3$R)), is_true()) 
        expect_equal(R_est_3$R["R_org", ], 0.505, tolerance = 0.001)
        expect_equal(R_est_3$R["R_link", ], 0.482, tolerance = 0.001)
        
        # original scale
        expect_equal(R_est_3$P$P_permut_org, 0.2, tolerance = 0.001)
        # link scale
        expect_equal(R_est_3$P$P_permut_link, 0.2, tolerance = 0.001)
        
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
R_est_1 <- rptPoisson(Egg ~ Treatment + (1|Container) + (1|Population), grname=c("Container", "Population"), data = BeetlesFemale,
        nboot=0, npermut=0)

test_that("rpt estimation works for two random effect, no boot, no permut, no parallelisation, logit link", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        # 1st random effect
        expect_equal(R_est_1$R["R_org", 1], 0.05716180, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", 1], 0.06018504, tolerance = 0.001)
        # 2nd random effect
        expect_equal(R_est_1$R["R_org", 2], 0.4881728, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", 2], 0.4668547, tolerance = 0.001)
})

test_that("LRTs works", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        expect_equal(R_est_1$P[1, "LRT_P"], 8.744869e-03, tolerance = 0.001)
        expect_equal(R_est_1$P[2, "LRT_P"], 1.187500e-16, tolerance = 0.001)
})


# run with one random effect, boot, no permut
R_est_2 <- rptPoisson(Egg ~ Treatment + (1|Container) + (1|Population), grname=c("Container", "Population"), data = BeetlesFemale,
        nboot=2, npermut=0)

test_that("rpt estimation works for two random effect, boot, no permut, no parallelisation, logit link", {
        
        # original scale
        # 1st random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org[1, "2.5%"]),  0.02544666, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org[1, "97.5%"]), 0.05610622, tolerance = 0.001)
        # 2nd random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org[2, "2.5%"]), 0.20203285 , tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org[2, "97.5%"]), 0.57533006, tolerance = 0.001)
        
        
        # link scale
        # 1st random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link[1, "2.5%"]),  0.02747585, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link[1, "97.5%"]),0.05960882, tolerance = 0.001)
        # 2nd random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link[2, "2.5%"]),  0.20581645, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link[2, "97.5%"]), 0.53812845, tolerance = 0.001)
        
})



# run with one random effect, no boot, permut
R_est_3 <- rptPoisson(Egg ~ Treatment + (1|Container) + (1|Population), grname=c("Container", "Population"), data = BeetlesFemale,
        nboot=0, npermut=5)

test_that("rpt estimation works for two random effect, no boot, permut, no parallelisation, logit link", {
        # original scale
        expect_equal(R_est_3$P$P_permut_org[1], 0.8, tolerance = 0.001)
        expect_equal(R_est_3$P$P_permut_org[2], 0.2, tolerance = 0.001)
        # link scale
        expect_equal(R_est_3$P$P_permut_link[1], 1, tolerance = 0.001)
        expect_equal(R_est_3$P$P_permut_link[2], 0.2, tolerance = 0.001)
        
})



# R_est_4 <- rptPoisson(Egg ~ Treatment + (1|Container), grname=c("Container"), data = BeetlesFemale,
#         nboot=5, npermut=5, parallel = TRUE)
