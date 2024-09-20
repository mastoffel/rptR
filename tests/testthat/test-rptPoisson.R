
context("rptPoisson")
library(tibble)
# Set a seed for reproducibility of the randomization 
suppressWarnings(RNGversion("3.5.0"))
set.seed(23)

# load data
data(BeetlesFemale)

########## checks for one random effect

# run with one random effect, no boot, no permut
R_est_1 <- rptPoisson(Egg ~ Treatment + (1|Container), grname=c("Container"), data = BeetlesFemale,
                      nboot=0, npermut=0)

test_that("rpt estimation works for one random effect, no boot, no permut, no parallelisation, logit link", {
        expect_true(is.numeric(unlist(R_est_1$R))) 
        expect_equal(R_est_1$R["R_org", ], 0.5449128, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", ], 0.5550295, tolerance = 0.001)
})

test_that("LRT works", {
        expect_true(is.numeric(unlist(R_est_1$R))) 
        expect_equal(R_est_1$P$LRT_P, 1.98e-46, tolerance = 0.001)
})


# run with one random effect, boot, no permut
R_est_2 <- rptPoisson(Egg ~ Treatment + (1|Container), grname=c("Container"), data = BeetlesFemale,
                      nboot=2, npermut=0)

test_that("rpt estimation works for one random effect, boot, no permut, no parallelisation, logit link", {
        
        expect_true(is.numeric(unlist(R_est_2$R))) 
        expect_equal(R_est_2$R["R_org", ], 0.5449128, tolerance = 0.001)
        expect_equal(R_est_2$R["R_link", ], 0.5550295, tolerance = 0.001)
        
        # original scale
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["2.5%"]),  0.4124189, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["97.5%"]), 0.553597, tolerance = 0.001)
        # link scale
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["2.5%"]), 0.424651, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["97.5%"]), 0.5656179, tolerance = 0.001)
        
})




# run with one random effect, no boot, permut
R_est_3 <- rptPoisson(Egg ~ Treatment + (1|Container), grname=c("Container"), data = BeetlesFemale,
                      nboot=0, npermut=5)

test_that("rpt estimation works for one random effect, no boot, permut, no parallelisation, logit link", {
        
        expect_true(is.numeric(unlist(R_est_3$R))) 
        expect_equal(R_est_3$R["R_org", ], 0.5449128, tolerance = 0.001)
        expect_equal(R_est_3$R["R_link", ], 0.5550295, tolerance = 0.001)
        
        # original scale
        expect_equal(R_est_3$P$P_permut_org, 0.2, tolerance = 0.001)
        # link scale
        expect_equal(R_est_3$P$P_permut_link, 0.2, tolerance = 0.001)
        
})






########## checks for two random effects

# run with two random effect, no boot, no permut
R_est_1 <- rptPoisson(Egg ~ Treatment + (1|Container) + (1|Population), grname=c("Container", "Population"), data = BeetlesFemale,
        nboot=0, npermut=0)

test_that("rpt estimation works for two random effect, no boot, no permut, no parallelisation, logit link", {
        expect_true(is.numeric(unlist(R_est_1$R))) 
        # 1st random effect
        expect_equal(R_est_1$R["R_org", 1], 0.03190105, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", 1],0.03799297, tolerance = 0.001)
        # 2nd random effect
        expect_equal(R_est_1$R["R_org", 2], 0.5041709, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", 2], 0.5195091, tolerance = 0.001)
})

test_that("LRTs works", {
        expect_true(is.numeric(unlist(R_est_1$R))) 
        expect_equal(R_est_1$P[1, "LRT_P"], 8.744869e-03, tolerance = 0.001)
        expect_equal(R_est_1$P[2, "LRT_P"], 1.187500e-16, tolerance = 0.001)
})


# check that random effect components plus overdispersion variance sum up to one
R_est_1 <- rptPoisson(Egg ~ Treatment + (1|Container) + (1|Population), 
        grname=c("Container", "Population", "Residual"), data = BeetlesFemale,
        nboot=0, npermut=0)

test_that("random effect components sum to up to one", {
        expect_equal(sum(R_est_1$R["R_link", ]), 1)
})


# run with two random effect, boot, no permut
R_est_2 <- rptPoisson(Egg ~ Treatment + (1|Container) + (1|Population), grname=c("Container", "Population"), data = BeetlesFemale,
        nboot=2, npermut=0)

test_that("rpt estimation works for two random effect, boot, no permut, no parallelisation, logit link", {
        
        # original scale
        # 1st random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org[1, "2.5%"]),  0.02431758, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org[1, "97.5%"]), 0.0256355, tolerance = 0.001)
        # 2nd random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org[2, "2.5%"]), 0.2346106 , tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org[2, "97.5%"]),0.5839979, tolerance = 0.001)
        
        
        # link scale
        # 1st random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link[1, "2.5%"]), 0.02642277, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link[1, "97.5%"]), 0.03203479, tolerance = 0.001)
        # 2nd random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link[2, "2.5%"]),  0.2449164, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link[2, "97.5%"]), 0.5944401, tolerance = 0.001)
        
})



# run with two random effect, no boot, permut
R_est_3 <- rptPoisson(Egg ~ Treatment + (1|Container) + (1|Population), grname=c("Container", "Population"), data = BeetlesFemale,
        nboot=0, npermut=5)

test_that("rpt estimation works for two random effect, no boot, permut, no parallelisation, logit link", {
        # original scale
        expect_equal(R_est_3$P$P_permut_org[1], 0.2, tolerance = 0.001)
        expect_equal(R_est_3$P$P_permut_org[2], 0.2, tolerance = 0.001)
        # link scale
        expect_equal(R_est_3$P$P_permut_link[1], 0.2, tolerance = 0.001)
        expect_equal(R_est_3$P$P_permut_link[2], 0.2, tolerance = 0.001)
        
})




# test that Ratio works

# run with one random effect, no boot, no permut
R_est_1 <- rptPoisson(Egg ~ Treatment + (1|Container), grname=c("Container"), data = BeetlesFemale,
        nboot=5, npermut=5, ratio = FALSE)

test_that("rpt estimation works for one random effect, no boot, no permut, no parallelisation, logit link", {
        expect_true(is.numeric(unlist(R_est_1$R))) 
        expect_true(is.na(R_est_1$R["R_org", ]))
        expect_equal(R_est_1$R["R_link", ], 0.3229846, tolerance = 0.001)
})

test_that("LRT works", {
        expect_true(is.numeric(unlist(R_est_1$R))) 
        expect_equal(R_est_1$P$LRT_P, 1.98e-46, tolerance = 0.001)
})

# run with two random effect, no boot, no permut
R_est_2 <- rptPoisson(Egg ~ Treatment + (1|Container) + (1|Population), grname=c("Container", "Population"), data = BeetlesFemale,
        nboot=0, npermut=0, ratio = FALSE)

test_that("rpt estimation works for two random effect, no boot, no permut, no parallelisation, logit link", {
        expect_true(is.numeric(unlist(R_est_1$R))) 
        # 1st random effect
        expect_true(is.na(R_est_2$R["R_org", 1]))
        expect_equal(R_est_2$R["R_link", 1],  0.02224455, tolerance = 0.001)
        # 2nd random effect
        expect_true(is.na(R_est_2$R["R_org", 2]))
        expect_equal(R_est_2$R["R_link", 2], 0.3041681, tolerance = 0.001)
})


# test that Residual and Overdispersion works

R_est_1 <-  rptPoisson(formula = Egg ~ Treatment + (1|Container) + (1|Habitat) ,
grname=c("Container", "Habitat", "Residual", "Overdispersion"), data = BeetlesFemale,
nboot=5, npermut=5, ratio = FALSE)

test_that("rpt estimation works for two random effects, estimation of Variance and Residual / Overdispersion", {
        expect_true(is.na(R_est_1$R["R_org", "Overdispersion"]))
        expect_true(is.na(R_est_1$R["R_org", "Residual"]))
        expect_equal(R_est_1$R["R_link", "Overdispersion"], 0.1004924, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", "Residual"],0.2570103, tolerance = 0.001)
})


R_est_1 <-  rptPoisson(formula = Egg ~ Treatment + (1|Container) + (1|Habitat) ,
        grname=c("Container", "Habitat", "Residual", "Overdispersion"), data = BeetlesFemale,
        nboot=0, npermut=0, ratio = TRUE)
R_est_2 <-  rptPoisson(formula = Egg ~ Treatment + (1|Container) + (1|Habitat) ,
        grname=c("Container", "Habitat"), data = BeetlesFemale,
        nboot=0, npermut=0, ratio = TRUE)


# test whether repeatabilities are equal for grouping factors independent of residual and overdispersion specification

test_that("repeatabilities are equal for grouping factors independent of residual and overdispersion specification", {
        expect_equal(R_est_1$R$Habitat, R_est_2$R$Habitat)
        expect_equal(R_est_1$R$Container, R_est_2$R$Container)
})


# check that grname sequence doesnt play a role
R_est_3 <-  rptPoisson(formula = Egg ~ Treatment + (1|Container) + (1|Habitat) ,
        grname=c("Habitat", "Container"), data = BeetlesFemale,
        nboot=0, npermut=0, ratio = TRUE)

test_that("Repeatabilities are equal for different order in grname argument", {
        expect_false(any(R_est_2$R$Container == R_est_3$R$Container) == FALSE)
        expect_false(any(R_est_2$R$Habitat == R_est_3$R$Habitat) == FALSE)
})

test_that("LRTs are equal for different different sequence in grname argument", {
        expect_equal(R_est_2$P["Container", "LRT_P"], R_est_3$P["Container", "LRT_P"])
        expect_equal(R_est_2$P["Habitat", "LRT_P"], R_est_3$P["Habitat", "LRT_P"])
})


R_est_4 <-  rptPoisson(formula = Egg ~ Treatment + (1|Habitat) + (1|Container),
        grname=c("Container", "Habitat"), data = BeetlesFemale,
        nboot=0, npermut=0, ratio = TRUE)

test_that("Repeatabilities are equal for different order in formula argument", {
        expect_false(any(R_est_2$R$Container == R_est_4$R$Container) == FALSE)
        expect_false(any(R_est_2$R$Habitat == R_est_4$R$Habitat) == FALSE)
})

test_that("LRTs are equal for different for different order in formula argument", {
        expect_equal(R_est_2$P["Container", "LRT_P"], R_est_4$P["Container", "LRT_P"])
        expect_equal(R_est_2$P["Habitat", "LRT_P"], R_est_4$P["Habitat", "LRT_P"])
})

# random slopes
R_est_5 <- rptPoisson(formula = Egg ~ Treatment + (1 + Treatment|Container),
        grname=c("Container"), data = BeetlesFemale,
        nboot=0, npermut=0, ratio = TRUE)

test_that("Random Slope repeatability point estimate is correct", {
        expect_equal(R_est_5$R["R_org", "Container"], 0.5539859, tolerance = 0.001)
})

test_that("Random slopes fitted without correlation (i.e. with grname occuring in multiple random effect terms) throw an error", {
        expect_error(rptPoisson(formula = Egg ~ (1|Treatment) + (0 + Treatment|Container),
                grname=c("Container"), data = BeetlesFemale,
                nboot=0, npermut=0, ratio = TRUE))
})

# run with one random effect, no boot, no permut
R_est_6 <- rptPoisson(Egg ~ Treatment + (1|Container), grname=c("Container"), data = tibble(BeetlesFemale),
                      nboot=0, npermut=0)

test_that("rptPoisson runs with a tibble as data", {
        expect_true(is.numeric(unlist(R_est_1$R))) 
})