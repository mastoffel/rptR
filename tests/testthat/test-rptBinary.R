
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
        expect_equal(R_est_1$R["R_org", ], 1.080313, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", ], 0.1879297, tolerance = 0.001)
})

test_that("LRT works", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        expect_equal(R_est_1$P$LRT_P, 8.656602e-15, tolerance = 0.001)
})



# run with one random effect, boot, no permut
R_est_2 <- rptBinary(Colour ~ (1|Population),  grname=c("Population"), data=BeetlesMale, nboot=2, npermut=0)

test_that("rpt estimation works for one random effect, boot, no permut, no parallelisation, logit link", {
        
        expect_that(is.numeric(unlist(R_est_2$R)), is_true()) 
        expect_equal(R_est_2$R["R_org", ], 1.080313, tolerance = 0.001)
        expect_equal(R_est_2$R["R_link", ], 0.1879297, tolerance = 0.001)
        
        # original scale
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["2.5%"]),1.072162, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["97.5%"]), 1.085152, tolerance = 0.001)
        # link scale
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["2.5%"]), 0.1818619, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["97.5%"]), 0.202833, tolerance = 0.001)
        
})



# run with one random effect, no boot, permut
R_est_3 <- rptBinary(Colour ~ (1|Population),  grname=c("Population"), data=BeetlesMale, nboot=0, npermut=5)

test_that("rpt estimation works for one random effect, no boot, permut, no parallelisation, logit link", {
        
        expect_that(is.numeric(unlist(R_est_3$R)), is_true()) 
        expect_equal(R_est_3$R["R_org", ], 1.080313, tolerance = 0.001)
        expect_equal(R_est_3$R["R_link", ], 0.1879297, tolerance = 0.001)
        
        # original scale
        expect_equal(R_est_3$P$P_permut_org, 0.2, tolerance = 0.001)
        # link scale
        expect_equal(R_est_3$P$P_permut_link, 0.2, tolerance = 0.001)
        
})





########## checks for two random effects

# run with one random effect, no boot, no permut
R_est_1 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Container", "Population"), 
                     data=BeetlesMale, nboot=0, npermut=0)

test_that("rpt estimation works for two random effect, no boot, no permut, no parallelisation, logit link", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        # 1st random effect
        expect_equal(R_est_1$R["R_org", 1], 0, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", 1], 0, tolerance = 0.001)
        # 2nd random effect
        expect_equal(R_est_1$R["R_org", 2], 1.080313, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", 2], 0.1879323, tolerance = 0.001)
})

test_that("LRTs works", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        expect_equal(R_est_1$P[1, "LRT_P"], 1, tolerance = 0.001)
        expect_equal(R_est_1$P[2, "LRT_P"], 5.81e-09, tolerance = 0.001)
})


# Variance components sum up to one
R_est_1 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Container", "Population", "Residual"), 
        data=BeetlesMale, nboot=0, npermut=0)

test_that("random effect components sum to up to one", {
        expect_equal(sum(R_est_1$R["R_link", ]), 1)
})

# run with one random effect, boot, no permut
R_est_2 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Container", "Population"), data=BeetlesMale, nboot=2, npermut=0)

test_that("rpt estimation works for two random effect, boot, no permut, no parallelisation, logit link", {
        
        # original scale
        # 1st random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org[1, "2.5%"]), 9.241147e-12 , tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org[1, "97.5%"]), 3.600875e-10, tolerance = 0.001)
        # 2nd random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org[2, "2.5%"]), 1.058988, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org[2, "97.5%"]), 1.103392, tolerance = 0.001)
        
        
        # link scale
        # 1st random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link[1, "2.5%"]),  1.520977e-11, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link[1, "97.5%"]), 5.927945e-10, tolerance = 0.001)
        # 2nd random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link[2, "2.5%"]),0.1392826, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link[2, "97.5%"]),0.1571744, tolerance = 0.001)
        
})



# run with two random effects, no boot, permut
R_est_3 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Container", "Population"), data=BeetlesMale, nboot=0, npermut=2)

test_that("rpt estimation works for two random effect, no boot, permut, no parallelisation, logit link", {
        # original scale
        expect_equal(R_est_3$P$P_permut_org[1], 1, tolerance = 0.001)
        expect_equal(R_est_3$P$P_permut_org[2], 0.5, tolerance = 0.001)
        # link scale
        expect_equal(R_est_3$P$P_permut_link[1], 1, tolerance = 0.001)
        expect_equal(R_est_3$P$P_permut_link[2], 0.5, tolerance = 0.001)
        
})


# variance estimation with boot and permut and residual/overdispersion
R_est_1 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Container", "Population", "Residual", "Overdispersion"), 
        data=BeetlesMale, nboot=2, npermut=2, ratio = FALSE)

test_that("Variance estimation works for two random effects, boot, permut, Residual and Overdispersion", {
        expect_equal(R_est_1$R$Container[2], 0, tolerance = 0.001)
        expect_equal(R_est_1$R$Population[2], 0.9458128 , tolerance = 0.001)
        expect_equal(R_est_1$R$Overdispersion[2], 0 , tolerance = 0.001)
        expect_equal(R_est_1$R$Residual[2], 4.086918, tolerance = 0.001)
})


# test whether repeatabilities are equal for grouping factors independent of residual and overdispersion specification

R_est_1 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Container", "Population", "Residual", "Overdispersion"), 
        data=BeetlesMale, nboot=0, npermut=0, ratio = TRUE)
R_est_2 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Container", "Population"), 
        data=BeetlesMale, nboot=0, npermut=0, ratio = TRUE)


test_that("Repeatabilities are equal for grouping factors independent of residual and overdispersion specification", {
        expect_false(any(R_est_1$R$Container == R_est_2$R$Container) == FALSE)
        expect_false(any(R_est_1$R$Population == R_est_2$R$Population) == FALSE)
})

# check that grname and formula sequence dont play a role 

R_est_3 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Population", "Container"), 
        data=BeetlesMale, nboot=0, npermut=0, ratio = TRUE)

test_that("Repeatabilities are equal for different order in grname argument", {
     expect_equal(R_est_2$R$Container, R_est_3$R$Container)
     expect_equal(R_est_2$R$Population, R_est_3$R$Population) 
})

test_that("LRTs are equal for different order in grname argument", {
        expect_equal(R_est_2$P["Container", "LRT_P"], R_est_3$P["Container", "LRT_P"])
        expect_equal(R_est_2$P["Population", "LRT_P"], R_est_3$P["Population", "LRT_P"])
})


R_est_4 <- rptBinary(Colour ~ (1|Population) + (1|Container),  grname=c("Container", "Population"), 
        data=BeetlesMale, nboot=0, npermut=0, ratio = TRUE)

test_that("Repeatabilities are equal for different order in formula argument", {
        expect_equal(R_est_2$R$Container, R_est_4$R$Container)
        expect_equal(R_est_2$R$Population, R_est_4$R$Population) 
})

test_that("LRTs are equal for different order in formula argument", {
        expect_equal(R_est_2$P["Container", "LRT_P"], R_est_4$P["Container", "LRT_P"])
        expect_equal(R_est_2$P["Population", "LRT_P"], R_est_4$P["Population", "LRT_P"])
})
