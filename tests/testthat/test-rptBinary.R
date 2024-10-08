
context("rptBinary")

# load data
data(BeetlesMale)

# Set a seed for reproducibility of the randomization 
suppressWarnings(RNGversion("3.5.0"))
set.seed(23)

########## checks for one random effect

# run with one random effect, no boot, no permut
R_est_1 <- rptBinary(Colour ~ (1|Population),  grname=c("Population"), data=BeetlesMale, nboot=0, npermut=0)

test_that("rpt estimation works for one random effect, no boot, no permut, no parallelisation, logit link", {
        expect_true(is.numeric(unlist(R_est_1$R))) 
        expect_equal(R_est_1$R["R_org", "Population"], 0.1853978, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", "Population"], 0.1879297, tolerance = 0.001)
})

test_that("LRT works", {
        expect_true(is.numeric(unlist(R_est_1$R))) 
        expect_equal(R_est_1$P$LRT_P, 8.656602e-15, tolerance = 0.001)
})



# run with one random effect, boot, no permut
R_est_2 <- rptBinary(Colour ~ (1|Population),  grname=c("Population"), data=BeetlesMale, nboot=2, npermut=0)

test_that("rpt estimation works for one random effect, boot, no permut, no parallelisation, logit link", {
        
        expect_true(is.numeric(unlist(R_est_2$R))) 
        expect_equal(R_est_2$R["R_org", "Population"], 0.1853978, tolerance = 0.001)
        expect_equal(R_est_2$R["R_link", "Population"], 0.1879296, tolerance = 0.001)
        
        # original scale
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["Population", "2.5%"]), 0.1448806, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["Population", "97.5%"]), 0.1660169, tolerance = 0.001)
        # link scale
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["Population", "2.5%"]), 0.1458783, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["Population", "97.5%"]), 0.1654149, tolerance = 0.001)
        
})



# run with one random effect, no boot, permut
R_est_3 <- rptBinary(Colour ~ (1|Population),  grname=c("Population"), data=BeetlesMale, nboot=0, npermut=2)

test_that("rpt estimation works for one random effect, no boot, permut, no parallelisation, logit link", {
        
        expect_true(is.numeric(unlist(R_est_3$R))) 
        expect_equal(R_est_3$R["R_org", "Population"], 0.1853978, tolerance = 0.001)
        expect_equal(R_est_3$R["R_link", "Population"], 0.1879296, tolerance = 0.001)
        
        # original scale
        expect_equal(R_est_3$P$P_permut_org, 0.5, tolerance = 0.001)
        # link scale
        expect_equal(R_est_3$P$P_permut_link, 0.5, tolerance = 0.001)
        
})





########## checks for two random effects

# run with one random effect, no boot, no permut
R_est_4 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Container", "Population", "Residual"), 
                     data=BeetlesMale, nboot=0, npermut=0)

test_that("rpt estimation works for two random effects plus residual, no boot, no permut, no parallelisation, logit link", {
        expect_true(is.numeric(unlist(R_est_4$R))) 
        # Container random effect
        expect_equal(R_est_4$R["R_org", "Container"], 0, tolerance = 0.001)
        expect_equal(R_est_4$R["R_link", "Container"], 0, tolerance = 0.001)
        # Population random effect
        expect_equal(R_est_4$R["R_org", "Population"], 0.1854004, tolerance = 0.001)
        expect_equal(R_est_4$R["R_link", "Population"], 0.1879323, tolerance = 0.001)
        # Residual variance
        #expect_equal(R_est_1$R["R_org", "Residual"], NA, tolerance = 0.001)
        expect_equal(R_est_4$R["R_link", "Residual"], 0.8120677, tolerance = 0.001)
})

test_that("LRTs works", {
        expect_true(is.numeric(unlist(R_est_4$R))) 
        # expect_equal(R_est_1$P[1, "LRT_P"], 1, tolerance = 0.001)
        expect_equal(R_est_4$P[1, "LRT_P"], 1, tolerance = 0.001) # previously 0.5
        expect_equal(R_est_4$P[2, "LRT_P"], 5.80928e-09, tolerance = 0.001)
})

test_that("random effect components sum to up to one", {
        expect_equal(sum(R_est_4$R["R_link", ]), 1)
})


# run with two random effect, boot, no permut
R_est_5 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Residual", "Container", "Population"), data=BeetlesMale, nboot=2, npermut=0)

test_that("rpt estimation is independent of the order in grname", {
        # Container random effect
        expect_equal(R_est_4$R["R_org", "Container"], R_est_5$R["R_org", "Container"], tolerance = 0.001)
        expect_equal(R_est_4$R["R_link", "Container"], R_est_5$R["R_link", "Container"], tolerance = 0.001)
        # Population random effect
        expect_equal(R_est_4$R["R_org", "Population"], R_est_5$R["R_org", "Population"], tolerance = 0.001)
        expect_equal(R_est_4$R["R_link", "Population"], R_est_5$R["R_link", "Population"], tolerance = 0.001)
        # Residual random effect
        expect_equal(R_est_4$R["R_org", "Residual"], R_est_5$R["R_org", "Residual"], tolerance = 0.001)
        expect_equal(R_est_4$R["R_link", "Residual"], R_est_5$R["R_link", "Residual"], tolerance = 0.001)
})

test_that("rpt estimation works for two random effect, boot, no permut, no parallelisation, logit link", {
        
        # original scale
        # Container random effect
        expect_equal(as.numeric(R_est_5$CI_emp$CI_org["Container", "2.5%"]), 0.000214565 , tolerance = 0.001)
        expect_equal(as.numeric(R_est_5$CI_emp$CI_org["Container", "97.5%"]),0.008368035, tolerance = 0.001)
        # Population random effect
        expect_equal(as.numeric(R_est_5$CI_emp$CI_org["Population", "2.5%"]),  0.091740174, tolerance = 0.001)
        expect_equal(as.numeric(R_est_5$CI_emp$CI_org["Population", "97.5%"]), 0.206301801, tolerance = 0.001)
        # Residual random effect
        #expect_equal(as.numeric(R_est_5$CI_emp$CI_org["Residual", "2.5%"]), NA, tolerance = 0.001)
        #expect_equal(as.numeric(R_est_5$CI_emp$CI_org["Residual", "97.5%"]), NA, tolerance = 0.001)
        
        
        # link scale
        # Container random effect
        expect_equal(as.numeric(R_est_5$CI_emp$CI_link["Container", "2.5%"]),  0.000215459, tolerance = 0.001)
        expect_equal(as.numeric(R_est_5$CI_emp$CI_link["Container", "97.5%"]), 0.0084029, tolerance = 0.001)
        # Population random effect
        expect_equal(as.numeric(R_est_5$CI_emp$CI_link["Population", "2.5%"]),  0.092182856, tolerance = 0.001)
        expect_equal(as.numeric(R_est_5$CI_emp$CI_link["Population", "97.5%"]),  0.2095188, tolerance = 0.001)
        # Residual random effect
        expect_equal(as.numeric(R_est_5$CI_emp$CI_link["Residual", "2.5%"]), 0.790265699, tolerance = 0.001)
        expect_equal(as.numeric(R_est_5$CI_emp$CI_link["Residual", "97.5%"]), 0.8994142, tolerance = 0.001)
        
})



# run with two random effects (no residual), no boot, permut
R_est_6 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Container", "Population"), data=BeetlesMale, nboot=0, npermut=2)

test_that("rpt estimation is independent of the presence of Residual in grname", {
        # Container random effect
        expect_equal(R_est_4$R["R_org", "Container"], R_est_6$R["R_org", "Container"], tolerance = 0.001)
        expect_equal(R_est_4$R["R_link", "Container"], R_est_6$R["R_link", "Container"], tolerance = 0.001)
        # Population random effect
        expect_equal(R_est_4$R["R_org", "Population"], R_est_6$R["R_org", "Population"], tolerance = 0.001)
        expect_equal(R_est_4$R["R_link", "Population"], R_est_6$R["R_link", "Population"], tolerance = 0.001)
})


test_that("rpt estimation works for two random effect, no boot, permut, no parallelisation, logit link", {
        # original scale
        expect_equal(R_est_6$P["Container", "P_permut_org"], 1, tolerance = 0.001)
        expect_equal(R_est_6$P["Population", "P_permut_org"], 0.5, tolerance = 0.001)
        # link scale
        expect_equal(R_est_6$P["Container", "P_permut_link"], 1, tolerance = 0.001)
        expect_equal(R_est_6$P["Population", "P_permut_link"], 0.5, tolerance = 0.001)
        
})


# variance estimation with boot and permut and residual/overdispersion and Fixed effect
R_est_7 <- rptBinary(Colour ~ Treatment + (1|Container) + (1|Population),  grname=c("Fixed", "Container", "Population", "Residual"), 
        data=BeetlesMale, nboot=2, npermut=2, ratio = FALSE)

test_that("Variance estimation works for two random effects, boot, permut, Residual and Fixed", {
        expect_equal(R_est_7$R["R_link", "Container"], 7.902152e-10, tolerance = 0.001)
        expect_equal(R_est_7$R["R_link", "Population"],  1.054672 , tolerance = 0.001)
        expect_equal(R_est_7$R["R_link", "Residual"], 4.086918, tolerance = 0.001)
        expect_equal(R_est_7$R["R_link", "Fixed"], 0.2438065, tolerance = 0.001)
})


# test whether repeatabilities are equal independent of other grouping factors

R_est_8 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Container", "Population", "Residual", "Fixed"), 
        data=BeetlesMale, nboot=0, npermut=0, ratio = TRUE, adjusted=FALSE)
R_est_9 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Fixed", "Container"), 
        data=BeetlesMale, nboot=0, npermut=0, ratio = TRUE, adjusted=FALSE)


test_that("Repeatabilities are equal for grouping factors independent of residual specification", {
        expect_false(any(R_est_8$R[,"Fixed"] == R_est_9$R[,"Fixed"]) == FALSE)
})

# check that grname and formula sequence dont play a role 

R_est_10 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Population", "Container"), 
        data=BeetlesMale, nboot=0, npermut=0, ratio = TRUE, adjusted=FALSE)

test_that("Repeatabilities are equal for different order in grname argument", {
     expect_equal(R_est_10$R[,"Container"], R_est_6$R[,"Container"])
     expect_equal(R_est_10$R[,"Population"], R_est_6$R[,"Population"]) 
})

test_that("LRTs are equal for different order in grname argument", {
        expect_equal(R_est_10$P["Container", "LRT_P"], R_est_6$P["Container", "LRT_P"])
        expect_equal(R_est_10$P["Population", "LRT_P"], R_est_6$P["Population", "LRT_P"])
})

R_est_11 <- rptBinary(Colour ~ (1|Population) + (1|Container),  grname=c("Container", "Population"), 
        data=BeetlesMale, nboot=0, npermut=0, ratio = TRUE)

test_that("Repeatabilities are equal for different order in formula argument", {
        expect_equal(R_est_10$R[,"Container"], R_est_11$R[,"Container"])
        expect_equal(R_est_10$R[,"Population"], R_est_11$R[,"Population"]) 
})

test_that("LRTs are equal for different order in formula argument", {
        expect_equal(R_est_10$P["Container", "LRT_P"], R_est_11$P["Container", "LRT_P"])
        expect_equal(R_est_10$P["Population", "LRT_P"], R_est_11$P["Population", "LRT_P"])
})

