
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
        expect_equal(R_est_1$R["R_org", "Population"], 0.1853978, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", "Population"], 0.1879297, tolerance = 0.001)
})

test_that("LRT works", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        expect_equal(R_est_1$P$LRT_P, 8.656602e-15, tolerance = 0.001)
})



# run with one random effect, boot, no permut
R_est_2 <- rptBinary(Colour ~ (1|Population),  grname=c("Population"), data=BeetlesMale, nboot=2, npermut=0)

test_that("rpt estimation works for one random effect, boot, no permut, no parallelisation, logit link", {
        
        expect_that(is.numeric(unlist(R_est_2$R)), is_true()) 
        expect_equal(R_est_2$R["R_org", "Population"], 0.1853978, tolerance = 0.001)
        expect_equal(R_est_2$R["R_link", "Population"], 0.1879296, tolerance = 0.001)
        
        # original scale
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["Population", "2.5%"]), 0.1797193, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["Population", "97.5%"]), 0.2007206, tolerance = 0.001)
        # link scale
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["Population", "2.5%"]), 0.1818618, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["Population", "97.5%"]), 0.2028332, tolerance = 0.001)
        
})



# run with one random effect, no boot, permut
R_est_3 <- rptBinary(Colour ~ (1|Population),  grname=c("Population"), data=BeetlesMale, nboot=0, npermut=2)

test_that("rpt estimation works for one random effect, no boot, permut, no parallelisation, logit link", {
        
        expect_that(is.numeric(unlist(R_est_3$R)), is_true()) 
        expect_equal(R_est_3$R["R_org", "Population"], 0.1853978, tolerance = 0.001)
        expect_equal(R_est_3$R["R_link", "Population"], 0.1879296, tolerance = 0.001)
        
        # original scale
        expect_equal(R_est_3$P$P_permut_org, 0.5, tolerance = 0.001)
        # link scale
        expect_equal(R_est_3$P$P_permut_link, 0.5, tolerance = 0.001)
        
})





########## checks for two random effects

# run with one random effect, no boot, no permut
R_est_1 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Container", "Population", "Residual"), 
                     data=BeetlesMale, nboot=0, npermut=0)

test_that("rpt estimation works for two random effects plus residual, no boot, no permut, no parallelisation, logit link", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        # Container random effect
        expect_equal(R_est_1$R["R_org", "Container"], 0, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", "Container"], 0, tolerance = 0.001)
        # Population random effect
        expect_equal(R_est_1$R["R_org", "Population"], 0.1854004, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", "Population"], 0.1879323, tolerance = 0.001)
        # Residual variance
        #expect_equal(R_est_1$R["R_org", "Residual"], NA, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", "Residual"], 0.8120677, tolerance = 0.001)
})

test_that("LRTs works", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        # expect_equal(R_est_1$P[1, "LRT_P"], 1, tolerance = 0.001)
        expect_equal(R_est_1$P[1, "LRT_P"], 0.5, tolerance = 0.001)
        expect_equal(R_est_1$P[2, "LRT_P"], 5.80928e-09, tolerance = 0.001)
})

test_that("random effect components sum to up to one", {
        expect_equal(sum(R_est_1$R["R_link", ]), 1)
})


# run with two random effect, boot, no permut
R_est_2 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Residual", "Container", "Population"), data=BeetlesMale, nboot=2, npermut=0)

test_that("rpt estimation is independent of the order in grname", {
        # Container random effect
        expect_equal(R_est_1$R["R_org", "Container"], R_est_2$R["R_org", "Container"], tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", "Container"], R_est_2$R["R_link", "Container"], tolerance = 0.001)
        # Population random effect
        expect_equal(R_est_1$R["R_org", "Population"], R_est_2$R["R_org", "Population"], tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", "Population"], R_est_2$R["R_link", "Population"], tolerance = 0.001)
        # Residual random effect
        expect_equal(R_est_1$R["R_org", "Residual"], R_est_2$R["R_org", "Residual"], tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", "Residual"], R_est_2$R["R_link", "Residual"], tolerance = 0.001)
})

test_that("rpt estimation works for two random effect, boot, no permut, no parallelisation, logit link", {
        
        # original scale
        # Container random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["Container", "2.5%"]), 0.000460106 , tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["Container", "97.5%"]), 0.01794413, tolerance = 0.001)
        # Population random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["Population", "2.5%"]), 0.08705837, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["Population", "97.5%"]), 0.226984, tolerance = 0.001)
        # Residual random effect
        #expect_equal(as.numeric(R_est_2$CI_emp$CI_org["Residual", "2.5%"]), NA, tolerance = 0.001)
        #expect_equal(as.numeric(R_est_2$CI_emp$CI_org["Residual", "97.5%"]), NA, tolerance = 0.001)
        
        
        # link scale
        # Container random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["Container", "2.5%"]),  0.0004719506, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["Container", "97.5%"]), 0.01840607, tolerance = 0.001)
        # Population random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["Population", "2.5%"]), 0.08733554, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["Population", "97.5%"]), 0.232777, tolerance = 0.001)
        # Residual random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["Residual", "2.5%"]), 0.748817, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["Residual", "97.5%"]), 0.9121925, tolerance = 0.001)
        
})



# run with two random effects (no residual), no boot, permut
R_est_3 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Container", "Overdispersion", "Population"), data=BeetlesMale, nboot=0, npermut=2)

test_that("rpt estimation is independent of the presence of Overdispersion/Residual in grname", {
        # Container random effect
        expect_equal(R_est_1$R["R_org", "Container"], R_est_3$R["R_org", "Container"], tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", "Container"], R_est_3$R["R_link", "Container"], tolerance = 0.001)
        # Population random effect
        expect_equal(R_est_1$R["R_org", "Population"], R_est_3$R["R_org", "Population"], tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", "Population"], R_est_3$R["R_link", "Population"], tolerance = 0.001)
})

test_that("rpt estimation works for overdispersion, boot, no permut, no parallelisation, logit link", {
        # original scale
        # Overdispersion random effect
        expect_equal(R_est_3$R["R_org", "Overdispersion"], 0, tolerance = 0.001)
        expect_equal(R_est_3$R["R_link", "Overdispersion"], 0, tolerance = 0.001)

})


test_that("rpt estimation works for two random effect, no boot, permut, no parallelisation, logit link", {
        # original scale
        expect_equal(R_est_3$P["Container", "P_permut_org"], 1, tolerance = 0.001)
        expect_equal(R_est_3$P["Population", "P_permut_org"], 0.5, tolerance = 0.001)
        # link scale
        expect_equal(R_est_3$P["Container", "P_permut_link"], 1, tolerance = 0.001)
        expect_equal(R_est_3$P["Population", "P_permut_link"], 0.5, tolerance = 0.001)
        
})


# variance estimation with boot and permut and residual/overdispersion and Fixed effect
R_est_1 <- rptBinary(Colour ~ Treatment + (1|Container) + (1|Population),  grname=c("Fixed", "Container", "Population", "Residual", "Overdispersion"), 
        data=BeetlesMale, nboot=2, npermut=2, ratio = FALSE)

test_that("Variance estimation works for two random effects, boot, permut, Residual, Overdispersion and Fixed", {
        expect_equal(R_est_1$R["R_link", "Container"], 7.902152e-10, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", "Population"],  1.054672 , tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", "Overdispersion"], 0 , tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", "Residual"], 4.086918, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", "Fixed"], 0.2438065, tolerance = 0.001)
})


# test whether repeatabilities are equal for overdispersion independent of other grouping factors

R_est_1 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Container", "Population", "Residual", "Overdispersion", "Fixed"), 
        data=BeetlesMale, nboot=0, npermut=0, ratio = TRUE, adjusted=FALSE)
R_est_2 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Overdispersion", "Fixed"), 
        data=BeetlesMale, nboot=0, npermut=0, ratio = TRUE, adjusted=FALSE)


test_that("Repeatabilities are equal for grouping factors independent of residual and overdispersion specification", {
        expect_false(any(R_est_1$R[,"Overdispersion"] == R_est_2$R[,"Overdispersion"]) == FALSE)
        expect_false(any(R_est_1$R[,"Fixed"] == R_est_2$R[,"Fixed"]) == FALSE)
})

# check that grname and formula sequence dont play a role 

R_est_3 <- rptBinary(Colour ~ (1|Container) + (1|Population),  grname=c("Population", "Container"), 
        data=BeetlesMale, nboot=0, npermut=0, ratio = TRUE, adjusted=FALSE)

test_that("Repeatabilities are equal for different order in grname argument", {
     expect_equal(R_est_1$R[,"Container"], R_est_3$R[,"Container"])
     expect_equal(R_est_1$R[,"Population"], R_est_3$R[,"Population"]) 
})

test_that("LRTs are equal for different order in grname argument", {
        expect_equal(R_est_1$P["Container", "LRT_P"], R_est_3$P["Container", "LRT_P"])
        expect_equal(R_est_1$P["Population", "LRT_P"], R_est_3$P["Population", "LRT_P"])
})


R_est_4 <- rptBinary(Colour ~ (1|Population) + (1|Container),  grname=c("Container", "Population"), 
        data=BeetlesMale, nboot=0, npermut=0, ratio = TRUE)

test_that("Repeatabilities are equal for different order in formula argument", {
        expect_equal(R_est_3$R[,"Container"], R_est_4$R[,"Container"])
        expect_equal(R_est_3$R[,"Population"], R_est_4$R[,"Population"]) 
})

test_that("LRTs are equal for different order in formula argument", {
        expect_equal(R_est_3$P["Container", "LRT_P"], R_est_4$P["Container", "LRT_P"])
        expect_equal(R_est_3$P["Population", "LRT_P"], R_est_4$P["Population", "LRT_P"])
})
