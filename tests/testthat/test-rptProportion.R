
context("rptProportion")

# Set a seed for reproducibility of the randomization 
set.seed(23)

# load data
data(BeetlesMale)

BeetlesMale$Dark <- BeetlesMale$Colour
BeetlesMale$Reddish <- (BeetlesMale$Colour-1)*-1
md <- aggregate(cbind(Dark, Reddish) ~ Population + Container, data=BeetlesMale, FUN=sum)

########## checks for one random effect

# run with one random effect, no boot, no permut
R_est_1 <- rptProportion(cbind(Dark, Reddish) ~ (1|Population), grname=c("Population"), data=md,
                         nboot=0, npermut=0)

test_that("rpt estimation works for one random effect, no boot, no permut, no parallelisation, logit link", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        expect_equal(R_est_1$R["R_org", ],0.186, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", ], 0.223, tolerance = 0.001)
})

test_that("LRT works", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        expect_equal(R_est_1$P$LRT_P, 5.81e-09, tolerance = 0.001)
})


# run with one random effect, boot, no permut
R_est_2 <- rptProportion(cbind(Dark, Reddish) ~ (1|Population), grname=c("Population"), data=md,
        nboot=2, npermut=0)

test_that("rpt estimation works for one random effect, boot, no permut, no parallelisation, logit link", {
        
        # original scale
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["2.5%"]), 0.1573669, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org["97.5%"]), 0.1911474, tolerance = 0.001)
        # link scale
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["2.5%"]), 0.1998974, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link["97.5%"]), 0.2258376, tolerance = 0.001)
        
})



# run with one random effect, no boot, permut
R_est_3 <- rptProportion(cbind(Dark, Reddish) ~ (1|Population), grname=c("Population"), data=md,
        nboot=0, npermut=5)

test_that("rpt estimation works for one random effect, no boot, permut, no parallelisation, logit link", {
        
        # original scale
        expect_equal(R_est_3$P$P_permut_org, 0.2, tolerance = 0.001)
        # link scale
        expect_equal(R_est_3$P$P_permut_link, 0.2, tolerance = 0.001)
        
})


# run with one random effect, no boot, no permut
# R_est_4 <- rptProportion(cbind(Dark, Reddish) ~ (1|Population), grname=c("Population"), data=md,
#         nboot=5, npermut=5, parallel = TRUE)



########## checks for two random effects

# run with one random effect, no boot, no permut
R_est_1 <- suppressWarnings(rptProportion(cbind(Dark, Reddish) ~ (1|Container) + (1|Population), grname=c("Container", "Population"), data = md,
        nboot=0, npermut=0))

test_that("rpt estimation works for two random effect, no boot, no permut, no parallelisation, logit link", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        # 1st random effect
        expect_equal(R_est_1$R["R_org", 1], 3.720300e-11, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", 1], 4.686765e-11 , tolerance = 0.001)
        # 2nd random effect
        expect_equal(R_est_1$R["R_org", 2], 0.1858069, tolerance = 0.001)
        expect_equal(R_est_1$R["R_link", 2], 0.2232977, tolerance = 0.001)
})

test_that("LRTs works", {
        expect_that(is.numeric(unlist(R_est_1$R)), is_true()) 
        expect_equal(R_est_1$P[1, "LRT_P"], 1, tolerance = 0.001)
        expect_equal(R_est_1$P[2, "LRT_P"], 5.81e-09, tolerance = 0.001)
})

# variance components sum up to 1
R_est_1 <- suppressWarnings(rptProportion(cbind(Dark, Reddish) ~ (1|Container) + (1|Population), 
        grname=c("Container", "Population", "Residual"), data = md,
        nboot=0, npermut=0))

test_that("random effect components sum to up to one", {
        expect_equal(sum(R_est_1$R["R_link", ]), 1)
})


# run with one random effect, boot, no permut
R_est_2 <- suppressWarnings(rptProportion(cbind(Dark, Reddish) ~ (1|Container) + (1|Population), grname=c("Container", "Population"), data = md,
        nboot=2, npermut=0))

test_that("rpt estimation works for two random effect, boot, no permut, no parallelisation, logit link", {
        
        # original scale
        # 1st random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org[1, "2.5%"]),  0, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org[1, "97.5%"]), 0, tolerance = 0.001)
        # 2nd random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org[2, "2.5%"]),  1.280478e-01 , tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_org[2, "97.5%"]), 1.814857e-01, tolerance = 0.001)
        
        
        # link scale
        # 1st random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link[1, "2.5%"]),  0, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link[1, "97.5%"]),0, tolerance = 0.001)
        # 2nd random effect
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link[2, "2.5%"]), 1.623094e-01, tolerance = 0.001)
        expect_equal(as.numeric(R_est_2$CI_emp$CI_link[2, "97.5%"]),  2.126430e-01, tolerance = 0.001)
        
})



# run with one random effect, no boot, permut
R_est_3 <- suppressWarnings(rptProportion(cbind(Dark, Reddish) ~ (1|Container) + (1|Population), grname=c("Container", "Population"), data = md,
        nboot=0, npermut=5))

test_that("rpt estimation works for two random effect, no boot, permut, no parallelisation, logit link", {
        # original scale
        expect_equal(R_est_3$P$P_permut_org[1], 1, tolerance = 0.001)
        expect_equal(R_est_3$P$P_permut_org[2], 0.2, tolerance = 0.001)
        # link scale
        expect_equal(R_est_3$P$P_permut_link[1], 1, tolerance = 0.001)
        expect_equal(R_est_3$P$P_permut_link[2], 0.2, tolerance = 0.001)
        
})

R_est_4 <- suppressWarnings(rptProportion(cbind(Dark, Reddish) ~ (1|Container) + (1|Population), 
        grname=c("Container", "Population", "Overdispersion", "Residual"), data = md,
        nboot=0, npermut=0))

test_that("repeatabilities are equal for grouping factors independent of residual and overdispersion specification", {
        expect_false(any(R_est_3$R$Container == R_est_4$R$Container) == FALSE)
        expect_false(any(R_est_3$R$Population == R_est_4$R$Population) == FALSE)
})



# check that estimates are independent from the order in grname or formula

R_est_5 <- suppressWarnings(rptProportion(cbind(Dark, Reddish) ~ (1|Container) + (1|Population), 
        grname=c("Population", "Container"), data = md,
        nboot=0, npermut=0))

test_that("repeatabilities are equal for different sequence in grname argument", {
        expect_false(any(R_est_3$R$Container == R_est_5$R$Container) == FALSE)
        expect_false(any(R_est_3$R$Population == R_est_5$R$Population) == FALSE)
})

test_that("LRTs are equal for different different sequence in grname argument", {
        expect_equal(R_est_3$P["Container", "LRT_P"], R_est_5$P["Container", "LRT_P"])
        expect_equal(R_est_3$P["Population", "LRT_P"], R_est_5$P["Population", "LRT_P"])
})




R_est_6 <- suppressWarnings(rptProportion(cbind(Dark, Reddish) ~ (1|Population) + (1|Container), 
        grname=c("Container", "Population"), data = md,
        nboot=0, npermut=0))

test_that("repeatabilities are equal for different sequence in formula argument", {
        expect_false(any(R_est_3$R$Container == R_est_6$R$Container) == FALSE)
        expect_false(any(R_est_3$R$Population == R_est_6$R$Population) == FALSE)
})

test_that("LRTs are equal for different order in formula argument", {
        expect_equal(R_est_3$P["Container", "LRT_P"], R_est_6$P["Container", "LRT_P"])
        expect_equal(R_est_3$P["Population", "LRT_P"], R_est_6$P["Population", "LRT_P"])
})
