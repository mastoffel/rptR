
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


