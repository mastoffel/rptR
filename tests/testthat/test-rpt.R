context("rpt")

# Set a seed for reproducibility of the randomization 
suppressWarnings(RNGversion("3.5.0"))
set.seed(23)

# load data
data(BeetlesBody)
data(BeetlesMale)
data(BeetlesFemale)

# prepare proportion data
BeetlesMale$Dark <- BeetlesMale$Colour
BeetlesMale$Reddish <- (BeetlesMale$Colour-1)*-1
md <- aggregate(cbind(Dark, Reddish) ~ Population + Container, data=BeetlesMale, FUN=sum)



test_that("wrapper function rpt works for all distributions without perm or boot and defaults", {
        expect_error(rpt(BodyL ~ (1|Population), grname="Population", data=BeetlesBody, 
                nboot=0, npermut=0, datatype = "Gaussian"), NA)
        expect_error(rpt(Egg ~ Treatment + (1|Container), grname=c("Container"), data = BeetlesFemale,
                nboot=0, npermut=0, datatype = "Poisson"), NA)
        expect_error(rpt(Colour ~ (1|Population), grname=c("Population"), 
                data=BeetlesMale, nboot=5, npermut=0, datatype = "Binary"), NA)
        expect_error(rpt(cbind(Dark, Reddish) ~ (1|Population), grname=c("Population"), data=md,
                nboot=0, npermut=0, datatype = "Proportion"), NA)
})


test_that("wrapper function rpt works for all distributions with perm and boot and defaults", {
        expect_error(rpt(BodyL ~ (1|Population), grname="Population", data=BeetlesBody, 
                nboot=2, npermut=2, datatype = "Gaussian"), NA)
        expect_error(rpt(Egg ~ Treatment + (1|Container), grname=c("Container"), data = BeetlesFemale,
                nboot=2, npermut=2, datatype = "Poisson"), NA)
        expect_error(rpt(Colour ~ (1|Population), grname=c("Population"), 
                data=BeetlesMale, nboot=2, npermut=2, datatype = "Binary"), NA)
        expect_error(rpt(cbind(Dark, Reddish) ~ (1|Population), grname=c("Population"), data=md,
                nboot=2, npermut=2, datatype = "Proportion"), NA)
})



