context("rpt")

# Set a seed for reproducibility of the randomization 
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


# test random slopes

#res1 <- rpt(BodyL ~ Treatment + Sex + (1|Population) + (0 + Treatment|Population),  grname=c("Population"), data=BeetlesBody, datatype="Gaussian", nboot=0, npermut=0)
#res2 <- rpt(BodyL ~ Treatment + Sex + (1|Population) + (-1 + Treatment|Population),  grname=c("Population"), data=BeetlesBody, datatype="Gaussian", nboot=0, npermut=0)
#res3 <- rpt(BodyL ~ Treatment + Sex + (1|Population) + (- 1 + Treatment|Population),  grname=c("Population"), data=BeetlesBody, datatype="Gaussian", nboot=0, npermut=0)
#res4 <- rpt(BodyL ~ Treatment + Sex + (0 + Treatment|Population) + (1|Population),  grname=c("Population"), data=BeetlesBody, datatype="Gaussian", nboot=0, npermut=0)
#res5 <- rpt(BodyL ~ Treatment + Sex + (-1 + Treatment|Population) + (1|Population),  grname=c("Population"), data=BeetlesBody, datatype="Gaussian", nboot=0, npermut=0)
#res6 <- rpt(BodyL ~ Treatment + Sex + (- 1 + Treatment|Population) + (1|Population) ,  grname=c("Population"), data=BeetlesBody, datatype="Gaussian", nboot=0, npermut=0)
