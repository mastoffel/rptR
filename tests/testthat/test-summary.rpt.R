context("summary.rpt")

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


# summary without bootstraps or permutations
rpt_gaussian <- rpt(BodyL ~ (1|Population), grname="Population", data=BeetlesBody, 
        nboot=0, npermut=0, datatype = "Gaussian")
rpt_poisson <- rpt(Egg ~ Treatment + (1|Container), grname=c("Container"), data = BeetlesFemale,
        nboot=0, npermut=0, datatype = "Poisson")
rpt_binary <- rpt(Colour ~ (1|Population), grname=c("Population"), 
        data=BeetlesMale, nboot=5, npermut=0, datatype = "Binary")
rpt_proportion <- rpt(cbind(Dark, Reddish) ~ (1|Population), grname=c("Population"), data=md,
                nboot=0, npermut=0, datatype = "Proportion")



test_that("wrapper function rpt works for all distributions without perm or boot and defaults", {
        expect_error(test <- summary(rpt_gaussian), NA)
        expect_error(rpt(Egg ~ Treatment + (1|Container), grname=c("Container"), data = BeetlesFemale,
                nboot=0, npermut=0, datatype = "Poisson"), NA)
        expect_error(rpt(Colour ~ (1|Population), grname=c("Population"), 
                data=BeetlesMale, nboot=5, npermut=0, datatype = "Binary"), NA)
        expect_error(rpt(cbind(Dark, Reddish) ~ (1|Population), grname=c("Population"), data=md,
                nboot=0, npermut=0, datatype = "Proportion"), NA)
})