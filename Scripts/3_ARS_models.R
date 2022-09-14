# # ------------------------------------------------------------------------

# ---- Parental age manuscript 
# ---- AMS
# ---- ARS models

# # ------------------------------------------------------------------------

rm(list=ls())


#~~ libraries

library(glmmTMB)
library(arm)
library(DHARMa)
library(performance)
library(data.table)

#~~ data 

ARSdata <- read.csv("Data/Sparks_Parental_Age_ARS_Data_2022.csv", stringsAsFactors = FALSE)


#~~ set factors etc

ARSdata$Age <- ARSdata$Year - ARSdata$BirthYear
ARSdata$AgeSq <- ARSdata$Age^2
ARSdata$logTQ <- log10(ARSdata$TQ)
ARSdata$HelpersF <- as.factor(ARSdata$HelpersF)
ARSdata$Siblings <- as.factor(ARSdata$Siblings)
ARSdata$MotherID<-as.factor(ARSdata$MotherID)
ARSdata$FatherID<-as.factor(ARSdata$FatherID)
ARSdata$DominantMaleID<-as.factor(ARSdata$DominantMaleID)
ARSdata$ID<-as.factor(ARSdata$ID)
ARSdata$BirthYearF<-as.factor(ARSdata$BirthYear)
ARSdata$Sex<-as.factor(ARSdata$Sex)
ARSdata$LifespanInt <- ARSdata$FinalYear - ARSdata$BirthYear



# # -----------------------------------------------------------------------


# ARS models: female offspring


# # -----------------------------------------------------------------------

#~~ model subset
#~~ NOTE: the dataset is already cropped to individuals who made it to at least 1 year,
#~~ who were not translocated and whose final year was < 2018 (last year of the pedigree) and birth year > 1994

ARSdataf <- droplevels(subset(ARSdata, Sex==0 & !is.na(MotherID) & !is.na(FatherID) & !is.na(DominantMaleID) & !is.na(logTQ) 
                              & !is.na(HelpersF) & !is.na(GroupSize)))


#~~ data

table(ARSdataf$ARS1yr)


#~~ model - run a variety of GLMMs in glmmTMB and see which model fits best

f.ars.zipoisson <- glmmTMB(ARS1yr~ rescale(Age) + rescale(AgeSq) + rescale(LifespanInt) + rescale(AgeMother) + rescale(AgeFather) + 
                            rescale(DominantMaleAge) +
                            rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                            Siblings +
                            (1|ID) + (1|Year),
                          data=ARSdataf,
                          ziformula=~1,
                          family=poisson)
summary(f.ars.zipoisson)

f.ars.zinbinom2 <- update(f.ars.zipoisson, family=nbinom2) # did not converge

f.ars.zinbinom1 <- update(f.ars.zipoisson, family=nbinom1) # did not converge

f.ars.poisson <- update(f.ars.zipoisson, ziformula=~0, family=poisson)

f.ars.nbinom2 <- update(f.ars.zipoisson, ziformula=~0, family=nbinom2) # did not converge

f.ars.nbinom1 <- update(f.ars.zipoisson, ziformula=~0, family=nbinom1) # did not converge


#~~ AIC comparisons - NOTE only poisson and ZI poisson models converged (AIC of poisson lower - just)

AIC(f.ars.zinbinom1, f.ars.zinbinom2, f.ars.poisson, f.ars.nbinom2, f.ars.nbinom1, f.ars.zipoisson)


#~~ check model residuals 

f.ars.poisson.resid <- simulateResiduals(fittedModel = f.ars.poisson, plot = T) # pass
f.ars.zipoisson.resid <- simulateResiduals(fittedModel = f.ars.zipoisson, plot = T) # pass


#~~ look at zero inflation component

summary(f.ars.zipoisson) # ZI intercept not sig
summary(f.ars.poisson) # model estimates very similar between both models

check_zeroinflation(f.ars.poisson) # ratio of observed to predicted zeroes is okay
check_zeroinflation(f.ars.zipoisson)


#~~ go with poisson model

f.ars.model <- glmmTMB(ARS1yr~ rescale(Age) + rescale(AgeSq) + rescale(LifespanInt) + rescale(AgeMother) + rescale(AgeFather) + 
                          rescale(DominantMaleAge) +
                          rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                          Siblings +
                          (1|ID) + (1|Year),
                        data=ARSdataf,
                        ziformula=~0,
                        family=poisson)
summary(f.ars.model)


#~~ model performance checks

check_overdispersion(f.ars.model)
check_zeroinflation(f.ars.model)
check_collinearity(f.ars.model)
f.ars.model.resid <- simulateResiduals(fittedModel = f.ars.model, plot = T)


#~~ save model output

ARSfmodelfixef<-as.data.frame((summary(f.ars.model)$coefficients$cond))
ARSfmodelfixef<-setDT(ARSfmodelfixef, keep.rownames = TRUE)[]

write.csv(ARSfmodelfixef, file="Models/ARS_female_model.csv")


#~~ significance of fixed effects

#~~ age2 - x=211.53  df=1  p<2.2e-16 

f.ars.model.2 <- update(f.ars.model, ~ . -rescale(AgeSq))

anova(f.ars.model, f.ars.model.2)


#~~ lifespan - x=8.4071  df=1  p=0.003738

f.ars.model.3 <- update(f.ars.model, ~ . -rescale(LifespanInt))

anova(f.ars.model, f.ars.model.3)


#~~ age of mother - x=0.13  df=1  p=0.7184

f.ars.model.4 <- update(f.ars.model, ~ . -rescale(AgeMother))

anova(f.ars.model, f.ars.model.4)


#~~ age of father - x=0.0158  df=1  p=0.9001

f.ars.model.5 <- update(f.ars.model, ~ . -rescale(AgeFather))

anova(f.ars.model, f.ars.model.5)


#~~ age of dominant male - x=1.4163  df=1  p=0.234

f.ars.model.6 <- update(f.ars.model, ~ . -rescale(DominantMaleAge))

anova(f.ars.model, f.ars.model.6)


#~~ birth year - x=0.0299  df=1  p=0.8627

f.ars.model.7 <- update(f.ars.model, ~ . -rescale(BirthYear))

anova(f.ars.model, f.ars.model.7)


#~~ tq - x=0.8776  df=1  p=0.3489

f.ars.model.8 <- update(f.ars.model, ~ . -rescale(logTQ))

anova(f.ars.model, f.ars.model.8)


#~~ group size - x=0.6339  df=1  p=0.4259

f.ars.model.9 <- update(f.ars.model, ~ . -rescale(GroupSize))

anova(f.ars.model, f.ars.model.9)


#~~ helpers - x=5.1086  df=1  p=0.02381

f.ars.model.10 <- update(f.ars.model, ~ . -HelpersF)

anova(f.ars.model, f.ars.model.10)


#~~ siblings - x=0.0817  df=1  p=0.775

f.ars.model.11 <- update(f.ars.model, ~ . -Siblings)

anova(f.ars.model, f.ars.model.11)



# # -----------------------------------------------------------------------


# ARS MODELS - MALES


# # -----------------------------------------------------------------------

#~~ model subset
#~~ NOTE: the dataset is already cropped to individuals who made it to at least 1 year,
#~~ who were not translocated and whose final year was < 2018 (last year of the pedigree) and birth year > 1994

ARSdatam <- droplevels(subset(ARSdata, Sex==1  & !is.na(MotherID) & !is.na(FatherID) & !is.na(DominantMaleID) & !is.na(logTQ) 
                              & !is.na(HelpersF) & !is.na(GroupSize)))


#~~ data

table(ARSdataf$ARS1yr)


#~~ model - run a variety of GLMMs in glmmTMB and see which model fits best

m.ars.zipoisson <- glmmTMB(ARS1yr~ rescale(Age) + rescale(AgeSq) + rescale(LifespanInt) + rescale(AgeMother) + rescale(AgeFather) +
                            rescale(DominantMaleAge) +
                            rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                            Siblings +
                            (1|ID) + (1|Year),
                          data=ARSdatam,
                          ziformula=~1,
                          family=poisson)
summary(m.ars.zipoisson)

m.ars.zinbinom2 <- update(m.ars.zipoisson, family=nbinom2)

m.ars.zinbinom1 <- update(m.ars.zipoisson, family=nbinom1)

m.ars.poisson <- update(m.ars.zipoisson, ziformula=~0, family=poisson)

m.ars.nbinom2 <- update(m.ars.zipoisson, ziformula=~0, family=nbinom2)

m.ars.nbinom1 <- update(m.ars.zipoisson, ziformula=~0, family=nbinom1)


#~~ AIC comparisons (poisson is the lowest, just)

AIC(m.ars.zinbinom1, m.ars.zinbinom2, m.ars.poisson, m.ars.nbinom2, m.ars.nbinom1, m.ars.zipoisson)


#~~ check model residuals - all pass

m.ars.zinbinom1.resid <- simulateResiduals(fittedModel = m.ars.zinbinom1, plot = T)

m.ars.nbinom1.resid <- simulateResiduals(fittedModel = m.ars.nbinom1, plot = T)

m.ars.poisson.resid <- simulateResiduals(fittedModel = m.ars.poisson, plot = T, n=1000)

m.ars.zipoisson.resid <- simulateResiduals(fittedModel = m.ars.zipoisson, plot = T)

m.ars.zinbinom2.resid <- simulateResiduals(fittedModel = m.ars.zinbinom2, plot = T)

m.ars.nbinom2.resid <- simulateResiduals(fittedModel = m.ars.nbinom2, plot = T)


#~~ conclusions the same for each model type

summary(m.ars.zipoisson)
summary(m.ars.zinbinom1)
summary(m.ars.zinbinom2)
summary(m.ars.nbinom1)
summary(m.ars.poisson)
summary(m.ars.nbinom2)


#~~ check model diagnostics for the poisson given it has the lowest AIC 

check_overdispersion(m.ars.poisson) # no overdispersion
check_zeroinflation(m.ars.poisson) # no zero inflation
check_collinearity(m.ars.poisson) 

m.ars.poisson.resid <- simulateResiduals(fittedModel = m.ars.poisson, plot = T, n=1000)
outliertest <- testOutliers(simulationOutput, type="bootstrap", nBoot=1000)


#~~ go with poisson model

m.ars.model <- glmmTMB(ARS1yr~ rescale(Age) + rescale(AgeSq) + rescale(LifespanInt) + rescale(AgeMother) + rescale(AgeFather) + 
                        rescale(DominantMaleAge) +
                        rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                        Siblings +
                        (1|ID) + (1|Year),
                      data=ARSdatam,
                      ziformula=~0,
                      family=poisson)
summary(m.ars.model)


#~~ model performance checks

check_overdispersion(m.ars.model)
check_zeroinflation(m.ars.model)
check_collinearity(m.ars.model)
m.ars.poisson.resid <- simulateResiduals(fittedModel = m.ars.poisson, plot = T, n=1000)


#~~ save model output

ARSmmodelfixef<-as.data.frame((summary(m.ars.model)$coefficients$cond))
ARSmmodelfixef<-setDT(ARSmmodelfixef, keep.rownames = TRUE)[]

write.csv(ARSmmodelfixef, file="Models/ARS_male_model.csv")


#~~ significance of fixed effects

#~~ age2 - x=213.39  df=1  p<2.2e-16 

m.ars.model.2 <- update(m.ars.model, ~ . -rescale(AgeSq))

anova(m.ars.model, m.ars.model.2)


#~~ lifespan - x=7.0173  df=1  p=0.008073

m.ars.model.3 <- update(m.ars.model, ~ . -rescale(LifespanInt))

anova(m.ars.model, m.ars.model.3)


#~~ age of mother - x=1.4612  df=1  p=0.2267

m.ars.model.4 <- update(m.ars.model, ~ . -rescale(AgeMother))

anova(m.ars.model, m.ars.model.4)


#~~ age of father - x=1.2558  df=1  p=0.2625

m.ars.model.5 <- update(m.ars.model, ~ . -rescale(AgeFather))

anova(m.ars.model, m.ars.model.5)


#~~ age of dominant male - x=2.032  df=1  p=0.154

m.ars.model.6 <- update(m.ars.model, ~ . -rescale(DominantMaleAge))

anova(m.ars.model, m.ars.model.6)


#~~ birth year - x=1.0621  df=1  p=0.3027

m.ars.model.7 <- update(m.ars.model, ~ . -rescale(BirthYear))

anova(m.ars.model, m.ars.model.7)


#~~ tq - x=1.035  df=1  p=0.309

m.ars.model.8 <- update(m.ars.model, ~ . -rescale(logTQ))

anova(m.ars.model, m.ars.model.8)


#~~ group size - x=2.8784  df=1  p=0.08977

m.ars.model.9 <- update(m.ars.model, ~ . -rescale(GroupSize))

anova(m.ars.model, m.ars.model.9)


#~~ helpers - x=0.005  df=1  p=0.9436

m.ars.model.10 <- update(m.ars.model, ~ . -HelpersF)

anova(m.ars.model, m.ars.model.10)


#~~ siblings - x=0.1828  df=1  p=0.6689

m.ars.model.11 <- update(m.ars.model, ~ . -Siblings)

anova(m.ars.model, m.ars.model.11)

