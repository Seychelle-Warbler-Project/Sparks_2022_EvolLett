# # ------------------------------------------------------------------------

# ---- Parental age manuscript
# ---- AMS
# ---- LRS models

# # ------------------------------------------------------------------------

rm(list=ls())


#~~ libraries

library(glmmTMB)
library(arm)
library(DHARMa)
library(data.table)
library(performance)
library(dplyr)


#~~  data

fitnessdata <- read.csv("Data/Sparks_Parental_Age_Lifespan_Data_2022.csv", stringsAsFactors = FALSE)


#~~ set factors etc

fitnessdata$MotherID <- as.factor(fitnessdata$MotherID)
fitnessdata$FatherID <- as.factor(fitnessdata$FatherID)
fitnessdata$DominantMaleID <- as.factor(fitnessdata$DominantMaleID)
fitnessdata$ID <- as.factor(fitnessdata$ID)
fitnessdata$BirthYearF <- as.factor(fitnessdata$BirthYear)
fitnessdata$Sex <- as.factor(fitnessdata$Sex)
fitnessdata$HelpersF <- as.factor(fitnessdata$HelpersF)
fitnessdata$TranslocatedF  <-  as.factor(fitnessdata$TranslocatedF)
fitnessdata$logTQ  <-  log10(fitnessdata$TQ)
fitnessdata$AgeMotherSq  <-  fitnessdata$AgeMother * fitnessdata$AgeMother
fitnessdata$AgeFatherSq  <-  fitnessdata$AgeFather * fitnessdata$AgeFather
fitnessdata$AgeSocFatherSq  <-  fitnessdata$DominantMaleAge * fitnessdata$DominantMaleAge
fitnessdata$Siblings  <-  as.factor(fitnessdata$Siblings)
fitnessdata$LifespanInt  <-  fitnessdata$FinalYear - fitnessdata$BirthYear


#~~ LRS data subset
#~~ NOTE: the dataset is already cropped to individuals who made it to at least 1 year
#~~ Lifetime reproductive success (LRS) measures are only available for birds that have completed their reproductive lifespan 
#~~ i.e. who are dead and whose final year was < 2018 (last complete year of the pedigree was 2017) and birth year > 1994
#~~ LRS is calculated irrespective of whether the individual became a dominant breeder or not
#~~ However, individuals who were translocated need to be removed from the dataset as they have incomplete LRS information

LRSdata <- droplevels(subset(fitnessdata, LifespanInt > 0 & !is.na(Sex) & FinalYear < 2018 & TranslocatedF==0))

str(LRSdata)

summary(LRSdata)


# # -----------------------------------------------------------------------

## LRS models: female offspring  

# # -----------------------------------------------------------------------

#~~ female model subset

LRSfdata <- droplevels(subset(LRSdata, Sex==0 & !is.na(GroupSize) & !is.na(AgeMother) & !is.na(AgeFather) & !is.na(DominantMaleAge) 
                              & !is.na(BirthYear) & !is.na(Siblings) & !is.na(HelpersF)))


#~~ data

table(LRSfdata$LRS1YEAR) 

hist(LRSfdata$LRS1YEAR)

summary(LRSfdata$LRS1YEAR)


#~~ model - run a variety of GLMMs in glmmTMB and see which model fits best

f.lrs.zipoisson <- glmmTMB(LRS1YEAR~ rescale(AgeMother) + rescale(AgeFather) + rescale(DominantMaleAge) +
                            rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                            Siblings + (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                          data=LRSfdata,
                          ziformula=~1,
                          family=poisson)

summary(f.lrs.zipoisson)

f.lrs.zinbinom2 <- update(f.lrs.zipoisson, family=nbinom2)

f.lrs.zinbinom1 <- update(f.lrs.zipoisson, family=nbinom1)

f.lrs.poisson <- update(f.lrs.zipoisson, ziformula=~0, family=poisson)

f.lrs.nbinom2 <- update(f.lrs.zipoisson, ziformula=~0, family=nbinom2)

f.lrs.nbinom1 <- update(f.lrs.zipoisson, ziformula=~0, family=nbinom1)

# all models converged


#~~ AIC comparisons
# lowest AIC is nbinom1 with zero inflation, then nbinom2 with zero inflation then zero inflated poisson

AIC(f.lrs.zinbinom1, f.lrs.zinbinom2, f.lrs.poisson, f.lrs.nbinom2, f.lrs.nbinom1, f.lrs.zipoisson)


#~~ check model residuals

f.lrs.zinbinom1.resid <- simulateResiduals(fittedModel = f.lrs.zinbinom1, plot = T) # fine

f.lrs.zinbinom2.resid <- simulateResiduals(fittedModel = f.lrs.zinbinom2, plot = T) # fine

f.lrs.poisson.resid <- simulateResiduals(fittedModel = f.lrs.poisson, plot = T) # not fine

f.lrs.nbinom2.resid <- simulateResiduals(fittedModel = f.lrs.nbinom2, plot = T) # not fine

f.lrs.nbinom1.resid <- simulateResiduals(fittedModel = f.lrs.nbinom1, plot = T) # not fine

f.lrs.zipoisson.resid <- simulateResiduals(fittedModel = f.lrs.zipoisson, plot = T) # fine


#~~ check zero inflation with 3 best models - ZI poisson accounts for observed zeroes best

check_zeroinflation(f.lrs.zinbinom1)
check_zeroinflation(f.lrs.zinbinom2)
check_zeroinflation(f.lrs.zipoisson)


#~~ go with poisson model
#~~ full model with linear and quadratic mother and father ages

lrs.model.f.1 <- glmmTMB(LRS1YEAR~ rescale(AgeMother) + rescale(AgeFather) + rescale(DominantMaleAge) +
                           rescale(AgeMotherSq) + rescale(AgeFatherSq) + rescale(AgeSocFatherSq) +
                            rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                            Siblings + (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                          data=LRSfdata,
                          ziformula=~1,
                          family=poisson)

summary(lrs.model.f.1)


#~~ remove non significant quadratic terms

#~~ drop age mother sq

lrs.model.f.2 <- update(lrs.model.f.1, ~ . -rescale(AgeMotherSq))

anova(lrs.model.f.1, lrs.model.f.2)

summary(lrs.model.f.2)


#~~ drop age dominant male sq

lrs.model.f.3 <- update(lrs.model.f.2, ~ . -rescale(AgeFatherSq))

anova(lrs.model.f.2, lrs.model.f.3)

summary(lrs.model.f.3)


#~~ drop age father sq

lrs.model.f.4 <- update(lrs.model.f.3, ~ . -rescale(AgeSocFatherSq))

anova(lrs.model.f.3, lrs.model.f.4)

summary(lrs.model.f.4)


#~~ base model

summary(lrs.model.f.4) # fine - parameter values > 10 are suspect


#~~ model performance checks

lrs.model.f.4.resid <- simulateResiduals(fittedModel = lrs.model.f.4, plot = T) 
check_collinearity(lrs.model.f.4)
check_zeroinflation(lrs.model.f.4)


#~~ save model output

LRSfmodelfixef<-as.data.frame((summary(lrs.model.f.4)$coefficients$cond))
LRSfmodelfixef<-setDT(LRSfmodelfixef, keep.rownames = TRUE)[]

write.csv(LRSfmodelfixef, file="Models/LRS_female_model.csv")


#~~ significance of fixed effects

#~~ age of mother - x=7.8169  df=1  p=0.005176

lrs.model.f.5 <- update(lrs.model.f.4, ~ . -rescale(AgeMother))

anova(lrs.model.f.4, lrs.model.f.5)


#~~ age of father - x=0.7979  df=1  p=0.3717

lrs.model.f.6 <- update(lrs.model.f.4, ~ . -rescale(AgeFather))

anova(lrs.model.f.4, lrs.model.f.6)


#~~ age of dominant male - x=0.2762  df=1  p=0.5992

lrs.model.f.7 <- update(lrs.model.f.4, ~ . -rescale(DominantMaleAge))

anova(lrs.model.f.4, lrs.model.f.7)


#~~ birth year - x=25.568  df=1  p=4.27e-07 

lrs.model.f.8 <- update(lrs.model.f.4, ~ . -rescale(BirthYear))

anova(lrs.model.f.4, lrs.model.f.8)


#~~ TQ - x=0.8922  df=1  p=0.3449

lrs.model.f.9 <- update(lrs.model.f.4, ~ . -rescale(logTQ))

anova(lrs.model.f.4, lrs.model.f.9)


#~~ group size - x=0.0461  df=1  p=0.83

lrs.model.f.10 <- update(lrs.model.f.4, ~ . -rescale(GroupSize))

anova(lrs.model.f.4, lrs.model.f.10)


#~~ helpers - x=10.055  df=1  p=0.001519

lrs.model.f.11 <- update(lrs.model.f.4, ~ . -HelpersF)

anova(lrs.model.f.4, lrs.model.f.11)


#~~ siblings - x=0.1668  df=1  p=0.683

lrs.model.f.12 <- update(lrs.model.f.4, ~ . -Siblings)

anova(lrs.model.f.4, lrs.model.f.12)


#~~ age mother sq - x=0.6464  df=1  p=0.4214

lrs.model.f.13 <- update(lrs.model.f.4, ~ . +rescale(AgeMotherSq))

anova(lrs.model.f.4, lrs.model.f.13)


#~~ age father sq - x=0.4561  df=1  p=0.4995

lrs.model.f.14 <- update(lrs.model.f.4, ~ . +rescale(AgeFatherSq))

anova(lrs.model.f.4, lrs.model.f.14)


#~~ dominant male age sq - x=0.4803  df=1  p=0.4883

lrs.model.f.15 <- update(lrs.model.f.4, ~ . +rescale(AgeSocFatherSq))

anova(lrs.model.f.4, lrs.model.f.15)


#~~ significance of interactions

#~~ TQ * mum age - x=0.6255  df=1  p=0.429

lrs.model.f.16 <- update(lrs.model.f.4, ~ . +rescale(logTQ)*rescale(AgeMother))

anova(lrs.model.f.4, lrs.model.f.16)


#~~ TQ * dad age - x=6.0852  df=1p=0.01363 

lrs.model.f.17 <- update(lrs.model.f.4, ~ . +rescale(logTQ)*rescale(AgeFather))

anova(lrs.model.f.4, lrs.model.f.17)

summary(lrs.model.f.17)


#~~ TQ * dominant male age - x=4.6984  df=1  p=0.03019

lrs.model.f.18 <- update(lrs.model.f.4, ~ . +rescale(logTQ)*rescale(DominantMaleAge))

anova(lrs.model.f.4, lrs.model.f.18)

summary(lrs.model.f.18)


#~~ group size * mum age - x=0.9509  df=1  p=0.3295

lrs.model.f.19 <- update(lrs.model.f.4, ~ . +rescale(GroupSize)*rescale(AgeMother))

anova(lrs.model.f.4, lrs.model.f.19)


#~~ group size * dad age - x=0.0455  df=1  p=0.8312

lrs.model.f.20 <- update(lrs.model.f.4, ~ . +rescale(GroupSize)*rescale(AgeFather))

anova(lrs.model.f.4, lrs.model.f.20)


#~~ group size * dominant male age - x=1.7727  df=1  p=0.1831

lrs.model.f.21 <- update(lrs.model.f.4, ~ . +rescale(GroupSize)*rescale(DominantMaleAge))

anova(lrs.model.f.4, lrs.model.f.21)


#~~ helpers * mum age - x=0.5297  df=1  p=0.4667

lrs.model.f.22 <- update(lrs.model.f.4, ~ . +HelpersF*rescale(AgeMother))

anova(lrs.model.f.4, lrs.model.f.22)


#~~ helpers * dad age - x=0.0125  df=1  p=0.9109

lrs.model.f.23 <- update(lrs.model.f.4, ~ . +HelpersF*rescale(AgeFather))

anova(lrs.model.f.4, lrs.model.f.23)


#~~ helpers * dominant male age - x=1.9956  df=1  p=0.1578

lrs.model.f.24 <- update(lrs.model.f.4, ~ . +HelpersF*rescale(DominantMaleAge))

anova(lrs.model.f.4, lrs.model.f.24)


#~~ is maternal age still significant with lifespan in model?

lrs.model.f.25 <- update(lrs.model.f.4, ~ . +rescale(LifespanInt))

summary(lrs.model.f.25)


#~~ lifespan - x=204.93  df=1  p<2.2e-16

anova(lrs.model.f.4, lrs.model.f.25)


#~~ maternal age - x=0.0567  df=1  p=0.8117

lrs.model.f.26 <- update(lrs.model.f.25, ~ . -rescale(AgeMother))

anova(lrs.model.f.26, lrs.model.f.25)


#~~ check results are consistent if age of father and age of dominant male included separately 

#~~ genetic parents

lrs.model.f.a <- glmmTMB(LRS1YEAR~ rescale(AgeMother) + rescale(AgeFather) + 
                           rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                           Siblings + (1|MotherID) + (1|FatherID) + (1|BirthYearF),
                         data=LRSfdata,
                         ziformula=~1,
                         family=poisson)
summary(lrs.model.f.a)


#~~ significance of maternal age - x=7.8941  df=1  p=0.00496

lrs.model.f.a1 <- update(lrs.model.f.a, ~ . -rescale(AgeMother))

anova(lrs.model.f.a1, lrs.model.f.a)


#~~ significance of paternal age - x=0.5728  df=1  p=0.4491

lrs.model.f.a2 <- update(lrs.model.f.a, ~ . -rescale(AgeFather))

anova(lrs.model.f.a2, lrs.model.f.a)


#~~ model with genetic mother and dominant male

lrs.model.f.b <- glmmTMB(LRS1YEAR~ rescale(AgeMother) + rescale(DominantMaleAge) + 
                           rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                           Siblings + (1|MotherID) + (1|DominantMaleID) + (1|BirthYearF),
                         data=LRSfdata,
                         ziformula=~1,
                         family=poisson)
summary(lrs.model.f.b)


#~~ significance of maternal age - x=7.9467  df=1  p=0.004818

lrs.model.f.b1 <- update(lrs.model.f.b, ~ . -rescale(AgeMother))

anova(lrs.model.f.b1, lrs.model.f.b)


#~~ significance of dominant male age - x=0.0816  df=1  p=0.7751

lrs.model.f.b2 <- update(lrs.model.f.b, ~ . -rescale(DominantMaleAge))

anova(lrs.model.f.b2, lrs.model.f.b)


#~~ check maternal age * paternal age interaction

#~~ maternal age * paternal age - x=0.3966  df=1  p=0.5288

lrs.model.f.c <- update(lrs.model.f.4, ~ . + rescale(AgeMother)*rescale(AgeFather))

summary(lrs.model.f.c)

anova(lrs.model.f.4, lrs.model.f.c)


#~~ maternal age * dominant male age - x=1.5031  df=1  p=0.2202

lrs.model.f.d <- update(lrs.model.f.4, ~ . + rescale(AgeMother)*rescale(DominantMaleAge))

summary(lrs.model.f.d)

anova(lrs.model.f.4, lrs.model.f.d)


# # -----------------------------------------------------------------------

## LRS models: female offspring

## Within-subject centering

# # -----------------------------------------------------------------------

#~~ set up variables

#~~ z-transform covariates

LRSfdata$cAgeMother <- rescale(LRSfdata$AgeMother)
LRSfdata$cAgeFather <- rescale(LRSfdata$AgeFather)
LRSfdata$cDominantMaleAge <- rescale(LRSfdata$DominantMaleAge)
LRSfdata$cBirthYear <- rescale(LRSfdata$BirthYear)
LRSfdata$clogTQ <- rescale(LRSfdata$logTQ)
LRSfdata$cGroupSize <- rescale(LRSfdata$GroupSize)


#~~ make mean maternal, paternal and dominant male age variables and deviation of each value from the mean

LRSfdata<-LRSfdata %>% group_by(MotherID) %>% mutate(cMeanAgeMother = mean(cAgeMother)) %>% ungroup()
LRSfdata$cDevMeanAgeMother<-LRSfdata$cAgeMother-LRSfdata$cMeanAgeMother

LRSfdata<-LRSfdata %>% group_by(FatherID) %>% mutate(cMeanAgeFather = mean(cAgeFather)) %>% ungroup()
LRSfdata$cDevMeanAgeFather<-LRSfdata$cAgeFather-LRSfdata$cMeanAgeFather

LRSfdata<-LRSfdata %>% group_by(DominantMaleID) %>% mutate(cMeanAgeSocMale = mean(cDominantMaleAge)) %>% ungroup()
LRSfdata$cDevMeanAgeSocMale<-LRSfdata$cDominantMaleAge-LRSfdata$cMeanAgeSocMale

summary(LRSfdata$cMeanAgeMother)
summary(LRSfdata$cMeanAgeSocMale)
summary(LRSfdata$cMeanAgeFather)


#~~ are within vs between parental ages significant?

lrs.model.f.wb <- glmmTMB(LRS1YEAR~cMeanAgeMother + cDevMeanAgeMother +
                            cMeanAgeFather + cDevMeanAgeFather + 
                            cMeanAgeSocMale + cDevMeanAgeSocMale +
                            cBirthYear + clogTQ + cGroupSize + HelpersF + Siblings +
                            (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                          data=LRSfdata,
                          ziformula=~1,
                          family=poisson)

summary(lrs.model.f.wb)


#~~ save model output

LRSfmodelwbfixef<-as.data.frame((summary(lrs.model.f.wb)$coefficients$cond))
LRSfmodelwbfixef<-setDT(LRSfmodelwbfixef, keep.rownames = TRUE)[]

write.csv(LRSfmodelwbfixef, file="Models/LRS_female_model_wnbn.csv")


#~~ significance of fixed effects

#~~ between maternal age effect - x=3.1969  df=1  p=0.07378 

lrs.model.f.wb.1 <- update(lrs.model.f.wb, ~ . -cMeanAgeMother)

anova(lrs.model.f.wb, lrs.model.f.wb.1)


#~~ within maternal age effect - x=9.6404  df=1  p=0.001903 

lrs.model.f.wb.2 <- update(lrs.model.f.wb, ~ . -cDevMeanAgeMother)

anova(lrs.model.f.wb, lrs.model.f.wb.2)


#~~ between paternal age effect - x=0.0669  df=1  p=0.7958

lrs.model.f.wb.3 <- update(lrs.model.f.wb, ~ . -cMeanAgeFather)

anova(lrs.model.f.wb, lrs.model.f.wb.3)


#~~ within paternal age effect - x=1.1033  df=1  p=0.2935

lrs.model.f.wb.4 <- update(lrs.model.f.wb, ~ . -cDevMeanAgeFather)

anova(lrs.model.f.wb, lrs.model.f.wb.4)


#~~ between dominant male age effect - x=0.6901  df=1  p=0.4061

lrs.model.f.wb.5 <- update(lrs.model.f.wb, ~ . -cMeanAgeSocMale)

anova(lrs.model.f.wb, lrs.model.f.wb.5)


#~~ within dominant male age effect - x=0.0067  df=1  p=0.935

lrs.model.f.wb.6 <- update(lrs.model.f.wb, ~ . -cDevMeanAgeSocMale)

anova(lrs.model.f.wb, lrs.model.f.wb.6)


#~~ birth year - x=23.569  df=1  p=1.205e-06

lrs.model.f.wb.7 <- update(lrs.model.f.wb, ~ . -cBirthYear)

anova(lrs.model.f.wb, lrs.model.f.wb.7)


#~~ territory quality - x=0.9269  df=1  p=0.3357

lrs.model.f.wb.8 <- update(lrs.model.f.wb, ~ . -clogTQ)

anova(lrs.model.f.wb, lrs.model.f.wb.8)


#~~ group size - x=6e-04  df=1  p=0.9806

lrs.model.f.wb.9 <- update(lrs.model.f.wb, ~ . -cGroupSize)

anova(lrs.model.f.wb, lrs.model.f.wb.9)


#~~ helpers - x=9.3356  df=1  p=0.002247

lrs.model.f.wb.10 <- update(lrs.model.f.wb, ~ . -HelpersF)

anova(lrs.model.f.wb, lrs.model.f.wb.10)


#~~ siblings - x=0.3325  df=1  p=0.5642

lrs.model.f.wb.11 <- update(lrs.model.f.wb, ~ . -Siblings)

anova(lrs.model.f.wb, lrs.model.f.wb.11)


#~~ are these within vs between slopes significantly different?

lrs.model.f.wb.b <- glmmTMB(LRS1YEAR~cMeanAgeMother + cAgeMother +
                              cMeanAgeFather + cAgeFather +
                              cMeanAgeSocMale + cDominantMaleAge +
                              cBirthYear + clogTQ + cGroupSize + HelpersF + Siblings +
                              (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                            data=LRSfdata,
                            ziformula=~1,
                            family=poisson)

summary(lrs.model.f.wb.b)


#~~ save model output

LRSfmodelwb.b.fixef<-as.data.frame((summary(lrs.model.f.wb.b)$coefficients$cond))
LRSfmodelwb.b.fixef<-setDT(LRSfmodelwb.b.fixef, keep.rownames = TRUE)[]

write.csv(LRSfmodelwb.b.fixef, file="Models/LRS_female_model_wnbn_b.csv")


#~~ significance of fixed effects

#~~ difference between within and between maternal age slopes - x=3.3436  df=1  p=0.06747 

lrs.model.f.wb.b2 <- update(lrs.model.f.wb.b, ~ . -cMeanAgeMother)

anova(lrs.model.f.wb.b, lrs.model.f.wb.b2)


#~~ within maternal age effect - x=9.6404  df=1  p=0.001903

lrs.model.f.wb.b3 <- update(lrs.model.f.wb.b, ~ . -cAgeMother)

anova(lrs.model.f.wb.b, lrs.model.f.wb.b3)


#~~ difference between within and between paternal age slopes - x=0.659  df=1  p=0.4169

lrs.model.f.wb.b4 <- update(lrs.model.f.wb.b, ~ . -cMeanAgeFather)

anova(lrs.model.f.wb.b, lrs.model.f.wb.b4)


#~~ within paternal age effect - x=1.1033  df=1  p=0.2935

lrs.model.f.wb.b5 <- update(lrs.model.f.wb.b, ~ . -cAgeFather)

anova(lrs.model.f.wb.b, lrs.model.f.wb.b5)


#~~ difference between within and between dominant male age slopes - x=0.2602  df=1  p=0.61

lrs.model.f.wb.b6 <- update(lrs.model.f.wb.b, ~ . -cMeanAgeSocMale)

anova(lrs.model.f.wb.b, lrs.model.f.wb.b6)


#~~ within dominant male age effect - x=0.0067  df=1  p=0.935

lrs.model.f.wb.b7 <- update(lrs.model.f.wb.b, ~ . -cDominantMaleAge)

anova(lrs.model.f.wb.b, lrs.model.f.wb.b7)


#~~ birth year - x=23.569  df=1  p=1.205e-06

lrs.model.f.wb.b8 <- update(lrs.model.f.wb.b, ~ . -cBirthYear)

anova(lrs.model.f.wb.b, lrs.model.f.wb.b8)


#~~ territory quality - x=0.9269  df=1  p=0.3357

lrs.model.f.wb.b9 <- update(lrs.model.f.wb.b, ~ . -clogTQ)

anova(lrs.model.f.wb.b, lrs.model.f.wb.b9)


#~~ group size - x=6e-04  df=1  p=0.9806

lrs.model.f.wb.b10 <- update(lrs.model.f.wb.b, ~ . -cGroupSize)

anova(lrs.model.f.wb.b, lrs.model.f.wb.b10)


#~~ helpers - x=9.3356     df=1   p=0.002247

lrs.model.f.wb.b11 <- update(lrs.model.f.wb.b, ~ . -HelpersF)

anova(lrs.model.f.wb.b, lrs.model.f.wb.b11)


#~~ siblings - x=0.3325  df=1  p=0.5642

lrs.model.f.wb.b12 <- update(lrs.model.f.wb.b, ~ . -Siblings)

anova(lrs.model.f.wb.b, lrs.model.f.wb.b12)


#~~ sanity check that maternal age estimates are similar if age of the father and age of the dominant male are not separated into a within vs between effect

lrs.model.f.wb.b1 <- glmmTMB(LRS1YEAR~cMeanAgeMother + cDevMeanAgeMother +
                               cAgeFather + cDominantMaleAge +
                               cBirthYear + clogTQ + cGroupSize + HelpersF + Siblings +
                               (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                             data=LRSfdata,
                             ziformula=~1,
                             family=poisson)

summary(lrs.model.f.wb.b1)

lrs.model.f.wb.b2 <- glmmTMB(LRS1YEAR~cMeanAgeMother + cAgeMother +
                               cAgeFather + cDominantMaleAge +
                               cBirthYear + clogTQ + cGroupSize + HelpersF + Siblings +
                               (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                             data=LRSfdata,
                             ziformula=~1,
                             family=poisson)

summary(lrs.model.f.wb.b2)


#~~ check model estimates similar without the outlier seen in the within maternal age graph

LRSfdata.check <- droplevels(subset(LRSfdata, DevMeanAgeMother > -6.9))

lrs.model.f.wb.check <- glmmTMB(LRS1YEAR~cMeanAgeMother + cDevMeanAgeMother +
                            cMeanAgeFather + cDevMeanAgeFather + 
                            cMeanAgeSocMale + cDevMeanAgeSocMale +
                            cBirthYear + clogTQ + cGroupSize + HelpersF + Siblings +
                            (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                          data=LRSfdata.check,
                          ziformula=~1,
                          family=poisson)

summary(lrs.model.f.wb.check)


#~~ save model output

LRSfmodelwbcheckfixef<-as.data.frame((summary(lrs.model.f.wb.check)$coefficients$cond))
LRSfmodelwbcheckfixef<-setDT(LRSfmodelwbcheckfixef, keep.rownames = TRUE)[]

write.csv(LRSfmodelwbcheckfixef, file="Models/LRS_female_model_wnbn_check.csv")


#~~ significance of fixed effects

#~~ between maternal age effect - x=3.079  df=1  p=0.07931

lrs.model.f.wb.check2 <- update(lrs.model.f.wb.check, ~ . -cMeanAgeMother)

anova(lrs.model.f.wb.check, lrs.model.f.wb.check2)


#~~ within maternal age effect - x=7.6947  df=1  p=0.005538 

lrs.model.f.wb.check3 <- update(lrs.model.f.wb.check, ~ . -cDevMeanAgeMother)

anova(lrs.model.f.wb.check, lrs.model.f.wb.check3)


#~~ between paternal age effect - x=0.0605  df=1  p=0.8057

lrs.model.f.wb.check4 <- update(lrs.model.f.wb.check, ~ . -cMeanAgeFather)

anova(lrs.model.f.wb.check, lrs.model.f.wb.check4)


#~~ within paternal age effect - x=1.1091  df=1  p=0.2923

lrs.model.f.wb.check5 <- update(lrs.model.f.wb.check, ~ . -cDevMeanAgeFather)

anova(lrs.model.f.wb.check, lrs.model.f.wb.check5)


#~~ between dominant male age effect - x=0.6885  df=1  p=0.4067

lrs.model.f.wb.check6 <- update(lrs.model.f.wb.check, ~ . -cMeanAgeSocMale)

anova(lrs.model.f.wb.check, lrs.model.f.wb.check6)


#~~ within dominant male age effect - x=0.0082  df=1  p=0.928

lrs.model.f.wb.check7 <- update(lrs.model.f.wb.check, ~ . -cDevMeanAgeSocMale)

anova(lrs.model.f.wb.check, lrs.model.f.wb.check7)


#~~ birth year - x=23.354  df=1  p=1.347e-06

lrs.model.f.wb.check8 <- update(lrs.model.f.wb.check, ~ . -cBirthYear)

anova(lrs.model.f.wb.check, lrs.model.f.wb.check8)


#~~ territory quality - x=0.8978  df=1  p=0.3434

lrs.model.f.wb.check9 <- update(lrs.model.f.wb.check, ~ . -clogTQ)

anova(lrs.model.f.wb.check, lrs.model.f.wb.check9)


#~~ group size - x=3e-04  df=1  p=0.986

lrs.model.f.wb.check10 <- update(lrs.model.f.wb.check, ~ . -cGroupSize)

anova(lrs.model.f.wb.check, lrs.model.f.wb.check10)


#~~ helpers - x=9.3469   df=1  p=0.002234

lrs.model.f.wb.check11 <- update(lrs.model.f.wb.check, ~ . -HelpersF)

anova(lrs.model.f.wb.check, lrs.model.f.wb.check11)


#~~ siblings - x=0.3384  df=1  p=0.5608

lrs.model.f.wb.check12 <- update(lrs.model.f.wb.check, ~ . -Siblings)

anova(lrs.model.f.wb.check, lrs.model.f.wb.check12)


# # -----------------------------------------------------------------------

## LRS models: male offspring  

# # -----------------------------------------------------------------------

#~~ male model subset

LRSmdata <- droplevels(subset(LRSdata, Sex==1 & !is.na(GroupSize) & !is.na(AgeMother) & !is.na(AgeFather) & !is.na(DominantMaleAge) 
                              & !is.na(BirthYear) & !is.na(Siblings) & !is.na(HelpersF)))


#~~ data

table(LRSmdata$LRS1YEAR) 

hist(LRSmdata$LRS1YEAR)

summary(LRSmdata$LRS1YEAR)


#~~ model - run a variety of GLMMs in glmmTMB and see which model fits best

m.lrs.zipoisson <- glmmTMB(LRS1YEAR~ rescale(AgeMother) + rescale(AgeFather) + rescale(DominantMaleAge) +
                             rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                             Siblings + (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                           data=LRSmdata,
                           ziformula=~1,
                           family=poisson)

summary(m.lrs.zipoisson)

m.lrs.zinbinom2 <- update(m.lrs.zipoisson, family=nbinom2)

m.lrs.zinbinom1 <- update(m.lrs.zipoisson, family=nbinom1)

m.lrs.poisson <- update(m.lrs.zipoisson, ziformula=~0, family=poisson)

m.lrs.nbinom2 <- update(m.lrs.zipoisson, ziformula=~0, family=nbinom2)

m.lrs.nbinom1 <- update(m.lrs.zipoisson, ziformula=~0, family=nbinom1)

# all models converged


#~~ AIC comparisons
# lowest AIC is nbinom1 with zero inflation

AIC(m.lrs.zinbinom1, m.lrs.zinbinom2, m.lrs.poisson, m.lrs.nbinom2, m.lrs.nbinom1, m.lrs.zipoisson)


#~~ check model diagnostics

m.lrs.zinbinom1.resid <- simulateResiduals(fittedModel = m.lrs.zinbinom1, plot = T) # fine

m.lrs.zinbinom2.resid <- simulateResiduals(fittedModel = m.lrs.zinbinom2, plot = T) # not fine

m.lrs.poisson.resid <- simulateResiduals(fittedModel = m.lrs.poisson, plot = T) # not fine

m.lrs.nbinom2.resid <- simulateResiduals(fittedModel = m.lrs.nbinom2, plot = T) # fine

m.lrs.nbinom1.resid <- simulateResiduals(fittedModel = m.lrs.nbinom1, plot = T) # fine

m.lrs.zipoisson.resid <- simulateResiduals(fittedModel = m.lrs.zipoisson, plot = T) # fine


#~~ check zero inflation with 4 best models - ZI poisson accounts for observed zeroes best

check_zeroinflation(m.lrs.zinbinom1)
check_zeroinflation(m.lrs.nbinom2)
check_zeroinflation(m.lrs.nbinom1)
check_zeroinflation(m.lrs.zipoisson)


#~~ full model with linear and quadratic mother and father ages

lrs.model.m.1 <- glmmTMB(LRS1YEAR~ rescale(AgeMother) + rescale(AgeFather) + rescale(DominantMaleAge) +
                           rescale(AgeMotherSq) + rescale(AgeFatherSq) + rescale(AgeSocFatherSq) +
                           rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                           Siblings + (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                         data=LRSmdata,
                         ziformula=~1,
                         family=poisson)

summary(lrs.model.m.1)


#~~ drop age mother sq

lrs.model.m.2 <- update(lrs.model.m.1, ~ . -rescale(AgeMotherSq))

anova(lrs.model.m.1,lrs.model.m.2)

summary(lrs.model.m.2)


#~~ drop dominant male age sq

lrs.model.m.3 <- update(lrs.model.m.2, ~ . -rescale(AgeSocFatherSq))

anova(lrs.model.m.3,lrs.model.m.2)

summary(lrs.model.m.3)


#~~ drop age father sq

lrs.model.m.4 <- update(lrs.model.m.3, ~ . -rescale(AgeFatherSq))

anova(lrs.model.m.3,lrs.model.m.4)


#~~ base model

summary(lrs.model.m.4) # fine - parameter values > 10 are suspect


#~~ model performance checks

lrs.model.m.4.resid <- simulateResiduals(fittedModel = lrs.model.m.4, plot = T) 
check_collinearity(lrs.model.m.4)
check_zeroinflation(lrs.model.m.4)


#~~ save output

LRSmmodelfixef<-as.data.frame((summary(lrs.model.m.4)$coefficients$cond))
LRSmmodelfixef<-setDT(LRSmmodelfixef, keep.rownames = TRUE)[]

write.csv(LRSmmodelfixef, file="Models/LRS_male_model.csv")


#~~ significance of fixed effects

#~~ age of mother - x=5.297  df=1  p=0.02136

lrs.model.m.5 <- update(lrs.model.m.4, ~ . -rescale(AgeMother))

anova(lrs.model.m.4, lrs.model.m.5)


#~~ age of father - x=2.2952  df=1  p=0.1298

lrs.model.m.6 <- update(lrs.model.m.4, ~ . -rescale(AgeFather))

anova(lrs.model.m.4, lrs.model.m.6)


#~~ dominant male age - x=0.0233  df=1  p=0.8788

lrs.model.m.7 <- update(lrs.model.m.4, ~ . -rescale(DominantMaleAge))

anova(lrs.model.m.4, lrs.model.m.7)


#~~ birth year - x=17.446  df=1  p=2.955e-05

lrs.model.m.8 <- update(lrs.model.m.4, ~ . -rescale(BirthYear))

anova(lrs.model.m.4, lrs.model.m.8)


#~~ TQ - x=1.4948  df=1  p=0.2215

lrs.model.m.9 <- update(lrs.model.m.4, ~ . -rescale(logTQ))

anova(lrs.model.m.4, lrs.model.m.9)


#~~ group size - x=0.7082  df=1  p=0.4

lrs.model.m.10 <- update(lrs.model.m.4, ~ . -rescale(GroupSize))

anova(lrs.model.m.4, lrs.model.m.10)


#~~ helpers - x=1.0689  df=1  p=0.3012

lrs.model.m.11 <- update(lrs.model.m.4, ~ . -HelpersF)

anova(lrs.model.m.4, lrs.model.m.11)


#~~ siblings - x=0.1558  df=1  p=0.6931

lrs.model.m.12 <- update(lrs.model.m.4, ~ . -Siblings)

anova(lrs.model.m.4, lrs.model.m.12)


#~~~ significance of quadratic parental effects

#~~ age mother sq - x=0.1865  df=1  p=0.6659

lrs.model.m.13 <- update(lrs.model.m.4, ~ . +rescale(AgeMotherSq))

anova(lrs.model.m.4, lrs.model.m.13)


#~~ age father sq - x=0.0779  df=1  p=0.7801

lrs.model.m.14 <- update(lrs.model.m.4, ~ . +rescale(AgeFatherSq))

anova(lrs.model.m.4, lrs.model.m.14)


#~~ dominant male age sq - x=0.0724  df=1  p=0.7879

lrs.model.m.15 <- update(lrs.model.m.4, ~ . +rescale(AgeSocFatherSq))

anova(lrs.model.m.4, lrs.model.m.15)


#~~~ significance of interactions

#~~ TQ * mum age - x=0.0816  df=1  p=0.7751

lrs.model.m.16 <- update(lrs.model.m.4, ~ . +rescale(logTQ)*rescale(AgeMother))

anova(lrs.model.m.4, lrs.model.m.16)


#~~ TQ * dad age - x=0.0086  df=1  p=0.9261

lrs.model.m.17 <- update(lrs.model.m.4, ~ . +rescale(logTQ)*rescale(AgeFather))

anova(lrs.model.m.4, lrs.model.m.17)


#~~ TQ * dominant male age - x=0.0097  df=1  p=0.9217

lrs.model.m.18 <- update(lrs.model.m.4, ~ . +rescale(logTQ)*rescale(DominantMaleAge))

anova(lrs.model.m.4, lrs.model.m.18)


#~~ group size * mum age - x=0.0902  df=1  p=0.7639

lrs.model.m.19 <- update(lrs.model.m.4, ~ . +rescale(GroupSize)*rescale(AgeMother))

anova(lrs.model.m.4, lrs.model.m.19)


#~~ group size * dad age - x=0.7135  df=1  p=0.3983

lrs.model.m.20 <- update(lrs.model.m.4, ~ . +rescale(GroupSize)*rescale(AgeFather))

anova(lrs.model.m.4, lrs.model.m.20)


#~~ group size * dominant male age - x=0.0858  df=1  p=0.7696

lrs.model.m.21 <- update(lrs.model.m.4, ~ . +rescale(GroupSize)*rescale(DominantMaleAge))

anova(lrs.model.m.4, lrs.model.m.21)


#~~ helpers * mum age - x=0.6121  df=1  p=0.434

lrs.model.m.22 <- update(lrs.model.m.4, ~ . +HelpersF*rescale(AgeMother))

anova(lrs.model.m.4, lrs.model.m.22)


#~~ helpers * dad age - x=0.7715  df=1  p=0.3798

lrs.model.m.23 <- update(lrs.model.m.4, ~ . +HelpersF*rescale(AgeFather))

anova(lrs.model.m.4, lrs.model.m.23)


#~~ helpers * dominant male age - x=0.0012  df=1  p=0.9721

lrs.model.m.24 <- update(lrs.model.m.4, ~ . +HelpersF*rescale(DominantMaleAge))

anova(lrs.model.m.4, lrs.model.m.24)


#~~ check results are consistent if age of father and age of dominant male included separately 

#~~ genetic parents

lrs.model.m.a <- glmmTMB(LRS1YEAR~ rescale(AgeMother) + rescale(AgeFather) + 
                           rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                           Siblings + (1|MotherID) + (1|FatherID) + (1|BirthYearF),
                         data=LRSmdata,
                         ziformula=~1,
                         family=poisson)
summary(lrs.model.m.a)


#~~ significance of maternal age - x=5.5545  df=1  p=0.01843

lrs.model.m.a1 <- update(lrs.model.m.a, ~ . -rescale(AgeMother))

anova(lrs.model.m.a1, lrs.model.m.a)


#~~ significance of paternal age - x=0.5728  df=1  p=0.4491

lrs.model.f.a2 <- update(lrs.model.f.a, ~ . -rescale(AgeFather))

anova(lrs.model.f.a2, lrs.model.f.a)


#~~ model with genetic mother and dominant male 

lrs.model.m.b <- glmmTMB(LRS1YEAR~ rescale(AgeMother) + rescale(DominantMaleAge) + 
                           rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                           Siblings + (1|MotherID) + (1|DominantMaleID) + (1|BirthYearF),
                         data=LRSmdata,
                         ziformula=~1,
                         family=poisson)
summary(lrs.model.m.b)


#~~ significance of maternal age - x=4.5819  df=1  p=0.03231

lrs.model.m.b1 <- update(lrs.model.m.b, ~ . -rescale(AgeMother))

anova(lrs.model.m.b, lrs.model.m.b1)


#~~ significance of dominant male age - x=0.0816  df=1  p=0.7751

lrs.model.m.b2 <- update(lrs.model.m.b, ~ . -rescale(DominantMaleAge))

anova(lrs.model.f.b2, lrs.model.f.b)


#~~ check maternal age * paternal age interaction

#~~ maternal age * paternal age - x=0.2073  df=1  p=0.6489

lrs.model.m.c <- update(lrs.model.m.4, ~ . + rescale(AgeMother)*rescale(AgeFather))

summary(lrs.model.m.c)

anova(lrs.model.m.4, lrs.model.m.c)


#~~ maternal age * dominant male age - x=0.0109  df=1  p=0.917

lrs.model.m.d <- update(lrs.model.m.4, ~ . + rescale(AgeMother)*rescale(DominantMaleAge))

summary(lrs.model.m.d)

anova(lrs.model.m.4, lrs.model.m.d)


# # -----------------------------------------------------------------------

## LRS models: male offspring

## Within-subject centering

# # -----------------------------------------------------------------------

#~~ set up variables

#~~ z-transform covariates

LRSmdata$cAgeMother <- rescale(LRSmdata$AgeMother)
LRSmdata$cAgeFather <- rescale(LRSmdata$AgeFather)
LRSmdata$cDominantMaleAge <- rescale(LRSmdata$DominantMaleAge)
LRSmdata$cBirthYear <- rescale(LRSmdata$BirthYear)
LRSmdata$clogTQ <- rescale(LRSmdata$logTQ)
LRSmdata$cGroupSize <- rescale(LRSmdata$GroupSize)


#~~ make mean maternal, paternal and dominant male age variables and deviation of each value from the mean

LRSmdata<-LRSmdata %>% group_by(MotherID) %>% mutate(cMeanAgeMother = mean(cAgeMother)) %>% ungroup()
LRSmdata$cDevMeanAgeMother<-LRSmdata$cAgeMother-LRSmdata$cMeanAgeMother

LRSmdata<-LRSmdata %>% group_by(FatherID) %>% mutate(cMeanAgeFather = mean(cAgeFather)) %>% ungroup()
LRSmdata$cDevMeanAgeFather<-LRSmdata$cAgeFather-LRSmdata$cMeanAgeFather

LRSmdata<-LRSmdata %>% group_by(DominantMaleID) %>% mutate(cMeanAgeSocMale = mean(cDominantMaleAge)) %>% ungroup()
LRSmdata$cDevMeanAgeSocMale<-LRSmdata$cDominantMaleAge-LRSmdata$cMeanAgeSocMale

summary(LRSmdata$cMeanAgeMother)
summary(LRSmdata$cMeanAgeSocMale)
summary(LRSmdata$cMeanAgeFather)


#~~ are within vs between parental ages significant?

lrs.model.m.wb <- glmmTMB(LRS1YEAR~cMeanAgeMother + cDevMeanAgeMother +
                            cMeanAgeFather + cDevMeanAgeFather + 
                            cMeanAgeSocMale + cDevMeanAgeSocMale +
                            cBirthYear + clogTQ + cGroupSize + HelpersF + Siblings +
                     (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                     data=LRSmdata,
                     ziformula=~1,
                     family=poisson)

summary(lrs.model.m.wb)


#~~ save model output

LRSmmodelwb<-as.data.frame((summary(lrs.model.m.wb)$coefficients$cond))
LRSmmodelwb<-setDT(LRSmmodelwb, keep.rownames = TRUE)[]

write.csv(LRSmmodelwb, file="Models/LRS_male_model_wnbn.csv")


#~~ significance of fixed effects

#~~ between maternal age effect - x=6.4238  df=1  p=0.01126

lrs.model.m.wb.1 <- update(lrs.model.m.wb, ~ . -cMeanAgeMother)

anova(lrs.model.m.wb, lrs.model.m.wb.1)


#~~ within maternal age effect - x=0.4023  df=1  p=0.5259

lrs.model.m.wb.2 <- update(lrs.model.m.wb, ~ . -cDevMeanAgeMother)

anova(lrs.model.m.wb, lrs.model.m.wb.2)


#~~ between paternal age effect - x=4.3834  df=1  p=0.03629 

lrs.model.m.wb.3 <- update(lrs.model.m.wb, ~ . -cMeanAgeFather)

anova(lrs.model.m.wb, lrs.model.m.wb.3)


#~~ within paternal age effect - x=0.067  df=1  p=0.7958

lrs.model.m.wb.4 <- update(lrs.model.m.wb, ~ . -cDevMeanAgeFather)

anova(lrs.model.m.wb, lrs.model.m.wb.4)


#~~ between dominant male age effect - x=0.4908  df=1  p=0.4836

lrs.model.m.wb.5 <- update(lrs.model.m.wb, ~ . -cMeanAgeSocMale)

anova(lrs.model.m.wb, lrs.model.m.wb.5)


#~~ within dominant male age effect - x=1.1294  df=1  p=0.2879

lrs.model.m.wb.6 <- update(lrs.model.m.wb, ~ . -cDevMeanAgeSocMale)

anova(lrs.model.m.wb, lrs.model.m.wb.6)


#~~ birth year - x=16.84  df=1  p=4.067e-05

lrs.model.m.wb.7 <- update(lrs.model.m.wb, ~ . -cBirthYear)

anova(lrs.model.m.wb, lrs.model.m.wb.7)


#~~ territory quality - x=1.2105  df=1  p=0.2712

lrs.model.m.wb.8 <- update(lrs.model.m.wb, ~ . -clogTQ)

anova(lrs.model.m.wb, lrs.model.m.wb.8)


#~~ group size - x=0.3044  df=1  p=0.5811

lrs.model.m.wb.9 <- update(lrs.model.m.wb, ~ . -cGroupSize)

anova(lrs.model.m.wb, lrs.model.m.wb.9)


#~~ helpers - x=1.6427  df=1  p=0.2

lrs.model.m.wb.10 <- update(lrs.model.m.wb, ~ . -HelpersF)

anova(lrs.model.m.wb, lrs.model.m.wb.10)


#~~ siblings - x=0.1175  df=1  p=0.7318

lrs.model.m.wb.11 <- update(lrs.model.m.wb, ~ . -Siblings)

anova(lrs.model.m.wb, lrs.model.m.wb.11)


#~~ are these within vs between slopes significantly different?

lrs.model.m.wb.b <- glmmTMB(LRS1YEAR~cMeanAgeMother + cAgeMother +
                              cMeanAgeFather + cAgeFather +
                              cMeanAgeSocMale + cDominantMaleAge +
                              cBirthYear + clogTQ + cGroupSize + HelpersF + Siblings +
                            (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                          data=LRSmdata,
                          ziformula=~1,
                          family=poisson)

summary(lrs.model.m.wb.b)


#~~ save model output

LRSmmodelwb.b<-as.data.frame((summary(lrs.model.m.wb.b)$coefficients$cond))
LRSmmodelwb.b<-setDT(LRSmmodelwb.b, keep.rownames = TRUE)[]

write.csv(LRSmmodelwb.b, file="Models/LRS_male_model_wnbn_b.csv")


#~~ significance of fixed effects

#~~ difference between within and between maternal age slopes - x=1.2186  df=1  p=0.2696

lrs.model.m.wb.b2 <- update(lrs.model.m.wb.b, ~ . -cMeanAgeMother)

anova(lrs.model.m.wb.b, lrs.model.m.wb.b2)


#~~ within maternal age effect - x=0.4023  df=1  p=0.5259

lrs.model.m.wb.b3 <- update(lrs.model.m.wb.b, ~ . -cAgeMother)

anova(lrs.model.m.wb.b, lrs.model.m.wb.b3)


#~~ difference between within and between paternal age slopes - x=2.4487  df=1  p=0.1176

lrs.model.m.wb.b4 <- update(lrs.model.m.wb.b, ~ . -cMeanAgeFather)

anova(lrs.model.m.wb.b, lrs.model.m.wb.b4)


#~~ within paternal age effect - x=0.067  df=1  p=0.7958

lrs.model.m.wb.b5 <- update(lrs.model.m.wb.b, ~ . -cAgeFather)

anova(lrs.model.m.wb.b, lrs.model.m.wb.b5)


#~~ difference between within and between dominant male age slopes - x=1.4927  df=1  p=0.2218

lrs.model.m.wb.b6 <- update(lrs.model.m.wb.b, ~ . -cMeanAgeSocMale)

anova(lrs.model.m.wb.b, lrs.model.m.wb.b6)


#~~ within dominant male age effect - x=1.1294  df=1  p=0.2879

lrs.model.m.wb.b7 <- update(lrs.model.m.wb.b, ~ . -cDominantMaleAge)

anova(lrs.model.m.wb.b, lrs.model.m.wb.b7)


#~~ birth year - x=16.84  df=1  p=4.067e-05

lrs.model.m.wb.b8 <- update(lrs.model.m.wb.b, ~ . -cBirthYear)

anova(lrs.model.m.wb.b, lrs.model.m.wb.b8)


#~~ territory quality - x=1.2105  df=1  p=0.2712

lrs.model.m.wb.b9 <- update(lrs.model.m.wb.b, ~ . -clogTQ)

anova(lrs.model.m.wb.b, lrs.model.m.wb.b9)


#~~ group size - x=0.3044  df=1  p=0.5811

lrs.model.m.wb.b10 <- update(lrs.model.m.wb.b, ~ . -cGroupSize)

anova(lrs.model.m.wb.b, lrs.model.m.wb.b10)


#~~ helpers - x=1.6427  df=1  p=0.2

lrs.model.m.wb.b11 <- update(lrs.model.m.wb.b, ~ . -HelpersF)

anova(lrs.model.m.wb.b, lrs.model.m.wb.b11)


#~~ siblings - x=0.1175  df=1  p=0.7318

lrs.model.m.wb.b12 <- update(lrs.model.m.wb.b, ~ . -Siblings)

anova(lrs.model.m.wb.b, lrs.model.m.wb.b12)


#~~ sanity check that maternal age estimates are similar if age of the father and age of the dominant male are not separated into a within and between effect

lrs.model.m.wb.b1 <- glmmTMB(LRS1YEAR~cMeanAgeMother + cDevMeanAgeMother +
                               cAgeFather + cDominantMaleAge +
                               cBirthYear + clogTQ + cGroupSize + HelpersF + Siblings +
                               (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                             data=LRSmdata,
                             ziformula=~1,
                             family=poisson)

summary(lrs.model.m.wb.b1)

lrs.model.m.wb.b2 <- glmmTMB(LRS1YEAR~cMeanAgeMother + cAgeMother +
                               cAgeFather + cDominantMaleAge +
                               cBirthYear + clogTQ + cGroupSize + HelpersF + Siblings +
                               (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                             data=LRSmdata,
                             ziformula=~1,
                             family=poisson)

summary(lrs.model.m.wb.b2)



