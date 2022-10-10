# # ------------------------------------------------------------------------

# ---- Parental age manuscript 
# ---- AMS
# ---- Lifespan models

# # ------------------------------------------------------------------------

rm(list=ls())


#~~ libraries

library(coxme)
library(arm)
library(survival)
library(glmmTMB)
library(DHARMa)
library(beepr)
library(performance)
library(data.table)
library(survminer)
library(dplyr)


#~~ data

fitnessdata <- read.csv("Data/Sparks_Parental_Age_Lifespan_Data_2022.csv", stringsAsFactors = FALSE)

#~~ NOTE: the dataset is already cropped to individuals who made it to at least 1 year


#~~ set factors

fitnessdata$MotherID<-as.factor(fitnessdata$MotherID)
fitnessdata$FatherID<-as.factor(fitnessdata$FatherID)
fitnessdata$DominantMaleID<-as.factor(fitnessdata$DominantMaleID)
fitnessdata$ID<-as.factor(fitnessdata$ID)
fitnessdata$BirthYearF<-as.factor(fitnessdata$BirthYear)
fitnessdata$Sex<-as.factor(fitnessdata$Sex)
fitnessdata$HelpersF <- as.factor(fitnessdata$HelpersF)
fitnessdata$TranslocatedF  <-  as.factor(fitnessdata$TranslocatedF)
fitnessdata$logTQ  <-  log10(fitnessdata$TQ)
fitnessdata$AgeMotherSq  <-  fitnessdata$AgeMother * fitnessdata$AgeMother
fitnessdata$AgeFatherSq  <-  fitnessdata$AgeFather * fitnessdata$AgeFather
fitnessdata$AgeSocFatherSq  <-  fitnessdata$DominantMaleAge * fitnessdata$DominantMaleAge
fitnessdata$Siblings  <-  as.factor(fitnessdata$Siblings)
fitnessdata$LifespanInt  <-  fitnessdata$FinalYear - fitnessdata$BirthYear


#~~ VIF function

##Variance inflation factors
vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}


# # -----------------------------------------------------------------------

## Lifespan GLMMs - females

# # -----------------------------------------------------------------------

#~~ lifespan GLMM data subset

# These GLMMs are only run on individuals who survived to at least one year, 
# and who were already dead (final year before 2019 and individuals not translocated) 

#~~ female model subset

lifespan.data.dead.f <- droplevels(subset(fitnessdata, LifespanInt > 0 & Sex==0 & FinalYear < 2019 & TranslocatedF==0
                                          & !is.na(GroupSize) & !is.na(AgeMother) & !is.na(AgeFather) & !is.na(DominantMaleAge) 
                                          & !is.na(BirthYear) & !is.na(Siblings) & !is.na(HelpersF)))


#~~ check all individuals dead

table(lifespan.data.dead.f$Event)


#~~ data

table(lifespan.data.dead.f$LifespanInt) 
hist(lifespan.data.dead.f$LifespanInt)
summary(lifespan.data.dead.f$LifespanInt)


#~~ model - run a variety of GLMMs in glmmTMB and see which model fits best

lifespan.f.poisson <- glmmTMB(LifespanInt~ rescale(AgeMother) + rescale(AgeFather) + rescale(DominantMaleAge) +
                             rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                             Siblings + (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                           data=lifespan.data.dead.f,
                           ziformula=~0,
                           family=poisson)
summary(lifespan.f.poisson)

lifespan.f.binom2 <- update(lifespan.f.poisson, family=nbinom2)

lifespan.f.binom1 <- update(lifespan.f.poisson, family=nbinom1)

# all models converged


#~~ AIC comparisons

AIC(lifespan.f.poisson, lifespan.f.binom2, lifespan.f.binom1)


#~~ check model residuals
# best fit seems to be poisson model
# slight deviation appears to be caused by birth year - but cannot fit birth year as a quadratic in the model 

lifespan.f.poisson.resid <- simulateResiduals(fittedModel = lifespan.f.poisson, plot = T)

plotResiduals(lifespan.f.poisson.resid, lifespan.data.dead.f$BirthYear)

lifespan.f.binom2.resid <- simulateResiduals(fittedModel = lifespan.f.binom2, plot = T)

lifespan.f.binom1.resid <- simulateResiduals(fittedModel = lifespan.f.binom1, plot = T)


#~~ go with poisson model
# full model with linear and quadratic mother and father ages

lifespan.f.glmm.model <- glmmTMB(LifespanInt~ rescale(AgeMother) + rescale(AgeFather) + rescale(DominantMaleAge) +
                         rescale(AgeMotherSq) + rescale(AgeFatherSq) + rescale(AgeSocFatherSq) +
                         rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                         Siblings + (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                       data=lifespan.data.dead.f,
                       ziformula=~0,
                       family=poisson)

summary(lifespan.f.glmm.model)


#~~ drop age mother sq

lifespan.f.glmm.model.2 <- update(lifespan.f.glmm.model, ~ . -rescale(AgeMotherSq))

anova(lifespan.f.glmm.model, lifespan.f.glmm.model.2)

summary(lifespan.f.glmm.model.2)


#~~ drop age father sq

lifespan.f.glmm.model.3 <- update(lifespan.f.glmm.model.2, ~ . -rescale(AgeFatherSq))

anova(lifespan.f.glmm.model.3, lifespan.f.glmm.model.2)

summary(lifespan.f.glmm.model.3)


#~~ drop dominant male age sq

lifespan.f.glmm.model.4 <- update(lifespan.f.glmm.model.3, ~ . -rescale(AgeSocFatherSq))

anova(lifespan.f.glmm.model.3, lifespan.f.glmm.model.4)


#~~ base model

summary(lifespan.f.glmm.model.4)

check_collinearity(lifespan.f.glmm.model.4)
check_overdispersion(lifespan.f.glmm.model.4)


#~~ save output

lifespanfmodelfixef<-as.data.frame((summary(lifespan.f.glmm.model.4)$coefficients$cond))
lifespanfmodelfixef<-setDT(lifespanfmodelfixef, keep.rownames = TRUE)[]

write.csv(lifespanfmodel, file="Models/Lifespan_female_model.csv")


#~~ significance of fixed effects

#~~ age of mum - x=10.569  df=1  p=0.00115

lifespan.f.glmm.model.5 <- update(lifespan.f.glmm.model.4, ~ . -rescale(AgeMother))

anova(lifespan.f.glmm.model.5, lifespan.f.glmm.model.4)


#~~ age of dad - x=0.2172  df=1  p=0.6412

lifespan.f.glmm.model.6 <- update(lifespan.f.glmm.model.4, ~ . -rescale(AgeFather))

anova(lifespan.f.glmm.model.6, lifespan.f.glmm.model.4)


#~~ age of dominant male - x=2.6763  df=1  p=0.1019

lifespan.f.glmm.model.7 <- update(lifespan.f.glmm.model.4, ~ . -rescale(DominantMaleAge))

anova(lifespan.f.glmm.model.7, lifespan.f.glmm.model.4)


#~~ birth year - x=22.852  df=1  p=1.749e-06

lifespan.f.glmm.model.8 <- update(lifespan.f.glmm.model.4, ~ . -rescale(BirthYear))

anova(lifespan.f.glmm.model.8, lifespan.f.glmm.model.4)


#~~ TQ - x=0.0542  df=1  p=0.816

lifespan.f.glmm.model.9 <- update(lifespan.f.glmm.model.4, ~ . -rescale(logTQ))

anova(lifespan.f.glmm.model.9, lifespan.f.glmm.model.4)


#~~ group size - x=1.9938  df=1  p=0.1579

lifespan.f.glmm.model.10 <- update(lifespan.f.glmm.model.4, ~ . -rescale(GroupSize))

anova(lifespan.f.glmm.model.10, lifespan.f.glmm.model.4)


#~~ helpers - x=2.7867  df=1  p=0.09505 

lifespan.f.glmm.model.11 <- update(lifespan.f.glmm.model.4, ~ . -HelpersF)

anova(lifespan.f.glmm.model.11, lifespan.f.glmm.model.4)


#~~ siblings - x=0.4389  df=1  p=0.5076

lifespan.f.glmm.model.12 <- update(lifespan.f.glmm.model.4, ~ . -Siblings)

anova(lifespan.f.glmm.model.12, lifespan.f.glmm.model.4)


#~~ age mother sq - x=0.0014  df=1  p=0.9702

lifespan.f.glmm.model.13 <- update(lifespan.f.glmm.model.4, ~ . +rescale(AgeMotherSq))

anova(lifespan.f.glmm.model.4, lifespan.f.glmm.model.13)


#~~ age father sq - x=0.4829  df=1  p=0.4871

lifespan.f.glmm.model.14 <- update(lifespan.f.glmm.model.4, ~ . +rescale(AgeFatherSq))

anova(lifespan.f.glmm.model.4, lifespan.f.glmm.model.14)


#~~ dominant male age sq - x=1.5468  df=1  p=0.2136

lifespan.f.glmm.model.15 <- update(lifespan.f.glmm.model.4, ~ . +rescale(AgeSocFatherSq))

anova(lifespan.f.glmm.model.4, lifespan.f.glmm.model.15)


#~~ significance of interactions

#~~ TQ * mum age - x=0.4069  df=1  p=0.5235

lifespan.f.glmm.model.16 <- update(lifespan.f.glmm.model.4, ~ . +rescale(logTQ)*rescale(AgeMother))

anova(lifespan.f.glmm.model.4, lifespan.f.glmm.model.16)


#~~ TQ * dad age - x=0.052  df=1  p=0.8196

lifespan.f.glmm.model.17 <- update(lifespan.f.glmm.model.4, ~ . +rescale(logTQ)*rescale(AgeFather))

anova(lifespan.f.glmm.model.4, lifespan.f.glmm.model.17)


#~~ TQ * social male age - x=0.121  df=1  p=0.728

lifespan.f.glmm.model.18 <- update(lifespan.f.glmm.model.4, ~ . +rescale(logTQ)*rescale(DominantMaleAge))

anova(lifespan.f.glmm.model.4, lifespan.f.glmm.model.18)


#~~ group size * mum age - x=2.6499  df=1  p=0.1036

lifespan.f.glmm.model.19 <- update(lifespan.f.glmm.model.4, ~ . +rescale(GroupSize)*rescale(AgeMother))

anova(lifespan.f.glmm.model.4, lifespan.f.glmm.model.19)


#~~ group size * dad age - x=0.1264  df=1  p=0.7222

lifespan.f.glmm.model.20 <- update(lifespan.f.glmm.model.4, ~ . +rescale(GroupSize)*rescale(AgeFather))

anova(lifespan.f.glmm.model.4, lifespan.f.glmm.model.20)


#~~ group size * social male age - x=1.407  df=1  p=0.2356

lifespan.f.glmm.model.21 <- update(lifespan.f.glmm.model.4, ~ . +rescale(GroupSize)*rescale(DominantMaleAge))

anova(lifespan.f.glmm.model.4, lifespan.f.glmm.model.21)


#~~ helpers * mum age - x=0.0016  df=1  p=0.9681

lifespan.f.glmm.model.22 <- update(lifespan.f.glmm.model.4, ~ . +HelpersF*rescale(AgeMother))

anova(lifespan.f.glmm.model.4, lifespan.f.glmm.model.22)


#~~ helpers * dad age - x=2.1633  df=1  p=0.1413

lifespan.f.glmm.model.23 <- update(lifespan.f.glmm.model.4, ~ . +HelpersF*rescale(AgeFather))

anova(lifespan.f.glmm.model.4, lifespan.f.glmm.model.23)


#~~ helpers * social male age - x=0.3257  df=1  p=0.5682

lifespan.f.glmm.model.24 <- update(lifespan.f.glmm.model.4, ~ . +HelpersF*rescale(DominantMaleAge))

anova(lifespan.f.glmm.model.4, lifespan.f.glmm.model.24)


#~~ check results are consistent if age of father and age of dominant male included separately 

#~~ genetic parents

lifespan.f.glmm.model.a <- glmmTMB(LifespanInt~ rescale(AgeMother) + rescale(AgeFather) +
                                     rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                                     Siblings + (1|MotherID) + (1|FatherID) + (1|BirthYearF),
                                   data=lifespan.data.dead.f,
                                   ziformula=~0,
                                   family=poisson)

summary(lifespan.f.glmm.model.a)


#~~ significance of maternal age - x=10.325  df=1  p=0.001312

lifespan.f.glmm.model.a1 <- update(lifespan.f.glmm.model.a, ~ . -rescale(AgeMother))

anova(lifespan.f.glmm.model.a1, lifespan.f.glmm.model.a)


#~~ significance of paternal age - x=0.1996  df=1  p=0.6551

lifespan.f.glmm.model.a2 <- update(lifespan.f.glmm.model.a, ~ . -rescale(AgeFather))

anova(lifespan.f.glmm.model.a2, lifespan.f.glmm.model.a)


#~~ model with genetic mother and dominant male

lifespan.f.glmm.model.b <- glmmTMB(LifespanInt~ rescale(AgeMother) + rescale(DominantMaleAge) +
                                     rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                                     Siblings + (1|MotherID) + (1|DominantMaleID) + (1|BirthYearF),
                                   data=lifespan.data.dead.f,
                                   ziformula=~0,
                                   family=poisson)

summary(lifespan.f.glmm.model.b)


#~~ significance of maternal age - x=14.342  df=1  p=0.0001524

lifespan.f.glmm.model.b1 <- update(lifespan.f.glmm.model.b, ~ . -rescale(AgeMother))

anova(lifespan.f.glmm.model.b1, lifespan.f.glmm.model.b)


#~~ significance of dominant male age - x=1.1464  df=1  p=0.2843

lifespan.f.glmm.model.b2 <- update(lifespan.f.glmm.model.b, ~ . -rescale(DominantMaleAge))

anova(lifespan.f.glmm.model.b2, lifespan.f.glmm.model.b)


#~~ check maternal age * paternal age interactions

#~~ maternal age * paternal age - x=0.1313  df=1  p=0.7171

lifespan.f.glmm.model.c <- glmmTMB(LifespanInt~ rescale(AgeMother) + rescale(AgeFather) + rescale(DominantMaleAge) +
                                     rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                                     Siblings + rescale(AgeMother)*rescale(AgeFather) + (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                                   data=lifespan.data.dead.f,
                                   ziformula=~0,
                                   family=poisson)

summary(lifespan.f.glmm.model.c)

anova(lifespan.f.glmm.model.4, lifespan.f.glmm.model.c)


# maternal age * dominant male age - x=2.7381  df=1  p=0.09798

lifespan.f.glmm.model.d <- glmmTMB(LifespanInt~ rescale(AgeMother) + rescale(AgeFather)  + rescale(DominantMaleAge) +
                                     rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                                     Siblings + rescale(AgeMother)*rescale(DominantMaleAge) + (1|MotherID) + (1|FatherID) + (1|DominantMaleID)  + (1|BirthYearF),
                                   data=lifespan.data.dead.f,
                                   ziformula=~0,
                                   family=poisson)

summary(lifespan.f.glmm.model.d)

anova(lifespan.f.glmm.model.4, lifespan.f.glmm.model.d)


# # -----------------------------------------------------------------------

## Lifespan GLMMs - females

## Within-subject centering

# # -----------------------------------------------------------------------

#~~ set up variables

#~~ z-transform covariates

lifespan.data.dead.f$cAgeMother <- rescale(lifespan.data.dead.f$AgeMother)
lifespan.data.dead.f$cAgeFather <- rescale(lifespan.data.dead.f$AgeFather)
lifespan.data.dead.f$cDominantMaleAge <- rescale(lifespan.data.dead.f$DominantMaleAge)
lifespan.data.dead.f$cBirthYear <- rescale(lifespan.data.dead.f$BirthYear)
lifespan.data.dead.f$clogTQ <- rescale(lifespan.data.dead.f$logTQ)
lifespan.data.dead.f$cGroupSize <- rescale(lifespan.data.dead.f$GroupSize)


#~~ make mean maternal, paternal and dominant male age variables and deviation of each value from the mean

lifespan.data.dead.f<-lifespan.data.dead.f %>% group_by(MotherID) %>% mutate(cMeanAgeMother = mean(cAgeMother)) %>% ungroup()
lifespan.data.dead.f$cDevMeanAgeMother<-lifespan.data.dead.f$cAgeMother-lifespan.data.dead.f$cMeanAgeMother

lifespan.data.dead.f<-lifespan.data.dead.f %>% group_by(FatherID) %>% mutate(cMeanAgeFather = mean(cAgeFather)) %>% ungroup()
lifespan.data.dead.f$cDevMeanAgeFather<-lifespan.data.dead.f$cAgeFather-lifespan.data.dead.f$cMeanAgeFather

lifespan.data.dead.f<-lifespan.data.dead.f %>% group_by(DominantMaleID) %>% mutate(cMeanAgeSocMale = mean(cDominantMaleAge)) %>% ungroup()
lifespan.data.dead.f$cDevMeanAgeSocMale<-lifespan.data.dead.f$cDominantMaleAge-lifespan.data.dead.f$cMeanAgeSocMale

lifespan.data.dead.f<-lifespan.data.dead.f %>% group_by(MotherID) %>% mutate(MeanAgeMother = mean(AgeMother)) %>% ungroup()
lifespan.data.dead.f$DevMeanAgeMother<-lifespan.data.dead.f$AgeMother-lifespan.data.dead.f$MeanAgeMother

summary(lifespan.data.dead.f$cMeanAgeMother)
summary(lifespan.data.dead.f$cMeanAgeSocMale)
summary(lifespan.data.dead.f$cMeanAgeFather)


#~~ are within vs between parental ages significant?

lifespan.model.f.wb <- glmmTMB(LifespanInt~cMeanAgeMother + cDevMeanAgeMother +
                            cMeanAgeFather + cDevMeanAgeFather + 
                            cMeanAgeSocMale + cDevMeanAgeSocMale +
                            cBirthYear + clogTQ + cGroupSize + HelpersF + Siblings +
                            (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                          data=lifespan.data.dead.f,
                          ziformula=~0,
                          family=poisson)

summary(lifespan.model.f.wb)


#~~ save model output

lifespanfmodelwb<-as.data.frame((summary(lifespan.model.f.wb)$coefficients$cond))
lifespanfmodelwb<-setDT(lifespanfmodelwb, keep.rownames = TRUE)[]

write.csv(lifespanfmodelwb, file="Models/Lifespan_female_model_wnbn.csv")


#~~ between maternal age effect - x=1.1397  df=1  p=0.2857

lifespan.model.f.wb.1 <- update(lifespan.model.f.wb, ~ . -cMeanAgeMother)

anova(lifespan.model.f.wb.1, lifespan.model.f.wb)


#~~ within maternal age effect - x=17.503  df=1  p=2.868e-05

lifespan.model.f.wb.2 <- update(lifespan.model.f.wb, ~ . -cDevMeanAgeMother)

anova(lifespan.model.f.wb.2, lifespan.model.f.wb)


#~~ between paternal age effect - x=0.0012  df=1  p=0.9726

lifespan.model.f.wb.3 <- update(lifespan.model.f.wb, ~ . -cMeanAgeFather)

anova(lifespan.model.f.wb.3, lifespan.model.f.wb)


#~~ within paternal age effect - x=1.3614  df=1  p=0.2433

lifespan.model.f.wb.4 <- update(lifespan.model.f.wb, ~ . -cDevMeanAgeFather)

anova(lifespan.model.f.wb.4, lifespan.model.f.wb)


#~~ between dominant male age effect - x=1.4372  df=1  p=0.2306

lifespan.model.f.wb.5 <- update(lifespan.model.f.wb, ~ . -cMeanAgeSocMale)

anova(lifespan.model.f.wb.5, lifespan.model.f.wb)


#~~ within dominant male age effect - x=2.9744  df=1  p=0.08459

lifespan.model.f.wb.6 <- update(lifespan.model.f.wb, ~ . -cDevMeanAgeSocMale)

anova(lifespan.model.f.wb.6, lifespan.model.f.wb)


#~~ birth year - x=19  df=1  p=1.307e-05

lifespan.model.f.wb.7 <- update(lifespan.model.f.wb, ~ . -cBirthYear)

anova(lifespan.model.f.wb.7, lifespan.model.f.wb)


#~~ territory quality - x=0.0753  df=1  p=0.7838

lifespan.model.f.wb.8 <- update(lifespan.model.f.wb, ~ . -clogTQ)

anova(lifespan.model.f.wb.8, lifespan.model.f.wb)


#~~ group size - x=1.8249  df=1  p=0.1767

lifespan.model.f.wb.9 <- update(lifespan.model.f.wb, ~ . -cGroupSize)

anova(lifespan.model.f.wb.9, lifespan.model.f.wb)


#~~ helpers - x=2.8057  df=1  p=0.09393

lifespan.model.f.wb.10 <- update(lifespan.model.f.wb, ~ . -HelpersF)

anova(lifespan.model.f.wb.10, lifespan.model.f.wb)


#~~ siblings - x=0.446  df=1  p=0.5042

lifespan.model.f.wb.11 <- update(lifespan.model.f.wb, ~ . -Siblings)

anova(lifespan.model.f.wb.11, lifespan.model.f.wb)


#~~ check model similar without outlier for females within maternal age graph

lifespanfdata.check <- droplevels(subset(lifespan.data.dead.f, DevMeanAgeMother > -6.9))

lifespan.model.f.wb.check <- glmmTMB(LifespanInt~cMeanAgeMother + cDevMeanAgeMother +
                                  cMeanAgeFather + cDevMeanAgeFather + 
                                  cMeanAgeSocMale + cDevMeanAgeSocMale +
                                  cBirthYear + clogTQ + cGroupSize + HelpersF + Siblings +
                                  (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                                data=lifespanfdata.check,
                                ziformula=~1,
                                family=poisson)

summary(lifespan.model.f.wb.check)


#~~ save model output

lifespanfmodelwbcheckfixef<-as.data.frame((summary(lifespan.model.f.wb.check)$coefficients$cond))
lifespanfmodelwbcheckfixef<-setDT(lifespanfmodelwbcheckfixef, keep.rownames = TRUE)[]

write.csv(lifespanfmodelwbcheckfixef, file="Models/Lifespan_female_model_wnbn_check.csv")


#~~ between maternal age effect - x=1.3093  df=1  p=0.2525

lifespan.model.f.wb.check1 <- update(lifespan.model.f.wb.check, ~ . -cMeanAgeMother)

anova(lifespan.model.f.wb.check1, lifespan.model.f.wb.check)


#~~ within maternal age effect - x=14.046  df=1  p=0.0001784

lifespan.model.f.wb.check2 <- update(lifespan.model.f.wb.check, ~ . -cDevMeanAgeMother)

anova(lifespan.model.f.wb.check2, lifespan.model.f.wb.check)


#~~ between paternal age effect - x=0.0111  df=1  p=0.9162

lifespan.model.f.wb.check3 <- update(lifespan.model.f.wb.check, ~ . -cMeanAgeFather)

anova(lifespan.model.f.wb.check3, lifespan.model.f.wb.check)


#~~ within paternal age effect - x=1.4603  df=1  p=0.2269

lifespan.model.f.wb.check4 <- update(lifespan.model.f.wb.check, ~ . -cDevMeanAgeFather)

anova(lifespan.model.f.wb.check4, lifespan.model.f.wb.check)


#~~ between dominant male age effect - x=1.3833  df=1  p=0.2395

lifespan.model.f.wb.check5 <- update(lifespan.model.f.wb.check, ~ . -cMeanAgeSocMale)

anova(lifespan.model.f.wb.check5, lifespan.model.f.wb.check)


#~~ within dominant male age effect - x=2.68  df=1  p=0.1016

lifespan.model.f.wb.check6 <- update(lifespan.model.f.wb.check, ~ . -cDevMeanAgeSocMale)

anova(lifespan.model.f.wb.check6, lifespan.model.f.wb.check)


#~~ birth year - x=18.46  df=1  p=1.735e-05

lifespan.model.f.wb.check7 <- update(lifespan.model.f.wb.check, ~ . -cBirthYear)

anova(lifespan.model.f.wb.check7, lifespan.model.f.wb.check)


#~~ territory quality - x=0.0669  df=1  p=0.7959

lifespan.model.f.wb.check8 <- update(lifespan.model.f.wb.check, ~ . -clogTQ)

anova(lifespan.model.f.wb.check8, lifespan.model.f.wb.check)


#~~ group size - x=1.7044  df=1  p=0.1917

lifespan.model.f.wb.check9 <- update(lifespan.model.f.wb.check, ~ . -cGroupSize)

anova(lifespan.model.f.wb.check9, lifespan.model.f.wb.check)


#~~ helpers - x=2.6777  df=1  p=0.1018

lifespan.model.f.wb.check10 <- update(lifespan.model.f.wb.check, ~ . -HelpersF)

anova(lifespan.model.f.wb.check10, lifespan.model.f.wb.check)


#~~ siblings - x=0.3753  df=1  p=0.5401

lifespan.model.f.wb.check11 <- update(lifespan.model.f.wb.check, ~ . -Siblings)

anova(lifespan.model.f.wb.check11, lifespan.model.f.wb.check)


#~~ are these within vs between slopes significantly different?

lifespan.model.f.wb.b <- glmmTMB(LifespanInt~cMeanAgeMother + cAgeMother +
                                   cMeanAgeFather + cAgeFather +
                                   cMeanAgeSocMale + cDominantMaleAge +
                                   cBirthYear + clogTQ + cGroupSize + HelpersF + Siblings +
                                 (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                               data=lifespan.data.dead.f,
                               ziformula=~0,
                               family=poisson)


summary(lifespan.model.f.wb.b)


#~~ save model output

lifespanfmodelwb.b.fixef<-as.data.frame((summary(lifespan.model.f.wb.b)$coefficients$cond))
lifespanfmodelwb.b.fixef<-setDT(lifespanfmodelwb.b.fixef, keep.rownames = TRUE)[]

write.csv(lifespanfmodelwb.b.fixef, file="Models/Lifespan_female_model_wnbn_b.csv")


#~~ significance of fixed effects

#~~ difference between within and between maternal age slopes - x=7.9391  df=1  p=0.004838 

lifespan.model.f.wb.b2 <- update(lifespan.model.f.wb.b, ~ . -cMeanAgeMother)

anova(lifespan.model.f.wb.b2, lifespan.model.f.wb.b)


#~~ within maternal age effect - x=17.503  df=1  p=2.868e-05

lifespan.model.f.wb.b3 <- update(lifespan.model.f.wb.b, ~ . -cAgeMother)

anova(lifespan.model.f.wb.b3, lifespan.model.f.wb.b)


#~~ difference between within and between paternal age slopes - x=0.9284  df=1  p=0.3353

lifespan.model.f.wb.b4 <- update(lifespan.model.f.wb.b, ~ . -cMeanAgeFather)

anova(lifespan.model.f.wb.b4, lifespan.model.f.wb.b)


#~~ within paternal age effect - x=1.3614  df=1  p=0.2433

lifespan.model.f.wb.b5 <- update(lifespan.model.f.wb.b, ~ . -cAgeFather)

anova(lifespan.model.f.wb.b5, lifespan.model.f.wb.b)


#~~ difference between within and between dominant male age slopes - x=0.6753  df=1  p=0.4112

lifespan.model.f.wb.b6 <- update(lifespan.model.f.wb.b, ~ . -cMeanAgeSocMale)

anova(lifespan.model.f.wb.b6, lifespan.model.f.wb.b)


#~~ within dominant male age effect - x=2.9744  df=1  p=0.08459

lifespan.model.f.wb.b7 <- update(lifespan.model.f.wb.b, ~ . -cDominantMaleAge)

anova(lifespan.model.f.wb.b7, lifespan.model.f.wb.b)


#~~ birth year - x=19  df=1  p=1.307e-05

lifespan.model.f.wb.b8 <- update(lifespan.model.f.wb.b, ~ . -cBirthYear)

anova(lifespan.model.f.wb.b8, lifespan.model.f.wb.b)


#~~ territory quality - x=0.0753  df=1  p=0.7838

lifespan.model.f.wb.b9 <- update(lifespan.model.f.wb.b, ~ . -clogTQ)

anova(lifespan.model.f.wb.b9, lifespan.model.f.wb.b)


#~~ group size - x=1.8249  df=1  p=0.1767

lifespan.model.f.wb.b10 <- update(lifespan.model.f.wb.b, ~ . -cGroupSize)

anova(lifespan.model.f.wb.b10, lifespan.model.f.wb.b)


#~~ helpers - x=2.8057  df=1  p=0.09393

lifespan.model.f.wb.b11 <- update(lifespan.model.f.wb.b, ~ . -HelpersF)

anova(lifespan.model.f.wb.b11, lifespan.model.f.wb.b)


#~~ siblings - x=0.446  df=1  p=0.5042

lifespan.model.f.wb.b12 <- update(lifespan.model.f.wb.b, ~ . -Siblings)

anova(lifespan.model.f.wb.b12, lifespan.model.f.wb.b)


# #  ----------------------------------------------------------------------

# Lifespan GLMMs - males

# # -----------------------------------------------------------------------

#~~ lifespan GLMM data subset

# These GLMMs are only run on individuals who survived to at least one year, 
# and who were already dead (final year before 2019 and individuals not translocated) 

#~~ male model subset

lifespan.data.dead.m <- droplevels(subset(fitnessdata, LifespanInt > 0 & Sex==1 & FinalYear < 2019 & TranslocatedF==0 
                                          & !is.na(GroupSize) & !is.na(AgeMother) & !is.na(AgeFather) & !is.na(DominantMaleAge) 
                                          & !is.na(BirthYear) & !is.na(Siblings) & !is.na(HelpersF)))

#~~ check all individuals dead

table(lifespan.data.dead.m$Event)


#~~ model - run a variety of GLMMs in glmmTMB and see which model fits best

lifespan.m.poisson <- glmmTMB(LifespanInt~ rescale(AgeMother) + rescale(AgeFather) + rescale(DominantMaleAge) +
                                rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                                Siblings + (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                              data=lifespan.data.dead.m,
                              ziformula=~0,
                              family=poisson)
summary(lifespan.m.poisson)

lifespan.m.binom2 <- update(lifespan.m.poisson, family=nbinom2)

lifespan.m.binom1 <- update(lifespan.m.poisson, family=nbinom1)

# all models converged


#~~ AIC comparisons

AIC(lifespan.m.poisson, lifespan.m.binom2, lifespan.m.binom1)


#~~ check model residuals

lifespan.m.poisson.resid <- simulateResiduals(fittedModel = lifespan.m.poisson, plot = T)

lifespan.m.binom2.resid <- simulateResiduals(fittedModel = lifespan.m.binom2, plot = T)

lifespan.m.binom1.resid <- simulateResiduals(fittedModel = lifespan.m.binom1, plot = T) # worst


#~~ go with poisson model
# full model with linear and quadratic mother and father ages

lifespan.m.glmm.model <- glmmTMB(LifespanInt~ rescale(AgeMother) + rescale(AgeFather) + rescale(DominantMaleAge) +
                                   rescale(AgeMotherSq) + rescale(AgeFatherSq) + rescale(AgeSocFatherSq) +
                                   rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                                   Siblings + (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                                 data=lifespan.data.dead.m,
                                 ziformula=~0,
                                 family=poisson)

summary(lifespan.m.glmm.model)


#~~ drop age father sq

lifespan.m.glmm.model.2 <- update(lifespan.m.glmm.model, ~ . -rescale(AgeFatherSq))

anova(lifespan.m.glmm.model, lifespan.m.glmm.model.2)

summary(lifespan.m.glmm.model.2)


#~~ drop dominant male age sq

lifespan.m.glmm.model.3 <- update(lifespan.m.glmm.model.2, ~ . -rescale(AgeSocFatherSq))

anova(lifespan.m.glmm.model.3, lifespan.m.glmm.model.2)

summary(lifespan.m.glmm.model.3)


#~~ drop age mother sq

lifespan.m.glmm.model.4 <- update(lifespan.m.glmm.model.3, ~ . -rescale(AgeMotherSq))

anova(lifespan.m.glmm.model.3, lifespan.m.glmm.model.4)


#~~ base model

summary(lifespan.m.glmm.model.4)

check_collinearity(lifespan.m.glmm.model.4)
check_overdispersion(lifespan.m.glmm.model.4)


#~~ save output

lifespanmmodel.fixef<-as.data.frame((summary(lifespan.m.glmm.model.4)$coefficients$cond))
lifespanmmodel.fixef<-setDT(lifespanmmodel.fixef, keep.rownames = TRUE)[]

write.csv(lifespanmmodel.fixef, file="Models/Lifespan_male_model.csv")


#~~ significance of fixed effects

#~~ age of mum - x=2.0983  df=1  p=0.1475

lifespan.m.glmm.model.5 <- update(lifespan.m.glmm.model.4, ~ . -rescale(AgeMother))

anova(lifespan.m.glmm.model.5, lifespan.m.glmm.model.4)


#~~ age of dad - x=0.545  df=1  p=0.4603

lifespan.m.glmm.model.6 <- update(lifespan.m.glmm.model.4, ~ . -rescale(AgeFather))

anova(lifespan.m.glmm.model.6, lifespan.m.glmm.model.4)


#~~ age of dominant male - x=2.3634  df=1  p=0.1242

lifespan.m.glmm.model.7 <- update(lifespan.m.glmm.model.4, ~ . -rescale(DominantMaleAge))

anova(lifespan.m.glmm.model.7, lifespan.m.glmm.model.4)


#~~ birth year - x=13.71  df=1  p=0.0002133

lifespan.m.glmm.model.8 <- update(lifespan.m.glmm.model.4, ~ . -rescale(BirthYear))

anova(lifespan.m.glmm.model.8, lifespan.m.glmm.model.4)


#~~ TQ - x=2.5977  df=1  p=0.107

lifespan.m.glmm.model.9 <- update(lifespan.m.glmm.model.4, ~ . -rescale(logTQ))

anova(lifespan.m.glmm.model.9, lifespan.m.glmm.model.4)


#~~ group size - x=0.8625  df=1  p=0.353

lifespan.m.glmm.model.10 <- update(lifespan.m.glmm.model.4, ~ . -rescale(GroupSize))

anova(lifespan.m.glmm.model.10, lifespan.m.glmm.model.4)


#~~ helpers - x=1.6313  df=1  p=0.2015

lifespan.m.glmm.model.11 <- update(lifespan.m.glmm.model.4, ~ . -HelpersF)

anova(lifespan.m.glmm.model.11, lifespan.m.glmm.model.4)


#~~ siblings - x=1.8346  df=1  p=0.1756

lifespan.m.glmm.model.12 <- update(lifespan.m.glmm.model.4, ~ . -Siblings)

anova(lifespan.m.glmm.model.12, lifespan.m.glmm.model.4)


#~~ age mother sq - x=2.5225  df=1  p=0.1122

lifespan.m.glmm.model.13 <- update(lifespan.m.glmm.model.4, ~ . +rescale(AgeMotherSq))

anova(lifespan.m.glmm.model.4, lifespan.m.glmm.model.13)


#~~ age father sq - x=0.8083  df=1  p=0.3686

lifespan.m.glmm.model.14 <- update(lifespan.m.glmm.model.4, ~ . +rescale(AgeFatherSq))

anova(lifespan.m.glmm.model.4, lifespan.m.glmm.model.14)


#~~ dominant male age sq - x=1.3471  df=1  p=0.2458

lifespan.m.glmm.model.15 <- update(lifespan.m.glmm.model.4, ~ . +rescale(AgeSocFatherSq))

anova(lifespan.m.glmm.model.4, lifespan.m.glmm.model.15)


#~~ significance of interactions

#~~ TQ* mum age - x=0.2903  df=1  p=0.5901

lifespan.m.glmm.model.16 <- update(lifespan.m.glmm.model.4, ~ . +rescale(logTQ)*rescale(AgeMother))

anova(lifespan.m.glmm.model.4, lifespan.m.glmm.model.16)


#~~ TQ*dad age - x=0.2954  df=1  p=0.5868

lifespan.m.glmm.model.17 <- update(lifespan.m.glmm.model.4, ~ . +rescale(logTQ)*rescale(AgeFather))

anova(lifespan.m.glmm.model.4, lifespan.m.glmm.model.17)


#~~ TQ* social male age - x=0.0948  df=1  p=0.7582

lifespan.m.glmm.model.18 <- update(lifespan.m.glmm.model.4, ~ . +rescale(logTQ)*rescale(DominantMaleAge))

anova(lifespan.m.glmm.model.4, lifespan.m.glmm.model.18)


#~~ group size * mum age - x=1.1644  df=1  p=0.2806

lifespan.m.glmm.model.19 <- update(lifespan.m.glmm.model.4, ~ . +rescale(GroupSize)*rescale(AgeMother))

anova(lifespan.m.glmm.model.4, lifespan.m.glmm.model.19)


#~~ group size * dad age - x=2.082  df=1  p=0.149

lifespan.m.glmm.model.20 <- update(lifespan.m.glmm.model.4, ~ . +rescale(GroupSize)*rescale(AgeFather))

anova(lifespan.m.glmm.model.4, lifespan.m.glmm.model.20)


#~~ group size * social male age - x=1.5628  df=1  p=0.2112

lifespan.m.glmm.model.21 <- update(lifespan.m.glmm.model.4, ~ . +rescale(GroupSize)*rescale(DominantMaleAge))

anova(lifespan.m.glmm.model.4, lifespan.m.glmm.model.21)


#~~ helpers*mum age - x=0.0565  df=1  p=0.8122

lifespan.m.glmm.model.22 <- update(lifespan.m.glmm.model.4, ~ . +HelpersF*rescale(AgeMother))

anova(lifespan.m.glmm.model.4, lifespan.m.glmm.model.22)


#~~ helpers*dad age - x=0.8383  df=1  p=0.3599

lifespan.m.glmm.model.23 <- update(lifespan.m.glmm.model.4, ~ . +HelpersF*rescale(AgeFather))

anova(lifespan.m.glmm.model.4, lifespan.m.glmm.model.23)


#~~ helpers*social male age - x=1.1215  df=1  p=0.2896

lifespan.m.glmm.model.24 <- update(lifespan.m.glmm.model.4, ~ . +HelpersF*rescale(DominantMaleAge))

anova(lifespan.m.glmm.model.4, lifespan.m.glmm.model.24)


#~~ check results are consistent if age of father and age of dominant male included separately 

#~~ genetic parents

lifespan.m.glmm.model.a <- glmmTMB(LifespanInt~ rescale(AgeMother) + rescale(AgeFather) +
                                     rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                                     Siblings + (1|MotherID) + (1|FatherID) + (1|BirthYearF),
                                   data=lifespan.data.dead.m,
                                   ziformula=~0,
                                   family=poisson)

summary(lifespan.m.glmm.model.a)


#~~ significance of maternal age - x=2.3337  df=1  p=0.1266

lifespan.m.glmm.model.a1 <- update(lifespan.m.glmm.model.a, ~ . -rescale(AgeMother))

anova(lifespan.m.glmm.model.a1, lifespan.m.glmm.model.a)


# significance of paternal age - x=0.0318  df=1  p=0.8584

lifespan.m.glmm.model.a2 <- update(lifespan.m.glmm.model.a, ~ . -rescale(AgeFather))

anova(lifespan.m.glmm.model.a2, lifespan.m.glmm.model.a)


#~~ model with genetic mother and dominant male

lifespan.m.glmm.model.b <- glmmTMB(LifespanInt~ rescale(AgeMother) + rescale(DominantMaleAge) +
                                     rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                                     Siblings + (1|MotherID) + (1|DominantMaleID) + (1|BirthYearF),
                                   data=lifespan.data.dead.m,
                                   ziformula=~0,
                                   family=poisson)


summary(lifespan.m.glmm.model.b)


#~~ significance of maternal age - x=1.7605  df=1  p=0.1846

lifespan.m.glmm.model.b1 <- update(lifespan.m.glmm.model.b, ~ . -rescale(AgeMother))

anova(lifespan.m.glmm.model.b1, lifespan.m.glmm.model.b)


#~~ significance of dominant male age - x=2.4305  df=1  p=0.119

lifespan.m.glmm.model.b2 <- update(lifespan.m.glmm.model.b, ~ . -rescale(DominantMaleAge))

anova(lifespan.m.glmm.model.b2, lifespan.m.glmm.model.b)


#~~ check maternal age by paternal age interactions

#~~ maternal age * paternal age - x=1.1442  df=1  p=0.2848

lifespan.m.glmm.model.c <- glmmTMB(LifespanInt~ rescale(AgeMother) + rescale(AgeFather) + rescale(DominantMaleAge) +
                                     rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                                     Siblings + rescale(AgeMother)*rescale(AgeFather) + (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                                   data=lifespan.data.dead.m,
                                   ziformula=~0,
                                   family=poisson)

summary(lifespan.m.glmm.model.c)


anova(lifespan.m.glmm.model.4, lifespan.m.glmm.model.c)


#~~ maternal age * dominant male age - x=0.0995  df=1  p=0.7525

lifespan.m.glmm.model.d <- glmmTMB(LifespanInt~ rescale(AgeMother) + rescale(AgeFather)  + rescale(DominantMaleAge) +
                                     rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                                     Siblings + rescale(AgeMother)*rescale(DominantMaleAge) + (1|MotherID) + (1|FatherID) + (1|DominantMaleID)  + (1|BirthYearF),
                                   data=lifespan.data.dead.m,
                                   ziformula=~0,
                                   family=poisson)

summary(lifespan.m.glmm.model.d)

anova(lifespan.m.glmm.model.4, lifespan.m.glmm.model.d)


# # -----------------------------------------------------------------------

## Lifespan GLMMs - males

## Within-subject centering

# # -----------------------------------------------------------------------

#~~ set up variables

#~~ z-transform covariates

lifespan.data.dead.m$cAgeMother <- rescale(lifespan.data.dead.m$AgeMother)
lifespan.data.dead.m$cAgeFather <- rescale(lifespan.data.dead.m$AgeFather)
lifespan.data.dead.m$cDominantMaleAge <- rescale(lifespan.data.dead.m$DominantMaleAge)
lifespan.data.dead.m$cBirthYear <- rescale(lifespan.data.dead.m$BirthYear)
lifespan.data.dead.m$clogTQ <- rescale(lifespan.data.dead.m$logTQ)
lifespan.data.dead.m$cGroupSize <- rescale(lifespan.data.dead.m$GroupSize)


#~~ make mean maternal, paternal and dominant male age variables and deviation of each value from the mean

lifespan.data.dead.m<-lifespan.data.dead.m %>% group_by(MotherID) %>% mutate(cMeanAgeMother = mean(cAgeMother)) %>% ungroup()
lifespan.data.dead.m$cDevMeanAgeMother<-lifespan.data.dead.m$cAgeMother-lifespan.data.dead.m$cMeanAgeMother

lifespan.data.dead.m<-lifespan.data.dead.m %>% group_by(FatherID) %>% mutate(cMeanAgeFather = mean(cAgeFather)) %>% ungroup()
lifespan.data.dead.m$cDevMeanAgeFather<-lifespan.data.dead.m$cAgeFather-lifespan.data.dead.m$cMeanAgeFather

lifespan.data.dead.m<-lifespan.data.dead.m %>% group_by(DominantMaleID) %>% mutate(cMeanAgeSocMale = mean(cDominantMaleAge)) %>% ungroup()
lifespan.data.dead.m$cDevMeanAgeSocMale<-lifespan.data.dead.m$cDominantMaleAge-lifespan.data.dead.m$cMeanAgeSocMale

summary(lifespan.data.dead.m$cMeanAgeMother)
summary(lifespan.data.dead.m$cMeanAgeSocMale)
summary(lifespan.data.dead.m$cMeanAgeFather)


#~~ are within vs between parental ages significant?

lifespan.model.m.wb <- glmmTMB(LifespanInt~cMeanAgeMother + cDevMeanAgeMother +
                                 cMeanAgeFather + cDevMeanAgeFather + 
                                 cMeanAgeSocMale + cDevMeanAgeSocMale +
                                 cBirthYear + clogTQ + cGroupSize + HelpersF + Siblings +
                                 (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                               data=lifespan.data.dead.m,
                               ziformula=~0,
                               family=poisson)

summary(lifespan.model.m.wb)


#~~ save model output

lifespanmmodelwb.fixef<-as.data.frame((summary(lifespan.model.m.wb)$coefficients$cond))
lifespanmmodelwb.fixef<-setDT(lifespanmmodelwb.fixef, keep.rownames = TRUE)[]

write.csv(lifespanmmodelwb.fixef, file="Models/Lifespan_male_model_wnbn.csv")


#~~ significance of fixed effects

#~~ between maternal age effect - x=1.7926  df=1  p=0.1806

lifespan.model.m.wb.1 <- update(lifespan.model.m.wb, ~ . -cMeanAgeMother)

anova(lifespan.model.m.wb.1, lifespan.model.m.wb)


#~~ within maternal age effect - x=0.9199  df=1  p=0.3375

lifespan.model.m.wb.2 <- update(lifespan.model.m.wb, ~ . -cDevMeanAgeMother)

anova(lifespan.model.m.wb.2, lifespan.model.m.wb)


#~~ between paternal age effect - x=0.189  df=1  p=0.6637

lifespan.model.m.wb.3 <- update(lifespan.model.m.wb, ~ . -cMeanAgeFather)

anova(lifespan.model.m.wb.3, lifespan.model.m.wb)


#~~ within paternal age effect - x=2.437  df=1  p=0.1185

lifespan.model.m.wb.4 <- update(lifespan.model.m.wb, ~ . -cDevMeanAgeFather)

anova(lifespan.model.m.wb.4, lifespan.model.m.wb)


#~~ between dominant male age effect - x=2.5701  df=1  p=0.1089

lifespan.model.m.wb.5 <- update(lifespan.model.m.wb, ~ . -cMeanAgeSocMale)

anova(lifespan.model.m.wb.5, lifespan.model.m.wb)


#~~ within dominant male age effect - x=0.6038  df=1  p=0.4371

lifespan.model.m.wb.6 <- update(lifespan.model.m.wb, ~ . -cDevMeanAgeSocMale)

anova(lifespan.model.m.wb.6, lifespan.model.m.wb)


#~~ birth year - x=12.252  df=1  p=0.0004647

lifespan.model.m.wb.7 <- update(lifespan.model.m.wb, ~ . -cBirthYear)

anova(lifespan.model.m.wb.7, lifespan.model.m.wb)


#~~ territory quality - x=2.6452  df=1  p=0.1039

lifespan.model.m.wb.8 <- update(lifespan.model.m.wb, ~ . -clogTQ)

anova(lifespan.model.m.wb.8, lifespan.model.m.wb)


#~~ group size - x=0.6543  df=1  p=0.4186

lifespan.model.m.wb.9 <- update(lifespan.model.m.wb, ~ . -cGroupSize)

anova(lifespan.model.m.wb.9, lifespan.model.m.wb)


#~~ helpers - x=1.9361  df=1  p=0.1641

lifespan.model.m.wb.10 <- update(lifespan.model.m.wb, ~ . -HelpersF)

anova(lifespan.model.m.wb.10, lifespan.model.m.wb)


#~~ siblings - x=1.9906  df=1  p=0.1583

lifespan.model.m.wb.11 <- update(lifespan.model.m.wb, ~ . -Siblings)

anova(lifespan.model.m.wb.11, lifespan.model.m.wb)


#~~ are these within vs between slopes significantly different?

lifespan.model.m.wb.b <- glmmTMB(LifespanInt~cMeanAgeMother + cAgeMother +
                                   cMeanAgeFather + cAgeFather +
                                   cMeanAgeSocMale + cDominantMaleAge +
                                   cBirthYear + clogTQ + cGroupSize + HelpersF + Siblings +
                                   (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                                 data=lifespan.data.dead.m,
                                 ziformula=~0,
                                 family=poisson)


summary(lifespan.model.m.wb.b)


#~~ save model output

lifespanmmodelwb.b.fixef<-as.data.frame((summary(lifespan.model.m.wb.b)$coefficients$cond))
lifespanmmodelwb.b.fixef<-setDT(lifespanmmodelwb.b.fixef, keep.rownames = TRUE)[]

write.csv(lifespanmmodelwb.b.fixef, file="Models/Lifespan_male_model_wnbn_b.csv")


#~~ significance of fixed effects

#~~ difference between within and between maternal age slopes - x=0.0686  df=1  p=0.7935

lifespan.model.m.wb.b2 <- update(lifespan.model.m.wb.b, ~ . -cMeanAgeMother)

anova(lifespan.model.m.wb.b2, lifespan.model.m.wb.b)


#~~ within maternal age effect - x=0.9199  df=1  p=0.3375

lifespan.model.m.wb.b3 <- update(lifespan.model.m.wb.b, ~ . -cAgeMother)

anova(lifespan.model.m.wb.b3, lifespan.model.m.wb.b)


#~~ difference between within and between paternal age slopes - x=2.5027  df=1  p=0.1137 

lifespan.model.m.wb.b4 <- update(lifespan.model.m.wb.b, ~ . -cMeanAgeFather)

anova(lifespan.model.m.wb.b4, lifespan.model.m.wb.b)


#~~ within paternal age effect - x=2.437  df=1  p=0.1185

lifespan.model.m.wb.b5 <- update(lifespan.model.m.wb.b, ~ . -cAgeFather)

anova(lifespan.model.m.wb.b5, lifespan.model.m.wb.b)


#~~ difference between within and between dominant male age slopes - x=0.1258  df=1  p=0.7229

lifespan.model.m.wb.b6 <- update(lifespan.model.m.wb.b, ~ . -cMeanAgeSocMale)

anova(lifespan.model.m.wb.b6, lifespan.model.m.wb.b)


#~~ within dominant male age effect - x=0.6038  df=1  p=0.4371

lifespan.model.m.wb.b7 <- update(lifespan.model.m.wb.b, ~ . -cDominantMaleAge)

anova(lifespan.model.m.wb.b7, lifespan.model.m.wb.b)


#~~ birth year - x=12.252  df=1  p=0.0004647

lifespan.model.m.wb.b8 <- update(lifespan.model.m.wb.b, ~ . -cBirthYear)

anova(lifespan.model.m.wb.b8, lifespan.model.m.wb.b)


#~~ territory quality - x=2.6452  df=1  p=0.1039

lifespan.model.m.wb.b9 <- update(lifespan.model.m.wb.b, ~ . -clogTQ)

anova(lifespan.model.m.wb.b9, lifespan.model.m.wb.b)


#~~ group size - x=0.6543  df=1  p=0.4186

lifespan.model.m.wb.b10 <- update(lifespan.model.m.wb.b, ~ . -cGroupSize)

anova(lifespan.model.m.wb.b10, lifespan.model.m.wb.b)


#~~ helpers - x=1.9361  df=1  p=0.1641

lifespan.model.m.wb.b11 <- update(lifespan.model.m.wb.b, ~ . -HelpersF)

anova(lifespan.model.m.wb.b11, lifespan.model.m.wb.b)


#~~ siblings - x=1.9906  df=1  p=0.1583

lifespan.model.m.wb.b12 <- update(lifespan.model.m.wb.b, ~ . -Siblings)

anova(lifespan.model.m.wb.b12, lifespan.model.m.wb.b)


# # -----------------------------------------------------------------------

## Lifespan coxme models - females

# # -----------------------------------------------------------------------

#~~ subset dataset for coxme analyses - all individuals who lived to one 

lifespan.data <- subset(fitnessdata, LifespanInt > 0)


#~~ female offspring subset
# check sample size

lifespan.data.fem <- subset(lifespan.data, Sex == 0 & !is.na(AgeMother) & !is.na(DominantMaleAge) 
                                  & !is.na(BirthYearF) & !is.na(logTQ)
                                  & !is.na(HelpersF) & !is.na(GroupSize) & !is.na(AgeFather))


# NOTE: females translocated or still alive are censored in the model

table(lifespan.data.fem$Event)


#~~ base model

lspanfmodela <- coxme(Surv(LifespanInt, Event) ~ AgeMother + AgeFather + DominantMaleAge + logTQ + GroupSize + HelpersF + Siblings +
                        (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF), data = lifespan.data.fem)

lspanfmodela

beep()

save(lspanfmodela, file="Models/Lifespan_coxme_model_females.RData")


#~~ check vif

vif.mer(lspanfmodela)


#~~ significance of fixed effects

#~~ mum age sq - x=0.008  df=1  p=0.9286

lspanfmodelb <- update(lspanfmodela, ~ . +AgeMotherSq)

anova(lspanfmodela, lspanfmodelb)


#~~ dad age sq - x=1.9541  df=1  p=0.1621

lspanfmodelb <- update(lspanfmodela, ~ . +AgeFatherSq)

anova(lspanfmodela, lspanfmodelc)


#~~ soc male age sq - x=0.016  df=1  p=0.8994

lspanfmodeld <- update(lspanfmodela, ~ . +AgeSocFatherSq)

anova(lspanfmodela, lspanfmodeld)


#~~ age mum - x=3.8709  df=1  p=0.04913

lspanfmodele <- update(lspanfmodela, ~ . -AgeMother)

anova(lspanfmodela, lspanfmodele)


#~~ age dad - x=0.013  df=1  p=0.9092

lspanfmodelf <- update(lspanfmodela, ~ . -AgeFather)

anova(lspanfmodela, lspanfmodelf)


#~~ soc male age - x=0.0175  df=1  p=0.8948

lspanfmodelg <- update(lspanfmodela, ~ . -DominantMaleAge)

anova(lspanfmodela, lspanfmodelg)


#~~ tq - x=4.707  df=1  p=0.03004

lspanfmodelh <- update(lspanfmodela, ~ . -logTQ)

anova(lspanfmodela, lspanfmodelh)


#~~ group size - x=2.3595  df=1  p=0.1245

lspanfmodeli <- update(lspanfmodela, ~ . -GroupSize)

anova(lspanfmodela, lspanfmodeli)


#~~ helpers - x=2.866  df=1  p=0.09047 

lspanfmodelj <- update(lspanfmodela, ~ . -HelpersF)

anova(lspanfmodela, lspanfmodelj)


#~~ siblings - x=0.0055  df=1  p=0.9408

lspanfmodelk <- update(lspanfmodela, ~ . -Siblings)

anova(lspanfmodela, lspanfmodelk)


#~~ TQ * mum age - x=0.0444  df=1  p=0.8331

lspanfmodell <- update(lspanfmodela, ~ . +logTQ*AgeMother)

anova(lspanfmodela, lspanfmodell)


#~~ TQ * dad age - x=0.1244  df=1  p=0.7243

lspanfmodelm <- update(lspanfmodela, ~ . +logTQ*AgeFather)

anova(lspanfmodela, lspanfmodelm)


#~~ TQ * social male age - x=0.0629  df=1  p=0.802

lspanfmodeln <- update(lspanfmodela, ~ . +logTQ*DominantMaleAge)

anova(lspanfmodela, lspanfmodeln)


#~~ group size * mum age - x=0.7484  df=1  p=0.387

lspanfmodelo <- update(lspanfmodela, ~ . +GroupSize*AgeMother)

anova(lspanfmodela, lspanfmodelo)


#~~ group size * dad age - x=1.3105  df=1  p=0.2523

lspanfmodelp <- update(lspanfmodela, ~ . +GroupSize*AgeFather)

anova(lspanfmodela, lspanfmodelp)


#~~ group size * social male age - x=0.0102  df=1  p=0.9194

lspanfmodelq <- update(lspanfmodela, ~ . +GroupSize*DominantMaleAge)

anova(lspanfmodela, lspanfmodelq)


#~~ helpers * mum age - x=0.0069  df=1  p=0.9339

lspanfmodelr <- update(lspanfmodela, ~ . +HelpersF*AgeMother)

anova(lspanfmodela, lspanfmodelr)


#~~ helpers * dad age - x=2.0422  df=1  p=0.153

lspanfmodels <- update(lspanfmodela, ~ . +HelpersF*AgeFather)

anova(lspanfmodela, lspanfmodels)


#~~ helpers * social male age - x=0.076  df=1  p=0.7828

lspanfmodelt <- update(lspanfmodela, ~ . +HelpersF*DominantMaleAge)

anova(lspanfmodela, lspanfmodelt)


#~~ maternal age * paternal age - x=0.4246  df=1  p=0.5147

lspanfmodelu <- update(lspanfmodela, ~ . +AgeMother*AgeFather)

anova(lspanfmodela, lspanfmodelu)


#~~ maternal age * dominant male age - x=2.2069  df=1  p=0.1374

lspanfmodelv <- update(lspanfmodela, ~ . +AgeMother*DominantMaleAge)

anova(lspanfmodela, lspanfmodelv)


#~~ check assumptions of final model 

coxme_fem <- coxph(Surv(LifespanInt, Event) ~ AgeMother + AgeFather + DominantMaleAge + logTQ + GroupSize + HelpersF + Siblings, data=lifespan.data.fem)

female.test.ph <- cox.zph(coxme_fem)
female.test.ph

ggcoxzph(female.test.ph)
ggcoxdiagnostics(coxme_fem, type="dfbeta", linear.predictions=FALSE, ggtheme=theme_bw())
ggcoxdiagnostics(coxme_fem, type="deviance", linear.predictions=FALSE, ggtheme=theme_bw())


# # -----------------------------------------------------------------------

## Lifespan coxme models - males

# # -----------------------------------------------------------------------

#~~ male offspring subset
# check sample size

lifespan.data.male <- subset(lifespan.data, Sex == 1 & !is.na(AgeMother) & !is.na(DominantMaleAge) 
                                  & !is.na(BirthYearF) & !is.na(logTQ)
                                  & !is.na(HelpersF) & !is.na(GroupSize) & !is.na(AgeFather))


# NOTE: males translocated or still alive are censored in the model

table(lifespan.data.male$Event)


#~~ base model

lspanmmodela <- coxme(Surv(LifespanInt, Event) ~ AgeMother + AgeFather + DominantMaleAge + logTQ + GroupSize + HelpersF + Siblings + 
                        (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF), data = lifespan.data.male)
lspanmmodela

beep()

save(lspanmmodela, file="Models/Lifespan_coxme_model_males.RData")


#~~ check vif

vif.mer(lspanmmodela)


#~~ significance of fixed effects

#~~ mum age sq - x=0.1045  df=1  p=0.7464

lspanmmodelb <- update(lspanmmodela, ~ . +AgeMotherSq)

anova(lspanmmodela, lspanmmodelb)


#~~ dad age sq - x=0.1515  df=1  p=0.6971

lspanmmodelc <- update(lspanmmodela, ~ . +AgeFatherSq)

anova(lspanmmodela, lspanmmodelc)


#~~ soc male age sq - x=0.8783  df=1  p=0.3487

lspanmmodeld <- update(lspanmmodela, ~ . +AgeSocFatherSq)

anova(lspanmmodela, lspanmmodeld)


#~~ age mum - x=0.2791  df=1  p=0.5973

lspanmmodele <- update(lspanmmodela, ~ . -AgeMother)

anova(lspanmmodela, lspanmmodele)


#~~ age dad - x=0.023  df=1  p=0.8795

lspanmmodelf <- update(lspanmmodela, ~ . -AgeFather)

anova(lspanmmodela, lspanmmodelf)


#~~ soc male age - x=3.6529  df=1  p=0.05597

lspanmmodelg <- update(lspanmmodela, ~ . -DominantMaleAge)

anova(lspanmmodela, lspanmmodelg)


#~~ tq - x=0.8909  df=1  p=0.3452

lspanmmodelh <- update(lspanmmodela, ~ . -logTQ)

anova(lspanmmodela, lspanmmodelh)


#~~ group size - x=3.651  df=1  p=0.05603

lspanmmodeli <- update(lspanmmodela, ~ . -GroupSize)

anova(lspanmmodela, lspanmmodeli)


#~~ helpers - x=3.1135  df=1  p=0.07765

lspanmmodelj <- update(lspanmmodela, ~ . -HelpersF)

anova(lspanmmodela, lspanmmodelj)


#~~ siblings - x=1.2858  df=1  p=0.2568

lspanmmodelk <- update(lspanmmodela, ~ . -Siblings)

anova(lspanmmodela, lspanmmodelk)


#~~ TQ * mum age - x=0.4195  df=1  p=0.5172

lspanmmodell <- update(lspanmmodela, ~ . +logTQ*AgeMother)

anova(lspanmmodela, lspanmmodell)


#~~ TQ * dad age - x=1.8421  df=1  p=0.1747

lspanmmodelm <- update(lspanmmodela, ~ . +logTQ*AgeFather)

anova(lspanmmodela, lspanmmodelm)


#~~ TQ * social male age - x=0.1465  df=1  p=0.7019

lspanmmodeln <- update(lspanmmodela, ~ . +logTQ*DominantMaleAge)

anova(lspanmmodela, lspanmmodeln)


#~~ group size * mum age - x=0.0315  df=1  p=0.8592

lspanmmodelo <- update(lspanmmodela, ~ . +GroupSize*AgeMother)

anova(lspanmmodela, lspanmmodelo)


#~~ group size * dad age - x=0.0111  df=1  p=0.9159

lspanmmodelp <- update(lspanmmodela, ~ . +GroupSize*AgeFather)

anova(lspanmmodela, lspanmmodelp)


#~~ group size * social male age - x=0.3812  df=1  p=0.537

lspanmmodelq <- update(lspanmmodela, ~ . +GroupSize*DominantMaleAge)

anova(lspanmmodela, lspanmmodelq)


#~~ helpers * mum age - x=0.7841  df=1  p=0.3759

lspanmmodelr <- update(lspanmmodela, ~ . +HelpersF*AgeMother)

anova(lspanmmodela, lspanmmodelr)


#~~ helpers * dad age - x=0.692  df=1  p=0.4055

lspanmmodels <- update(lspanmmodela, ~ . +HelpersF*AgeFather)

anova(lspanmmodela, lspanmmodels)


#~~ helpers * social male age - x=0.0196  df=1  p=0.8887

lspanmmodelt <- update(lspanmmodela, ~ . +HelpersF*DominantMaleAge)

anova(lspanmmodela, lspanmmodelt)


#~~ maternal * paternal age - x=0.1313  df=1  p=0.7171

lspanmmodelu <- update(lspanmmodela, ~ . +AgeMother*AgeFather)

anova(lspanmmodela, lspanmmodelu)


#~~ maternal * dominant male age - x=0.5391  df=1  p=0.4628

lspanmmodelv <- update(lspanmmodela, ~ . +AgeMother*DominantMaleAge)

anova(lspanmmodela, lspanmmodelv)


#~~ check assumptions of final model 

coxme_male <- coxph(Surv(LifespanInt, Event) ~ AgeMother + AgeFather + DominantMaleAge + logTQ + GroupSize + HelpersF + Siblings, data=lifespan.data.male)

male.test.ph <- cox.zph(coxme_male)
male.test.ph

ggcoxzph(male.test.ph)
ggcoxdiagnostics(coxme_male, type="dfbeta", linear.predictions=FALSE, ggtheme=theme_bw())
ggcoxdiagnostics(coxme_male, type="deviance", linear.predictions=FALSE, ggtheme=theme_bw())


# # -----------------------------------------------------------------------

# Reproductive lifespan models - females

# # -----------------------------------------------------------------------

#~~ Reproductive lifespan data subset
#~~ NOTE: the dataset is already cropped to individuals who made it to at least 1 year
#~~ Reproductive lifespan measures are only available for birds that have completed their reproductive lifespan 
#~~ i.e. who are dead and whose final year was < 2018 (last complete year of the pedigree was 2017) and birth year > 1994
#~~ Reproductive lifespans are calculated irrespective of whether the individual became a dominant breeder or not
#~~ However, individuals who were translocated need to be removed from the dataset as they have incomplete lifetime data

reprolifespandata <- droplevels(subset(fitnessdata, LifespanInt > 0 & !is.na(Sex) & FinalYear < 2018 & TranslocatedF==0))


#~~ female model subset

reprolifespan.data.f <- droplevels(subset(reprolifespandata, Sex==0 & !is.na(GroupSize) & !is.na(AgeMother) & !is.na(AgeFather) & !is.na(DominantMaleAge) 
                              & !is.na(BirthYear) & !is.na(Siblings) & !is.na(HelpersF)))

#~~ data

table(reprolifespan.data.f$ReproLifespan)


#~~ model - run a variety of GLMMs in glmmTMB and see which model fits best

rls.f.zipoisson <- glmmTMB(ReproLifespan~ rescale(AgeMother) + rescale(AgeFather) + rescale(DominantMaleAge) +
                            rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                            Siblings + (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                          data=reprolifespan.data.f,
                          ziformula=~1,
                          family=poisson)
summary(rls.f.zipoisson)

rls.f.zinbinom2 <- update(rls.f.zipoisson, family=nbinom2)

rls.f.zinbinom1 <- update(rls.f.zipoisson, family=nbinom1)

rls.f.poisson <- update(rls.f.zipoisson, ziformula=~0, family=poisson)

rls.f.nbinom2 <- update(rls.f.zipoisson, ziformula=~0, family=nbinom2)

rls.f.nbinom1 <- update(rls.f.zipoisson, ziformula=~0, family=nbinom1)

# all models converged


#~~ AIC comparisons

AIC(rls.f.zipoisson, rls.f.zinbinom2, rls.f.zinbinom1, rls.f.poisson, rls.f.nbinom2, rls.f.nbinom1)


#~~ check model residuals

rls.f.zipoisson.resid <- simulateResiduals(fittedModel = rls.f.zipoisson, plot = T) # fine

rls.f.zinbinom2.resid <- simulateResiduals(fittedModel = rls.f.zinbinom2, plot = T) # can't use this one

rls.f.zinbinom1.resid <- simulateResiduals(fittedModel = rls.f.zinbinom1, plot = T) # can't use this one

rls.f.poisson.resid <- simulateResiduals(fittedModel = rls.f.poisson, plot = T) # can't use this one

rls.f.nbinom2.resid <- simulateResiduals(fittedModel = rls.f.nbinom2, plot = T) # fine

rls.f.nbinom1.resid <- simulateResiduals(fittedModel = rls.f.nbinom1, plot = T) # fine


#~~ check zero inflation with 3 best models - ZI poisson accounts for observed zeroes best

check_zeroinflation(rls.f.zipoisson)
check_zeroinflation(rls.f.nbinom2)
check_zeroinflation(rls.f.nbinom1)


#~~ go with zero inflated poisson model
#~~ base model

rls.f.model <- glmmTMB(ReproLifespan~ rescale(AgeMother) + rescale(AgeFather) + rescale(DominantMaleAge) +
                             rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                             Siblings + (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                           data=reprolifespan.data.f,
                           ziformula=~1,
                           family=poisson)
summary(rls.f.model)


#~~ model checks

check_collinearity(rls.f.model, component="count")
check_zeroinflation(rls.f.model)
rls.f.model.resid <- simulateResiduals(fittedModel = rls.f.model, plot = T)


#~~ save output

rls.f.model.fixef<-as.data.frame((summary(rls.f.model)$coefficients$cond))
rls.f.model.fixef<-setDT(rls.f.model.fixef, keep.rownames = TRUE)[]

write.csv(rls.f.model.fixef, file="Models/Repro_lifespan_female_model.csv")


#~~ significance of fixed effects

#~~ maternal age - x=9.9316  df=1  p=0.001625

rls.f.model2 <- update(rls.f.model, ~ . -rescale(AgeMother))

anova(rls.f.model, rls.f.model2)


#~~ paternal age - x=0.5753  df=1  p=0.4481

rls.f.model3 <- update(rls.f.model, ~ . -rescale(AgeFather))

anova(rls.f.model, rls.f.model3)


#~~ dominant male age - x=0.0976  df=1  p=0.7547

rls.f.model4 <- update(rls.f.model, ~ . -rescale(DominantMaleAge))

anova(rls.f.model, rls.f.model4)


#~~ birth year - x=32.013  df=1  p=1.531e-08

rls.f.model5 <- update(rls.f.model, ~ . -rescale(BirthYear))

anova(rls.f.model, rls.f.model5)


#~~ tq - x=1.7209  df=1  p=0.1896

rls.f.model6 <- update(rls.f.model, ~ . -rescale(logTQ))

anova(rls.f.model, rls.f.model6)


#~~ group size - x=3.9775  df=1  p=0.04611

rls.f.model7 <- update(rls.f.model, ~ . -rescale(GroupSize))

anova(rls.f.model, rls.f.model7)


#~~ helpers - x=3.5902  df=1  p=0.05812

rls.f.model8 <- update(rls.f.model, ~ . -HelpersF)

anova(rls.f.model, rls.f.model8)


#~~ sibling - x=1.0021  df=1  p=0.3168

rls.f.model9 <- update(rls.f.model, ~ . -Siblings)

anova(rls.f.model, rls.f.model9)


# # -----------------------------------------------------------------------

# Reproductive lifespan models - males

# # -----------------------------------------------------------------------

#~~ male model subset

reprolifespan.data.m <- droplevels(subset(reprolifespandata, Sex==1 & !is.na(GroupSize) & !is.na(AgeMother) & !is.na(AgeFather) & !is.na(DominantMaleAge) 
                                          & !is.na(BirthYear) & !is.na(Siblings) & !is.na(HelpersF)))

#~~ data

table(reprolifespan.data.m$ReproLifespan)


#~~ model - run a variety of GLMMs in glmmTMB and see which model fits best

rls.m.zipoisson <- glmmTMB(ReproLifespan~ rescale(AgeMother) + rescale(AgeFather) + rescale(DominantMaleAge) +
                             rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                             Siblings + (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                           data=reprolifespan.data.m,
                           ziformula=~1,
                           family=poisson)
summary(rls.m.zipoisson)

rls.m.zinbinom2 <- update(rls.m.zipoisson, family=nbinom2)

rls.m.zinbinom1 <- update(rls.m.zipoisson, family=nbinom1)

rls.m.poisson <- update(rls.m.zipoisson, ziformula=~0, family=poisson)

rls.m.nbinom2 <- update(rls.m.zipoisson, ziformula=~0, family=nbinom2)

rls.m.nbinom1 <- update(rls.m.zipoisson, ziformula=~0, family=nbinom1)

# all models converged


#~~ AIC comparisons

AIC(rls.m.zipoisson, rls.m.zinbinom2, rls.m.zinbinom1, rls.m.poisson, rls.m.nbinom2, rls.m.nbinom1)


#~~ check model residuals

rls.m.zipoisson.resid <- simulateResiduals(fittedModel = rls.m.zipoisson, plot = T) # fine

rls.m.zinbinom2.resid <- simulateResiduals(fittedModel = rls.m.zinbinom2, plot = T) # fine

rls.m.zinbinom1.resid <- simulateResiduals(fittedModel = rls.m.zinbinom1, plot = T) # fine

rls.m.poisson.resid <- simulateResiduals(fittedModel = rls.m.poisson, plot = T) # not fine

rls.m.nbinom2.resid <- simulateResiduals(fittedModel = rls.m.nbinom2, plot = T) # fine

rls.m.nbinom1.resid <- simulateResiduals(fittedModel = rls.m.nbinom1, plot = T) # fine


#~~ check zero inflation with best models - ZI poisson accounts for observed zeroes best

check_zeroinflation(rls.m.zipoisson)
check_zeroinflation(rls.m.zinbinom2)
check_zeroinflation(rls.m.zinbinom1)
check_zeroinflation(rls.m.nbinom2)
check_zeroinflation(rls.m.nbinom1)

                    
#~~ go with zero inflated poisson model
#~~ base model

rls.m.model <- glmmTMB(ReproLifespan~ rescale(AgeMother) + rescale(AgeFather) + rescale(DominantMaleAge) +
                         rescale(BirthYear) + rescale(logTQ) + rescale(GroupSize) + HelpersF +
                         Siblings + (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                       data=reprolifespan.data.m,
                       ziformula=~1,
                       family=poisson)
summary(rls.m.model)


#~~ model checks

check_collinearity(rls.m.model, component="count")
check_zeroinflation(rls.m.model)
rls.m.model.resid <- simulateResiduals(fittedModel = rls.m.model, plot = T)


#~~ save output

rls.m.model.fixef<-as.data.frame((summary(rls.m.model)$coefficients$cond))
rls.m.model.fixef<-setDT(rls.m.model.fixef, keep.rownames = TRUE)[]

write.csv(rls.m.model.fixef, file="Models/Repro_lifespan_male_model.csv")


#~~ maternal age - x=3.5711  df=1  p=0.05879

rls.m.model2 <- update(rls.m.model, ~ . -rescale(AgeMother))

anova(rls.m.model, rls.m.model2)


#~~ paternal age - x=0.2946  df=1  p=0.5873

rls.m.model3 <- update(rls.m.model, ~ . -rescale(AgeFather))

anova(rls.m.model, rls.m.model3)


#~~ dominant male age - x=0.5857  df=1  p=0.4441

rls.m.model4 <- update(rls.m.model, ~ . -rescale(DominantMaleAge))

anova(rls.m.model, rls.m.model4)


#~~ birth year - x=21.23  df=1  p==4.073e-06 

rls.m.model5 <- update(rls.m.model, ~ . -rescale(BirthYear))

anova(rls.m.model, rls.m.model5)


#~~ tq - x=1.1926  df=1  p=0.2748

rls.m.model6 <- update(rls.m.model, ~ . -rescale(logTQ))

anova(rls.m.model, rls.m.model6)

#~~ group size - x=0.1761  df=1  p=0.6747

rls.m.model7 <- update(rls.m.model, ~ . -rescale(GroupSize))

anova(rls.m.model, rls.m.model7)


#~~ helpers - x=1.0022  df=1  p=0.3168

rls.m.model8 <- update(rls.m.model, ~ . -HelpersF)

anova(rls.m.model, rls.m.model8)


#~~ sibling - x=0.7113  df=1  p=0.399

rls.m.model9 <- update(rls.m.model, ~ . -Siblings)

anova(rls.m.model, rls.m.model9)