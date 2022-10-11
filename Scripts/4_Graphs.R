# # ------------------------------------------------------------------------

# ---- Parental age manuscript
# ---- AMS
# ---- Graphs

# # ------------------------------------------------------------------------

rm(list=ls())


#~~ libraries

library(ggplot2)
library(cowplot)
library(survival)
library(survminer)
library(arm)
library(car)
library(glmmTMB)
library(dplyr)
library(interactions)


#~~ data

fitnessdata <- read.csv("Data/Sparks_Parental_Age_Lifespan_Data_2022.csv", stringsAsFactors = FALSE)


#~~ set factors

fitnessdata$MotherID<-as.factor(fitnessdata$MotherID)
fitnessdata$FatherID<-as.factor(fitnessdata$FatherID)
fitnessdata$DominantMaleID<-as.factor(fitnessdata$DominantMaleID)
fitnessdata$ID<-as.factor(fitnessdata$ID)
fitnessdata$BirthYearF<-as.factor(fitnessdata$BirthYear)
fitnessdata$Sex<-as.factor(fitnessdata$Sex)
fitnessdata$AgeMotherInt <- as.integer(fitnessdata$AgeMother)
fitnessdata$AgeClassMother<-car::recode(fitnessdata$AgeMotherInt,"0:3=0;4:6 =1; 7:9=2; 10:20=3")
fitnessdata$AgeClassMother<-as.factor(fitnessdata$AgeClassMother)
fitnessdata$HelpersF <- as.factor(fitnessdata$HelpersF)
fitnessdata$TranslocatedF  <-  as.factor(fitnessdata$TranslocatedF)
fitnessdata$logTQ <- log10(fitnessdata$TQ)
fitnessdata$AgeMotherSq <- fitnessdata$AgeMother * fitnessdata$AgeMother
fitnessdata$AgeFatherSq <- fitnessdata$AgeFather * fitnessdata$AgeFather
fitnessdata$AgeSocFatherSq <- fitnessdata$DominantMaleAge * fitnessdata$DominantMaleAge
fitnessdata$Siblings <- as.factor(fitnessdata$Siblings)
fitnessdata$LifespanInt <- fitnessdata$FinalYear - fitnessdata$BirthYear


#~~ summary function

##To get a table with means sd se and ci:
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


# # -----------------------------------------------------------------------

# Figure S1 - histograms of LRS and lifespan

# # -----------------------------------------------------------------------

#~~ summary stats

lifespanstats <- subset(fitnessdata, LifespanInt > 0 & !is.na(Sex) & Event==1)

lifespanbysex <- summarySE(lifespanstats, measurevar = "LifespanInt", groupvars = "Sex")
lifespanbysex

lrsstats <- subset(fitnessdata, LifespanInt > 0 & !is.na(Sex) & FinalYear < 2018 & TranslocatedF==0)

lrsbysex <- summarySE(lrsstats, measurevar = "LRS1YEAR", groupvars = "Sex")
lrsbysex


#~~ graphs

lifespanf <- subset(lifespanstats, Sex==0)
lifespanm <- subset(lifespanstats, Sex==1)
lrsf <- subset(lrsstats, Sex==0)
lrsm <- subset(lrsstats, Sex==1)

fig.s1a <- ggplot(lifespanf, aes(x=LifespanInt)) +
  geom_histogram(binwidth=1, alpha=0.8, colour="black", fill="darkorange") +
  xlab("Female lifespan") + ylab("Count") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey1", size=12),
        axis.title.y=element_text(colour="grey1", size=12), 
        axis.line=element_line(colour="grey6"), 
        text = element_text(size=12, colour="grey6")) 

fig.s1b <- ggplot(lifespanm, aes(x=LifespanInt)) +
  geom_histogram(binwidth=1, alpha=0.6, colour="black", fill="#999999") +
  xlab("Male lifespan") + ylab("Count") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey1", size=12),
        axis.title.y=element_text(colour="grey1", size=12), 
        axis.line=element_line(colour="grey6"), 
        text = element_text(size=12, colour="grey6")) 

fig.s1c <- ggplot(lrsf, aes(x=LRS1YEAR)) +
  geom_histogram(binwidth=1, alpha=0.8, colour="black", fill="darkorange") +
  xlab("Female LRS") + ylab("Count") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey1", size=12),
        axis.title.y=element_text(colour="grey1", size=12), 
        axis.line=element_line(colour="grey6"), 
        text = element_text(size=12, colour="grey6")) + 
  scale_x_continuous(labels=c(0, 3, 6, 9), breaks=c(0,3,6,9))

fig.s1d <- ggplot(lrsm, aes(x=LRS1YEAR)) +
  geom_histogram(binwidth=1, alpha=0.6, colour="black", fill="#999999") +
  xlab("Male LRS") + ylab("Count") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey1", size=12),
        axis.title.y=element_text(colour="grey1", size=12), 
        axis.line=element_line(colour="grey6"), 
        text = element_text(size=12, colour="grey6")) 

lrslspanhists <- plot_grid(fig.s1a, fig.s1b, fig.s1c, fig.s1d, labels="AUTO", ncol=2, nrow=2)

tiff("Graphs/FigureS1_Histograms.tif", width=3000, height=2500, res=400)

lrslspanhists

dev.off()


# # -----------------------------------------------------------------------

# Figure S2 - correlation between maternal age and paternal age 

# # -----------------------------------------------------------------------

#~~ correlation stats

cor.test(fitnessdata$AgeMother, fitnessdata$AgeFather)
cor.test(fitnessdata$AgeMother, fitnessdata$DominantMaleAge)
cor.test(fitnessdata$AgeFather, fitnessdata$DominantMaleAge)


#~~ figure

fig.s2<-ggplot(fitnessdata ,aes(x=AgeMother, y=AgeFather, color=Sex)) + geom_point(alpha=0.5) +
  xlab("Age of mother (years)") + ylab("Age of father (years)") + labs(colour="Offspring sex") + theme_classic() +
  scale_color_manual(labels = c("Female", "Male"), values=c("darkorange", "#999999")) +
  theme(axis.title.x=element_text(colour="grey1", size=13), 
        axis.title.y=element_text(colour="grey1", size=13), 
        axis.line=element_line(colour="grey60"),
        text = element_text(size=15, colour="grey26"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
  
fig.s2

tiff("Graphs/FigureS2_ParAgeCorrelations.tif", width=2500, height=2000, res=400)

fig.s2

dev.off()


# # -----------------------------------------------------------------------

# Figure S5 - correlation between lifespan and LRS 

# # -----------------------------------------------------------------------

#~~ correlation stats

cor.test(lrsstats$LifespanInt, lrsstats$LRS1YEAR)


#~~ figure

fig.s5 <- ggplot(lrsstats,aes(x=LifespanInt,y=LRS1YEAR, color=Sex)) +
  geom_jitter(alpha=0.5, width=0.35, height=0.35) +
  xlab("Lifespan (years)") + ylab("Lifetime reproductive success") + theme_classic() +
  scale_color_manual(labels = c("Female", "Male"), values=c("darkorange", "#999999")) +
  theme(axis.title.x=element_text(colour="grey16", size=13), 
        axis.title.y=element_text(colour="grey16", size=13), 
        axis.line=element_line(colour="grey60"),
        text = element_text(size=15, colour="grey26"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  scale_y_continuous(breaks=c(0,2,4,6,8,10))

fig.s5

tiff("Graphs/FigureS5_LRSLifespanCorr.tif", width=2500, height=2000, res=400)

fig.s5

dev.off()


# # -----------------------------------------------------------------------

# Figure 1 - survival plot for female offspring split by maternal age

# # -----------------------------------------------------------------------

#~~ lifespan female dataset
# same dataset included for the coxme models - so individuals still alive or translocated are censored

lifespan.coxme.data.fem <- subset(fitnessdata, LifespanInt > 0 & Sex == 0 & !is.na(AgeMother) & !is.na(DominantMaleAge) 
                            & !is.na(BirthYearF) & !is.na(logTQ)
                            & !is.na(HelpersF) & !is.na(GroupSize) & !is.na(AgeFather))


#~~ model

fem.coxme.model <- survfit(Surv(LifespanInt, Event) ~ AgeClassMother, data = lifespan.coxme.data.fem)


#~~ graph

fig.1 <- ggsurvplot(fem.coxme.model, data = lifespan.coxme.data.fem,
                       conf.int = FALSE,
                       pval = FALSE,
                       fun = "pct",
                       risk.table = FALSE,
                       font.y=12,
                       font.x=12,
                       font.tickslab=11,
                       size = 1,
                       linetype = "strata",
                       palette = c("#e66101", "#fdb863", 
                                   "#b2abd2",  "#5e3c99"),
                       legend = "bottom",
                       xlab="Time (years)",
                       ylab="Female offspring survival probability (%)",
                       legend.title = "Maternal age (years)",
                       legend.labs = c("0-3","4-6", "7-9","10+"),
                       font.legend=11,
                       censor.shape="")

fig.1

pdf("Graphs/Figure1_SurvPlot.pdf", width=5.25, height=5.5, onefile=FALSE)

print(fig.1)

dev.off()


# # -----------------------------------------------------------------------

# Figure S3 - within-individual maternal age effects on female offspring lifespan

# # -----------------------------------------------------------------------

#~~ data subset (see main lifespan script for more info)

lifespan.data.dead.f <- droplevels(subset(fitnessdata, LifespanInt > 0 & Sex==0 & FinalYear < 2019 & TranslocatedF==0 
                                          & !is.na(GroupSize) & !is.na(AgeMother) & !is.na(AgeFather) & !is.na(DominantMaleAge) 
                                          & !is.na(BirthYear) & !is.na(Siblings) & !is.na(HelpersF)))


#~~ z-transform covariates

lifespan.data.dead.f$cAgeMother <- rescale(lifespan.data.dead.f$AgeMother)
lifespan.data.dead.f$cAgeFather <- rescale(lifespan.data.dead.f$AgeFather)
lifespan.data.dead.f$cDominantMaleAge <- rescale(lifespan.data.dead.f$DominantMaleAge)
lifespan.data.dead.f$cBirthYear <- rescale(lifespan.data.dead.f$BirthYear)
lifespan.data.dead.f$clogTQ <- rescale(lifespan.data.dead.f$logTQ)
lifespan.data.dead.f$cGroupSize <- rescale(lifespan.data.dead.f$GroupSize)


#~~ make mean maternal, paternal and dominant male age variables and deviation of each value from the mean

lifespan.data.dead.f<-lifespan.data.dead.f %>% dplyr::group_by(MotherID) %>% dplyr::mutate(cMeanAgeMother = mean(cAgeMother)) %>% dplyr::ungroup()
lifespan.data.dead.f$cDevMeanAgeMother<-lifespan.data.dead.f$cAgeMother-lifespan.data.dead.f$cMeanAgeMother

lifespan.data.dead.f<-lifespan.data.dead.f %>% dplyr::group_by(FatherID) %>% dplyr::mutate(cMeanAgeFather = mean(cAgeFather)) %>% dplyr::ungroup()
lifespan.data.dead.f$cDevMeanAgeFather<-lifespan.data.dead.f$cAgeFather-lifespan.data.dead.f$cMeanAgeFather

lifespan.data.dead.f<-lifespan.data.dead.f %>% dplyr::group_by(DominantMaleID) %>% dplyr::mutate(cMeanAgeSocMale = mean(cDominantMaleAge)) %>% dplyr::ungroup()
lifespan.data.dead.f$cDevMeanAgeSocMale<-lifespan.data.dead.f$cDominantMaleAge-lifespan.data.dead.f$cMeanAgeSocMale

lifespan.data.dead.f<-lifespan.data.dead.f %>% dplyr::group_by(MotherID) %>% dplyr::mutate(MeanAgeMother = mean(AgeMother)) %>% dplyr::ungroup()
lifespan.data.dead.f$DevMeanAgeMother<-lifespan.data.dead.f$AgeMother-lifespan.data.dead.f$MeanAgeMother


#~~ model

lifespan.model.f.wb <- glmmTMB(LifespanInt~cMeanAgeMother + cDevMeanAgeMother +
                                 cMeanAgeFather + cDevMeanAgeFather + 
                                 cMeanAgeSocMale + cDevMeanAgeSocMale +
                                 cBirthYear + clogTQ + cGroupSize + HelpersF + Siblings +
                                 (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                               data=lifespan.data.dead.f,
                               ziformula=~0,
                               family=poisson)

summary(lifespan.model.f.wb)


#~~ make predict dataframe
#~~ all continuous variables set to mean values except for within-individual maternal age effects
#~~ estimated for offspring with no helpers or siblings in the natal nest

predDF.f1 <- data.frame(cMeanAgeMother = mean(lifespan.data.dead.f$cMeanAgeMother),
                        cDevMeanAgeMother = unique(sort(lifespan.data.dead.f$cDevMeanAgeMother)), 
                        cMeanAgeFather = mean(lifespan.data.dead.f$cMeanAgeFather),
                        cDevMeanAgeFather = mean(lifespan.data.dead.f$cDevMeanAgeFather),
                        cMeanAgeSocMale = mean(lifespan.data.dead.f$cMeanAgeSocMale),
                        cDevMeanAgeSocMale = mean(lifespan.data.dead.f$cDevMeanAgeSocMale),
                        cBirthYear = mean(lifespan.data.dead.f$cBirthYear),
                        clogTQ = mean(lifespan.data.dead.f$clogTQ),
                        cGroupSize = mean(lifespan.data.dead.f$cGroupSize),
                        HelpersF = '0',
                        Siblings='0',  
                        MotherID=NA,
                        FatherID=NA,
                        DominantMaleID=NA,
                        BirthYearF=NA)


#~~ add to dataframe delta age mother in years (not z-transformed) for plotting

DevAgeMothercorrect2 <- lifespan.data.dead.f %>%
  dplyr::select(DevMeanAgeMother, cDevMeanAgeMother) %>%
  unique()

predDF.f1 <- merge(predDF.f1, DevAgeMothercorrect2, by="cDevMeanAgeMother")


#~~ add model predictions and standard errors to the predict dataframe

predLRSf1 <- predict(lifespan.model.f.wb, predDF.f1, type='response', re.form=NA, se.fit=T)

predDF.f1$pred<-predLRSf1$fit
predDF.f1$upper<-predLRSf1$fit+predLRSf1$se.fit
predDF.f1$lower<-predLRSf1$fit-predLRSf1$se.fit


#~~ graph 

lifespanfgraph <- ggplot(lifespan.data.dead.f, aes(x=DevMeanAgeMother, y=LifespanInt)) +  
  geom_jitter(shape=19, alpha=0.5, colour="grey4", width=0, height=0.15) +
  xlab("Delta age mother (years)") + ylab("Female offspring lifespan") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=13)) +
  theme(axis.title.y=element_text(colour="grey6", size=13)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=13, colour="grey26")) + 
  geom_line(data=predDF.f1, aes(x=DevMeanAgeMother, y=pred), colour="darkorange", size=0.8) + 
  geom_line(data=predDF.f1, aes(x=DevMeanAgeMother, y=upper), colour="darkorange", linetype="longdash", size=0.4) +
  geom_line(data=predDF.f1, aes(x=DevMeanAgeMother, y=lower), colour="darkorange", linetype="longdash", size=0.4) +
  scale_x_continuous(breaks=c(-6,-4,-2,0,2,4,6)) +
  scale_y_continuous(breaks=c(0,4,8,12,16,20)) 

lifespanfgraph

tiff("Graphs/FigureS3_LifespanFemWnMatAgePlot.tif", width=2100, height=1800, res=400)

lifespanfgraph

dev.off()


# # -----------------------------------------------------------------------

# Figure 2 - maternal age and offspring LRS
# Figure S4 - within and between individual parental age effects on offspring LRS

# # -----------------------------------------------------------------------

#~~ data subset (see main LRS script for more info)

LRSdata <- droplevels(subset(fitnessdata, LifespanInt > 0 & !is.na(Sex) & FinalYear < 2018 & TranslocatedF==0)) 

summary(LRSdata)


#~~ first make female graphs: figure 2a and figure s4a

#~~ female model subset

LRSfdata <- droplevels(subset(LRSdata, Sex==0 & !is.na(GroupSize) & !is.na(AgeMother) & !is.na(AgeFather) & !is.na(DominantMaleAge) 
                              & !is.na(BirthYear) & !is.na(Siblings) & !is.na(HelpersF)))


#~~ z-transform covariates

LRSfdata$cAgeMother <- rescale(LRSfdata$AgeMother)
LRSfdata$cAgeFather <- rescale(LRSfdata$AgeFather)
LRSfdata$cDominantMaleAge <- rescale(LRSfdata$DominantMaleAge)
LRSfdata$cBirthYear <- rescale(LRSfdata$BirthYear)
LRSfdata$clogTQ <- rescale(LRSfdata$logTQ)
LRSfdata$cGroupSize <- rescale(LRSfdata$GroupSize)


#~~ model with cross-sectional parental age effects on female offspring LRS

lrs.model.f.1 <- glmmTMB(LRS1YEAR~ cAgeMother + cAgeFather + cDominantMaleAge +
                           cBirthYear + clogTQ + cGroupSize + HelpersF +
                           Siblings + (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                         data=LRSfdata,
                         ziformula=~1,
                         family=poisson)
summary(lrs.model.f.1)


#~~ make dataframe - for cross-sectional maternal age effects on female offspring LRS
#~~ all continuous variables set to mean values except for maternal age
#~~ estimated for offspring with no helpers or siblings in the natal nest

predDF.f <- data.frame(cAgeFather = mean(LRSfdata$cAgeFather),
                        cDominantMaleAge = mean(LRSfdata$cDominantMaleAge), 
                        cBirthYear = mean(LRSfdata$cBirthYear),
                        clogTQ = mean(LRSfdata$clogTQ),
                        cGroupSize = mean(LRSfdata$cGroupSize),
                        HelpersF = '0',
                        Siblings='0',  
                        cAgeMother = unique(sort(LRSfdata$cAgeMother)),
                        MotherID=NA,
                        FatherID=NA,
                        DominantMaleID=NA,
                        BirthYearF=NA)


#~~ add to dataframe maternal age in years (not z-transformed) for plotting

AgeMothercorrect1 <- LRSfdata %>%
  dplyr::select(AgeMother, cAgeMother) %>%
  unique()

predDF.f <- merge(predDF.f, AgeMothercorrect1, by="cAgeMother")


#~~ add model predictions and standard errors to the predict dataframe

predLRSf <- predict(lrs.model.f.1, predDF.f, type='response', re.form=NA, se.fit=T)

predDF.f$pred<-predLRSf$fit
predDF.f$upper<-predLRSf$fit+predLRSf$se.fit
predDF.f$lower<-predLRSf$fit-predLRSf$se.fit


#~~ graph: Figure 2A

fig.2a <- ggplot(LRSfdata, aes(x=AgeMother, y=LRS1YEAR)) +    
  geom_jitter(shape=19, alpha=0.5, colour="grey4", width=0, height=0.15) +
  xlab("Maternal age (years)") + ylab("Female offspring LRS") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=13)) +
  theme(axis.title.y=element_text(colour="grey6", size=13)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=13, colour="grey26")) + 
  geom_line(data=predDF.f, aes(x=AgeMother, y=pred), colour="darkorange", size=0.8) + 
  geom_line(data=predDF.f, aes(x=AgeMother, y=upper), colour="darkorange", linetype="longdash", size=0.4) +
  geom_line(data=predDF.f, aes(x=AgeMother, y=lower), colour="darkorange", linetype="longdash", size=0.4) +
  scale_x_continuous(breaks=c(1,3,5,7,9,11,13)) +
  scale_y_continuous(limits=c(0,10), breaks=c(0,2,4,6,8,10)) 

fig.2a


#~~ make mean maternal, paternal and dominant male age variables and deviation of each value from the mean

LRSfdata<-LRSfdata %>% dplyr::group_by(MotherID) %>% dplyr::mutate(cMeanAgeMother = mean(cAgeMother)) %>% dplyr::ungroup()
LRSfdata$cDevMeanAgeMother<-LRSfdata$cAgeMother-LRSfdata$cMeanAgeMother

LRSfdata<-LRSfdata %>% dplyr::group_by(FatherID) %>% dplyr::mutate(cMeanAgeFather = mean(cAgeFather)) %>% dplyr::ungroup()
LRSfdata$cDevMeanAgeFather<-LRSfdata$cAgeFather-LRSfdata$cMeanAgeFather

LRSfdata<-LRSfdata %>% dplyr::group_by(DominantMaleID) %>% dplyr::mutate(cMeanAgeSocMale = mean(cDominantMaleAge)) %>% dplyr::ungroup()
LRSfdata$cDevMeanAgeSocMale<-LRSfdata$cDominantMaleAge-LRSfdata$cMeanAgeSocMale

LRSfdata<-LRSfdata %>% dplyr::group_by(MotherID) %>% dplyr::mutate(MeanAgeMother = mean(AgeMother)) %>% dplyr::ungroup()
LRSfdata$DevMeanAgeMother<-LRSfdata$AgeMother-LRSfdata$MeanAgeMother


#~~ within-subject centering model for females

lrs.model.f.wb <- glmmTMB(LRS1YEAR~cMeanAgeMother + cDevMeanAgeMother +
                            cMeanAgeFather + cDevMeanAgeFather + 
                            cMeanAgeSocMale + cDevMeanAgeSocMale +
                            cBirthYear + clogTQ + cGroupSize + HelpersF + Siblings +
                            (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                          data=LRSfdata,
                          ziformula=~1,
                          family=poisson)

summary(lrs.model.f.wb)


#~~ make dataframe - for within-individual maternal age effects on female offspring LRS
#~~ all continuous variables set to mean values except for within-individual maternal age effects
#~~ estimated for offspring with no helpers or siblings in the natal nest

predDF.f2 <- data.frame(cMeanAgeMother = mean(LRSfdata$cMeanAgeMother),
                        cDevMeanAgeMother = unique(sort(LRSfdata$cDevMeanAgeMother)), 
                        cMeanAgeFather = mean(LRSfdata$cMeanAgeFather),
                        cDevMeanAgeFather = mean(LRSfdata$cDevMeanAgeFather),
                        cMeanAgeSocMale = mean(LRSfdata$cMeanAgeSocMale),
                        cDevMeanAgeSocMale = mean(LRSfdata$cDevMeanAgeSocMale),
                       cBirthYear = mean(LRSfdata$cBirthYear),
                       clogTQ = mean(LRSfdata$clogTQ),
                       cGroupSize = mean(LRSfdata$cGroupSize),
                       HelpersF = '0',
                       Siblings='0',  
                       MotherID=NA,
                       FatherID=NA,
                       DominantMaleID=NA,
                       BirthYearF=NA)


#~~ add to dataframe within-individual maternal age in years (not z-transformed) for plotting

DevAgeMothercorrect <- LRSfdata %>%
  dplyr::select(DevMeanAgeMother, cDevMeanAgeMother) %>%
  unique()

predDF.f2 <- merge(predDF.f2, DevAgeMothercorrect, by="cDevMeanAgeMother")


#~~ add model predictions and standard errors to the predict dataframe

predLRSf2 <- predict(lrs.model.f.wb, predDF.f2, type='response', re.form=NA, se.fit=T)

predDF.f2$pred<-predLRSf2$fit
predDF.f2$upper<-predLRSf2$fit+predLRSf2$se.fit
predDF.f2$lower<-predLRSf2$fit-predLRSf2$se.fit


#~~ graph: Figure S4A

fig.s4a <- ggplot(LRSfdata, aes(x=DevMeanAgeMother, y=LRS1YEAR)) +  
  geom_jitter(shape=19, alpha=0.5, colour="grey4", width=0, height=0.15) +
  xlab("Delta age mother (years)") + ylab("Female offspring LRS") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=13)) +
  theme(axis.title.y=element_text(colour="grey6", size=13)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=13, colour="grey26")) + 
  geom_line(data=predDF.f2, aes(x=DevMeanAgeMother, y=pred), colour="darkorange", size=0.8) + 
  geom_line(data=predDF.f2, aes(x=DevMeanAgeMother, y=upper), colour="darkorange", linetype="longdash", size=0.4) +
  geom_line(data=predDF.f2, aes(x=DevMeanAgeMother, y=lower), colour="darkorange", linetype="longdash", size=0.4) +
  scale_x_continuous(breaks=c(-6,-4,-2,0,2,4,6)) +
  scale_y_continuous(breaks=c(0,2,4,6,8,10)) 

fig.s4a


#~~ now for the male graphs: figure 2b and figure s4b and figure s4c

#~~ male model subset

LRSmdata <- droplevels(subset(LRSdata, Sex==1 & !is.na(GroupSize) & !is.na(AgeMother) & !is.na(AgeFather) & !is.na(DominantMaleAge) 
                              & !is.na(BirthYear) & !is.na(Siblings) & !is.na(HelpersF)))

#~~ z-transform covariates

LRSmdata$cAgeMother <- rescale(LRSmdata$AgeMother)
LRSmdata$cAgeFather <- rescale(LRSmdata$AgeFather)
LRSmdata$cDominantMaleAge <- rescale(LRSmdata$DominantMaleAge)
LRSmdata$cBirthYear <- rescale(LRSmdata$BirthYear)
LRSmdata$clogTQ <- rescale(LRSmdata$logTQ)
LRSmdata$cGroupSize <- rescale(LRSmdata$GroupSize)


#~~ model with cross-sectional parental age effects on male offspring LRS

lrs.model.m.1 <- glmmTMB(LRS1YEAR~ cAgeMother + cAgeFather + cDominantMaleAge +
                           cBirthYear + clogTQ + cGroupSize + HelpersF +
                           Siblings + (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                         data=LRSmdata,
                         ziformula=~1,
                         family=poisson)
summary(lrs.model.m.1)


#~~ make dataframe - for cross-sectional maternal age effects on male offspring LRS
#~~ all continuous variables set to mean values except for maternal age
#~~ estimated for offspring with no helpers or siblings in the natal nest

predDF.m <- data.frame(cAgeFather = mean(LRSmdata$cAgeFather),
                        cDominantMaleAge = mean(LRSmdata$cDominantMaleAge), 
                        cBirthYear = mean(LRSmdata$cBirthYear),
                        clogTQ = mean(LRSmdata$clogTQ),
                        cGroupSize = mean(LRSmdata$cGroupSize),
                        HelpersF = '0',
                        Siblings='0',  
                        cAgeMother = unique(sort(LRSmdata$cAgeMother)),
                        MotherID=NA,
                        FatherID=NA,
                        DominantMaleID=NA,
                        BirthYearF=NA)


#~~ add to dataframe maternal age in years (not z-transformed) for plotting

AgeMothercorrect2 <- LRSmdata %>%
  dplyr::select(AgeMother, cAgeMother) %>%
  unique()

predDF.m <- merge(predDF.m, AgeMothercorrect2, by="cAgeMother")


#~~ add model predictions and standard errors to the predict dataframe

predLRSm <- predict(lrs.model.m.1, predDF.m, type='response', re.form=NA, se.fit=T)

predDF.m$pred<-predLRSm$fit
predDF.m$upper<-predLRSm$fit+predLRSm$se.fit
predDF.m$lower<-predLRSm$fit-predLRSm$se.fit


#~~ graph: Figure 2B

fig.2b <- ggplot(LRSmdata, aes(x=AgeMother, y=LRS1YEAR)) +  
  geom_jitter(shape=19, alpha=0.5, colour="grey4", width=0, height=0.15) +
  xlab("Maternal age (years)") + ylab("Male offspring LRS") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=13)) +
  theme(axis.title.y=element_text(colour="grey6", size=13)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=13, colour="grey26")) + 
  geom_line(data=predDF.m, aes(x=AgeMother, y=pred), colour="darkorange", size=0.8) + 
  geom_line(data=predDF.m, aes(x=AgeMother, y=upper), colour="darkorange", linetype="longdash", size=0.4) +
  geom_line(data=predDF.m, aes(x=AgeMother, y=lower), colour="darkorange", linetype="longdash", size=0.4) +
  scale_x_continuous(breaks=c(1,3,5,7,9,11,13)) +
  scale_y_continuous(breaks=c(0,2,4,6,8,10)) 

fig.2b


#~~ make mean maternal, paternal and dominant male age variables and deviation of each value from the mean

LRSmdata<-LRSmdata %>% dplyr::group_by(MotherID) %>% dplyr::mutate(cMeanAgeMother = mean(cAgeMother)) %>% dplyr::ungroup()
LRSmdata$cDevMeanAgeMother<-LRSmdata$cAgeMother-LRSmdata$cMeanAgeMother

LRSmdata<-LRSmdata %>% dplyr::group_by(FatherID) %>% dplyr::mutate(cMeanAgeFather = mean(cAgeFather)) %>% dplyr::ungroup()
LRSmdata$cDevMeanAgeFather<-LRSmdata$cAgeFather-LRSmdata$cMeanAgeFather

LRSmdata<-LRSmdata %>% dplyr::group_by(DominantMaleID) %>% dplyr::mutate(cMeanAgeSocMale = mean(cDominantMaleAge)) %>% dplyr::ungroup()
LRSmdata$cDevMeanAgeSocMale<-LRSmdata$cDominantMaleAge-LRSmdata$cMeanAgeSocMale

LRSmdata<-LRSmdata %>% dplyr::group_by(MotherID) %>% dplyr::mutate(MeanAgeMother = mean(AgeMother)) %>% dplyr::ungroup()

LRSmdata<-LRSmdata %>% dplyr::group_by(FatherID) %>% dplyr::mutate(MeanAgeFather = mean(AgeFather)) %>% dplyr::ungroup()


#~~ within-subject centering model for males

lrs.model.m.wb <- glmmTMB(LRS1YEAR~cMeanAgeMother + cDevMeanAgeMother +
                            cMeanAgeFather + cDevMeanAgeFather + 
                            cMeanAgeSocMale + cDevMeanAgeSocMale +
                            cBirthYear + clogTQ + cGroupSize + HelpersF + Siblings +
                            (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                          data=LRSmdata,
                          ziformula=~1,
                          family=poisson)

summary(lrs.model.m.wb)


#~~ make dataframe - for between-individual maternal age effects on male offspring LRS
#~~ all continuous variables set to mean values except for between-individual maternal age effects
#~~ estimated for offspring with no helpers or siblings in the natal nest

predDF.m2 <- data.frame(cMeanAgeMother = unique(sort(LRSmdata$cMeanAgeMother)),
                        cDevMeanAgeMother = mean(LRSmdata$cDevMeanAgeMother), 
                        cMeanAgeFather = mean(LRSmdata$cMeanAgeFather),
                        cDevMeanAgeFather = mean(LRSmdata$cDevMeanAgeFather),
                        cMeanAgeSocMale = mean(LRSmdata$cMeanAgeSocMale),
                        cDevMeanAgeSocMale = mean(LRSmdata$cDevMeanAgeSocMale),
                        cBirthYear = mean(LRSmdata$cBirthYear),
                        clogTQ = mean(LRSmdata$clogTQ),
                        cGroupSize = mean(LRSmdata$cGroupSize),
                        HelpersF = '0',
                        Siblings='0',  
                        MotherID=NA,
                        FatherID=NA,
                        DominantMaleID=NA,
                        BirthYearF=NA)


#~~ add to dataframe between-individual maternal age in years (not z-transformed) for plotting

MeanAgeMothercorrect2 <- LRSmdata %>%
  dplyr::select(MeanAgeMother, cMeanAgeMother) %>%
  unique()

predDF.m2 <- merge(predDF.m2, MeanAgeMothercorrect2, by="cMeanAgeMother")


#~~ add model predictions and standard errors to the predict dataframe

predLRSm2 <- predict(lrs.model.m.wb, predDF.m2, type='response', re.form=NA, se.fit=T)

predDF.m2$pred<-predLRSm2$fit
predDF.m2$upper<-predLRSm2$fit+predLRSm2$se.fit
predDF.m2$lower<-predLRSm2$fit-predLRSm2$se.fit


#~~ graph: Figure S4B

fig.s4b <- ggplot(LRSmdata, aes(x=MeanAgeMother, y=LRS1YEAR)) +  
  geom_jitter(shape=19, alpha=0.5, colour="grey4", width=0, height=0.15) +
  xlab("Mean maternal age (years)") + ylab("Male offspring LRS") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=13)) +
  theme(axis.title.y=element_text(colour="grey6", size=13)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=13, colour="grey26")) + 
  geom_line(data=predDF.m2, aes(x=MeanAgeMother, y=pred), colour="darkorange", size=0.8) + 
  geom_line(data=predDF.m2, aes(x=MeanAgeMother, y=upper), colour="darkorange", linetype="longdash", size=0.4) +
  geom_line(data=predDF.m2, aes(x=MeanAgeMother, y=lower), colour="darkorange", linetype="longdash", size=0.4) +
  scale_x_continuous(breaks=c(1,3,5,7,9,11,13)) +
  scale_y_continuous(breaks=c(0,2,4,6,8,10)) 

fig.s4b


#~~ make dataframe - for between-individual paternal age effects on male offspring LRS
#~~ all continuous variables set to mean values except for between-individual paternal age effects
#~~ estimated for offspring with no helpers or siblings in the natal nest


predDF.m3 <- data.frame(cMeanAgeMother = mean(LRSmdata$cMeanAgeMother),
                        cDevMeanAgeMother = mean(LRSmdata$cDevMeanAgeMother), 
                        cMeanAgeFather = unique(sort(LRSmdata$cMeanAgeFather)),
                        cDevMeanAgeFather = mean(LRSmdata$cDevMeanAgeFather),
                        cMeanAgeSocMale = mean(LRSmdata$cMeanAgeSocMale),
                        cDevMeanAgeSocMale = mean(LRSmdata$cDevMeanAgeSocMale),
                        cBirthYear = mean(LRSmdata$cBirthYear),
                        clogTQ = mean(LRSmdata$clogTQ),
                        cGroupSize = mean(LRSmdata$cGroupSize),
                        HelpersF = '0',
                        Siblings='0',  
                        MotherID=NA,
                        FatherID=NA,
                        DominantMaleID=NA,
                        BirthYearF=NA)


#~~ add to dataframe between-individual paternal age in years (not z-transformed) for plotting

MeanAgeFathercorrect <- LRSmdata %>%
  dplyr::select(MeanAgeFather, cMeanAgeFather) %>%
  unique()

predDF.m3 <- merge(predDF.m3, MeanAgeFathercorrect, by="cMeanAgeFather")


#~~ add model predictions and standard errors to the predict dataframe

predLRSm3 <- predict(lrs.model.m.wb, predDF.m3, type='response', re.form=NA, se.fit=T)

predDF.m3$pred<-predLRSm3$fit
predDF.m3$upper<-predLRSm3$fit+predLRSm3$se.fit
predDF.m3$lower<-predLRSm3$fit-predLRSm3$se.fit


#~~ graph: Figure S4C

fig.s4c <- ggplot(LRSmdata, aes(x=MeanAgeFather, y=LRS1YEAR)) +  
  geom_jitter(shape=19, alpha=0.5, colour="grey4", width=0, height=0.15) +
  xlab("Mean paternal age (years)") + ylab("Male offspring LRS") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=13)) +
  theme(axis.title.y=element_text(colour="grey6", size=13)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=13, colour="grey26")) + 
  geom_line(data=predDF.m3, aes(x=MeanAgeFather, y=pred), colour="darkorange", size=0.8) + 
  geom_line(data=predDF.m3, aes(x=MeanAgeFather, y=upper), colour="darkorange", linetype="longdash", size=0.4) +
  geom_line(data=predDF.m3, aes(x=MeanAgeFather, y=lower), colour="darkorange", linetype="longdash", size=0.4) +
  scale_x_continuous(breaks=c(1,3,5,7,9,11,13)) +
  scale_y_continuous(breaks=c(0,2,4,6,8,10)) 

fig.s4c

#~~ compile figures for figure 2

fig.2 <- plot_grid(fig.2a, fig.2b, ncol=2, nrow=1, labels="AUTO")

pdf("Graphs/Figure2_LRSGraphs.pdf", width=8.4, height=4, onefile=FALSE)

print(fig.2)

dev.off()


#~~ compile figures for figure s4

fig.s4 <- plot_grid(fig.s4a, fig.s4b, fig.s4c, labels="AUTO", ncol=3, nrow=1)

tiff("Graphs/FigureS4_LRSwnbnGraphs.tif", width=4200, height=1500, res=400)

fig.s4

dev.off()


# # -----------------------------------------------------------------------

# Figure S6-S9: TQ*paternal age effects on female LRS

# # -----------------------------------------------------------------------

#~~ run the model on non rescaled variables to help visualise

#~~ recode birth year to help model converge

LRSfdata$BYRecoded <- LRSfdata$BirthYear-1994

LRSfdatacheck <- subset(LRSfdata, AgeFather > 12)


#~~ model

lrs.model.int <- glmmTMB(LRS1YEAR~ AgeMother + AgeFather + DominantMaleAge +
                           logTQ + GroupSize + HelpersF + Siblings + BYRecoded + AgeFather*logTQ +
                           (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                         data=LRSfdata,
                         ziformula=~1,
                         family=poisson)


summary(lrs.model.int)


#~~ figure s6

fig.s6 <- interact_plot(lrs.model.int, pred=AgeFather, modx=logTQ,
              data=LRSfdata, interval=TRUE, 
              outcome.scale="response", int.type="confidence", 
              int.width = 0.95, plot.points = TRUE, x.label="Age of genetic father",
              y.label="Female offspring LRS", legend.main="Territory quality") 

fig.s6

tiff("Graphs/FigureS6_TQ_parage_interaction_SD.tif", width=2500, height=1500, res=400)

fig.s6

dev.off()


#~~ figure s7

fig.s7 <- interact_plot(lrs.model.int, pred=AgeFather, modx=logTQ, modx.values="terciles",
              data=LRSfdata, interval=TRUE, 
              outcome.scale="response", int.type="confidence", 
              int.width = 0.95, plot.points = TRUE, facet.modx=TRUE, x.label="Age of genetic father",
              y.label="Female offspring LRS") 


fig.s7

tiff("Graphs/FigureS7_TQ_parage_interaction_terciles.tif", width=3000, height=1500, res=400)

fig.s7

dev.off()


#~~ check dominant male effects are similar

lrs.model.int2 <- glmmTMB(LRS1YEAR~ AgeMother + AgeFather + DominantMaleAge +
                           logTQ + GroupSize + HelpersF + Siblings + BYRecoded + DominantMaleAge*logTQ +
                           (1|MotherID) + (1|FatherID) + (1|DominantMaleID) + (1|BirthYearF),
                         data=LRSfdata,
                         ziformula=~1,
                         family=poisson)


summary(lrs.model.int2)


#~~ figure s8 

fig.s8 <- interact_plot(lrs.model.int2, pred=DominantMaleAge, modx=logTQ,
              data=LRSfdata, interval=TRUE, 
              outcome.scale="response", int.type="confidence",
              int.width = 0.95, plot.points = TRUE, x.label="Age of dominant male",
              y.label="Female offspring LRS", legend.main="Territory quality") 

tiff("Graphs/FigureS8_TQ_socmaleage_interaction_SD.tif", width=2500, height=1500, res=400)

fig.s8

dev.off()

LRSfdataolderdominantmales <- subset(LRSfdata, DominantMaleAge>10)


#~~ figure s9

fig.s9 <- interact_plot(lrs.model.int2, pred=DominantMaleAge, modx=logTQ, modx.values="terciles",
              data=LRSfdata, interval=TRUE, 
              outcome.scale="response", int.type="confidence", 
              int.width = 0.95, plot.points = TRUE, facet.modx=TRUE, x.label="Age of dominant male",
              y.label="Female offspring LRS") 

tiff("Graphs/FigureS9_TQ_socmaleage_interaction_terciles.tif", width=3000, height=1500, res=400)

fig.s9

dev.off()