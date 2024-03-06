########################################################################
### R script for statistical analysis of barley grain isotope values ###
###   Hald et al. 2024 Journal of Archaeological Science: Reports    ###
### Title: "Farming during turbulent times: Agriculture, food crops, ### 
###   and manuring practices in Bronze Age to Viking Age Denmark"    ###
###   This R code runs all analyses presented in the manuscript      ###
###      Authors: A Styring, amy.styring@arch.ox.ac.uk               ###
########################################################################

# load the libraries into R
library(readxl)
library(dplyr)
library(plyr)
library(nlme)
library(car)

# load data
data <- read_excel("Supplementary Table 3.xlsx")
names(data)[names(data) == "Hordeum vulgare"] <- "BarleyType"
names(data)[names(data) == "WeightMean (mg)"] <- "WeightMean"

#####################################################################
###	2.2 Isotope analysis                                        	###
#####################################################################

# Select Iron Age data only
IA<-data[data$Centre_cal2sigma>-500 & data$Centre_cal2sigma<750,]

# Summarise d15N values
min(IA$normd15N)
max(IA$normd15N)
mean(IA$normd15N)
sd(IA$normd15N)
nrow(IA)

# Check normality of d15N values
qqPlot(IA$normd15N) ## There are 2 high d15N values
shapiro.test(IA$normd15N) ## p = 0.0001, so not normally distributed.

# Check normality of d15N values when removing the outliers > 14 per mil
batch <- IA[IA$normd15N < 14,]
nrow(batch)
shapiro.test(batch$normd15N) ## p = 0.092

# Summarise d13C values
min(IA$normd13C)
max(IA$normd13C)
mean(IA$normd13C)
sd(IA$normd13C)

## Check normality of d13C values
qqPlot(IA$normd13C) 
shapiro.test(IA$normd13C) ## p = 0.154

# Summarise D13C values
min(IA$D13C)
max(IA$D13C)
mean(IA$D13C)
sd(IA$D13C)

# Summarise data by context
summarydat <- ddply(data, c("SampleID","`Hordeum vulgare`", "Site"), 
                 function(x) c(d15N=mean(x$normd15N), 
                               sdN=sd(x$normd15N), diffd15N=max(x$normd15N)-min(x$normd15N), d13C=mean(x$normd13C), 
                               sdC=sd(x$normd13C), diffd13C=max(x$normd13C)-min(x$normd13C), Delta13C=mean(x$D13C),
                               DateMidpoint=mean(x$Centre_cal2sigma), n=nrow(x)))

# Summarise d15N values by context
min(summarydat$diffd15N[summarydat$diffd15N!=0])
max(summarydat$diffd15N)
min(summarydat$sdN, na.rm=T)
max(summarydat$sdN, na.rm=T)
mean(summarydat$sdN, na.rm=T)

# Summarise d13C values by context
min(summarydat$diffd13C[summarydat$diffd13C!=0])
max(summarydat$diffd13C)
min(summarydat$sdC, na.rm=T)
max(summarydat$sdC, na.rm=T)
mean(summarydat$sdC, na.rm=T)

##################################################################################
## Perform Linear Mixed Model to test the effect of fixed variables on d15N values

# Select Iron Age data only
IA<-data[data$Centre_cal2sigma>-500 & data$Centre_cal2sigma<750,]

# Remove the two outliers
batch <- IA[IA$normd15N < 14,]

# Select data only from barley grains from contexts for which there is a 
# directly determined radiocarbon age (Relation2isotop is "Primary" or "Secondary")
batch$Relation2isotop[is.na(batch$Relation2isotop)] <- "Unknown"
batch<-batch[batch$Relation2isotop=="Primary" | batch$Relation2isotop=="Secondary",]

# Choose a model by AIC

interceptOnly <- gls(normd15N ~ 1, data=batch, method="ML")

randomInterceptOnly <- lme(normd15N ~ 1, data=batch, random=~1|Site, method="ML")

interceptDate <- gls(normd15N ~ Centre_cal2sigma, data=batch, method="ML")

randomInterceptDate <- lme(normd15N ~ Centre_cal2sigma, data=batch, random=~1|Site, method="ML")
# Lowest AIC so likely to be the best model
# No relationship between date and d15N

randomInterceptDateRegion <- lme(normd15N ~ Centre_cal2sigma + Region, data=batch, random=~1|Site, method="ML")

randomInterceptDateBarley <- lme(normd15N ~ Centre_cal2sigma + BarleyType, data=batch, random=~1|Site, method="ML")

anova(interceptOnly, randomInterceptOnly, interceptDate, randomInterceptDate, randomInterceptDateRegion, randomInterceptDateBarley)

# model found: randomInterceptDate
# normd15N ~ Centre_cal2sigma, data=batch, random=~1|Site
summary(randomInterceptDate)
intervals(randomInterceptDate)

##################################################################################
## Perform Linear Mixed Model to test the effect of fixed variables on D13C values

# Choose a model by AIC

interceptOnly <- gls(D13C ~ 1, data=batch, method="ML")

randomInterceptOnly <- lme(D13C ~ 1, data=batch, random=~1|Site, method="ML")
# Lowest AIC so likely to be the best model

interceptDate <- gls(D13C ~ Centre_cal2sigma, data=batch, method="ML")

randomInterceptDate <- lme(D13C ~ Centre_cal2sigma, data=batch, random=~1|Site, method="ML")
# No relationship between date and D13C

randomInterceptDateRegion <- lme(D13C ~ Centre_cal2sigma + Region, data=batch, random=~1|Site, method="ML")

randomInterceptDateBarley <- lme(D13C ~ Centre_cal2sigma + BarleyType, data=batch, random=~1|Site, method="ML")

anova(interceptOnly, randomInterceptOnly, interceptDate, randomInterceptDate, randomInterceptDateRegion, randomInterceptDateBarley)

# model: randomInterceptDate
# D13C ~ Centre_cal2sigma, data=batch, random=~1|Site
summary(randomInterceptDate)
intervals(randomInterceptDate)

#######################################################################################
## Perform Linear Mixed Model to test the effect of fixed variables on estimated weight

# Select only hulled barley
batch<-batch[batch$BarleyType=="hulled",]
# Select only complete barley grains whose XYZ dimensions could be measured
batch<-batch[batch$WeightMean !="NA",]
batch$WeightMean<-as.numeric(batch$WeightMean)

# Choose a model by AIC

interceptOnly <- gls(WeightMean ~ 1, data=batch, method="ML")
# Lowest AIC so likely to be best model

randomInterceptOnly <- lme(WeightMean ~ 1, data=batch, random=~1|Site, method="ML")

interceptDate <- gls(WeightMean ~ Centre_cal2sigma, data=batch, method="ML")

randomInterceptDate <- lme(WeightMean ~ Centre_cal2sigma, data=batch, random=~1|Site, method="ML")
summary(randomInterceptDate) ## No relationship between date and estimate weight

randomInterceptDateRegion <- lme(WeightMean ~ Centre_cal2sigma + Region, data=batch, random=~1|Site, method="ML")

anova(interceptOnly, randomInterceptOnly, interceptDate, randomInterceptDate, randomInterceptDateRegion)

# model: randomInterceptDate
# WeightMean ~ Centre_cal2sigma, data=batch, random=~1|Site
summary(randomInterceptDate)
intervals(randomInterceptDate)

#####################################################################
###	4.2 Soil enrichment practices                                 ###
#####################################################################

# Calculating mean d15N value of wild herbivore bone collagen from Denmark
# load data
DIANA_Fauna <- read_excel("Supplementary Table 4.xlsx")

wild<-DIANA_Fauna[DIANA_Fauna$TYPE=="WildHerbivore",]
wild<-wild[complete.cases(wild$TYPE),]

deer<-wild[wild$TAXON=="Red deer" | wild$TAXON=="Roe deer",]
names(deer)[names(deer) == "SI_d15N"] <- "d15N"
names(deer)[names(deer) == "SI_d13C"] <- "d13C"
names(deer)[names(deer) == "TAXON_LATIN"] <- "Species"
deer<-deer[,c(14,2,4)]

# load data
MaringRiede<-read_excel("Maring_Riede_2019.xlsx")
deer2<-MaringRiede[MaringRiede$Species=="Cervus elaphus" | MaringRiede$Species=="Capreolus capreolus",]
deer2<-deer2[,c(3,8,9)]

dat<-rbind(deer,deer2)

# Summarise deer bone collagen d15N values
mean(dat$d15N)
model<-lm(d15N ~ 1, dat)
confint(model)
max(dat$d15N)
min(dat$d15N)
