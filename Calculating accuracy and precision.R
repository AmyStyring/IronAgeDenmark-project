#######################################################################
###     R script for calculating accuracy and precision of data     ### 
###                 in Supplementary Table 1                        ###
###  Hald et al. 2024 Journal of Archaeological Science: Reports    ###
###    Title: "Farming during turbulent times: Agriculture, food    ### 
###          crops, and manuring practices in Bronze Age to         ###
###                       Viking Age Denmark"                       ###
###   This R code runs all analyses presented in the manuscript     ###
###       Authors: A Styring, amy.styring@arch.ox.ac.uk             ###
#######################################################################

# load the libraries into R
library(readxl)
library(dplyr)
library(multiway)

###################################################################
####				EXTRACTING REFERENCE MATERIALS FROM RUNFILES				###
###################################################################

# load data
data1 <- read_excel("Supplementary Table 1.xlsx", sheet = "Data-Session 1")
data2 <- read_excel("Supplementary Table 1.xlsx", sheet = "Data-Session 2")
data3 <- read_excel("Supplementary Table 1.xlsx", sheet = "Data-Session 3")
data4 <- read_excel("Supplementary Table 1.xlsx", sheet = "Data-Session 4")
data5 <- read_excel("Supplementary Table 1.xlsx", sheet = "Data-Session 5")
data6 <- read_excel("Supplementary Table 1.xlsx", sheet = "Data-Session 6 (N only)")
data7 <- read_excel("Supplementary Table 1.xlsx", sheet = "Data-Session 7 (C only)")
data8 <- read_excel("Supplementary Table 1.xlsx", sheet = "Data-Session 8")

data<-rbind(data1, data2, data3, data4, data5, data6, data7, data8)
RawStandards <- data[grep("LEU|COW|SEAL|P2|SALANINE", data$ID), ]
RawStandards <- RawStandards[ ,c(2,15,19,20)]
names(RawStandards) <- c("ID", "Runfile","normd13C","normd15N")

###################################################################
####					 EXTRACTING REPLICATE SAMPLE MEASUREMENTS					###
###################################################################
data$Runfile <- as.character(data$`Runfile date`)

# Replicate samples are denoted by DA and DB after the Sample ID
RepCA<-data[grep("DA$", data$ID), ]
RepCA <- RepCA[order(as.character(RepCA$ID)), ]
RepCA <- RepCA[ ,c(2,15, 19, 20)]
names(RepCA) <- c("ID", "Runfile","normd13C_DuplA","normd15N_DuplA")
RepCA$ID<-gsub('DA$','',RepCA$ID)

RepCB<-data[grep("DB$", data$ID), ]
RepCB <- RepCB[order(as.character(RepCB$ID)), ]
RepCB <- RepCB[ ,c(2,19, 20)]
names(RepCB) <- c("ID", "normd13C_DuplB","normd15N_DuplB")
RepCB$ID<-gsub('DB$','',RepCB$ID)

RepCar<-merge(RepCA, RepCB, "ID")

###################   Difference between modern duplicates     #################

RepCar$Diffd13C<-RepCar$normd13C_DuplA-RepCar$normd13C_DuplB
RepCar$Diffd13C<-abs(RepCar$Diffd13C)

RepCar$Diffd15N<-RepCar$normd15N_DuplA-RepCar$normd15N_DuplB
RepCar$Diffd15N<-abs(RepCar$Diffd15N)

################################################################################        
##                    CALCULATING ACCURACY AND PRECISION                      ##
##                                  CARBON 						                        ##
################################################################################

################### Mean and Stdev for all analytical sessions #################

Cmean<-aggregate(RawStandards$normd13C, list(RawStandards$ID), mean)
Csd<-aggregate(RawStandards$normd13C, list(RawStandards$ID), sd)
N<-count(RawStandards, "ID")

all.standards<-cbind(N,Cmean, Csd)
all.standards<-all.standards[,c(3,2,4,6)]
names(all.standards)<-c("RM","Number", "d13Cmean","d13Csd")
all.standards$srm<-(all.standards$Number-1)*(all.standards$d13Csd^2)

################################################################################        
##                    Check and calibration standards                         ##
##          Pooled standard deviation of each standard (Ssrm) 						    ##
##                and the degrees of freedom (dfsrm)				   						    ##
################################################################################
dfsrm<-sum(all.standards$Number)-nrow(all.standards)
Ssrm<-sqrt(sum(all.standards$srm)/dfsrm)

################################################################################        
##                            Check standards                                 ##
################################################################################
checkC<-subset(all.standards,all.standards$RM=="P2"|all.standards$RM=="LEU")

CheckS.1<--28.19#### P2
CheckS.2<--28.23 #### LEU

y="P2" #### change if using different check standards - this should reflect the name of CheckS.1
fun1<-function(x,y) if(x==y) {CheckS.1} else {CheckS.2}
checkC$known<-mapply(fun1, checkC$RM, y)
CheckS.1sd<-0.14 #### P2
CheckS.2sd<-0.07 #### LEU

fun1<-function(x,y) if(x==y) {CheckS.1sd} else {CheckS.2sd}
checkC$knownsd<-mapply(fun1, checkC$RM, y)

checkC$Diff_measured_known<-checkC$d13Cmean-checkC$known
################################################################################
# RMS bias = the root mean square of the difference between the observed mean	 #
# and the known values of standard reference materials which are treated as    #
# unknowns during analysis (Szpak's "check standards").                 			 #
# u_cref = the root mean square of the known standard deviations of the RMs    #
# used as check standards.                                                		 # 
################################################################################

RMSbias<-sqrt(sumsq(checkC$Diff_measured_known)/nrow(checkC))
u_cref<-sqrt(sumsq(checkC$knownsd)/nrow(checkC))

################################################################################		 
##                                ACCURACY	                                  ##
################################################################################

x<-list(RMSbias, u_cref)
u_bias<-sqrt(sumsq(x))

################################################################################		 
##                               SAMPLE REPLICATES	                          ##
################################################################################
RepCar$Sd<-apply(subset(RepCar,select = c("normd13C_DuplA","normd13C_DuplB")),1,sd)

RepCar$Mean<-apply(subset(RepCar,select = c("normd13C_DuplA","normd13C_DuplB")),1,mean)

RepCar$Number<-2## this needs to change if you run a replicate more than twice in a run

RepCar$RepSsrm<-1*(RepCar$Sd^2)
dfrep<-sum(RepCar$Number)-nrow(RepCar)
Srep<-sqrt((sum(RepCar$RepSsrm))/dfrep)

################################################################################
##                                 PRECISION                                  ##
## Random errors within the laboratory = precision ( Menditto et al 2007)     ##
## or repeatability (Carter and Fry 2013) - summed standard deviation of      ##
## all repeated measurements during relevant analytical sessions - including  ##
##  check and calibration RMs (Ssrum) and the replicates (Srep)	              ##
################################################################################

uRw<-sqrt((Ssrm^2)+(Srep^2)/2)

################################################################################		 
##                      Standard uncertainty (Szpak 2017)	                    ##
################################################################################

y<-list(u_bias, uRw)
Uc<-sqrt(sumsq(y))

################################################################################		 
##                      VALUES TO REPORT (cf. Szpak 2017)                     ##
################################################################################
print(uRw) ## precision (u(Rw))
print(u_bias) ## accuracy, or systematic error (u(bias))
print(Uc) ## total analytical uncertainty (Uc)

################################################################################        
##                    CALCULATING ACCURACY AND PRECISION                      ##
##                                NITROGEN 						                        ##
################################################################################

################### Mean and Stdev for all analytical sessions #################

Nmean<-aggregate(RawStandards$normd15N, list(RawStandards$ID), mean)
Nsd<-aggregate(RawStandards$normd15N, list(RawStandards$ID), sd)
N<-count(RawStandards, "ID")

all.standards<-cbind(N,Nmean, Nsd)
all.standards<-all.standards[,c(3,2,4,6)]
names(all.standards)<-c("RM","Number", "d15Nmean","d15Nsd")
all.standards$srm<-(all.standards$Number-1)*(all.standards$d15Nsd^2)

################################################################################        
##                    Check and calibration standards                         ##
##          Pooled standard deviation of each standard (Ssrm) 						    ##
##                and the degrees of freedom (dfsrm)				   						    ##
################################################################################
dfsrm<-sum(all.standards$Number)-nrow(all.standards)
Ssrm<-sqrt(sum(all.standards$srm)/dfsrm)

################################################################################        
##                            Check standards                                 ##
################################################################################
checkN<-subset(all.standards,all.standards$RM=="P2"|all.standards$RM=="LEU")

CheckS.1<--1.57#### P2
CheckS.2<-6.35 #### LEU

y="P2" #### this should reflect the name of CheckS.1
fun1<-function(x,y) if(x==y) {CheckS.1} else {CheckS.2}
checkN$known<-mapply(fun1, checkC$RM, y)
CheckS.1sd<-0.19 #### P2
CheckS.2sd<-0.14 #### LEU

fun1<-function(x,y) if(x==y) {CheckS.1sd} else {CheckS.2sd}
checkN$knownsd<-mapply(fun1, checkN$RM, y)

checkN$Diff_measured_known<-checkN$d15Nmean-checkN$known
################################################################################
# RMS bias = the root mean square of the difference between the observed mean	 #
# and the known values of standard reference materials which are treated as    #
# unknowns during analysis (Szpak's "check standards").                 			 #
# u_cref = the root mean square of the known standard deviations of the RMs    #
# used as check standards.                                                		 # 
################################################################################

RMSbias<-sqrt(sumsq(checkN$Diff_measured_known)/nrow(checkN))
u_cref<-sqrt(sumsq(checkN$knownsd)/nrow(checkN))

################################################################################		 
##                                ACCURACY	                                  ##
################################################################################

x<-list(RMSbias, u_cref)
u_bias<-sqrt(sumsq(x))

################################################################################		 
##                               SAMPLE REPLICATES	                          ##
################################################################################
RepCar$Sd<-apply(subset(RepCar,select = c("normd15N_DuplA","normd15N_DuplB")),1,sd)

RepCar$Mean<-apply(subset(RepCar,select = c("normd15N_DuplA","normd15N_DuplB")),1,mean)

RepCar$Number<-2## this needs to change if you run a replicate more than twice in a run

RepCar$RepSsrm<-1*(RepCar$Sd^2)
dfrep<-sum(RepCar$Number)-nrow(RepCar)
Srep<-sqrt((sum(RepCar$RepSsrm))/dfrep)

################################################################################
##                                 PRECISION                                  ##
## Random errors within the laboratory = precision ( Menditto et al 2007)     ##
## or repeatability (Carter and Fry 2013) - summed standard deviation of      ##
## all repeated measurements during relevant analytical sessions - including  ##
##  check and calibration RMs (Ssrum) and the replicates (Srep)	              ##
################################################################################

uRw<-sqrt((Ssrm^2)+(Srep^2)/2)

################################################################################		 
##                      Standard uncertainty (Szpak 2017)	                    ##
################################################################################

y<-list(u_bias, uRw)
Uc<-sqrt(sumsq(y))

################################################################################		 
##                      VALUES TO REPORT (cf. Szpak 2017)                     ##
################################################################################
print(uRw) ## precision (u(Rw))
print(u_bias) ## accuracy, or systematic error (u(bias))
print(Uc) ## total analytical uncertainty (Uc)