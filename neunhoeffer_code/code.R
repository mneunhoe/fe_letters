
# Runs with RStudio Version 0.99.902 and R version 3.3.1 (2016-06-21) -- "Bug in Your Hair"


#######################################################################################
### Check if necessary packages are installed, install additional packages and load ###
#######################################################################################

# Packages needed for this analysis
list.of.packages <- c("rstudioapi","blockTools", "nbpMatching","dplyr", "plyr","sp","rgdal","RColorBrewer","ri","parallel","experiment")


# Check if all packages are installed. Install new packages if necessary.
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)

# Load all packages
sapply(list.of.packages, require, character.only = TRUE)

# Remove all objects in the Environment
rm(list=ls())

# Set Working Directory to directory with the necessary files. Works only with RStudio!
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Check if Working Directory was set right.
getwd()


# Source auxilliary functions
source("functions/functions.R")
source("functions/FRTCI_function_library.R")

# Turn off plotting device (if active)
if(!dev.cur()==1) dev.off()


# Make Paper if TRUE plots will be saved to directory and number of simulations are high (computationally demanding!). 
# Set TRUE for reproducing results in the paper
MAKE_PAPER <- FALSE


######################################
### Load the pre-treatment dataset ###
######################################

pre_data <- read.csv2("data/pre_data.csv",stringsAsFactors = F)

# Select the communities with <= 1500 eligible voters in 2011 and make sample_data and rest_data
sample_data <- pre_data[which(pre_data$eligible_voters_2011<=1500),]
rest_data <-  pre_data[-which(pre_data$eligible_voters_2011<=1500),]


# block and match small communities according to 
sample_block <- block(sample_data, groups = "electoral_district", n.tr = 2, id.vars = c("community_long_ID")
                      , block.vars = c("eligible_voters_2011", "turnout_2011","green_percent_2011","proportion_male","proportion_german","proportion_under18","proportion_18to64","proportion_over65")
                      , algorithm="optGreedy"
                      , distance = "mahalanobis"
                      , level.two = FALSE, verbose = TRUE)

# Assign one community in each block to treatment. Set Seed to get reproducible results
sample_assignment <- assignment(sample_block,seed=100290)

# Make new directory in working directory for the assignment files
dir.create(paste(getwd(),"/assignment_files",sep=""))

# Temporarily set Working Directory to assignment_files folder
setwd(paste(getwd(),"/assignment_files",sep=""))

# Make .csv files of sample_assignment
outCSV(sample_assignment)

# Make one assignment_list
files <-  list.files(pattern="Group")
assignment_list <-   do.call(rbind, lapply(files, function(x) read.csv(x, stringsAsFactors = FALSE)))

# Order assignment_list according to distance
assignment_list <- assignment_list[order(assignment_list$Distance) , ]

# Add electoral_district for all blocks
assigned.distr1 <- rep(NA,length(assignment_list$Treatment.1))
assigned.distr2 <- rep(NA,length(assignment_list$Treatment.2))

suppressWarnings(
  for (i in 1:length(assignment_list$Treatment.1)){
  assigned.distr1[i] <- sample_data$electoral_district[sample_data$community_long_ID==assignment_list$Treatment.1[i]]
  }
) 

suppressWarnings(
  for (i in 1:length(assignment_list$Treatment.2)){
  assigned.distr2[i] <- sample_data$electoral_district[sample_data$community_long_ID==assignment_list$Treatment.2[i]]
  }
)

assigned.district <- ifelse(is.na(assigned.distr1),assigned.distr2,assigned.distr1 )

assignment_list$electoral_district <- assigned.district

# Give useful names to variables
names(assignment_list) <- c("Rank_in_electoral_district","Treatment","Control","Multivariate_Distance","electoral_district")

# Make assignment_list.csv file
write.csv(assignment_list,paste(Sys.Date(),"_final_assignment_list.csv",sep=""))

################################
### End of random assignment ###
################################

# Reset working directory to source
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Check working directory
getwd()

#############################
### Make figures in paper ###
#############################

# Make graphs directory
if ( MAKE_PAPER ) {
dir.create(paste(getwd(),"/graphs",sep=""))
}

################
### Figure 1 ###
################

# Figure 1 in paper. Map with location of the 185 sample communities.

# Read map data from gadm.org "GADM database of Global Administrative Areas"
gadm <- readRDS("data/DEU_adm4.rds")

bw <- gadm[gadm$NAME_1=="Baden-Württemberg",]

rm(gadm)

# Clean up community names to same format
bw@data$NAME_4[bw@data$NAME_4=="Widdern"] <- "Widdern, Stadt"
bw@data$NAME_4[bw@data$NAME_4=="Hettingen"] <- "Hettingen, Stadt"
bw@data$NAME_4[bw@data$NAME_4=="Langenburg"] <- "Langenburg, Stadt"

# Set Colors for the map
myColours <- rep("white", length(bw@data$NAME_4))
myColours[bw@data$NAME_4%in%sample_data$community_name] <- gray.colors(3)[1]
myColours[64] <- "white" #pick right Altdorf
myColours[116] <- "white" #pick right Altheim
myColours[370] <- "white" #pick right Dürnau
myColours[403] <- "white" #pick right Talheim

# Make Figure 1 from Paper. Map with the location of the 185 sample communities
# Uncomment next line to save map to graphs directory
if ( MAKE_PAPER ) {
  pdf("graphs/map_sample.pdf",width=9, height=8)
}

par(
  family = "serif",
  oma = c(1,1,1,1),
  mar = c(5,5,5,5),
  adj = 0.5
)
plot(bw,col = myColours, border = 'gray80')
title("The Location of the 185 small communities \n in Baden-Württemberg")
legend("bottom", inset=c(0,-0.1), legend=c("Community with <= 1,500 eligible voters in 2011","Community with > 1,500  eligible voters in 2011"), fill=c(gray.colors(3)[1],"white"),border="gray80",bty="n",cex=0.8,xpd=TRUE)

if ( MAKE_PAPER ) {
  dev.off()
}

################
### Figure 2 ###
################


##############################################################
### Load final dataset for the 82 experimental communities ###
##############################################################

plot_data <- read.csv2("data/plot_data.csv",stringsAsFactors = F)


# Figure 2 from paper. Map for treatment and control communities
sample <- myColours==gray.colors(3)[1]

trt_map <- which(sample)[bw@data$NAME_4[sample]%in%as.character(plot_data$community_name)[plot_data$treatment_indicator_Z==1]]
con_map <- which(sample)[bw@data$NAME_4[sample]%in%as.character(plot_data$community_name)[plot_data$treatment_indicator_Z==0]]

myColours <- rep("white", length(bw@data$NAME_4))
myColours[trt_map] <- gray.colors(3)[1]
myColours[con_map] <- gray.colors(3)[2]

# Uncomment next line to save map to graphs directory
if ( MAKE_PAPER ) {
  pdf("graphs/map_trt_con.pdf",width=9, height=8)
}
par(
  family = "serif",
  oma = c(1,1,1,1),
  mar = c(5,5,5,5),
  adj = 0.5
)
plot(bw,col = myColours, border = 'gray80')
title("The Location of the Treatment and Control Communities")
legend("bottom", inset=c(0,-0.1), legend=c("Treatment","Control"), fill=c(gray.colors(3)[1],gray.colors(3)[2]),border="gray80",bty="n",cex=0.8,xpd=TRUE)

if ( MAKE_PAPER ) {
  dev.off()
}


################
### Figure 3 ###
################

# Parameters for Power Calculations

assumed.ate <- seq(0.001,0.02,0.002)
n.pairs <- 2:100


# Power Plot

power <- vector("list",length(assumed.ate))
names(power) <- paste(assumed.ate)
for(i in 1:length(assumed.ate)){
  power[[i]] <- power.function(assumed.ate = assumed.ate[i],var.assumed.diff = 0.001,n.pairs = n.pairs)
  names(power[[i]]) <- paste(n.pairs,"pairs")
}

col <- rev(grey.colors(length(power)))

if ( MAKE_PAPER ) {
col <- rep("black",10)
}
if ( MAKE_PAPER ) {
  pdf("graphs/power_curves.pdf",width=9, height=8)
}
par(
  family = "serif",
  oma = c(1,1,1,1),
  mar = c(5,5,5,2),
  adj = 0.5
)
plot(n.pairs,power[[1]],type="l"
     ,ylim=c(0,1),col=col[1]
     ,ylab=expression(paste("Power of the Experiment 1-",beta))
     ,xlab="Number of Matched Pairs m"
     ,main=bquote(atop(bold("Power of the Experiment for different hypothetical ATE"),
                       "and " ~ alpha ~ " = 0.05")))
for(i in 2:length(power)){
  lines(n.pairs,power[[i]],col=col[i])
}
abline(h=0.8,lty=2)
abline(v=41,lty=2)
label <- paste("ATE",names(power))
text(75,power[[1]][74]+0.03, label[1],cex = 0.7,col=col[1]) 
text(70,power[[2]][69]+0.03, label[2],cex = 0.7,col=col[2],srt=8) 
text(65,power[[3]][64]+0.03, label[3],cex = 0.7,col=col[3],srt=13) 
text(60,power[[4]][59]+0.03, label[4],cex = 0.7,col=col[4],srt=18) 
text(55,power[[5]][54]+0.03, label[5],cex = 0.7,col=col[5],srt=23) 
text(50,power[[6]][49]+0.03, label[6],cex = 0.7,col=col[6],srt=28) 
text(45,power[[7]][44]+0.03, label[7],cex = 0.7,col=col[7],srt=30) 
text(40,power[[8]][39]+0.03, label[8],cex = 0.7,col=col[8],srt=32) 
text(35,power[[9]][34]+0.03, label[9],cex = 0.7,col=col[9],srt=35) 
text(30,power[[10]][29]+0.03, label[10],cex = 0.7,col=col[10],srt=38) 

if ( MAKE_PAPER ) {
  dev.off()
}

# Power Calculation with more general functional form
estimation_data <- read.csv2("data/estimation_data.csv")

grp <- as.numeric(as.factor(plot_data$community_ID))
power.est <- ATEcluster(Y=turnout_2011,Z=Z,data=estimation_data,grp=grp,match = pair)



dif.k <- power.est$diff
w.k <- aggregate(estimation_data$eligible_voters_2011,by=list(estimation_data$pair),sum)$x

n <- sum(estimation_data$eligible_voters_2011)

# Define sensible search ATE for power = 0.8
search.ate <- seq(0.014,0.018,0.0005)

power.neq <- power.function.neq(assumed.ate = search.ate,n.pairs = 41,w.k=w.k,dif.k=dif.k,n=n)
names(power.neq) <- paste(search.ate, "ATE")

# Display power for search ATE
round(power.neq,2)==0.8

# Display minimal detectable effect
round(power.neq,2)[round(power.neq,2)==0.8]


################
### Figure 4 ###
################

# Figure 4 in paper. Covariate Balance in treatment and control communities

tmp <- data.frame(plot_data[,c(4,7,9,10,11,13,14:16,18)])
tmp <- data.frame(tmp[,1:2],tmp[,3:9]*100,tmp[,10])


colnames(tmp) <- c("district","Number of Eligible Voters 2011","Turnout in 2011","Vote Share for \n Green Party in 2011","Share of \n Male Population","Share of \n German Citizens","Share of Population \n under 18 years","Share of Population \n between 18 and 65 years","Share of Population \n over 65 years","Z")
tmp_plot_data <- as.matrix(tmp)
tmp_plot_data[,1] <- as.numeric(as.factor(tmp_plot_data[,1]))
tmp_plot_data <- data.frame(rbind(unname(cbind(tmp_plot_data[,c(1,10)],rep(colnames(tmp_plot_data)[1],82)))
                                  ,unname(cbind(tmp_plot_data[,c(2,10)],rep(colnames(tmp_plot_data)[2],82)))
                                  ,unname(cbind(tmp_plot_data[,c(3,10)],rep(colnames(tmp_plot_data)[3],82)))
                                  ,unname(cbind(tmp_plot_data[,c(4,10)],rep(colnames(tmp_plot_data)[4],82)))
                                  ,unname(cbind(tmp_plot_data[,c(5,10)],rep(colnames(tmp_plot_data)[5],82)))
                                  ,unname(cbind(tmp_plot_data[,c(6,10)],rep(colnames(tmp_plot_data)[6],82)))
                                  ,unname(cbind(tmp_plot_data[,c(7,10)],rep(colnames(tmp_plot_data)[7],82)))
                                  ,unname(cbind(tmp_plot_data[,c(8,10)],rep(colnames(tmp_plot_data)[8],82)))
                                  ,unname(cbind(tmp_plot_data[,c(9,10)],rep(colnames(tmp_plot_data)[9],82)))))


names(tmp_plot_data) <- c("value","Z","Variable")

tmp_plot_data$value <- as.numeric(as.character(tmp_plot_data$value))
tmp_plot_data$Z <- as.numeric(as.character(tmp_plot_data$Z))
tmp_plot_data$Variable <- as.character(tmp_plot_data$Variable)

number <- rep(NA,738)

for(i in 1:9){
  number[(i*82-81):(i*82)] <- rep(0+i,82)
}

number <- rev(number)
tmp_plot_data$number <- number

# Uncomment next line to make plot to graphs directory.
if ( MAKE_PAPER ) {
pdf("graphs/balance.pdf",width=9, height=8)
}
par(
  family = "serif",
  oma = c(1,2,1,1),
  mar = c(5,10,5,2),
  adj = 0.5
)
boxplot(value~Z + number, data = tmp_plot_data[165:738,]
        , at = c(1, 2, 5, 6, 9, 10, 13, 14, 17, 18, 21, 22, 25, 26)
        #, xaxt='n'
        , yaxt='n'
        , ylim = c(0,100)
        , xlab = "in %"
        , col = c(gray.colors(3)[3], gray.colors(3)[2])
        , horizontal = T)

abline(h=19.5,lty=2)
axis(side=2, at=c(25.5, 21.5,17.5,13.5,9.5,5.5,1.5), las=2,labels=unique(tmp_plot_data$Variable)[3:9], line=0, lwd=0)
axis(side=2, at=19.5, las=2,labels="- - - - - - - - - - - - - - - - - - - - -", line=0, lwd=0)
axis(side=2, at=23.5,outer=T,labels=" Previous \n Elections ",cex.axis=0.8, line=-1.5,lwd=0)
axis(side=2, at=9.5,outer=T,labels=" Socio-Demographic Variables ",cex.axis=0.8, line=-1,lwd=0)

title('Balance of Covariates \n in Treatment and Control Groups')

# Add lables to the boxplots
text(c(80, 80), c(26,25), c('Treatment', 'Control'),adj=0,cex=0.8)

if ( MAKE_PAPER ) {
  dev.off()
}


####################################################
### Randomization Inference using the ri package ###
####################################################
rm(list=setdiff(ls(), "MAKE_PAPER"))
source("functions/functions.R")
source("functions/FRTCI_function_library.R")
estimation_data <- read.csv2("data/estimation_data.csv",stringsAsFactors = F)
plot_data <- read.csv2("data/plot_data.csv")
# Set seed to make results reproducible
set.seed(100290)

#####################################
# Make long data to use ri package
#####################################


N <- sum(estimation_data$eligible_voters_16)
n <- estimation_data$eligible_voters_16
Z <- estimation_data$Z
green <- estimation_data$green_abs_16
voted <- estimation_data$voters_16
cluster <- 1:82
block <- estimation_data$pair

long_ri <- data.frame(NULL)

long_Z <- rep(Z,n)
long_cluster <- rep(cluster,n)
long_block <- rep(block,n)
long_voted <- NULL
for(i in 1:length(n)){
  long_voted <- append(long_voted,c(rep(0,n[i]-voted[i]),rep(1,voted[i])))
}
long_green <-NULL
for(i in 1:length(n)){
  long_green <- append(long_green,c(rep(0,n[i]-green[i]),rep(1,green[i])))
}

long_ri <- data.frame(clustvar=long_cluster,blockvar=long_block,Z=long_Z,voted=long_voted,green=long_green)

###################################
# Set up variables for ri package #
###################################

###############################
# First, randomization checks #
###############################

blockvar <- estimation_data$pair
clustvar <- 1:82
Z <- estimation_data$Z

probs <- genprobexact(Z,blockvar,clustvar)

numiter <- 100 # Number of possible randomizations is 2 199 023 255 552, for reproduction of results set to 5000 simulations

if ( MAKE_PAPER ) {
  numiter <- 5000
}

perms <- genperms(Z,blockvar,clustvar,maxiter=numiter) 

covs <- as.matrix(plot_data[,c(4,6,9:11,13:16)])
covs <- cbind(as.factor(covs[,1]),as.numeric(covs[,2]),as.numeric(covs[,3]),as.numeric(covs[,4]),as.numeric(covs[,5]),as.numeric(covs[,6]),as.numeric(covs[,7]),as.numeric(covs[,8]),as.numeric(covs[,9]))

##########################
# Calculate F-Statistics #
##########################

Fstat <- summary(lm(Z~covs))$fstatistic[1]   # F-statistic from actual data

Fstatstore <- rep(NA,numiter)
for (i in 1:numiter) {
  Fstatstore[i] <- summary(lm(perms[,i]~covs))$fstatistic[1]   # F-statistic under the null of random assignment of Z
}

################
# Show p-value #
################

mean(Fstatstore >= Fstat)


############################################
# Prepare Data for randomization inference #
############################################

Z <- long_ri$Z
blockvar <- long_ri$block
clustvar <- long_ri$clustvar

probs <- genprobexact(Z,blockvar,clustvar)


ptm <- proc.time()
perms <- genperms(Z,blockvar,clustvar,maxiter=numiter) 
proc.time()-ptm

###################################
# P-values for the ATE on turnout #
###################################


# Generate full schedule of potential outcomes with ATE = 0
Ys <- genouts(Y=long_ri$voted,Z,ate=0)

# Observed difference-in-means ATE
ate <- estate(Y=long_ri$voted,Z,prob = probs)

# Generate Randomization Distribution
ptm <- proc.time()
distout <- gendist(Ys,perms,prob=probs)
proc.time()-ptm

# Density Line for Plot
d <- density(distout)

# Save results for plot
plot <- dispdist(distout,ate,display.plot = F)

# Display Results
plot 

# Difference-in-totals ATE
ateHT <- estate(Y=long_ri$voted,Z,prob=probs,HT=T)

# Difference-in-totals Randomization Distribution
ptm <- proc.time()
distoutHT <- gendist(Ys,perms,prob=probs,HT=T)
proc.time()-ptm

# Density line for plot
dHT <- density(distoutHT)

# Save Results
plotHT <- dispdist(distoutHT,ateHT,display.plot = F)

# Display Results
plotHT

# Set limits for plot
xlim <- c(-max(c(abs(min(distoutHT)),max(distoutHT))),max(c(abs(min(distoutHT)),max(distoutHT))))

################
### Figure 5 ###
################

if ( MAKE_PAPER ) {
pdf("graphs/p_plot_turnout.pdf",width = 10,height = 7)
  adapt <- 40
}

adapt <- 0
col <- grey.colors(3)
par(
  family = "serif",
  oma = c(1,1,1,1),
  mar = c(10,4,10,4),
  adj = 0.5
)
par(mfrow=c(1,2))

hist(distout, freq = F, ylab="",xlab = "Simulated ATE", main = paste("Difference-in-means"), 
     breaks = (length(distout)^0.5)-adapt, lwd = 1,col=col[3],border="white",xlim=xlim,cex.main=0.9,axes=F)
axis(side=1)
boxplot(c(plot$quantile[1],quantile(distout,.25)
          ,0,quantile(distout,.75),plot$quantile[2])
        ,horizontal=T,add=T,boxwex=12,frame=F,at=4.5
        ,col=rgb(red=1,green=1,blue=1,alpha=0.7,maxColorValue = 1))
lines(c(ate, ate), c(-1000, 1000), col = col[1], lwd = 2, 
      lty = 2)
lines(d)

mtext(text=paste("observed \n ATE = ",round(ate,4)),side=3,at=ate,cex=0.7)
text(ate-0.0035 ,38,paste("p =", plot$'two.tailed.p.value'),srt=270,col=col[1],cex=0.7)




hist(distoutHT, freq = F, ylab="",xlab = "Simulated ATE", main = paste("Difference-in-totals"), 
     breaks = length(distoutHT)^0.5, lwd = 1,col=col[3],border="white",xlim=xlim,cex.main=0.9,axes=F)
axis(side=1)
boxplot(c(plotHT$quantile[1],quantile(distoutHT,.25)
          ,0,quantile(distoutHT,.75),plotHT$quantile[2])
        ,horizontal=T,add=T,boxwex=4,frame=F,at=2
        ,col=rgb(red=1,green=1,blue=1,alpha=0.7,maxColorValue = 1))

lines(c(ateHT, ateHT), c(-1000, 1000), col = col[1], lwd = 2, 
      lty = 2)
lines(dHT)

mtext(text=paste("observed \n ATE = ",round(ateHT,4)),side=3,at=ateHT,cex=0.7)
text(ateHT+0.0035,13,paste("p =", plotHT$'two.tailed.p.value'),srt=270,col=col[1],cex=0.7)

title("The ATE on turnout"
      ,outer=TRUE,line=-3,cex.main=1.5)
mtext(text="Randomization Distributions from 5,000 simulations of the ATE under the sharp null of no treatment effect \n for the difference-in-means estimator and the difference-in-totals estimator."
      , side=1,outer=T,line=-3,adj=0)

if ( MAKE_PAPER ) {
  dev.off()
}


############################################
# P-values for the ATE on green vote share #
############################################

# Generate full schedule of potential outcomes with ATE = 0
Ys_green <- genouts(Y=long_ri$green,Z,ate=0)

# Observed difference-in-means ATE on Green vote share
ate_green <- estate(Y=long_ri$green,Z,prob = probs)

# Generate Randomization Distribution
ptm <- proc.time()
distout_green <- gendist(Ys_green,perms,prob=probs)
proc.time()-ptm

# Density line for plot
d_green <- density(distout_green)

# Save Results for Plot
plot_green <- dispdist(distout_green,ate_green,display.plot = F)

# Display results
plot_green

# Difference-in-totals ATE on Green vote share
ateHT_green <- estate(Y=long_ri$green,Z,prob=probs,HT=T)

# Generate Randomization Distribution
ptm <- proc.time()
distoutHT_green <- gendist(Ys_green,perms,prob=probs,HT=T)
proc.time()-ptm

# Density line for plot
dHT_green <- density(distoutHT_green)

# Save Results for plot
plotHT_green <- dispdist(distoutHT_green,ateHT_green,display.plot = F)

# Display Results
plotHT_green

# Limit x-axis in plot
xlim_green <- c(-max(c(abs(min(distoutHT_green)),max(distoutHT_green))),max(c(abs(min(distoutHT_green)),max(distoutHT_green))))


################
### Figure 7 ###
################

if ( MAKE_PAPER ) {
  pdf("graphs/p_plot_green.pdf",width = 10,height = 7)
}
col <- grey.colors(3)
par(
  family = "serif",
  oma = c(1,1,1,1),
  mar = c(10,4,10,4),
  adj = 0.5
)
par(mfrow=c(1,2))

hist(distout_green, freq = F, ylab="",xlab = "Simulated ATE", main = paste("Difference-in-means"), 
     breaks = (length(distout_green)^0.5)-adapt, lwd = 1,col=col[3],border="white",xlim=xlim_green,cex.main=0.9,axes=F)
axis(side=1)
boxplot(c(plot_green$quantile[1],quantile(distout_green,.25)
          ,0,quantile(distout_green,.75),plot_green$quantile[2])
        ,horizontal=T,add=T,boxwex=20,frame=F,at=6.5
        ,col=rgb(red=1,green=1,blue=1,alpha=0.7,maxColorValue = 1))
lines(c(ate_green, ate_green), c(-1000, 1000), col = col[1], lwd = 2, 
      lty = 2)
lines(d_green)

mtext(text=paste("observed \n ATE = ",round(ate_green,4)),side=3,at=ate_green,cex=0.7)
text(ate_green-0.0025,60,paste("p =", plot_green$'two.tailed.p.value'),srt=270,col=col[1],cex=0.7)

hist(distoutHT_green, freq = F, ylab="",xlab = "Simulated ATE", main = paste("Difference-in-totals"), 
     breaks = length(distoutHT_green)^0.5, lwd = 1,col=col[3],border="white",xlim=xlim_green,cex.main=0.9,axes=F)
axis(side=1)
boxplot(c(plotHT_green$quantile[1],quantile(distoutHT_green,.25)
          ,0,quantile(distoutHT_green,.75),plotHT_green$quantile[2])
        ,horizontal=T,add=T,boxwex=15,frame=F,at=5
        ,col=rgb(red=1,green=1,blue=1,alpha=0.7,maxColorValue = 1))

lines(c(ateHT_green, ateHT_green), c(-1000, 1000), col = col[1], lwd = 2, 
      lty = 2)
lines(dHT_green)

mtext(text=paste("observed \n ATE = ",round(ateHT_green,4)),side=3,at=ateHT_green,cex=0.7)
text(ateHT_green+0.0025,40,paste("p =", plotHT_green$'two.tailed.p.value'),srt=270,col=col[1],cex=0.7)

title("The ATE on green vote share"
      ,outer=TRUE,line=-3,cex.main=1.5)
mtext(text="Randomization Distributions from 5,000 simulations of the ATE under the sharp null of no treatment effect \n for the difference-in-means estimator and the difference-in-totals estimator."
      , side=1,outer=T,line=-3,adj=0)

if ( MAKE_PAPER ) {
  dev.off()
}

#######################################################
# CIs assuming observed ATE is the true ATE - Turnout #
#######################################################

Ys_CI <- genouts(Y=long_ri$voted,Z,ate=ate)

ptm <- proc.time()
distout_CI <- gendist(Ys_CI,perms,prob=probs)
proc.time()-ptm

d_CI <- density(distout_CI)

plot_CI <- dispdist(distout_CI,ate,display.plot = F)
plot_CI

YsHT_CI <- genouts(Y=long_ri$voted,Z,ate=ateHT)


ptm <- proc.time()
distoutHT_CI <- gendist(YsHT_CI,perms,prob=probs,HT=T)
proc.time()-ptm

dHT_CI <- density(distoutHT_CI)
plotHT_CI <- dispdist(distoutHT_CI,ateHT,display.plot = F)
plotHT_CI

################
### Figure 6 ###
################

if ( MAKE_PAPER ) {
  pdf("graphs/p_plot_turnout_CI.pdf",width = 10,height = 7)
}

col <- grey.colors(3)
par(
  family = "serif",
  oma = c(1,1,1,1),
  mar = c(10,4,10,4),
  adj = 0.5
)
par(mfrow=c(1,2))

hist(distout_CI, freq = F, ylab="",xlab = "Simulated ATE", main = paste("Difference-in-means"), 
     breaks = (length(distout_CI)^0.5), lwd = 1,col=col[3],border="white",cex.main=0.9,axes=F)
axis(side=1)
boxplot(c(plot_CI$quantile[1],quantile(distout_CI,.25)
          ,ate,quantile(distout_CI,.75),plot_CI$quantile[2])
        ,horizontal=T,add=T,boxwex=15,frame=F,at=5
        ,col=rgb(red=1,green=1,blue=1,alpha=0.7,maxColorValue = 1))
lines(c(ate, ate), c(-1000, 1000), col = col[1], lwd = 2, 
      lty = 2)
lines(d_CI)

abline(v=0)

mtext(text=paste("observed \n ATE = ",round(ate,4)),side=3,at=ate,cex=0.7)



hist(distoutHT_CI, freq = F, ylab="",xlab = "Simulated ATE", main = paste("Difference-in-totals"), 
     breaks = length(distoutHT_CI)^0.5, lwd = 1,col=col[3],border="white",cex.main=0.9,axes=F)
axis(side=1)
boxplot(c(plotHT_CI$quantile[1],quantile(distoutHT_CI,.25)
          ,ateHT,quantile(distoutHT_CI,.75),plotHT_CI$quantile[2])
        ,horizontal=T,add=T,boxwex=5,frame=F,at=1.5
        ,col=rgb(red=1,green=1,blue=1,alpha=0.7,maxColorValue = 1))

lines(c(ateHT, ateHT), c(-1000, 1000), col = col[1], lwd = 2, 
      lty = 2)
lines(dHT_CI)

abline(v=0)

mtext(text=paste("observed \n ATE = ",round(ateHT,4)),side=3,at=ateHT,cex=0.7)

title("Confidence Intervals for the ATE on turnout"
      ,outer=TRUE,line=-3,cex.main=1.5)
mtext(text="Randomization Distributions from 5,000 simulations of the ATE under the assumption that the observed ATE is the true ATE \n for the difference-in-means estimator and the difference-in-totals estimator."
      , side=1,outer=T,line=-3,adj=0)

if ( MAKE_PAPER ) {
  dev.off()
}


######################################
# CI for the ATE on green vote share #
######################################

Ys_green_CI <- genouts(Y=long_ri$green,Z,ate=ate_green)



ptm <- proc.time()
distout_green_CI <- gendist(Ys_green_CI,perms,prob=probs)
proc.time()-ptm

d_green_CI <- density(distout_green_CI)

plot_green_CI <- dispdist(distout_green_CI,ate_green,display.plot = F)
plot_green_CI

YsHT_green_CI <- genouts(Y=long_ri$green,Z,ate=ateHT_green)

ptm <- proc.time()
distoutHT_green_CI <- gendist(YsHT_green_CI,perms,prob=probs,HT=T)
proc.time()-ptm

dHT_green_CI <- density(distoutHT_green_CI)
plotHT_green_CI <- dispdist(distoutHT_green_CI,ateHT_green,display.plot = F)
plotHT_green_CI

################
### Figure 8 ###
################

if ( MAKE_PAPER ) {
  pdf("graphs/p_plot_green_CI.pdf",width = 10,height = 7)
}
col <- grey.colors(3)
par(
  family = "serif",
  oma = c(1,1,1,1),
  mar = c(10,4,10,4),
  adj = 0.5
)
par(mfrow=c(1,2))

hist(distout_green_CI, freq = F, ylab="",xlab = "Simulated ATE", main = paste("Difference-in-means"), 
     breaks = (length(distout_green_CI)^0.5), lwd = 1,col=col[3],border="white",cex.main=0.9,axes=F)
axis(side=1)
boxplot(c(plot_green_CI$quantile[1],quantile(distout_green_CI,.25)
          ,ate_green,quantile(distout_green_CI,.75),plot_green_CI$quantile[2])
        ,horizontal=T,add=T,boxwex=20,frame=F,at=7
        ,col=rgb(red=1,green=1,blue=1,alpha=0.7,maxColorValue = 1))

lines(c(ate_green, ate_green), c(-1000, 1000), col = col[1], lwd = 2, 
      lty = 2)
lines(d_green_CI)

abline(v=0)

mtext(text=paste("observed \n ATE = ",round(ate_green,4)),side=3,at=ate_green,cex=0.7)


hist(distoutHT_green_CI, freq = F, ylab="",xlab = "Simulated ATE", main = paste("Difference-in-totals"), 
     breaks = length(distoutHT_green_CI)^0.5, lwd = 1,col=col[3],border="white",cex.main=0.9,axes=F)
axis(side=1)
boxplot(c(plotHT_green_CI$quantile[1],quantile(distoutHT_green_CI,.25)
          ,ateHT_green,quantile(distoutHT_green_CI,.75),plotHT_green_CI$quantile[2])
        ,horizontal=T,add=T,boxwex=15,frame=F,at=5
        ,col=rgb(red=1,green=1,blue=1,alpha=0.7,maxColorValue = 1))
lines(c(ateHT_green, ateHT_green), c(-1000, 1000), col = col[1], lwd = 2, 
      lty = 2)
lines(dHT_green_CI)

abline(v=0)

mtext(text=paste("observed \n ATE = ",round(ateHT_green,4)),side=3,at=ateHT_green,cex=0.7)

title("Confidence Intervals for the ATE on green vote share"
      ,outer=TRUE,line=-3,cex.main=1.5)
mtext(text="Randomization Distributions from 5,000 simulations of the ATE under the assumption that the observed ATE is the true ATE \n for the difference-in-means estimator and the difference-in-totals estimator."
      , side=1,outer=T,line=-3,adj=0)

if ( MAKE_PAPER ) {
  dev.off()
}

###################################################
### Checking robustness with experiment package ###
###################################################


long_data <- read.csv2("data/long_data.csv")


ate.robust <- ATEcluster(Y=voted,Z=Z,data=long_data,grp=cluster,match = block)


# Display ATE on turnout according to Imai et al. 2009
ate.robust$est 

# Display 95%-CI according to Imai et al. 2009
ate.robust$est+qt(c(.025, .975), df=ate.robust$m-1)*sqrt(ate.robust$var)




ate.robust.green <- ATEcluster(Y=green,Z=Z,data=long_data,grp=cluster,match = block)

# Display ATE on turnout according to Imai et al. 2009
ate.robust.green$est 

# Display 95%-CI according to Imai et al. 2009
ate.robust.green$est+qt(c(.025, .975), df=ate.robust.green$m-1)*sqrt(ate.robust.green$var)



###############################################
# Looking for heterogenous treatment effects ##
###############################################

B <- 100
grid.size <- 51

if ( MAKE_PAPER ) {
B <- 1000   # number of trials for permutation test
grid.size <- 301 # how many points to search in grid.
}

blockvar1 <- estimation_data$pair
clustvar1 <- 1:82
Z1 <- estimation_data$Z
Y1 <- as.numeric(estimation_data$turnout_16)


cat( "Testing for idiosyncratic treatment effect variation on turnout \n" )
ptm <- proc.time()
het.test.turnout <-  FRTCI( Y1, Z1, B=B, grid.size = grid.size,
              get.z.star = get.z.star1 )
time <- proc.time()-ptm

het.test.turnout 

Y1 <- estimation_data$green_abs_16/estimation_data$eligible_voters_16

cat( "Testing for idiosyncratic treatment effect variation on Green vote share \n" )
ptm <- proc.time()
het.test.green <-  FRTCI( Y1, Z1, B=B, grid.size = grid.size,
                            get.z.star = get.z.star1 )
time <- proc.time()-ptm

het.test.green
