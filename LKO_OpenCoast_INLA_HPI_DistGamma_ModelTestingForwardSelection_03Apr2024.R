rm(list = ls())

##############################
# Load the data and packages #
##############################
setwd("C:/Users/MidwoodJ/Desktop/LKO_Open_Coast_INLA/Analysis - Reboot 2023")
#setwd('/Users/highstat/appldfat/HighlandStatistdfs/Courses/Data')
df <- read.csv("C:/Users/MidwoodJ/Desktop/LKO_Open_Coast_INLA/Analysis - Reboot 2023/Data/LKO_OC_WorkingDataset_wSlope_N214_02Nov2023.csv")
head(df)
names(df)
str(df)
df$fG3  <- factor(df$Group.3)
df$fG4  <- factor(df$Group.4)
df$fLgStrahelr <- factor(df$RiverStrahler_Large)
#Load support files and packages
library(lattice)
library(ggplot2)
library(rgdal)
library(sp)
library(fields)
library(gstat)
library(ggmap)
library(reshape)
library(raster)
require(dismo)
library(INLA)
library(maps)
library(maptools)
library(mapdata)
library(rgeos)
library(brinla)
library(mgcv)
source("HighstatLibV11.R")

################
# Housekeeping #
################
# To make the numerical estimation process
# easier, we standardize all continous 
# covariates.
df$fG3 <- df$fG3.JDM
df$fG4 <- df$fG4.JDM

df.2 <- df %>% dplyr::select(Ibi,AdjIbi,Hpi,AdjHpi,
                             M1LAT,M1LONG,
                             RiverDist_Large,ClosestDist,Temp,Mean.Fetch,WetlandDist_ALL,MarinaDist_ALL,RiverDist_ALL,Slope, ## Oct2023 - added slope
                             INLA_dist,fG3,fG4,
                             Flag)

df.2$Slope.std          <- MyStd(df.2$Slope) ## Oct2023 - added slope
df.2$Temp.std           <- MyStd(df.2$Temp) 
df.2$Fetch.std          <- MyStd(df.2$Mean.Fetch)
df.2$MarinaDist.std     <- MyStd(df.2$MarinaDist_ALL)
df.2$RiverDist.std      <- MyStd(df.2$RiverDist_ALL)
df.2$WetlandDist.std    <- MyStd(df.2$WetlandDist_ALL)
df.2$ClosestDist.std    <- MyStd(df.2$ClosestDist) 
#df.2$Lat.std            <- MyStd(df.2$M1LAT)
#df.2$Long.std           <- MyStd(df.2$M1LONG)
#df.2$Ibi.std            <- MyStd(df.2$Ibi)


hist(df.2$Hpi,breaks=40)
df.2$HPI.Beta<-((df.2$Hpi/100)*(nrow(df.2)-1)+0.5)/nrow(df.2) # scaling the IBI
hist(df.2$HPI.Beta,breaks=40)
df.2$HPI.Gamma<-df.2$Hpi+0.025
hist(df.2$HPI.Gamma,breaks=40)

#df.2$log.HPI<-log(df.2$Hpi)
df.2$INLA_dist<-df.2$INLA_dist/1000 ## make smaller

# Re-order the data by year (allows for easier plotting)
I1  <- order(df.2$INLA_dist)
df2 <- df.2[I1,]

#df2 <- df2[complete.cases(df2), ]
names(df2)

##########################
## Next we define knots ##     
##########################
Knots <- seq(from = 0, to = 205, by = 2) # 
Knots 
length(Knots) #That isn't much

# Now we have to allocate each year to a knot.
# There is a function fInterval that can do that.
Intervals <- findInterval(df2$INLA_dist,
                          c(0, Knots, 205))
Intervals

# plot knots
range(df2$INLA_dist)
par(mar = c(5,2,2,2), cex.lab = 1.5)
plot(0,0, type = "n",
     ylim = c(0,0.5),
     xlim = c(0,205),
     axes = FALSE,
     xlab = "Distance (km)",
     ylab = "")
for (i in 1:length(Knots)) {
	segments(Knots[i],0,Knots[i], 0.1, lty = 2)
}    
axis(1)
text(df2$INLA_dist, rep(0, nrow(df2)), "|", cex = 1.5)
## SAVE

#########################
## Define spatial mesh ##
#########################
#1. Make a mesh.
mesh1d <- inla.mesh.1d(loc = Knots)
# 2. Define the weighting factors a_ik (also called the projector matrix).
A4     <- inla.spde.make.A(mesh = mesh1d, 
                           loc = df2$INLA_dist)               
dim(A4)
head(A4, 10)

#spde <- inla.spde2.matern(mesh1d,alpha=2)
# 3. define spde (matern correlation structure)
spde <- inla.spde2.pcmatern( ## model for the precision
  mesh=mesh1d, ## mesh supplied
  alpha=2, ## smoothness parameter
  prior.range = c(4, 0.05), ## P(range < 1) = 0.01 - prob. range <50
  prior.sigma = c(1.5, 0.05)) ## P(sigma > 1) = 0.5

# range # diff(range(mesh$loc[, 1]))/2
# sigma # sd(df$y)/10

#######################
## Make global stack ##
#######################
# 5. Make a stack.
X <- data.frame(
  fG2 = model.matrix(~fG3,data = df2)[,2],
  fG3 = model.matrix(~fG3,data = df2)[,3], 
  fG4 = model.matrix(~fG3,data = df2)[,4],
  fG5 = model.matrix(~fG3,data = df2)[,5],
  fG6 = model.matrix(~fG3,data = df2)[,6],
  fG7 = model.matrix(~fG3,data = df2)[,7],
  MarinaDist.std = model.matrix(~MarinaDist.std,data = df2)[,2],
  RiverDist.std = model.matrix(~RiverDist.std,data = df2)[,2],
  Fetch.std = model.matrix(~Fetch.std,data = df2)[,2],
  Temp.std = model.matrix(~Temp.std,data = df2)[,2],
  Slope.std = model.matrix(~Slope.std,data = df2)[,2],
  ClosestDist.std = model.matrix(~ClosestDist.std,data = df2)[,2],
  WetlandDist.std = model.matrix(~WetlandDist.std,data = df2)[,2]
)
head(X)

StackFit <- inla.stack(
             tag = "Fit",
             data = list(HPI.Gamma = df2$HPI.Gamma),  
	         A = list(1, 1,A4),                  
	         effects = list(  
	              Intercept = rep(1,nrow(df2)),  # defining intercept 
	              X    = X,
	              INLA_dist = 1:length(Knots)))

# 6. Specify the model formula in terms of the 
# Use the PC prior
U <- 1
hyper.prec = list(theta1 = list(prior = "pc.prec", param = c(U, 0.05)))



################
## NULL Model ##
################
f.null1 <- HPI.Gamma ~ -1 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7 +
  f(INLA_dist, model = spde) 
I.null1 <- inla(f.null1, 
                control.compute = list(config = TRUE,dic = TRUE,waic=TRUE),
                control.predictor = list(A = inla.stack.A(StackFit),
                                         compute = TRUE, 
                                         quantiles = c(0.025, 0.975)),
                family = "Gamma",
                data = inla.stack.data(StackFit))
summary(I.null1)

f.null2 <- HPI.Gamma ~ -1 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7 
I.null2 <- inla(f.null2, 
                control.compute = list(config = TRUE,dic = TRUE,waic=TRUE),
                control.predictor = list(A = inla.stack.A(StackFit),
                                         compute = TRUE, 
                                         quantiles = c(0.025, 0.975)),
                family = "Gamma",
                data = inla.stack.data(StackFit))
summary(I.null2)

#########################################################
## INLA Step Function Forward Selection w Spatial Term ##
#########################################################
#https://rdrr.io/cran/INLAutils/man/INLAstep.html
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/INLAutils/INLAutils_0.0.5.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
library(INLAutils)

result.marina <- INLAstep(fam1 = "Gamma", 
                          df2,
                          spde=spde,
                          in_stack = StackFit,
                          invariant = "0 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7 + f(INLA_dist, model = spde)",
                          direction = 'forwards',
                          include = 20:22, # marina = 22, river = 23, wetland = 24, closest = 25
                          y = 'HPI.Gamma',
                          y2 = 'HPI.Gamma',
                          powerl = 1,
                          inter = 2, ## number of levels of interactions
                          thresh = 2)
summary(result.marina$best_model)
#autoplot(result.marina$best_model, which = c(1, 5), CI = TRUE)

result.river <- INLAstep(fam1 = "Gamma", 
                         df2,
                         spde=spde,
                         in_stack = StackFit,
                         invariant = "0 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7 + f(INLA_dist, model = spde)",
                         direction = 'forwards',
                         include = c(20:21,23), ## marina = 22, river = 23, wetland = 24, closest = 25
                         y = 'HPI.Gamma',
                         y2 = 'HPI.Gamma',
                         powerl = 1,
                         inter = 2, ## number of levels of interactions
                         thresh = 2)
summary(result.river$best_model)
#autoplot(result.river$best_model, which = c(1, 5), CI = TRUE)

result.wetland <- INLAstep(fam1 = "Gamma", 
                           df2,
                           spde=spde,
                           in_stack = StackFit,
                           invariant = "0 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7 + f(INLA_dist, model = spde)",
                           direction = 'forwards',
                           include = c(20:21,24), ## marina = 22, river = 23, wetland = 24, closest = 25
                           y = 'HPI.Gamma',
                           y2 = 'HPI.Gamma',
                           powerl = 1,
                           inter = 2, ## number of levels of interactions
                           thresh = 2)
summary(result.wetland$best_model)
#autoplot(result.wetland$best_model, which = c(1, 5), CI = TRUE)

result.closest <- INLAstep(fam1 = "Gamma", 
                           df2,
                           spde=spde,
                           in_stack = StackFit,
                           invariant = "0 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7 + f(INLA_dist, model = spde)",
                           direction = 'forwards',
                           include = c(20:21,25), ## marina = 22, river = 23, wetland = 24, closest = 25
                           y = 'HPI.Gamma',
                           y2 = 'HPI.Gamma',
                           powerl = 1,
                           inter = 2, ## number of levels of interactions
                           thresh = 2)
summary(result.closest$best_model)
#autoplot(result.closest$best_model, which = c(1, 5), CI = TRUE)

##########################################################
## INLA Step Function Forward Selection NO Spatial Term ##
##########################################################
result.marina <- INLAstep(fam1 = "Gamma", 
                          df2,
                          spde=spde,
                          in_stack = StackFit,
                          invariant = "0 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7",
                          direction = 'forwards',
                          include = 20:23, # marina = 23, river = 24, wetland = 25, closest = 26
                          y = 'HPI.Gamma',
                          y2 = 'HPI.Gamma',
                          powerl = 1,
                          inter = 2, ## number of levels of interactions
                          thresh = 2)
summary(result.marina$best_model)
#autoplot(result.marina$best_model, which = c(1, 5), CI = TRUE)

result.river <- INLAstep(fam1 = "Gamma", 
                         df2,
                         spde=spde,
                         in_stack = StackFit,
                         invariant = "0 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7",
                         direction = 'forwards',
                         include = c(20:22,24), ## marina = 23, river = 24, wetland = 25, closest = 26
                         y = 'HPI.Gamma',
                         y2 = 'HPI.Gamma',
                         powerl = 1,
                         inter = 2, ## number of levels of interactions
                         thresh = 2)
summary(result.river$best_model)
#autoplot(result.river$best_model, which = c(1, 5), CI = TRUE)

result.wetland <- INLAstep(fam1 = "Gamma", 
                           df2,
                           spde=spde,
                           in_stack = StackFit,
                           invariant = "0 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7",
                           direction = 'forwards',
                           include = c(20:22,25), ## marina = 23, river = 24, wetland = 25, closest = 26
                           y = 'HPI.Gamma',
                           y2 = 'HPI.Gamma',
                           powerl = 1,
                           inter = 2, ## number of levels of interactions
                           thresh = 2)
summary(result.wetland$best_model)
#autoplot(result.wetland$best_model, which = c(1, 5), CI = TRUE)

result.closest <- INLAstep(fam1 = "Gamma", 
                           df2,
                           spde=spde,
                           in_stack = StackFit,
                           invariant = "0 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7",
                           direction = 'forwards',
                           include = c(20:22,26), ## marina = 23, river = 24, wetland = 25, closest = 26
                           y = 'HPI.Gamma',
                           y2 = 'HPI.Gamma',
                           powerl = 1,
                           inter = 2, ## number of levels of interactions
                           thresh = 2)
summary(result.closest$best_model)
#autoplot(result.closest$best_model, which = c(1, 5), CI = TRUE)

###########################
## Try to run best model ##
###########################
f1a <- HPI.Gamma ~ -1 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7 +
  Slope.std:Fetch.std + f(INLA_dist, model = spde) 
f1a <- HPI.Gamma ~ -1 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7 

I1a <- inla(f1a, 
            control.compute = list(config = TRUE,dic = TRUE,waic=TRUE),
            control.predictor = list(A = inla.stack.A(StackFit),
                                     compute = TRUE, 
                                     quantiles = c(0.025, 0.975)),
            family = "Gamma",
            data = inla.stack.data(StackFit))
summary(I1a)


######################
## Model validation ##
######################
IBest <- I1a ## by changing this, you should be able to produce all the figures required to test model fit.
summary(IBest) 
# OUTPUT FOR THE BETAS
Beta1 <- IBest$summary.fixed[, c("mean",
                                 "sd",
                                 "0.025quant",  
                                 "0.975quant")] 
print(Beta1, digits = 3)

# 8. Inspect the results.
## Pearson Residuals for Gamma
mu <- IBest$summary.fitted.values[1:214,"mean"]
phi <- IBest$summary.hyperpar[,"mean"]
Vary <- mu^2/phi
E1 <- (df2$HPI.Gamma - mu) /sqrt(Vary) # Pearson residuals

Fit.inla <- IBest$summary.fitted.values[1:214, "mean"]
plot(E1) # 
plot(acf(E1)) # all good

# And now we can apply our usual tricks
# for model validation
par(mfrow=c(2,2))
par(mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = Fit.inla[1:214], 
     y = E1, 
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lty = 2)
## SAVE
## structure, 

par(mar = c(5, 5, 2, 2), cex.lab = 1.5, pty = "s")
plot(x = Fit.inla[1:214], 
     y = df2$HPI.Gamma, 
     xlim = c(0, 40),
     ylim = c(0, 60),
     xlab = "Fitted values",
     ylab = "Observed HPI.Gamma")
## meh fit...
## SAVE

#Normality (this is the least important assumption)
hist(E1, main = "Normality", breaks=10)
# meh

## SAVE AS: HPI_BestModel_N214_FigureS3_06Apr2024


#Plot residuals versus covariates
MyVar <- c("Fetch.std", "Slope.std") ## include only if treatments are in the model... can change these
df2$E1 <- E1
MyMultipanel.ggp2(Z = df2, 
                  varx = MyVar, 
                  vary = "E1", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
## SAVE AS: HPI_BestModel_N214_Covariates_FigureS3_06Apr2024


# Plot with treatments 
p <- ggplot()
p <- p + geom_point(data = df2, 
                    aes(y = E1, x = INLA_dist),
                    shape = 16, 
                    size = 2)

p <- p + geom_hline(yintercept = 0)      
p <- p + xlab("Distance (km)") + 
  ylab("Residuals")
p <- p + theme(text = element_text(size=15))
p <- p + facet_wrap( ~ fG3, scales = "fixed")
p ## SAVE - lower than usual scores at the eastern sites

## SAVE AS: HPI_BestModel_N214_FigureS4_06Apr2024

# Make a variogram of the residuals
df2$MyOnes <- rep(1, nrow(df2))
MyData <- data.frame(E1    = E1, 
                     Distance  = df2$INLA_dist, 
                     Dummy = df2$MyOnes)
coordinates(MyData) <- c("Distance", "Dummy")
V1 <- variogram(E1 ~ 1, MyData, cutoff = 205)
plot(V1)
# not a great plot for this...

# spatial plot
xyplot(M1LAT~ M1LONG,
       aspect = "iso",
       col = 1, 
       pch = 16,
       data = df2)

# Option 1:
MyCex <- 3 * abs(E1) / max(E1) + 0.5
Sign <- as.numeric(E1 >=0) + 1
MyPch <- c(1, 3)[Sign]
MyCol <- c("black", "red")[Sign]
xyplot(M1LAT~ M1LONG,
       data = df2,
       cex = MyCex,
       pch = MyPch,
       col = MyCol,
       aspect = "iso",
       xlab = list(label = "Latitude", cex = 1.5),
       ylab = list(label = "Longitude", cex = 1.5)
) # 
# SAVE AS: HPI_BestModel_N214_Spatial_FigureS5_06Apr2024



# And this is the trend itself, with 95% CI
Yearsm   <- IBest$summary.random$INLA_dist
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = Yearsm$"ID",
     y = Yearsm$"mean",
     xlab = "Knot",
     ylab = "Trend",
     type = "l",
     ylim = c(min(Yearsm$"0.025quant"), max(Yearsm$"0.975quant")))

lines(Yearsm$"ID", Yearsm$"0.025quant", lty = 2)
lines(Yearsm$"ID", Yearsm$"0.975quant", lty = 2)
## SAVE
# SAVE AS: FishPA_BestModel_N214_RangeTrendPlot_05Apr2024
#####################################
## Predict model output (plotting) ##
#####################################
X <- data.frame(
  fG2 = model.matrix(~fG3,data = df2)[,2],
  fG3 = model.matrix(~fG3,data = df2)[,3], 
  fG4 = model.matrix(~fG3,data = df2)[,4],
  fG5 = model.matrix(~fG3,data = df2)[,5],
  fG6 = model.matrix(~fG3,data = df2)[,6],
  fG7 = model.matrix(~fG3,data = df2)[,7],
  MarinaDist.std = model.matrix(~MarinaDist.std,data = df2)[,2],
  RiverDist.std = model.matrix(~RiverDist.std,data = df2)[,2],
  Fetch.std = model.matrix(~Fetch.std,data = df2)[,2],
  Temp.std = model.matrix(~Temp.std,data = df2)[,2],
  Slope.std = model.matrix(~Slope.std,data = df2)[,2],
  ClosestDist.std = model.matrix(~ClosestDist.std,data = df2)[,2],
  WetlandDist.std = model.matrix(~WetlandDist.std,data = df2)[,2]
)
head(X)

# Stack for fitting model
StackFit.model <- inla.stack(
  tag = "Fit",
  data = list(HPI.Gamma = df2$HPI.Gamma),  
  A = list(1, 1,A4),                  
  effects = list(  
    Intercept = rep(1,nrow(df2)),  # defining intercept 
    X    = X, # covariates
    INLA_dist = 1:length(Knots)) # spatial component
  )

# Build Stack 
# stack for predicting year trend
#range(df2$Temp.std)
range(df2$Fetch.std) #-2.562259  1.757419
range(df2$Slope.std) # -6.522487  1.111866

Slope.x <- seq(from = -6.6, to = 1.2, length = 25) # Fetch values to make prediction for
Fetch.x <- rep(0,length(Slope.x))  # MarinaDist values to make prediction for

#RiverDist.x <- seq(from = -.76, to = 2.65, length = 25) # MarinaDist values to make prediction for
#lope.x <- rep(0,length(RiverDist.x))  # MarinaDist values to make prediction for

Intercept.pred <-rep(1, length(Slope.x))

StackFit.pred <- inla.stack(
  tag = "Pred.Marina",                   # - Name (nametag) of the stack
  data = list(HPI.Gamma = NA),                 # - dependent variable - set to NA to make prediction
  A = list(1, 1, 1),                      # 
  effects = list(                      # - model componetns
    Intercept = Intercept.pred, # - intercept - vector of 1s
    Fetch.std = Fetch.x,
    Slope.std = Slope.x
  )
)

# Combine stacks
StackCombo<- inla.stack(StackFit.model, StackFit.pred) # stack for spatial model



##########################
## plotting interaction ##
##########################
X <- data.frame(
  fG2 = model.matrix(~fG3,data = df2)[,2],
  fG3 = model.matrix(~fG3,data = df2)[,3], 
  fG4 = model.matrix(~fG3,data = df2)[,4],
  fG5 = model.matrix(~fG3,data = df2)[,5],
  fG6 = model.matrix(~fG3,data = df2)[,6],
  fG7 = model.matrix(~fG3,data = df2)[,7],
  MarinaDist.std = model.matrix(~MarinaDist.std,data = df2)[,2],
  RiverDist.std = model.matrix(~RiverDist.std,data = df2)[,2],
  Fetch.std = model.matrix(~Fetch.std,data = df2)[,2],
  Temp.std = model.matrix(~Temp.std,data = df2)[,2],
  Slope.std = model.matrix(~Slope.std,data = df2)[,2],
  ClosestDist.std = model.matrix(~ClosestDist.std,data = df2)[,2],
  WetlandDist.std = model.matrix(~WetlandDist.std,data = df2)[,2]
)
head(X)

# Stack for fitting model
StackFit.model <- inla.stack(
  tag = "Fit",
  data = list(HPI.Gamma = df2$HPI.Gamma),  
  A = list(1, 1,A4),                  
  effects = list(  
    Intercept = rep(1,nrow(df2)),  # defining intercept 
    X    = X, # covariates
    INLA_dist = 1:length(Knots)) # spatial component
)

# Build Stack 
# stack for predicting year trend
range(df2$Slope.std)
range(df2$Fetch.std)

Slope.base <- seq(from = -6.5, to = 1.1, length = 25) # Fetch values to make prediction for
Fetch.base <- seq(from = -2.5, to = 1.7, length = 25) # MarinaDist values to make prediction for

Slope.x <- rep(Slope.base,25) # rep to combine
Fetch.base2 <- rep(Fetch.base,25)# Rep to combine
Fetch.order <- order(Fetch.base2)
Fetch.x <- Fetch.base2[Fetch.order]
X.inter<- data.frame(
  Slope.std = Slope.x,
  Fetch.std = Fetch.x
)

Intercept.pred <-rep(1, length(Fetch.x))


StackFit.pred <- inla.stack(
  tag = "Pred.Marina",                   # - Name (nametag) of the stack
  data = list(HPI.Gamma = NA),                 # - dependent variable - set to NA to make prediction
  A = list(1, 1),                      # 
  effects = list(                      # - model componetns
    Intercept = Intercept.pred, # - intercept - vector of 1s
    X.inter)
)

# Combine stacks
StackCombo<- inla.stack(StackFit.model, StackFit.pred) # stack for spatial model

I.pred <- inla(f1a, 
                 control.compute = list(config = TRUE,dic = TRUE,waic = TRUE),
                 control.predictor = list(A = inla.stack.A(StackCombo),
                                          compute = TRUE, 
                                          link = 1),
                 quantiles = c(0.025, 0.975), 
                 family = "Gamma",
                 data = inla.stack.data(StackCombo))

Index <- inla.stack.index(StackCombo, tag = "Pred.Marina")$data
modI4 = I.pred
pred.glm <- as.data.frame(list(
  x = Fetch.x * sd(df2$Mean.Fetch) + mean(df2$Mean.Fetch),
  y = Slope.x * sd(df2$Slope) + mean(df2$Slope), 
  mean = modI4$summary.fitted.values[Index, "mean"],
  LCI = modI4$summary.fitted.values[Index, "0.025quant"],
  UCI = modI4$summary.fitted.values[Index, "0.975quant"],
  type = "GLM"
))
#pred.glm$IBI<-(((pred.glm$mean*nrow(df2))-0.5)/(nrow(df2)-1))*100 # convert back to IBI score
#pred.glm$LCI.fix<-(((pred.glm$LCI*nrow(df2))-0.5)/(nrow(df2)-1))*100 # convert back to IBI score
#pred.glm$UCI.fix<-(((pred.glm$UCI*nrow(df2))-0.5)/(nrow(df2)-1))*100 # convert back to IBI score

library(akima)
library(rgl)
library(plotly)
> x=runif(1000)
> y=runif(1000)
> z=rnorm(1000)
s<-interp(pred.glm$x,pred.glm$y,pred.glm$mean,extrap=T)
dim(s$z)
surface3d(s$x,s$y,s$z)

pred.glm$Fetch<-pred.glm$x
pred.glm$Slope<-pred.glm$y


p <- plot_ly(pred.glm, x = ~x, y = ~y, z = ~mean, type = 'scatter3d', mode = 'lines',
             opacity = 1, line = list(width = 6, reverscale = FALSE))
p

xlabel <- list(title = "Fetch (km)")
ylabel <- list(title = "Slope")
p <- plot_ly(x = s$x, y = s$y, z = s$z) %>% 
  add_surface() %>% 
  layout(scene=list(xaxis =list(titlefont=list(size=30),
                                tickfont = list(size=17),
                                title = "Fetch (km)"), 
                    yaxis = list(titlefont=list(size=30),
                                 tickfont = list(size=17),
                                 title = "Slope"), 
                    zaxis = list(titlefont=list(size=30),
                                 tickfont = list(size=17),
                                 title="Habitat Productivity Index")))
p

## SAVE AS: HPI_BestModel_N214_CovariateInteractions_Figure_31Jan2024
png("HPI_BestModel_N214_CovariateInteractions_Figure_03Apr2024.png",
    width = 2400, height = 2400,units="px",res=300)
p
dev.off()


p <- ggplot()
p <- p + geom_point(data = df2, 
                    aes(y = Hpi, x = Mean.Fetch),
                    shape = 1, 
                    size = 1)
p <- p + xlab("Mean.Fetch (m)") + ylab("HPI")
p <- p + theme(text = element_text(size=25))
p <- p + geom_line(data = pred.glm, 
                   aes(x = x, y = mean), 
                   colour = "black")
p <- p + geom_ribbon(data = pred.glm, 
                     aes(x = x, 
                         ymax = UCI, 
                         ymin = LCI),
                     alpha = 0.4)
p



#########################
## Linear combinations ##
#########################
## see https://becarioprecario.bitbucket.io/inla-gitbook/ch-INLAfeatures.html

df2$Intercept <- 1
lcI2 <- inla.make.lincomb(Intercept = 1, fG2 = -1) # intercept vs fG2
lcI3 <- inla.make.lincomb(Intercept = 1, fG3 = -1) # intercept vs fG3
lcI4 <- inla.make.lincomb(Intercept = 1, fG4 = -1) # intercept vs fG4
lcI5 <- inla.make.lincomb(Intercept = 1, fG5= -1) # intercept vs fG5
lcI6 <- inla.make.lincomb(Intercept = 1, fG6= -1) # intercept vs fG5
lcI7 <- inla.make.lincomb(Intercept = 1, fG7= -1) # intercept vs fG5
lc23 <- inla.make.lincomb(fG2 = 1, fG3 = -1) # fg2 vs fg3
lc24 <- inla.make.lincomb(fG2 = 1, fG4 = -1) # fg2 vs fg4
lc25 <- inla.make.lincomb(fG2 = 1, fG5= -1) # fg2 vs fg5
lc26 <- inla.make.lincomb(fG2 = 1, fG6= -1) # fg2 vs fg5
lc27 <- inla.make.lincomb(fG2 = 1, fG7= -1) # fg2 vs fg5
lc34 <- inla.make.lincomb(fG3 = 1, fG4 = -1) # fg3 vs fg4
lc35 <- inla.make.lincomb(fG3 = 1, fG5= -1) # fg3 vs fg5
lc36 <- inla.make.lincomb(fG3 = 1, fG6 = -1) # fg3 vs fg4
lc37 <- inla.make.lincomb(fG3 = 1, fG7= -1) # fg3 vs fg5
lc45 <- inla.make.lincomb(fG4 = 1, fG5= -1) # fg4 vs fg5
lc46 <- inla.make.lincomb(fG4 = 1, fG6= -1) # fg4 vs fg5
lc47 <- inla.make.lincomb(fG4 = 1, fG7= -1) # fg4 vs fg5
lc56 <- inla.make.lincomb(fG5 = 1, fG6= -1) # fg4 vs fg5
lc57 <- inla.make.lincomb(fG5 = 1, fG7= -1) # fg4 vs fg5
lc67 <- inla.make.lincomb(fG6 = 1, fG7= -1) # fg4 vs fg5
names(lcI2) = "lcI2"
names(lcI3) = "lcI3"
names(lcI4) = "lcI4"
names(lcI5) = "lcI5"
names(lcI6) = "lcI6"
names(lcI7) = "lcI7"
names(lc23) = "lc23"
names(lc24) = "lc24"
names(lc25) = "lc25"
names(lc26) = "lc26"
names(lc27) = "lc27"
names(lc34) = "lc34"
names(lc35) = "lc35"
names(lc36) = "lc36"
names(lc37) = "lc37"
names(lc45) = "lc45"
names(lc46) = "lc46"
names(lc47) = "lc47"
names(lc56) = "lc56"
names(lc57) = "lc57"
names(lc67) = "lc67"


lc <- c(lcI2, lcI3, lcI4, lcI5,lcI6,lcI7, lc23, lc24, lc25, lc26, lc27, lc34, lc35,lc36,lc37, lc45, lc46, lc47,lc56,lc57,lc67)

I.lc <- inla(f1a, 
            control.compute = list(config = TRUE,dic = TRUE,waic = TRUE),
            control.predictor = list(A = inla.stack.A(StackFit),
                                     compute = TRUE, 
                                     link = 1),
            quantiles = c(0.025, 0.975), 
            family = "Gamma",
            lincomb = lc,
            #control.inla = list(lincomb.derived.only = F),
            data = inla.stack.data(StackFit))


summary(I.lc)
HyperPar.I.lc.pc <- bri.hyperpar.summary(I.lc)
round(HyperPar.I.lc.pc, digits = 3)

LC <- I.lc$summary.lincomb.derived[,c("mean", "0.025quant", "0.975quant")]
round(LC, digits = 3)

plot(LC)

#########################
## Boxplotish Figure ? ##
#########################
lcpredI = inla.make.lincomb(Intercept = 1, fG2 = 0)
lcpred1 = inla.make.lincomb(Intercept = 1, fG2 = 1)
lcpred2 = inla.make.lincomb(Intercept = 1, fG3 = 1)
lcpred3 = inla.make.lincomb(Intercept = 1, fG4 = 1)
lcpred4 = inla.make.lincomb(Intercept = 1, fG5 = 1)
lcpred5 = inla.make.lincomb(Intercept = 1, fG6 = 1)
lcpred6 = inla.make.lincomb(Intercept = 1, fG7 = 1)
names(lcpredI) = "lcpredI"
names(lcpred1) = "lcpred1"
names(lcpred2) = "lcpred2"
names(lcpred3) = "lcpred3"
names(lcpred4) = "lcpred4"
names(lcpred5) = "lcpred5"
names(lcpred6) = "lcpred6"

lcpred<-c(lcpredI,lcpred1,lcpred2,lcpred3,lcpred4,lcpred5,lcpred6)

I.lc.plot <- inla(f1a, 
            control.compute = list(config = TRUE,dic = TRUE,waic = TRUE),
            control.predictor = list(A = inla.stack.A(StackFit),
                                     compute = TRUE, 
                                     link = 1),
            quantiles = c(0.025, 0.975), 
            family = "Gamma",
            lincomb = lcpred,
            #control.inla = list(lincomb.derived.only = F),
            data = inla.stack.data(StackFit))

summary(I.lc.plot)
HyperPar.I.lc.plot.pc <- bri.hyperpar.summary(I.lc.plot)
round(HyperPar.I.lc.plot.pc, digits = 3)

LC <- I.lc.plot$summary.lincomb.derived[,c("mean", "0.025quant", "0.975quant")]
round(LC, digits = 3)
LC$Treatment<-c("Gravel Beach","Groyne","Mixed Shoreline","Mixed Beach","Revetment","Surcharged Revetment","Sand Beach")
LC$Upper<-LC$'0.975quant'
LC$Lower<-LC$'0.025quant'

p <- ggplot()
p <- p + geom_point(data = LC,
                    aes(x = Treatment, 
                        y = mean)
)
p <- p + geom_errorbar(data = LC,
                       aes(x = Treatment, 
                           ymax = Upper, 
                           ymin = Lower), 
                       width=0.2)
p <- p + xlab("Treatment") + ylab("Total Catch")
p <- p + theme(text = element_text(size = 15)) 
p <- p + theme(legend.position="none") 
p
## SAVE AS: HPI_BestModel_N214_FigureS2_06Apr2024


#I9c$marginals.lincomb$lcpred1
plot(I.lc.plot$marginals.lincomb.derived$lcpredI,type="l",xlim=c(-2,4),
     ylim=c(0,2.1),lwd=2,lty=2,col="grey",
     xlab="Posterior Distribution",ylab="Probability",main="HPI")
lines(I.lc.plot$marginals.lincomb.derived$lcpred1,lwd=2,col="black")
lines(I.lc.plot$marginals.lincomb.derived$lcpred2,lwd=2,col="blue")
lines(I.lc.plot$marginals.lincomb.derived$lcpred3,lwd=2,lty=2,col="blue")
lines(I.lc.plot$marginals.lincomb.derived$lcpred4,lwd=2,lty=1,col="brown")
lines(I.lc.plot$marginals.lincomb.derived$lcpred5,lwd=2,lty=1,col="red")
lines(I.lc.plot$marginals.lincomb.derived$lcpred6,lwd=2,lty=2,col="yellow")


png(file="HPI_BestModel_N214_Lincomb_Figure_05Apr2024.png", res = 300, width = 9, height = 8, units = 'in')
par(mar=c(5,6,4,1))
plot(I.lc.plot$marginals.lincomb.derived$lcpredI,type="l",xlim=c(-2,4),
     ylim=c(0,2.7),lwd=3,lty=2,col="grey23",
     xlab="Posterior Distribution",ylab="Probability",main="Habitat Productivity Index",cex.lab=2, cex.axis=1.5, cex.main=2, cex.sub=1.5)
lines(I.lc.plot$marginals.lincomb.derived$lcpred1,lwd=3,col="black")
lines(I.lc.plot$marginals.lincomb.derived$lcpred2,lwd=3,col="cadetblue4")
lines(I.lc.plot$marginals.lincomb.derived$lcpred3,lwd=3,lty=2,col="blue")
lines(I.lc.plot$marginals.lincomb.derived$lcpred4,lwd=3,lty=1,col="brown")
lines(I.lc.plot$marginals.lincomb.derived$lcpred5,lwd=3,lty=1,col="red")
lines(I.lc.plot$marginals.lincomb.derived$lcpred6,lwd=3,lty=2,col="darkgoldenrod2")
legend(-2,2.5, # places a legend at the appropriate place
       title="Shoreline Type",
       c("Gravel Beach - c", "Groyne - abc","Mixed Shoreline - a","Mixed Beach - bc","Revetment - ab","Surcharged Revetment - bc","Sand Beach - abc"), # puts text in the legend
       lty=c(2,1,1,2,1,1,2), # gives the legend appropriate symbols (lines)
       lwd=c(3.5),
       col=c("grey23","black","cadetblue4","blue","brown","red","darkgoldenrod2"),
       cex=1.25) # gives the legend lines the correct color and width
dev.off()