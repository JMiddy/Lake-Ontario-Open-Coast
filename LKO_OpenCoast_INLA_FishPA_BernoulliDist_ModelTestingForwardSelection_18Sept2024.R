rm(list = ls())

##############################
# Load the data and packages #
##############################
setwd("C:/Users/MidwoodJ/Documents/Projects/LKO/LKO Open Coast/Analysis - Reboot 2023")
#setwd('/Users/highstat/appldfat/HighlandStatistdfs/Courses/Data')
df <- read.csv("C:/Users/MidwoodJ/Documents/Projects/LKO/LKO Open Coast/Analysis - Reboot 2023/Data/LKO_OC_WorkingDataset_wSlope_N214_02Nov2023.csv")
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

df.2 <- df %>% dplyr::select(Ibi,AdjIbi,Hpi,AdjHpi,Numbers,
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

hist(df.2$AdjIbi,breaks=40) # definitely zero-infalted

df.2$AdjIbi01<-ifelse(df.2$AdjIbi==0,0,1)
df.2$AdjIbiModel<-ifelse(df.2$AdjIbi>0,df.2$AdjIbi,NA)
hist(df.2$AdjIbiModel,breaks=40) # definitely zero-infalted

df.2$AdjHpi01<-ifelse(df.2$AdjHpi==0,0,1)
df.2$AdjHpiModel<-ifelse(df.2$AdjHpi>0,df.2$AdjHpi,NA)
hist(df.2$AdjHpiModel,breaks=40) # definitely zero-infalted
#df.2$log.HPI<-log(df.2$Hpi)
df.2$INLA_dist<-df.2$INLA_dist/1000 ## make smaller

# Re-order the data by year (allows for easier plotting)
I1  <- order(df.2$INLA_dist)
df2 <- df.2[I1,]

#df2 <- df2[complete.cases(df2), ]
names(df2)

#df2<-subset(df2,Flag!="Add") ## drop four sites with >60 catch

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

#######################
## Make global stack ##
#######################
# 5. Make a stack.
Xm <- model.matrix(~fG3 + WetlandDist.std + RiverDist.std + Fetch.std + Temp.std + ClosestDist.std + MarinaDist.std + Slope.std, data=df2) ## add slope Oct2023
N <- nrow(df2)

X01 <- data.frame(
  Intercept = rep(1,N),
  fG2 = Xm[,2],
  fG3 = Xm[,3], 
  fG4 = Xm[,4],
  fG5 = Xm[,5],
  fG6 = Xm[,6],
  fG7 = Xm[,7],
  WetlandDist.std =  Xm[,8],
  RiverDist.std = Xm[,9],
  Fetch.std = Xm[,10],
  Temp.std = Xm[,11],
  ClosestDist.std = Xm[,12],
  MarinaDist.std = Xm[,13],
  Slope.std = Xm[,14]
)

StackFit <- inla.stack(
  tag = "Fit",
  data = list(AdjIbi01 = df2$AdjIbi01),  
  A = list(1, 1,A4),                  
  effects = list(  
    Intercept = rep(1,nrow(df2)),  # defining intercept 
    X    = X01,
    INLA_dist = 1:length(Knots)))

# 6. Specify the model formula in terms of the 
# Use the PC prior
U <- 1
hyper.prec = list(theta1 = list(prior = "pc.prec", param = c(U, 0.05)))



################
## NULL Model ##
################
f.null1 <- AdjIbi01 ~ -1 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7 +
  f(INLA_dist, model = spde) 
I.null1 <- inla(f.null1, 
                control.compute = list(config = TRUE,dic = TRUE,waic=TRUE),
                control.predictor = list(A = inla.stack.A(StackFit),
                                         compute = TRUE, 
                                         quantiles = c(0.025, 0.975)),
                family = "binomial",
                data = inla.stack.data(StackFit))
summary(I.null1)

f.null2 <- AdjIbi01 ~ -1 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7 
I.null2 <- inla(f.null2, 
                control.compute = list(config = TRUE,dic = TRUE,waic=TRUE),
                control.predictor = list(A = inla.stack.A(StackFit),
                                         compute = TRUE, 
                                         quantiles = c(0.025, 0.975)),
                family = "binomial",
                data = inla.stack.data(StackFit))
summary(I.null2)

#########################################################
## INLA Step Function Forward Selection w Spatial Term ##
#########################################################
#https://rdrr.io/cran/INLAutils/man/INLAstep.html
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/INLAutils/INLAutils_0.0.5.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
library(INLAutils)

result.marina <- INLAstep(fam1 = "binomial", 
                          df2,
                          spde=spde,
                          in_stack = StackFit,
                          invariant = "0 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7 + f(INLA_dist, model = spde)",
                          direction = 'forwards',
                          include = 20:23, # marina = 23, river = 24, wetland = 25, closest = 26
                          y = 'AdjIbi01',
                          y2 = 'AdjIbi01',
                          powerl = 1,
                          inter = 2, ## number of levels of interactions
                          thresh = 2)
summary(result.marina$best_model)
#autoplot(result.marina$best_model, which = c(1, 5), CI = TRUE)

result.river <- INLAstep(fam1 = "binomial", 
                         df2,
                         spde=spde,
                         in_stack = StackFit,
                         invariant = "0 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7 + f(INLA_dist, model = spde)",
                         direction = 'forwards',
                         include = c(20:22,24), ## marina = 23, river = 24, wetland = 25, closest = 26
                         y = 'AdjIbi01',
                         y2 = 'AdjIbi01',
                         powerl = 1,
                         inter = 2, ## number of levels of interactions
                         thresh = 2)
summary(result.river$best_model)
#autoplot(result.river$best_model, which = c(1, 5), CI = TRUE)

result.wetland <- INLAstep(fam1 = "binomial", 
                           df2,
                           spde=spde,
                           in_stack = StackFit,
                           invariant = "0 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7 + f(INLA_dist, model = spde)",
                           direction = 'forwards',
                           include = c(20:22,25), ## marina = 23, river = 24, wetland = 25, closest = 26
                           y = 'AdjIbi01',
                           y2 = 'AdjIbi01',
                           powerl = 1,
                           inter = 2, ## number of levels of interactions
                           thresh = 2)
summary(result.wetland$best_model)
#autoplot(result.wetland$best_model, which = c(1, 5), CI = TRUE)

result.closest <- INLAstep(fam1 = "binomial", 
                           df2,
                           spde=spde,
                           in_stack = StackFit,
                           invariant = "0 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7 + f(INLA_dist, model = spde)",
                           direction = 'forwards',
                           include = c(20:22,26), ## marina = 23, river = 24, wetland = 25, closest = 26
                           y = 'AdjIbi01',
                           y2 = 'AdjIbi01',
                           powerl = 1,
                           inter = 2, ## number of levels of interactions
                           thresh = 2)
summary(result.closest$best_model)
#autoplot(result.closest$best_model, which = c(1, 5), CI = TRUE)

##########################################################
## INLA Step Function Forward Selection NO Spatial Term ##
##########################################################
result.marina <- INLAstep(fam1 = "binomial", 
                          df2,
                          spde=spde,
                          in_stack = StackFit,
                          invariant = "0 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7",
                          direction = 'forwards',
                          include = 20:23, # marina = 23, river = 24, wetland = 25, closest = 26
                          y = 'AdjIbi01',
                          y2 = 'AdjIbi01',
                          powerl = 1,
                          inter = 2, ## number of levels of interactions
                          thresh = 2)
summary(result.marina$best_model)
#autoplot(result.marina$best_model, which = c(1, 5), CI = TRUE)

result.river <- INLAstep(fam1 = "binomial", 
                         df2,
                         spde=spde,
                         in_stack = StackFit,
                         invariant = "0 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7",
                         direction = 'forwards',
                         include = c(20:22,24), ## marina = 23, river = 24, wetland = 25, closest = 26
                         y = 'AdjIbi01',
                         y2 = 'AdjIbi01',
                         powerl = 1,
                         inter = 2, ## number of levels of interactions
                         thresh = 2)
summary(result.river$best_model)
#autoplot(result.river$best_model, which = c(1, 5), CI = TRUE)

result.wetland <- INLAstep(fam1 = "binomial", 
                           df2,
                           spde=spde,
                           in_stack = StackFit,
                           invariant = "0 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7",
                           direction = 'forwards',
                           include = c(20:22,25), ## marina = 23, river = 24, wetland = 25, closest = 26
                           y = 'AdjIbi01',
                           y2 = 'AdjIbi01',
                           powerl = 1,
                           inter = 2, ## number of levels of interactions
                           thresh = 2)
summary(result.wetland$best_model)
#autoplot(result.wetland$best_model, which = c(1, 5), CI = TRUE)

result.closest <- INLAstep(fam1 = "binomial", 
                           df2,
                           spde=spde,
                           in_stack = StackFit,
                           invariant = "0 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7",
                           direction = 'forwards',
                           include = c(20:22,26), ## marina = 23, river = 24, wetland = 25, closest = 26
                           y = 'AdjIbi01',
                           y2 = 'AdjIbi01',
                           powerl = 1,
                           inter = 2, ## number of levels of interactions
                           thresh = 2)
summary(result.closest$best_model)
#autoplot(result.closest$best_model, which = c(1, 5), CI = TRUE)

###########################
## Try to run best model ##
###########################
f1a <- AdjIbi01 ~ -1 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7 +
  MarinaDist.std + f(INLA_dist, model = spde) 
#f1a <- AdjIbi01 ~ -1 + Intercept + fG2 + fG3 + fG4 + fG5 + fG6 + fG7 +
#  MarinaDist.std * Fetch.std + f(INLA_dist, model = spde) 

I1a <- inla(f1a, 
            control.compute = list(config = TRUE,dic = TRUE,waic=TRUE),
            control.predictor = list(A = inla.stack.A(StackFit),
                                     compute = TRUE, 
                                     quantiles = c(0.025, 0.975)),
            family = "binomial",
            data = inla.stack.data(StackFit))
summary(I1a)


###################################
## Explore Output for best model ##
###################################
IBest <- I1a ## by changing this, you should be able to produce all the figures required to test model fit.
summary(IBest) # 

# OUTPUT FOR THE BETAS
Beta1 <- IBest$summary.fixed[, c("mean",
                                 "sd",
                                 "0.025quant",  
                                 "0.975quant")] 
print(Beta1, digits = 3)
# SAVE

# 8. Inspect the results.

#  Pearson residuals for Bernoulli model
Fit.inla <- IBest$summary.fitted.values[1:214, "mean"]
Pi <- IBest$summary.fitted.values[1:214, "mean"]
ExpY <- Pi
VarY <- (1 - Pi) * Pi
E1 <- (df2$AdjIbi01 - ExpY) / sqrt(VarY)


SpatField.w <- inla.spde2.result(inla = IBest,
                                 name = "INLA_dist",
                                 spde = spde,
                                 do.transfer = TRUE)

Kappa <- inla.emarginal(function(x) x, 
                        SpatField.w$marginals.kappa[[1]] )

Sigma_u <- inla.emarginal(function(x) sqrt(x), 
                          SpatField.w$marginals.variance.nominal[[1]] )

Range <- inla.emarginal(function(x) x, 
                        SpatField.w$marginals.range.nominal[[1]] )

Kappa # 0.125
Sigma_u #1.34
Range #39.6      #Distance at which the correlation diminishes


######################
## Model validation ##
######################
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
## 

par(mar = c(5, 5, 2, 2), cex.lab = 1.5, pty = "s")
plot(x = Fit.inla[1:214], 
     y = df2$AdjIbi01, 
     xlim = c(0, 1),
     ylim = c(0, 1),
     xlab = "Fitted values",
     ylab = "Observed Numbers")
## 
## SAVE

#Normality (this is the least important assumption)
hist(E1, main = "Normality", breaks=10)
# fine

## SAVE AS: FishPA_BestModel_N214_FigureS3_05Apr2024

#Plot residuals versus covariates
MyVar <- c( "MarinaDist.std")
df2$E1 <- E1
MyMultipanel.ggp2(Z = df2, 
                  varx = MyVar, 
                  vary = "E1", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
## SAVE AS: FishPA_BestModel_N214_Covariates_Marina_FigureS3_05Apr2024

MyVar <- c( "Slope.std")
df2$E1 <- E1
MyMultipanel.ggp2(Z = df2, 
                  varx = MyVar, 
                  vary = "E1", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
## SAVE AS: FishPA_BestModel_N214_Covariates_Slope_FigureS3_05Apr2024


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
p ## SAVE AS: FishPA_BestModel_N214_FigureS4_05Apr2024


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

### add lko map
library(pacman)
p_load(rgeos)
p_load(plotKML)
p_load(sp)
p_load(fetchR)
p_load(sf)
p_load(sfheaders)
p_load(rgdal)
filename <- "C:/Users/MidwoodJ/Documents/Projects/LKO/Electrofishing/Regional IBI Comparisons/GIS/GLAF_LKO_Water.shp"
LKO_poly <- st_read(filename)
polygon_list <- st_geometry(LKO_poly)
#LKO_poly.sp = SpatialPolygons(polygon_list,
#                              proj4string = CRS("+proj=utm +zone=17 +datum=WGS84"))
plot(LKO_poly[3], col = "lightgrey")
points(M1LAT~ M1LONG,data=df2,col="red")

library(sp)
library(lattice)
library(latticeExtra)
library(rgdal)
filename <- "C:/Users/MidwoodJ/Documents/Projects/LKO/Electrofishing/Regional IBI Comparisons/GIS/GLAF_LKO_Water.shp"
polygon_data <- readOGR(dsn = filename)
poly_df <- as(polygon_data, "data.frame")
panel_with_basemap <- function(...) {
  panel.xyplot(...)           # plot points (your xyplot data)
  sp.polygons(poly_df, col = "lightgrey")  # add the basemap
}

plot(poly_df[2])
points

# Option 1:
par(mfrow=c(1,2))
MyCex <- 3 * abs(E1) / max(E1) + 0.5
Sign <- as.numeric(E1 >=0) + 1
MyPch <- c(1, 3)[Sign]
MyCol <- c("black", "red")[Sign]
xyplot(M1LAT~ M1LONG,
       data = df2,
       panel = panel_with_basemap,
       cex = MyCex,
       pch = MyPch,
       col = MyCol,
       aspect = "iso",
       xlab = list(label = "Latitude", cex = 1.5),
       ylab = list(label = "Longitude", cex = 1.5)
) # hard to see, but looks like all good
# SAVE AS: FishPA_BestModel_N214_Spatial_FigureS5_05Apr2024

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
# And of course now you need to apply model validation

#####################################
## Predict model output (plotting) ##
#####################################
X <- data.frame(
  Intercept = rep(1,N),
  fG2 = Xm[,2],
  fG3 = Xm[,3], 
  fG4 = Xm[,4],
  fG5 = Xm[,5],
  fG6 = Xm[,6],
  fG7 = Xm[,7],
  WetlandDist.std =  Xm[,8],
  RiverDist.std = Xm[,9],
  Fetch.std = Xm[,10],
  Temp.std = Xm[,11],
  ClosestDist.std = Xm[,12],
  MarinaDist.std = Xm[,13],
  Slope.std = Xm[,14]
)

head(X)

# Stack for fitting model
StackFit.model <- inla.stack(
  tag = "Fit",
  data = list(AdjIbi01 = df2$AdjIbi01),  
  A = list(1, 1,A4),                  
  effects = list(  
    Intercept = rep(1,nrow(df2)),  # defining intercept 
    X    = X, # covariates
    INLA_dist = 1:length(Knots)) # spatial component
)

# Build Stack 
# stack for predicting year trend
range(df2$MarinaDist.std) #-1.065782  2.723567
range(df2$Slope.std) # -6.522487  1.111866

## Marina Distance Plot
MarinaDist.x <- seq(from = -1.1, to = 2.8, length = 25) # MarinaDist values to make prediction for
Slope.x <- rep(0,length(MarinaDist.x))  # Slope values to make prediction for

Intercept.pred <-rep(1, length(MarinaDist.x))

StackFit.pred <- inla.stack(
  tag = "Pred.Marina",                   # - Name (nametag) of the stack
  data = list(AdjIbi01 = NA),                 # - dependent variable - set to NA to make prediction
  A = list(1, 1, 1),                      # 
  effects = list(                      # - model componetns
    Intercept = Intercept.pred, # - intercept - vector of 1s
    MarinaDist.std = MarinaDist.x,
    Slope.std = Slope.x
  )
)

# Combine stacks
StackCombo<- inla.stack(StackFit.model, StackFit.pred) # stack for spatial model

I.pred <- inla(f1a, 
               control.compute = list(config = TRUE,dic = TRUE,waic = TRUE),
               control.predictor = list(A = inla.stack.A(StackCombo),
                                        compute = TRUE, 
                                        link = 1),
               quantiles = c(0.025, 0.975), 
               family = "binomial",
               data = inla.stack.data(StackCombo))



Index <- inla.stack.index(StackCombo, tag = "Pred.Marina")$data
modI4 = I.pred
pred.glm <- as.data.frame(list(
  x = MarinaDist.x * sd(df2$MarinaDist_ALL) + mean(df2$MarinaDist_ALL), ## convert back to normal units?
  mean = modI4$summary.fitted.values[Index, "mean"],
  LCI = modI4$summary.fitted.values[Index, "0.025quant"],
  UCI = modI4$summary.fitted.values[Index, "0.975quant"],
  type = "GLM"
))

par(mar=c(5,5,4,1))
p <- ggplot()
p <- p + geom_point(data = df2, 
                    aes(y = AdjIbi01, x = MarinaDist_ALL/1000),
                    shape = 1, 
                    size = 3)
p <- p + xlab("Harbour Distance (km)") + ylab("Nearshore Presence/Absence")
p <- p + theme_bw(base_size = 25)
p <- p + geom_line(data = pred.glm, 
                   aes(x = x/1000, y = mean), 
                   colour = "black")
p <- p + geom_ribbon(data = pred.glm, 
                     aes(x = x/1000, 
                         ymax = UCI, 
                         ymin = LCI),
                     alpha = 0.4)
p
## SAVE AS: FishPA_BestModel_N214_Covariate_Marina_Figure_05Apr2024
png("FishPA_BestModel_N214_Covariate_Marina_Figure_18Sept2024.png",
    width = 2400, height = 2400,units="px",res=300)
p
dev.off()



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
             family = "binomial",
             lincomb = lc,
             #control.inla = list(lincomb.derived.only = F),
             data = inla.stack.data(StackFit))


summary(I.lc)
HyperPar.I.lc.pc <- bri.hyperpar.summary(I.lc)
round(HyperPar.I.lc.pc, digits = 3)

LC <- I.lc$summary.lincomb.derived[,c("mean", "0.025quant", "0.975quant")]
round(LC, digits = 3)
## SAVE

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
                  family = "binomial",
                  lincomb = lcpred,
                  #control.inla = list(lincomb.derived.only = F),
                  data = inla.stack.data(StackFit))


summary(I.lc.plot)
HyperPar.I.lc.plot.pc <- bri.hyperpar.summary(I.lc.plot)
round(HyperPar.I.lc.plot.pc, digits = 3)

LC <- I.lc.plot$summary.lincomb.derived[,c("mean", "0.025quant", "0.975quant")]
round(LC, digits = 3)
LC$Treatment<-c("Gravel Beach","Groyne","Mixed shoreline","Mixed Beach","Revetment","Surcharged Revetment","Sand Beach")
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
p <- p + xlab("Treatment") + ylab("Presence / Absence")
p <- p + theme(text = element_text(size = 25)) 
p <- p + theme(legend.position="none") 
p
## SAVE AS: FishPA_BestModel_N214_FigureS2_05Apr2024


#I9c$marginals.lincomb$lcpred1
plot(I.lc.plot$marginals.lincomb.derived$lcpredI,type="l",xlim=c(-4,4),
     ylim=c(0,2.1),lwd=2,lty=2,col="grey",
     xlab="Posterior Distribution",ylab="Probability",main="Probability of Non-Pelagic Fish")
lines(I.lc.plot$marginals.lincomb.derived$lcpred1,lwd=2,col="black")
lines(I.lc.plot$marginals.lincomb.derived$lcpred2,lwd=2,col="blue")
lines(I.lc.plot$marginals.lincomb.derived$lcpred3,lwd=2,lty=2,col="blue")
lines(I.lc.plot$marginals.lincomb.derived$lcpred4,lwd=2,lty=1,col="brown")
lines(I.lc.plot$marginals.lincomb.derived$lcpred5,lwd=2,lty=1,col="red")
lines(I.lc.plot$marginals.lincomb.derived$lcpred6,lwd=2,lty=2,col="yellow")


png(file="FishPA_BestModel_N214_Lincomb_Figure_03Apr2024.png", res = 300, width = 9.5, height = 8, units = 'in')
par(mar=c(5,6,4,1))
plot(I.lc.plot$marginals.lincomb.derived$lcpredI,type="l",xlim=c(-4,4),
     ylim=c(0,2.1),lwd=3,lty=2,col="grey23",
     xlab="Posterior Distribution",ylab="Probability",main="Probability of Non-Pelagic Fish",cex.lab=2, cex.axis=1.5, cex.main=2, cex.sub=1.5)
lines(I.lc.plot$marginals.lincomb.derived$lcpred1,lwd=3,col="black")
lines(I.lc.plot$marginals.lincomb.derived$lcpred2,lwd=3,col="cadetblue4")
lines(I.lc.plot$marginals.lincomb.derived$lcpred3,lwd=3,lty=2,col="blue")
lines(I.lc.plot$marginals.lincomb.derived$lcpred4,lwd=3,lty=1,col="brown")
lines(I.lc.plot$marginals.lincomb.derived$lcpred5,lwd=3,lty=1,col="red")
lines(I.lc.plot$marginals.lincomb.derived$lcpred6,lwd=3,lty=2,col="darkgoldenrod2")
legend(-4,2, # places a legend at the appropriate place
       title="Shoreline Type",
       c("Gravel Beach - ab", "Groyne - ab","Mixed Shoreline - a","Mixed Beach - ab","Revetment - ab","Surcharged Revetment - b","Sand Beach - ab"), # puts text in the legend
       lty=c(2,1,1,2,1,1,2), # gives the legend appropriate symbols (lines)
       lwd=c(3.5),
       col=c("grey23","black","cadetblue4","blue","brown","red","darkgoldenrod2"),
       cex=1.25) # gives the legend lines the correct color and width
dev.off()

