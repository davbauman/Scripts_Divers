library(gstat)
library(lattice)
library(sp)

# Importation du fichier contenant les coordonnées spatiales des placeaux
# et les variables environnementales associées (transformées par Cox-Box)

mydata <- read.table("topography.txt", h = T, sep = "\t", row.names = 1)
colnames(mydata)

# Sélectionner seulement la variable d'intérêt et les coord. spa. pour ne pas éliminer
# de lignes inutilement en éliminant les NA:

mydata = mydata[, c(1,2,8)]
mydata <- na.omit(mydata)

# Pour Treignes uniquement:
# *************************
piquets <- read.table("Piquets.txt", h = T, sep = "\t", row.names = 1)
colnames(piquets) <- c("x", "y")
coordinates(piquets) <- ~ x + y

# Convert the (full and subset) data frame(s), and the interpolation
# grid, to sp classes. They can later be recovered as data.frames with 
# "as.data.frame()". And create the dataset positions to be predicted.

grd <- expand.grid(x = seq(0, 417, by = 3), y = seq(0, 84, by = 3)) 
# x = seq(z/2, y, by = z) serait la formule générale

# seq() va du min au max par pas de 'by', mais le premier pas = "point 0", il faut
# donc ajouter un pas de 3 au max pour y arriver, sinon on s'arrêtera à max-3.

coordinates(grd) <- ~ x + y
gridded(grd) <- TRUE
plot(grd, cex = 0.5)
points(mydata, pch = 1, col = 'red', cex = 0.7)
points(piquets, pch = 1,col = 'black', cex = 0.7)

coordinates(mydata) <- ~ x + y

# Compute the omnidirectional spatial variogram of the wanted variable and
# fit a variogram model to it.

# First we examine the variogram cloud to see all the point-
# pairs; then we compute the experimental variogram. Then we specify a
# variogram model and its initial parameters (determined by eye looking
# at the experimental variogram); Finally we fit it.

# Create gstat objects (objects that hold all the information necessary 
# for geostatistical prediction

g <- gstat(id = "Flow.direction", formula = Flow.direction ~ 1, data = mydata)

# Experimental variogam
# Argument 'cutoff = ' précise jusqu'à quelle distance considérer les paires de
# points pour l'estimation du variogramme.
# Argument 'width' précise, au sein de la gamme allant de 0 à cutoff, les intervalles
# réguliers de distance pour lesquelles les points seront regroupés.

v <- variogram(g)
plot(v, pl = T)

# ,boundaries=c(0,25,50,75,100,125,150,175,200,225,250,275,300)

# Estimate variogram model form and parameters by eye:
# ****************************************************

m.var1 <- vgm(1.5 , 'Sph', 40, 0)  
# Arg1= psil; Arg2= Mod.; Arg3= Height at which equilibrium is reached; Arg4= Nugget (if any)

plot(v, pl=T, model=m.var1)

# fit model parameters by weighted least-squares:
# ***********************************************

(v.fit <- fit.variogram(v, m.var1))
plot(v, pl = T, model = v.fit, as.table = TRUE)

# Use the fitted model of regionalization to interpolate the variable
# with OK on the prediction grid.

g <- gstat(g, id = "Flow.direction", model = v.fit)
p <- predict(g, model = v.fit, newdata = grd) 
 
# remplacer 'grd' par un objet contenant uniquement les coord. des points pour 
# lesquels voulons des valeurs de la variable, au cas où ne voulons pas toute 
# la grille mais qu'une fenêtre par ex.

# RGB Representation of the kriging values
# ****************************************

# Back transformation of the Box-Cox values:

p2 <- as.data.frame(p)
lambda <- 0.46
p2[, 3] <- ((p2[, 3] * lambda)+1)^(1/lambda)
coordinates(p2) <- ~x+y
gridded(p2) <- TRUE

myPal<- colorRampPalette(c("goldenrod", "darkred"))
image(p, col=myPal(30))
contour(p, add=TRUE, drawlabels=FALSE, col='black')
points(mydata, pch=4, col = 2, cex=0.5)
points(piquets, pch=4, col = 1, cex=0.5)
pts <- list("sp.points", mydata, pch = 4, col = "black", cex=0.5)
spplot(p2, zcol="Flow.direction.pred", col.regions= myPal(30), cuts=19, sp.layout=list(pts), 
       contour=TRUE, labels=FALSE, pretty=TRUE, col='brown')

# Enregistrement des résultats 
# ****************************

z=data.frame(p)
write.table(z, "Results.Flow.Direction.txt", append = F, sep = "\t")

# Cross-validation
# ****************

source("ck_plotfns.R")
# plot.kresults(p, "Slope", mydata, "Slope")

# Cross-validate the OK predictions of the target-variable and compute
# the diagnostic measures.

cv.o <- krige.cv(Flow.direction ~ 1, mydata, model = v.fit, nfold = nrow(mydata))
summary(cv.o)
res <- as.data.frame(cv.o)$residual
# RMSE (precision) --> Ideally small value
sqrt(mean(res^2))
# Mean Error (of the residuals): the bias --> Ideally 0
mean(res)
# Mean Squared Deviation Ratio (MSDR) : If < 1 --> predictions are less
# variables than reality (: to be expected, since kriging is a smoothing procedure).
mean(res^2/as.data.frame(cv.o)$var1.var)
rm(res)

# Pour connaître l'ampleur absolue de l'erreur sur les prédictions, il faut
# retransformer en valeurs non transformées et en faire de même pour l'erreur.
# L'erreur est à additionner aux valeurs prédites (et pas à considérer à part).
# (v. p. 27 dans Rossiter 2009).















######################################################################
######################################################################
######################################################################



# Modelling a bivariate co-regionalisation (co-kriging):
########################################################

# Display a scatterplot of the co-variable vs. the target variable and compute
# the correlation coefficient. Describe the feature-space relation.

attach(as.data.frame(mydata))
xyplot(Carbo ~ pH_H2O, pch = 20, cex = 1.2, col = "blue", ylab = "Carbo", 
   xlab = "pH_H2O")
cor(pH_H2O, Carbo)
sum(is.na(pH_H2O))
(cor <- cor(pH_H2O, Carbo, use = "complete"))   # 'complete' --> cor que pour non NA.
cor^2   # R2

# Is it sufficient for a successful co-kriging?
# We model the omnidirectional spatial variogram of the covariable. We do not
# intend to interpolate it; however, we need to know its spatial structure
# for co-kriging.
# Si nécessaire, faire un subset de la covariable pour éviter les NA.
# var.co <- subset(as.data.frame(mydata), !is.na(om), c(x, y, "nom de var"))

var.co <- subset(as.data.frame(mydata), !is.na(Carbo), c(x, y, Carbo))
coordinates(var.co) <- ~ x + y
v.co <- variogram(Carbo ~ 1, var.co)
plot(v.co, pl=T)
# model by eye
m.co <- vgm(.2, "Sph", 40, .4)
# fit
(m.co.fit <- fit.variogram(v.co, m.co))
plot(v.co, pl=T, model=m.co.fit)

# Compare variogram structure to target variable

var.target <- subset(as.data.frame(mydata), !is.na(pH_H2O), c(x, y, pH_H2O))
coordinates(var.target) <- ~ x + y
v.tar <- variogram(pH_H2O ~ 1, var.target)
plot(v.tar, pl=T)
# model by eye
m.tar <- vgm(3.2e+07, "Sph", 40, 6e+07)
# fit
(m.tar.fit <- fit.variogram(v.tar, m.tar))
plot(v.tar, pl=T, model=m.tar.fit)

# Comparons les 'ranges'
m.co.fit$range[2]; m.tar.fit$range[2]

round(m.co.fit$psill[1]/sum(m.co.fit$psill),2)
round(m.tar.fit$psill[1]/sum(m.tar.fit$psill),2)

# The linear model of co-regionalisation requires the ranges of the target
# and co-variable be the same; here we see they are not too different.

# Building a data structure to model co-regionalisation:
# ******************************************************

# Now comes the hardest part of co-kriging: modelling the co-regionalisation.
# We have to fit models to both the direct and cross-variograms simultaneously,
# and these models must lead to a positive definite cokriging system. The 
# easiest way to ensure this is to fit a linear model of co-regionalisation: 
# all models (direct and cross) have the SAME SHAPE AND RANGE,
# BUT MAY HAVE DIFFERENT PARTIAL SILLS AND NUGGETS.

# VOIR PAGE 31 de Rossiter pour explications step-by-step de ce que fait gstat.

# First frame (of the list) = target var., second frame = covariable.

(g <- gstat(NULL, id = "pH_H2O", form = pH_H2O ~ 1, data = var.target))
(g <- gstat(g, id = "Carbo", form = Carbo ~ 1, data=var.co))

# Compute and display the two direct variograms and one cross-variogram.

v.cross <- variogram(g)
str(v.cross)
plot(v.cross, pl=T)

# Fitting a linear model of co-regionalisation:
# *********************************************

# The next step is to add variogram models to the gstat object and then fit these.
# We use the fitted model for the target value as a starting point for all
# three variogram models. To use the linear model of co-regionalisation, 
# --> single range and structure.
# By filling all the frames with one model, these conditions are automat. met.

(g <- gstat(g, id = "pH_H2O", model = m.tar.fit, fill.all=T))

# The gstat object contains now both the data and the models.

# We fit the the three variograms together, ensuring they lead to a positive-
# definite co-kriging system.

(g <- fit.lmc(v.cross, g))
plot(variogram(g), model=g$model)

# We can compare the 2 direct variograms with the parameters defined for each
# one separately (pour voir à quel point ajuster le modèle de co-régionalisation
# à ajusté le fit à chaque variable séparément).

str(m.co.fit)
str(g, max.level = 1)

str(g$data, max.level = 1)
str(g$model, max.level = 1)
str(g$data$pH_H2O, max.level = 1)
str(g$model$pH_H2O, max.level = 1)

# Comparaison
# ***********
# Covariate
g$model$Carbo$psill - m.co.fit$psill
sum(g$model$Carbo$psill) - sum(m.co.fit$psill)
sum(g$model$Carbo$psill)   # Total sill

# Target variable
g$model$pH_H2O$psill - m.tar.fit$psill
sum(g$model$pH_H2O$psill) - sum(m.tar.fit$psill)
sum(g$model$pH_H2O$psill)   # Total sill

# Les deux dernières lignes (pour covariate et target) sont censées être proches.

##### CO-KRIGING WITH ONE CO-VARIABLE
#####################################

# The wrapper method krige that was used for OK (§5) can only be used for 
# univariate kriging; here we must use the predict.gstat method. This takes a 
# gstat object as the ?rst argument and the prediction points data frame as 
# the second argument.

# interpolate
k.c <- predict.gstat(g, grd)
str(k.c)
# summarize predictions and their errors
summary(k.c$pH_H2O.pred); summary(k.c$pH_H2O.var)

# Display the predictions and their errors as maps:
plot.kresults(k.c, "Carbo", var.co, var.target, "pH_H2O", 
   "CK with Carbo covariable, ")

# Cross-Validation
# ################

cv.c <- gstat.cv(g)
str(cv.c)
summary(cv.c$residual)
sqrt(mean(cv.c$residual^2))   # RMSE (mean error)
mean(cv.c$residual)   # MSDR (biais)
mean(cv.c$residual^2/cv.c$ltpb.var)


