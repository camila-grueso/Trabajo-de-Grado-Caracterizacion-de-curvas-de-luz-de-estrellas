###############################################################################
#######  Modelo de Segmentacion de series temporales con aproximacion  ########
######           Bayesiana, para series con efecto funcional            #######
#####                 Para multiples Series temporales                   ######
###############################################################################

install.packages("dplyr")
library(dplyr)
library(readr)
library(zoo)

source("~/Creacion de bases/bases_gen.R")
source("~/Creacion de bases/BayesSegmentation_functions.R")

curvas <- read.csv("cp1i.csv")
curvas <- as.list(curvas)

### Suavizacion de las curvas de luz plegadas
curvas_S <- list()

for (i in 1:length(curvas)) {
  serie <- as.ts(unlist(curvas[[i]]))
  suavizada <- rollmean(serie, k=20,fill = NA )
  magnitud_suavizada <- na.approx(suavizada,rule=2)
  curvas_S[[i]] <- as.vector(magnitud_suavizada)
}


curvas_S[[1]]

plot.ts(curvas_S[[1]])

# creacion de la F matriz y sus funciones segun el tamaño de la serie
t <- (1:length(curvas_S[[1]]))
Fmatrix<-CreationBases(t,7,7,10)[[1]]
Fmatrix
dim(Fmatrix)
ti <- replicate(67, 1:length(curvasL[,1]), simplify = FALSE)
m1 <- na.omit(curvas_S[[1]])
m2 <- na.omit(curvas_S[[2]])

# Calculo de la F matriz con series de diferentes tamaños
#FM <- lapply(ti, function(x) CreationBases(x,7,7,10))
#FM[[40]]$F.Bases


### Funcion Metropolis H, para estimar betagamma landar

# Se establecen los palametros para realizar la funcion de metropolis hasting
itertot <- 100000
burnin <- 35000
lec1 <- 50
lec2 <- 50
lec1i <- rep(50, 67)
lec2i <- rep(50, 67)
nbseginit <- 5 
nbfuncinit <- 6
nbtochangegamma <- 1 
nbtochanger <- 1 
Pi <- c(1,rep(0.01,length(m1)))   
eta <- c(1,rep(0.01,dim(Fmatrix)[2]-1))

#######------ Calculo del Metropolis Hasting para una serie -----------####

MH <- segmentation_bias_MH(m1, itertot, burnin,lec1,lec2, Fmatrix, nbseginit,nbtochangegamma,nbfuncinit,nbtochanger,Pi,eta,printiter=FALSE) 

#####---- Calculo del Metropolis Hasting para varias series de tiempo ----####
MMH <- lapply(curvas_S, function(x) segmentation_bias_MH(x, itertot, burnin,lec1,lec2, Fmatrix, nbseginit,nbtochangegamma,nbfuncinit,nbtochanger,Pi,eta,printiter=FALSE))

# begin 12:08 pm 23 de junio 2025 metropoli - 5:26 pm
MMH[[1]][5]

######--- Seleccion de los puntos de cambio y las funciones para una sola serie ---#### 

# se obtiene la estimacion de Beta gamma y lambda R
ML <- as.numeric(unlist(MMH[[4]][1]))
ML2 <- as.numeric(unlist(MMH[[4]][2]))

PC <- unique(round(ML/(itertot-burnin),7)) # codigo para elegir la probabilidad posterior para validacion cruzada threshold
quantile(PC, 0.9)

PCF <- unique(round(ML2/(itertot-burnin),7)) # codigo para elegir la probabilidad posterior para validacion cruzada threshold
quantile(PCF, 0.9)

threshold <-0.018  # umbral de las probabilidades posteriores para los puntos de cambio
str(breakpoints)
#Seleccion de los puntos de cambio
breakpoints <- c(which(ML/(itertot-burnin) > threshold))-1
u <- as.vector(as.numeric(unlist(breakpoints[1])))
nbselections <- ML[[1]][u+1]  #Numero de Puntos de cambio seleccionados
cbind(breakpoints,nbselections)
gammahat <- rep(0,length(m1)) #gamma estimado 
for (i in (breakpoints+1)){gammahat[i] <- 1}


PC <- unique(round(ML2/(itertot-burnin),7)) # codigo para elegir la probabilidad posterior para validacion cruzada threshold
quantile(PC, 0.9)  
# Umbral de las probabilidas posteriores para las funciones          
threshold <- 0.0
#Seleccion de las funciones 
basefunctions <- c(which(ML2/(itertot-burnin) > threshold))
nbselections2 <- ML2[basefunctions] # numero de funciones seleccionadas
cbind(basefunctions,nbselections2,colnames(Fmatrix)[basefunctions])
rhat <- rep(0,dim(Fmatrix)[2])  # r estimado
for (i in basefunctions[-1]){rhat[i] <- 1}

#######---- Seleccion de los puntos de cambio y las funciones para varias series ----####

# Se obtiene la estimacion de beta gamma y lambda r en una lista para todas las series 
ML <- list()
ML2 <- list()
for (i in 1:length(MMH)){
  
  ML[i] = MMH[[i]][1]
  ML2[i] = MMH[[i]][2]}

# Seleccion de los punto de cambio con el umbral de la probabilidad posteriori para todas las series
breakpointsi <- lapply(ML, function(x) c(which((x)/(itertot-burnin) > 0.025078))-1)
# coordenadas de los puntos de cambio para cada serie
nbselectionsi <- mapply(function(x,y) return(x[y+1]), x = ML, y = breakpointsi)
ML[1]
nbselectionsi[1]

# Estimacion de gamma para todas las series 
gammahati <- replicate(89, numeric(365), simplify = FALSE)
for (i in 1:89) {
   for (j in (as.numeric(unlist(breakpointsi[i])) + 1)){gammahati[[i]][j] = 1}}

as.numeric(unlist(breakpointsi[1])) + 1
gammahati[1]

# Seleccion de las funciones para varias series de tiepo con el umbral de la probabilidad posterior
basefunctionsi  <- lapply(ML2, function(x)  c(which(x/(itertot-burnin) > 0.014935)))
# coordenadas de las funciones en todas las series
nbselections2i  <- mapply(function(x,y) return(x[y]), x = ML2, y = basefunctionsi)
# Nombre de las funciones seleccionadas 
namesfun <- lapply(basefunctionsi, function(x)  colnames(Fmatrix)[x])

## Estimacion de r para todas las series 
rhati <- replicate(89, numeric(dim(Fmatrix)[2]), simplify = FALSE)
for (i in 1:89) {
  for (j in (as.numeric(unlist(basefunctionsi[i]))[-1])){rhati[[i]][j] = 1}}

as.numeric(unlist(basefunctionsi[1]))[-1]
rhati[2]

#######----- Calculo del algoritmo de muestreo de gibbs ------#######

# Once the breakpoints and functions have beeen selected, betagamma, lambdar et sigma2 can be estimated, 

# Parametros de la funcion que realiza el muestreo de gibbs 
priorminsigma2 <- 0.001
priormaxsigma2 <- 5

# calculo de la funcion muestreo de gibbs para una sola serie 
estim <- estimation_moy_biais(curvasL[,1], itertot, burnin, lec1, lec2, Fmatrix, unlist(gammahati[1]),unlist(rhati[1]),priorminsigma2,priormaxsigma2,printiter=FALSE)

# Calculo de la funcion muestreo de gibbs para varias series 
resultados <- mapply(
  estimation_moy_biais,
  serie = curvas_S,
  lec1 = lec1i,
  lec2 = lec2i,
  gammahat = gammahati,
  rhat = rhati,
  MoreArgs = list(nbiter = itertot, nburn = burnin,Fmatrix = Fmatrix, priorminsigma2 = priorminsigma2, priormaxsigma2 = priormaxsigma2,printiter=FALSE),
  SIMPLIFY = FALSE)

str(rhati)
estimi[[1]][6]
str(estimi)
resultados[[38]][4]
####---- Reconstruccion de la serie temporal con los puntos de cambio y las funciones ----####

#######---- Reconstruccion de la curva de luz una sola serie -----#######

# vector donde almacenan los puntos de cambios 
muest <- rep(0,length(unlist(breakpointsi[1])))
muest[1] <- unlist(resultados[[1]][4])[1] # almacena la media de los segmentos

# Reconstrucion de la serie de acuerdo a los puntos de cambio y la media presente en los segmentos
reconstructionmu <- rep(0,dim(Fmatrix)[1])
if (length(breakpoints)>1){
  for (i in 2:length(breakpoints)){muest[i] <- unlist(estim[[4]])[i]+muest[i-1]}
  compt <- 1
  for (i in 2:length(c(breakpoints,dim(Fmatrix)[1]))){
    reconstructionmu[compt:(c(breakpoints,dim(Fmatrix)[1])[i])] <- muest[i-1]
    compt <- compt + (c(breakpoints,dim(Fmatrix)[1])[i]-c(breakpoints,dim(Fmatrix)[1])[i-1])
  }
}


# Reconstrucion de la serie de acuerdo con las funciones presentes en cada segmento 
reconstructionf <- rep(0,dim(Fmatrix)[1])
if (length(basefunctions)>1){
  for (i in 2:length(basefunctions)){
    reconstructionf <- reconstructionf + estim[[5]][i-1]*Fmatrix[,basefunctions[i]]
  }   
}

# reconstrucion total = reconstrucion de la media mas la reconstrucion de las funciones 
reconstructiontot <- reconstructionmu + reconstructionf

# grafico reconstrucion de las funciones 
plot.ts(reconstructionf)
plot.ts(m1)


##########----- Reconstruccion de la curva de luz de varias series ------#######

# Se guarda un vector de los puntos de cambio de todas las series en una lista
muesti <-list()
for (i in 1:length(breakpointsi)){
muesti[i] <- rep(0,length(unlist(breakpointsi[i])))
 for (j in 1:length(unlist(breakpointsi[i]))) {
     muesti[[i]][1] <- unlist(resultados[[i]][4])[1] }} # almacena la media de los segmentos de cada serie


# Reconstrucion de la serie de acuerdo a los puntos de cambio y la media presente en los segmentos
reconstructionmui <- replicate(89, numeric(dim(Fmatrix)[1]), simplify = FALSE)
compti <- list()

for (i in  1:length(breakpointsi)) {

if (length(unlist(breakpointsi[i]))>1){
  for (j in 2:length(unlist(breakpointsi[i]))){muesti[[i]][j] <-  unlist(resultados[[i]][4])[j]+muesti[[i]][j-1]}
  compti[i] <- 1
  for (j in 2:length(c(unlist(breakpointsi[i]),dim(Fmatrix)[1]))){
    reconstructionmui[[i]][unlist(compti[i]):(c(unlist(breakpointsi[[i]]),dim(Fmatrix)[1])[j])] <- muesti[[i]][j-1]
    compti[i] <- unlist(compti[i]) + (c(unlist(breakpointsi[[i]]),dim(Fmatrix)[1])[j]-c(unlist(breakpointsi[i]),dim(Fmatrix)[1])[j-1])
  }}}

# Reconstrucion de cada serie de acuerdo con las funciones presentes en cada segmento 
reconstructionfi <- replicate(89, numeric(dim(Fmatrix)[1]), simplify = FALSE)

for(i in 1:length(basefunctionsi)) {
  
if (length(unlist(basefunctionsi[i]))>1){
  for (j in 2:length(unlist(basefunctionsi[i]))){
    reconstructionfi[[i]] <- unlist(reconstructionfi[i]) + unlist(resultados[[i]][5])[j-1]*Fmatrix[,unlist(basefunctionsi[[i]])[j]]
  }}}

# reconstrucion total igual a la reconstrucion de la media de los puntos de cambio y las funciones para cada serie
reconstructiontoti <- unlist(reconstructionmui[2]) + unlist(reconstructionfi[2])
