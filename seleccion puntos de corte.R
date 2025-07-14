###################################################################################
#####   Codigo para hacer validacion cruzada al punto de corte de seleccion   #####
####                                                                           ####
###################################################################################

install.packages("Metrics")  
library(Metrics) 

source("~/Creacion de bases/BayesSegmentation_functions.R")

# Seleccion de los punto de cambio con el umbral de la probabilidad posteriori para todas las series
breakpointsi <- lapply(ML, function(x) c(which((x)/(itertot-burnin) > 0.01734))-1)

# Seleccion de las funciones para varias series de tiepo con el umbral de la probabilidad posterior
basefunctionsi  <- lapply(ML2, function(x)  c(which(x/(itertot-burnin) > 0.005036)))

# coordenadas de los puntos de cambio para cada serie
nbselectionsi <- mapply(function(x,y) return(x[y+1]), x = ML, y = breakpointsi)

# Estimacion de gamma para todas las series 
gammahati <- replicate(89, numeric(365), simplify = FALSE)
for (i in 1:89) {
  for (j in (as.numeric(unlist(breakpointsi[i])) + 1)){gammahati[[i]][j] = 1}}

# coordenadas de las funciones en todas las series
nbselections2i  <- mapply(function(x,y) return(x[y]), x = ML2, y = basefunctionsi)
# Nombre de las funciones seleccionadas 
namesfun <- lapply(basefunctionsi, function(x)  colnames(Fmatrix)[x])

## Estimacion de r para todas las series 
rhati <- replicate(89, numeric(dim(Fmatrix)[2]), simplify = FALSE)
for (i in 1:89) {
  for (j in (as.numeric(unlist(basefunctionsi[i]))[-1])){rhati[[i]][j] = 1}}



# Calculo de la funcion muestreo de gibbs para varias series 
resultados <- mapply(
  estimation_moy_biais,
  serie = cp1i,
  lec1 = lec1i,
  lec2 = lec2i,
  gammahat = gammahati,
  rhat = rhati,
  MoreArgs = list(nbiter = itertot, nburn = burnin,Fmatrix = Fmatrix, priorminsigma2 = priorminsigma2, priormaxsigma2 = priormaxsigma2,printiter=FALSE),
  SIMPLIFY = FALSE)



######---------    Funcion para reconstruir la curva de luz    --------############
Reconstrucion <- function(breakpointsi,resultados,Fmatrix,basefunctionsi){
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


# Se calculan los Y estimados o predichos, se reconstruye la curva de luz
reconstructiontoti <- list()
for (i in 1:length(reconstructionmui)){
  reconstructiontoti[[i]] <- unlist(reconstructionmui[i]) + unlist(reconstructionfi[i])
}
return(reconstructiontoti )}
#


#######------ Funcion para calcular la media de los errores cuadrados medios -------########
MediaECM <- function(curvasLL,reconstructiontoti){
# Calcular MSE
ErrorCM <- list()
for (i in 1:length(curvasLL)){
  ErrorCM[[i]] <- mse(unlist(curvasLL[i]), unlist(reconstructiontoti[i]))
}
ErrorCMM <- mean(as.vector(unlist(ErrorCM))) 
rmse <- sqrt(ErrorCMM)
 return (c(ErrorCMM,rmse))}
#
################################################################################
is.na(reconstructiontoti)

reconstrucionto4 <- Reconstrucion(breakpointsi,resultados,Fmatrix,basefunctionsi)

MediaECM(cp1i, reconstrucionto4); MediaECM(curvasLL, reconstrucionto1)
str(reconstructiontoti)
