########     Codigo Plegado en phase de las curvas de luz           #########
#######         y Seleccion de la muestra aleatoria > 365            ########
#############################################################################

library(readr)
library(dplyr)
# Cargar los datos de la magnitud de la estrella 
Datos <-read_csv("DatosStar_EB_Con")

# sub base de datos donde se almacenan los dias julianos de las mediciones de cada estrella 
Datosmag1 = Datos[,1:9]

# Muestra aleatoria de las estrellas que se van a analizar 
set.seed(202524)

muestra <- sort(sample(1:47735, 150))

# Data frame de las estrellas seleccionadas 
df <- data.frame()
for (i in 1:length(muestra)){
t <- Datosmag1[Datosmag1$Estrella==muestra[i],]
 df <- rbind(df,t)  
}

# Este codigo es para enumerar nuevamente a las estrellas 
na <- is.na(df[,1])
na[FALSE,] = 0

Estrella = 1
for (i in 1:(nrow(df))) {
  
  if (na[i] == 0){
    df$Estrella[i] = Estrella
  }
  
  if (na[i] == 1){
    df$Estrella[i] = Estrella + 1
    Estrella = Estrella + 1 
  }
}

### Transformacion de los datos pre procesamiento 
df$P = round(as.numeric(df$P),5)
# Eliminar el signo igual
df$To <- gsub("=", "", df$To)
# Convertir a numérico
df$To <- as.numeric(df$To)
df$Mag_Mean <- as.numeric(df$Mag_Mean)
df$Tipo <- as.factor(df$Tipo)
str(df)
summary(df)

89 * (30490/(135675 + 30490))

na.omit(df)
BaseEb <- data.frame()
for (i in 1:150) {
  BaseEb[i,] = df[df$Estrella==1,][1,]
}
table(tipoE)

# formula para hallar la fase  
fase = (((df$HJD - df$To) /df$P) %% 1) 

plot(df$fase[540:900],df$Magn_Inf[540:900])

df$fase = round(fase,6)

######## Lista del calculo de la phase para cada curva de luz #######
fasei <- list(0)
for (i in 1:length(tiempi)){
  for (j in 1:length(tiempi[[i]])) {
    fasei[[i]][j] = ((((tiempi[[i]][j]) - t[i]) /p[i]) %% 1) + 1
  }
}


####################################################################

# Seleccionar las fases con datos superiores a 365 dias julianos 
faselis <- list()
Datosfase = df[,c(10,2,4)]
Datosfase = na.omit(Datosfase)

for (i in 1:length(muestra)){
  faselis[i] = Datosfase[Datosfase$Estrella==i,]
}

fase365 <- list()
for (i in 1:length(faselis)) {
  U <- na.omit(unlist(faselis[i]))
  if (length(U) >= 365) {
    fase365[i] = faselis[i]
  }}

fase365 <- Filter(Negate(is.null), fase365)

#Convertir la lista de datos en un data frame 
phases <- bind_rows(lapply(fase365, function(x) as.data.frame(t(x))))

phases = phases[,-(366:12717)]  # eliminar los datos superiores a 365 dias
phasesL <- as.data.frame(t(phases))

# pasar el data.frame a lista 
phasesLL <- list()
for (i in 1:length(phases[,1])){
  m1 <-  phasesL[i]
  phasesLL[i] = m1}


# Seleccionar las magnitudes con datos superiores a 365

Datosmag = na.omit( df[,c(2,4,10)])

# Seleccionar las estrellas de la muestra aleatoria 
muestrastar <- list()
for (i in 1:length(muestra)){
  muestrastar[i] = Datosmag[Datosmag$Estrella==i,]
}

# Seleccionar las estrellas con datos superiores a 365 dias julianos 
muestrastar365 <- list()

for (i in 1:length(muestrastar)) {
  U <- na.omit(unlist(muestrastar[i]))
  if (length(U) >= 365) {
    muestrastar365[i] = muestrastar[i]
  }}

# Eliminacion de datos nulos 
muestrastar365= na.omit(muestrastar365)
muestrastar365 <- Filter(Negate(is.null), muestrastar365)

#Convertir la lista de datos en un data frame 
curvas <- bind_rows(lapply(muestrastar365, function(x) as.data.frame(t(x))))

curvas = curvas[,-(366:12718)]  # eliminar los datos superiores a 365 dias
curvasL <- as.data.frame(t(curvas))

# pasar el data.frame a lista 
curvasLL <- list()
for (i in 1:length(curvas[,1])){
  m1 <-  curvasL[i]
  curvasLL[i] = m1
}


# Ordenar las magnitudes en funcion de la phase para capturar la forma de la curva de luz 
cp1i <- list()
ordenfase <- list()
for (i in 1:length(phasesLL)){
  ordenfase[[i]] <- order(as.vector(unlist(phasesLL[i])))
}
for (i in 1:length(curvasLL)) {
    cp1i[[i]] <- curvasLL[[i]][as.vector(unlist(ordenfase[[i]]))]
}
