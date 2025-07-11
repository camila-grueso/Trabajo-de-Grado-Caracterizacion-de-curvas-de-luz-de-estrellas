#------------------------------------------------------------------------------#
#       Algoritmo cesta de mercado con las funciones curvas de luz          ####
#------------------------------------------------------------------------------#

wd="C:\\Users\\camila\\Documents\\10 Semestre\\Trabajo de grado\\Resultados\\"  # Ruta al Directorio de trabajo
setwd(wd)                                                 # Configuración al Directorio de trabajo                                     

install.packages("arules")                                # Librería para el análisis de reglas de asociación.
library("arules")
install.packages("arulesViz")
library(arulesViz);

# Lectura en formato data.frame
curvas <- read.csv("cp1i.csv")
curvas <- as.list(curvas)


funciones=read.table("TransacionFunciones.csv",header=T,sep=";")  # Lectura de Datos en R
View(funciones)
str(funciones)

Tr=arules::read.transactions("clipboard",format="basket",header=TRUE,
                             sep=",",cols =NULL ,rm.duplicates = TRUE)

?read.transactions
View(Tr)
str(Tr)
summary(Tr)  
colnames(Tr)        


## 1.2. Análisis Exploratorio de transacciones (longitud, items freceuntes)   ####
#------------------------------------------------------------------------------#

# Longitud de las transacciones
long.Tr = size(Tr)  # Tamaño de las transacciones
long.Tr=long.Tr[-90]
windows()
M=matrix(c(1,1,1,1,1,1,2,2,2), ncol=3,byrow=TRUE)
layout(M)
barplot(table(long.Tr[-90]),freq=FALSE, col="skyblue2",xlim=c(0,27),xlab="Número de funciones en la curva de luz",main="")
boxplot(long.Tr[-90], col="skyblue2",range=3,ylim=c(0,27),horizontal = TRUE)

quantile(long.Tr[-90],probs=(0:10)/10)

# Productos y su frecuencia de ventas
itemFrequency(Tr,type="absolute")   # Frecuencia/soporte de todos los items

n.item = 20   # n items (productos o combinaciones de productos) más vendidos
item_MF = sort(itemFrequency(Tr,type="absolute"),decreasing = TRUE)[1:n.item]

windows()
par(oma=c(2,3,1,2))
barplot(item_MF,horiz=TRUE,xlab="Número de Transacciones",las=1,cex.names = 0.8,
        main=paste(n.item, "Items más frecuentes"), col="skyblue2")

item_MFt <- as.matrix(item_MF)
# Extraccion de itemsets frecuentes
minsup = 0.033 

# itemset frecuentes
itemset=arules::apriori(Tr,parameter=list(support=minsup,
                                          minlen=1, maxlen=20, target = "frequent itemset"))
summary(itemset)
inspect(sort(itemset))
size(itemset)
# representación grafica de los 20 itemset más frecuentes segun longitud de itemset
L=list()      #Almacenamiento de los k-itemset, k=1,2,3,4

windows()
par(mfrow=c(2,2))
par(oma=c(2,6,1,2))
for (k in 2:3){
  h=itemset[size(itemset)==k]
  L[[k]]=as(sort(h, by = "support", decreasing = TRUE),Class = "data.frame")
  top_20_itemset=as(sort(h, by = "support", decreasing = TRUE)[1:min(10,length(h))]
                    ,Class = "data.frame")  
  barplot(top_20_itemset$count,horiz=TRUE,xlab="Número de curvas de luz",las=1,cex.names = 0.8,
          names.arg = top_20_itemset$items,main=paste("Mas Frecuente",k,"-itemset"),col = "skyblue2")
}

#####------ Filtrado de itemsets frecuentes, involucrando un producto especifico, suponiendo interés en el producto whole milk  #####

itemsets_filtrado <- arules::subset(itemset,
                                    subset = items %in% "whole milk")
inspect(sort(itemsets_filtrado, by = "support", decreasing = TRUE)[1:20])

# suponiendo interés en en la venta conjunta de whole milk and rolls/buns

itemsets_filtrado <- arules::subset(itemset,
                                    subset = items %ain% c("whole milk","rolls/buns"))
inspect(sort(itemsets_filtrado, by = "support", decreasing = TRUE)[1:min(20,length(itemsets_filtrado))])

########################################################################################
## 1.3. generación de Reglas de Asociación - Apriori                          ####
#------------------------------------------------------------------------------#
minsup = 0.033  # Al menos 2 veces por día, 30 transacciones en el mes,  
minconf= 0.5

Reglas=arules::apriori(Tr,parameter=list(support=minsup, minlen=1,maxlen=10,
                                         confidence=minconf, target = "rules"))
print(Reglas)      # Número de Reglas
summary(Reglas)    # Estadísticos resumen de las Reglas - general
inspect(Reglas)    # Medidas de asociación de las Reglas - Todas
quality(Reglas)

Reglas = sort(Reglas,by="confidence",decreasing=TRUE)  # Reordenameinto descendente por confianza
inspect(head(Reglas,20,by="confidence")) # Evaluando las 20 reglas de mayor confianza

## 1.4. Visualización y filtro de las Reglas                                           ####
#------------------------------------------------------------------------------#
windows()
plot(Reglas,measure=c("support","confidence"),shading="lift")     # grafico de dispersion - Soporte vs confianza de todas las reglas

inspect(Reglas[quality(Reglas)$confidence > 0.75]& quality(Reglas)$confidence<0.99 & quality(Reglas)$suport<0.05)
inspect(Reglas[quality(Reglas)$support > 10/nrow(Tr)])
inspect(Reglas[quality(Reglas)$support > 44/nrow(Tr) & quality(Reglas)$confidence > 0.75 ])

windows()
plot(Reglas, measure=c("support","confidence"),shading="order", control=list(main ="Diagrama de dispersión para las 86 reglas "))

windows()
plot(head(Reglas,15), method="graph", control=list(type="items"),col="blue", main="Red 15 Reglas con mas confianza") # Grafico de red 20 reglas con mayor confianza
inspect(head(Reglas,15))


windows()
plot(head(Reglas,15), method="paracoord", control=list(type="items")) # Grafico de red 20 reglas con mayor confianza

windows()
plot(head(Reglas,15), method="grouped matrix")    # Grafico de matriz de 100 regla con mayor confianza

windows()
sel = plot(Reglas, measure=c("support","confidence"), shading="lift", engine="interactive");
inspect(sel)

windows()
plot(head(Reglas,15),method = "graph",engine = "htmlwidget")

windows()
plot(head(Reglas,15), method = "graph", control = list(type = "items"))

windows()
plot.ts(unlist(curvas[[5]]), ylab="Magnitud en banda I", main = "Curva de luz plegada en fase")


##################################################
########## Graficos F matriz ###########
windows()
par(mfrow=c(2,3))

plot.ts(Fmatrix_1_$`cos(2pi1posx)`, lwd=2, main = "cos(2pi1posx)", ylab="F matrix", col="blue",cex.main = 2.5, cex.lab = 2.5,cex.axis = 2)
plot.ts(Fmatrix_1_$`sin(2pi*1posx)`,lwd=2, main = "sin(2pi1posx)", ylab="F matrix",col="blue",cex.main = 2.5, cex.lab = 2.5,cex.axis = 2)
plot.ts(Fmatrix_1_$Haar0.14 ,lwd=2, main = "Haar0.14", ylab="F matrix",col="blue",cex.main = 2.5, cex.lab = 2.5,cex.axis = 2)
plot.ts(Fmatrix_1_$bspl3.7,lwd=2, main = "bspl3.7", ylab="F matrix",col="blue",cex.main = 2.5, cex.lab = 2.5,cex.axis = 2)
plot.ts(Fmatrix_1_$x,lwd=2, main = "X", ylab="F matrix",col="blue",cex.main = 2.5, cex.lab = 2.5,cex.axis = 2)
plot.ts(Fmatrix_1_$xx,lwd=2, main = "XX", ylab="F matrix",col="blue",cex.main = 2.5, cex.lab = 2.5,cex.axis = 2)


windows()
par(mfrow=c(2,2))

plot.ts(Fmatrix_1_$`cos(2pi1posx)`, lwd=2, main = "cos(2pi1posx)", ylab="F matrix", col="blue",cex.main = 1.5, cex.lab = 2.5,cex.axis = 2)
plot.ts(Fmatrix_1_$`cos(2pi2posx)`,lwd=2, main = "cos(2pi2posx)", ylab="F matrix",col="blue",cex.main = 1.5, cex.lab = 2.5,cex.axis = 2)
plot.ts(Fmatrix_1_$Haar0.14 ,lwd=2, main = "Haar", ylab="F matrix",col="blue",cex.main = 1.5, cex.lab = 1.5,cex.axis = 2)
lines(Fmatrix_1_$Haar0.429, col = "red", lwd=2)
lines(Fmatrix_1_$Haar0.937, col = "green3", lwd=2)
lines(Fmatrix_1_$Haar1, col ="#FFA54F", lwd=2)
lines(Fmatrix_1_$Haar0.929, col = "#EEEE00", lwd=2)
lines(Fmatrix_1_$Haar0.046, col ="#00FF7F", lwd=2)
lines(Fmatrix_1_$Haar0.007, col ="#6959CD")
lines(Fmatrix_1_$Haar0.601, col = "#00F5FF", lwd=2)
abline(h=0,col = "black")

legend("topright", 
       legend=c("Haar0.429", "Haar0.937","Haar0.14","Haar1","Haar0.929","Haar0.046","Haar0.007","Haar0.601"), 
       col=c( "red","green3","blue", "#FFA54F", "#EEEE00","#00FF7F","#6959CD","#00F5FF"), 
       lwd=2, 
       cex=0.8, 
       title="Funciones")

plot.ts(Fmatrix_1_$bspl3.7,lwd=2, main = "B splines", ylab="F matrix",col="blue",cex.main = 1.5, cex.lab = 2.5,cex.axis = 2)
lines(Fmatrix_1_$bspl3.1, col="red",lwd=2)
lines(Fmatrix_1_$bspl3.3, col ="green3",lwd=2)
lines(Fmatrix_1_$bspl3.11, col = "#EEA2AD", lwd=2)
abline(h=0,col = "black")
legend("topright", 
       legend=c("bspl3.1", "bspl3.3","bspl3.7","bspl3.11"), 
       col=c( "red","green3","blue", "#EEA2AD"), 
       lwd=2, 
       cex=0.8, 
       title="Funciones")

star <- rep("C",16)
nc <- rep("NC",73)
starc <- c(star,nc)

barplot(table(starc))

prop <- prop.table(table(starc))
windows()
bp <- barplot(table(starc), 
              col="skyblue2", 
              main="Tipos de estrellas variables binarias eclipsantes",
              ylab="Número de estrellas", xlab ="Tipo de estrella")
grid(nx=NA, ny=NULL, col="gray", lty="dotted")
# Añadir texto de proporción encima de cada barra
text(bp, table(starc) -4, labels=paste0(round(prop*100,1),"% - ",table(starc)))

