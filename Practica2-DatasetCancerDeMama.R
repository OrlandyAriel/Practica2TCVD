
#Carga de librerías
install.packages("https://cran.r-project.org/src/contrib/nortest_1.0-4.tar.gz", repos=NULL)
library(caret)
library(rattle)
library(lda)
library(readr)
library(C50)
library(ggplot2)
library(gmodels)
library(corrplot)
library(psych)
library(GGally)
library(class)
library(arules)
library(grid)
library(gridExtra)
library(normtest)
library(nortest)
#Carga datos
datosCancer<- read_csv("/datos/breast-cancer-wisconsin.data.csv")
# Se cambia el nómbre a los atributos para facilitar su lectura en los distintos gráficos
names(datosCancer)<-c ("v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10")
summary(datosCancer)
#class = 2 venigno, 4 maligno
# Cargo los datos en memoría para poder trabajar más rapido con ellos 
attach(datosCancer)
#información del dataset
names(datosCancer)
dim(datosCancer)
summary(datosCancer)
sapply(datosCancer, function(x) class(x))
table(datosCancer$v10)
head(datosCancer)

# Limpieza de los datos
# ver los nulos, vacios y cero (en caso que no sea normal)
# Estadísticas de valores nulos, vacios y ceros
colSums(is.na(datosCancer))
colSums(datosCancer== "")
colSums(datosCancer== 0) ## en este caso no se puede hacer nada puesto que los ceros aquí son validos al ser variables númericas
# elimino las variables sin importancia, en este caso es el id de la fila
datosCancer$v0 <- NULL
# elimino los datos repetido
# podría datosCancer <- unique(datosCancer)

datos <- datosCancer
#Convierto la variable 6, puesto que se que es numerica y R la interpreta como factor, a una variable numerica.
datos$v6 <- as.numeric(as.factor(datos$v6))
# normalizo los datos
normalizacion <- function(x){
  return ((x-min(x))/(max(x)-min(x)))
}
  # en el proceso se pierden formato por tanto se realizan unas conversiones para poder mantener el formato de las celdas a numeric
datos <- as.data.frame(lapply(datos,normalizacion))
datos<-format(datos, digits=3, nsmall =2)

datos <- as.data.frame(sapply(datos,as.numeric))


## normalidad y homogeniedad

benigno <- subset(datos, v10 == 0)
maligno <- subset(datos, v10 == 1) 
hist(benigno$v1)
hist(datos$v2, breaks = 20, main = "", xlab = "Cell size uniformity ", border = "darkred")
hist(datos$v3, breaks = 20, main = "", xlab = "Cell shape uniformity",border = "blue")
# normalidad Geary
normlGeary <- function(x){
  for (i in 1:ncol(x)) {
    if(0.05 < geary.norm.test(x[,i],nrepl=2000)$p.value){
      cat(names(x[i]))
      cat(" Es normal\n")
    }
    
  }
}
datosPrecioPuntosNorm <- datosPrecioPuntosNorm[rnum,]

normlGeary(benigno[-10])

## seleccionadno los datos en benigno y maligno 
benigno2 <- benigno[benigno$v2 < 0.2 & benigno$v3 <0.2,]
maligno2 <- maligno[maligno$v2 > 0.2 & maligno$v3 > 0.2,]

porctajeBenigno<- (length(benigno2$v1) /length(benigno$v1)) * 100
porctajeMaligno<- (length(maligno2$v1) /length(maligno$v1)) * 100

# matriz decorrelación

cor(x = datos, method = "pearson")
# puedo ver gráficamente y al menos visualmente no parecen tener correlación lineal 
multi.hist(x = datos, dcol = c("blue", "red"), dlty = c("dotted", "solid"),  main = "")

#ver la significancia de la correlación entre las variables v2 y v10 que parecen ser las que mas tienen entre ellas
cor.test(x = datos$v2, y = datos$v10, alternative = "two.sided", conf.level = 0.95, method = "pearson")

#ver la significancia de la correlación entre las variables v3 y v10 que  tambien es alta
cor.test(x = datos$v3, y = datos$v10, alternative = "two.sided", conf.level = 0.95, method = "pearson")

#ver la significancia de la correlación entre las variables v7 y v10 que  tambien es alta, aunque la más baja de entre de v2,v3 y v7
cor.test(x = beni$v7, y = datos$v10, alternative = "two.sided", conf.level = 0.95, method = "pearson")

#veamos una mátriz de correlación, donde podemos ver claramente los que tienen la correlación más fuerte
corrplot(corr = cor(x = benigno[c(3:8)], method = "pearson"), method = "number")
#mátriz de correlación con grafica en la diagonal
pairs.panels(x = datos, ellipses = FALSE, lm = TRUE, method = "pearson")
#mátriz de correlación con grafica
ggpairs(datos, lower = list(continuous = "smooth"), 
        diag = list(continuous = "bar"), axisLabels = "none")
# Test de shapiro
shapiro.test(datos$v2)
shapiro.test(datos$v3)

# regresión lineal simple
lm(formula = v2 ~ v10, data = datos) 
#Cálculo del modelo de regresión lineal simple 
cor.test(x = datos$v2, y = datos$v10, method = "spearman")
cor.test(x = datos$v2, y = datos$v10, method = "pearson")
# en ambos casos los coeficientes de correlación son significativos.
# correlación
modelo_lineal <- lm(v2 ~ v10, datos)
# lm() devuelve el valor de la variable y para x=0 (intersección) junto 
# con la pendiente de la recta.
# Para ver la información del modelo se requiere summary().
summary(modelo_lineal)

#intervalos de confianza 
confint(modelo_lineal)
## Resultados matriz filtro por v2 y v3
porctajeBenigno<- (length(benigno2$v1) /length(benigno$v1)) * 100
porctajeMaligno<- (length(maligno2$v1) /length(maligno$v1)) * 100

matrizFiltroV2V3 <- matrix(nc = 2, nr = 2)

matrizFiltroV2V3[1, 1] <- "Porcentaje benigno con v2 y v3 < 0.2 "
matrizFiltroV2V3[2, 1] <- "Porcentaje Maligno con v2 y v3 > 0.2"
matrizFiltroV2V3[1, 2] <- porctajeBenigno
matrizFiltroV2V3[2, 2] <- porctajeMaligno
print(matrizFiltroV2V3)
######################------KNN---------#########################
summary(datos)

datosParaKNN <-datos
datosParaKNN$v10[datosParaKNN$v10 == 0] <- "Benigno"
datosParaKNN$v10[datosParaKNN$v10 == 1] <- "Maligno"
str(datosParaKNN)
#ver cantidad de casos de cada tipo de cancer
table(datosParaKNN$v10)

# ver en porcentaje
round(prop.table(table(datos$v10))* 100, digits = 1 )

# preparación datos para entrenamiento y pruebas
# se hará el siguiente reparto 70% entrenamientos 30% pruebas
particion <- sample(2, nrow(datosParaKNN), replace =  TRUE, prob = c(0.7, 0.3))

trainDatos <- datosParaKNN[particion == 1,]
testDatos <- datosParaKNN[particion == 2,]
# guardo las los datos separados sin la clase
trainDatosNoTag <- trainDatos[-10]
testDatosNoTag <- testDatos[-10]
# compruebo que los datos se han guadado como deben
dim(trainDatos)
dim(trainDatosNoTag)
dim(testDatos)
dim(testDatosNoTag)

# guardo las clases separadas 
trainDatosTag <- trainDatos$v10
testDatosTag <- testDatos$v10

# se pasa a realizar el entrenamiento

testPredicion <- knn(train = trainDatosNoTag, test = testDatosNoTag, cl = trainDatosTag, k = 3, prob = TRUE)
summary(testPredicion)

# evaluación del modelo obtenido, matriz de confusión
matrizConfusion <-CrossTable(x = testDatosTag, y = testPredicion, prop.chisq =  FALSE)

