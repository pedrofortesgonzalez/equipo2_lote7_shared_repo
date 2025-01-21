# 1-PRE-PROCESADO --------------------------------------------------------------
## 1.1. WEnV -------------------------------------------------------------------
rm(list=ls()) # resetear WEnv
path <- "/Users/pedrofortesgonzalez/Desktop/MASTER_UNIR/11_ALGORITMOS_IA/actividades/Act3_grupal/Material complementario"
setwd(path)

# Librerías y random seed
set.seed(1999)
library(tidyverse) # contiene ggplot2, dplyr, tidyr, readr, etc.
library(caret)       # funciones para ML
library(randomForest)# " para random forests
library(klaR)        # " para visualizar regularización
library(glmnet)      # funciones para calcular regularización
library(pROC)        # representación de curvas ROC
library(stats)       # Contiene funciones básicas para análisis multivariado, como PCA (prcomp), MDS (cmdscale) y escalado (scale).
library(Rtsne)       # Implementación del método t-SNE, ideal para explorar relaciones no lineales en los datos y visualizarlos en espacios de baja dimensión.
library(FNN)         # Proporciona herramientas para cálculos rápidos de vecinos más cercanos (KNN), necesarias para t-SNE, UMAP, Isomap, entre otros.
library(plotly)      # Permite crear gráficos interactivos en 2D y 3D, facilitando la exploración de los datos proyectados o clusters.
library(factoextra)  # Herramientas para calcular y visualizar resultados de clustering (k-means, clustering jerárquico) y análisis multivariado.
library(cluster)     # Incluye métodos de clustering como k-means, clustering jerárquico (agnes), y divisivo (DIANA).

## Estructurar df---------------------------------------------------------------
clases <- read.csv("./1_data/classes.csv", header = FALSE, sep = ";", col.names = c("Muestra", "Clase"), row.names = 1)
columnas <- read_lines("1_data/column_names.txt")
gen_exp <-read.csv("1_data/gene_expression.csv", header = FALSE, sep = ";", col.names = (columnas))
gen_exp <- scale(gen_exp, center = TRUE, scale=TRUE)  # Normalizamos datos de expresión génica
df <- cbind(clases, gen_exp)
df$Clase <- as.factor(df$Clase) # convertir variable clase a factor

## Procesar NAs-----------------------------------------------------------------
any(sum(is.na(df))) # Comprobamos si hay valores NA
na_columns <- colnames(gen_exp)[colSums(is.na(gen_exp)) > 0] # Tras ver que sí, las guardo en un vector 
na_columns # y las visualizo (son solo 3 de 501, podemos eliminarlas)
df <- df[, !colnames(df) %in% na_columns] # Y las eliminamos ya que no aportan información y tenemos muchas otras aún por analizar
gen_exp <- gen_exp[, !colnames(gen_exp) %in% na_columns] # Las eliminamos también en gen_exp (si no, en los algoritmos NSL R nos dará errores)
any(sum(is.na(df))) #Comprobamos que ya no hay columnas con todo 0
write.csv(df,"1_data/df1_completo.csv") # Saco otra copia del dataframe para backup

## Dimensiones del df-----------------------------------------------------------
dim(df)

# Tenemos un df con alta dimensionalidad: n_cols >> n_rows / 10




# 2.AP. NO SUPERVISADO----------------------------------------------------------
# Non Supervised Learning
# Aquí utilizaremos los genes sin etiquetas (que están en el df: gen_exp)

## 2.1.REDUCIR DIMENSIONES------------------------------------------------------
### PCA-------------------------------------------------------------------------
#### Cálculos PCA --------------------------------------------------------------
pca.results <- prcomp(gen_exp, center = TRUE, scale = FALSE) # Calcular componentes principales
pca.df <- data.frame(pca.results$x) # Resultado de las componentes principales
varianzas <- pca.results$sdev^2 # Varianza (cuadrado de la desviación típica)
total.varianza <- sum(varianzas) # Total de la varianza de los datos
varianza.explicada <- varianzas / total.varianza # Varianza explicada por cada componente principal
varianza.acumulada <- cumsum(varianza.explicada) # Calcular la varianza acumulada
n.pc <- min(which(varianza.acumulada > 0.9)) # Tomar nº de PCs que explican el 90% de la varianza

# Imprimir resultados
print(paste("Número de componentes principales que explican el 90% de la varianza:", n.pc))
# Necesitamos 155 PCs para explicar el 90% de la varianza... PCA no parece que funcione demasiado bien

###Verificación de cuántas variables explican el porcentaje y la variabilidad de los genes
varianza.acumulada <- cumsum(varianza.explicada) # Calcular la varianza acumulada
print(varianza.acumulada) # Mostrar la proporción acumulada


#### Gráfica PCA ---------------------------------------------------------------
# Etiquetas de los ejes del gráfico
x_label <- paste0(paste('PC1', round(varianza.explicada[1] * 100, 2)), '%')
y_label <- paste0(paste('PC2', round(varianza.explicada[2] * 100, 2)), '%')

# Representación gráfica de las primeras dos componentes principales respecto a los datos
ggplot(pca.df, aes(x=PC1, y=PC2, color=clases$Clase)) +
  geom_point(size=3) +
  scale_color_manual(values=c('red', 'blue', 'green', 'orange', 'purple')) +
  labs(title='PCA RNA-seq genes', x=x_label, y=y_label, color='Grupo') +
  theme_classic() +
  theme(panel.grid.major = element_line(color="black"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title = element_text(hjust = 0.5))



### t-SNE-----------------------------------------------------------------------
#### Cálculos t-SNE-------------------------------------------------------------
# Reducción de dimensionalidad con t-SNE en 2D
tsne_2d <- Rtsne(X = gen_exp, perplexity = 15, dims = 2, check_duplicates = FALSE, theta = 0.2, pca = TRUE)
tsne_result_2d <- data.frame(tsne_2d$Y)
# Reducción de dimensionalidad con t-SNE en 3D
tsne_3d <- Rtsne(X = gen_exp, perplexity = 15, dims = 3, check_duplicates = FALSE, theta = 0.2, pca = TRUE)
tsne_result_3d <- data.frame(tsne_3d$Y)
# Número de vecinos más cercanos a considerar
k <- 10

# Función para calcular la tasa de conservación de los k-vecinos más cercanos
conservation_rate <- function(original_data, reduced_data, k) {
  original_nn <- get.knnx(data = original_data, query = original_data, k = k)
  reduced_nn <- get.knnx(data = reduced_data, query = reduced_data, k = k)
  overlap_count <- sapply(1:nrow(original_data), function(i) {
    length(intersect(original_nn$nn.index[i, ], reduced_nn$nn.index[i, ]))
    })
  mean(overlap_count) / k
}

###Cálculo de nivel de conservacion de los de los plots

##Calculo la tasa de conservación para dos configuraciones: 2D y 3D, utilizando un 
##valor de **k** para definir cuántos vecinos considerar. Imprimo el resultado de la 
##tasa de conservación para cada caso, lo que me indica cuán bien se conservan las relaciones 
##de proximidad de los puntos en el espacio reducido respecto al espacio original. 
### Valor cercano a 1 sugiere una conservación alta
### Valores más bajos indican una menor preservación de la estructura en la reducción dimensional.

# Calcular la tasa de conservación en 2D
rate_2d <- conservation_rate(original_data = gen_exp, reduced_data = tsne_result_2d, k = k)
print(paste("La tasa de conservación de los", k, "vecinos más cercanos en 2D es:", rate_2d))
# Calcular la tasa de conservación en 3D
rate_3d <- conservation_rate(original_data = gen_exp, reduced_data = tsne_result_3d, k = k)
print(paste("La tasa de conservación de los", k, "vecinos más cercanos en 3D es:", rate_3d))


#### Gráfica t-SNE -------------------------------------------------------------
# Graficar los resultados de t-SNE en 2D
ggplot(tsne_result_2d, aes(x = X1, y = X2, color = clases$Clase)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Método t-SNE (2D) RNA-seq genes", x = "Dim 1", y = "Dim 2", color = "Grupo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "black"), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "grey95"), 
        plot.title = element_text(hjust = 0.5))




## 2.2. CLUSTERING--------------------------------------------------------------
### Jerárquico: DIANA ----------------------------------------------------------
##### Implementación -----------------------------------------------------------
diana_euclidean <- diana(gen_exp, metric = "euclidean", stand = FALSE)
diana_manhattan <- diana(gen_exp, metric = "manhattan", stand = FALSE)
colors <- rainbow(5) # definimos paleta de colores
clust_diana_euclidean <- fviz_dend(diana_euclidean, cex = 0.5, k = 5,
                                   palette = colors, main = 'Euclidean',
                                   xlab = "Índice de Observaciones",
                                   ylab = "Distancia") + theme_classic()

clust_diana_manhattan <- fviz_dend(diana_manhattan, cex = 0.5, k = 5,
                                   palette = colors, main = 'Manhattan',
                                   xlab = "Índice de Observaciones",
                                   ylab = "Distancia") + theme_classic()

grid.arrange(clust_diana_euclidean, clust_diana_manhattan, nrow = 2)

##### Análisis del resultado----------------------------------------------------
####Al realizar el análisis de clustering jerárquico utilizando **DIANA**, los dendrogramas 
##generados no presentan una visualización clara de los grupos más pequeños. Esto podría 
##deberse a la alta dimensionalidad de los datos, ya que cuando se manejan numerosos genes 
##o muestras, el dendrograma se vuelve denso y complejo, dificultando su interpretación.
##Las diferencias entre los grupos más pequeños no son fácilmente distinguibles, lo que
##podría ser resultado de que las distancias entre los puntos de datos no están suficientemente 
##diferenciadas, lo que impide que los grupos más pequeños se destaquen adecuadamente. Por tanto
##no es una buena técnica a usar con este tipo de muestras.



### No jerárquico: K.means------------------------------------------------------
##### Planteamiento-------------------------------------------------------------
##Primeramente, se determina el número óptimo de clusters mediante el método del **codo** 
##(WSS: within-cluster sum of squares), utilizando la función `fviz_nbclust`. Este método 
##ayuda a identificar el número adecuado de clusters observando la caída en la variabilidad 
##dentro de los grupos a medida que se aumentan los clusters. Luego, se realiza el clustering 
## **K-means** con 4 clusters definidos manualmente (`centers = 4`) por la regla del codo

##### Cálculos y plot ----------------------------------------------------------
# Número óptimo de clusters
fviz_nbclust(gen_exp, kmeans, method = "wss") +
  ggtitle("Número óptimo de clusters", subtitle = "") +  theme_classic()

# Realizar el clustering k-means
kmeans.result <- kmeans(gen_exp, centers = 4, iter.max = 100, nstart = 25)

# Visualizar el clustering y las etiquetas
fviz_cluster(kmeans.result, gen_exp, xlab = '', ylab = '', geom = "point", labelsize = 4, 
             ggtitle("Cluster K-means de RNA-seq genes") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5)))

##El análisis de K-means con 4 clusters muestra que tres de ellos se superponen, mientras 
##que uno se separa claramente. Aunque la regla del codo indica que 4 clusters es el número 
##óptimo, la superposición de los clusters sugiere que los datos no se agrupan de manera 
##clara en 4 grupos, lo que podría indicar que este número no es el más adecuado. Al probar 
##con solo 2 clusters, la separación es más evidente, con la mayoría de los datos concentrados 
##en el primer cluster y otro grupo distinto en el segundo. Sin embargo, al aumentar el número 
##de clusters (3, 4, etc.), la superposición aumenta dentro del primer cluster, lo que sugiere
##que la mayoría de la expresión génica se concentra en un grupo principal, con valores similares 
##o cercanos entre sí.

# Realizar el clustering k-means
kmeans.result <- kmeans(gen_exp, centers = 2, iter.max = 100, nstart = 25)

# Visualizar el clustering y las etiquetas
fviz_cluster(kmeans.result, gen_exp, xlab = '', ylab = '', geom = "point", labelsize = 4, 
             ggtitle("Cluster K-means de RNA-seq genes") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5)))



# 3.REGULARIZACIÓN: LASSO-------------------------------------------------------
## Tenemos 801 filas y 498 columnas = alta dimensionalidad.
## Debemos escoger un método para elegir columnas ya que, si no, tendremos mucha info (y por tanto ruido) en nuestro análisis
## Probamos a aplicar LASSO para quedarnos con el subconjunto óptimo de variables

## Implementación --------------------------------------------------------------
x <- as.matrix(df[,3:498])
y <- as.factor(df$Clase)
lasso.result <- cv.glmnet(x,y, family = "multinomial", alpha = 1)
lasso.coef <- coef(lasso.result, s = "lambda.min")

## Convertir los coeficientes a matrices, y estas a dataframes
coef_AGH_df <- as.data.frame(as.matrix(lasso.coef[["AGH"]]))
coef_CFB_df <- as.data.frame(as.matrix(lasso.coef[["CFB"]]))
coef_CGC_df <- as.data.frame(as.matrix(lasso.coef[["CGC"]]))
coef_CHC_df <- as.data.frame(as.matrix(lasso.coef[["CHC"]]))
coef_HPB_df <- as.data.frame(as.matrix(lasso.coef[["HPB"]]))
## Asegurarse de que los nombres de las variables estén en las columnas
colnames(coef_AGH_df) <- "AGH"
colnames(coef_CFB_df) <- "CFB"
colnames(coef_CGC_df) <- "CGC"
colnames(coef_CHC_df) <- "CHC"
colnames(coef_HPB_df) <- "HPB"

## Combinar los dataframes
df_final <- cbind(coef_AGH_df, coef_CFB_df, coef_CGC_df, coef_CHC_df, coef_HPB_df)
df_final <- df_final[2:497, ] #quitar el intercept
#Quitamos variables donde para todos los niveles el valor sea 0 
df_final_2 <- df_final[rowSums(df_final) != 0, ]  
#recuperamos el nombre de aquellas variables que hemos seleccionado con LASSO
variables_elegidas <- rownames(df_final_2)
#Del df original, solo nos quedamos con la variable categorica y con las variables numericas seleccionadas
df_analisis <- df %>% dplyr :: select(Clase, all_of(variables_elegidas))
dim(df_analisis)
write.csv(df,"1_data/df2_lasso.csv") # Saco una copia del dataframe para backup

## Resultado -------------------------------------------------------------------
## Nos hemos quedado con 801 filas y 104 columnas. Ahora podremos hacer 
## el análisis en un dataframe sin alta dimensionalidad





# 4.AP. SUPERVISADO-------------------------------------------------------------
## 4.0.TRAIN / TEST-------------------------------------------------------------
#Separar en conjunto de entrenamiento y conjunto de prueba
index_train <- createDataPartition(df_analisis$Clase, p = 0.8, list = FALSE)
df_train <- df_analisis[index_train, ]
df_test <- df_analisis[-index_train, ]
# Para los metodos de analisis discriminante, hay que aplicar una formula:
formula  <- as.formula(paste("Clase ~ ",paste(variables_elegidas, collapse = "+")))
# Guardamos clases reales en una variable aparte para construir matrices de confusión
clases.reales <- df_test$Clase



## 4.1.RDA----------------------------------------------------------------------
rda.result <- rda(formula, data = df_train) # entrenamos RDA
rda.prediction <- predict(rda.result, newdata = df_test) # predecimos en test dataset



## 4.2.SVM----------------------------------------------------------------------
#Como el RDA ha ido bien > parece que los datos tienen relaciones lineales > usamos la variable lineal 
svm.result <- train(Clase ~., data = df_train, method = "svmLinear",
                    trControl = trainControl(method = "cv", number = 5),
                    prob.model = TRUE, tuneGrid = expand.grid(C =seq(1,20,0.5)))

svm.prediction <- predict(svm.result, newdata = df_test)
plot(svm.result) # C con el Accuracy más alto = 7.5



## 4.3.RANDOM FOREST------------------------------------------------------------
#Es mas costoso computacionalmente, pero puede merecer la pena usar un metodo mas robusto para compararlo con los otros dos mas simples
randomforest.result <- train(Clase ~., data = df_train, method = "rf",
                             trControl = trainControl(method = "cv", number = 5),
                             tuneLength = 30)

plot(randomforest.result)
randomforest.prediction <- predict(randomforest.result, newdata = df_test)

varImp(randomforest.result) #variables mas importantes en el modelo
varImpPlot(randomforest.result$finalModel)




## 4.4 EVALUAR AP. SUPERVISADO--------------------------------------------------
### Matrices de confusión
### Las obtengo para construir la matriz conjunta
#### RDA
probabilidades.rda <- rda.prediction$posterior
clases.rda <- rda.prediction$class
matriz.rda <- confusionMatrix(clases.rda, clases.reales, mode = "everything") # matriz de confusión RDA
t_m_rda <- t_m_rf <- t(matriz.rda$table) # matriz traspuesta (para mejor visualización al verlas todas juntas)
colnames(t_m_rda) <- paste("RDA", colnames(t_m_rda), sep="_") #añado título del método a cada col. de la conf.matrix

#### SVM LINEAL
probabilidades.svm <- predict(svm.result, newdata = df_test, type = "prob")
matriz.svm <- confusionMatrix(svm.prediction, clases.reales, mode = "everything")
t_m_svm <- t(matriz.svm$table) # matriz traspuesta
colnames(t_m_svm) <- paste("SVM", colnames(t_m_svm), sep="_") # añado título del método a cada col.

#### RF
probabilidades.rf <- predict(randomforest.result, newdata = df_test, type = "prob")
matriz.rf <- confusionMatrix(randomforest.prediction, clases.reales, mode = "everything")
t_m_rf <- t(matriz.rf$table) # matriz traspuesta
colnames(t_m_rf) <- paste("RF", colnames(t_m_rf), sep="_") # añado título del método a cada col.

#### matriz de confusión conjunta: construir
cm <- cbind(t_m_rda, t_m_svm, t_m_rf) # Combino las matrices traspuestas
rownames(cm) <- paste("Ref", rownames(matriz.rf$table), sep="_") # renombro filas con prefijo "Ref" (referencia)


#### Matriz de confusión conjunta-----------------------------------------------
cm # Visualizar la matriz conjunta



#### Obtener performance metrics de las matrices de confusión --> construir tabla resumen (performance metrics + AUCs)
matrices = list(matriz.rda$byClass, matriz.svm$byClass, matriz.rf$byClass) # lista de matrices
counts = c(146,300,141,136,78) # vector de conteos / clase (sacado de summary(df$Clase))
pmetrix <- c("Sensitivity", "Specificity", "Precision", "F1") # cols que nos interesan
model_names <- c("RDA", "SVM", "RF")

# Como las clases están desbalanceadas y tengo 5 métricas por modelo (una por cada clase)
# Haré una media ponderada de cada métrica para que dicha media sea representativa de la métrica en los conteos de cada clase
# Primero definimos la función para calcular la media ponderada para cada clase
weighted_mean_func <- function(metrics, counts) {
  sum(metrics * counts) / sum(counts)
}
# Ahora filtramos las columnas y calculamos su media ponderada
analisis <- sapply(matrices, function(matrix) {
  filtered_matrix <- matrix[, pmetrix] # Para cada matriz, dejar solo cols en pmetrix
  sapply(1:ncol(filtered_matrix), function(i) { # una vez filtradas, calcular media ponderada para cada col de la matriz
    weighted_mean_func(filtered_matrix[, i], counts)
  })})
analisis



###### Curvas ROC---------------------------------------------------------------
table(df$Clase) #Este ejemplo tratamos un problema multiclase (>2 clases)

#No se pueden usar curvas ROC ni PR como tal > solo cuando hay dos clases y aquí tenemos 5
#Pero si que podemos usar multiclass.roc: calcula las curvas ROC para cada clase y luego genera un promedio.
#No se puede graficar, pero se puede ver el area bajo la curva y se interpreta igual que las curvas ROC y PR
# Calculamos las multiclass.roc para cada método
curva_rda <- multiclass.roc(df_test$Clase, as.matrix(probabilidades.rda))
curva_svm <- multiclass.roc(df_test$Clase, as.matrix(probabilidades.svm))
curva_rf <- multiclass.roc(df_test$Clase, as.matrix(probabilidades.rf))

# Y guardamos resultados en el vector aucs
aucs <- c(curva_rda$auc,curva_svm$auc, curva_rf$auc)

## Perf.metrics: unificar todas en dataframe
# Convertir analisis a dframe para representar todos los datos en él
analisis_df <- as.data.frame(analisis) 
analisis_df <- rbind(analisis_df, aucs)
colnames(analisis_df) <- model_names
rownames(analisis_df) <- c(pmetrix, "AUC")



### RESUMEN: P.Metrics------------------------------------------------------
analisis_df


# 5.PREGUNTAS SOBRE LA ACTIVIDAD -----------------------------------------------
## 1. Procesamiento de los datos (0,5 puntos):
### ¿Qué método habéis escogido para llevar a cabo la imputación de los datos? Razonad vuestra respuesta. (0,3 puntos).
### ¿Habéis llevado a cabo algún otro tipo de procesamiento? Razonad vuestra respuesta. (0,2 puntos).


## 2.Métodos no supervisados (1 punto):
### ¿Cuál es el motivo por el cual habéis seleccionado estas técnicas de reducción de dimensionalidad? (0,3 puntos).
### ¿Cuál es el motivo por el cual habéis seleccionado estas técnicas de clusterización? (0,3 puntos).
### En ambos casos, ¿qué aspectos positivos y negativos tienen cada una? (0,2 puntos).
### En el caso de la clusterización, ¿podéis afirmar con certeza que los clústeres generados son los mejores posibles? Razonad vuestra respuesta. (0,2 puntos).


## 3.Métodos supervisados (1,75 puntos):
### ¿Cuál es el motivo por el cual habéis seleccionado ambas técnicas de aprendizaje supervisado? ¿Cuál ha dado mejores resultados a la hora de clasificar las muestras? Razonad vuestra respuesta (1 punto).
### ¿Habéis considerado oportuno implementar algún método de reducción de dimensionalidad para procesar los datos antes de implementarlos en dichas técnicas? ¿Por qué? (0,5 puntos).
### ¿Qué aspectos positivos y negativos tienen cada una de las técnicas que habéis escogido? (0,25 puntos).

## 4.De estas cuatro opciones, ¿qué tipo de arquitectura de deep learning sería la más adecuada para procesar datos de expresión génica? Razonad vuestra respuesta (0,25 puntos).
### a) Red de perceptrones (multiperceptron layers).
### b) Redes convolucionales.
### c) Redes recurrentes.
### d) Redes de grafos.

### La respuesta correcta es la opción A); la red de perceptrones es la arquitectura deep learning más utilizada en el ámbito y es más adecuada a la estructura y naturaleza de los datos. Estos están separados por tabulaciones y  componen una matriz de datos de expresión génica (vectores) y metadatos (que pueden codificarse para que el algoritmo lo comprenda). 
### El resto de arquitecturas no son adecuadas:
### Redes convolucionales: Destaca por el procesado de estructuras espaciales (imágenes o datos 2D/3D), los datos de expresión génica carecen de esta estructura.
### Redes recurrentes: Destaca por el procesado de datos secuenciales o temporales, puesto que aprenden en cada secuencia de procesado. En los datos de expresión génica generales no se cuenta con esta relación (a menos que sean estudios longitudinales).
### Redes de grafos: Destaca por el procesado de datos estructurados entre nodos, en el ámbito de la expresión génica podría visualizarse las relaciones entre genes, siendo estos los nodos y las relaciones las conexiones. Sin embargo, sólamente con datos de expresión génica no puede emplearse esta arquitectura.
