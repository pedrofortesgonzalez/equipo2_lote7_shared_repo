#######################################
# Cargar librerias #
#######################################

library(stats)       # Contiene funciones básicas para análisis multivariado, como PCA (prcomp), MDS (cmdscale) y escalado (scale).
library(ggplot2)     # Paquete para crear visualizaciones gráficas, útil para representar resultados en 2D o 3D tras la reducción de dimensionalidad.
library(Rtsne)       # Implementación del método t-SNE, ideal para explorar relaciones no lineales en los datos y visualizarlos en espacios de baja dimensión.
library(FNN)         # Proporciona herramientas para cálculos rápidos de vecinos más cercanos (KNN), necesarias para t-SNE, UMAP, Isomap, entre otros.
library(plotly)      # Permite crear gráficos interactivos en 2D y 3D, facilitando la exploración de los datos proyectados o clusters.
library(factoextra)  # Herramientas para calcular y visualizar resultados de clustering (k-means, clustering jerárquico) y análisis multivariado.
library(cluster)     # Incluye métodos de clustering como k-means, clustering jerárquico (agnes), y divisivo (DIANA).


#######################################
# cargar datos #
#######################################

#Base de datos tratada
df_genes <- read.csv("df.csv")

labels <- read.csv(file.choose(), header = TRUE) # Cargar archivo de classes (con etiquetas de cáncer)

# Extraer las etiquetas
class_labels <- labels$cancer

# Validar que las etiquetas contengan la columna 'cancer'
if (!"cancer" %in% colnames(labels)) {
  stop("El archivo de etiquetas no contiene una columna llamada 'cáncer'. Verifica los datos.")
}

# Extraer las etiquetas
class_labels <- labels$cancer

# Crear un nuevo dataframe a partir de la tercera columna de df
df_genes_sin_etiquetas <- df_genes[, 3:ncol(df_genes)]

# Validar que todas las columnas de df_genes_sin_etiquetas sean numéricas
if (!all(sapply(df_genes_sin_etiquetas, is.numeric))) {
  stop("El dataframe de genes contiene valores no numéricos. Revisa y corrige los datos.")
}

##Escalar datos

df_genes_sin_etiquetas_scale <- scale(df_genes_sin_etiquetas)  # Normalización z-score


#######################################
# Algoritmos de reducción de dimensionalidad #
#######################################

#######################################
# Analisis de componentes principales (PCA) #
#######################################

# Calcular componentes principales con la función prcomp
pca.results <- prcomp(df_genes_sin_etiquetas_scale, center = TRUE, scale = FALSE)

# Resultado de las componentes principales
pca.df <- data.frame(pca.results$x)

# Varianza (cuadrado de la desviación típica)
varianzas <- pca.results$sdev^2

# Total de la varianza de los datos
total.varianza <- sum(varianzas)

# Varianza explicada por cada componente principal
varianza.explicada <- varianzas / total.varianza

# Calcular la varianza acumulada
varianza.acumulada <- cumsum(varianza.explicada)

# Tomar el número de componentes principales que explican el 90% de la varianza
n.pc <- min(which(varianza.acumulada > 0.9))

# Imprimir resultados
print(paste("Número de componentes principales que explican el 90% de la varianza:", n.pc))

###Verificación de cuántas variables explican el porcentaje y la variabilidad de los genes

# Calcular la varianza acumulada
varianza.acumulada <- cumsum(varianza.explicada)

# Mostrar la proporción acumulada
print(varianza.acumulada)

#######################################
# Grafica (PCA) #
#######################################

# Etiquetas de los ejes del gráfico
x_label <- paste0(paste('PC1', round(varianza.explicada[1] * 100, 2)), '%')
y_label <- paste0(paste('PC2', round(varianza.explicada[2] * 100, 2)), '%')

# Representación gráfica de las primeras dos componentes principales respecto a los datos
ggplot(pca.df, aes(x=PC1, y=PC2, color=labels$cancer)) +
  geom_point(size=3) +
  scale_color_manual(values=c('red', 'blue', 'green', 'orange', 'purple')) +
  labs(title='PCA RNA-seq genes', x=x_label, y=y_label, color='Grupo') +
  theme_classic() +
  theme(panel.grid.major = element_line(color="black"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title = element_text(hjust = 0.5))


#######################################
# t-SNE #
#######################################

# Ingresar la semilla para hacer el algoritmo replicable
set.seed(1234)

# Reducción de dimensionalidad con t-SNE en 2D
tsne_2d <- Rtsne(X = df_genes_sin_etiquetas_scale, perplexity = 15, dims = 2, check_duplicates = FALSE, theta = 0.2, pca = TRUE)
tsne_result_2d <- data.frame(tsne_2d$Y)

# Reducción de dimensionalidad con t-SNE en 3D
tsne_3d <- Rtsne(X = df_genes_sin_etiquetas_scale, perplexity = 15, dims = 3, check_duplicates = FALSE, theta = 0.2, pca = TRUE)
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
##de proximidad de los puntos en el espacio reducido respecto al espacio original. Un valor
##cercano a 1 sugiere una conservación alta, mientras que valores más bajos indican una menor 
##preservación de la estructura en la reducción dimensional.

# Calcular la tasa de conservación en 2D
rate_2d <- conservation_rate(original_data = df_genes_sin_etiquetas_scale, reduced_data = tsne_result_2d, k = k)
print(paste("La tasa de conservación de los", k, "vecinos más cercanos en 2D es:", rate_2d))

# Calcular la tasa de conservación en 3D
rate_3d <- conservation_rate(original_data = df_genes_sin_etiquetas_scale, reduced_data = tsne_result_3d, k = k)
print(paste("La tasa de conservación de los", k, "vecinos más cercanos en 3D es:", rate_3d))


# Graficar los resultados de t-SNE en 2D
ggplot(tsne_result_2d, aes(x = X1, y = X2, color = labels$cancer)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Método t-SNE (2D) RNA-seq genes", x = "Dim 1", y = "Dim 2", color = "Grupo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "black"), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "grey95"), 
        plot.title = element_text(hjust = 0.5))


#######################################
# Clustering Jerárquico (Método DIANA) #
#######################################

# Implementación del clustering divisivo

diana_euclidean <- diana(df_genes_sin_etiquetas_scale, metric = "euclidean", stand = FALSE)
diana_manhattan <- diana(df_genes_sin_etiquetas_scale, metric = "manhattan", stand = FALSE)

colors <- rainbow(5)
clust_diana_euclidean <- fviz_dend(diana_euclidean, 
                                   cex = 0.5, 
                                   k = 5,
                                   palette = colors, 
                                   main = 'Euclidean',
                                   xlab = "Índice de Observaciones",
                                   ylab = "Distancia") + 
  theme_classic()

colors <- rainbow(5)
clust_diana_manhattan <- fviz_dend(diana_manhattan, 
                                   cex = 0.5, 
                                   k = 5,
                                   palette = colors, 
                                   main = 'Manhattan',
                                   xlab = "Índice de Observaciones",
                                   ylab = "Distancia") + 
  theme_classic()

grid.arrange(clust_diana_euclidean, clust_diana_manhattan, nrow = 2)

####Al realizar el análisis de clustering jerárquico utilizando **DIANA**, los dendrogramas 
##generados no presentan una visualización clara de los grupos más pequeños. Esto podría 
##deberse a la alta dimensionalidad de los datos, ya que cuando se manejan numerosos genes 
##o muestras, el dendrograma se vuelve denso y complejo, dificultando su interpretación.
##Las diferencias entre los grupos más pequeños no son fácilmente distinguibles, lo que
##podría ser resultado de que las distancias entre los puntos de datos no están suficientemente 
##diferenciadas, lo que impide que los grupos más pequeños se destaquen adecuadamente. Por tanto
##no es una buena técnica a usar con este tipo de muestras.

#######################################
# clustering no jerárquico: K.means #
#######################################

##Primeramente, se determina el número óptimo de clusters mediante el método del **codo** 
##(WSS: within-cluster sum of squares), utilizando la función `fviz_nbclust`. Este método 
##ayuda a identificar el número adecuado de clusters observando la caída en la variabilidad 
##dentro de los grupos a medida que se aumentan los clusters. Luego, se realiza el clustering 
## **K-means** con 4 clusters definidos manualmente (`centers = 4`) por la regla del codo


# Número óptimo de clusters
fviz_nbclust(df_genes_sin_etiquetas_scale, kmeans, method = "wss") +
  ggtitle("Número óptimo de clusters", subtitle = "") +
  theme_classic()

# Realizar el clustering k-means
kmeans.result <- kmeans(df_genes_sin_etiquetas_scale, centers = 4, iter.max = 100, nstart = 25)

# Visualizar el clustering y las etiquetas
fviz_cluster(kmeans.result, df_genes_sin_etiquetas_scale, 
             xlab = '', ylab = '', 
             geom = "point", 
             labelsize = 4, 
             ggtitle("Cluster K-means de RNA-seq genes") +
               theme_minimal() +
               theme(plot.title = element_text(hjust = 0.5, margin = margin(b = -10))))

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
kmeans.result <- kmeans(df_genes_sin_etiquetas_scale, centers = 2, iter.max = 100, nstart = 25)

# Visualizar el clustering y las etiquetas
fviz_cluster(kmeans.result, df_genes_sin_etiquetas_scale, 
             xlab = '', ylab = '', 
             geom = "point", 
             labelsize = 4, 
             ggtitle("Cluster K-means de RNA-seq genes") +
               theme_minimal() +
               theme(plot.title = element_text(hjust = 0.5, margin = margin(b = -10))))

