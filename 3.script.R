#librerias
library(dplyr)
library(gtsummary)
library(stats)
library(factoextra)
library(MASS)
library(ggplot2)

#dataset
df <- read.csv("Dataset expresión genes.csv")
df_genes <- df %>% select(starts_with("AQ"))


#PCA
set.seed(1234)
pca.result <- prcomp(df_genes, center = TRUE, scale = TRUE )
View(pca.result$x)

#Graficos PCA ¿cuantos Componentes principales elegimos?
#Si miramos los eigenvalues:

eigenvalues <- get_eigenvalue(pca.result)
eigenvalues #parece que hay que coger hasta la dimension 5 > miramos grafico
fviz_eig(pca.result, addlabels = TRUE) #hasta la dimension 4 ¿?
fviz_pca_var(pca.result, col.var = "cos2", gradient.cols = c("blue","yellow", "red"), repel = TRUE)

#Clusterizacion
kmeans <- kmeans(t(df_genes), centers = 2)
summary(as.factor(kmeans$cluster))

fviz_pca_var(pca.result, col.var = kmeans$cluster, gradient.cols = c("blue","green", "red"), legend.title ="Cluster", repel = TRUE) 

#Nombre cluster verde > metabolismo energetico
#Nombre cluster azul > respuesta inmune e inflamatoria
#Nombre cluster rojo > metabolismo celular y respuesta inmune 

# comentario Clusters:
## el cluster rojo contiene genes tanto del metabolismo celular como de la respuesta inmune. Debemos mejorar la clusterización (tema 4 de Algoritmos)

fviz_contrib(pca.result, choice ="var", axes = 2, top = 20)

fviz_cluster(kmeans, df_genes)

fviz_pca_ind(pca.result, col.ind = "cos2", gradient.cols = c("blue","yellow", "red"), repel = TRUE)

#codificar una nueva variable
df$metastasisnosi <- as.factor(ifelse( df$extension == "metastasico", "metastasis", "no metastasis"))

df_pca <- as.data.frame(pca.result$x)

ggplot(df_pca, aes(x =PC1, y = PC4, color = df$metastasisnosi))+geom_point(size = 3)


#Tabla descriptiva en terciles

#Regresion logistica

