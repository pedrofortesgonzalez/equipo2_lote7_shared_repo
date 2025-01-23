#librerias
library(dplyr)
library(gtsummary)
library(stats)
library(factoextra)
library(MASS)
library(ggplot2)
library(glmnet)

#dataset
df <- read.csv("1_data/Dataset expresión genes.csv")
df_genes <- df %>% dplyr::select(starts_with("AQ"))
df_nogenes <- df %>% dplyr::select(-starts_with("AQ"))


#PCA
set.seed(1234)
pca.result <- prcomp(df_genes, center = TRUE, scale = TRUE )
View(pca.result$x)

#Graficos PCA ¿cuantos Componentes principales elegimos?
#Si miramos los eigenvalues:

eigenvalues <- get_eigenvalue(pca.result)
eigenvalues #parece que hay que coger hasta la dimension 5 > miramos grafico
fviz_eig(pca.result, addlabels = TRUE) #hasta la dimension 4 ¿?
fviz_pca_var(pca.result, col.var = "cos2", gradient.cols = c("blue","yellow", "red"), repel = TRUE) #variables segun calidad cos2

#Clusterizacion
kmeans <- kmeans(t(df_genes), centers = 2)
summary(as.factor(kmeans$cluster))

fviz_pca_var(pca.result, col.var = kmeans$cluster, gradient.cols = c("blue","green", "red"), legend.title ="Cluster", repel = TRUE) 

#Nombre cluster verde > metabolismo energetico
#Nombre cluster azul > respuesta inmune e inflamatoria
#Nombre cluster rojo > metabolismo celular y respuesta inmune 

# comentario Clusters:
## el cluster rojo contiene genes tanto del metabolismo celular como de la respuesta inmune. Debemos mejorar la clusterización (tema 4 de Algoritmos)

fviz_contrib(pca.result, choice ="var", axes = 2, top = 20) #variables que mas influyen a la dimension 2

kmeans2 <- kmeans(df_genes, centers = 3)
fviz_cluster(kmeans2, df_genes) #clusterizacion pacientets

#individuos en las dos primeras dimensiones
fviz_pca_ind(pca.result, col.ind = "cos2", gradient.cols = c("blue","yellow", "red"), repel = TRUE)

#codificar una nueva variable
df$metastasisnosi <- as.factor(ifelse( df$extension == "metastasico", "metastasis", "no metastasis"))

df_pca <- as.data.frame(pca.result$x)

df_pca <- df_pca[1:5] #solo estamos analizando los primeros cinco componentes principales

#alguna pareja de componentes permite separar bien los tipos de pacientes ¿?
ggplot(df_pca, aes(x =PC1, y = PC4, color = df$metastasisnosi))+geom_point(size = 3)


#Tabla descriptiva en terciles

#Regresion logistica
df_regresion <- cbind(df_pca, df_nogenes)
y <- as.factor(df$metastasisnosi)

#Hay que convertir df_regresion todo a valores numericos

variables <- colnames(df_regresion)

variables_numericas <- colnames(df_regresion[ ,sapply(df_regresion, is.numeric)]) #variables numericas
variables_no_numericas <- variables[!variables %in% variables_numericas] #variables no numericas 
variables_no_numericas_no_sexo <- variables_no_numericas[variables_no_numericas != "sexo"] #variables de tipo si/no

for (i in variables_no_numericas_no_sexo) {
  
  df_regresion[[i]] <- ifelse (df_regresion[[i]] == "si", 1, 0

  )
  
}

#ahora manualmente se haria el sexo

df_regresion$sexo <- ifelse(df_regresion$sexo == "mujer", 1, 0)

df_regresion <- df_regresion %>% dplyr::select(-X,-id)


x <- as.matrix(df_regresion)

lasso.result <- cv.glmnet(x,y,alpha = 1, familiy = "binomial")





