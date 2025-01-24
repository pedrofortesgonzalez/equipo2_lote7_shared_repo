###############################################################################
#                               LIBRERIAS
###############################################################################

library(dplyr)
library(gtsummary)
library(stats)
library(factoextra)
library(MASS)
library(ggplot2)
library(glmnet)
library(gridExtra)

###############################################################################
#                               DATASET
###############################################################################

df <- read.csv("1_data/Dataset expresión genes.csv")
df_genes <- df %>% dplyr::select(starts_with("AQ"))
df_nogenes <- df %>% dplyr::select(-starts_with("AQ"))


###############################################################################
#                               PCA
###############################################################################

set.seed(1234)
pca.result <- prcomp(df_genes, center = TRUE, scale = TRUE )


#####Graficos PCA ¿cuantos Componentes principales elegimos?
#Si miramos los eigenvalues:

eigenvalues <- get_eigenvalue(pca.result)
eigenvalues #parece que hay que coger hasta la dimension 5 (70% varianza)> miramos grafico
fviz_eig(pca.result, addlabels = TRUE) #hasta la dimension 4 ¿?
#Si aplicaramos la regla del codo, en la dimension 2-3 ya se aplana la grafica, pero no llegamos al 70% de varianza

###############################################################################
#                               GRAFICOS PCA
###############################################################################

fviz_pca_var(pca.result, col.var = "cos2", gradient.cols = c("blue","yellow", "red"), repel = TRUE, axes = c(1,2),
             title = "Cos2 de variables en PC1 y PC2") #variables segun calidad cos2
#cos2 indica que tan bien esta representada las variables por los componentes principales
# (como de bien "caben" las variables en el espacio de menor dimension) > parece que casi todas las variables se ajustan bien al epacio de menor dimension,
#menos algunos genes como ADIPOQ o NOX5

#Clusterizacion
kmeans <- kmeans(t(df_genes), centers = 2)

clusters <- list()

for (i in 1:4) {
  
  cluster.plot <- fviz_pca_var(pca.result, col.var = kmeans$cluster, gradient.cols = c("blue","green", "red"), axes = c(i,i+1), 
                               legend.title ="Clusterizacion", repel = TRUE) 
  clusters[[i]] <- cluster.plot
}

todos.clusters <- grid.arrange(grobs = clusters, nrow = 3)
todos.clusters


#Nombre cluster verde > metabolismo energetico
#Nombre cluster azul > respuesta inmune e inflamatoria
#Nombre cluster rojo > metabolismo celular y respuesta inmune 

# comentario Clusters:
## el cluster rojo contiene genes tanto del metabolismo celular como de la respuesta inmune. Debemos mejorar la clusterización (tema 4 de Algoritmos)
#Es dificil hacer la clusterizacion, casi todos los genes parece que podrian pertencer al mismo bloque, 
#que se relacionan de manera negativa con l adimension 1 con bastante fuerza. 
#tambien hay otro conjunto pequeño de genes que se relacionna de manera positiva con dim1, pero son los de cos2 bajo. 
#realmente solo sale un grupo de genes muy heterogeneo y dificl de agrupar, y luego los dos genes que no se ajustn bien al nuevo espacio de menor dimension
#El resto de comparaciones no muestran esa tendencia tan clara como con la dimension 1


#Para darle nombre a los componentes, se puede ver la contribucion de las variables a las dimensiones

graficos <- list()
for (i in 1:5) {
  grafico <- fviz_contrib(pca.result, choice ="var", axes = i, top = 5) 
  graficos[[i]] <- grafico
}

graficos_varianza <- grid.arrange(grobs = graficos, nrow = 3)
#Componente 1 > inflamacion, respuesta inmune, estres celular, metabolismo (pero la relacion con las variables es negativa, recordar)
#Podria ser algo asi como disminucion de la respuesta imune, estres celular

#Componente 2: alteraciones en la inmunorregulacion (aumento y disminucion de la expresion)

#Componente 3: alteraciones metabolismo energetico y señalizacion en estres celular

#Componente 4: alteraciones inmunidad innata y homeostasis

#Componente 5:alteraciones en inflamacion y regeneracion celular


###Graficos pacientes
kmeans2 <- kmeans(df_genes, centers = 3)
fviz_cluster(kmeans2, df_genes) #clusterizacion pacientes > hay muchos solapamientos
#probablemente los pacientes, aunque tienen enfermedades/rasgos diferentes, como los genes
#pueden estar implicados en diferentes rutas metabolicas y procesos tumorales, puede haber pacientes
#con diferente enfermedad pero que tengan un patron genico similar o coincidente en  algunos puntos


#individuos en las dos primeras dimensiones
fviz_pca_ind(pca.result, col.ind = "cos2", gradient.cols = c("blue","yellow", "red"), repel = TRUE)
#muchos pacientes parecen ajustarse bien al nuevo espacio reducido, pero algunos como el 41, 50,10,34,14 no. 

#codificar una nueva variable > variable para categorizar a los pacientes segun su tumor
df$metastasisnosi <- as.factor(ifelse( df$extension == "metastasico", "metastasis", "no_metastasis"))

df_pca <- as.data.frame(pca.result$x)

df_pca <- df_pca[1:5] #solo estamos analizando los primeros cinco componentes principales

#alguna pareja de componentes permite separar bien los tipos de pacientes ¿?

componentes <- c("PC1", "PC2", "PC3", "PC4", "PC5")
componentes.plots <- list()

for (i in 1:(length(componentes) - 1))  {
  grafico <- ggplot(df_pca, aes_string(x =componentes[i], y = componentes[i+1], color = df$metastasisnosi))+geom_point(size = 3)
  componentes.plots[[i]] <- grafico
}

grid.arrange(grobs = componentes.plots, nrow = 3)

#No parece que se  puedan separar bien los pacientes > igual no es el mejor metodo para la reduccion
#de la dimensionalidad ¿poner algo mas?


###############################################################################
#                               TABLA DESCRIPTIVA TERCILES
###############################################################################


###############################################################################
#                               REGRESION  LOGISTICA
###############################################################################
#Primero aplicamos LASSO
df_regresion <- cbind(df_pca, df_nogenes)
#Necesitamos que todas las variables sean numericas
variables <- colnames(df_regresion)
variables_numericas <- colnames(df_regresion[ ,sapply(df_regresion, is.numeric)]) #variables numericas

#Escalado variables numericas:
df_numerico <- df_regresion %>% dplyr::select(variables_numericas)
df_numerico <- scale(df_numerico)
df_regresion <- df_regresion %>% dplyr::select(-variables_numericas)
df_regresion <- cbind(df_regresion, df_numerico)

variables_no_numericas <- variables[!variables %in% variables_numericas] #variables no numericas 
variables_no_numericas_no_sexo <- variables_no_numericas[variables_no_numericas != "sexo"] #variables de tipo si/no

for (i in variables_no_numericas_no_sexo) {
  
  df_regresion[[i]] <- ifelse (df_regresion[[i]] == "si", 1, 0

  )
  
}

#ahora manualmente se haria el sexo

df_regresion$sexo <- ifelse(df_regresion$sexo == "mujer", 1, 0)

df_regresion <- df_regresion %>% dplyr::select(-X,-id)

y <-  as.factor(df$metastasisnosi)
x <- as.matrix(df_regresion)

grid <- 10^seq(10,-2, length=100)

lasso.result <- glmnet :: cv.glmnet(x, y, family = "binomial", lambda = grid, alpha = 1)
lambda_min <- lasso.result$lambda.min
lasso.coef <- coef(lasso.result, s = "lambda.min")
lasso.coef <- as.data.frame(as.matrix(lasso.coef))
lasso.coef$s1 <- round(lasso.coef$s1, 5)

View(lasso.coef)








