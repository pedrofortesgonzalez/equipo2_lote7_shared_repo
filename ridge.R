#librerias
library(dplyr)
library(caret)
library(randomForest)
library(klaR)
library(glmnet)
library(pROC)

#Base de datos tratada
df <- read.csv("1_data/df.csv")
dim(df) #Tenemos muchas variables > podemos aplicar LASSO para quedarnos con el subconjunto de variables mas optimo

################################################################################
#                               RIDGE
################################################################################

x <- as.matrix(df[ ,3:499])
y <- as.factor(df$Clase)
ridge.result <- cv.glmnet(x,y, family = "multinomial", alpha = 0)
ridge.coef <- coef(ridge.result, s = "lambda.min")

# Convertir los coeficientes a matrices
coeff_AGH <- as.matrix(ridge.coef[["AGH"]])  
coeff_CFB <- as.matrix(ridge.coef[["CFB"]])  
coeff_CGC <- as.matrix(ridge.coef[["CGC"]])  
coeff_CHC <- as.matrix(ridge.coef[["CHC"]])  
coeff_HPB <- as.matrix(ridge.coef[["HPB"]])  

# Convertir las matrices a dataframes, conservando los nombres de las variables
coef_AGH_df <- as.data.frame(coeff_AGH)
coef_CFB_df <- as.data.frame(coeff_CFB)
coef_CGC_df <- as.data.frame(coeff_CGC)
coef_CHC_df <- as.data.frame(coeff_CHC)
coef_HPB_df <- as.data.frame(coeff_HPB)

# Asegurarse de que los nombres de las variables estén en las columnas
colnames(coef_AGH_df) <- "AGH"
colnames(coef_CFB_df) <- "CFB"
colnames(coef_CGC_df) <- "CGC"
colnames(coef_CHC_df) <- "CHC"
colnames(coef_HPB_df) <- "HPB"

# Combinar los dataframes
df_final <- cbind(coef_AGH_df, coef_CFB_df, coef_CGC_df, coef_CHC_df, coef_HPB_df)

df_final <- df_final[2:199, ] #quitar el intercept

#Eliminamos las variables cuya contribucion a todas las clases sea baja
df_final_2 <- df_final[rowSums(abs(df_final)) > abs(0.1) ]
dim(df_final_2) #71 variables me parece mas razonable
variables_elegidas <- rownames(df_final_2) #variables elegidas
#Del df original, solo nos quedamos con la variable categorica y con las variables numericas seleccionadas
df_analisis <- df %>% dplyr :: select(Clase, all_of(variables_elegidas))
dim(df_analisis) 

################################################################################
#                 DIVISION EN ENTRENAMIENTO Y PRUEBA
################################################################################

#convertir a variable categorica
df_analisis$Clase <- as.factor(df_analisis$Clase)

#Separar en conjunto de entrenamiento y conjunto de prueba
index_train <- createDataPartition(df_analisis$Clase, p = 0.8, list = FALSE)
df_train <- df_analisis[index_train, ]
df_test <- df_analisis[-index_train, ]

#Para los metodos de analisis del discriminante,  hay que aplicar una formula:
formula  <- as.formula(paste("Clase ~ ",paste(variables_elegidas, collapse = "+")))


################################################################################
#                               RDA
################################################################################

rda.result <- rda(formula, data = df_train)
rda.prediction <- predict(rda.result, newdata = df_test)
rda.clases <- rda.prediction$class
clases.reales <- df_test$Clase
probabilidades.rda <-  rda.prediction$posterior
matriz.rda <- confusionMatrix(rda.clases, clases.reales)
matriz.rda

################################################################################
#                               SVM
################################################################################

#Como el RDA ha ido bien > parece que los datos tienen relaciones lineales > usamos la variable lineal 

svm.result <- train(Clase ~., data = df_train, method = "svmLinear",
                    trControl = trainControl(method = "cv", number = 5),
                    prob.model = TRUE, tuneGrid = expand.grid(C =seq(1,20,0.5)))

plot(svm.result) #para ir mirando el valor apropiado de C


svm.prediction <- predict(svm.result, newdata = df_test)
svm.probabilidades <- predict(svm.result, newdata = df_test, type = "prob")
matriz.svm <- confusionMatrix(svm.prediction, clases.reales)
matriz.svm

################################################################################
#                               RANDOM FOREST
################################################################################
#Es mas costoso computacionalmente, pero puede merecer la pena usar un metodo mas robusto para compararlo con los otros dos mas simples
randomforest.result <- train(Clase ~., data = df_train, method = "rf",
                             trControl = trainControl(method = "cv", number = 5),
                             tuneLength = 30)

plot(randomforest.result)

randomforest.prediction <- predict(randomforest.result, newdata = df_test)
rf.probabilidades <- predict(randomforest.result, newdata = df_test, type = "prob")
rf.matriz <- confusionMatrix(randomforest.prediction,  clases.reales)
rf.matriz

varImp(randomforest.result) #variables mas importantes en el modelo
varImpPlot(randomforest.result$finalModel)

################################################################################
#                               CURVAS ROC
################################################################################
#Este ejemplo se trata de un problema multiclase:
table(df$Clase)

#No se pueden usar curvas ROC ni PR como tal > solo cuando hay dos clases
#Pero si que podemos usar multiclass.roc
#Esta funcion calcula las curvas ROC para cada clase y luego genera un promedio.
#No se puede graficar, pero se puede ver el area bajo la curva y se interpreta igual que las curvas ROC y PR

curva_rda <- multiclass.roc(df_test$Clase, as.matrix(probabilidades.rda))
curva_svm <- multiclass.roc(df_test$Clase, as.matrix(svm.probabilidades))
curva_rf <- multiclass.roc(df_test$Clase, as.matrix(rf.probabilidades))


metodos <- c("RDA", "SVM","RF")
areas <- c(curva_rda$auc,curva_svm$auc, curva_rf$auc)

analisis <- data.frame(modelos = metodos,  AUC = areas) #Pedro, por si quieres usar este df de base para incorporar el resto de parametros





