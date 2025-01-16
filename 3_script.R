library(readr)
#Estructurar el Dataframe
clases <- read.csv("1_data/classes.csv", header = FALSE, sep = ";", col.names = c("Muestra", "Clase"), row.names = 1)
columnas <- read_lines("1_data/column_names.txt")
gen_exp <-read.csv("1_data/gene_expression.csv", header = FALSE, sep = ";", col.names = (columnas))
datos <- cbind(clases, gen_exp)
#Imputación
sum(is.na(datos)) #En principio no hay NAs
datos$Clase <- as.factor(datos$Clase)
str(datos)
ceros <- colnames(datos[,2:501])[colSums(datos[,2:501]) ==0] #Que genes dan 0
datos_imputados <- datos[, !colnames(datos) %in% ceros] #Quitamos ceros
datos_imputados <- datos[, colSums(datos[,2:501]) != 0] #Quitamos ceros (más simple)
datos #Ha removido 3 variables

#######################Aporte Cris#############################################
#al ejecutar esto
any(colSums(datos_imputados[ ,sapply(datos_imputados,is.numeric)]) == 0) ##sigue habiendo columnas con 0
genes <- gen_exp[ ,colSums(gen_exp) != 0] #sobre df con todo numeros, quitamos las columnas con todo 0
genes <- scale(genes) #escalamos tambien
df <- cbind(clases, genes)
any(is.na(df)) #no datos faltantes
any(colSums(df[ ,sapply(df,is.numeric)]) == 0) #ya no hay columnas con todo 0
dim(df) #ahora si que se han eliminado las variables con todo 0
write.csv(df,"1_data/df.csv")

