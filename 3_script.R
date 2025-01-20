# 1-PRE-PROCESADO --------------------------------------------------------------
# 1.1.Set WEnV -----------------------------------------------------------------
rm(list=ls()) # resetear WEnv
path <- "/Users/pedrofortesgonzalez/Desktop/MASTER_UNIR/11_ALGORITMOS_IA/actividades/Act3_grupal/Material complementario"
setwd(path)

# Librerías y random-seed-------------------------------------------------------
set.seed(1999)
library(readr)
library(dplyr)
library(caret)
library(randomForest)
library(klaR)
library(glmnet)
library(pROC)

# Estructurar el Dataframe a partir de Material_complementario------------------
clases <- read.csv("./1_data/classes.csv", header = FALSE, sep = ";", col.names = c("Muestra", "Clase"), row.names = 1)
columnas <- read_lines("1_data/column_names.txt")
gen_exp <-read.csv("1_data/gene_expression.csv", header = FALSE, sep = ";", col.names = (columnas))
df <- cbind(clases, gen_exp)
df$Clase <- as.factor(df$Clase) # convertir variable clase a factor

# Procesar Valores faltantes----------------------------------------------------
any(sum(is.na(df))) #En principio no hay NAs
any(colSums(df[ ,sapply(df,is.numeric)]) == 0) # Pero sí que hay columnas con todo 0s
ceros <- colnames(df[,2:501])[colSums(df[,2:501]) ==0] # Las almacenamos aquí por si queremos verlas
df <- df[, !colnames(df) %in% ceros] # Y las eliminamos ya que no aportan información y tenemos muchas otras aún por analizar
any(colSums(df[ ,sapply(df,is.numeric)]) == 0) #Comprobamos que ya no hay columnas con todo 0
write.csv(df,"1_data/1_df.csv") # Y lo exportamos (por si acaso)

# Dimensiones del df------------------------------------------------------------
dim(df)

## Tenemos 801 filas y 498 columnas = alta dimensionalidad.
## Debemos escoger un método para elegir columnas ya que, si no, tendremos mucha info (y por tanto ruido) en nuestro análisis
## Probamos a aplicar LASSO para quedarnos con el subconjunto óptimo de variables

