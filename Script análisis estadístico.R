#En primer lugar cargamos las librerías que vamos a utilizar

library(SummarizedExperiment)    #Para trabajar con esta clase
library(metabolomicsWorkbenchR)  #Para descargar los datos de Metabolomics Workbench
library(POMA)                    #Para llevar a cabo el pretratamiento de los datos
library(ggplot2)                 #Para gráficos
library(ggraph)                  #Para gráficos
library(plotly)                  #Para gráficos
library(factoextra)              #Para gráficos de cargas del loading plot
library(ComplexHeatmap)          #Para crear un agrupamiento con heatmap

#Descargamos los datos en como el objeto se de la clase SummarizedExperiment desde la base de datos Metabolomics Workbench

se <- do_query(
  context = 'study',
  input_item = 'study_id',       
  input_value = 'ST000910',      #Indicamos que busque el id del estudio en cuestión
  output_item = 'SummarizedExperiment'
)

#Comprobamos que el objeto se es de la clase correcta y vemos sus atributos

class(se)
rowData(se)                      #Información de las filas (dataframe)
colData(se)                      #Información de las columnas (dataframe)
metadata(se)                     #Metadatos relevantes del estudio
head(assay(se))                  #matriz de resultados

#Guardamos el objeto en formato .Rda

save(se, file = "SummarizedExperiment.rda")

#Seguimos el workflow de pre procesado de datos de POMA:
#Primero comprobamos la existencia de los valore nulos y, en caso de haberlos, los eliminamos con PomaImpute

imputed <- se %>% 
  PomaImpute(method = "knn", zeros_as_na = TRUE, remove_na = TRUE, cutoff = 20)

#Ahora normalizamos los datos para una mejor visualización

normalized <- imputed %>% 
  PomaNorm(method = "log_pareto")

#Hacemos un boxplot usando POMA con los datos normalizados para visualizarlos

PomaBoxplots(normalized, x = "features")

#Visualizamos las muestras para detectar outliers con PomaOutliers()

PomaOutliers(normalized, outcome = "Diagnosis")

#Finalmente guardamos únicamente los datos que no son diagnosticados como outliers

pre_processed <- PomaOutliers(normalized, outcome = "Diagnosis")$data
pre_processed                    #Así tenemos los datos ya pre procesados. Sólo hemos perdido una muestra

#Preparamos los datos para llevar a cabo el PCA

matrix <- assay(pre_processed)                        #Extraemos la matriz del objeto de datos preprocesados
fac <- colData(pre_processed)                         #Extraemos la información de las columnas del objeto de datos preprocesados
metabolites <- rowData(se)                            #Extraemos la información de las filas del objeto de datos preprocesados

#Por alguna razón no se han guardado los nombres de las filas en la matriz
#Lo incluimos a partir del dataframe que acabamos de obtener

rownames(matrix) <- metabolites$metabolite_id 

#Para hacer las pruebas con los componentes principales que queramos hacemos una función en la que se puedan especificar:

plotPCA <- function ( X, scale=FALSE, firstPC=1, secondPC=2) #Opción de poner la escala y seleccionar los PC 
{
  pcX<-prcomp(X, scale=scale)                         #Creamos la matriz de componentes principales con prcomp()
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)     #Calculamos el porcentaje de varianza explicado para todos los PC
  xlab<-c(paste("PC",firstPC,loads[firstPC],"%"))     #String para el eje x con el porcentaje de varianza explicada
  ylab<-c(paste("PC",secondPC,loads[secondPC],"%"))   #String para el eje y con el porcentaje de varianza explicada
  
  #Creamos con ggplot el scatter plot para los PC indicados con distinguiendo en función de los factores diagnosis (color) y sex (forma)
  
  scatter <- ggplot(pcX$x, aes(x = pcX$x[,firstPC], y = pcX$x[,secondPC], color = fac$Diagnosis, shape = fac$Sex)) +
    geom_point(size = 2) +                            #Tamaño puntos
    scale_color_viridis_d() +                         #Escala de colores
    theme_minimal() +                                 #Tema
    xlab(xlab) +                                      #Texto eje x
    ylab(ylab) +                                      #Texto eje y
    labs(shape="Sexo", color="Diagnosis")             #Nombres para la leyenda
  
  #Creamos el diagrama de cargas con la siguiente función
  
  loadings <- fviz_pca_var(pcX, axes = c(firstPC, secondPC), col.var="contrib", repel=TRUE)
  
  #Imprimimos los gráficos
  
  print(scatter)
  print(loadings)
}

#Aplicamos la función a nuestros datos para los dos primeros PCs con la matriz transpuesta

plotPCA(t(matrix), scale=TRUE, firstPC=1, secondPC=2)

#Aplicamos la función a nuestros datos para los dos siguientes PCs con la matriz transpuesta

plotPCA(t(matrix), scale=TRUE, firstPC=3, secondPC=4)

#Pasamos a crear el heatmap:
#Primero nombramos a las columnas de la matriz por el factor diagnosis

colnames(matrix) <- fac$Diagnosis

#Creamos el heatmap

Heatmap(t(matrix), 
        name = " ",                                                    #Leyenda
        column_title = "Metabolitos", row_title = "Muestras",          #Nombre de los ejes
        row_names_gp = gpar(fontsize = 6)                              #Tamaño del texto
)


