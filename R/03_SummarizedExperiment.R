# Notas clase 11/02/2026

## Bioconductor
#Todos los paquetes de bioconductor tienen una **"vignettes"** que son documentos en donde prueban los paquetes y su documentación, además de que diario lo evaluan sus paquetes en varios sistemas operativos para saber su funcionalidad y que no esten fallando.

#En la página de [Bioconductor](https://www.bioconductor.org/) puedes revisar los paquetes que hay actualmente, habíendo 4 tipos:
#Los paquetes no son exclusivos de un solo tipo, pueden entre conectarse o tener más de una funcionalidad.
#Los paquetes se instalan de diferente manera entre Bioconductor y de Cran, para instalar un paquete de bioconductor, debes de usar la función de **BiocManager::install**, con el no solo puede instalar paquete de Bioconductor, si no también de Cran.
# Cada paquete es un directorio.
#En Bioconductor es que la gran mayoría de tus funciones tengan ejemplos de uso que se puedan ejecutar, para asegurar su uso, además que tienen un control de versiones, donde los autores dicen que es lo que se fue actualizando.
# A bioconductor lo actualizan dos veces al año, es decir, cada 6 meses en abril y octubre, R se actualiza una vez al año en el mes de abril.

## Lets build our first SummarizedExperiment object
library("SummarizedExperiment")
## ?SummarizedExperiment

## De los ejemplos en la ayuda oficial

## Creamos los datos para nuestro objeto de tipo SummarizedExperiment
## para 200 genes a lo largo de 6 muestras
nrows <- 200
ncols <- 6
## Números al azar de cuentas
set.seed(20210223)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
## Información de nuestros genes
rowRanges <- GRanges(
  rep(c("chr1", "chr2"), c(50, 150)),
  IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
  strand = sample(c("+", "-"), 200, TRUE),
  feature_id = sprintf("ID%03d", 1:200)
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))
## Información de nuestras muestras
colData <- DataFrame(
  Treatment = rep(c("ChIP", "Input"), 3),
  row.names = LETTERS[1:6]
)
## Juntamos ahora toda la información en un solo objeto de R
rse <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  rowRanges = rowRanges,
  colData = colData
)

## Exploremos el objeto resultante
rse


## Comandos para explorar nuestros datos

## 1) Número de genes y muestras
# dim(rse)
## 2) IDs de nuestros genes y muestras
# dimnames(rse)
## 3) Nombres de tablas de cuentas que tenemos (RPKM, CPM, counts, logcounts, etc)
# assayNames(rse)
## 4) El inicio de nuestra tabla de cuentas
# head(assay(rse))
## 5) Información de los genes en un objeto de Bioconductor
# rowRanges(rse)
## 6) Tabla con información de los genes
# rowData(rse) # es idéntico a 'mcols(rowRanges(rse))'
## 7)Tabla con información de las muestras
#colData(rse)

## Ejercicio. Crear un boxplot para el gene_3 por treatment

boxplot(assay(rse)["gene_3", ] ~ colData(rse)$Treatment)

## Descarguemos unos datos de spatialLIBD
sce_layer <- spatialLIBD::fetch_data("sce_layer")
sce_layer

## Revisemos el tamaño de este objeto
lobstr::obj_size(sce_layer)

# Al igual que nuestro objeto rse podemos usar iSEE::iSEE() para explorar los datos.
iSEE::iSEE(sce_layer)
