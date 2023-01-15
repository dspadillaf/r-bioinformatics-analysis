

#Limpia informaci?n actual
rm(list=ls(all=TRUE))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affy")
library(affy)
BiocManager::install("GEOquery")
BiocManager::install("CLL")
BiocManager::install("Rcpp")
BiocManager::install("affydata")
BiocManager::install("vsn")
BiocManager::install("gcrma")
BiocManager::install("latticeExtra")
BiocManager::install("RColorBrewer")
BiocManager::install("oligo")
BiocManager::install("pd.e.coli.2")


#library corrplot y Matrix del CRAN
install.packages("Matrix")
install.packages("corrplot")
install.packages("pheatmap")


#1.Cargar datos
#2.Verificar calidad / Filtro genes y muestras
#3. Normalizar

#1.Cargar datos
library(GEOquery)
gset <- getGEO("GSE127710", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
class(gset)
gset
# set parameters and draw the plot

#dev.new(width=4+dim(gset)[[2]]/5, height=6)
#par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
#title <- paste ("GSE127710", '/', annotation(gset), " selected samples", sep ='')
boxplot(exprs(gset), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
class(exprs(gset))
library(corrplot)
corrplot(cor(exprs(gset)))
hist(exprs(gset)[,1])
write.table(exprs(gset),"tabla1.csv")
###############################

#en librería affy con Readaffy se crea un affybatch
#Importaci?n de datos con formato Affymetrix
#library(affy)
#myAB=ReadAffy()#sin argumentos lee todos los archivos CEL del directorio (en caso de que tengan archivos ya bajados del GEO den NCBI)
#myAB=ReadAffy(filenames=c("GSM770862.CEL","GSM770863.CEL"),full.names=TRUE)

################################
             
#con read.celfiles
setwd("~/2022/cursos/Genómica/Talleres/Rarchivos")
library(oligo)
library(pd.e.coli.2)
?read.celfiles

myAB2=read.celfiles("GSM770862.CEL","GSM770863.CEL")
myAB2
class(myAB2)
#verificar calidad

##############################
#read.table("tabla.txt")  (ver taller 1)
##########################

#ejemplo con datos de leucemia linfocitica cronica del paquete CLL, el cual ya es de tipo Affybatch
#2.Verificar calidad
library(CLL)
data(CLLbatch)
class(CLLbatch)
CLLbatch
#acceder a toda la informacion de anotacion contenida en el objeto
sampleNames(CLLbatch)
str(CLLbatch)

#informacion adicional en otra tabla de datos (estado del paciente, si enfermo o control)
data(disease)
head(disease)
#se puede incluir en la parte anotacion del objeto


#Determinar y controlar la calidad

library(lattice)
library(Matrix)
library(genefilter)
library(gcrma)

#Calidad con descripcion de datos:

class(CLLbatch)
corrplot(cor(exprs(CLLbatch)))
corrplot.mixed(cor(exprs(CLLbatch)), lower="number",upper="circle",number.cex=0.6,tl.cex=0.6)
boxplot(exprs(CLLbatch),las=2)
hist(exprs(CLLbatch[,1]))

#DENDROGRAMA ES OPCIONAL, PORQUE MISMA INFO QUE CORRELOGRAMA
dd= dist2(log2(exprs(CLLbatch)))
View(dd)
diag(dd)=0
dd.row<-as.dendrogram(hclust(as.dist(dd)))
row.ord<-order.dendrogram(dd.row) #organizacion de filas y columnas para el heatmap que refleja las distancias entre muestras

library(RColorBrewer)
library(latticeExtra)
legend=list(top=list(fun=dendrogramGrob,args=list(x=dd.row,side="top")))
lp=levelplot(dd[row.ord,row.ord],xlab="", ylab="",legend=legend, las=3)
lp

#quitar alguna muestra?
#ejemplo:
exprsdef<-exprs(CLLbatch)[,-2]#elimina columna 2

#faltaria filtrar genes si es necesario (por bajos conteos, por ejemplo)

#3. Normalizacion

#rma de la libreria affy o de la libreria oligo


#################################
class(CLLbatch)
#convertir a expressionset y hacer normalizacion
tmp <- new ("ExpressionSet", phenoData = phenoData(CLLbatch) ,
            featureData = featureData(CLLbatch), experimentData =
            experimentData(CLLbatch), annotation =
            annotation(CLLbatch),  assayData= assayData(CLLbatch))
class(tmp)
class(CLLbatch)
#eset <- rma(CLLbatch)#expression set
########################

library(vsn)#m?todo vsn de Huber 
citation("vsn")

tmpmat<-exprs(tmp)
class(tmpmat)
dim(tmpmat)
tmpmat<-tmpmat[1:20000,]
tmpnorm = justvsn(tmpmat)
meanSdPlot(tmpmat[1:20000,], ranks=TRUE)
meanSdPlot(tmpnorm[1:20000,], ranks=TRUE)



data("kidney")
class(kidney)
dim(kidney)
boxplot(exprs(kidney))
par(mfrow=c(1,2))
hist(exprs(kidney)[,1],main="green")
hist(exprs(kidney)[,2],main="red")
par(mfrow=c(1,1))
meanSdPlot(kidney, ranks=TRUE)

#se observa que sd aumenta mucho en los valores altos
nkid2 = justvsn(kidney)
meanSdPlot(nkid2, ranks=TRUE)
#se ha estabilizado la varianza

par(mfrow=c(1,1))
boxplot(exprs(nkid2))
par(mfrow=c(1,2))
hist(exprs(nkid2)[,1],main="green")
hist(exprs(nkid2)[,2],main="red")
#Tambien se ha mejorado la asimetria


par(mfrow=c(1,1))

#Filtro de genes haciendo en crudos y visualizando en datos filtrados
select = (0==rowSums(exprs(kidney)<=0))#filtro genes
plot(log2(exprs(kidney)[select, ]), main = "a) raw", pch = ".", asp=1)
plot(exprs(nkid2), main = "b) vsn", pch = ".", asp=1,col=ifelse(select, "black", "orange")) # en naranja los genes que no est?n en el filtro

#genes con valores de cero expresion, que se podrian eliminar
R=exprs(kidney)
Rfilt=exprs(kidney)[select, ]
Rfiltnorm=exprs(nkid2)[select, ]
#write.table

class(R)
plot(R, pch=".")
abline(a=0,b=1,col="blue")
#correlacion entre dos muestras
#para varias se revisa con corrplot

#ejemplo con otros datos RNAseq-RPKM
#en carpeta Rarchivos, editar linea siguiente seg?n cada computador
setwd("~/2022/cursos/Genómica/Talleres/Rarchivos")
RNA<-read.table("RNAseq_Athaliana.txt",h=T)
View(RNA)#datos de conteos ya normalizados por RPKM
dim(RNA)
names(RNA)
#par(mfrow=c(1,1))
boxplot(RNA[,2:26],las=2)
gennames<-RNA$GenMicroAT
#RNAdef<-RNA[,2:26]
#sin embargo muchos at?picos, hay que aplicar vsn
hist(RNA[,2])
class(RNA)
RNA<-as.matrix(RNA[,2:26])
#library(vsn)
RNA2 = justvsn(RNA)
boxplot(RNA2,las=2)
meanSdPlot(RNA, ranks=TRUE)
meanSdPlot(RNA2, ranks=TRUE)
library(pheatmap)
pheatmap(RNA2)
#ahora estan mucho mas comparables y con varianza estabilizada
correl<-cor(RNA2)
View(correl)
symnum(correl)
?symnum
#library(corrplot)
corrplot(correl)


#cuantos grupos de muestras parece haber?
#Tissues were collected at three different times post-inoculation (6, 12 and 24 hours)

#Para uso de otras librer?as Bioconductor convertir a ExpressionSet
class(RNA2)
as.matrix(RNA2)->matriz
eset <- ExpressionSet(assayData=RNA2)
class(eset)

#Taller para datos espec?ficos de affymetrix
#http://bioinf.wehi.edu.au/marray/ibc2004/lab2/lab2.html



################Taller individual#############

#Sobre los datos GSE5752 o GSE4116 (o datos de su escogencia) realice los siguientes pasos:
#1. Consulte y resuma la descripcion del experimento en NCBI (GEO datasets) y describalo brevemente
#2. Haga estadísticas descriptivas y de calidad  (boxplot, heatmap, correlacion entre replicas etc.); eliminar las muestras problematicas; eliminar genes que considere poco informativos.
#3. Normalice los datos con el metodo vsn y verifique estabilizacion de la varianza con meansdplot
#4. Proponga una pregunta biologica para el an?lisis de estos datos
