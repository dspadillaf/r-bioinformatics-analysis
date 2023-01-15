#STAR - Protocols TCGA Analysis
# https://www.costalab.org/wp-content/uploads/2020/11/R_class_D3.html#2_TCGA_data

#BiocManager::install("TCGAbiolinks")
BiocManager::install("TCGAbiolinks")
BiocManager::install("glmnet")
BiocManager::install("caret")
BiocManager::install("survminer")
BiocManager::install("gProfileR")
BiocManager::install("airway")
install.packages("rlang")
install.packages("recipes")


library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")

#TCGA data

GDCprojects = getGDCprojects() 
head(GDCprojects[c("project_id", "name")])
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")

#Descargar datos de conteos par BRCA (Breast Cancer)
query_TCGA = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  access ="open",
  )
#Se verifica los sample para homogenizar los tipo de muestras
brca_res = getResults(query_TCGA)
colnames(brca_res)
head(brca_res$sample_type) 
summary(factor(brca_res$sample_type))

#Se modifica la query dependiendo de esto
query_TCGA = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  access ="open",
  sample.type = "Primary Tumor",
  data.type="Gene Expression Quantification" # tener en cuenta que esto es necesario para que funcione bien el query
)

#Se descargan los datos
setwd("J:/INVESTIGACION/TCGA/BRCA")
GDCdownload(query = query_TCGA)

#Se preparan los datos
tcga_data = GDCprepare(query_TCGA)
    
#En esta libreria la columna 1106 tiene valores NA por lo que se eliminan

tcga_data2 = tcga_data[,-1106]
#Se cambia el nombre para más adelante convertir el objeto en edgeR
head(tcga_data2)
names(tcga_data2@assays)[1] <- "counts"
head(tcga_data2)



#Filtrar datos

#Filtrar para mantener solo genes codificantes de proteinas

library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
ensembl
out <- getBM(attributes=c("ensembl_gene_id_version", "gene_biotype"), 
             filters="ensembl_gene_id_version", values=rownames(tcga_data2), mart=ensembl)
out <- out[match(rownames(tcga_data2), out$ensembl_gene_id_version),] # now in correct order/length.

out_filtrado = out[(out$gene_biotype=="protein_coding" & !is.na(out$gene_biotype)),]

dds_filtrado = tcga_data2[out_filtrado$ensembl_gene_id_version, ]
dim(tcga_data2)
dim(dds_filtrado)

#Se crea el objeto edgeR
library(SummarizedExperiment)
library("edgeR")

d <- SE2DGEList(dds_filtrado)

##Generar la lista de ID Ensembl con su correspondiente GenID o Gene Name
mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
write.csv(listAttributes(mart),"list.csv")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id_version","hgnc_symbol","external_gene_name"),
  values = rownames(dds_filtrado))

write.csv(annotLookup,"list.csv")
dim(dds_filtrado)
dim(annotLookup)
View(annotLookup)


##SUMARIZACIÓN##

do.map=function(a,probe){
  #Filtra por plataforma y deja solo la columna MSU6 (nombres de fila son probeid)
  probeB=probe[names(probe)%in%c("ensembl_gene_id_version","external_gene_name")]#head(probeB)
  #Convierte la matriz en un data.frame
  probeC=as.data.frame(probeB,row.names=probeB$ensembl_gene_id_version)#head(probeC) 
  #Combinacion de data.frames
  b=new.cbind(probeC,as.data.frame(a)) #head(b)    
  #Elimina columnas donde estaban los nombres de probes y de genes para acelerar la conversi?n a matriz
  c=b[,-c(1,2)]#head(c)
  #Convierte de data.frame a matrix para poder agregar nombres de fila repetidos
  c=as.matrix(c)#head(c)
  #Convierte los datos de las columnas a valores numericos
  c=apply(c,2,as.numeric)#head(c)
  #Agrega nombres de fila
  rownames(c)<-b[,2]#head(c)  
  #Para genes que aparecen en m?s de un probe, toma el valor m?ximo entre los probes para cada ensayo
  c=aggregate(c, list(rownames(c)), max)#head(c)
  #Agrega nombres de fila
  rownames(c)<-c[,1]#head(c)
  #Elimina columna donde estaban los nombres de fila
  c=c[,-1]#head(c)  
  return(c)
}

#Funci?n para combinar tablas por nombres de filas
new.cbind <- function(...){
  input <- eval(substitute(list(...), env = parent.frame()))  
  names.orig <- NULL
  nrows <- numeric()
  for (i in 1:length(input))
  {
    nrows[i] <- nrow(input[[i]])
    names.orig <- c(names.orig, colnames(input[[i]])) 
  }  
  idx <- (1:length(input))[order(nrows, decreasing=T)]
  x <- NULL
  for (i in 1:length(input))
  {
    x <- c(x, rownames(input[[idx[i]]]))
  }  
  r <- data.frame(row.names=unique(x))
  for (i in 1:length(input))
  {
    r <- cbind(r, data.frame(input[[i]][match(rownames(r), rownames(input[[i]])),]))
  }  
  colnames(r) <- names.orig  
  return(r)
}

#Reemplaza probes por genes
E=do.map(d$counts, annotLookup)#El n?mero indica la longitud del nombre del gen
dim(d)
dim(E)
View(E)

dds_filtrado <- na.omit(E)
dim(dds_filtrado)

  #sin filtrar
dim(dds_filtrado)

  #filtro expresion > 10000
dds_filtrado=dds_filtrado[1e5>apply(dds_filtrado,1,max),] # eliminar genes con lecturas muy altas
dim(dds_filtrado)

  #filtro con sumatoria de expresion < 10
dds_filtrado=dds_filtrado[rowSums(dds_filtrado) >= 10,]
dim(dds_filtrado)

  #filtro expresion 
dds_filtrado=dds_filtrado[5<apply(dds_filtrado,1,var),] # eliminar genes con varianza menor a 5
dim(dds_filtrado)


#Filtrado por columnas para solo seleecionar triple negativo
#er_status_by_ihc/breast_carcinoma_estrogen_receptor_status = Negative
#pr_status_by_ihc/breast_carcinoma_progesterone_receptor_status  =Negative
#her2_status_by_ihc/lab_proc_her2_neu_immunohistochemistry_receptor_status

tnbc_barcode <- as.matrix(read.table("TNBC.txt"))
#Obtener las columnas que tienen el ID del paciente contenido
colnames(dds_filtrado)[grepl(paste(tnbc_barcode,collapse = "|"),colnames(d))]
d_tnbc <- dds_filtrado[,colnames(d)[grepl(paste(tnbc_barcode,collapse = "|"),colnames(d))]]

d_tnbc
dim(d_tnbc)
#VER Datos

par(mar=c(7,4,2,1), cex.axis=0.4, cex.main = 0.8)
boxplot(d_tnbc, boxwex=0.7, notch=T, outline=FALSE, las=2, main="Box Plot - Muestras TNBC")

#Filtrar columnas(samples) con baja correlación
library(corrplot)

corrplot(cor(d_tnbc), tl.cex = 0.1)

  #https://stackoverflow.com/questions/33809379/eliminate-highly-correlated-columns-but-conserving-my-non-numerical-columns
  #Aplicando filtro

#Filtrado de TNBC con alta correlación
library(caret)
cor_tnbc <- (cor(d_tnbc))
index_tnbc <- findCorrelation(cor_tnbc, .80)
index_tnbc
keep_tnbc <- colnames(cor_tnbc)[index_tnbc]
keep_tnbc
d_tnbc_filtrado <- d_tnbc[,keep_tnbc]

dim(d_tnbc)
dim(d_tnbc_filtrado)

write.csv(d_tnbc, "tnbc_counts.csv")
write.table(d_tnbc_filtrado, "tnbc_correlacionalta_counts.csv")
boxplot(d_tnbc_filtrado, boxwex=0.7, notch=T, outline=FALSE, las=2, main="Box Plot - Muestras TNBC cor > 0.8")
corrplot(cor(d_tnbc_filtrado), tl.cex = 0.1)

###Selección aleatoria de 100 muestras de Breast Cancer que no sean TNBC.###
#Obtener las columnas que NO están contenidos en el archivo de TNBC

d_brca <- dds_filtrado[,colnames(d)[!grepl(paste(tnbc_barcode,collapse = "|"),colnames(d))]] # Selecciona las columnas que no estan en el archivo TNBC
dim(dds_filtrado)
dim(d_tnbc)
dim(d_brca)

              d_brca_100 <- d_brca[, sample(ncol(d_brca), size=100)] # Selecciona aleatoriamente 100 muestras
              dim(d_brca_100)
              write.table(d_brca_100$counts, "brca_100_aleatorio_counts.csv")
              
              d_brca_100 <- read.table("brca_100_aleatorio_counts.csv")
              boxplot(d_brca_100, boxwex=0.7, notch=T, outline=FALSE, las=2, main="Box Plot - 100 Muestras Aleatorias no-TNBC")
              corrplot(cor(d_brca_100), tl.cex = 0.1)
              cor_brca <- (cor(d_brca_100))
              index_brca <- findCorrelation(cor_brca, .80)
              index_brca
              keep_brca <- colnames(cor_brca)[index_brca]
              keep_brca
              d_brca_filtrado <- d_brca_100[,keep_brca]
              
              boxplot(d_brca_filtrado, boxwex=0.7, notch=T, outline=FALSE, las=2, main="Box Plot - Muestras Aleatorias no-TNBC cor > 0.8")
              corrplot(cor(d_brca_filtrado), tl.cex = 0.1)
              dim(d_brca_100)
              dim(d_brca_filtrado)
              write.table(d_brca_filtrado, "brca_correlacionalta_counts.csv")
              
              tnbc_final <- as.data.frame(read.table("tnbc_correlacionalta_counts.csv"))
              brca_final <- as.data.frame(read.table("brca_correlacionalta_counts.csv"))

####FILTRO DE MAS CORRELACIONADOS DE LOS 990 INICIALES
        d_brca_100 <- as.matrix(d_brca)
        dim(d_brca_100)
        cor_brca <- (cor(d_brca_100))
        index_brca <- findCorrelation(cor_brca, .80)
        index_brca
        keep_brca <- colnames(cor_brca)[index_brca]
        keep_brca
        d_brca_filtrado <- d_brca_100[,keep_brca]
        dim(d_brca_100)
        dim(d_brca_filtrado)
        
        ###DE LOS 719 RESTANTES SE SELECCIONAN 100 AL AZAR
        d_brca_100 <- d_brca_filtrado[, sample(ncol(d_brca_filtrado), size=100)] # Selecciona aleatoriamente 100 muestras
        dim(d_brca_100)
        write.table(d_brca_100, "brca_100_aleatorio_counts.csv")
        
        d_brca_100 <- read.table("brca_100_aleatorio_counts.csv")
        boxplot(d_brca_100, boxwex=0.7, notch=T, outline=FALSE, las=2, main="Box Plot - 100 Muestras Aleatorias no-TNBC")
        corrplot(cor(d_brca_100), tl.cex = 0.1)
        cor_brca <- (cor(d_brca_100))
        index_brca <- findCorrelation(cor_brca, .80)
        index_brca
        keep_brca <- colnames(cor_brca)[index_brca]
        keep_brca
        d_brca_filtrado <- d_brca_100[,keep_brca]
        
        boxplot(d_brca_filtrado, boxwex=0.7, notch=T, outline=FALSE, las=2, main="Box Plot - Muestras Aleatorias no-TNBC cor > 0.8")
        corrplot(cor(d_brca_filtrado), tl.cex = 0.1)
        dim(d_brca_100)
        dim(d_brca_filtrado)
        write.table(d_brca_filtrado, "brca_correlacionalta_counts.csv")
        
        tnbc_final <- as.data.frame(read.table("tnbc_correlacionalta_counts.csv"))
        brca_final <- as.data.frame(read.table("brca_correlacionalta_counts.csv"))
        
        

##Se unen nuevamente los datos de conteos ya filtrados##
dim(tnbc_final)
dim(brca_final)
dataCounts = cbind(tnbc_final,brca_final)
dim(dataCounts)
colnames(dataCounts) %in% colnames(tnbc_final)
colnames(dataCounts) %in% colnames(brca_final)

##Se crea el objeto para análisis edgeR Wilcoxon rank-sum test
n_tnbc = ncol(tnbc_final)
n_notnbc = ncol(brca_final)
dataCondition <- c(rep("TNBC",n_tnbc),rep("noTNBC",n_notnbc))
dataCondition <- as.factor(dataCondition)
data_final <- DGEList(counts=dataCounts,group=dataCondition)

colnames(data_final[,data_final$samples$group == "TNBC"]) == colnames(tnbc_final)
colnames(data_final[,data_final$samples$group == "noTNBC"]) == colnames(brca_final)

##Remove rows consistently have zero or very low counts
keep_final <- filterByExpr(data_final)
data_final <- data_final[keep_final,keep.lib.sizes=FALSE]
dim(data_final)

##Normalización
data_final <- calcNormFactors(data_final,method="TMM")
count_norm=cpm(data_final)
count_norm<-as.data.frame(count_norm)
write.csv(count_norm,file="final_norm_count.csv")

par(mar=c(7,4,2,1), cex.axis=0.4, cex.main = 0.8)
boxplot(count_norm, boxwex=0.7, notch=T, outline=FALSE, las=2, main="Box Plot - Normalizado")
corrplot(cor(count_norm), tl.cex = 0.1)

###Wilcoxon rank-sum test Método 1

##Run the Wilcoxon rank-sum test for each gene

pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),dataCondition)
  p=wilcox.test(gene~dataCondition, data)$p.value
  return(p)
})
fdr=p.adjust(pvalues,method = "fdr")

##Calculate the fold-change for each gene

conditionsLevel<-levels(as.factor(dataCondition))
dataCon1=count_norm[,c(which(dataCondition==conditionsLevel[1]))]
dataCon2=count_norm[,c(which(dataCondition==conditionsLevel[2]))]
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))

##Exportar datos
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
fdrThres=0.05
write.table(outRst[outRst$FDR<fdrThres,], file="WilcoxonTest.rst.tsv",sep="\t", quote=F,row.names = T,col.names = T)

###Wilcoxon rank-sum test Método 2
#BiocManager::install("matrixTests")
library("matrixTests")

conditionsLevel<-levels(dataCondition)
dataMem1=count_norm[,c(which(dataCondition==conditionsLevel[1]))]
dataMem2=count_norm[,c(which(dataCondition==conditionsLevel[2]))]
pvalue<-row_wilcoxon_twosample(dataMem1,dataMem2)$pvalue
fdr=p.adjust(pvalue,method = "BH")
fdrcutoff=0.000000001
outputRst=na.omit(as.data.frame(cbind(row.names(count_norm)[which(fdr<fdrcutoff)],pvalue[which(fdr<fdrcutoff)],fdr[which(fdr<fdrcutoff)]),stringsAsFactors=F))
write.table(outputRst,file="WilcoxonTest.rst.tsv",sep="\t", quote=F,row.names = F,col.names = c("Gene","p-value","FDR"))

    ###edgeR
    design <- model.matrix(~dataCondition)
    data_final_edgeR <- estimateDisp(data_final,design) # Considering robust=T when you think your data has potential outlier issue.
    #perform quasi-likelihood F-tests:
    fit <- glmQLFit(data_final_edgeR,design)
    qlf <- glmQLFTest(fit,coef=2)
    fdrcutoff=0.05
    res<-topTags(qlf, n=nrow(dataCounts), adjust.method = "BH", sort.by = "PValue", p.value = fdrcutoff)
    write.table(res, file="edgeR.rst.tsv",sep="\t", quote=F,row.names = T,col.names = T)
    

##Importar datos nuevamente
DEG_Final <- read.table("WilcoxonTest.rst.tsv")
dim(DEG_Final)
DEG_Final

##ORDERNAR GENES CON DE MENOR A MAYOR FDR 

DEG_Final_ordenado <- DEG_Final
colnames(DEG_Final_ordenado) <- DEG_Final_ordenado[1,]
DEG_Final_ordenado <- DEG_Final_ordenado[-1,]
rownames(DEG_Final_ordenado) <- DEG_Final_ordenado[,1]
DEG_Final_ordenado <- DEG_Final_ordenado[,-1]
DEG_Final_ordenado <- DEG_Final_ordenado[order(as.numeric(DEG_Final_ordenado$FDR)),]
write.csv(DEG_Final_ordenado, file = "DEG_Final_ordenado.csv")
dim(DEG_Final_ordenado)

##Filtrar Genes DEG en Matrix
E <- na.omit(count_norm[as.matrix(DEG_Final[1]),])
dim(E)
write.csv(E, file="DEG_norm_final.csv")


##Visualziar 

  ###OPCION 1
  ### Set a color palette
  heat.colors <- brewer.pal(6, "YlOrRd")
  
  ### Run pheatmap
  library(pheatmap)
  pheatmap(E, color = heat.colors, cluster_rows = T, show_rownames=F, border_color=NA, fontsize = 1, scale="row",
           fontsize_row = 1, height=20)

  ###OPCION 2
  highly_variable_lcpm <- as.matrix(E)
  dim(highly_variable_lcpm)
  
  ## Get some nicer colours
  mypalette <- brewer.pal(11,"RdYlBu")
  morecols <- colorRampPalette(mypalette)
  # Set up colour vector for celltype variable
  col.cell <- c("purple","orange")[data_final$samples$group]
  
  # Plot the heatmap
  
  heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="DEG",ColSideColors=col.cell,scale="row")
  
#https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02648-4
#https://rpubs.com/LiYumei/806213

################################################################
# 2. Cálculo de la matriz de similitud.
################################################################

#Una vez conformada la matriz de expresi?n se calcula la similitud por pares de genes.
#A continuaci?n se presentan dos medidas de similitud que capturan la dependencia entre perfiles de expresi?n.

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2.1 Coeficiente de Correlaci?n de Pearson
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#El Coeficiente de correlaci?n de Pearson es la medida m?s usada.
#Ventajas: 
#- Facilidad en la implementaci?n e interpretaci?n.
#Desventajas: 
#- No refleja las asociaciones no lineales entre perfiles, se ve afectado por observaciones extremas y puede tener valores muy altos que dificultan la selecci?n del umbral.

#Calcula el CCP:
#install.packages("minet")
#install.packages("igraph")
#install.packages("infotheo")
library(igraph) #Análisis y manipulación de redes
library(infotheo)#del CRAN

S<-abs(cor(t(E),use =  "pairwise.complete.obs"))
write.table(S, file="matrix_correlacion_ccp.tsv")
#cor(E): calcula correl entre columnas

View(S)

par(mar=c(7,4,2,1), cex.axis=0.4, cex.main = 0.8)
#Revisa histograma de similitudes:
hist(S,xlim=c(0,1), xlab = "Correlación", main = "Histograma de Correlaciones")

#Ejemplo para un par de genes 1:
##Posición de los genes en la matriz de expresi?n:
gen1=which(rownames(S)=="CPA3")
gen2=which(rownames(S)=="MS4A2")
##Perfiles de expresión de los genes:
matplot(t(E[c(gen1,gen2),]), type="l",xlab="Muestras",ylab="Nivel de expresión",las=1, pch = 1:4)
legend("topright", legend = rownames(E)[c(gen1,gen2)], lty = 1:2, col = 1:2)
##Similitud entre los dos genes
S[gen1,gen2]
##Niveles de expresi?n del gen 1 Vs Niveles de expresi?n del gen 2
plot(t(E[gen1,]),t(E[gen2,]),xlab=rownames(E)[gen1],ylab=rownames(E)[gen2],pch=16,col="blue")


#Ejemplo para un par de genes 2:
##Posición de los genes en la matriz de expresi?n:
gen3=which(rownames(S)=="CKAP2L")
gen4=which(rownames(S)=="BUB1")
##Perfiles de expresión de los genes:
matplot(t(E[c(gen3,gen4),]), type="l",xlab="Muestras",ylab="Nivel de expresión",las=1)
legend("topright", legend = rownames(E)[c(gen3,gen4)], lty = 1:2, col = 1:2)
##Similitud entre los dos genes
S[gen3,gen4]
##Niveles de expresi?n del gen 1 Vs Niveles de expresi?n del gen 2
plot(t(E[gen3,]),t(E[gen4,]),xlab=rownames(E)[gen3],ylab=rownames(E)[gen4],pch=16,col="blue")


#Ejemplo para un par de genes 3:
##Posición de los genes en la matriz de expresi?n:
gen5=which(rownames(S)=="FOSL1")
gen6=which(rownames(S)=="DUSP7")
##Perfiles de expresión de los genes:
matplot(t(E[c(gen5,gen6),]), type="l",xlab="Muestras",ylab="Nivel de expresi?n",las=1)
legend("topright", legend = rownames(E)[c(gen5,gen6)], lty = 1:2, col = 1:2)
##Similitud entre los dos genes
S[gen5,gen6]
##Niveles de expresi?n del gen 1 Vs Niveles de expresi?n del gen 2
plot(t(E[gen5,]),t(E[gen6,]),xlab=rownames(E)[gen5],ylab=rownames(E)[gen6],pch=16,col="blue")

#Ejemplo para un par de genes 4:
##Posición de los genes en la matriz de expresi?n:
gen7=which(rownames(S)=="CCNA2")
gen8=which(rownames(S)=="MAD2L1")
##Perfiles de expresión de los genes:
matplot(t(E[c(gen7,gen8),]), type="l",xlab="Muestras",ylab="Nivel de expresi?n",las=1)
legend("topright", legend = rownames(E)[c(gen7,gen8)], lty = 1:2, col = 1:2)
##Similitud entre los dos genes
S[gen7,gen8]
##Niveles de expresi?n del gen 1 Vs Niveles de expresi?n del gen 2
plot(t(E[gen7,]),t(E[gen8,]),xlab=rownames(E)[gen7],ylab=rownames(E)[gen8],pch=16,col="blue")

#Ejemplo para un par de genes 5:
##Posición de los genes en la matriz de expresi?n:
gen9=which(rownames(S)=="CDC20")
gen10=which(rownames(S)=="KIF2C")
##Perfiles de expresión de los genes:
matplot(t(E[c(gen9,gen10),]), type="l",xlab="Muestras",ylab="Nivel de expresi?n",las=1)
legend("topright", legend = rownames(E)[c(gen9,gen10)], lty = 1:2, col = 1:2)
##Similitud entre los dos genes
S[gen9,gen10]
##Niveles de expresi?n del gen 1 Vs Niveles de expresi?n del gen 2
plot(t(E[gen9,]),t(E[gen10,]),xlab=rownames(E)[gen9],ylab=rownames(E)[gen10],pch=16,col="blue")




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2.2 Coeficiente de Información Mutua (CIM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#El coeficiente de información mutua (CIM) es una medida basada en la teor?a de la informaci?n
#Ventajas: 
#- Detecta cualquier tipo de dependencia entre perfiles de expresi?n.
#Desventajas:
#- Requiere discretizar los datos de expresi?n y es de dif?cil implementaci?n. 

#Se usa la libreria minet para calcular la informaci?n mutua entre dos variables (Meyer et al., 2008).

library(minet)

#Los niveles de expresión de cada gen se discretizan en valores enteros seg?n dos par?metros:
# - N?mero de intervalos (nbins): Se puede usar la regla de Sturges: nb = 1+log2(p) (Scott and Sain, 2005).
# - Tipo de discretización: Intervalos de igual amplitud "globalequalwidth" o intervalos de igual frecuencia "equalfreq".
matplot(t(E), type="l",xlab="Muestras",ylab="Nivel de expresi?n",ylim=c(0,15),las=1)

library(infotheo)

Ed=discretize(t(E), disc="globalequalwidth",nbins=ceiling(1+log(ncol(E)))) 

matplot(t(E), type="l",xlab="Muestras",ylab="Nivel de expresi?n",ylim=c(0,15),las=1)
matplot(Ed, type="l",xlab="Muestras",ylab="Nivel de expresi?n discretizado",las=1, ylim=c(1,max(Ed)))


#Ejemplo para un par de genes:
##Antes de discretizar:
matplot(t(E[c(gen1,gen2),]), type="l",xlab="Muestras",ylab="Nivel de expresi?n",ylim=c(0,15),las=1)
##Despu?s de discretizar:
matplot(Ed[,c(gen1,gen2)], type="l",xlab="Muestras",ylab="Nivel de expresi?n discretizado",las=1, ylim=c(1,max(Ed)))

#Calcula la informaci?n m?tua entre pares de genes:
##Se prefiere el estimador de entrop?a de Hausser (2006).  
IM=mutinformation(Ed,method= "shrink")

#La informaci?n m?tua var?a entre [0,8). Se aplica una transformaci?n para calcular el CIM en el intervalo [0,1] (Dionisio et al. 2004):
S1<-sqrt(1-exp(-2*IM))
S1[which(is.na(S1))]<-0

#Revisa histograma de similitudes:
hist(S1,xlim=c(0,1))

#Guarda matriz de similitud:
write.table(S1,"S-IM.txt",sep="\t",col.names=TRUE)


################################################################
# 3. Definición de un umbral de similitud
################################################################

#A partir de la matriz de similitud se selecciona el umbral de similitud (tao).
#El umbral de similitud determina las aristas finales de la red de acuerdo con la funci?n de adyacencia (Zhang & Horvath, 2005).
#La adyacencia entre dos genes i y j viene dada por:
# if (s[i,j]>=tao) {a[i,j]=1} else (a[i,j]=0)
#El umbral se selecciona por el m?todo de Elo et al. (2007).
#El m?todo compara el coeficiente de agrupamiento de la red (Co) con el coeficiente de agrupamiento esperado para una red aleatoria (Cr), a distintos valores de tao.

#Se usa la libreria igraph para calcular los coeficientes de agrupamiento (Csardi et.al. 2006):
library(igraph)

#N?mero de genes:
#S=S1#para hacerlo con IM
n=nrow(S)

#Crea un vector que guarda los coeficientes de agrupamiento locales (Ci) por cada valor de tao:
C=matrix(nrow=n, ncol=100)

#Crea un vector que guarda el grado de nodo de cada gen (ki) por cada valor de tao:
K=matrix(nrow=n, ncol=100)

#Crea un vector de umbrales a ser evaluados:
ltaos=seq(0.01,0.99,by=0.01)

#Calcula el grado de nodo (ki) y el coeficiente de agrupamiento local (Ci) por cada valor de tao:
for(tao in ltaos){
  print(tao)
  ##Matriz de adyacencia:
  A=matrix(0,nrow=n,ncol=n)    
  ##Completa la matriz de adyacencia usando la funci?n de adyacencia:
  for(i in 1:n){  
    A[which(S[,i]>=tao),i]<-1
    A[which(S[,i]<tao),i]<-0
  }
  ##Transforma la matriz A en un objeto igraph:
  A=graph.adjacency(A,mode="undirected",diag=FALSE)
  ##Calcula el Ci de los nodos:
  Cv=transitivity(A,type="local")
  ##Calcula el ki de los nodos:
  Kv=degree(A,loops=FALSE)
  ##Guarda Ci y ki en los vectores C y K respectivamente:
  K[,round(tao*100,0)]<-Kv
  C[,round(tao*100,0)]<-Cv
}

#Calcula el coeficiente de agrupamiento de la red (Co) y el coeficiente de agrupamiento esperado para una red aleatoria (Cr), a distintos valores de tao:
##Define vectores que guardan los coeficientes:
Cr=Co=rep(0,100)
##Para cada valor de tao:
for(i in round(ltaos*100,0)){  
  gn<-which(K[,i]>=1)#Posici?n de los genes conectados en la red
  kn=length(gn)#N?mero de nodos en la red
  k1=1/kn*sum(K[gn,i])#Variable en ecuaci?n 3 (Ver Elo et.al. 2007)
  k2=1/kn*sum(K[gn,i]^2)#Variable en ecuaci?n 3 (Ver Elo et.al. 2007)
  Co[i]=((k2-k1)^2)/(kn*k1^3) #Coeficiente de agrupamiento esperado para una red aleatoria
  if(kn==0){Co[i]=0}#Si no hay nodos conectados: Co=0
  gn<-which(K[,i]>1)#Posici?n de los genes con k1>1.
  kn=length(gn)#N?mero de genes con m?s de una arista en la red
  Cr[i]=1/kn*sum(C[gn,i])#Coeficiente de agrupamiento observado en la red.
  if(kn==0){Cr[i]=0}#Si no hay nodos conectados: Cr=0
}

#Grafica la curva |Cr-Co|:
plot(ltaos,abs(Cr-Co)[ltaos*100],t="l",xlab="Threshold",ylab="|C-Co|")

#Indentifica el primer m?ximo local en la curva |Cr-Co|
##En ocasiones se requiere suavizar la curva para observar el crecimiento continuo:
dif=runmed(abs(Cr-Co),k=3,endrule="constant")[1:100]  
plot(ltaos,dif[ltaos*100],t="l",xlab="Threshold",ylab="|C-Co|")
#(tao=identify(ltaos,dif[ltaos*100],n=1)/100)

#Identificar el máximo local

xplot <- as.vector(ltaos)
yplot <- as.vector(dif[ltaos*100])

plot(xplot, yplot,xlab="Threshold",ylab="|C-Co|")
xmax = xplot[which.max(yplot)]
ymax = yplot[which.max(yplot)]
xmax
ymax
## Add a point
p <- c(xmax, ymax)
points(t(p), pch=16)

## At axes
ps <- diag(2)*p  # get points at axes
points(ps, col="red", pch=c("|", "-"), cex=1:2)

#El máximo se fijo en 0.90 por el resultado del plot
tao = xmax-0.05

################################################################
# 4. Cálculo de la matriz de adyacencia y visualizaci?n de la red.
################################################################

#Finalmente se obtiene la matriz de adyacencia de la red al umbral seleccionado.

#Define matriz de adyacencia:
A=matrix(0,nrow=n,ncol=n)    

#Completa la matriz de adyacencia usando la funci?n de adyacencia:
for(i in 1:n){  
  A[which(S[,i]>=tao),i]<-1
  A[which(S[,i]<tao),i]<-0
}

#Agrega nombres a filas y columnas

rownames(E)
colnames(A)<-rownames(A)<-rownames(E)


#Convierte la diagonal en ceros (red no dirigida):
diag(A)<-0

#Elimina nodos no conectados:
A=A[which(K[,round(tao*100,0)]>0),which(K[,round(tao*100,0)]>0)]

class(A)
dim(A)

#Crea objeto igraph:
A=graph.adjacency(A,mode="undirected",add.colnames=NULL,diag=FALSE)

plot(A)

plot.igraph(A,layout=layout_with_kk,vertex.color = 'lightblue',vertex.size=2,edge.color="darkgreen",vertex.label.font=1,vertex.label.cex=0.6 )
#tydigraph
#cytoskape
#gephi

#Revisa algunas propiedades de la red:
##N?mero de nodos:
length(V(A))
##Aristas:
E(A)

#Guarda listado de aristas:
write.graph(A,"edges.txt",format="ncol")

#Para visualizar, analiza y manipular la red se recomienda el programa Cytoscape (Shannon et.al. 2003)
#En Cytoscape se carga el archivo "edges.txt" desde: File-->Import-->Network-->File

#Extrae clusters:
clusters=clusters(A, mode=c("weak"))

#Revisa el tama?o de los clusters:
clusters$csize

transitivity(A, type="global") 

#A network diameter is the longest geodesic distance (length of the shortest path between two nodes) in the network. 
diameter(A)

#degree
deg <- degree(A, mode="all")
barplot(table(deg),las=2)
plot(A, vertex.size=deg)
plot.igraph(A,layout=layout_with_kk,vertex.color = 'lightblue',vertex.size=deg,edge.color="darkgreen",vertex.label.font=1,vertex.label.cex=0.6 )

hist(deg, breaks=1:vcount(A)-1, main="Histogram of node degree")


deg.dist <- degree_distribution(A, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency")

cd<-centr_degree(A, mode="in", normalized=T)
cd$centralization
cd$res


eigc<-eigen_centrality(A, directed=F, weights=NA)
#The hubs and authorities algorithm developed by Jon Kleinberg was initially used to examine web pages. Hubs were expected to contain catalogs with a large number of outgoing links
hs <- hub_score(A, weights=NA)$vector
hs
sort(hs,decreasing = TRUE)
sort(deg,decreasing = TRUE)

#centr_betw(A, directed=F, normalized=T)
#closeness(A, mode="all", weights=NA) 
#centr_clo(A, mode="all", normalized=T) 



################################################################
# Referencias
################################################################

#Hausser J (2006) Improving Entropy Estimation and the Inference of Genetic Regulatory Networks. 33.
#Dionisio A, Menezes R, Mendes D (2004) Mutual information: a measure of dependency for nonlinear time series. Physica A 344:326?329. doi: 10.1016/j.physa.2004.06.144
#Csardi G, Nepusz T (2006) The igraph software package for complex network research. InterJournal 1:1695.
#Elo LL, J?rvenp?? H, Oresic M, et al. (2007) Systematic construction of gene coexpression networks with applications to human T helper cell differentiation process. Bioinformatics 23:2096?2103.
#Meyer PE, Lafitte F, Bontempi G (2008) minet: A R/Bioconductor Package for Inferring Large Transcriptional Networks Using Mutual Information. BMC Bioinformatics 9:461.
#Scott D, Sain S (2005) Multidimensional Density Estimation. In: Rao CR, Wegman EJ, Solka JL (eds) Handb. Stat. Data Min. Data Vis. Elsevier Science, Amsterdam, pp 229?261
#Zhang B, Horvath S (2005) A general framework for weighted gene co-expression network analysis. Stat Appl Genet Mol Biol 4:Article17
#Irizarry, RA, Hobbs, B, Collin, F, Beazer-Barclay, YD, Antonellis, KJ, Scherf, U, Speed, TP (2003) Exploration, Normalization, and Summaries of High Density Oligonucleotide Array Probe Level Data. Biostatistics .Vol. 4, Number 2: 249-264
#Weston, D. J., Karve, A. a., Gunter, L. E., Jawdy, S. S., Yang, X., Allen, S. M. & Wullschleger, S. D. (2011). Comparative physiology and transcriptional networks underlying the heat shock response in Populus trichocarpa, Arabidopsis thaliana and Glycine max., Plant, cell & environment 34(9): 1488-506. 
#Shannon P,Markiel A, Ozier O, Baliga NS,Wang JT, Ramage D, Amin N, Schwikowski B, Ideker T. 2003. Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome Research 13:2498?2504


#Taller individual
#Sobre un conjunto de datos de expresión
#1. Realice un filtro de genes (genes más activos o DEGs o genes de una categoría de anotación,...)
#3. Construya la red con Información mutua o valor absoluto de Pearson seleccionando el umbral con base en el percentil 90 o el coeficiente de clustering
