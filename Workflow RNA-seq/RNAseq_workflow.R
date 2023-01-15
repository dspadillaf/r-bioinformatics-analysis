#Taller 3 de entregar estadística genómica
rm(list=ls(all=TRUE))
setwd("G:/Universidad Nacional Maestria/2do Semestre/Estadistica Genomica/Taller 3")

#1. Cargar datos
#BiocManager::install("farms")
#BiocManager::install("farms")


  #GEO
  library(GEOquery)
  gset <- getGEO("GSE61697", GSEMatrix =TRUE, getGPL=FALSE)
  if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  write.AnnotatedDataFrame(phenoData(gset),"result/GSE61697_pheno.txt")
  
  #Affy
  library(affy)
  celfilesdir <- "data/GSE61697/celfiles/"
  celfileslist <- list.celfiles(celfilesdir)
  affyset <- ReadAffy(filenames = celfileslist, celfile.path = celfilesdir)
  
  #Obtener los nombres de accesion y phenodata apartir de ExpressionSet
    library(magrittr)
    library(purrr)
    library(stringr)
    pd <- pData(gset)
    pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)
    pd$cel_file
  
  #Cambiar los nombres de los samples al codigo de accesion
    sampleNames(affyset)[sampleNames(affyset)==pd$cel_file]  <- pd$geo_accession[pd$cel_file==sampleNames(affyset)]
    sampleNames(affyset)
  
  #Agregar phenodata
    phenoData(affyset) <- phenoData(gset)[sampleNames(affyset),]
    phenoData(affyset)

ex<-exprs(gset)    

#2. probar calidad
    
    #Affyset
    par(mar=c(7,4,2,1))
    boxplot(exprs(affyset), boxwex=0.7, notch=T, outline=FALSE, las=2)
    class(exprs(gset))
    library(corrplot)
    corrplot(cor(exprs(gset)))
    hist(exprs(gset)[,1])
    write.table(exprs(gset),"result/taller3_crudos.csv")

  #Sumarization
    library(farms)
    eset <- qFarms(affyset, laplacian=TRUE)
    exprs(eset)
    
    #Filtrado de genes y pruebas
    INIs <- INIcalls(eset)
    summary(INIs)    
    
    I_data <- getI_Eset(INIs)
    I_data
    exprs(I_data)
    
    gset <- I_data
      
#Normal
par(mar=c(7,4,2,1))
boxplot(exprs(gset), boxwex=0.7, notch=T, outline=FALSE, las=2)
class(exprs(gset))
library(corrplot)
corrplot(cor(exprs(gset)))
hist(exprs(gset)[,1])
write.table(exprs(gset),"result/RNA_exp_norm.csv")


#Log
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex_log <- log2(ex) }

par(mar=c(7,4,2,1))
title <- paste ("GSE61697", "/", annotation(gset), sep ="")
boxplot(ex_log, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
  
#2. Filtrar datos de baja expresion y anomalos

varLabels(gset)
table(gset$`cell population:ch1`)

#Filtrar quitando las columnas anomalas
filtro<-colnames(gset)[gset$"geo_accession"!="GSM1511427" & gset$"geo_accession"!="GSM1511432" & gset$"geo_accession"!="GSM1511438" & gset$"geo_accession"!="GSM1511439"]
gset.filtrado <-gset[,filtro]
table(gset.filtrado$geo_accession)

filtrado <- exprs(gset.filtrado)

  #Graficar nuevo descartando
  qx <- as.numeric(quantile(filtrado, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { filtrado[which(filtrado <= 0)] <- NaN
  ex_log_filtrado <- log2(filtrado) }
  
  par(mar=c(7,4,2,1))
  title <- paste ("GSE61697", "/", annotation(gset), sep ="")
  boxplot(ex_log_filtrado, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)

  #Filtro sugerido
  e.mat <- 2**(exprs(gset.filtrado))
  
  # look at mean, sd, & cv for each gene across arrays
  gene.mean <- apply(e.mat,1,mean)
  gene.sd <- apply(e.mat,1,sd)
  gene.cv <- gene.sd/gene.mean
  
  # make plots
  library(geneplotter); library(RColorBrewer)
  blues.ramp <- colorRampPalette(brewer.pal(9,"Blues")[-1])
  dCol <- densCols(log(gene.mean),log(gene.sd),
                   colramp=blues.ramp)
  par(mfrow=c(1,1))
  plot(gene.mean,gene.sd,log='xy',col=dCol,pch=16,cex=0.1)
  abline(v=100,lwd=3,col='red')
  hist(log(gene.cv),main=NA)
  abline(v=log(.7),lwd=3,col='red')
  
  # filter: keep genes with cv between .7 and 10,
  # and where 20% of samples had exprs. > 100
  library(genefilter)
  f1 <- pOverA(0.20,100)
  f2 <- cv(0.7,10)
  
  #Medida de presión superior a 200 en al menos 5 muestras
  f3 <- kOverA(5, 200)

  #Aplicar filtros
  ffun_combined <- filterfun(f1, f2)
  t.fil <- genefilter(e.mat, ffun_combined)
  small.eset <- log2(e.mat[t.fil,])
  
  ffun_combined <- filterfun(f3)
  t.fil <- genefilter(small.eset, ffun_combined)
  small.eset <-small.eset[t.fil,]
  
  # apply filter, and put expression back on log scale
  dim(e.mat)
  dim(small.eset)
  
  # find the gene names
  gn.keep <- rownames(e.mat)
  gn.keep

#Baja expresión
Rfilt=small.eset[rowSums(small.eset[])>4,]
write.table(Rfilt,"result/RNA_exp_fil.csv")
boxplot(Rfilt, boxwex=0.7, notch=T, outline=FALSE, las=2)


##3. Normalización
matrix_fil <- read.table("result/RNA_exp_fil.csv", header=TRUE ,sep =" ", row.names =  1)
RNA_exp_fil <-new("ExpressionSet", exprs=as.matrix(matrix_fil))

library(vsn)#m?todo vsn de Huber 
citation("vsn")

RNA_exp_norm<-exprs(RNA_exp_fil)
class(RNA_exp_norm)
dim(RNA_exp_norm)
RNA_exp_norm = justvsn(RNA_exp_norm)

#Guardar resultado
write.table(RNA_exp_norm,"result/RNA_exp_norm.csv")

#Comprobar normalización

  #Affy
  RNA_exp_fil<-exprs(affyset)
  RNA_exp_norm <- exprs(gset)
  plotMDS(d, col = as.numeric(group))
  
meanSdPlot(RNA_exp_fil, ranks=TRUE)
meanSdPlot(RNA_exp_norm, ranks=TRUE)

#Nota: se esperaria que cambiara la varianza a una escala muy pequeña, pero se debe preguntar si se sigue filtrando o no

boxplot(RNA_exp_norm,las=2)
par(mfrow=c(1,2))
hist(RNA_exp_norm[,3],main="green")
hist(RNA_exp_norm[,4],main="red")

##Determinar control de la normalización
par(mfrow=c(1,1))
corrplot(cor(RNA_exp_norm))
corrplot.mixed(cor(RNA_exp_norm), lower="number",upper="circle",number.cex=0.6,tl.cex=0.6)
boxplot(RNA_exp_norm,las=2)
hist(RNA_exp_norm[,1])

gset$platform_id
##3.5. probes to genID
  require("biomaRt")
  mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  annotLookup <- getBM(
    mart=mart,
    attributes=c(
      "affy_hg_u133_plus_2",
      "ensembl_gene_id",
      "gene_biotype",
      "external_gene_name"),
    filter = "affy_hg_u133_plus_2",
    values = rownames(exprs(gset)), uniqueRows=TRUE)

  annotLookup
  
  #Guardar con Ensembl ID
  indicesLookup <- match(rownames(gset), annotLookup$affy_hg_u133_plus_2)
  head(annotLookup[indicesLookup, "ensembl_gene_id"])  
  dftmp <- data.frame(rownames(RNA_exp_norm), annotLookup[indicesLookup, c("affy_hg_u133_plus_2", "ensembl_gene_id")])

  #Comprobacion
  head(dftmp, 20)
  table(dftmp[,1] == dftmp[,2])
  
  rownames(RNA_exp_norm) <- paste(annotLookup[indicesLookup, "external_gene_name"])
  RNA_exp_norm
  
##4. MA plot
inaive = seq(1, 6)
itscm = seq(7, 12)
itcm =seq(13, 17)
item =seq(18, 23)

M_naive_tscm= RNA_exp_norm[,inaive]-RNA_exp_norm[,itscm] #log ratios= fold change
A_naive_tscm=(RNA_exp_norm[,inaive]+RNA_exp_norm[,itscm])/2 

M_naive_tcm= RNA_exp_norm[,inaive]-RNA_exp_norm[,itcm] #log ratios= fold change
A_naive_tcm=(RNA_exp_norm[,inaive]+RNA_exp_norm[,itcm])/2 

M_naive_tem= RNA_exp_norm[,inaive]-RNA_exp_norm[,item] #log ratios= fold change
A_naive_tem=(RNA_exp_norm[,inaive]+RNA_exp_norm[,item])/2 

M_tscm_tcm= RNA_exp_norm[,itscm]-RNA_exp_norm[,itcm] #log ratios= fold change
A_tscm_tcm=(RNA_exp_norm[,itscm]+RNA_exp_norm[,itcm])/2 

M_tscm_tem= RNA_exp_norm[,itscm]-RNA_exp_norm[,item] #log ratios= fold change
A_tscm_tem=(RNA_exp_norm[,itscm]+RNA_exp_norm[,item])/2 

M_tcm_tem= RNA_exp_norm[,itcm]-RNA_exp_norm[,item] #log ratios= fold change
A_tcm_tem=(RNA_exp_norm[,itcm]+RNA_exp_norm[,item])/2 

  #Affy
  gset$`cell population:ch1`
  
  M_naive_tscm= (exprs(gset)[, sampleNames(gset)[gset$`cell population:ch1`=="Naive T cell"]]) - (exprs(gset)[, gset$`cell population:ch1`=="Stem cell memory T cell"])
  A_naive_tscm=((exprs(gset)[, sampleNames(gset)[gset$`cell population:ch1`=="Naive T cell"]]) + (exprs(gset)[, gset$`cell population:ch1`=="Stem cell memory T cell"]))/2 
  
  M_naive_tcm= (exprs(gset)[, sampleNames(gset)[gset$`cell population:ch1`=="Naive T cell"]]) - (exprs(gset)[, gset$`cell population:ch1`=="Central memory T cell"]) #log ratios= fold change
  A_naive_tcm=((exprs(gset)[, sampleNames(gset)[gset$`cell population:ch1`=="Naive T cell"]]) + (exprs(gset)[, gset$`cell population:ch1`=="Central memory T cell"]))/2 
  
  M_naive_tem= (exprs(gset)[, sampleNames(gset)[gset$`cell population:ch1`=="Naive T cell"]]) - (exprs(gset)[, gset$`cell population:ch1`=="Effector memory T cell"]) #log ratios= fold change
  A_naive_tem=((exprs(gset)[, sampleNames(gset)[gset$`cell population:ch1`=="Naive T cell"]]) + (exprs(gset)[, gset$`cell population:ch1`=="Effector memory T cell"]))/2 
  
  #Affy MAPlot
  plot(rowMeans(A_naive_tscm), rowMeans(M_naive_tscm), main=" ", xlab="A", ylab="M", pch=".")
  abline(h=0, col="red")
  abline(h=2, col="green")
  abline(h=-2, col="green")
  
  plot(rowMeans(A_naive_tem), rowMeans(M_naive_tem), main=" ", xlab="A", ylab="M", pch=".")
  abline(h=0, col="red")
  abline(h=2, col="green")
  abline(h=-2, col="green")
  
  plot(rowMeans(A_naive_tcm), rowMeans(M_naive_tcm), main=" ", xlab="A", ylab="M", pch=".")
  abline(h=0, col="red")
  abline(h=2, col="green")
  abline(h=-2, col="green")

  
dim(M)#relacion entre dos condiciones=fold change
dim(A)#promedio de expresi?n


plot(rowMeans(A_naive_tscm), rowMeans(M_naive_tscm), main=" ", xlab="A", ylab="M", pch=".")
abline(h=0, col="red")
abline(h=2, col="green")
abline(h=-2, col="green")

smoothScatter(rowMeans(A), rowMeans(M), main=" ", xlab="A", ylab="M", pch=".")
abline(h=0, col="red")
abline(h=2, col="green")
abline(h=-2, col="green")

library(pheatmap)
pheatmap(RNA_exp_norm)

##5. Prueba preliminar DEG
library(Matrix)
library(genefilter)

#naive 1-6 ; tscm 7-11; tcm 12-16; tem 17-20
naive_tscm=as.matrix(RNA_exp_norm[,1:11])
f = c(0,0,0,0,0,0,1,1,1,1,1)#2 condiciones
f = as.factor(f)

r2 = rowttests(naive_tscm, f)
r2
head(r2)


p.adjust(r2$p.value, method="bonf")->Pajustados 
length(Pajustados)
Pajustados

  #Affy prueba
  library(Matrix)
  library(genefilter)
  gset$`cell population:ch1`
  naive_tscm=as.matrix(RNA_exp_norm[,sampleNames(gset)[gset$`cell population:ch1`=="Naive T cell" | gset$`cell population:ch1`=="Stem cell memory T cell"]])
  naive_tcm=as.matrix(RNA_exp_norm[,sampleNames(gset)[gset$`cell population:ch1`=="Naive T cell" | gset$`cell population:ch1`=="Central memory T cell"]])
  naive_tem=as.matrix(RNA_exp_norm[,sampleNames(gset)[gset$`cell population:ch1`=="Naive T cell" | gset$`cell population:ch1`=="Effector memory T cell"]])
  
  
  f = c(0,0,0,0,0,0,1,1,1,1,1,1)#2 condiciones
  f = as.factor(f)
  
  r2 = rowttests(naive_tscm, f)
  r2
  head(r2)
  
  
  p.adjust(r2$p.value, method="bonf")->Pajustados 
  length(Pajustados)
  Pajustados

##6. Analisis SAM

  library(splines)
  #BiocManager::install("siggenes")
  library(siggenes)
  
  sam.out_1 <- sam(naive_tscm, f)
  sam.out_2 <- sam(naive_tcm, f)
  sam.out_3 <- sam(naive_tem, f)
  
  sam.out_1
  sam.out_2
  sam.out_3
  
  #Naive vs tscm
  print(sam.out_1, seq(1.6,2.4,3.2))
  plot(sam.out_1,1.6)
  plot(sam.out_1,2.4)  
  plot(sam.out_1,3.2)  
  sum.sam.out_1 <- summary(sam.out_1, 2.4)
  sum.sam.out_1 #Tabla DE
  write.csv(sum.sam.out_1@mat.sig,"result/DEG_naive_tscm_sam.csv")
  list.siggenes(sam.out_1, 2.4) 
  
  #Naive vs tcm
  print(sam.out_2, seq(1.9,3.7))
  plot(sam.out_2,1.9)
  plot(sam.out_2,3.7)
  sum.sam.out_2 <- summary(sam.out_2, 1.9)
  sum.sam.out_2 #Tabla DE
  write.csv(sum.sam.out_2@mat.sig,"result/DEG_naive_tcm_sam.csv")
  list.siggenes(sam.out_2, 1.9) 
  
  #Naive vs tem
  print(sam.out_3, seq(2.5,5.0,7.4))
  plot(sam.out_3,2.5)
  plot(sam.out_3,5.0)  
  plot(sam.out_3,7.4)
  sum.sam.out_3 <- summary(sam.out_3, 2.5)
  sum.sam.out_3 #Tabla DE
  write.csv(sum.sam.out_3@mat.sig,"result/DEG_naive_tem_sam.csv")
  list.siggenes(sam.out_1, 2.5) 
  
  #Herramienta buscar
  findDelta(sam.out, fdr = 0.04)
  findDelta(sam.out, genes = 500)  
  
  
##7. Análisis LIMA
  library("limma")
  #BiocManager::install("statmod")
  library("statmod")
  
  #Web: https://medium.com/biosyntax/using-limma-to-find-differentially-expressed-genes-1c12c905e3b1
  
  # Create a design matrix
  gset$`cell population:ch1`
  design <- model.matrix(~ 0+factor(c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4)))
  colnames(design) <- c("Naive", "Tscm", "Tcm", "Tem")

  # Contrast matrix
  cont_matrix <- makeContrasts(naivevstscm = Naive-Tscm, naivevstcm = Naive-Tcm, naivevstem = Naive-Tem, levels=design)

  #Analysis
  
  # Fit the expression matrix to a linear model
  fit <- lmFit(RNA_exp_norm, design)
  
  # Compute contrast
  
  fit_contrast <-contrasts.fit(fit, cont_matrix)
  fit_contrast_1 <- contrasts.fit(fit, cont_matrix[,1])
  fit_contrast_2 <- contrasts.fit(fit, cont_matrix[,2])
  fit_contrast_3 <- contrasts.fit(fit, cont_matrix[,3])
  fit_contrast_1 #Naive vs Tscm
  fit_contrast_2 #Naive vs Tcm
  fit_contrast_3 #Naive vs Tem
  
  # Bayes statistics of differential expression
  # *There are several options to tweak!*
  fit_contrast_1 <- eBayes(fit_contrast_1)
  fit_contrast_2 <- eBayes(fit_contrast_2)
  fit_contrast_3 <- eBayes(fit_contrast_3)

  # Generate a vocalno plot to visualize differential expression
  volcanoplot(fit_contrast_1)
  volcanoplot(fit_contrast_2)
  volcanoplot(fit_contrast_3)
  
  # Generate a list of top 100 differentially expressed genes
  top_genes_1 <- topTable(fit_contrast_1, adjust = "BH")
  top_genes_2 <- topTable(fit_contrast_2, adjust = "BH")
  top_genes_3 <- topTable(fit_contrast_3, adjust = "BH")
  
  write.csv(top_genes_1,"result/DEG_naive_tscm_limma.csv")
  write.csv(top_genes_2,"result/DEG_naive_tcm_limma.csv")
  write.csv(top_genes_3,"result/DEG_naive_tem_limma.csv")
  
  # Summary of results (number of differentially expressed genes)
  result <- decideTests(fit_contrast)
  summary(result)
  
  plotMD(fit_contrast_1)
  abline(h=0, col="red")
  abline(h=2, col="green")
  abline(h=-2, col="green")
  
  plotMD(fit_contrast_2)
  abline(h=0, col="red")
  abline(h=2, col="green")
  abline(h=-2, col="green")
  
  plotMD(fit_contrast_3)
  abline(h=0, col="red")
  abline(h=2, col="green")
  abline(h=-2, col="green")
  