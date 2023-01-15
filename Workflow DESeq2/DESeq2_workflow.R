## Taller DESeq2

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")


##### paquetes ####
require(DESeq2)
require(ggplot2)
require(corrplot)

library(apeglm)


##### lectura de datos ####
#setwd
setwd("~/2022/cursos/Genómica/Talleres/Taller3A-DESeq2")

Ferrets=read.delim("GSE147507_RawReadCounts_Ferret.tsv",stringsAsFactors = F) # lectura
Ferrets=Ferrets[,c(1,2:5,6:9,10:13,16:19)] # individuos de interés
Ferrets=Ferrets[rowSums(Ferrets[,-1])!=0,] # sondas con lecturas

##### renombrar genes Mustela con nombres de genes humano ####
Mustela=read.csv("Mustela.csv",stringsAsFactors = F) # diccionario
Mustela=data.frame(X=Mustela$gene_id, gene=Mustela$gene_name,stringsAsFactors = F) # nombres GS
Mustela=Mustela[Mustela$X %in% Ferrets$X,] # sondas de inter?s
Mustela=unique.data.frame(Mustela) # remover duplicados
Mustela[sample(dim(Mustela)[1],5),]

Reads=cbind(Ferrets[order(Ferrets$X),2:17],Mustela$gene[order(Mustela$X)]) # agregar GS
colnames(Reads)[17]<-"GS"
Reads=Reads[Reads$GS!="",] # eliminar lecturas sin GS
reads=aggregate(Reads[,-17], by=list(Reads$GS), max) # agregar por genes
colnames(reads)[1]<-"GS"

E=reads[,-1]#matriz de expresi?n
colnames(E)<-c('Ctl.d1.r1','Ctl.d1.r2','Covid.d1.r1','Covid.d1.r2',
               'Ctl.d3.r1','Ctl.d3.r2','Covid.d3.r1','Covid.d3.r2',
               'Ctl.d7.r1','Ctl.d7.r2','Covid.d7.r1','Covid.d7.r2',
               'Ctl.d14.r1','Ctl.d14.r2','Covid.d14.r1','Covid.d14.r2')
rownames(E)<-reads[,1]
head(E)
boxplot(E,las=2)

# filtro
modE=E[1e5>apply(E,1,max),] # eliminar genes con lecturas muy altas
modE=modE[0<apply(modE[,1:12],1,var),] # eliminar genes sin varianza

##### DESeq2 ####
coldata <- data.frame(
  time=factor(c(1,1,1,1,3,3,3,3,7,7,7,7,14,14,14,14)),
  condition=factor(rep(rep(c('h','d'),ea=2),4))
)
rownames(coldata)= colnames(E)
coldata$condition <- relevel(coldata$condition, ref = 'h')
#ref es la condicion control


## ignorando el tiempo 
#modelo para cada gen con la condici?n covid vs control
dds0 <- DESeqDataSetFromMatrix(
  countData = modE[,1:12],
  colData = coldata[1:12,],
  design= ~ condition)

########################
########NORMALIZACION 

dds0 <- estimateSizeFactors(dds0)
sizeFactors(dds0)

head(counts(dds0,normalized=TRUE))
boxplot(counts(dds0,normalized=TRUE),las=2)

dds0 <- DESeq(dds0, test="LRT", reduced=~1)#prueba de verosimilitud 
dds0 <- DESeq(dds0, test="Wald")#prueba de Wald solo para dos condiciones
res0 <- results(dds0)

resultsNames(dds0) # coeficientes
plotDispEsts(dds0) # estimaciones dispersiones

#una vez ajustado el modelo se aplica el metodo de shrinkage
resLFC0 <- lfcShrink(dds0, coef="condition_d_vs_h", type="apeglm")
#funcion del paquete "apeglm"
resLFC0 # LFC reducidos

plotMA(resLFC0)

Padj0=data.frame(p=res0$padj[res0$padj<0.95])
ggplot(Padj0, aes(x=p)) + geom_histogram(breaks = seq(0,0.95,0.05)) +
  labs(y='frecuencia') + scale_x_continuous(breaks = seq(0,0.95,0.05))

plotCounts(dds0, gene=which.min(res0$padj), intgroup="condition")

# diseño completo con factor tiempo y condicion
dds <- DESeqDataSetFromMatrix(
  countData = modE[,1:12],
  colData = coldata[1:12,],
  design= ~ time*condition)



dds <- DESeq(dds, test="LRT", reduced=~1)
res <- results(dds)
resultsNames(dds) # coeficientes

Padj=data.frame(p=res$padj[res$padj<0.95])
ggplot(Padj, aes(x=p)) + geom_histogram(breaks = seq(0,0.95,0.05)) +
  labs(y='frecuencia') + scale_x_continuous(breaks = seq(0,0.95,0.05))

plotMA(res, ylim=c(-4,4))

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

# 9 genes con menor valor p
padj=res$padj
genes=rownames(modE)[order(padj)]
par(mfrow=c(3,3), mar=c(2,2,1,1))
for(i in 1:9){
  plotCounts(dds, gene=genes[i], intgroup=c("time"),
             pch=20 ,col=c(1,1,2,2,1,1,2,2,1,1,2,2,1,2))
}

par(mfrow = c(1,1))

rld <- rlog(dds, blind=FALSE) # transformación log regularizado
boxplot(assay(rld), pch=20, las=2)

corrplot(cor(assay(rld)), cl.lim = c(0.95,1), is.corr = F, tl.col = "black",
         col = colorRampPalette(c("blue","white","black"))(200))

resOrdered <- res[order(res$padj),]
resSig <- subset(res, padj < 0.1)
resSig

write.csv(as.data.frame(resSig), 
          file="time7.conditiond_results.csv")

#Taller individual
#Realice un an?lisis de expresi?n diferencial para unos datos de RNAseq (tabla de conteos) de su selecci?n.
