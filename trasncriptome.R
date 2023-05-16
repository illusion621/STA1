library(DESeq2)
library(tidyverse)
library(GOplot)
library(clusterProfiler)
library(ggplot2)
library(AnnotationDbi)
library(AnnotationHub)
library(org.Mm.eg.db)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(ggThemeAssist)
library(ggplot2)
library(ggrepel)
library(ggthemes)
#coutns转换
a=read.csv("shuisu.csv",header=T,row.names=1)
counts=a
counts$gene=rownames(counts) #设置新一列gene，然后这列等于counts的行名
counts=separate(data = counts, col = gene, into = c("ensemble", "SYMBOL","function1"), sep = "-",convert = T) #把counts数据的gene这一列按照”-“符号劈开，然后分为三列，分别是ensemble，symbol，function1
counts$ensemble=NULL #去掉esemble这一列
counts$function1=NULL #去掉function1这一列
ID=bitr(counts$SYMBOL,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = "org.Mm.eg.db") #设置一个新的数据框，将counts数据symbol这一列拿出来从symbolgene 向 ensemble 基因名转换
counts=inner_join(counts,ID,by="SYMBOL") #将转换好的ID，包含symbol和ensemble两列，与counts数据根据symbol名合并
counts$SYMBOL=NULL #消去symbol一列
counts=counts[!duplicated(counts[,dim(counts)[2]]),] #根据最后一列，也就是counts列消去重复值
rownames(counts)=counts$ENSEMBL #设置行名为ensemble一列
counts$ENSEMBL=NULL #消去ensemble一列


eff_length=read.csv("gene_length2.csv",header = T) #读取小鼠基因长度文件
rownames(eff_length)=eff_length$gene_id #格式，行名
colnames(eff_length)=c("gene_id","eff_length") #格式 列名
rownames(eff_length)=do.call(rbind,strsplit(as.character(eff_length$gene_id),"\\."))[,1] #去掉genelength的小数点
feature_ids=rownames(counts) #将ID拿出来
counts=counts[feature_ids %in% rownames(eff_length),] #判断样本文件基因ID是否都在基因注释文件里
mm=match(rownames(counts),rownames(eff_length)) #将两文件基因名位置匹配
eff_length=eff_length[mm,] #在小鼠基因长度文件中筛选出有counts基因的数据

x=counts / eff_length$eff_length #tpm算法
counts_tpm=t(t(x) / colSums(x) ) * 1e6 #tpm算法
colSums(counts_tpm) #看是否标准化为tpm 值都为1e6

counts_tpm=as.data.frame(counts_tpm) #将tpm文件转化为矩阵
counts_tpm$ENSEMBL=rownames(counts_tpm) #ensemble为行名

DI=bitr(counts_tpm$ENSEMBL,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Mm.eg.db") #跟上面一样，这次是把ensemble转化为symbol ID
counts_tpm=inner_join(counts_tpm,DI,by="ENSEMBL")
counts_tpm$ENSEMBL=NULL
counts_tpm=counts_tpm[!duplicated(counts_tpm[,dim(counts_tpm)[2]]),]
rownames(counts_tpm)=counts_tpm$SYMBOL
counts_tpm$SYMBOL=NULL

#差异分析
color = colorRampPalette(rev(brewer.pal(n=11, name='RdYlGn')))(round(100))
a3=a[,c(1,2,3,4,11,12,13)]
group3=read.csv("group3.csv",stringsAsFactors = T)
dds3=DESeqDataSetFromMatrix(countData = a3,
                            colData = group3,
                            design = ~dex)
dds3=DESeq(dds3)
res3=results(dds3)
res3=data.frame(res3)

res3$gene=rownames(res3)
res3=separate(data = res3, col = gene, into = c("ensemble", "SYMBOL","function1"), sep = "-")
res3$function1=NULL
res3=res3[!duplicated(res3[,dim(res3)[2]]),]
rownames(res3)=res3$SYMBOL
res3=res3[c(rownames(coldss)),]
res3=res3[abs(res3$log2FoldChange)>0.5,]
p3=counts_tpm[,c(1,2,3,4,11,12,13)]
p3=p3[rownames(res3),]
pheatmap(p3,cluster_cols = F,cluster_rows = F,scale="row",
         color = colorRampPalette(rev(brewer.pal(n=11, name='RdYlGn')))(round(1000)),
         main = "DEGheatmap",annotation_col=rowdss,annotation_colors = annocolors,annotation_row = coldss)

#第四组 C1 C2 C3 VS G1 G2 G3
a4=a[,c(5,6,7,11,12,13)]
group4=read.csv("group4.csv",stringsAsFactors = T)
dds4=DESeqDataSetFromMatrix(countData = a4,
                            colData = group4,
                            design = ~dex)
dds4=DESeq(dds4)
res4=results(dds4)
res4=data.frame(res4)

res4$gene=rownames(res4)
res4=separate(data = res4, col = gene, into = c("ensemble", "SYMBOL","function1"), sep = "-")
res4$function1=NULL
res4=res4[!duplicated(res4[,dim(res4)[2]]),]
rownames(res4)=res4$SYMBOL
res4=res4[c(rownames(coldss)),]
res4=res4[abs(res4$log2FoldChange)>0.5,]
p4=counts_tpm[,c(5,6,7,11,12,13)]
p4=p4[rownames(res4),]


#水苏DSS PPAR热图
col511=read.csv("PPAR.csv",header = T,row.names = 1) 
row511=read.csv("511fenzu.csv",header = T,row.names = 1)
PPAR=counts_tpm[rownames(counts_tpm)%in%rownames(col511),]
PPAR=PPAR[,colnames(PPAR)%in%rownames(row511)]
PPAR=PPAR[-16,]
pheatmap(PPAR,cluster_cols = F,cluster_rows = F,scale="row",
         color = colorRampPalette(rev(brewer.pal(n=11, name='RdYlGn')))(round(1000)),
         main = "PPARheatmap",annotation_col=row511,annotation_colors = annocolors,annotation_row = col511)


#水苏DSS TH17/TREG热图
col511=read.csv("th17treg.csv",header = T,row.names = 1) 
row511=read.csv("511fenzu.csv",header = T,row.names = 1)
TT=counts_tpm[rownames(counts_tpm)%in%rownames(col511),]
TT=TT,colnames(TT)%in%rownames(row511)]
pheatmap(TT,cluster_cols = F,cluster_rows = F,scale="row",
         color = colorRampPalette(rev(brewer.pal(n=11, name='RdYlGn')))(round(1000)),
         main = "TH17/TREGheatmap",annotation_col=row511,annotation_colors = annocolors,annotation_row = col511)




