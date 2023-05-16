micro_matrix=read.csv("ss_fina.csv",header = T,sep="\t",row.names = 1,comment.char = "!")
NC=micro_matrix[,c(1:3)]
DSS=micro_matrix[,c(4:6)]
SS-L=micro_matrix[,c(7:9)]
SS-M=micro_matrix[,c(10:12)]
SS-H=micro_matrix[,c(13:15)]

NC$rowMeans=apply(NC,1,function(x)mean(as.numeric(x),na.rm=T))
NC=NC[apply(NC[4],1,function(x) all(x!=0)),]
NC$rowMeans=NULL

SS-L$rowMeans=apply(SS-L,1,function(x)mean(as.numeric(x),na.rm=T))
SS-L=SS-L[apply(SS-L[4],1,function(x) all(x!=0)),]
SS-L$rowMeans=NULL

SS-M$rowMeans=apply(SS-M,1,function(x)mean(as.numeric(x),na.rm=T))
SS-M=SS-M[apply(SS-M[4],1,function(x) all(x!=0)),]
SS-M$rowMeans=NULL

SS-H$rowMeans=apply(SS-H,1,function(x)mean(as.numeric(x),na.rm=T))
SS-H=SS-H[apply(SS-H[4],1,function(x) all(x!=0)),]
SS-H$rowMeans=NULL

DSS$rowMeans=apply(DSS,1,function(x)mean(as.numeric(x),na.rm=T))
DSS=DSS[apply(DSS[4],1,function(x) all(x!=0)),]
DSS$rowMeans=NULL


library(VennDiagram)
List_ID=list("NC"=rownames(NC),
             "DSS"=rownames(DSS),
             "SS-L"=rownames(SS-L),
             "SS-M"=rownames(SS-M),
             "SS-H"=rownames(SS-H)
)

venn.plot=venn.diagram(x=List_ID,filename="Venn_up.pdf",
                       col="black",lwd=2,fontface="bold",
                       fill = c("cornflowerblue", "turquoise", "red","green","pink"))


pca=micro_matrix
pca=pca[,dim(2)]
pca$rowMeans=apply(pca,1,function(x)mean(as.numeric(x),na.rm=T))
pca=pca[order(pca$rowMeans,decreasing = T),]
pca=pca[c(1:1126),]
pca$rowMeans=NULL
library(vegan)
library(ggplot2)
library(ggrepel)
library(ape)
pca_df=as.data.frame(t(pca))
exprfl=read.csv("DEGfenlei.csv",header=T,sep = ",",row.names = 1,stringsAsFactors = F)

distance <- vegdist(pca_df,method = "bray")
df.pcoa<-pcoa(distance,correction = "none")
df.pcoa$vectors
df.pcoa$values

df.plot<-data.frame(df.pcoa$vectors)
head(df.plot)
library(ggplot2)
x_label<-round(df.pcoa$values$Relative_eig[1]*100,2)
y_label<-round(df.pcoa$values$Relative_eig[2]*100,2)
x_label
y_label
ggplot(data=df.plot,aes(x=Axis.1,y=Axis.2))+
  geom_point()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=paste0("PCoA1 ",x_label,"%"),
       y=paste0("PCoA2 ",y_label,"%"))


df.plot$group=exprfl$type
ggplot(data=df.plot,aes(x=Axis.1,y=Axis.2,
                        color=group,shape=group))+
  geom_point(size=5)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=paste0("PCoA1 ",x_label,"%"),
       y=paste0("PCoA2 ",y_label,"%"))+
  stat_ellipse(data=df.plot,
               geom = "polygon",
               aes(fill=group),
               alpha=0.3)+
  scale_fill_manual(values = c("#e31a1c","#1f78b4","#3AFEFE"))

#堆叠图
library(tidyr)
library(dplyr)
library(tibble)
phylum=read.csv("p.csv",header = T)


ggplot(phylum,aes(sample,proportion,fill = type)) + 
  geom_bar(position = "stack",stat = "identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
guides(fill=guide_legend(ncol=1))

#热图 #Family
library(pheatmap)
family=read.csv("family.csv",header = T,row.names = 1)
family = family[apply(family, 1, function(x) sd(x)!=0),] 
family_type=read.csv("family_type.csv",header = T,row.names = 1)
pheatmap(family,cluster_cols = F,cluster_rows=T,scale="row",
         show_colnames = T,show_rownames = T,
         annotation_row = family_type)

#热图 genus
genus=read.csv("genus.csv",header = T,row.names = 1)
genus = genus[apply(genus, 1, function(x) sd(x)!=0),] 
genus_type=read.csv("genus_type.csv",header = T,row.names = 1)
pheatmap(genus,cluster_cols = F,cluster_rows=T,scale="row",
         show_colnames = T,show_rownames = T,
         annotation_row = genus_type)

#lefes
BiocManager::install("microeco")
library(microeco)
library(ggplot2)
library(phyloseq)

otu_p <- as.data.frame(read.table("OTU.csv", header = TRUE, sep = ",", row.names = 1))



tax_p <- as.data.frame(read.table("TAX.csv", header = TRUE, sep = ",", row.names = 1))


sample_p=as.data.frame(read.table("sample.csv", header = TRUE, sep = ",", row.names = 1))


dataset<- microtable$new(sample_table = sample_p,
                         otu_table = otu_p, 
                         tax_table = tax_p)

lefse <- trans_diff$new(dataset = dataset, 
                        method = "lefse", 
                        group = "Group", 
                        alpha = 0.05, 
                        lefse_subgroup = NULL,
                        p_adjust_method = "none")

head(lefse$res_diff)

lefse$plot_diff_bar(use_number = 1:30, 
                    width = 0.8, 
                    group_order = c("NC", "DSS", "SS-L","SS-M","SS-H")) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg()


library(ggtree)
lefse$plot_diff_cladogram(use_taxa_num = 200, 
                          use_feature_num = 50, 
                          clade_label_level = 5, 
                          group_order = c("NC","SS-M", "DSS","SS-L","SS-H"))
