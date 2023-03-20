library(ape) 
library(adegenet) 
library(ggplot2) 

pca=1100 #number of PCs 
rep=100 #number of replicates 
nc=20 #max number of clusters 

seq<-read.dna("HpGP_snps.fasta", format="fasta")#input fasta file
gen<-DNAbin2genind(seq) #convert to genind
for(i in 1:rep){ 
  grp <- find.clusters(gen, max.n.clust = nc, n.pca = pca, choose.n.clust = FALSE) #iterate find.clusters
  if(i==1){ 
    k<-grp$Kstat 
  } 
  else{ 
    k<-rbind(k, grp$Kstat) 
  } 
} 
p<-as.data.frame(k) 
colnames(p)<- c(1:nc) 
rownames(p)<-c() 
write.table(p, "kmeans_all.txt", sep="\t", quote=FALSE) 
m<-c() 
se<-c() 
for(i in 1:nc){ 
  v<-p[,i] 
  m<-append(m, mean(as.numeric(v))) #mean value for each k
  se<-append(se, 3*sd(v)) #standard error for each k
   
} 
k_means<-data.frame(m, se) #convert in dataframe
pdf("pylori_kmeans_100rep.pdf") #plot in pdf
ggplot(k_means, aes(x=1:nc,y=m)) + geom_line() +  geom_errorbar(width=.1, aes(ymin=m-se, ymax=m+se)) + geom_point() + ggtitle("Value of BIC versus number of clusters") + ylab("Value of BIC") + xlab("Number of Clusters") + theme_minimal() 
dev.off()