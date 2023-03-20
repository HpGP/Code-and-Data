library(adegenet)
library(ape)
seq<-read.dna("pylori_snps.fasta", format="fasta")
py<-DNAbin2genind(seq)
rm(seq)
grp<-find.clusters(py, n.pca=56, n.clust=6) #find genetic clusters
write.table(grp$grp, "kmeans_17g.txt", sep="\t", quote = FALSE) #find genetic clusters
dc <- dapc(py, grp$grp, n.pca=56, n.da=6)
DAPC_colors<-c("#EB3D4F","#C93342","#7F7F7F","#000000","#8B8000","#D62FAE","#AACDE3","#DB5ABC","#0096FF","#FFD966","#556B2F ","#E6A5D6","#4169e1","#7B4BA0","#008041","#FAB778","#d5b7ed") #colors
pdf("DAPC_hp_17g_final.pdf")
scatter(dc, scree.da=FALSE, bg="white", pch=20, cstar=0, col=DAPC_colors, cex=1, clab=.5) #DAPC plot
dev.off()
pdf("membership_prob_k17_all_final.pdf",21,7) #plot membership probabilities
compoplot(dc, col=DAPC_colors, legend=FALSE) 
dev.off() 

#DAPC without outlier groups
seq<-read.dna("pylori_snps_17g_noAFSh.fasta", format="fasta")
py<-DNAbin2genind(seq)
rm(seq)
grp<-read.table("1kmeans_17k_noAfSh.txt", sep="\t", h=TRUE)
dc <- dapc(py, grp$x, n.pca=56, n.da=6)
DAPC_colors<-c("#EB3D4F","#C93342","#7F7F7F","#8B8000","#D62FAE","#AACDE3","#DB5ABC","#0096FF","#FFD966","#556B2F ","#E6A5D6","#4169e1","#7B4BA0","#008041","#FAB778")
pdf("DAPC_hp_17g_noAfSh_final.pdf")
scatter(dc, scree.da=FALSE, bg="white", pch=20, cstar=0, col=DAPC_colors, cex=1, clab=.5)
dev.off()
