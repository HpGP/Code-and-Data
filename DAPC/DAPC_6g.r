library(adegenet)
library(ape)
seq<-read.dna("pylori_snps.fasta", format="fasta")
py<-DNAbin2genind(seq)
rm(seq)
grp<-find.clusters(py, n.pca=56, n.clust=6) #find genetic clusters
write.table(grp$grp, "kmeans_6g.txt", sep="\t", quote = FALSE)
dc <- dapc(py, grp$grp, n.pca=56, n.da=6)
DAPC_colors<-c("#008041","#7F7F7F","#DB5ABC","#000000","#EB3D4F","#d5b7ed") #colors
pdf("DAPC_hp_6g_final.pdf")
scatter(dc, scree.da=FALSE, bg="white", pch=20, cstar=0, col=DAPC_colors, cex=1, clab=.5) #DAPC plot
dev.off()
pdf("membership_prob_k6_all_final.pdf",21,7) 
compoplot(dc, col=DAPC_colors, legend=FALSE) #membership probabilities plot
dev.off()

#DAPC without outlier groups
seq<-read.dna("pylori_snps_16g_noAFSh.fasta", format="fasta")
py<-DNAbin2genind(seq)
rm(seq)
grp<-read.table("kmeans_6g_noAfSh.txt", sep="\t", h=TRUE)
dc <- dapc(py, grp$x, n.pca=56, n.da=6)
DAPC_colors<-c("#008041","#7F7F7F","#DB5ABC","#EB3D4F")
pdf("DAPC_hp_6g_noAfSh_final.pdf")
scatter(dc, scree.da=FALSE, bg="white", pch=20, cstar=0, col=DAPC_colors, cex=1, clab=.5)
dev.off()