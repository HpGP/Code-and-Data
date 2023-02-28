##################################
#### Script : FineSTRUCTURE   ####
##################################

# Modified script from
# https://people.maths.bris.ac.uk/~madjl/finestructure/finestructureR.html

# Note: It is required to have the FineSTRUCTURE v2 software previously installed.
# https://people.maths.bris.ac.uk/~madjl/finestructure/finestructure.html

# Read in the R functions, which also calls the needed packages
library("gplots")
source("FinestructureLibrary.R")

# Load output files from FineSTRUCTURE main command
chunkfile<-"Hp_linked_hap.chunkcounts.out"
mcmcfile<-"Hp_linked_hap_mcmc.xml"
treefile<-"Hp_linked_hap_tree.xml"

#Extract additional files from fineSTRUCTURE
mappopchunkfile<-"Hp.EMlinked.mapstate.csv" # population-by-population chunkcount file for the populations used in the MAP (i.e tree)
command <- paste("fs fs -X -Y -e X2",chunkfile,treefile,mappopchunkfile)
system(paste("fs fs -X -Y -e X2",chunkfile,treefile,mappopchunkfile))

meancoincidencefile<-"Hp.EMlinked.meancoincidence.csv" # pairwise coincidence, .i.e. proportion of MCMC files where individuals are found in the same 
command <- paste("fs fs -X -Y -e meancoincidence",chunkfile,mcmcfile,meancoincidencefile)
system(paste("fs fs -X -Y -e meancoincidence",chunkfile,mcmcfile,meancoincidencefile))

#Read in the chunkcout file
dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=1)) # read in the pairwise coincidence 

#Read in the MCMC files
mcmcxml<-xmlTreeParse(mcmcfile) # read into xml format
mcmcdata<-as.data.frame.myres(mcmcxml) # convert this into a data frame

#Read and prepare the tree files
treexml<-xmlTreeParse(treefile) # read the tree as xml format
ttree<-extractTree(treexml) # extract the tree into ape's phylo format
ttree$node.label<-NULL # remove the labels of internal nodes in order not to plot them
tdend<-myapetodend(ttree,factor=1) # convert to dendrogram format

# the MAP state
mapstate<-extractValue(treexml,"Pop") # map state as a finestructure clustering
mapstatelist<-popAsList(mapstate) # and as a list of individuals in populations
popnames<-lapply(mapstatelist,NameSummary) # population names IN A REVERSIBLE FORMAT (I.E LOSSLESS)

popnamesplot<-lapply(mapstatelist,NameMoreSummary) # a nicer summary of the populations
names(popnames)<-popnamesplot # for nicety only
names(popnamesplot)<-popnamesplot # for nicety only

popdend<-makemydend(tdend,mapstatelist) # use NameSummary to make popdend
popdend<-fixMidpointMembers(popdend) # needed for obscure dendrogram reasons

popdendclear<-makemydend(tdend,mapstatelist,"NameMoreSummary") # use NameMoreSummary to make popdend
popdendclear<-fixMidpointMembers(popdendclear) # needed for obscure dendrogram reasons

# Pairwise coincidences 
fullorder<-labels(tdend) # the order according to the tree
mcmcmatrixraw<-as.matrix(read.csv(meancoincidencefile,row.names=1)) # read in the pairwise coincidence file we created earlier
mcmcmatrix<-mcmcmatrixraw[fullorder,fullorder] 
mapstatematrix<-groupingAsMatrix(mapstatelist)[fullorder,fullorder] # map state for reference

# Coancestry matrix
datamatrix<-dataraw[fullorder,fullorder] 
datamatrix

ids_list <-read.csv (file="ids_colours_list.csv",header=FALSE,sep = ",",check.names = FALSE)
ids_list

colnames(ids_list)

ids_list %>% 
  mutate(alias = deframe(df2)[cliente])

# Save coancestry matrix with IND ids
write.csv(datamatrix, "datamatrix.csv")

###### HACER CAMBIO DE IDS DE FORMA MANUAL
#leer matrix
InputFileData <- "datamatrix_R.csv"

df <-read.table(file=InputFileData,header=TRUE,sep = ",",check.names = FALSE)
row.names(df)

#names_pop <- df$X
names_pop <- df[,1]
df[,1] <- NULL
#df$X <- NULL
rownames(df) = names_pop
df2 <- as.matrix(df)

# make some colours
some.colors<-MakeColorYRP() # these are yellow-red-purple
some.colorsEnd<-MakeColorYRP(final=c(0.2,0.2,0.2)) # as above, but with a dark grey final for capped values

tmatmax<-200 # cap the heatmap
tmpmat<-df2 
tmpmat[tmpmat>tmatmax]<-tmatmax  


#el color de los labels de col
cols <- rep('#828282', ncol(df2))

#el color de los labels de row
cols_rows <- rep('#828282', nrow(df2))


#---- Colores ----
list_color <-read.table(file="list_colour_2022.csv",header=TRUE,sep = ",")
list_color

colours_list_Pop <- rep('black', nrow(list_color))
colours_list_Pop
colours_list_Pop[list_color$POP %in% "HpAfrica2"] <- '#000000'
colours_list_Pop[list_color$POP %in% "hspAfrica1SAfrica"] <- '#D9D9D9' 
colours_list_Pop[list_color$POP %in% "hpNEAfrica"] <- '#008041' 
colours_list_Pop[list_color$POP %in% "hspAfrica1WAfrica"] <- '#7F7F7F' 
colours_list_Pop[list_color$POP %in% "hspAfrica1NAmerica"] <- '#947145' 
colours_list_Pop[list_color$POP %in% "hspAfrica1MiscAmerica"] <- '#B08958' 
colours_list_Pop[list_color$POP %in% "hspSWEuropeLATAM"] <- '#0096FF' 
colours_list_Pop[list_color$POP %in% "hspSWEuropeChile"] <- '#aacde3' 
colours_list_Pop[list_color$POP %in% "hspSWEurope"] <- '#e6a5d6' 
colours_list_Pop[list_color$POP %in% "hspEurasia"] <- '#DB5ABC' 
colours_list_Pop[list_color$POP %in% "hspNEurope"] <- '#d62fae'
colours_list_Pop[list_color$POP %in% "HpAsia2"] <- '#FFD966'
colours_list_Pop[list_color$POP %in% "hspUraI"] <- '#FAB778'
colours_list_Pop[list_color$POP %in% "hspSahul"] <- '#D5B7ED'
colours_list_Pop[list_color$POP %in% "hpNorthAsia"] <- '#EB3D4F'
colours_list_Pop[list_color$POP %in% "hspIndigenous"] <- '#7B4BA0'
colours_list_Pop[list_color$POP %in% "hspEAsia"] <- '#C93342'

#---- colours_list_Pop (DAPC) ----
list_color <-read.table(file="list_colour_2022.csv",header=TRUE,sep = ",")#(file="list_color_1020Hp.csv",header=TRUE,sep = ",")
list_color

colours_list_DAPC <- rep('black', nrow(list_color))
colours_list_DAPC
colours_list_DAPC[list_color$DAPC %in% "1"] <- '#008041' 
colours_list_DAPC[list_color$DAPC %in% "2"] <- '#7F7F7F' 
colours_list_DAPC[list_color$DAPC %in% "3"] <- '#DB5ABC' 
colours_list_DAPC[list_color$DAPC %in% "4"] <- '#000000' 
colours_list_DAPC[list_color$DAPC %in% "5"] <- '#EB3D4F' 
colours_list_DAPC[list_color$DAPC %in% "6"] <- '#D5B7ED' 

rlab=as.matrix(t(cbind(colours_list_Pop,colours_list_DAPC)))
clab=as.matrix((cbind(colours_list_Pop,colours_list_DAPC)))


pdf("Draft005.pdf",height=48,width=48)

p<-heatmap.3( tmpmat, #la funcion que se utiliznado usualmente es heatmap.2
             scale = "none", #scale data per row, column, or none
             Rowv=tdend, #nota importante tdend se obtiene corriendo primero el codigoR_Hp.r de fineSTRUCTURE
             Colv=tdend, 
             dendrogram = "both", #"both","row","column","none"
             trace="none",
             key=TRUE,
             keysize = 0.5,
             cexRow=0.3, #0.1 
             cexCol = 0.3, #0.1
             col=some.colorsEnd, #dendrogram color palette
             #symkey=FALSE,
             #xlab = "Identifier",
             #ylab = "Rows",
             #colCol = cols,
             RowSideColorsSize=2,
             ColSideColorsSize=2,
             ColSideColors=clab,#colours_list_Pop,
             RowSideColors=rlab#colours_list_Pop#colours_list_cagPAI, #colours_list_Fs#colours_list_ISMEJ # 
             #colCol = cols,
             #colRow = cols_rows,
             #ColSideColors=colours_list,
             #lmat=rbind(c(0,5,4,0,0), c(0,3,2,1,0)), 
             #lhei=c(2,5),
             #lwid=c(1,1,4,0.25,1)
             )


dev.off()


#---- Dendogram figure ----
library(dendextend)
IDS <- list_color$ID_FILES
IDS <- list_color$Strain_name
pdf(file="Hp_FullDendrogram_Colour_Version21Nov2022.pdf",height=8,width=150) #24 normal orginal
par(mar=c(6,0,2,0),mfrow=c(1,1))
tdend %>% set("labels", IDS ) %>%  set("labels_col", colours_list_Pop) %>% set("labels_cex", 1) %>% plot #Colour

tdend %>% set("labels", IDS)  %>% set("labels_cex", 0.5) %>% plot #RAW

#tdend %>% rect.dendrogram(k=24, border = 8, lty = 5, lwd = 1)
tdend %>% rect.dendrogram(k=18, border = 8, lty = 5, lwd = 1)
#tdend %>% rect.dendrogram(k=18, border = 8, lty = 5, lwd = 1)
#tdend %>% rect.dendrogram(k=13, border = 8, lty = 5, lwd = 1) ##This 

abline(h = 9, lty = 2,col="orange",lwd = 2) #18
abline(h = 10, lty = 2,col="red",lwd = 2) #13
abline(h = 11, lty = 2,col="blue",lwd = 2) #23

#the_bars =cbind(DAPC=colours_list,ISMEJ=colours_list_ISMEJ,new=colours_list_Fs,cagPAI=colours_list_cagPAI)
the_bars =cbind(#Level10_K13=colours_list_Group13,
                #subregion=colours_list_subregion, #NCBI
                #Level11_K18=colours_list_k18, #NCBI
                #Level9_K18=colours_list_Group18,
                #Level8_K23=colours_list_Group23, 
                #Level7_K24=colours_list_Group24, 
                #Level6_K31=colours_list_Group31, 
                #Previous=colours_list_Pop,
                Pop=colours_list_Pop)

colored_bars(colors = the_bars, dend = tdend)

dev.off()

#---- PCA Principal Components Analysis ----
data_dataraw <-data.frame(dataraw)#####extra
write.csv(data_dataraw, file = "data_dataraw.csv")#####extra


pcares<-mypca(dataraw)
pcares$values

saveRDS(pcares, file="pcares.RData")

# For figuring out how many PCs are important; see Lawson & Falush 2012
# You need packages GPArotation and paran
tmap<-optimalMap(dataraw)
warnings()
thorn<-optimalHorn(dataraw)
thorn
c(tmap,thorn) # 11 and 5. Horn typically underestimates, Map is usually better
pcapops<-getPopIndices(rownames(dataraw),mapstatelist)
pcapops
data_pcapops <-data.frame(pcapops)#####extra
write.csv(data_pcapops, file = "data_pcapops.csv")#####extra

pcanames<-rownames(dataraw)
pcanames

write.csv(pcanames, file = "data_pcanames.csv")#####extra

rcols<-rainbow(max(pcapops))
pcapops

data_rcols<-data.frame(rcols) #####extra
write.csv(data_rcols, file = "data_rcols.csv")#####extra


que <-data.frame(rcols[pcapops])#####extra
Colors_pop_file = read.delim("Colors.txt", header=FALSE)#####extra

Colors_pop <- Colors_pop_file$V1 #####extra

###PCANAMES
pcanames_file = read.csv("data_pcanames3.csv")#####extra
pcanames2 <- pcanames_file$x #####extra

class(pcanames)
pcanames
class(pcanames2)
pcanames2<-as.character(pcanames2)
class(pcanames2)
pcanames2

###PCANAMES_NOW
pcanames_file_NOW = read.csv("data_pcanames2_NOW.csv")#####extra
pcanames2_NOW <- pcanames_file_NOW$x #####extra

class(pcanames2_NOW)
pcanames2_NOW<-as.numeric(pcanames2_NOW)
class(pcanames2_NOW)



###RCOLS
file_rcols2 = read.csv("rcols_2020.csv")#####extra
#file_rcols2[,1] <- NULL #####extra
rcols2 <- file_rcols2$rcols #####extra

class(rcols)
rcols
class(rcols2)
rcols2
rcols2<-as.character(rcols2)
class(rcols2)

###PCAPOPS

#pcapops2 = read.csv("data_pcapops_IDS.csv") #####extra
#names <- as.character(pcapops2[,1]) #####extra
#row.names(pcapops2)<- names #####extra
#pcapops2[,1] <- NULL #####extra
#pcapops2 <- unclass(pcapops2) #####extra



pcares$vectors[,i]

pdf("PCA.pdf",height=48,width=36)
par(mfrow=c(4,3))
for(i in 1:10) for(j in (i+1):11) {
  plot(pcares$vectors[,i],
       pcares$vectors[,j],
       #col=rcols[pcapops], #original
       #col=rcols2[pcapops], #####extra
       col=Colors_pop, #####extra
       xlab=paste("PC",i),
       ylab=paste("PC",j),
       main=paste("PC",i,"vs",j),
       # pch = c(16, 17, 18)[as.numeric(Species)]
       pch=16, #orignal
  )
      # pch= pcanames2_NOW )
  #text(pcares$vectors[,i],pcares$vectors[,j],labels=pcanames,col=rcols[pcapops],cex=0.5,pos=1) #original
  text(pcares$vectors[,i],pcares$vectors[,j],labels=pcanames2,cex=0.5,pos=1) #####extra
}
dev.off()



