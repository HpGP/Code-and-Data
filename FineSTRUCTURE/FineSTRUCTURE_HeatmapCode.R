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

# Save coancestry matrix with IND ids
write.csv(datamatrix, "datamatrix.csv")

# Note
# Replace the INDXXX IDs with the IDs of each strain according to the list_colour_2022.csv file, 
# using an external program (such as Excel using the VLOOKUP function).

#Read the matrix with ids names
InputFileData <- "datamatrix_R.csv"
datamatrix_R <-read.table(file=InputFileData,header=TRUE,sep = ",",check.names = FALSE)

names_pop <- datamatrix_R[,1]
datamatrix_R[,1] <- NULL
rownames(datamatrix_R) = names_pop
datamatrix_R <- as.matrix(datamatrix_R)

# make some colours
some.colors<-MakeColorYRP() # these are yellow-red-purple
some.colorsEnd<-MakeColorYRP(final=c(0.2,0.2,0.2)) # as above, but with a dark grey final for capped values

tmatmax<-200 # cap the heatmap
tmpmat<-datamatrix_R 
tmpmat[tmpmat>tmatmax]<-tmatmax  

#Colours
cols <- rep('#828282', ncol(datamatrix_R))
cols_rows <- rep('#828282', nrow(datamatrix_R))
list_color <-read.table(file="list_colour.csv",header=TRUE,sep = ",")

colours_list_Pop <- rep('black', nrow(list_color))
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

colours_list_DAPC <- rep('black', nrow(list_color))
colours_list_DAPC
colours_list_DAPC[list_color$DAPC %in% "1"] <- '#008041' 
colours_list_DAPC[list_color$DAPC %in% "2"] <- '#7F7F7F' 
colours_list_DAPC[list_color$DAPC %in% "3"] <- '#DB5ABC' 
colours_list_DAPC[list_color$DAPC %in% "4"] <- '#000000' 
colours_list_DAPC[list_color$DAPC %in% "5"] <- '#EB3D4F' 
colours_list_DAPC[list_color$DAPC %in% "6"] <- '#D5B7ED' 

rlab=as.matrix(t(cbind(colours_list_Pop)))
clab=as.matrix((cbind(colours_list_Pop,colours_list_DAPC)))

#Making heatmap
pdf("Fs_coancestry_heatmap.pdf",height=48,width=48)

library(downloader)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

p<-heatmap.3( tmpmat, 
             scale = "none",
             Rowv=tdend,
             Colv=tdend, 
             dendrogram = "both", 
             trace="none",
             key=TRUE,
             keysize = 0.5,
             cexRow=0.3,  
             cexCol = 0.3, 
             col=some.colorsEnd, #dendrogram color palette
             colCol = cols,
             colRow = cols_rows,
             RowSideColorsSize=2,
             ColSideColorsSize=2,
             ColSideColors=clab,
             RowSideColors=rlab
             
             )

dev.off()
