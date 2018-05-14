library(psych)
library(DoubletDecon)

######################
#                    #
# Correlation Test 1 #
#                    #
# Test with 10% of   #
# each group removed #
#                    #
######################

#Grimes without doublet removal
# expt_name="GrimesICGS"
# location="/Users/car6wn/Documents/Projects/cellHarmony/"
# rawDataFile=paste0(location,"exp.Grimes_HaromizeReference_filt.txt")
# expressionFile=read.table(rawDataFile, sep="\t", header=T, row.names=1)

#Grimes MF without doublet removal
expt_name="GrimesMF"
location="/Users/car6wn/Documents/Projects/cellHarmony/"
rawDataFile=paste0(location,"Mm-Grimes-Nature-MarkerFinder-cellHarmony-reference-v2.txt")
expressionFile=read.table(rawDataFile, sep="\t", header=T, row.names=1)

#Melanoma without doublet removal
expt_name="MelanomaICGS"
location="/Users/car6wn/Documents/Projects/cellHarmony/"
rawDataFile=paste0(location,"Human-Melanoma-Guide3-cellHarmony-reference.txt")
expressionFile=read.table(rawDataFile, sep="\t", header=T, row.names=1)

#Melanoma MF without doublet removal
expt_name="MelanomaMF"
location="/Users/car6wn/Documents/Projects/cellHarmony/"
rawDataFile=paste0(location,"Melanoma-MarkerFinder-cellHarmony-reference.txt")
expressionFile=read.table(rawDataFile, sep="\t", header=T, row.names=1)

#Melanoma with doublet removal
expt_name="MelanomaICGS_DD"
location="/Users/car6wn/Documents/Projects/cellHarmony/"
rawDataFile=paste0(location,"Human-Melanoma-Guide3-cellHarmony-reference.txt")
groupsFile=paste0(location,"groups.Human-Melanoma-Guide3-cellHarmony-reference.txt")
expressionFile=read.table(rawDataFile, sep="\t", header=T, row.names=1)
removeCC=F
species="mmu"
write=F
recluster="doublets_decon"
isADoublet=1
filename="DD_and_CH_test"
corrCutoff=1.15
DD=Main_Doublet_Decon(rawDataFile, groupsFile, filename, removeCC, species, corrCutoff, write, recluster, isADoublet)
expressionFile=expressionFile[,colnames(expressionFile) %in% row.names(DD$Final_nondoublets_groups)]

#Melanoma with doublet removal MF
expt_name="MelanomaMF_DD"
location="/Users/car6wn/Documents/Projects/cellHarmony/"
rawDataFile=paste0(location,"Melanoma-MarkerFinder-cellHarmony-reference.txt")
groupsFile=paste0(location,"groups.Human-Melanoma-Guide3-cellHarmony-reference.txt")
expressionFile=read.table(rawDataFile, sep="\t", header=T, row.names=1)
removeCC=F
species="mmu"
write=F
recluster="doublets_decon"
isADoublet=1
filename="DD_and_CH_test"
corrCutoff=1.15
DD=Main_Doublet_Decon(rawDataFile, groupsFile, filename, removeCC, species, corrCutoff, write, recluster, isADoublet)
expressionFile=expressionFile[,colnames(expressionFile) %in% row.names(DD$Final_nondoublets_groups)]

#Grimes, new groups
expt_name="GrimesNew"
location="/Users/car6wn/Documents/Projects/cellHarmony/BoneMarrow/"
rawDataFile=paste0(location,"exp.Grimes_HaromizeReference_filt_NEW.txt")
expressionFile=read.table(rawDataFile, sep="\t", header=T, row.names=1)


#Testing
ntests=100
ncor=length(seq(0,1,0.05))

fullCorTable=cor(expressionFile)
groups=unique(as.numeric(expressionFile[1,2:ncol(expressionFile)]))
#for(grp in 1:length(groups)){

#print(grp)
kappa_table=as.data.frame(matrix(nrow=ncor, ncol=ntests))
placement_table=as.data.frame(matrix(nrow=ncor, ncol=ntests))
newExpressionFile=expressionFile[,2:ncol(expressionFile)] #removes the cluster

for(set in 1:ntests){
  
  print(set)
  files=Create_Testing_Set(newExpressionFile)
  #a=files$testingExpression
  #b=expressionFile[,which(expressionFile[1,]==grp)]
  #testingExpression=cbind(files$testingExpression, expressionFile[,which(expressionFile[1,]==grp)])
  testingExpression=files$testingExpression
  trainingExpression=files$trainingExpression
  
  sapply(1:ncor, wrapper_CH2, trainingExpression, testingExpression, fullCorTable)
}
row.names(kappa_table)=c(seq(0,1,0.05))
row.names(placement_table)=c(seq(0,1,0.05))
#assign(paste0(expt_name, "_kappa_table_clus_", grp), kappa_table)
assign(paste0(expt_name, "_1_kappa_table_clus"), kappa_table)
assign(paste0(expt_name, "_1_placement_table_clus"), placement_table)
#}

MelanomaMF_DD_1_kappa_table_clus[is.na(MelanomaMF_DD_1_kappa_table_clus)]<-1
full_table=(MelanomaMF_DD_1_kappa_table_clus+MelanomaMF_DD_1_placement_table_clus)/2
par(mfrow=c(2,3))

plot(MelanomaMF_DD_1_kappa_table_clus[,1], ylim=c(0,1), type="l", xaxt = "n", main="Kappa Table", ylab="Cohens K", xlab="correlation cutoff")
axis(1, at=1:ncor, labels=seq(0,1,0.05))
apply(MelanomaMF_DD_1_kappa_table_clus, 2, points, type="l")
plot(MelanomaMF_DD_1_placement_table_clus[,1], ylim=c(0,1), type="l", xaxt = "n", main="Placement Table", ylab="% Assigned", xlab="correlation cutoff")
axis(1, at=1:ncor, labels=seq(0,1,0.05))
apply(MelanomaMF_DD_1_placement_table_clus, 2, points, type="l")
plot(full_table[,1], ylim=c(0,1), type="l", xaxt = "n", main="Combined Table", ylab="Error Value", xlab="correlation cutoff")
axis(1, at=1:ncor, labels=seq(0,1,0.05))
apply(full_table, 2, points, type="l")

mean_kappa=apply(MelanomaMF_DD_1_kappa_table_clus, 1, mean)
sd_kappa=apply(MelanomaMF_DD_1_kappa_table_clus, 1, sd)
plot(mean_kappa, ylim=c(0,1), type="l", xaxt = "n", main="Kappa Table", ylab="Cohens K", xlab="correlation cutoff")
axis(1, at=1:ncor, labels=seq(0,1,0.05))
arrows(1:21, mean_kappa-sd_kappa, 1:21, mean_kappa+sd_kappa, length=0.02, angle=90, code=3)
mean_placement=apply(MelanomaMF_DD_1_placement_table_clus, 1, mean)
sd_placement=apply(MelanomaMF_DD_1_placement_table_clus, 1, sd)
plot(mean_placement, ylim=c(0,1), type="l", xaxt = "n", main="Placement Table", ylab="% Assigned", xlab="correlation cutoff")
axis(1, at=1:ncor, labels=seq(0,1,0.05))
arrows(1:21, mean_placement-sd_placement, 1:21, mean_placement+sd_placement, length=0.02, angle=90, code=3)
mean_full=apply(full_table, 1, mean)
sd_full=apply(full_table, 1, sd)
plot(mean_full, ylim=c(0,1), type="l", xaxt = "n", main="Combined Table", ylab="Error Value", xlab="correlation cutoff")
axis(1, at=1:ncor, labels=seq(0,1,0.05))
arrows(1:21, mean_full-sd_full, 1:21, mean_full+sd_full, length=0.02, angle=90, code=3)





# par(mfrow=c(ceiling(length(groups)/3),3))
# for(index in 1:length(groups)){
#   plot(get(paste0(expt_name, "_kappa_table_clus_", index))[,1], ylim=c(0,1), type="l", xaxt = "n", main=paste0("Cluster ", index), ylab="Cohens K", xlab="SDs below the max correlation")
#   axis(1, at=1:ncor, labels=seq(0,2,0.1))
#   apply(get(paste0(expt_name, "_kappa_table_clus_", index)), 2, points, type="l")
#   
# }

######################
#                    #
# Correlation Test 2 #
#                    #
# Test with 10% of   #
# each group removed #
# and one full group #
#                    #
######################

#Grimes without doublet removal
expt_name="GrimesICGS"
location="/Users/car6wn/Documents/Projects/cellHarmony/"
rawDataFile=paste0(location,"exp.Grimes_HaromizeReference_filt.txt")
expressionFile=read.table(rawDataFile, sep="\t", header=T, row.names=1)

#Grimes MF without doublet removal
expt_name="GrimesMF"
location="/Users/car6wn/Documents/Projects/cellHarmony/"
rawDataFile=paste0(location,"Mm-Grimes-Nature-MarkerFinder-cellHarmony-reference-v2.txt")
expressionFile=read.table(rawDataFile, sep="\t", header=T, row.names=1)

#Melanoma without doublet removal
expt_name="MelanomaICGS"
location="/Users/car6wn/Documents/Projects/cellHarmony/"
rawDataFile=paste0(location,"Human-Melanoma-Guide3-cellHarmony-reference.txt")
expressionFile=read.table(rawDataFile, sep="\t", header=T, row.names=1)

#Melanoma MF without doublet removal
expt_name="MelanomaMF"
location="/Users/car6wn/Documents/Projects/cellHarmony/"
rawDataFile=paste0(location,"Melanoma-MarkerFinder-cellHarmony-reference.txt")
expressionFile=read.table(rawDataFile, sep="\t", header=T, row.names=1)

#Melanoma with doublet removal
expt_name="MelanomaICGS_DD"
location="/Users/car6wn/Documents/Projects/cellHarmony/"
rawDataFile=paste0(location,"Human-Melanoma-Guide3-cellHarmony-reference.txt")
groupsFile=paste0(location,"groups.Human-Melanoma-Guide3-cellHarmony-reference.txt")
expressionFile=read.table(rawDataFile, sep="\t", header=T, row.names=1)
removeCC=F
species="mmu"
write=F
recluster="doublets_decon"
isADoublet=1
filename="DD_and_CH_test"
corrCutoff=1.15
DD=Main_Doublet_Decon(rawDataFile, groupsFile, filename, removeCC, species, corrCutoff, write, recluster, isADoublet)
expressionFile=expressionFile[,colnames(expressionFile) %in% row.names(DD$Final_nondoublets_groups)]

expt_name="MelanomaMF_DD"
location="/Users/car6wn/Documents/Projects/cellHarmony/"
rawDataFile=paste0(location,"Melanoma-MarkerFinder-cellHarmony-reference.txt")
groupsFile=paste0(location,"groups.Human-Melanoma-Guide3-cellHarmony-reference.txt")
expressionFile=read.table(rawDataFile, sep="\t", header=T, row.names=1)
removeCC=F
species="mmu"
write=F
recluster="doublets_decon"
isADoublet=1
filename="DD_and_CH_test"
corrCutoff=1.15
DD=Main_Doublet_Decon(rawDataFile, groupsFile, filename, removeCC, species, corrCutoff, write, recluster, isADoublet)
expressionFile=expressionFile[,colnames(expressionFile) %in% row.names(DD$Final_nondoublets_groups)]

#Grimes, new groups
expt_name="GrimesNew"
location="/Users/car6wn/Documents/Projects/cellHarmony/BoneMarrow/"
rawDataFile=paste0(location,"exp.Grimes_HaromizeReference_filt_NEW.txt")
expressionFile=read.table(rawDataFile, sep="\t", header=T, row.names=1)

#Testing
ntests=10
ncor=length(c(seq(0,1,0.05)))

fullCorTable=cor(expressionFile)
groups=unique(as.numeric(expressionFile[1,2:ncol(expressionFile)]))
for(grp in 1:length(groups)){
  
  print(grp)
  kappa_table=as.data.frame(matrix(nrow=ncor, ncol=ntests))
  newExpressionFile=expressionFile[,which(expressionFile[1,]!=grp)] #removes the cluster
  
  for(set in 1:ntests){
    
    print(set)
    files=Create_Testing_Set(newExpressionFile)
    a=files$testingExpression
    b=expressionFile[,which(expressionFile[1,]==grp)]
    testingExpression=cbind(files$testingExpression, expressionFile[,which(expressionFile[1,]==grp)])
    trainingExpression=files$trainingExpression
    
    sapply(1:ncor, wrapper_CH, trainingExpression, testingExpression, fullCorTable)
  }
  row.names(kappa_table)=c(seq(0,1,0.05))
  #row.names(kappa_table)=0.7
  assign(paste0(expt_name, "_kappa_table_clus_", grp), kappa_table)
}

par(mfrow=c(ceiling(length(groups)/3),3))
maxes=NULL
for(index in 1:length(groups)){
  mean_full=apply(get(paste0(expt_name, "_kappa_table_clus_", index)), 1, mean)
  maxes=c(maxes, max(mean_full))
  plot(get(paste0(expt_name, "_kappa_table_clus_", index))[,1], ylim=c(0,1), type="l", xaxt = "n", main=paste0("Cluster ", index), ylab="Cohens K", xlab="SDs below the max correlation")
  axis(1, at=1:ncor, labels=seq(0,1,0.05))
  apply(get(paste0(expt_name, "_kappa_table_clus_", index)), 2, points, type="l")
}

par(mfrow=c(2,3))
expt_name=c("GrimesICGS", "GrimesMF", "MelanomaICGS", "MelanomaMF", "MelanomaICGS_DD", "MelanomaMF_DD")
means2=NULL
for(name in 1:length(expt_name)){
  temp=NULL
  for(index in 1:length(groups)){
    apply(get(paste0(expt_name[name], "_kappa_table_clus_", index)), 1, mean)
    temp=cbind(temp, apply(get(paste0(expt_name[name], "_kappa_table_clus_", index)), 1, mean))
  }
  means=apply(temp, 1, mean)
  plot(means, ylim=c(0,1), type="l", xaxt="n", main=paste0("All Cluster Average - ", expt_name[name]), ylab="Cohens K", xlab="Correlation Cutoff")
  axis(1, at=1:ncor, labels=seq(0,1,0.05))
  means2=c(means2,max(means))
}
par(mfrow=c(1,1))
boxplot(means2, main="Average Cohens K", ylim=c(0,1))
# par(mfrow=c(1,1))
# kitty=c(unlist(get(paste0(expt_name, "_kappa_table_clus_", index))[1:9]))
# plot(kitty, ylim=c(0,1))

#####################
#                   #
#     Functions     #
#                   #
#####################

#For one cluster out
wrapper_CH <- function(index, trainingExpression, testingExpression, fullCorTable){
  
  cors=c(seq(0,1,0.05))
  #print("Start CH")
  files2=CellHarmony2(trainingExpression, testingExpression, cors[index], fullCorTable)
  #print("End CH")
  #test for errors across correlation cutoffs
  #a=files2$infoTable[which(files2$infoTable[,5]!=grp),5:6] #remove any cluster 9 cells
  a=files2$infoTable[,5:6] 
  b=cbind(t(files2$outlierExpression[1,]), rep(grp, ncol(files2$outlierExpression)))
  colnames(b)=colnames(a)
  c=rbind(a,b)
  d=cohen.kappa(c)
  kappa_table[index,set]<<-d$kappa
}

#For all clusters used
wrapper_CH2 <- function(index, trainingExpression, testingExpression, fullCorTable){
  
  cors=c(seq(0,1,0.05))
  #print("Start CH")
  files2=CellHarmony2(trainingExpression, testingExpression, cors[index], fullCorTable)
  #print("End CH")
  #test for errors across correlation cutoffs
  #a=files2$infoTable[which(files2$infoTable[,5]!=grp),5:6] #remove any cluster 9 cells
  a=files2$infoTable[,5:6]
  b=cbind(t(files2$outlierExpression[1,]), rep(0, ncol(files2$outlierExpression)))
  colnames(b)=colnames(a)
  c=rbind(a,b)
  #d=cohen.kappa(c)
  #print(paste0(set,":", index))
  if(length(which(is.na(a)))/2!=ncol(testingExpression)){
    d=cohen.kappa(a)
    kappa_table[index,set]<<-d$kappa
  }
  e=1-(ncol(files2$outlierExpression)/ncol(testingExpression))
  placement_table[index,set]<<-e
}

#Inputs: expressionFile: ICGS expression with cell groups as the first row and gene groups as the first column
#Outputs: trainingExpression: Expression file containing training of cells, cell groups as the first row and gene groups as the first column
#         testingExpression: Expression file containing testing cells, gene groups as the first column (no cell groups)
Create_Testing_Set <- function(expressionFile){
  
  #pull out 30% of cells from each cluster as testing set
  cellsToPull=NULL
  groups=unique(as.numeric(expressionFile[1,2:ncol(expressionFile)]))
  for(grp in 1:length(groups)){
    temp=expressionFile[,which(expressionFile[1,]==grp)]
    cellsToPull=c(cellsToPull,sample(colnames(temp),ceiling((ncol(temp)/30)),replace=FALSE)) #need to round up
  }
  testingExpression=expressionFile[,which(colnames(expressionFile) %in% cellsToPull)]
  trainingExpression=expressionFile[,which(!(colnames(expressionFile) %in% cellsToPull))]
  
  return(list(trainingExpression=trainingExpression, testingExpression=testingExpression))
  
}

#Inputs: trainingExpression: Expression file containing training of cells, cell groups as the first row and gene groups as the first column
#        testingExpression: Expression file containing testing cells, gene groups as the first column (no cell groups)
#        multi: Multiplier for testing various correlation cutoffs for minimum correlation to training cell to be matched
#Outputs: outlierExpression: Expression file containing cells that do not meet the minimum correlation for classification
#         trainingExpression: New expression file with testing cells matched to training cells, cell groups as the first row and gene groups as the first column
CellHarmony2 <- function(trainingExpression, testingExpression, multi, fullCorTable){
  
  #create data frame for cells that don't classify
  outlierExpression=as.data.frame(matrix(ncol=0, nrow=nrow(testingExpression)))
  infoTable=as.data.frame(matrix(ncol=6, nrow=ncol(testingExpression)))
  
  #for each cell pulled out, find which cell it is most correlated to.
  #correlationTable=cor(trainingExpression[2:nrow(trainingExpression),], testingExpression[2:nrow(testingExpression),])
  #correlationTable3=cor(trainingExpression, testingExpression)
  #correlationTable=fullCorTable[row.names(fullCorTable) %in% colnames(testingExpression), colnames(fullCorTable) %in% colnames(trainingExpression)]
  
  correlationTable=fullCorTable[colnames(fullCorTable) %in% colnames(trainingExpression), row.names(fullCorTable) %in% colnames(testingExpression)]
  maxCorrelation=max(correlationTable)
  temp2=sort(c(correlationTable))
  temp3=temp2[floor(length(temp2)*0.025):(length(temp2)-floor(length(temp2)*0.025))] #cutoff the top and bottom 2.5% of correlations to remove outliers for sd
  temp4=sd(temp3)
  #cutoff=maxCorrelation-(multi*temp4)
  cutoff=multi
  
  #if max correlation is below a threshold (test this range), save it to a different table
  #otherwise, place it in the position to the right of the cell in a new expression file
  #do this for all of the cells
  for(cell in 1:ncol(testingExpression)){
    cell_name=colnames(testingExpression)[cell]
    #print(paste0(cell, " of ", ncol(testingExpression)))
    correct_column=correlationTable[,which(colnames(correlationTable)==colnames(testingExpression)[cell])]
    localMaxCorrelation=max(correct_column)
    infoTable[cell,1]=localMaxCorrelation
    if(localMaxCorrelation >= cutoff){
      #These are the cells that we will put back in
      if(length(which(correct_column==localMaxCorrelation))>1){
        #if there is more than one cell that is most correlated
        closestCellName=row.names(correlationTable)[as.numeric(which(correct_column==localMaxCorrelation))]+1
        infoTable[cell,2]=closestCellName
      }else{
        #closestCellCol=which(correlationTable[,cell]==localMaxCorrelation)+1
        closestCellName=row.names(correlationTable)[as.numeric(which(correct_column==localMaxCorrelation))]
        infoTable[cell,2]=closestCellName
      }
      closestCellCol=which(colnames(trainingExpression)==closestCellName)
      infoTable[cell,3]=closestCellCol
      newLocation=which(colnames(trainingExpression)==closestCellName)
      infoTable[cell,4]=newLocation
      temp=testingExpression[,cell]
      infoTable[cell,5]=temp[1]
      temp[1]=trainingExpression[,newLocation][1]
      infoTable[cell,6]=temp[1]
      temp2=as.data.frame(temp)
      colnames(temp2)=colnames(testingExpression)[cell]
      if(newLocation==ncol(trainingExpression)){ #if you are attaching to the last cell
        trainingExpression=cbind(trainingExpression, temp2)
      }else{
        trainingExpression=cbind(trainingExpression[1:newLocation],temp2, trainingExpression[(newLocation+1):ncol(trainingExpression)]) #TODO: fix colnames
      }
    }else{
      #These are the cells we spit out as not in the original
      outlierExpression=cbind(outlierExpression, testingExpression[,cell])
    }
  }
  colnames(infoTable)=c("max_correaltion", "closestCellCol", "closesCellName", "newLocation", "old_group", "new_group")
  row.names(infoTable)=colnames(testingExpression)
  
  return(list(outlierExpression=outlierExpression, trainingExpression=trainingExpression, infoTable=infoTable))
}
