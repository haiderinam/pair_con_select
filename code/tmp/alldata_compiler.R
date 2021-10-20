alldata_compiler=function(alldata,nameposctrl1,nameposctrl2,namegene,mtn,nameposctrl1mt,nameposctrl2mt){
#Sorting by whether user wants to search by specific mutations or the general format overall
if(mtn=='N') {
  Positive_Ctrl1=as.numeric(!grepl("NaN|^[0]",do.call("$",list(alldata,nameposctrl1)),ignore.case = TRUE)) #Searching for 'NaN' or '0' in desired gene. A positive hit returns 0 and vise versa.
  Positive_Ctrl2=as.numeric(!grepl("NaN|^[0]",do.call("$",list(alldata,nameposctrl2)),ignore.case = TRUE))
} else {
  Positive_Ctrl1=as.numeric(grepl(nameposctrl1mt,do.call("$",list(alldata,nameposctrl1)),ignore.case = TRUE))
  Positive_Ctrl2=as.numeric(grepl(nameposctrl2mt,do.call("$",list(alldata,nameposctrl2)),ignore.case = TRUE))
}
genex=as.numeric(!grepl(paste(c("NaN","0"), collapse = "|"),do.call("$",list(alldata,namegene)))) #Searching for 'NaN' or '0' in desired gene. A positive hit returns 0 and vise versa.

alldata_comp=cbind(alldata[,c(1,2)],Positive_Ctrl1,Positive_Ctrl2,genex) #Data frame with all the data compiled. I add in random array to this later

genex_replication_prop=sum(genex)/length(genex) #Calculating Replication Proportion

#Creating Random Array (Negative Control)
rndmarray=rbinom(length(Positive_Ctrl1),1,genex_replication_prop) #change alldata            #Decided to just call random array in the for-loop
alldata_comp=cbind(alldata_comp,rndmarray) #Adding in random array to compiled gene data

# max subsample size is whichever integer is the minimum of: number of genex in sample vs. number of NOTgenex or PC1 or PC2 in sample*genexreplication proportion. For example we had 165 PC1s, therefore 340-175=165 PC1s in our data. Since we want to sample 7.35% positive hits and 92.65% negative hits, a sample containing too many positive hits can impose an upper bound on our subsampling size.
# maxsubsamplesize=min(c(genex_replication_prop*length(alldata_comp$Positive_Ctrl1),floor((length(alldata_comp$Positive_Ctrl1)-sum(alldata_comp$Positive_Ctrl1))*genex_replication_prop),floor((length(alldata_comp$Positive_Ctrl1)-sum(alldata_comp$Positive_Ctrl2))*genex_replication_prop),floor((length(alldata_comp$Positive_Ctrl1)-sum(alldata_comp$genex))*genex_replication_prop)))
# nsubsamples=maxsubsamplesize
return(list(genex_replication_prop,alldata_comp))
}