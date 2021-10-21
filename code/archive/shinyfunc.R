shinyfunc=function(cohort_size,goi_abundance,pc1_abundance,pc2_abundance,pc1_goi,pc1_pc2){
  library(plotly)
  library(knitr)
  library(tictoc)
  library(workflowr)
  library(VennDiagram)
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  library(devtools)
  library(ggsignif)

  source("code/contab_maker.R")
  source("code/alldata_compiler.R")
  source("code/mut_excl_genes_generator.R")
  source("code/mut_excl_genes_datapoints.R")
  source("code/simresults_generator.R")

  cleanup=theme_bw() +
    theme(plot.title = element_text(hjust=.5),
          panel.grid.major = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"))

  ###Shiny Function
  # cohort_size=500
  # goi_abundance=25
  # pc1_abundance=100
  # pc2_abundance=100
  # pc1_goi=20
  # pc1_pc2=2
  # # rbind(c(cohort_size-pc1_abundance-goi_abundance+pc1_goi,pc1_abundance-pc1_goi),c(goi_abundance-pc1_goi,pc1_goi))
  # fisher.test(rbind(c(cohort_size-pc1_abundance-goi_abundance+pc1_goi,pc1_abundance-pc1_goi),c(goi_abundance-pc1_goi,pc1_goi)))

  fet_pc1_goi=fisher.test(rbind(c(pc1_goi,pc1_abundance-pc1_goi),c(goi_abundance-pc1_goi,cohort_size-pc1_abundance-goi_abundance+pc1_goi)))
  fet_pc1_pc2=fisher.test(rbind(c(pc1_pc2,pc1_abundance-pc1_pc2),c(pc2_abundance-pc1_pc2,cohort_size-pc1_abundance-pc2_abundance+pc1_pc2)))
  gene_pair_1=unlist(mut_excl_genes_generator(cohort_size,goi_abundance,fet_pc1_goi$estimate,fet_pc1_pc2$estimate)[1])
  gene_pair_2=unlist(mut_excl_genes_generator(cohort_size,goi_abundance,fet_pc1_goi$estimate,fet_pc1_pc2$estimate)[2])
  alldata_1=mut_excl_genes_datapoints(gene_pair_1)
  alldata_2=mut_excl_genes_datapoints(gene_pair_2)
      # contab_maker(alldata_1$gene1,alldata_1$gene2,alldata_1)
      # contab_maker(alldata_2$gene1,alldata_2$gene2,alldata_2)
  alldata_comp_1=alldata_compiler(alldata_1,"gene2","gene3","gene1",'N',"N/A","N/A")[[2]]
  genex_replication_prop_1=alldata_compiler(alldata_1,"gene2","gene3","gene1",'N',"N/A","N/A")[[1]]
  alldata_comp_2=alldata_compiler(alldata_2,"gene2","gene3","gene1",'N',"N/A","N/A")[[2]]
  genex_replication_prop_2=alldata_compiler(alldata_2,"gene2","gene3","gene1",'N',"N/A","N/A")[[1]]

  #Number of simulations are the minimum of: subsampling with replacement combinations of incidence of gene pair 1, pair 2, and 10000. 10000 because I don't want to do more than 1k simulations per experiment.
  nsims=min(c(choose(sum(alldata_comp_1$genex)+sum(alldata_comp_1$genex)-1,sum(alldata_comp_1$genex)),
              choose(sum(alldata_comp_2$genex)+sum(alldata_comp_2$genex)-1,sum(alldata_comp_2$genex)),1000))
  # nsubsamples=sum(alldata_comp_2$genex) #Verify that this is what it is.

  # Max subsample size is whichever integer is the minimum of:
  # number of genex in sample vs. number of NOTgenex or PC1 or PC2 in sample*genexreplication proportion. For example we had 165 PC1s, therefore 340-175=165 PC1s in our data. Since we want to sample 7.35% positive hits and 92.65% negative hits, a sample containing too many positive hits can impose an upper bound on our subsampling size.
  alldata_comp=alldata_comp_1
  genex_replication_prop=genex_replication_prop_1
  maxsubsamplesize_1=min(c(genex_replication_prop*length(alldata_comp$Positive_Ctrl1),
                           floor((length(alldata_comp$Positive_Ctrl1)-sum(alldata_comp$Positive_Ctrl1))*genex_replication_prop),
                           floor((length(alldata_comp$Positive_Ctrl1)-sum(alldata_comp$Positive_Ctrl2))*genex_replication_prop),
                           floor((length(alldata_comp$Positive_Ctrl1)-sum(alldata_comp$genex))*genex_replication_prop)))
  alldata_comp=alldata_comp_2
  genex_replication_prop=genex_replication_prop_2
  maxsubsamplesize_2=min(c(genex_replication_prop*length(alldata_comp$Positive_Ctrl1),
                           floor((length(alldata_comp$Positive_Ctrl1)-sum(alldata_comp$Positive_Ctrl1))*genex_replication_prop),
                           floor((length(alldata_comp$Positive_Ctrl1)-sum(alldata_comp$Positive_Ctrl2))*genex_replication_prop),
                           floor((length(alldata_comp$Positive_Ctrl1)-sum(alldata_comp$genex))*genex_replication_prop)))

  nsubsamples=min(c(maxsubsamplesize_1,maxsubsamplesize_2)) #minimum of either of the subsamples.
  nexperiments=7


  #####If Cohort size is greater than gene1total, run simresults_generator for that simulation####

  # if(cohort_size_vals[k]>gene1_total_vals[j]){

    simresults_pair1=simresults_generator(alldata_comp_1,7,nsims,nsubsamples,genex_replication_prop)
    simresults_pair1$gene_pair=1

    simresults_pair2=simresults_generator(alldata_comp_2,7,nsims,nsubsamples,genex_replication_prop)
    simresults_pair2$gene_pair=2

    simresults=rbind(simresults_pair1,simresults_pair2)
    simresults=as.data.frame(simresults, stringsAsFactors = F)

    # simresults$true_OR=true_or_vals[i]
    # simresults$gene1_total=gene1_total_vals[j]
    # simresults$cohort_size=cohort_size_vals[k]
    # simresults$or_pair2=or_pair2[l]
    simresults$true_OR=fet_pc1_goi$estimate
    simresults$gene1_total=goi_abundance
    simresults$cohort_size=cohort_size
    simresults$or_pair2=fet_pc1_pc2$estimate

    # colnames(simresults_compiled)=names(simresults)
    # simresults_compiled[[i+j+k+l-3]]=simresults
    # simresults_compiled[[ct]]=simresults
    # simresults_compiled=rbind(simresults_compiled,simresults)

  # }
    simresults
}
