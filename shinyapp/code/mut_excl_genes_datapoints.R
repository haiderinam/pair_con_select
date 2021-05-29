mut_excl_genes_datapoints=function(gene_pair){
  a=gene_pair[1]
  b=gene_pair[2]
  c=gene_pair[3]
  d=gene_pair[4]
  columns=c(1:(a+b+c+d))
  all_sim_data=data.frame(columns)
  all_sim_data$gene1=0 #gene 1 is like alk
  #Non mutually exclusive alk hits
  all_sim_data$gene1[sample(c(1:(a+c)),a,replace = F,prob = NULL)]=1 #note that this is sampling a number of events
  #Mutually exclusive alk hits
  all_sim_data$gene1[sample(c((a+c+1):(a+b+c+d)),b,replace = F,prob = NULL)]=1 #note that this is sampling a number of events
  all_sim_data$gene2=0 #gene 2 is like braf
  all_sim_data$gene2[c(1:(a+c))]=1
  all_sim_data$gene3=0
  all_sim_data$gene3[floor(c((a+c-5):(a+b+c+d)))]=1
  all_sim_data
}
