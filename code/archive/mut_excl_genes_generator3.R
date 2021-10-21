mut_excl_genes_generator3=function(cohort_size,incidence,or_pair1,or_pair2,d_pair1,d_pair2){
  # cohort_size=500
  # incidence=25 #aka incidence of the gene of interest
  # or_pair1=.9
  # or_pair2=.1
  # d_pair1=50
  # d_pair2=50


  c1=cohort_size-d_pair1-incidence
  a1=or_pair1*incidence/(or_pair1+(d_pair1/c1))
  b1=incidence-a1

  incidence_pc1=a1+c1
  b2=cohort_size-incidence_pc1-d_pair2
  a2=or_pair2*incidence_pc1/(or_pair2+(d_pair2/b2))
  c2=incidence_pc1-a2

  gene_pair_1=(as.numeric(c(a1,b1,c1,d_pair1))/cohort_size)
  gene_pair_2=(as.numeric(c(a2,b2,c2,d_pair2))/cohort_size)
  return(list(gene_pair_1,gene_pair_2))
}
