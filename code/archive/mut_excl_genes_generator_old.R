mut_excl_genes_generator=function(cohort_size,incidence,or_pair1,or_pair2){
  # or_pair1=.9
  # or_pair2=.01
  #
  # incidence=40
  # cohort_size=1000

  or_p=or_pair2
  or_q=or_pair1
  q00=0.05
  p00=0.05
  ###PC1 vs GOI pair###
  p_goi=incidence/cohort_size

  q10=(1-q00)-p_goi
  q11=(p_goi*(1/q00)*or_q*((1-q00)-p_goi))/(1+(1/q00)*or_q*((1-q00)-p_goi))
  q01=(p_goi)/(1+(1/q00)*or_q*((1-q00)-p_goi))

  ###PC1 vs PC2 pair###
  p_pc1=q11+q10

  p01=(1-p00)-p_pc1
  p11=(p_pc1*(1/p00)*or_p*((1-p00)-p_pc1))/(1+((1/p00)*or_p*((1-p00)-p_pc1)))
  p10=(p_pc1)/(1+(1/p00)*or_p*((1-p00)-p_pc1))

  # gene_pair_1=(as.numeric(c(q11,q10,q01,q00)))
  # gene_pair_2=(as.numeric(c(p11,p10,p01,p00)))
  gene_pair_1=(as.numeric(c(q11,q01,q10,q00)))
  gene_pair_2=(as.numeric(c(p11,p10,p01,p00)))
  return(list(gene_pair_1,gene_pair_2))
}

# mut_excl_genes_generator2(500,12,.1,.05)
# mut_excl_genes_generator2(500,36,.1,.05)
# mut_excl_genes_generator(500,12,.1,.05)
# mut_excl_genes_generator(500,36,.1,.05)

