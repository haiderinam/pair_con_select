mut_excl_genes_generator=function(cohort_size,incidence,or_pair1,or_pair2){
  # source("quadratic_solver.R")
  ###Scenario 1
  d_c=1 ### odds of d in c. Probability, therefore, is x/(1+x) for scenario 1
  a=or_pair1*incidence/(1+or_pair1) ###Fix rounding
  b=incidence-a ### Fix rounding
  c=(cohort_size-incidence)*d_c/(1+d_c)
  d=c*d_c
  ###Scenario 2
    ###Makeshift solution to prevent quadratic errors:
    c2=(cohort_size-a-c)/2
    d2=c2 ##or d2=c2*d_c
    a2=(a+c)*or_pair2/(1+or_pair2)
    b2=a+c-a2 ##look at pic taken on 12/13 for map of this
  # c2=c #aka positive control 1 and 2 are the same size
  # b2=c #aka positive control 1 and 2 are the same size
  # if((b2+c2)>=cohort_size){
  #   return("Positive Control1 too big")
  # }
  # if(((-(cohort_size-b2-c2))^2)-(4*1*(b2*c2*or_pair2))<0){ #b2-4ac. Note that a,b,c are coefficients of the quadratic
  #   return("imaginary roots for a2")
  # }
  #
  # a2=quadratic_solver(1,b2+c2-cohort_size,(c2*b2*or_pair2))
  # a2=min(a2)
  # d2=cohort_size-b2-c2
  # gene_pair_1=round(signif(as.numeric(c(a,b,c,d)),1))
  # gene_pair_2=round(signif(as.numeric(c(a2,b2,c2,d2)),1))

  # gene_pair_1=round(as.numeric(c(a,b,c,d)))
  # gene_pair_2=round(as.numeric(c(a2,b2,c2,d2)))

  gene_pair_1=(as.numeric(c(a,b,c,d))/cohort_size)
  gene_pair_2=(as.numeric(c(a2,b2,c2,d2))/cohort_size)
  return(list(gene_pair_1,gene_pair_2))
}
