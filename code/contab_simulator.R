contab_simulator=function(contab,nsims,cohort_size){
  # Another function that takes in a contingency table, and resamples it based on its probabilities on a multinomial distribution. Number of resampling simulations are an input of the function
  contabs_resampled=rmultinom(nsims, #number of datasets to generate
                              cohort_size, #sample size of the total population
                              c(contab[1,1],
                                contab[1,2],
                                contab[2,1],
                                contab[2,2])) #probability of positive event
  contabs_resampled_t=t(contabs_resampled)
  colnames(contabs_resampled_t)=c("p11","p10","p01","p00")
  contabs_resampled_t
}
