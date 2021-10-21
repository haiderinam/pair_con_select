contab_downsampler<- function(contab,prob) {
  p11=contab[1,1]
  p10=contab[1,2]
  p01=contab[2,1]
  p00=contab[2,2]
  contab_downsampled=contab ###Just making a dummy contingency table
  contab_downsampled[1,1]=prob*p11/(p10+p11)
  contab_downsampled[1,2]=prob*p10/(p10+p11)
  contab_downsampled[2,1]=(1-prob)*p01/(p00+p01)
  contab_downsampled[2,2]=(1-prob)*p00/(p00+p01)
  contab_downsampled
}
