contab_maker<- function(gene1,gene2,arrayname) {
  gene1_gene2_nn=nrow(arrayname[gene1==0 & gene2==0,]) #N N
  gene1_gene2_ny=nrow(arrayname[gene1==0 & gene2!=0,]) #N Y
  gene1_gene2_yn=nrow(arrayname[gene1!=0 & gene2==0,]) #Y N
  gene1_gene2_yy=nrow(arrayname[gene1!=0 & gene2!=0,]) #Y Y
  
  con_tab_gene1_gene2=rbind(c(gene1_gene2_nn,gene1_gene2_ny),c(gene1_gene2_yn,gene1_gene2_yy))
  con_tab_gene1_gene2
}
# 
# #Contingency Tables:
# con_tab_pctrl1_genex=contab_maker(alldata_comp$Positive_Ctrl1,alldata_comp$genex,alldata_comp)
# con_tab_pctrl2_genex=contab_maker(alldata_comp$Positive_Ctrl2,alldata_comp$genex,alldata_comp)
# con_tab_pctrl1_pctrl2=contab_maker(alldata_comp$Positive_Ctrl1,alldata_comp$Positive_Ctrl2,alldata_comp)
# con_tab_pctrl1_rndmarray=contab_maker(alldata_comp$Positive_Ctrl1,alldata_comp$rndmarray,alldata_comp)
# con_tab_pctrl2_rndmarray=contab_maker(alldata_comp$Positive_Ctrl2,alldata_comp$rndmarray,alldata_comp)
# 
# #Fisher's Exact Test
# p_p1_genex=fisher.test(con_tab_pctrl1_genex)
# p_p2_genex=fisher.test(con_tab_pctrl2_genex)
# p_p1_p2=fisher.test(con_tab_pctrl1_pctrl2)
# p_p1_rndm=fisher.test(con_tab_pctrl1_rndmarray)
# p_p2_rndm=fisher.test(con_tab_pctrl2_rndmarray)