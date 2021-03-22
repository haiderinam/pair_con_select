simresults_generator=function(alldata_comp,nexperiments,nsims,nsubsamples,genex_replication_prop){

  all_genex     =alldata_comp[alldata_comp$genex!=0,] #These will be used later for sampling
  all_not_genex =alldata_comp[alldata_comp$genex==0,]
  all_pc1       =alldata_comp[alldata_comp$Positive_Ctrl1!=0,]
  all_not_pc1   =alldata_comp[alldata_comp$Positive_Ctrl1==0,]
  all_pc2       =alldata_comp[alldata_comp$Positive_Ctrl2!=0,]
  all_not_pc2   =alldata_comp[alldata_comp$Positive_Ctrl2==0,]
  all_rndm      =alldata_comp[alldata_comp$rndmarray!=0,]
  all_not_rndm  =alldata_comp[alldata_comp$rndmarray==0,]

  ct=1
  simresults<-matrix(nrow=nsims*1*nexperiments,ncol=7)
  # subsamplenumber=10
  subsamplenumber=nsubsamples
  # for (subsamplenumber in 1:nsubsamples){
    simresults[c(ct:((ct+(nsims*nexperiments))-1)),2]=subsamplenumber #Subsample Number: This updates the next n rows with the subsample number that the loop is on. n is calculated by number of experiments * number of simulations
    subsamplenumbernegativehits=round(subsamplenumber/genex_replication_prop)-subsamplenumber #This determines the number of non-hits for the sample
    for (simnumber in 1:nsims){
      alldata_comp$rndmarray=rbinom(length(alldata_comp$Positive_Ctrl1),1,genex_replication_prop) #Creates random array. Note: we're updating random array at each simulation within each subsample size.

      simresults[c(ct:((ct+nexperiments)-1)),3]=simnumber #Simulation Number. This updates the next 6 rows with the simulation number that the loop is on.

      sample_genex     =all_genex[sample(length(all_genex$genex),subsamplenumber,replace=T,prob=NULL),]
      sample_not_genex =all_not_genex[sample(length(all_not_genex$genex),subsamplenumbernegativehits,replace=T,prob=NULL),]
      sample_pc1       =all_pc1[sample(length(all_pc1$Positive_Ctrl1),subsamplenumber,replace=T,prob=NULL),]
      sample_not_pc1   =all_not_pc1[sample(length(all_not_pc1$Positive_Ctrl1),subsamplenumbernegativehits,replace=T,prob=NULL),]
      sample_pc2       =all_pc2[sample(length(all_pc2$Positive_Ctrl2),subsamplenumber,replace=T,prob=NULL),]
      sample_not_pc2   =all_not_pc2[sample(length(all_not_pc2$Positive_Ctrl2),subsamplenumbernegativehits,replace=T,prob=NULL),]
      sample_rndm      =all_rndm[sample(length(all_rndm$rndmarray),subsamplenumber,replace=T,prob=NULL),]
      sample_not_rndm  =all_not_rndm[sample(length(all_not_rndm$rndmarray),subsamplenumbernegativehits,replace=T,prob=NULL),]

      #Combining n (subsample number) hits and m (subsamplenegativenumber) non-hits
      sample_genex_comb =rbind(sample_genex,sample_not_genex) #Dunno if this is the right way to do it
      sample_pc1_comb   =rbind(sample_pc1,sample_not_pc1)
      sample_pc2_comb   =rbind(sample_pc2,sample_not_pc2)
      sample_rndm_comb  =rbind(sample_rndm,sample_not_rndm)

      # alldata_comp_sample<-as.data.frame(cbind(sample_pc1_comb,sample_pc2_comb,sample_genex_comb,sample_rndm_comb))  #Creating array with samples so that contab_maker can use them

      #contingency tables
      con_tab_sample_pctrl1_genex=contab_maker(sample_pc1_comb$Positive_Ctrl1,sample_pc1_comb$genex,sample_pc1_comb) #Double check to see if this is the right way to do it
      ###
      con_tab_sample_genex_pctrl1=contab_maker(sample_genex_comb$Positive_Ctrl1,sample_genex_comb$genex,sample_pc1_comb) #change this to genex_comb
      p_s_p1_genex=fisher.test(con_tab_sample_pctrl1_genex)

      p_s_genex_p1=fisher.test(con_tab_sample_genex_pctrl1)
      # p_s_p1_genex$p.value
      # p_s_genex_p1$p.value
      ###
      con_tab_sample_pctrl2_genex=contab_maker(sample_pc2_comb$Positive_Ctrl2,sample_pc2_comb$genex,sample_pc2_comb)
      con_tab_sample_pctrl1_pctrl2=contab_maker(sample_pc1_comb$Positive_Ctrl1,sample_pc1_comb$Positive_Ctrl2,sample_pc1_comb) #note how sampling pc1 vs pc2 and then finding pc2 vs pc1 has drastically different p-values. Add code that looks at which of the samples has a smaller size and then chooses to sample for 7% of that first
      con_tab_sample_pctrl1_pctrl2=contab_maker(sample_pc2_comb$Positive_Ctrl2,sample_pc2_comb$Positive_Ctrl1,sample_pc2_comb)
      con_tab_sample_pctrl1_rndmarray=contab_maker(sample_pc1_comb$Positive_Ctrl1,sample_pc1_comb$rndmarray,sample_pc1_comb)
      con_tab_sample_pctrl2_rndmarray=contab_maker(sample_pc2_comb$Positive_Ctrl2,sample_pc2_comb$rndmarray,sample_pc2_comb)
      con_tab_sample_genex_rndmarray=contab_maker(sample_genex_comb$genex,sample_genex_comb$rndmarray,sample_genex_comb)

      #Fishers exact test
      p_s_p1_genex=fisher.test(con_tab_sample_pctrl1_genex) #p_s_ stands p-value, sample. Can name these better in the future
      p_s_p2_genex=fisher.test(con_tab_sample_pctrl2_genex) #Could check if pc1_genex and genex_pc1 give the same p-value just to verify our test
      p_s_p1_p2=fisher.test(con_tab_sample_pctrl1_pctrl2)
      p_s_p1_rndm=fisher.test(con_tab_sample_pctrl1_rndmarray)
      p_s_p2_rndm=fisher.test(con_tab_sample_pctrl2_rndmarray)
      p_s_genex_rndm=fisher.test(con_tab_sample_genex_rndmarray)

      #Creating df that has all the simulation data so far. This will be used by the simresults df later.
      # alldatasamplepvals<-as.data.frame(cbind(p_s_p1_genex$p.value,p_s_p2_genex$p.value,p_s_p1_p2$p.value,p_s_p1_rndm$p.value,p_s_p2_rndm$p.value,p_s_genex_rndm$p.value))
      # colnames(alldatasamplepvals)=c("p_s_p1_genex","p_s_p2_genex","p_s_p1_p2","p_s_p1_rndm","p_s_p2_rndm","p_s_genex_rndm")
      # alldatasampleORvals<-as.data.frame(cbind(p_s_p1_genex$estimate,p_s_p2_genex$estimate,p_s_p1_p2$estimate,p_s_p1_rndm$estimate,p_s_p2_rndm$estimate,p_s_genex_rndm$estimate))
      # colnames(alldatasampleORvals)=c("p_s_p1_genex","p_s_p2_genex","p_s_p1_p2","p_s_p1_rndm","p_s_p2_rndm","p_s_genex_rndm")
      alldatasamplepvals<-as.data.frame(cbind(p_s_p1_genex$p.value,p_s_genex_p1$p.value,p_s_p2_genex$p.value,p_s_p1_p2$p.value,p_s_p1_rndm$p.value,p_s_p2_rndm$p.value,p_s_genex_rndm$p.value))
    colnames(alldatasamplepvals)=c("p_s_p1_genex","p_s_genex_p1","p_s_p2_genex","p_s_p1_p2","p_s_p1_rndm","p_s_p2_rndm","p_s_genex_rndm")
      alldatasampleORvals<-as.data.frame(cbind(p_s_p1_genex$estimate,p_s_genex_p1$estimate,p_s_p2_genex$estimate,p_s_p1_p2$estimate,p_s_p1_rndm$estimate,p_s_p2_rndm$estimate,p_s_genex_rndm$estimate))
      colnames(alldatasampleORvals)=c("p_s_p1_genex","p_s_genex_p1","p_s_p2_genex","p_s_p1_p2","p_s_p1_rndm","p_s_p2_rndm","p_s_genex_rndm")
      # for (expnumber in 1:nexperiments) {
      #  Ideally, all the code above that assumes there are 3 genes+ a random array would be fed in to this for-loop that just does all the pairwise combinations automatically
      # }

      for (expnumber in 1:nexperiments) {
        simresults[ct,1]=ct #Total Count
        simresults[ct,4]=expnumber #Experiment Number
        simresults[ct,5]=colnames(alldatasamplepvals[expnumber]) #Experiment Name
        simresults[ct,6]=alldatasamplepvals[1,expnumber] #Experiment P-value. Trying to convert it into numerics
        simresults[ct,7]=as.numeric(alldatasampleORvals[1,expnumber]) #Experiment OR-value

        ct=ct+1
      }
    }
  # }
  #Renaming columns of simresults and transforming simresults
  colnames(simresults)=c("totCt","subsample_size","sim_num","exp_num","exp_name","p_val","OR_val")
  simresults=as.data.frame(simresults, stringsAsFactors = F)

  simresults[c("totCt","subsample_size","sim_num","exp_num","p_val","OR_val")]=lapply(simresults[c("totCt","subsample_size","sim_num","exp_num","p_val","OR_val")], as.character)
  simresults[c("totCt","subsample_size","sim_num","exp_num","p_val","OR_val")]=lapply(simresults[c("totCt","subsample_size","sim_num","exp_num","p_val","OR_val")], as.numeric)

  return(simresults)
  }
