library(shiny)
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
library(shinydashboard)
library(shinycssloaders)

source("code/contab_maker.R")
source("code/alldata_compiler.R")
source("code/quadratic_solver.R")
source("code/mut_excl_genes_generator.R")
source("code/mut_excl_genes_datapoints.R")
source("code/simresults_generator.R")
source("code/shinyfunc.R")
######################Cleanup for GGPlot2#########################################
cleanup=theme_bw() +
  theme(plot.title = element_text(hjust=.5),
        panel.grid.major = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"))


# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {
  ###Intro Box
  output$intro<-renderText({
    "This is the pairwise comparison for conditional selection demo.
This method determines whether a gene of interest, GOI, is as mutually exclusive with a positive control gene as two mutually exclusive positive control genes are with each other.
This analysis requires starting off with count data on the abundance of variants in a cohort. Once the Apply Button is clicked, bootstrapping simulations of resampled positive controls are performed in real time.
<br> <br> To begin, please enter the abundances of a gene of interest and two positive control genes in the sidebar and click 'Apply'."
  })
  ###Warnings box
  output$warning_info=renderInfoBox({
    if((input$goi_abundance+input$pc1_abundance+input$pc2_abundance)>input$cohort_size){
      warnings_text=paste("Error: The cohort size is smaller than the sum of the abundances of the genes.")
      warnings_color="red"
      warnings_status="Error"
    }
    else if(input$pc1_goi>input$goi_abundance){
      warnings_text=paste("Error: The overlap between the gene of interest and the positive control gene exceed their abundance.")
      warnings_color="red"
      warnings_status="Error"
    }
    else{
      warnings_text=""
      warnings_color="green"
      warnings_status="No Warnings"
    }
    infoBox(title='Warnings',
            value =warnings_status,
            subtitle=warnings_text,
            color = warnings_color,
            icon=icon("exclamation"))
  })
  ###GOI vs PC1 Venn Diagram
  output$vennDiagram_pc1goi <- renderPlot({
    draw.pairwise.venn(area1 = as.numeric(input$pc1_abundance),
                       area2 = as.numeric(input$goi_abundance),
                       cross.area = as.numeric(input$pc1_goi),
                       category = c("Positive Control 1","Gene of Interest"),
                       fill = c("#FB8C62","#66C1A4"),
                       ext.text = TRUE,
                       ext.percent = c(0.2,0.2,0.2),
                       ext.length = 0.3,
                       label.col = rep("gray10",3),
                       lwd = 1,
                       cex = 1.5,
                       fontface = rep("bold",3),
                       fontfamily = rep("sans",3),
                       cat.cex = 1.5,
                       cat.fontface = rep("plain",2),
                       cat.fontfamily = rep("sans",2),
                       cat.pos = c(0, 0),
                       print.mode = c("percent","raw")
    )
  })

  ###PC1 vs PC2 Venn Diagram
  output$vennDiagram_pc1pc2 <- renderPlot({
    draw.pairwise.venn(area1 = as.numeric(input$pc1_abundance),
                       area2 = as.numeric(input$pc2_abundance),
                       cross.area = as.numeric(input$pc1_pc2),
                       category = c("Positive Control 1","Positive Control 2"),
                       fill = c("#FB8C62","#E689C2"),
                       ext.text = TRUE,
                       ext.percent = c(0.2,0.2,0.2),
                       ext.length = 0.3,
                       label.col = rep("gray10",3),
                       lwd = 1,
                       cex = 1.5,
                       fontface = rep("bold",3),
                       fontfamily = rep("sans",3),
                       cat.cex = 1.5,
                       cat.fontface = rep("plain",2),
                       cat.fontfamily = rep("sans",2),
                       cat.pos = c(0, 0),
                       print.mode = c("percent","raw")
    )
  })
  fet_pc1vsgoi=reactive({
    fisher.test(rbind(c(input$pc1_goi,input$pc1_abundance-input$pc1_goi),c(input$goi_abundance-input$pc1_goi,input$cohort_size-input$pc1_abundance-input$goi_abundance+input$pc1_goi)))
  })
  fet_pc1vspc2=reactive({
    fisher.test(rbind(c(input$pc1_pc2,input$pc1_abundance-input$pc1_pc2),c(input$pc2_abundance-input$pc1_pc2,input$cohort_size-input$pc1_abundance-input$pc2_abundance+input$pc1_pc2)))
  })

  ###PC1GOI box
  output$pc1goi_info=renderInfoBox({
    if(fet_pc1vsgoi()$p.value<.05){
      warnings_status=("Mutually Exclusive Genes")
      warnings_text=paste("Gene of interest is mutually exclusive with positive control 1. P value is ",signif(fet_pc1vsgoi()$p.value,3), "Pairwise Comparisons can help interpret the strength of this conditional selection when compared to a positive control gene pair with strong mutual exclusivity")
      warnings_color=("green")
    }
    else{
      warnings_status=("GOI not mutually exclusive with PC1 based on FET")
      warnings_text=paste("With a p-value of",signif(fet_pc1vsgoi()$p.value,3),"and an odds ratio of",signif(fet_pc1vsgoi()$estimate,3),", GOI is not mutually exclusive with positive control 1. Pairwise comparisons will enable statistical power to corroborate this finding. This may be especially relevant if the abundance of the gene of interest is rare.")
      warnings_color=("yellow")
      ###Add conditions for co-occurrence here
    }
    infoBox(title='Gene of Interest Gene Pair',
            value =warnings_status,
            subtitle=warnings_text,
            color = warnings_color,
            icon=icon("info"))
  })
  ###PC1PC2 box
  output$pc1pc2_info=renderInfoBox({
    if(fet_pc1vspc2()$p.value<.05){
      warnings_status=("Mutually Exclusive Positive Controls")
      warnings_text=paste("You have a reasonably mutually exclusive positive control gene pair. P value is ",signif(fet_pc1vspc2()$p.value,3), "and the odds ratio is:",signif(fet_pc1vspc2()$estimate,3), "These genes are a good positive control for pairwise comparisons.")
      warnings_color=("green")
    }
    else{
      warnings_status=("Non Mutually Exclusive Positive Controls")
      warnings_text=paste("With a P-value of",signif(fet_pc1vspc2()$p.value,3),"and an odds ratio of",signif(fet_pc1vspc2()$estimate,3),"this pair of positive controls does not have a statistically significant mutual exclusivity. Adjust your interpretation from the pairwise comparisons accordingly or consider using a different gene pair with a higher level of mutual exclusivity")
      warnings_color=("yellow")
      ###Add conditions for co-occurrence here
    }
    infoBox(title='Positive Controls Gene Pair',
            value =warnings_status,
            subtitle=warnings_text,
            color = warnings_color,
            icon=icon("info"))
  })
  ###Abundance counter info box
  output$abundance_info=renderInfoBox({
    if(input$goi_abundance<=10){
      warnings_text=paste("Attention! There is not enough statistical power to complete this analysis. You currently have", input$goi_abundance, "instances of your gene of interest. Approximately ",10-input$goi_abundance, "more events are needed to complete this analysis. Based on the frequency of the gene of interest of ",signif(input$goi_abundance*100/input$cohort_size,3),"% an additional,", (10-input$goi_abundance)*(input$cohort_size)/input$goi_abundance, "patients need to be screened to complete this analysis.")
      warnings_status="Warning"
      warnings_color="red"
    }
    else{
      warnings_text="The abundance of the gene of interest is high enough to get enough statistical power."
      warnings_status="No Warning"
      warnings_color="green"
    }
    infoBox(title='Gene Abundance Tracker',
            value =warnings_status,
            subtitle=warnings_text,
            color = warnings_color,
            icon=icon("info"))
  })

  ###Simresults generating function
  simresults_reactive=reactive({
    input$update
    # simresults=shinyfunc(500,25,100,100,20,5)
    simresults=isolate(shinyfunc(input$cohort_size,input$goi_abundance,input$pc1_abundance,input$pc2_abundance,input$pc1_goi,input$pc1_pc2))
    isolate(simresults)
  })
  ###Simresults plotting
  output$plotlyplot=renderPlot({
    input$update
    simresults=simresults_reactive()
    isolate(ggplot(simresults,aes(x=factor(gene1_total),y=log10(p_val),fill=factor(gene_pair)))+
      geom_boxplot(position = position_dodge(1))+
      facet_wrap(~true_OR,ncol=4)+
      # geom_text(data = label.df, label = "*",size=15,color="black")+
      scale_fill_brewer(palette = "Set2",name="Gene Pair",labels=c("Positive Control 1 vs GOI", "Positive Control 1 vs 2"),direction = -1)+
      cleanup+
      scale_y_continuous(name="Log(P-Value from Fisher's Test)",
                         limits = c(-3,0))+
      scale_x_discrete(name="Abundance of GOI in cohort")+
      theme(plot.title = element_text(hjust=.5),
            text = element_text(size=11,face="bold"),
            axis.title = element_text(face="bold",size="11",color="black"),
            axis.text=element_text(face="bold",size="11",color="black")))
  })
  ###KS Results BOX (maybe make this reactive so that it can be used by multiple boxes)
  # output$ksresult=renderText({
  #   input$update
  #   simresults=simresults_reactive()
  #   simresults_filtered_pair1=simresults%>%
  #     filter((exp_num==2&gene_pair==1))
  #   simresults_filtered_pair2=simresults%>%
  #     filter((exp_num==2&gene_pair==2))
  #   ks_result=ks.test(simresults_filtered_pair1$p_val,simresults_filtered_pair2$p_val,alternative = "less")$p.value
  #
  #   if(input$goi_abundance<=10)
  #     warning_text=paste("Attention! There is not enough statistical power to complete this analysis. You currently have", input$goi_abundance, "instances of your gene of interest. Approximately ",10-input$goi_abundance, "more events are needed to complete this analysis. Based on the frequency of the gene of interest of ",signif(input$goi_abundance*100/input$cohort_size,3),"% an additional,", (10-input$goi_abundance)*(input$cohort_size)/input$goi_abundance, "patients need to be screened to complete this analysis.")
  #   else
  #     warning_text="The abundance of the gene of interest is high enough to get enough statistical power."
  #   if(ks_result<=.05)
  #       isolate(return(paste(warning_text,"Result: There is a significant difference between the distributions of the two gene pairs. A one-tailed KS-Test yielded a p-value of:",signif(ks_result,3))))
  #   else
  #     isolate(paste(warning_text,"There is no significant difference between these two distributions. A one-tailed KS-Test yielded a p-value of:",signif(ks_result,3)))
  # })

  ###Text box asking users to make sure to click apply
  output$take_action=renderText({
    "Please make sure you click Apply to refresh these results"
  }
  )
  ###Simresults generating function
  # ksresults_reactive=reactive({
  #   input$update
  #   simresults=simresults_reactive()
  #   simresults=isolate(shinyfunc(input$cohort_size,input$goi_abundance,input$pc1_abundance,input$pc2_abundance,input$pc1_goi,input$pc1_pc2))
  #   isolate(simresults)
  # })
  ###Infobox about KS results
  output$ks_info=renderInfoBox({
    input$update
    simresults=simresults_reactive()
    simresults_filtered_pair1=simresults%>%
      filter((exp_num==2&gene_pair==1))
    simresults_filtered_pair2=simresults%>%
      filter((exp_num==2&gene_pair==2))
    ks_result=isolate(ks.test(simresults_filtered_pair1$p_val,simresults_filtered_pair2$p_val,alternative = "less")$p.value)
    # ks_result=isolate(ks.test(simresults_filtered_pair1$p_val,simresults_filtered_pair2$p_val,alternative = "two.sided")$p.value)
    # ks_result=isolate(ks.test(simresults_filtered_pair1$p_val,simresults_filtered_pair2$p_val)$p.value)
    # ks.test()

    if(ks_result<=.005){
      warnings_text=paste("Result: There is a significant difference between the distributions of the two gene pairs. A one-tailed KS-Test yielded a p-value of:",signif(ks_result,3))
      warnings_status="Confirmation of a lack of conditional selection"
      warnings_color="green"
    }
    else{
      warnings_text=paste("There is no significant difference between the distribution of gene of interest vs positive control 1 and positive control 1 vs positive control 2. A one-tailed KS-Test yielded a p-value of:",signif(ks_result,3))
      warnings_status="Confirmation of positive conditional selection"
      warnings_color="red"
    }
    isolate(infoBox(title='Pairwise Distributions',
            value =warnings_status,
            subtitle=warnings_text,
            color = warnings_color,
            icon=icon("poll")))
    })
  ks_debug=renderText({
    return(ks_result$p.value())
  })
})
