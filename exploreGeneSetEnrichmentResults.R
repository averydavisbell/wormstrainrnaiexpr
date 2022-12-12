#! /usr/bin/env/ Rscript
# Pull in, plot gene set enrichment results from Wormbase. One-off, potentially
# by Avery Davis Bell, begun 2022.01.11
require(data.table, quietly = T)
require(ggplot2, quietly = T)

#### Functions ####
gotermplot<-function(pdata, facetby = "Strain", mytitle = ""){
  # Makes gotermplot. NOT flexible, for hard-coded strains currently. GO terms ordered by median highest significance [but NAs removed - so being high in 1 strain is enough to have it be in top]
  # In: pdata, data to plot. Subset of allgsea as created this script. All rows will be plotted so subset before hand and with facet
  #     facetby, column name of pdata to facet on
  #     mytitle, title of plot
  # Out: ggplot2 object
  pdata<-copy(pdata)
  if(facetby=="Strain"){
    pdata$Strain<-factor(pdata$Strain, levels = c("N2", "JU1088", "EG4348", "CB4856", "QX1211")) # ref strain first then ordered by RNAi sensitivity
  }
  # Order GO terms by mean significance
  pdata$GOTerm<-factor(pdata$GOTerm, levels = levels(reorder(pdata$GOTerm, pdata$negLog10Q, function(x) median(x, na.rm = T)))) 
  
  # plot
  plt<-ggplot(pdata, aes(GOTerm, negLog10Q)) + 
    geom_bar(stat = "identity", aes(fill=GOTerm)) + coord_flip() +
    facet_grid(~eval(as.name(facetby))) + ggtitle(mytitle) + ylab("-log10(Q value)") + xlab("") +
    theme(legend.position = "none", axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          strip.text.x = element_text(size = 15),
          title = element_text(size = 16)) 
  return(plt)
}


#### Get in data ####
mydir<-"~/Dropbox (GaTech)/BioSci-Paaby-Ext/DifferentialExpressionRNAi/20220105_ADB_RNAiRNASeqSalmonStrainSpec_lrts_ws276_20210121cendr/enrichment/20220111_initial/"
myodir<-"~/Dropbox (GaTech)/BioSci-Paaby-Ext/DifferentialExpressionRNAi/20220105_ADB_RNAiRNASeqSalmonStrainSpec_lrts_ws276_20210121cendr/enrichment/20220111_initial/alltogether/"

compsets<-data.table(Strain = rep(c("N2", "JU1088", "EG4348", "CB4856", "QX1211"), each = 6),
                     Comparison = rep(rep(c("PAR1vsCTR", "POS1vsCTR"), each = 3), 5),
                     Direction =rep(c("sig", "up", "down"), 10))
compsets[,resfile:=paste0(mydir,
                          Strain, "_", Comparison, "/", Direction, "_results.txt")]
allgsea<-rbindlist(lapply(1:nrow(compsets), function(x){
  if(file.exists(compsets[x, resfile])){
    res<-fread(compsets[x, resfile], sep = "\t")
    # split out GO string
    res[,GOTerm:=unlist(lapply(strsplit(res$Term, split = " ", fixed = T), function(y) paste(y[1:(length(y) -1)], collapse = " ")))]
    
    # add metadata
    res[,`:=`(Strain = compsets[x, Strain],
              Comparison = compsets[x, Comparison],
              Direction = compsets[x, ifelse(Direction=="sig", "All DE", Direction)])]
    setcolorder(res, c("Strain", "Comparison", "Direction"))
    return(res)
  }else{
    return(NULL)
  }
}))
allgsea[,negLog10Q:=(-1*log10(`Q value`))]

#**manually shorten any super-long GO terms. NOT elegant.
#
allgsea[GOTerm=="oxidoreductase activity acting on paired donors with incorporation or reduction of molecular oxygen reduced flavin or flavoprotein as one donor and incorporation of one atom of oxygen",
        GOTerm:="oxidoreductase activity...flavin/flavoprotein..."]
allgsea[GOTerm=="biological process involved in interspecies interaction between organisms",
        GOTerm:="biol. process involved in interspecies interaction"]
allgsea[GOTerm=="SCF-dependent proteasomal ubiquitin-dependent protein catabolic process",
        GOTerm:="SCF-dependent proteasomal ubiq.-dependent prot. catabolic proc."]


# save combined results table!
write.table(allgsea, paste0(myodir, "allstrainstreatmentsdirstogether_GEAresults.txt"),
            sep = "\t", quote = F, row.names = F)

#### Plot ####
# Per comparison and direction
par1up<-gotermplot(pdata = allgsea[Comparison=="PAR1vsCTR" & Direction=="up",], facetby = "Strain",
                   mytitle = "PAR1 upregulated genes")
par1down<-gotermplot(pdata = allgsea[Comparison=="PAR1vsCTR" & Direction=="down",], facetby = "Strain",
                     mytitle = "PAR1 downregulated genes")
pos1up<-gotermplot(pdata = allgsea[Comparison=="POS1vsCTR" & Direction=="up",], facetby = "Strain",
                   mytitle = "POS1 upregulated genes")
pos1down<-gotermplot(pdata = allgsea[Comparison=="POS1vsCTR" & Direction=="down",], facetby = "Strain",
                             mytitle = "POS1 downregulated genes")
pdf(paste0(myodir, "gotermbarcharts_comparisondirection_strainfacets.pdf"), 12, 11)
print(par1up)
print(par1down)
print(pos1up)
print(pos1down)
invisible(dev.off())

# *** go back to manually shorten after min # plots generated

#### Summarize # times different terms show up? Or?
