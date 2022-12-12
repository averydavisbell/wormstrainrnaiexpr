#! /usr/bin/env/ Rscript
# Explore genes that are 'off' in one or more strains, strain/RNAi treatment data. Extension/update of initial offgenes_straintreatDE_deseq2.R, which I want to preserve.
# by Avery Davis Bell, begun 2021.12.30
require(data.table, quietly = T)
require(argparser, quietly = T)
require(DESeq2, quietly = T)
require(ashr, quietly = T)
require(ggplot2, quietly = T)
require(beeswarm, quietly = T)
require(RColorBrewer, quietly = T)

#### Functions ####
summoffgenes<-function(dt, colscheck, thresh, idcol = "gene_id"){
  # IDs rows with values in colscheck <= thresh; returns the idcol of these + number of colscheck they were under threshold in + ID of colscheck they were under threshold in
  # In: dt, data.table with data to check. must have columns idcol, all in colscheck
  #     colscheck, vector of column names in dt to check for values being <= thresh
  #     thresh, threshold for any column to be <= for row to be output
  #     idcol, column containing IDs - will output this + summary info
  # Out: data.table with columns idcols, nunder (# cols under threshold), whichunder (text of column names under threshold, comma separated)
  
  percol<-rbindlist(lapply(colscheck, function(x){
    dt[get(x) <= thresh, .(ID = get(idcol), belthresh = x)]
  }))
  setkey(percol, ID)
  out<-percol[,.(nunder = .N, whichunder = paste(belthresh, collapse = ",")), by = ID]
  
  setnames(out, "ID", idcol)
  setkeyv(out, eval(idcol))
  return(out)
}

summnonnagenes<-function(dt, colscheck, idcol = "gene_id"){
  # IDs rows with values in colscheck that are NOT NA; returns the idcol of these + number of colscheck they were not NA in + ID of colscheck they were not NA in
  # In: dt, data.table with data to check. must have columns idcol, all in colscheck
  #     colscheck, vector of column names in dt to check for values being not NA
  #     idcol, column containing IDs - will output this + summary info
  # Out: data.table with columns idcols, nnonna (# cols non NA), whichnonna (text of column names with not NAs, comma separated)
  
  percol<-rbindlist(lapply(colscheck, function(x){
    dt[!is.na(get(x)), .(ID = get(idcol), notna = x)]
  }))
  setkey(percol, ID)
  out<-percol[,.(nnonna = .N, whichnonna = paste(notna, collapse = ",")), by = ID]
  
  setnames(out, "ID", idcol)
  setkeyv(out, eval(idcol))
  return(out)
}

all1vecinother<-function(vecsearch, vectobesearched){
  # For element in first vector, do all of the (comma-separated) components of that character show up in matched element of second vector?
  #   i.e. if vecsearch[1] is 'a,b,c' and vectobesearched[1] is 'a_vs_b,b_vs_c' will be T
  # In: vecsearch, character vector. Each element is comma-separated string, each will be split out and searched against vectobesearched
  #     vectobesearched, character vector, matched to vecsearch. Each element will be searched for the elements in matching vecsearch position
  # Out: vecsearch-length logical vector: T = all vecsearch in vectobesearched at that element; F if not; ***NA if vectobesearched is NA
  
  lsearch<-strsplit(vecsearch, ",", fixed = T)
  mymatch<-sapply(1:length(lsearch), function(x){
    ifelse(is.na(vectobesearched[x]), NA,
           sapply(lsearch[[x]], grepl, vectobesearched[x]))
  })
  return(mymatch)
}

formatnocovgs<-function(nocovgf, gs, ids, swap = NULL){
  # Reads in genes' info, converts to be one row per gene in gs...see ins and outs!
  # In: nocovgf, Path to list of genes flagged as not having coverage in these strains - output of genelistsfromcovprop.R. Columns gene_id (as in input dds's), SampleID (strains), meannormval (value from filtering).
  #             One row for each strain-low coverage gene pair.
  #     gs, vector of gene_ids to get presence/absence info from nocovgf for
  #     ids, which SampleIDs in nocovgf are of interest - one output column for each of these!
  #     swap, OPTIONAL data.table for ID swapping. Columns inbed - sampleID in nocovgf; inout - how to name sample ID in output. Ie for cases where the isotype representative sequenced by CeNDR is different than strain from that isotype we've sequenced
  # Out: data.table with one row per entry in gs. Columns gene_id (keyed by this), one for each in ids (named by <id>.nocov) containing T or F: 
  #               was this strain/gene pair in nocovgf?
  
  nocovs<-fread(nocovgf, header = T)
  
  # Swap sample IDs if needed
  if(!is.null(swap)){
    for(i in 1:nrow(swap)){
      nocovs[SampleID==swap[i, inbed], SampleID:=swap[i, inout]]
    }
  }
  
  out<-data.table(gene_id = gs,
              do.call(cbind, lapply(ids, function(x){
    gs%in%nocovs[SampleID==x, gene_id]
  })))
  setnames(out, c("gene_id", paste(ids, "nocov", sep = ".")))
  setkey(out, gene_id)
  return(out)
}

summperstrain<-function(offs.zero.rnai.cov, strn){
  # Summarizes # of 'off' genes per strain, various category breakdowns
  # In: offs.zero.rnai.cov, off gene data.table with one row per gene. Must have columns -
  #           nStrainsZeroExp, # strains with MEAN expression overall at 0
  #           whichStrainsZeroExp, Comma-separated string specifying which strains have 0 expression
  #           nnocov, # of strains with no coverage at this gene (as defined however input gene no-coverage list was created)
  #           anyZeroExpIsNoCov,  T, F, or NA: are ANY of the strains called zero expression also no coverage at this gene? NA if no no-coverage strains at this gene.
  #           AllZeroExpIsNoCov,  T, F, or NA: are ALL of the strains called zero expression also no coverage at this gene? NA if no no-coverage strains at this gene.
  #           <strn>.nocov, One column per strain, T or F - was this gene flagged as no coverage (in input no coverage list) for this strain?
  #           nRNAiDE, # significant RNAi DE comparisons (one = POS1 in a strain, e.g.)
  #           whichRNAiDE, Comma-separated string specifying which RNAi comparisons have DE
  #     strn, strain ID to summarize - as in column names of offs.zero.rnai.cov, whichStrainsZeroExp
  # Out: named vector of # genes in different category. ***SEE comments when these get defined for full breakdown***
  
  # Narrow to genes 'off' in this strain
  strgns<-offs.zero.rnai.cov[grepl(strn, whichStrainsZeroExp),]
  
  out<-c("n" = nrow(strgns), # genes 'off' in this strain
         "nStrainOnly" = strgns[, sum(nStrainsZeroExp == 1)], # genes 'off' in ONLY this strain
         "nStrainPlus" = strgns[, sum(nStrainsZeroExp > 1)], # genes 'off' in at least one other strain
         "nCovInStrain" = strgns[, sum(!get(paste0(strn, ".nocov")))], # off and not called no-cov in this strain
         "nNoCovInStrain" = strgns[, sum(get(paste0(strn, ".nocov")))], # off and called no-cov in this strain
         "nStrainOnly_CovInStrain" = strgns[, sum(nStrainsZeroExp == 1 & !get(paste0(strn, ".nocov")))], # off ONLY this, not no-cov
         "nStrainOnly_NoCovInStrain" = strgns[, sum(nStrainsZeroExp == 1 & get(paste0(strn, ".nocov")))], # off ONLY this, called no-cov
         "nStrainPlus_CovInStrain" = strgns[, sum(nStrainsZeroExp > 1 & !get(paste0(strn, ".nocov")))], # off this and more strains, not no-cov
         "nStrainPlus_NoCovInStrain" = strgns[, sum(nStrainsZeroExp > 1 & get(paste0(strn, ".nocov")))], # off this and more strains, called no-cov
         "nRNAi" = strgns[, sum(!is.na(nRNAiDE))], # off in this strain AND has RNAi effect in some strain
         "nRNAi_CovInStrain" = strgns[, sum(!is.na(nRNAiDE) & !get(paste0(strn, ".nocov")))], # off + RNAi effect + not no-cov in this strain
         "nRNAi_NoCovInStrain" = strgns[, sum(!is.na(nRNAiDE) & get(paste0(strn, ".nocov")))], # off + RNAi effect + called no-cov in this strain
         "nRNAi_nStrainOnly" = strgns[, sum(!is.na(nRNAiDE) & nStrainsZeroExp == 1)], # off ONLY here + RNAi effect
         "nRNAi_nStrainPlus" = strgns[, sum(!is.na(nRNAiDE) & nStrainsZeroExp > 1)] # off this and more strains + RNAi effect
         )
  
  return(out)
}

#### Arguments & Input ####
p<-arg_parser("Explore genes that are 'off' in one or more strains, strain/RNAi treatment data", 
              name = "offgenes_straintretDE_deseq2_dnacov.R", hide.opts = TRUE)
p<-add_argument(p, "--dds",
                help = "DESeq2 strain/treatment interaction analysis object (from differentialexpr_straintreat_salmon_deseq2.R) path.
                Must be named dds (internally).
                Must provide BOTH this AND --ddsgrp.",
                default = NA)
p<-add_argument(p, "--ddsgrp",
                help = "DESeq2 dds object for GROUP model object (from differentialexpr_straintreat_salmon_deseq2.R) path.
                Must be named dds_grp (internally).
                Must provide BOTH this AND --dds.",
                default = NA)
p<-add_argument(p, "--outdir",
                help = "Output directory. Created if it doesn't exist.",
                default = "out")
p<-add_argument(p, "--outstem",
                help = "Output filestem.",
                default = "out")
p<-add_argument(p, "--alpha",
                help = "Alpha p-value threshold for FDR-like filtering *from LRT test of strain*.",
                type = "numeric",
                default = 0.1)
p<-add_argument(p, "--lfcthresh",
                help = "Log2 fold change threshold for summarizing RNAi/treatment results. Not used for filtering/multiple hypothesis testing correction, just for categorizing results passing alpha threshold only for RNAi. Default (0.5849625) corresponds to 1.5x fold change.",
                type = "numeric",
                default = 0.5849625)
p<-add_argument(p, "--nocovgenes",
                help = "Path to list of genes flagged as not having coverage in these strains - output of genelistsfromcovprop.R. Columns gene_id (as in input dds's), SampleID (strains), meannormval (value from filtering).
                One row for each strain-low coverage gene pair.",
                type = "character")
p<-parse_args(p)

#### Get DE results information to use ####
# Output directory
if(!dir.exists(p$outdir)){dir.create(p$outdir, recursive = T)}
setwd(p$outdir)

# Load
load(p$dds)
load(p$ddsgrp)

# LRT test for CTR samples - this is NEW. 
dds_ctr<-dds[,colData(dds)$Treatment=="CTR"]
design(dds_ctr)<- ~Strain
dds_ctr<-DESeq(dds_ctr, test="LRT", reduced=~1)
res_ctr_lrt<-as.data.table(results(dds_ctr))[,.(baseMean, stat, pvalue, padj)]
res_ctr_lrt[,`:=`(gene_id = rownames(dds_ctr))]
setcolorder(res_ctr_lrt, c("gene_id"))
setkey(res_ctr_lrt, gene_id)

# # LRT test for all samples to get overall Treatment p-value
# dds_rnai<-copy(dds)
# dds_rnai<-DESeq(dds_rnai, test = "LRT", reduced = ~Strain) # removes Treatment + Strain:Treatment

# Interaction results
res.int.shrink<-lapply(resultsNames(dds), function(x){
  res<-results(dds, name = x, alpha = p$alpha)
  out<-as.data.table(lfcShrink(dds, res = res, type = 'ashr'))
  out[,`:=`(gene_id=rownames(res), gene_name=rowData(dds)$locus, biotype = rowData(dds)$biotype)]
  setcolorder(out, c("gene_id", "gene_name", "biotype"))
  setkey(out, gene_id)
  return(out)
})
names(res.int.shrink)<-resultsNames(dds)

# Contrast results: among strain in CTR condition; among CTR & RNAi conditions within strain
## Build contrast pairs. Results will be grp1/grp2
strains<-levels(colData(dds_grp)$Strain)
conpairs_ctr<-rbindlist(lapply(1:(length(strains) - 1), function(x){
  rbindlist(lapply((x+1):length(strains), function(y){
    data.table(grp1 =  paste0(strains[y], ".CTR"), grp2 = paste0(strains[x], ".CTR"))
  }))
}))
conpairs_intreat<-rbindlist(lapply(strains, function(x){
  data.table(grp1 = c(paste0(x, ".PAR1"), paste0(x, ".POS1")),
             grp2 = c(paste0(x, ".CTR"), paste0(x, ".CTR")))
}))

## Get result for each contrast pairs
res.ctrstrains<-lapply(1:nrow(conpairs_ctr), function(x){
  res<-results(dds_grp, contrast = c("Group", conpairs_ctr[x, grp1], conpairs_ctr[x, grp2]), alpha = p$alpha)
  out<-as.data.table(lfcShrink(dds, res = res, type = 'ashr'))
  out[,`:=`(gene_id=rownames(res), gene_name=rowData(dds)$locus, biotype = rowData(dds)$biotype)]
  setcolorder(out, c("gene_id", "gene_name", "biotype"))
  setkey(out, gene_id)
  return(out)
})
names(res.ctrstrains)<-conpairs_ctr[,paste(grp1, "vs", grp2, sep = "_")]

res.treatinstrain<-lapply(1:nrow(conpairs_intreat), function(x){
  res<-results(dds_grp, contrast = c("Group", conpairs_intreat[x, grp1], conpairs_intreat[x, grp2]), alpha = p$alpha)
  out<-as.data.table(lfcShrink(dds, res = res, type = 'ashr'))
  out[,`:=`(gene_id=rownames(res), gene_name=rowData(dds)$locus, biotype = rowData(dds)$biotype)]
  setcolorder(out, c("gene_id", "gene_name", "biotype"))
  setkey(out, gene_id)
  return(out)
})
names(res.treatinstrain)<-conpairs_intreat[,paste(grp1, "vs", grp2, sep = "_")]

#### Get mean expression WITHIN EACH GROUP for easy cross referencing ####
grpmeanexp<-as.data.table(cbind(rownames(dds_grp),
                                do.call(cbind, lapply(levels(colData(dds_grp)$Group), function(x){
                                  rowMeans(counts(dds_grp, normalized = T)[,rownames(colData(dds_grp))[colData(dds_grp)$Group==x]])
                                }))))
setnames(grpmeanexp, c("gene_id", levels(colData(dds_grp)$Group)))
setkey(grpmeanexp, gene_id)
setcolorder(grpmeanexp, c("gene_id", paste(strains, rep(c("CTR", "PAR1", "POS1"), each = 5), sep = "."))) # put CTR first
# Also, within strains *as a whole*
grpmeanexp<-as.data.table(cbind(grpmeanexp, do.call(cbind,
                                                    lapply(strains, function(x){
                                                      rowMeans(counts(dds_grp, normalized = T)[,rownames(colData(dds_grp))[colData(dds_grp)$Strain==x]])
                                                    }))))
setnames(grpmeanexp, paste0("V", c(1:length(strains))), paste(strains, "all", sep = "."))
setkey(grpmeanexp, gene_id)

#### Identify potential 'off' genes with expression in one or more strains from LRT p value. NO LFC here. ####
destrain.ids<-res_ctr_lrt[padj < p$alpha, gene_id]
destrain<-data.table(gene_id = destrain.ids, 
                     res.ctrstrains[[1]][destrain.ids,.(gene_name, biotype)],
                     do.call(cbind, lapply(res.ctrstrains, function(x){
                       x[destrain.ids, ifelse(padj<p$alpha, log2FoldChange, NA)]
                     })),
                     grpmeanexp[destrain.ids, .SD, .SDcols = -1])
setkey(destrain, gene_id)

# Genes that are 'off' based on expression cutoffs - find
offs<-lapply(c(0, 1, 5), function(x){
  gs<-summoffgenes(destrain, colscheck = paste0(strains, ".all"), thresh = x, idcol = "gene_id")
  out<-destrain[gs]
  setcolorder(out, c("gene_id", "gene_name", "biotype", "nunder", "whichunder"))
  return(out)
})
names(offs)<-c("zeroexp", "onexp", "fivexp")

## Summarize by which strains might have no expression (& save)
offs.summ<-lapply(names(offs), function(x){
  setkey(offs[[x]], whichunder)
  out<-offs[[x]][, .(nstrains = nunder[1], ngenes = .N),by = whichunder]
  out<-out[order(nstrains, decreasing = F)]
})
names(offs.summ)<-c("zeroexp", "onexp","fivexp")
write.table(offs.summ$zeroexp, "zeroexponeormorestrains_Nsummary.txt", sep ="\t", quote = F, row.names = F)
write.table(offs.summ$onexp, "oneexponeormorestrains_Nsummary.txt", sep ="\t", quote = F, row.names = F)
write.table(offs.summ$fivexp, "fiveexponeormorestrains_Nsummary.txt", sep ="\t", quote = F, row.names = F)

#### Add RNAi results for these genes ####
# Determine significance via LRT dropping Treatment! NO NOT IMPLEMENTING YET
    # PUTTING A HOLD ON THIS
# Add log2FCs for all within-strain treatment comparisons. As previously, NA if not significant
# ##### xxxxxx DO NOT DELETE, COPY & CHANGE THAT COMMENT OUT OLD
offs.zero.rnai<-data.table(offs$zeroexp,
                           do.call(cbind, lapply(res.treatinstrain, function(x){
                             x[offs$zeroexp$gene_id, ifelse(padj<p$alpha & abs(log2FoldChange) > p$lfcthresh, log2FoldChange, NA)]
                           }))) 
setkey(offs.zero.rnai, gene_id)
summ.rnaide<-summnonnagenes(offs.zero.rnai, names(res.treatinstrain), "gene_id")
offs.zero.rnai<-summ.rnaide[offs.zero.rnai]
setcolorder(offs.zero.rnai, c("gene_id", "gene_name", "biotype", "nunder", "whichunder", "nnonna", "whichnonna"))
setnames(offs.zero.rnai, c("nunder", "whichunder", "nnonna", "whichnonna"),
         c("nStrainsZeroExp", "whichStrainsZeroExp", "nRNAiDE", "whichRNAiDE"))


#### Add DNA coverage/alignment information for these genes ####
gcov<-formatnocovgs(nocovgf = p$nocovgenes, gs = offs.zero.rnai$gene_id, ids = strains,
                    swap = data.table(inbed = "EG4349", inout = "EG4348")) ### SHOULD STANDARDIZE to swap isotype for strain
gcov[,nnocov:=rowSums(gcov[,.SD,.SDcols = -1])] # Add # strains with nocov
setcolorder(gcov, c("gene_id", "nnocov"))
offs.zero.rnai.cov<-gcov[offs.zero.rnai]
setcolorder(offs.zero.rnai.cov, c("gene_id", "gene_name", "biotype", "nStrainsZeroExp", "whichStrainsZeroExp"))
## Add if 1) all no-exp strains are no-cov from DNA; 2) ANY no-exp strains are no-cov from DNA
offs.zero.rnai.cov[nnocov>0, zeroExpIsNoCovCt:=0]
for(str in strains){
  offs.zero.rnai.cov[get(paste0(str, ".nocov"))==T & grepl(str, whichStrainsZeroExp), zeroExpIsNoCovCt:=zeroExpIsNoCovCt+1]
}
offs.zero.rnai.cov[, anyZeroExpIsNoCov:=ifelse(zeroExpIsNoCovCt>0, T, F)]
offs.zero.rnai.cov[, allZeroExpIsNoCov:=ifelse(nStrainsZeroExp==zeroExpIsNoCovCt, T, F)]
offs.zero.rnai.cov[, zeroExpIsNoCovCt:=NULL] # remove intermediate column
setcolorder(offs.zero.rnai.cov, c(names(offs.zero.rnai.cov)[1:6], c("anyZeroExpIsNoCov", "allZeroExpIsNoCov")))

# Save at this point!
write.table(offs.zero.rnai.cov, paste0(p$outdir, "/", p$outstem, "_zeroexp_rnai_nocov_genes.txt"), 
            quote = F, row.names = F, sep = "\t")

#### Summarize! ####
# Number of off genes in different categories (coverage, RNAi DE, etc) per strain
perstrsumm<-as.data.table(do.call(cbind, 
                                  lapply(strains, summperstrain, offs.zero.rnai.cov = offs.zero.rnai.cov)),
                          keep.rownames = T)
setnames(perstrsumm, c("category", strains))
write.table(perstrsumm, paste0(p$outdir, "/", p$outstem, "_zeroexp_rnai_nocov_genes_perstrain_nsummary.txt"),
            row.names = F, sep = "\t", quote = F)

# Numbers etc. All combos??
# per STRAIN
# per nocov category?
# per strain combo??

#### Session Info for reproducibility ####
cat("....offgenes_straintreatDE_deseq2_dnacov.R complete....\n")
cat("Session Information:\n")
sessionInfo()

