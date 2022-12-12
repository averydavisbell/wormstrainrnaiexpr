#! /usr/bin/env/ Rscript
# Overlapping genes flagged one way - written for as differentially expressed across strain via LRT - with no-coverage calls (from genelistsfromcovprop.R)
# by Avery Davis Bell, begun 2021.12.13
require(data.table, quietly = T)
require(argparser, quietly = T)
require(ggplot2, quietly = T)
require(ggVennDiagram, quietly = T)
require(eulerr, quietly = T)

#### Functions ####
formatpergnocovgs<-function(nocovgf, mygenes){
  # Reads in no-coverage gene info from nocovgf; retains only genes in mygenes; formats to have one row per no-coverage gene (with info on # strains it's no-cov in)
  # In: nocovgf, Path to list of genes flagged as not having coverage in strains of interest - output of genelistsfromcovprop.R. Presence in ANY strain will be enough for gene to be included in 'nocov' set - want to match 'hits' input with this carefully.
  #                   Columns gene_id (as in other inputs), SampleID (strains), meannormval (value from filtering)
  #     mygenes, vector of gene IDs of interest. Genes will be kept if they're in this set. E.g. all the genes used for a DE analysis
  # Out: data.table with one row per no-coverage gene. Columns: gene_id, nstrains (# rows this gene showed up in), whichstrains (comma-separated Strain IDs for strains that were nocov here)
  init<-fread(p$nocovgenes)[gene_id%in%mygenes, ]
  setkey(init, gene_id)
  nocovs<-init[,.(nstrains = .N, whichstrains = paste(SampleID, collapse =",")), by = gene_id]
  setkey(nocovs, gene_id)
  return(nocovs)
}

summnocovs<-function(nocovs, ntot){
  # Summarizes strain breakdown of nocov genes
  # In: nocovs, data.table with columngs gene_id, nstrains, whichstrains. Output of formatpergnocovgs()
  #     ntot, total # genes that could have been nocov - used for denominator for some proportions
  # Out: list of data.tables:
  #     $strainsumm, data.table with per-strain summaries of how that strain contributes no-coverage genes to total set. Columns:
  #             strain, strain ID
  #             n, total # no-coverage genes that were no-coverage in this strain
  #             n.unique, # no-coverage genes that were no-coverage *only* in this strain
  #             n.withothers, # no-coverage genes that were no-coverage in this strain and at least one other strain
  #             p.nocov, .nocov denotes proportion of total no-coverage gene set. Other parts of following names (and this one) as above.
  #             p.uniq.nocov,
  #             p.withothers.nocov,
  #             p.total, .total detotes proportion of total gene set (ntot). Other parts of following names (and this one) as above.
  #             p.uniq.total,
  #             p.withothers.total,
  #     $whstrains, data.table with column whichstrains, combination of strains; 
  #         N - # nocov genes in that combination
  #         p.nocov - proportion of nocov genes (in nocovs) with this combination
  #         p.total - proportion of ntot genes with this combination
  
  mystrains<-nocovs[nstrains==1, sort(unique(whichstrains))]
  whstrains<-as.data.table(nocovs[,sort(table(whichstrains), decreasing = T)])
  whstrains[,`:=`(p.nocov = N/sum(N),
                  p.total = N/ntot)]
  
  strainsumm<-rbindlist(lapply(mystrains, function(x){
    data.table(strain = x,
               n = nocovs[grepl(x, whichstrains), .N],
               n.unique = nocovs[grepl(x, whichstrains), sum(nstrains==1)],
               n.withothers = nocovs[grepl(x, whichstrains), sum(nstrains>1)])
  }))
  strainsumm[,`:=`(p.nocov = n/whstrains[,sum(N)],
                   p.uniq.nocov = n.unique/whstrains[,sum(N)],
                   p.withothers.nocov = n.withothers/whstrains[,sum(N)],
                   p.total = n/ntot,
                   p.uniq.total = n.unique/ntot,
                   p.withothers.total = n.withothers/ntot)]
  
  return(list(strainsumm = strainsumm, whstrains = whstrains))     
}

#### Parse arguments ####
p<-arg_parser("Overlapping genes flagged one way - written for as differentially expressed across strain via LRT - with no-coverage calls (from genelistsfromcovprop.R)", 
              name = "de_dnacov_overlap.R", hide.opts = TRUE)
p<-add_argument(p, "--genelist",
                help = "Path to no-header list of all genes that were considered to derive hits (e.g. included in differential expression analysis thanks to having enough coverage). 
                All of these should have also been included in coverage analysis (not checked!). Used as background/total set.",
                type = "character")
p<-add_argument(p, "--hitgenes",
                help = "Path to file containing genes that are hits to overlap with no coverage (e.g. differentially expressed genes).
                Must have column gene_id (format as in other inputs); may have other columns that will be included in output but not otherwise used here.",
                type = "character")
p<-add_argument(p, "--nocovgenes",
                help = "Path to list of genes flagged as not having coverage in strains of interest - output of genelistsfromcovprop.R. Presence in ANY strain will be enough for gene to be included in 'nocov' set - want to match 'hits' input with this carefully.
                Columns gene_id (as in other inputs), SampleID (strains), meannormval (value from filtering).
                One row for each strain-low coverage gene pair.",
                type = "character")
p<-add_argument(p, "--out",
                help = "Base out name - absolute path (include directory if not running from within desired output directory).",
                type = "character")

p<-parse_args(p)

#### Do overlaps ####
# Read in inputs
mygenes<-fread(p$genelist, header = F)$V1
hitgenes<-fread(p$hitgenes)
setkey(hitgenes, gene_id)
## Get no-cov genes; format by gene; Only keep non-covered genes that are in list of interest
nocovs<-formatpergnocovgs(nocovgf = p$nocovgenes, mygenes = mygenes)

# Summarize this nocov set w/r/t strain
nocovsumm<-summnocovs(nocovs, ntot = length(mygenes))
write.table(nocovsumm$strainsumm, paste0(p$out, "_nocovgenewhichstrainsummary.txt"), sep = "\t",
            quote = F, row.names = F)
write.table(nocovsumm$whstrains, paste0(p$out, "_nocovgenewhichstraincombosummary.txt"), sep = "\t",
            quote = F, row.names = F)

# Annotate hitgenes with no-coverage status & save
hitgenes<-nocovs[hitgenes]
setnames(hitgenes, c("nstrains", "whichstrains"), c("nocov_nstrains", "nocov_whichstrains"))
hitgenes[is.na(nocov_nstrains), nocov_nstrains:=0]
write.table(hitgenes, gzfile(paste0(p$out, "_inputhitsnocovannotated.txt.gz")), sep = "\t",
            quote = F, row.names = F)

# Save breakdown of # strains no-cov for genes that overlap; which strains
whstrainsumm<-as.data.table(hitgenes[,sort(table(nocov_whichstrains), decreasing = T)])
write.table(whstrainsumm, paste0(p$out, "_hitnocovgenewhichstrainsummary.txt"), sep = "\t",
            quote = F, row.names = F)
nstrainsumm<-as.data.table(hitgenes[,table(nocov_nstrains)])
write.table(nstrainsumm, paste0(p$out, "_hitnocovgenenstrainsummary.txt"), sep = "\t",
            quote = F, row.names = F)

# Summarize overlap with numbers
ovsumm<-hitgenes[,.(nHits = length(gene_id),
                    nNoCov = nrow(nocovs),
                    nHitsNoCov = sum(nocov_nstrains > 0),
                    nHitsNoCov1Strain = sum(nocov_nstrains == 1),
                    nHitsNoCovMultStrains = sum(nocov_nstrains > 1),
                    pAllGenesNoCov = nrow(nocovs)/length(mygenes),
                    pHitGenesNoCov = sum(nocov_nstrains > 0) / length(gene_id),
                    pAllGenesHits = length(gene_id)/length(mygenes),
                    pNoCovGenesHits = sum(nocov_nstrains > 0)/nrow(nocovs))]
## Add p value for hypergeometric test
ovsumm[,pvalueHypGeom := phyper(q = nHitsNoCov - 1, # overlap; going for p >= this
                                m = nHits, # number of total successes in one set
                                n = length(mygenes) - nHits, # number of non-successes in one set
                                k = nNoCov, # number successes in the other set
                                lower.tail = F)]
## Save
write.table(ovsumm, paste0(p$out, "_hitsnocov_overlapsummary.txt"), sep = "\t",
            quote = F, row.names = F)

#### Summary plots ####
# Venn & venn-like
## Venn diagram without scaling; nice because it has numbers but looks off
myvenn<-ggVennDiagram(list(all_genes = mygenes, hits = hitgenes$gene_id, NoDNACoverage = nocovs$gene_id),
                      category.names = c("All genes", "Hits", "No DNA coverage genes"), set_size = 5) + 
  scale_fill_gradient(low="white", high = "gray30") 
pdf(paste0(p$out, "_hitsnocov_overlapunscaledvenn.pdf"), 7.5, 6)
print(myvenn)
invisible(dev.off())

## Euler diagram with scaling; no numbers but looks good
pdf(paste0(p$out, "_hitsnocov_overlapeulerplot.pdf"), 7.5, 6)
plot(euler(list(all_genes = mygenes, hits = hitgenes$gene_id, NoDNACoverage = nocovs$gene_id)))
invisible(dev.off())

# proportion plot w/ error bars? prop all genes vs. prop DE genes that are no-cov in at least 1 strain
## Generate data for plot in appropriate forma
pdata<-data.table(descrip = c("All genes", "Hits", "Not hits"),
                  nTot = c(length(mygenes), ovsumm[,nHits], length(mygenes) - ovsumm[,nHits]),
                  nNoOv = c(ovsumm[,nNoCov], ovsumm[,nHitsNoCov], ovsumm[,nNoCov - nHitsNoCov])) # nNoOv is number OVERLAP/SUBSET
pdata<-data.table(pdata, rbindlist(lapply(1:nrow(pdata), function(x){
  res<-pdata[x, binom.test(nNoOv, nTot)]
  return(data.table(p = res$estimate, lowCI = res$conf.int[1], highCI = res$conf.int[2]))
})))
## Make and save plot
bp<-ggplot(pdata, aes(descrip, p)) + geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lowCI, ymax = highCI), width = 0.05) + 
  xlab("Gene set") + ylab("Proportion of genes that are no-DNA-coverage") +
  theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
        axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14))
pdf(paste0(p$out, "_propnocovgenesplot.pdf"), 7, 5.5)
print(bp)
invisible(dev.off())

#### Session Info for reproducibility ####
cat("....de_dnacov_overlap.R complete....\n")
cat("Session Information:\n")
sessionInfo()
