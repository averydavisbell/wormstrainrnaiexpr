#! /usr/bin/env/ Rscript
# DESeq2 differential expression analysis with interactions from salmon data
# Script data import, some structure from differentialexpr_straintreat_salmon_deseq2.R
# by Avery Davis Bell, begun 2021.12.13
require(data.table, quietly = T)
require(argparser, quietly = T)
require(DESeq2, quietly = T)
require(tximport, quietly = T) # installed via BiocManager::install("tximport")
require(ashr, quietly = T)
require(ggplot2, quietly = T)
require(beeswarm, quietly = T)
require(RColorBrewer, quietly = T)
require(pheatmap, quietly = T)
require(ggVennDiagram, quietly = T)

#### Functions ####
tximportdata<-function(sampinfo, exampquantsf, tx2genef, tmpdir = paste0("tmp",sample.int(10e03, 1))){
  # Gets in salmon RNA quantification data via tximport
  # In: sampinfo, data.table of sample information for samples to read in. Must include column SampleID
  #     exampquantsf, example filepath to salmon quant.sf (or quant.sf.gz) RNA quantifiaction file for one sample. Transcripts in name-sorted order. Where each Sample ID goes, needs to have _sampid_ (e.g. path/to/file/_sampid__genecounts.txt.gz).
  #               For any other differences in filepath, include * for interpolation.
  #     tx2genef, Path to file mapping transcripts to genes. Two columns (transcript ID, gene ID), no header.
  #     tmpdir, directory that will be created for temporary salmon files. Created and deleted by function.
  # Out: list of - 
  #   $missingts, record of transcripts that weren't in 1 or more input files and had to be added (non-reference strains). Columns:
  #      transcript_id, fromgene (gene this comes from), nSampMissingIn (# samples that didn't have this - multiple of strains that didn't have it),
  #      ntranscriptsthisgene, # transcripts this gene has, in service of:
  #      genemissing, T or F - T if this gene only has this one transcript
  #      MissingIn, sample IDs that had transcript missing
  #   $txi, tximport output for this data. To be passed directly to DESeqDataSetFromTximport
  #   $tx2gene, transcript ID to gene ID mapping data.table. Two columns: transcript_id, gene_id
  
  # Subfunctions [from previous script!]
  getalltrs<-function(sfile, tx2gene, ofile){
    # Reads salmon file, adds trancsripts that weren't included, save which ones these were
    # Writes OUT to ofile
    
    s<-fread(sfile)
    setkey(s, Name)
    s<-s[tx2gene$tr,]
    
    missingts<-tx2gene[tr%in%s[is.na(Length), Name],]
    if(nrow(missingts)>0){
      s[which(Name%in%missingts$tr), `:=`(Length = 10, EffectiveLength = 10, TPM = 0, NumReads = 0)]
      
      write.table(s, ofile, sep = "\t", row.names = F, quote = F)
      return(list(newfile = ofile, missingts = missingts, tx2gene = tx2gene))
    }else{
      return(NULL)
    }
  }
  
  # Get all files
  myfiles<-sapply(sampinfo$SampleID, function(x){
    Sys.glob(gsub("_sampid_", x, exampquantsf))
  })
  
  # Make sure all files have all transcripts (they don't necessarily if strain-specific/different strains!)
  ## Get in transcripts
  tx2gene<-fread(tx2genef, header = F)
  setnames(tx2gene, c("tr","g"))
  setkey(tx2gene, "tr")
  ## New files with all transcripts are temporarily generated: Write TEMPORARY salmon files that have ALL transcripts, zeroing out those that don't exist in alt strains
  if(!dir.exists(tmpdir)){dir.create(tmpdir, recursive = T)}
  newfs<-lapply(names(myfiles), function(x) getalltrs(myfiles[[x]], tx2gene = tx2gene, ofile = paste0(tmpdir, "/", x, "_quant.sf")))
  ## Update to read these files where needed
  files.use<-sapply(1:length(myfiles), function(x) ifelse(is.null(newfs[[x]]), myfiles[[x]], newfs[[x]]$newfile))
  ## Format missing transcripts to save/return
  missingts<-rbindlist(lapply(1:length(myfiles), function(x){
    if(!is.null(newfs[[x]])){
      out<-newfs[[x]]$missingts
      out[, Samp:=names(myfiles)[x]]
      setnames(out, c("tr", "g"), c("transcript_id", "gene_id"))
    }else{
      return(NULL)
    }
  }))
  setkey(missingts, "transcript_id")
  missingts<-missingts[,.(gene_id = gene_id[1], nSampsMissingIn = .N, MissingIn = paste(Samp, collapse = ",")), by = "transcript_id"]
  ### Record if this is only transcript from gene
  setkey(tx2gene, g)
  tperg<-tx2gene[, .N, by = g]
  setkey(missingts, gene_id)
  missingts<-tperg[missingts]
  setnames(missingts, c("g", "N"), c("fromgene", "ntranscriptsthisgene"))
  missingts[,genemissing:=(ntranscriptsthisgene==1)]
  setcolorder(missingts, c("transcript_id", "fromgene", "nSampsMissingIn", "ntranscriptsthisgene", "genemissing"))
  
  # Read with tximport 
  setkey(tx2gene, "tr")
  names(files.use)<-names(myfiles)
  txi.salmon <- tximport(files.use, type = "salmon", tx2gene = tx2gene)
  
  # Delete TEMPORARY salmon files
  unlink(tmpdir, recursive = T)
  
  # Return
  setnames(tx2gene, c("transcript_id", "gene_id"))
  return(list(missingts = missingts, txi = txi.salmon))
}

getlocusmetadata<-function(gene_ids, gff3file){
  # Gets the locus names for each input gene_id *in gene_id order*. no_locus_name for those without locus.
  # In: gene_ids, character vector of gene_ids (from GTF) to match. These should match the "Name" information in last column of provided GFF3
  #     gff3file, path to *genes only* gff3 file containing info on all gene_ids
  # Out: data.table with one row per input gene_id. Columns: Name, locus, biotype, chr.gff, start.gff, end.gff, strand.gff [start, end, strand names are RESERVED by some DESeq2 need]
  
  # Subfunctions
  formatgff<-function(gff, namesget = c("Name","locus", "biotype")){ # originally written in formatgff3_vits_20200803.R
    # Formats GFF3 to have useful columns only; only one column per gene
    # Input: data.table, no column names, of GFF information. *Needs to be pre-Tested for each entry to be GENE ONLY*
    #         namesget, vector of names of information to get from last column of GFF
    # Output: data.table with columns names Name, llocus, biotype, chr, start, end, strand 
    
    # subfunctions
    procinfo<-function(oneinfo, namesget){
      # Breaks down last column of GFF. In: one string of this information.
      # namesget is vector of names of information from info to get. e.g. ID, Name, locus, sequence_name, etc
      eachinfo<-strsplit(oneinfo,";")[[1]]
      eachinfo<-rbindlist(lapply(strsplit(eachinfo,"="), as.list))
      setnames(eachinfo,c("name", "value"))
      out<-data.table(matrix(eachinfo[match(namesget, name), value], nrow = 1))
      setnames(out, namesget)
      return(out)
    }
    
    setnames(gff, c("chr", "source", "type", "start", "end", "exclude1", "strand", "exclude2", "info"))
    gff<-data.table(gff, gff[, procinfo(oneinfo = info, namesget = namesget),by = 1:nrow(gff)])
    return(gff[,c(namesget, "chr", "start", "end", "strand"), with = F])
  }
  
  gffinfo<-formatgff(fread(gff3file, skip = 8, header=F)) 
  gffinfo<-merge(gffinfo, data.frame(gene_ids), by.x = "Name", by.y = "gene_ids", all.y = T)
  setnames(gffinfo, c("chr", "start", "end", "strand"), paste(c("chr", "start", "end", "strand"), "gff", sep = "."))
  
  gffinfo[is.na(locus), locus:="no_locus_name"]
  return(gffinfo)
}

plotPCA_givePCs<-function (object, intgroup = "condition", ntop = 500, returnData = FALSE, xpc = 1, ypc = 2, mytitle = "") {
  # This is DESeq2's plotPCA function (accessed via ESeq2:::plotPCA.DESeqTransform) with added functionality to plot other than first 1 or 2 PCs
  # In: see ?plotPCA except for:
  #   xpc: integer. Which PC to plot on X axis.
  #   ypc: integer. Which PC to plot on Y axis
  #   mytitle: character. Title for plot.
  
  xpcname<-paste0("PC", xpc)
  ypcname<-paste0("PC", ypc)
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(xpcname = pca$x[, xpc], ypcname = pca$x[, ypc], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[xpc:ypc]
    return(d)
  }
  ggplot(data = d, aes_string(x = "xpcname", y = "ypcname", color = "group")) + 
    geom_point(size = 3) + 
    xlab(paste0(xpcname, ": ", round(percentVar[xpc] * 100), "% variance")) + 
    ylab(paste0(ypcname, ": ", round(percentVar[ypc] * 100), "% variance")) + 
    coord_fixed() + 
    ggtitle(mytitle) + 
    theme(plot.title = element_text(hjust = 0.5))
}

plotPCA_givePCs_givegenes<-function (object, intgroup = "condition", genevec, returnData = FALSE, xpc = 1, ypc = 2, mytitle = "") {
  # This is DESeq2's plotPCA function (accessed via ESeq2:::plotPCA.DESeqTransform) with added functionality to plot other than first 1 or 2 PCs
  # instead of most variable genes (ntop) takes VECTOR of gene names
  # In: see ?plotPCA except for:
  #   xpc: integer. Which PC to plot on X axis.
  #   ypc: integer. Which PC to plot on Y axis
  #   mytitle: character. Title for plot.
  
  xpcname<-paste0("PC", xpc)
  ypcname<-paste0("PC", ypc)
  
  pca <- prcomp(t(assay(object)[genevec, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(xpcname = pca$x[, xpc], ypcname = pca$x[, ypc], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[xpc:ypc]
    return(d)
  }
  ggplot(data = d, aes_string(x = "xpcname", y = "ypcname", color = "group")) + 
    geom_point(size = 3) + 
    xlab(paste0(xpcname, ": ", round(percentVar[xpc] * 100), "% variance")) + 
    ylab(paste0(ypcname, ": ", round(percentVar[ypc] * 100), "% variance")) + 
    coord_fixed() + 
    ggtitle(mytitle) + 
    theme(plot.title = element_text(hjust = 0.5))
}

plotPCA_givePCs_colshape<-function (object, colgroup = "Strain", shapegroup = "Treatment", ntop = 500,
                                    returnData = FALSE, xpc = 1, ypc = 2, mytitle = "") {
  # This is DESeq2's plotPCA function (accessed via ESeq2:::plotPCA.DESeqTransform) with added functionality to plot other than first 1 or 2 PCs
  # also expects to plot color as one variable, shape as another
  # In: see ?plotPCA except for:
  #   colgroup: what to color by
  #   shapegroup: what to shape by
  #   xpc: integer. Which PC to plot on X axis.
  #   ypc: integer. Which PC to plot on Y axis
  #   mytitle: character. Title for plot.
  
  xpcname<-paste0("PC", xpc)
  ypcname<-paste0("PC", ypc)
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(c(colgroup, shapegroup) %in% names(colData(object)))) {
    stop("the argument colgroup and shapegroup should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, c(colgroup, shapegroup), 
                                               drop = FALSE])
  
  d <- data.frame(xpcname = pca$x[, xpc], ypcname = pca$x[, ypc], 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[xpc:ypc]
    return(d)
  }
  ggplot(data = d, aes_string(x = "xpcname", y = "ypcname", color = colgroup, shape = shapegroup)) + 
    geom_point(size = 3) + 
    xlab(paste0(xpcname, ": ", round(percentVar[xpc] * 100), "% variance")) + 
    ylab(paste0(ypcname, ": ", round(percentVar[ypc] * 100), "% variance")) + 
    coord_fixed() + 
    ggtitle(mytitle) + 
    theme(plot.title = element_text(hjust = 0.5))
}

plotPCA_givePCs_colshape_givegenes<-function (object, colgroup = "Strain", shapegroup = "Treatment", genevec, returnData = FALSE, xpc = 1, ypc = 2, mytitle = "") {
  # This is DESeq2's plotPCA function (accessed via ESeq2:::plotPCA.DESeqTransform) with added functionality to plot other than first 1 or 2 PCs
  # instead of most variable genes (ntop) takes VECTOR of gene names
  # In: see ?plotPCA except for:
  #   xpc: integer. Which PC to plot on X axis.
  #   ypc: integer. Which PC to plot on Y axis
  #   mytitle: character. Title for plot.
  
  xpcname<-paste0("PC", xpc)
  ypcname<-paste0("PC", ypc)
  
  pca <- prcomp(t(assay(object)[rownames(vsd)%in%genevec, ])) # want to include any of the genes that ARE there even if some aren't
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(c(colgroup, shapegroup) %in% names(colData(object)))) {
    stop("the argument colgroup and shapegroup should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, c(colgroup, shapegroup), 
                                               drop = FALSE])
  
  d <- data.frame(xpcname = pca$x[, xpc], ypcname = pca$x[, ypc], 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[xpc:ypc]
    return(d)
  }
  ggplot(data = d, aes_string(x = "xpcname", y = "ypcname", color = colgroup, shape = shapegroup)) + 
    geom_point(size = 3) + 
    xlab(paste0(xpcname, ": ", round(percentVar[xpc] * 100), "% variance")) + 
    ylab(paste0(ypcname, ": ", round(percentVar[ypc] * 100), "% variance")) + 
    coord_fixed() + 
    ggtitle(mytitle) + 
    theme(plot.title = element_text(hjust = 0.5))
}

heatmapeucdist<-function(object, usetoname){
  # Gets Euclidean distance between samples for all genes included in object; plots these as heatmaps
  # In: object, DESeq2 object to use; must have function assay(). Usually vst-transformed counts.
  #     usetoname. Vector of column names of colData(object) to add to sample names for plotting (e.g. Strain and Treatment)
  # Out: none. Heatmap created as side effect
  
  # Get distances
  sampleDists <- dist(t(assay(object)), method = "euclidean")
  sampleDistMatrix <- as.matrix(sampleDists)
  mynames<-apply(as.data.table(colData(object)[, usetoname]), 1, function(x) paste(unlist(x), collapse = "-"))
  rownames(sampleDistMatrix) <- mynames
  colnames(sampleDistMatrix) <- mynames
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  # plot
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
}

getlrtres<-function(dds_lrt, descrip, myalpha){
  # Pulls relevant colulmns of results() for a DESeq2 likelihood-ratio test object
  # In: dds_lrt, DESeq2 likelihood-ratio test analysis output
  #     descrip, description of results - for output
  #     myalpha, p-value threshold for counting genes as signif
  # Out: list of data.tables:
  #     $res, results data.table, one row per gene, keyed by gene id. Columns:
  #         gene_id, gene ID (from rownames of dds_lrt)
  #         gene_name, gene locus name
  #         biotype, gene biotype
  #         baseMean, base mean expression of gene
  #         stat, likelihood ratio test test statistic
  #         pvalue, unadjusted p-value
  #         padj, genome-wide adjusted p-value (NA for genes with too few reads, outliers)
  #     $summ, summary of results. One row total. Columns:
  #         description, descrip
  #         nSig, number of genes with padj < myalpha
  #         nTotal, total # genes in input (same for any test)
  #         nNonNA, total # genes with adjusted p-values calculated (those with too few reads for this test, outliers are excluded) (can vary across tests)
  #         pSig_ofAll, proportion of all genes in input that are significant
  #         pSig_ofNonNA, proportion of genes with adjusted p-values calculated that are significant
  
  res<-as.data.table(results(dds_lrt, alpha = myalpha))[,.(baseMean, stat, pvalue, padj)]
  res[,`:=`(gene_id=rownames(dds_lrt), gene_name=rowData(dds_lrt)$locus, biotype = rowData(dds_lrt)$biotype)]
  setcolorder(res, c("gene_id", "gene_name", "biotype"))
  setkey(res, gene_id)

  summ<-res[, .(description = descrip, nSig = sum(padj < myalpha, na.rm = T), nTotal = nrow(res), nNonNA = sum(!is.na(padj)))]
  summ[, `:=`(pSig_ofAll = nSig/nTotal, pSig_ofNonNA = nSig/nNonNA)]
  
  return(list(res = res, summ = summ))
}

getcontrastres<-function(dds_grp, colname = "Group", grp1, grp2, descrip, myalpha, mylfcthresh){
  # Pulls ashr-shrunken results for the  contrast between grp1 and grp2, nicely annotated.
  #     Also computes one-row summary of # genes DE, in different directions, etc from these results
  # In: dds_grp, DESeq2 object with results for contrasting grp1 and grp2
  #     colname, name of column in dds_grp sample info that contains grp1 and grp2 designations
  #     grp1, name of numerator group/class to compare
  #     grp2, name of denominator group/class to comparte
  #     descrip, description of contrast for output summary
  #     myalpha, p-value threshold below which genes are considered significant
  #     mylfcthresh, absolute value log2FoldChange threshold for genes to be considered significant (when LFCs included) - just used for post-hoc filtering NOT significance adjusting
  # Out:  list of data.tables:
  #     $res, results data.table, one row per gene, keyed by gene id. Columns:
  #         gene_id, gene ID (from rownames of dds_lrt)
  #         gene_name, gene locus name
  #         biotype, gene biotype
  #         baseMean, base mean expression of gene
  #         log2FoldChange, ashr-shrunken log2FoldChange
  #         lfcSE, standard error on log2FoldChange
  #         pvalue, unadjusted p-value
  #         padj, genome-wide adjusted p value
  #     $summ, summary of results. One row total. Columns:
  #         description, descrip
  #         nSig_p, # genes significant at padj<myalpha
  #         nSig_pLogFC, # genes significant at padj<myalpha AND abs(log 2 fold change) exceeds mylfcthresh
  #         nUp_p, # upregulated genes (p value threshold)
  #         nDown_p, # downregulated genes (p value threshold)
  #         nUp_pLogFC, # upregulated genes, p-value and log2FC threshold
  #         nDown_pLogFC, # downregulated genes, p-value and log2FC threshold
  #         nTotal,total # genes in input (same for any test)
  #         nNonNA, total # genes with adjusted p-values calculated (those with too few reads for this test, outliers are excluded) (can vary across tests)
  #         pSig_p_ofAll, proportion of all genes significant with just p-value threshold
  #         pSig_p_ofNonNA, proportion of non-NA genes significant with just p-value threshold
  #         pSig_pLogFC_ofAll, proportion of all genes significant & passing log2FC threshold
  #         pSig_pLogFC_ofNonNA, proportion of non-NA genes significant & passing log2FC threshold
  #         pUp_p_ofAll, proportion of all genes that are upregulated (p value threshold)
  #         pDown_p_ofAll, proportion of all genes that are downregulated (p value threshold)
  #         pUp_pLogFC_ofAll, proportion of all genes that are upregulated (p value & log2FC thresholds)
  #         pDown_pLogFC_ofAll, proportion of all genes that are downregulated (p value & log2FC thresholds)
  
  # Pull results
  res1<-results(dds_grp, contrast = c(colname, grp1, grp2), alpha = myalpha)
  res<-as.data.table(lfcShrink(dds_grp, res = res1, type = 'ashr'))
  res[,`:=`(gene_id=rownames(res1), gene_name=rowData(dds_grp)$locus, biotype = rowData(dds_grp)$biotype)]
  setcolorder(res, c("gene_id", "gene_name", "biotype"))
  setkey(res, gene_id)
  
  # Summarize
  summ<-res[, .(description = descrip, nSig_p = sum(padj < myalpha, na.rm = T), 
                nSig_pLogFC = sum((padj < myalpha) & (abs(log2FoldChange) > mylfcthresh), na.rm = T),
                nUp_p = sum(padj < myalpha & log2FoldChange > 0, na.rm = T),
                nDown_p = sum(padj < myalpha & log2FoldChange < 0, na.rm = T),
                nUp_pLogFC =  sum((padj < myalpha) & (log2FoldChange > mylfcthresh), na.rm = T),
                nDown_pLogFC = sum((padj < myalpha) & (log2FoldChange < (-1*mylfcthresh)), na.rm = T),
                nTotal = nrow(res), nNonNA = sum(!is.na(padj)))]
  
  summ[, `:=`(pSig_p_ofAll = nSig_p/nTotal, pSig_p_ofNonNA = nSig_p/nNonNA,
              pSig_pLogFC_ofAll = nSig_pLogFC/nTotal, pSig_pLogFC_ofNonNA = nSig_pLogFC/nNonNA,
              pUp_p_ofAll = nUp_p/nTotal,
              pDown_p_ofAll = nDown_p/nTotal,
              pUp_pLogFC_ofAll = nUp_pLogFC/nTotal,
              pDown_pLogFC_ofAll = nDown_pLogFC/nTotal)] # not adding p of not NA here -just too many columns
  
  # Return
  return(list(res = res, summ = summ))
}

devenn<-function(res.list, categ.names = names(res.list), myalpha = 0.1, mylfcthresh = log2(1.5),
                 mytitle = ""){
  # Makes Venn diagrams showing number of DE genes (and upregulated, downregulated) for input reults data tables list
  # In: res.list, list of data.tables, one per Venn circle to draw. Must have columns gene_id, log2FoldChange, padj
  #     categ.names, name of elements in list/for labeling Venn diagram
  #     myalpha, p-value threshold below which genes are considered significant
  #     mylfcthresh, absolute value log2FoldChange threshold for genes to be considered significant 
  #     mytitle, main title/description for all plots, i.e. describing the comparison producing these DE results. also used in output summary table
  # Out: List of:
  #       $venns, ggplot venn diagram objects. Itself a list of elements sig, up, down (venn diagrams of DE genes, upregulated genes, downregulated genes respectively)
  #       $novs, data.table summarizing the number of genes DE in each input and how this breaks down in terms of unique vs. shared across categories (i.e. strains) 
  #         - denominator for each row is number DE in that category, not overall. Columns:
  #       test, mytitle input
  #       direction, differentially expressed, upregulated, or downregulated
  #       categ, categ.name this row describes (e.g. strain)
  #       n, number of hits/elements in this category
  #       n.unique, # hits/elements ONLY in this category
  #       n.shared, # hits/elements in this category and at least one more
  #       n.shared.1other, # hits/elements in this category and only one more
  #       n.shared.multiple, # hits/elements in this category and more than one more
  #       p.unique, proportion of hits that are unique (n is denominator for all p-columns)
  #       p.shared, proportion of hits that are shared
  #       p.shared.1other, proportion of hits shared by only one other category
  #       p.shared.multiple, proportion of hits shared my more than one other category
  #       
  
  # Subfunctions
  getnoverlap<-function(sigs, categ.names){
    # Gets number of elements that are unique vs. shared across lists in sigs
    # In: sigs, list of elements to check sharing for (e.g. gene IDs)
    #     categ.names, names of the sigs list (for output)
    # Out: data.table with one row per categ.name. Columns:
    #       categ, categ.name this row describes
    #       n, number of hits/elements in this category
    #       n.unique, # hits/elements ONLY in this category
    #       n.shared, # hits/elements in this category and at least one more
    #       n.shared.1other, # hits/elements in this category and only one more
    #       n.shared.multiple, # hits/elements in this category and more than one more
    #       p.unique, proportion of hits that are unique (n is denominator for all p-columns)
    #       p.shared, proportion of hits that are shared
    #       p.shared.1other, proportion of hits shared by only one other category
    #       p.shared.multiple, proportion of hits shared my more than one other category
    
    # Get counts, categories per ID
    perobs<-data.table(id = unlist(sigs),
                       inwhich = unlist(lapply(1:length(categ.names), function(x) rep(categ.names[x], times = length(sigs[[x]])))))
    setkey(perobs, id)
    perid<-perobs[,.(n = .N, inwhich = list(inwhich)), by = id]
    
    # For each categ.name, get sharing breakdown
    out<-rbindlist(lapply(categ.names, function(x){
      perid[,test:=sapply(1:nrow(perid), function(y) x%in%perid[y, inwhich][[1]])]
      # Number
      out<-data.table(categ = x,
                      n = perid[,sum(test)],
                      n.unique = perid[n==1, sum(test)],
                      n.shared = perid[n>1, sum(test)],
                      n.shared.1other = perid[n==2, sum(test)],
                      n.shared.multiple = perid[n>2, sum(test)])
      # Proportion OF hits in this category (of the n in this same row) that're shared
      out[,`:=`(p.unique = n.unique/n,
                p.shared = n.shared/n,
                p.shared.1other = n.shared.1other/n,
                p.shared.multiple = n.shared.multiple/n)]
      return(out)
    }))
    
    return(out)
  }
  
  # Generate gene name lists for Venns
  sigs<-lapply(res.list, function(x) x[padj < myalpha & abs(log2FoldChange) > mylfcthresh, gene_id])
  ups<-lapply(res.list, function(x) x[padj < myalpha & log2FoldChange > mylfcthresh, gene_id])
  downs<-lapply(res.list, function(x) x[padj < myalpha & log2FoldChange < (-1*mylfcthresh), gene_id])
  
  # Make Venn diagrams
  v.sig<-ggVennDiagram(sigs, categ.names, set_size = 5) + scale_fill_gradient(low="white", high = "gray30") + 
    ggtitle(mytitle, subtitle = "Differentially expressed genes") +
    theme(title = element_text(size = 17))
  v.up<-ggVennDiagram(ups, categ.names, set_size = 5) + scale_fill_gradient(low="white", high = "gray30") + 
    ggtitle(mytitle, subtitle = "Upregulated genes") +
    theme(title = element_text(size = 17))
  v.down<-ggVennDiagram(downs, categ.names, set_size = 5) + scale_fill_gradient(low="white", high = "gray30") + 
    ggtitle(mytitle, subtitle = "Downregulated genes") +
    theme(title = element_text(size = 17))
  
  # Numerical overlap summary 
  novs<-rbind(data.table(test = mytitle, direction = "differentially expressed", getnoverlap(sigs, categ.names)),
               data.table(test = mytitle, direction = "upregulated", getnoverlap(ups, categ.names)),
               data.table(test = mytitle, direction = "downregulated", getnoverlap(downs, categ.names)))
  
  # Return
  return(list(venns = list(sig = v.sig, up = v.up, down = v.down), 
              novs = novs))
}

#### Parse arguments ####
p<-arg_parser("DESeq2 differential expression analysis with interactions from salmon-quantified RNA-seq data; likelihood ratio test to get overall effects of terms", 
              name = "diffexpr_lrt_straintreat_salmon_deseq2.R", hide.opts = TRUE)
p<-add_argument(p, "--sampinfo",
                help = "Path to sample information file. Must include column SampleID and any columns that are used in modeldesign.",
                type = "character")
p<-add_argument(p, "--baseoutname",
                help = "Base name for all output files",
                type = "character",
                default = "out")
p<-add_argument(p, "--outdir",
                help = "Outer output directory. Sub-directories will be created internally.",
                type = "character")
p<-add_argument(p, "--exampquantsf",
                help = "example filepath to salmon quant.sf (or quant.sf.gz) RNA quantifiaction file for one sample. Transcripts in name-sorted order. Where each Sample ID goes, needs to have _sampid_ (e.g. path/to/file/_sampid__genecounts.txt.gz).
                For any other differences in filepath, include * for interpolation.",
                type = "character")
p<-add_argument(p, "--tx2genef",
                help = "Path to file mapping transcripts to genes. Two columns (transcript ID, gene ID), no header.",
                type = "character")
p<-add_argument(p, "--refcategoryinfo",
                help = "Path to matrix describing the reference level for each factor in the model. Columns 'colname', 'reflevel'. 'colname' must have one entry for every column of sampinfo used in the modeldesign.",
                type = "character")
p<-add_argument(p, "--modeldesign",
                help = "Quote-wrapped model formula for DESeq2, e.g. ~ Batch + Strain + Treatment + Strain:Treatment. Must include an interaction term for the analyses performed in this script.
                 Formula MUST include spaces between all terms including ~. The interaction term is used to generate progressively reduced models for likelihood-ratio tests (interaction dropped first, then the terms making up the interaction separately)
                and to create a secondary model to get condition-specific results.",
                type = "character")
p<-add_argument(p, "--genegff",
                help = "Path to *genes only* gff3 file containing info on all gene_ids present in input counts file",
                type = "character")
p<-add_argument(p, "--alpha",
                help = "Alpha p-value threshold for FDR-like filtering.",
                type = "numeric",
                default = 0.1)
p<-add_argument(p, "--lfcthresh",
                help = "Log2 fold change threshold for summarizing among-group/pairwise comparison results. Not used for LRT tests or for filtering/multiple hypothesis testing correction, just for categorizing results passing alpha threshold. 
                Default (0.5849625) corresponds to 1.5x fold change.",
                type = "numeric",
                default = 0.5849625)
p<-parse_args(p)

#### Prep data & run DESeq2 analysis ####
sampinfo<-fread(p$sampinfo)
reflevels<-fread(p$refcategoryinfo, header = T)
cat("....Preparing data and running initial DESeq2 analyses....\n")
# Import gene expression data via tximport
txidata<-tximportdata(sampinfo, exampquantsf = p$exampquantsf, tx2genef = p$tx2genef)

# Get sample metadata in DESeq format
coldata<-data.frame(copy(sampinfo)[, .SD, .SDcols=-1])
rownames(coldata)<-sampinfo$SampleID
## Convert to factor; Re-level as appropriate based on input
for(x in 1:nrow(reflevels)){
  coldata[,reflevels[x, colname]]<-factor(coldata[,reflevels[x, colname]])
  coldata[,reflevels[x, colname]]<-relevel(coldata[,reflevels[x, colname]], ref = reflevels[x, reflevel])
}

# DESeq2 format
dds<-DESeqDataSetFromTximport(txidata$txi,
                              colData = coldata,
                              design = as.formula(p$modeldesign))
## Remove genes with too few reads; those missing in some samples (due to strain-specific transcriptome)
g.few<-which(rowSums(counts(dds)) < 10)
g.miss<-which(rownames(dds) %in% txidata$missingts[genemissing==T, fromgene])
dds<-dds[-(unique(c(g.few, g.miss))), ]
cat(paste(nrow(dds), "genes remain after removal of", length(g.few), "genes with < 10 reads across samples, including",
          sum(g.miss%in%g.few),"genes with < 10 reads across samples AND missing in some samples' input;",
          sum(!g.miss%in%g.few), "genes missing in some samples' input but having enough reads.", "\n"))

## Add gene-level metadata: nice-format locus name
geneinfo<-getlocusmetadata(rownames(dds), p$genegff)
mcols(dds)<-DataFrame(mcols(dds), geneinfo)

## Run group/contrast model & save this 
splitmod<-strsplit(substr(p$modeldesign, 2, nchar(p$modeldesign)), " ")[[1]]
forgrp<-strsplit(splitmod[grep(":", splitmod)],":")[[1]] # finds interaction term & splits out components. These are for group, also are MAIN EFFECTS
dds_grp<-copy(dds)
dds_grp$Group<-factor(paste(colData(dds)[,forgrp[1]],colData(dds)[,forgrp[2]],sep=".")) # ref levels don't matter as will only look within-group
pregrp<-splitmod[-which(splitmod%in%c(forgrp, paste(forgrp, collapse = ":"), "~", "+"))]
if(length(pregrp)>0){ # model should be should be non-group stuff + group
  design(dds_grp)<-as.formula(paste("~", paste(pregrp, collapse = " + "),
                                    "+", "Group"))
}else{ # model should just be ~Group
  design(dds_grp)<- ~ Group
}

cat(paste(".........Model used for determining within-group effects:", paste(design(dds_grp), collapse=""), "\n",
          "Group is", paste(forgrp, collapse = "."), "\n")) ## print out what design used for group
cat("........Running group model....\n")
dds_grp<-DESeq(dds_grp)
save(dds_grp, file = paste0(p$outdir, "/", p$baseoutname, "_dds_group.RData"), ascii = F)

## LIKELIHOOD RATIO TESTS: all samples ##
### Drop interaction term
if(length(pregrp)>0){ # model has terms before those that were in interaction; should be those then main effects
  nointeraction<-as.formula(paste("~", paste(pregrp, collapse = " + "),
                             "+", paste(forgrp, collapse = " + ")))
}else{ # model should just be main effects
  nointeraction<-as.formula(paste("~", paste(forgrp, collapse = " + ")))
}
cat("........Running LRT for interaction effect (full model vs. model with interaction dropped)....\n")
dds_int<-DESeq(dds, test = "LRT", reduced = nointeraction) # do comparison. full is full, reduced is nointeraction
save(dds_int, file = paste0(p$outdir, "/", p$baseoutname, "_dds_LRT_interaction.RData"), ascii = F)

### Drop main terms (that were in interaction term)
if(length(pregrp)>0){
  onlymain1<-as.formula(paste("~", paste(pregrp, collapse = " + "),
                   "+", forgrp[1]))
  onlymain2<-as.formula(paste("~", paste(pregrp, collapse = " + "),
                              "+", forgrp[2]))
}else{
  onlymain1<-as.formula(paste("~", forgrp[1]))
  onlymain2<-as.formula(paste("~", forgrp[2]))
}
main1<-forgrp[1] # just for housekeeping ease
main2<-forgrp[2]

design(dds)<-nointeraction
cat(paste0("........Running LRT for ", main1, " effect (model with interaction dropped vs. model with other main effect dropped)....\n"))
dds_main1<-DESeq(dds, test = "LRT", reduced = onlymain2) # effect of main 1 is for DROPPING this, retaining only main2 
save(dds_main1,  file = paste0(p$outdir, "/", p$baseoutname, "_dds_LRT_", main1, ".RData"), ascii = F)
cat(paste0("........Running LRT for ", main2, " effect (model with interaction dropped vs. model with other main effect dropped)....\n"))
dds_main2<-DESeq(dds, test = "LRT", reduced = onlymain1) # effect of main 2 is for DROPPING this, retaining only main1 
save(dds_main2,  file = paste0(p$outdir, "/", p$baseoutname, "_dds_LRT_", main2, ".RData"), ascii = F)

## LIKELIHOD RATIO TESTS: just within reference-level samples (i.e. effect of Strain in CTR condition) ##
# After some looking, seems likely need to sample subset to do this
### Subset data
dds_main1refsonly<-dds[,as.character(colData(dds)[,main1])==reflevels[colname!=main2, reflevel]] # has only N2 samples in strain/treatment scheme
colData(dds_main1refsonly)<-droplevels(colData(dds_main1refsonly))
design(dds_main1refsonly)<-onlymain2 # Design should include only main 2 (main 1 only one level included) plus any extras e.g. Batch included in input design

dds_main2refsonly<-dds[,as.character(colData(dds)[,main2])==reflevels[colname!=main1, reflevel]] # has only CTR samples in strain/treatment scheme
colData(dds_main2refsonly)<-droplevels(colData(dds_main2refsonly))
design(dds_main2refsonly)<-onlymain1  # Design should include only main 1 (main 2 only one level included) plus any extras e.g. Batch included in input design

### Run analyses
cat((paste0("........Running LRT for ", main1, " effect in ", reflevels[colname!=main1, reflevel], " samples only....\n"))) # e.g. Strain effect in CTR only
dds_main2refsonly<-DESeq(dds_main2refsonly, test = "LRT", reduced = ~1)
save(dds_main2refsonly, file = paste0(p$outdir, "/", p$baseoutname, "_dds_LRT_", main1, "EffectIn", reflevels[colname!=main1, reflevel], "Only.RData"), ascii = F)

cat((paste0("........Running LRT for ", main2, " effect in ", reflevels[colname!=main2, reflevel], " samples only....\n"))) # e.g. Treatment effect in N2 only
dds_main1refsonly<-DESeq(dds_main1refsonly, test = "LRT", reduced = ~1)
save(dds_main1refsonly, file = paste0(p$outdir, "/", p$baseoutname, "_dds_LRT_", main2, "EffectIn", reflevels[colname!=main2, reflevel], "Only.RData"), ascii = F)

#### Overview/QC plots ####
cat("....Generating QC/overview plots....\n")
# PC plots (showing all combinations of columns of interest)
vsd<-vst(dds, blind = T)
pcadir<-paste0(p$outdir, "/pcaplots")
if(!dir.exists(pcadir)){dir.create(pcadir, recursive = T)}
for(x in reflevels$colname){
  pdf(paste0(pcadir, "/", p$baseoutname, "_pcaplot_", x, ".pdf"), 7, 5.5)
  plot1<-plotPCA_givePCs(vsd, intgroup = x, ntop = 500, returnData = F, xpc = 1, ypc = 2,
                         mytitle = paste(x, "(PCA Plot, variance stabilizing transformed, top 500 most variable genes)"))
  print(plot1)
  plot2<-plotPCA_givePCs(vsd, intgroup = x, ntop = 500, returnData = F, xpc = 2, ypc = 3,
                         mytitle = paste(x, "(PCA Plot, variance stabilizing transformed, top 500 most variable genes)"))
  print(plot2)
  plot3<-plotPCA_givePCs(vsd, intgroup = x, ntop = 500, returnData = F, xpc = 3, ypc = 4,
                         mytitle = paste(x, "(PCA Plot, variance stabilizing transformed, top 500 most variable genes)"))
  print(plot3)
  invisible(dev.off())
}
## with interaction columns plotted as shape/color on same plot
colshape<-forgrp
pdf(paste0(pcadir, "/", p$baseoutname, "_pcaplot_", paste(colshape, collapse = "_and_"), ".pdf"), 7, 5.5)
plotPCA_givePCs_colshape(vsd, colgroup = colshape[1], shapegroup = colshape[2], ntop = 500, xpc = 1, ypc = 2,
                         mytitle = "PCA Plot, variance stabilizing transformed, top 500 most variable genes")
plotPCA_givePCs_colshape(vsd, colgroup = colshape[1], shapegroup = colshape[2], ntop = 500, xpc = 2, ypc = 3,
                         mytitle = "PCA Plot, variance stabilizing transformed, top 500 most variable genes")
plotPCA_givePCs_colshape(vsd, colgroup = colshape[1], shapegroup = colshape[2], ntop = 500, xpc = 3, ypc = 4,
                         mytitle = "PCA Plot, variance stabilizing transformed, top 500 most variable genes")
invisible(dev.off())

# Heatmaps based on Euclidean distance from samples from all genes
heatdir<-paste0(p$outdir, "/heatplots")
if(!dir.exists(heatdir)){dir.create(heatdir, recursive = T)}
pdf(paste0(heatdir, "/", p$baseoutname, "_eucdistvstheatmap.pdf"), 8, 8)
heatmapeucdist(vsd, usetoname = colshape)
invisible(dev.off())


#### General summaries of results: # significant at various levels, terms ####
cat("....Generating summaries of number of DE genes....\n")
dedir<-paste0(p$outdir, "/diffexpgenes/")
if(!dir.exists(dedir)){dir.create(dedir, recursive = T)}

# LRT results 
lrt.res<-list(interaction = getlrtres(dds_int, descrip = paste("interaction - ", paste(forgrp, collapse =":")), myalpha = p$alpha),
              main1 = getlrtres(dds_main1, descrip = paste(main1, "from full model all samples"), myalpha = p$alpha),
              main2 = getlrtres(dds_main2, descrip = paste(main2, "from full model all samples"), myalpha = p$alpha),
              main1inref2 = getlrtres(dds_main2refsonly, paste(main1, "effect in", reflevels[colname!=main1, reflevel], "samples only"), myalpha = p$alpha),
              main2inref1 = getlrtres(dds_main1refsonly, paste(main2, "effect in", reflevels[colname!=main2, reflevel], "samples only"), myalpha = p$alpha))
lrt.summ<-rbindlist(lapply(lrt.res, function(x) x$summ))
write.table(lrt.summ, paste0(dedir, p$baseoutname, "_degenes_LRTs_numsummary.txt"), sep = "\t", quote = F, row.names = F)

# Contrast/Group results (all relevant ones) (using ashr logfoldchange shrinkage)
## Get all pairs to contrast (within-strain, within-treatment; denominator is reference level when possible)
m1s<-levels(colData(dds_grp)[,main1])
m2s<-levels(colData(dds_grp)[,main2])

conpairs_inm1<-rbindlist(lapply(m2s, function(m2){
  rbindlist(lapply(1:(length(m1s) - 1), function(x){
    rbindlist(lapply((x+1):length(m1s), function(y){
      data.table(inwhich = m2,
                 name1 = m1s[y],
                 name2 = m1s[x],
                 grp1 = paste(m1s[y], m2, sep = "."), 
                 grp2 = paste(m1s[x], m2, sep = "."))
    }))
  }))
}))
conpairs_inm1[,compareAmong := main1]
setcolorder(conpairs_inm1, "compareAmong")

conpairs_inm2<-rbindlist(lapply(m1s, function(m1){
  rbindlist(lapply(1:(length(m2s) - 1), function(x){
    rbindlist(lapply((x+1):length(m2s), function(y){
      data.table(inwhich = m1,
                 name1 = m2s[y],
                 name2 = m2s[x],
                 grp1 = paste(m1, m2s[y], sep = "."), grp2 = paste(m1, m2s[x], sep = "."))
    }))
  }))
}))
conpairs_inm2[,compareAmong := main2]
setcolorder(conpairs_inm2, "compareAmong")

conpairs<-rbind(conpairs_inm1, conpairs_inm2)

## Pull contrast results
con.res<-lapply(1:nrow(conpairs), function(x){
  getcontrastres(dds_grp, colname = "Group", 
                 grp1 = conpairs[x, grp1],
                 grp2 = conpairs[x, grp2],
                 descrip = paste(conpairs[x, name1], "vs", conpairs[x, name2], "in", conpairs[x, inwhich]),
                 myalpha = p$alpha,
                 mylfcthresh = p$lfcthresh)
})
names(con.res)<-conpairs[,paste(compareAmong, name1, "vs", name2, "in", inwhich, sep = "_")]
con.summ<-rbindlist(lapply(con.res, function(x) x$summ))
write.table(con.summ, paste0(dedir, p$baseoutname, "_degenes_pairwise_numsummary.txt"), sep = "\t", quote = F, row.names = F)

#### Overlap of results across terms/categories ####
cat("....Generating summaries of number of DE gene overlaps across categories....\n")
ovdir<-paste0(dedir, "overlaps/")
dir.create(ovdir)
dir.create(paste0(ovdir, main1))
dir.create(paste0(ovdir, main2))

# Among main term 1
tocomp<-unique(conpairs[compareAmong==main1, .(name1, name2)])
main1compovs<-lapply(1:nrow(tocomp), function(x){
  # Set up inputs
  res.list<-lapply(con.res[conpairs[name1==tocomp[x, name1] & name2 == tocomp[x, name2], paste(compareAmong, name1, "vs", name2, "in", inwhich, sep = "_")]],
                   function(y) y$res)
  categ.names<-conpairs[name1==tocomp[x, name1] & name2 == tocomp[x, name2], inwhich]
  
  # Generate overlaps
  ovs<-devenn(res.list, categ.names, myalpha = p$alpha, mylfcthresh = p$lfcthresh, 
              mytitle = paste(tocomp[x, name1], "vs.", tocomp[x, name2]))
  
  # Save Venn diagrams in appropriate directory
  pdf(paste0(ovdir, main1, "/", p$baseoutname, paste("", tocomp[x, name1], "vs", tocomp[x, name2], sep = "_"), 
             "_DEoverlapmultiple", main2, ".pdf"), 7, 5.5)
  invisible(lapply(ovs$venns, print))
  invisible(dev.off())
  
  return(ovs$novs)
})
# Combine & save overlap numerical summaries
write.table(rbindlist(main1compovs), paste0(ovdir, main1, "/", p$baseoutname, "_numoverlapsacross", main1, "_comparisons", ".txt"),
            sep = "\t", quote = F, row.names = F)

# Among main term 2
tocomp<-unique(conpairs[compareAmong==main2, .(name1, name2)])
main2compovs<-lapply(1:nrow(tocomp), function(x){
  # Set up inputs
  res.list<-lapply(con.res[conpairs[name1==tocomp[x, name1] & name2 == tocomp[x, name2], paste(compareAmong, name1, "vs", name2, "in", inwhich, sep = "_")]],
                   function(y) y$res)
  categ.names<-conpairs[name1==tocomp[x, name1] & name2 == tocomp[x, name2], inwhich]
  
  # Generate overlaps
  ovs<-devenn(res.list, categ.names, myalpha = p$alpha, mylfcthresh = p$lfcthresh, 
              mytitle = paste(tocomp[x, name1], "vs.", tocomp[x, name2]))
  
  # Save Venn diagrams in appropriate directory
  pdf(paste0(ovdir, main2, "/", p$baseoutname, paste("", tocomp[x, name1], "vs", tocomp[x, name2], sep = "_"), 
             "_DEoverlapmultiple", main1, ".pdf"), 7, 5.5)
  invisible(lapply(ovs$venns, print))
  invisible(dev.off())
  
  return(ovs$novs)
})
# Combine & save overlap numerical summaries
write.table(rbindlist(main2compovs), paste0(ovdir, main2, "/", p$baseoutname, "_numoverlapsacross", main2, "_comparisons", ".txt"),
            sep = "\t", quote = F, row.names = F)


#### Save DE genes for easy use in downstream categorization applications ####
cat("....Saving tables of differentially expressed genes....\n")
# LRT DE genes
## Set up clear names for saving out
goodnames.lrtres<-data.table(inlist = names(lrt.res),
                             forout = c(paste0("_interaction", main1, main2, "_sigLRTDEgenes.txt.gz"),
                                        paste0("_", main1, "_allsamplesfullmodel_sigLRTDEgenes.txt.gz"),
                                        paste0("_", main2, "_allsamplesfullmodel_sigLRTDEgenes.txt.gz"),
                                        paste0("_", main1, "EffectIn", reflevels[colname!=main1, reflevel], "Only_sigLRTDEGenes.txt.gz"),
                                        paste0("_", main2, "EffectIn", reflevels[colname!=main2, reflevel], "Only_sigLRTDEGenes.txt.gz")))
## Save out
lrtdedir<-paste0(dedir, "lrt_sigde/")
dir.create(lrtdedir)
invisible(lapply(1:nrow(goodnames.lrtres), function(x){
  sigs<-lrt.res[[goodnames.lrtres[x, inlist]]]$res[padj<p$alpha, ]
  write.table(sigs, gzfile(paste0(lrtdedir, p$baseoutname, goodnames.lrtres[x, forout])), 
              sep = "\t", quote = F, row.names = F)
}))

# Pairwise DE genes
pairdedir<-paste0(dedir, "pair_sigde/")
dir.create(pairdedir)
invisible(lapply(names(con.res), function(x){
  sigs<-con.res[[x]]$res[padj < p$alpha & abs(log2FoldChange) > p$lfcthresh,]
  write.table(sigs, gzfile(paste0(pairdedir, p$baseoutname, "_", x, "_sigPairwiseDE.txt.gz")),
              sep = "\t", quote = F, row.names = F)
}))

cat("....diffexp_lrt_straintreat_salmon_deseq2.R complete! Session information:....\n")
sessionInfo()
