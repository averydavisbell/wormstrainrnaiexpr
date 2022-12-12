#! /usr/bin/env/ Rscript
# Computes coverage *per gene* from the coverage of the merged exons of that gene.
#   Also, Plots coverage across genes for multiple samples; various analysis on this. On output of mosdepthmergedexons.nf workflow processes
#   Built on/from exploregenecoverage.R
# by Avery Davis Bell, begun 2022.03.07
require(data.table, quietly = T)
require(argparser, quietly = T)
require(ggplot2, quietly = T)

#### Functions ####
freadannexonbed<-function(exonbedf, genegff, genelist = NA){
  # Reads in multi-sample/strain coverage-per-merged exon file (e.g. output of combinedpbeds process, mosdepthmergedexons.nf workflow).
  #   Gets mean coverage across gene (sum exon coverage * exon length / sum exon lengths), annotates with info from gff
  # In: exonbedf, Columns chr, start, end, name, <1 per sample containing mean coverage over exon for that sample>. **One row per MERGED EXON, so multiple rows per name (name = gene ID)
  #     genegff, genegff, Path to *genes only* gff3 file containing info on all genes in name column of --genecov input
  #     genelist, if NA, all genes read in and processed; if not NA, path to (no-header) list of gene IDs (as in name column of exonbedf) to restrict analyses to (e.g. expressed genes, protein-coding genes, etc)
  # Out: data.table with one row per gene (in genelist if provided), keyed by gene_id. with all columns of exonbedf (name now gene_id), plus, before sample columns, columns:
  #       locus, sequence_name, biotype, chr, start.gff, end.gff, strand - from GFF
  #       nMergedExons - number merged exons provided for this gene
  #       length.exons - from summing across merged exons
         
  # Subfunctions
  formatgff<-function(gff, namesget = c("Name","locus", "sequence_name", "biotype")){ # originally written in formatgff3_vits_20200803.R
    # Formats GFF3 to have useful columns only; only one column per gene
    # Input: data.table, no column names, of GFF information. *Needs to be pre-Tested for each entry to be GENE ONLY*
    #         namesget, vector of names of information to get from last column of GFF
    # Output: data.table with columns names Name, locus, biotype, chr, start, end, strand
    
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
  
  # Read in appropriate genes
  excov<-fread(exonbedf, header = T)
  if(!is.na(genelist)){
    gs<-fread(genelist, header = F)$V1
    cat(paste("....Restricting to", length(gs), "input genes from --genelist....\n"))
    excov<-excov[name%in%gs,]
  }
  setnames(excov, "name", "gene_id")
  setkey(excov, gene_id)
  
  # Compute per-gene coverage
  samps<-names(excov)[5:ncol(excov)]
  excov[,exonLn:=end - start]
  gcov<-excov[, lapply(samps, function(x) sum(get(x)*exonLn) / sum(exonLn)), by=gene_id]
  setnames(gcov, c("gene_id", samps))
  ## Keep length.exons, nMergedExons
  gcov<-excov[, .(nMergedExons = .N, length.exons = sum(exonLn)), by = gene_id][gcov]
  
  # Annotate with GFF information & return
  gffinfo<-formatgff(fread(genegff, skip = 8, header=F))[Name%in%gcov$gene_id]
  setnames(gffinfo, "Name", "gene_id")
  setkey(gffinfo, gene_id)

  return(gffinfo[gcov])
}

persampcovplots<-function(gcov, samps, datname ="Coverage (normalized to mean gene)", sname = "Strain",
                          rowdescrip = "genes"){
  # copied from exploregenecoverage.R, updated mean->median
  # Makes lots of plots summarizing gcov
  # In: gcov, data.table with columns including samps (one per character name in samps) that contain data to be faceted/plotted
  #           ***data assumed to all be POSITIVE numbers
  #     samps, vector of names of sample columns in gcov
  #     datname, name of data in the sample columns for plotting - used as axis labels etc
  #     sname, name of samples category. e.g. Strain or Sample. Used for plotting - axis labels etc
  #     rowdescrip, description of data in each row. Used for plotting. 'Number of' <this> will be histogram y axis, for example.
  # Out:  list of lists of ggplots. In each list, some subset of the following: 
  #         $plain, all data, no axis restrictions
  #         $logten, log10'ed data (0s/negatives excluded)
  #         $to5xmedian, axis restricted to 5x of *median* value across all values
  #         $tomedian, axis restricted to *median* value across all values
  #   $viols, violin plots. 
  #   $hists, histograms (faceted)
  #   $lines, line plots of ordered data (faceted)
  
  # Reformat data to repeat for each sample, ggplot2-style
  pdata<-data.table()
  for(s in samps){
    onedata<-gcov[, c(names(gcov)[!names(gcov)%in%samps], s), with = F]
    setnames(onedata, s, "dat")
    onedata[,samp:=s]
    pdata<-rbind(pdata, onedata)
  }
  
  # Violin plots
  v1<-ggplot(pdata, aes(dat, samp)) + 
    geom_violin() + xlab(datname) + ylab(sname) +
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14))
  ## log scale ***omits 0s, which are of interest, so treat with caution *******
  vlog<-v1 + scale_x_log10() +  xlab(paste(datname, "(log10, 0s excluded)")) # *****omits 0s, which are of interest, so treat with caution *******
  ## Restricted to 0-5x median coverage [violin plot is RE-COMPUTED HERE with new data - ones exceeding this excluded, so starts from scratch]
  vto5x<-ggplot(pdata, aes(dat, samp)) + 
    geom_violin(scale = "count") + xlab(paste(datname, "(to 5x median)")) + ylab(sname) + # ** scaled to COUNT here so areas could be different
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14)) +
    xlim(0, 5*pdata[,median(dat)])
  ## Restricted to median coverage or less
  vtomedian<-ggplot(pdata, aes(dat, samp)) + 
    geom_violin(scale = "count") + xlab(paste(datname, "(to median coverage)")) + ylab(sname) + # ** scaled to COUNT here so areas could be different
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14)) +
    xlim(0, pdata[,median(dat)])
  
  # Histograms
  h1<-ggplot(pdata, aes(dat)) + geom_histogram(bins = 100) +
    facet_grid(samp~.)+  xlab(datname) + ylab(paste("Number of", rowdescrip)) +
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          strip.text.y = element_text(size = 14))
  ## log scale ***omits 0s, which are of interest, so treat with caution *******
  hlog<-h1 + scale_x_log10() +  xlab(paste(datname, "(log10, 0s excluded)"))
  ## Restricted to 0-5x median coverage 
  hto5x<-h1 + xlim(0, 5*pdata[,median(dat)]) + xlab(paste(datname, "(to 5x median coverage)"))
  ## Restricted to just median coverage or less
  htomedian<-h1 + xlim(0, pdata[,median(dat)]) + xlab(paste(datname, "(to median coverage)"))
  
  # ordered line plots. nol sure yet if this is useful.
  pdata<-pdata[order(dat),]
  l1<-ggplot(pdata, aes(x = as.numeric(row.names(pdata)), y = dat)) + geom_line() +
    facet_grid(samp~.) + xlab(paste(rowdescrip, "(ascending)")) + ylab(datname) +
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          strip.text.y = element_text(size = 14))
  ## Restricted to 0-5x median coverage 
  lto5x<-l1 + ylim(0, 5*pdata[,median(dat)]) + xlab(paste(datname, "(to 5x median coverage)"))
  ## Restricted to just median coverage or less
  ltomedian<-l1 + ylim(0, pdata[,median(dat)]) + xlab(paste(datname, "(to median coverage)"))
  
  # RETURN
  return(list(viols = list(plain = v1, logten = vlog, to5xmedian = vto5x, tomedian = vtomedian),
              hists = list(plain = h1, logten = hlog, to5xmedian = hto5x, tomedian = htomedian),
              lines = list(plain = l1, to5xmedian = lto5x, tomedian = ltomedian)))
}

#### Arguments & input ####
p<-arg_parser("Computes coverage *per gene* from the coverage of merged exons of that gene. Also: Plots coverage across genes for multiple samples; various analysis on this. On output of mosdepthgenesmergedexons.nf workflow processes", 
              name = "exploregenecoverage_fromexons.R", hide.opts = TRUE)
p<-add_argument(p, "--exoncov",
                help = "Path to multi-sample/strain coverage-per-merged exon file (e.g. output of combinedpbeds process, mosdepthmergedexons.nf workflow).
                Columns chr, start, end, name, <1 per sample containing mean coverage over exon for that sample>.
                **One row per MERGED EXON, so multiple rows per name (name = gene ID)",
                type = "character")
p<-add_argument(p, "--covsumm",
                help = "Path to multi-sample/strain mosdepth summary file (e.g. output of combinedpsumms process, mosdepthmergedexons.nf workflow).
                Columns SampleID, chrom (chrom_region is for just merged exons on that chromosome; total is genome-wide), length, bases, mean, min, max",
                type = "character")
p<-add_argument(p, "--genegff",
                help = "Path to *genes only* gff3 file containing info on all genes in name column of --genecov input",
                type = "character")
p<-add_argument(p, "--outstem",
                help = "Output filestem. Include preceding path if don't want outputs in current directory.",
                type = "character",
                default = "out")
p<-add_argument(p, "--genelist",
                help = "OPTIONAL path to (no-header) list of gene IDs (as in name column of --exoncov) to restrict analyses to (e.g. expressed genes, protein-coding genes, etc)",
                type = "character",
                default = NA)
p<-parse_args(p)


out<-file.path(p$outstem)

#### Get data set up ####
cov.summ<-fread(p$covsumm, header = T)
samps<-cov.summ[,unique(SampleID)]
# Read in per exon coverage & get per gene coverage annotated with gene info
gcov<-freadannexonbed(exonbedf = p$exoncov, genegff = p$genegff, genelist = p$genelist)
## Save
write.table(gcov, gzfile(paste0(out, "_genecoveragefromexons_raw.txt.gz")), quote = F, row.names = F, sep = "\t")

# Normalize to median gene coverage
gcov.mednorm<-data.table(gcov[, names(gcov)[!names(gcov)%in%samps], with = F], 
                          gcov[,lapply(samps, function(x) get(x)/median(get(x)))])
setnames(gcov.mednorm, names(gcov))
write.table(gcov.mednorm, gzfile(paste0(out, "_genecoveragefromexons_mednorm.txt.gz")), quote = F, row.names = F, sep = "\t")

## Organized data.table info so can loop/lapply through plots etc. largely copied form exploregenecoverage.R!
dat.info<-data.table(datname = c("gcov", "gcov.mednorm"),
                     fdescrip = c("rawcoverage", "medgenenormcoverage"),
                     ldescrip = c("Raw coverage", "Coverage normalized to median gene's"))
dat.list<-list(gcov, gcov.mednorm)
names(dat.list)<-dat.info$datname

#### Analyses -  copied form exploregenecoverage.R!####
# Generate & save all within-strain plots
pdfht<-5 + length(samps)*0.5 # scale output pdf height to # samples
invisible(lapply(1:nrow(dat.info), function(x){
  plts<-persampcovplots(gcov = dat.list[[dat.info[x, datname]]], samps = samps,
                        datname = dat.info[x, ldescrip], 
                        sname = "Strain",
                        rowdescrip = "genes")
  pdf(paste(out, dat.info[x, fdescrip], "genecovviolinplots.pdf", sep = "_"), 8, pdfht)
  lapply(plts$viols, print)
  invisible(dev.off())
  
  pdf(paste(out, dat.info[x, fdescrip], "genecovhists.pdf", sep = "_"), 8, pdfht)
  lapply(plts$hists, print)
  invisible(dev.off())
  
  pdf(paste(out, dat.info[x, fdescrip], "genecovrankedlines.pdf", sep = "_"), 8, pdfht)
  lapply(plts$lines, print)
  invisible(dev.off())
  
  return(NULL)
}))

# Plot strain vs. strain** worth figuring out!!
invisible(lapply(1:nrow(dat.info), function(x){
  dat<-dat.list[[dat.info[x, datname]]]
  pdf(paste(out, dat.info[x, fdescrip], "strainVstrainplots.pdf", sep = "_"), pdfht, pdfht)
  ## All
  plot(dat[,samps,with = F], cex = 0.2, pch = 16)
  ## Axes to 5x global median
  plot(dat[,samps,with = F], cex = 0.2, pch = 16, xlim = c(0, median(as.matrix(dat[,samps,with = F]))*5),
       ylim = c(0, median(as.matrix(dat[,samps,with = F]))*5))
  ## Axes to global median
  plot(dat[,samps,with = F], cex = 0.2, pch = 16, xlim = c(0, median(as.matrix(dat[,samps,with = F]))),
       ylim = c(0, median(as.matrix(dat[,samps,with = F]))))
  invisible(dev.off())
  return(NULL)
}))

# Number genes various thresholds of proportion median coverage
ns<-sapply(seq(0, 0.5, 0.05), function(x) unlist(gcov.mednorm[,lapply(samps, function(y) sum(get(y)<=x))])) # Number genes
ps<-apply(ns, 1:2, function(x) x/nrow(gcov.mednorm)) # Proportion genes
## Format nicely & save
ns<-as.data.table(ns)
setnames(ns, paste0("nUnderPropMedianCov", seq(0, 0.5, 0.05)))
ns[,Strain:=samps]
setcolorder(ns, "Strain")
write.table(ns, paste0(out, "_ngenesunderpropofmediancoverage.txt"), sep = "\t", quote = F, row.names = F)

ps<-as.data.table(ps)
setnames(ps, paste0("propUnderPropMedianCov", seq(0, 0.5, 0.05)))
ps[,Strain:=samps]
setcolorder(ps, "Strain")
write.table(ps, paste0(out, "_propgenesunderpropofmediancoverage.txt"), sep = "\t", quote = F, row.names = F)

#### Record session info ####
cat("....exploregenecoverage_fromexons.R complete....\n")
cat("Session Information:\n")
sessionInfo()