#! /usr/bin/env/ Rscript
# Explores breaking genes into truly no coverage etc, along with overlaps with off genes DE. (written for 5 strains x 3 conditions RNA-seq data)
# Also pulls in per-gene divergence metrics as output by nucdivcendr_geneswindows_allandasestrains.R 

# by Avery Davis Bell, begun 2023.01.03
require(data.table, quietly = T)
require(argparser, quietly = T)
require(ggplot2, quietly = T)

#### Functions ####
covhistogram<-function(gcov, samps, datname ="Coverage", sname = "Strain",
                       rowdescrip = "genes", vline = F){
  # Makes lots of plots summarizing gcov
  # In: gcov, data.table with columns including samps (one per character name in samps) that contain data to be faceted/plotted
  #           ***data assumed to all be POSITIVE numbers
  #     samps, vector of names of sample columns in gcov
  #     datname, name of data in the sample columns for plotting - used as axis labels etc
  #     sname, name of samples category. e.g. Strain or Sample. Used for plotting - axis labels etc
  #     rowdescrip, description of data in each row. Used for plotting. 'Number of' <this> will be histogram y axis, for example.
  #     vline, F or where to draw a vertical line
  # Out: histogram
  
  # Reformat data to repeat for each sample, ggplot2-style
  pdata<-data.table()
  for(s in samps){
    onedata<-gcov[, c(names(gcov)[!names(gcov)%in%samps], s), with = F]
    setnames(onedata, s, "dat")
    onedata[,samp:=s]
    pdata<-rbind(pdata, onedata)
  }
  pdata[,samp:=factor(samp, levels = samps)]
  
  # Make hist
  h1<-ggplot(pdata, aes(dat)) + geom_histogram(bins = 100) +
    facet_grid(samp~.)+  xlab(datname) + ylab(paste("Number of", rowdescrip)) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          strip.text.y = element_text(size = 14))
  if(vline!=F){
    h1<-h1 + geom_vline(xintercept = vline, col = "blue", lty = "dashed")
  }
  
  return(h1)
}

# This was repurposed from other scripts
flexscat<-function(pdat, mycolors, colcol = "offInStrain", xcol = "rawcov", ycol = "pSegSites",
                   mytitle = "", mysubt = "", xlabel = "DNA sequencing coverage (raw x)",
                   ylabel = "Proportion segregating sites vs. N2",
                   labcolby = "", outlineby = NA, laboutlineby = NA, eqaxes = F, axmin = NA, axmax = NA,
                   annotlines = F){
  # Makes a scatter plot of DE vs ASE colored by regulatory pattern (more flexible than this in fact)
  # In: pdat, data with one row per gene/point to plot. Columns xcol value, ycol value, colcol value 
  #     mycolors, colors for bars: vector of length of categories in data; also used as ORDER for bars!! Names are categories, values are colors
  #         *also used to level factor AND SORT DATA
  #         Must have all values that are in colcol column of pdat
  #     colcol, name of column with data to color by
  #     xcol, character name of column for x axis values, typically DE
  #     ycol, character name of column for y axis values, typically ASE
  #     mytitle, title for plot
  #     mysubt, subtitle for plot
  #     xlabel, x axis label for plot
  #     ylabel, y axis label for plot
  #     labcolby, title of color legend label
  #     outlineby, optional column of res containing TRUE and FALSE as only values to have TRUE values outlined in black (NA to not separate points by outline)
  #     laboutlineby, What to label outline legend. If NA and outlineby provided, outlineby used
  #     eqaxes, F for auto-axes; T for make equal axes OR provide axis limits with next 2 options
  #     axmin, minimum value for both X and y axes. If NA, this will be minimum value of data in either xcol or ycol.
  #     axmax, maximum value for both x and y axes, if NA, this will be maximum value of data in either xcol or ycol
  #     annotlines, T or F: draw x = 0, y = 0 lines on plot?
  # Out: ggplot
  
  # Format data
  pdat<-copy(pdat) # don't want to mess with original data
  ## Make color col factor of desired order
  pdat[, colbycol:=factor(get(colcol), levels = names(mycolors))]
  pdat<-pdat[order(get(colcol))]
  
  # nexcl<-plt.data[,sum((get(xcol) < axmin | get(ycol) < axmin | 
  #                         get(xcol) > axmax | get(ycol) > axmax), na.rm = T)]
  
  # Make plot
  plt<-ggplot(pdat, aes(eval(as.name(xcol)), eval(as.name(ycol)))) + geom_point(aes(fill = colbycol), pch = 21, alpha = 0.4, stroke = 0.001) +
    scale_fill_manual(values = mycolors) + labs(fill = labcolby) +
    ggtitle(mytitle, subtitle = mysubt) + xlab(xlabel) + ylab(ylabel) +
    guides(fill = guide_legend(override.aes = list(size = 5))) +
    theme_bw() + theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14), 
                       axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
                       title = element_text(size = 17), legend.text = element_text(size = 13),
                       plot.subtitle = element_text(size = 15), strip.text.x = element_text(size = 15),
                       strip.text.y = element_text(size = 15),
                       legend.position = "bottom")
  
  if(!is.na(outlineby)){
    laboutlineby<-ifelse(is.na(laboutlineby), outlineby, laboutlineby)
    plt<-plt + geom_point(aes(stroke = 0.7, color = eval(as.name(outlineby))), pch = 21, alpha = 0.6, fill = NA) +
      scale_color_manual(values = c("TRUE" = 'black', "FALSE" = NA)) + labs(color = laboutlineby) +
      guides(color = guide_legend(override.aes = list(size = 5)))
  }
  
  if(annotlines){
    plt<-plt + geom_vline(xintercept = 0, col = "darkgray", lty = "dashed") + geom_hline(yintercept = 0, col = "darkgray", lty = "dashed")
  }
  
  if(eqaxes){
    ## Axis values
    if(is.na(axmin)){
      axmin<-min(c(pdat[, min(get(xcol), na.rm = T)], pdat[, min(get(ycol), na.rm = T)]))
    }
    if(is.na(axmax)){
      axmax<-max(c(pdat[, max(get(xcol), na.rm = T)], pdat[, max(get(ycol), na.rm = T)]))
    }
    plt<-plt + xlim(c(axmin, axmax)) + ylim(c(axmin, axmax)) 
  }
  
  return(plt)
}


#### Arguments & inputs ####
p<-arg_parser("Computes coverage *per gene* from the coverage of merged exons of that gene. Also: Plots coverage across genes for multiple samples; various analysis on this. On output of mosdepthgenesmergedexons.nf workflow processes", 
              name = "exploregenecoverage_fromexons.R", hide.opts = TRUE)
p<-add_argument(p, "--baseoutname",
                help = "Base name for all output files",
                type = "character",
                default = "out")
p<-add_argument(p, "--outdir",
                help = "Outer output directory. Sub-directories will be created internally.",
                type = "character")
p<-add_argument(p, "--strains",
                help = "Strains to process - in other inputs. Either comma-separated (no spaces) list or path to no-header file with one line per strain.
                Must match how strains are named in input files. Strains will be plotted/leveled in this order.",
                default = "N2,JU1088,EG4348,CB4856,QX1211")
p<-add_argument(p, "--isos",
                help = "Isotypes of the strains provided in --strains *in same order*. Either comma-separated (no spaces) list or path to no-header file with one line per strain.",
                default = "N2,JU1088,EG4349,CB4856,QX1211")
p<-add_argument(p, "--rawcov",
                help = "Path to file containing at least gene_id, raw coverage for each gene of interest in each strain of interest (i.e. _genecoveragefromexons_raw.txt.gz output of exploregenecoverage_fromexons.R)",
                type = "character")
p<-add_argument(p, "--offgenes",
                help = "Path to *_zeroexp_rnai_nocov_genes.txt file",
                type = "character")
p<-add_argument(p, "--nucdiv",
                help = "per-gene divergence metrics as output by nucdivcendr_geneswindows_allandasestrains.R ",
                type = "character")

# Parse arguments
cat("....Parsing arguments....\n")
p<-parse_args(p)

# Output directory
if(!dir.exists(p$outdir)){dir.create(p$outdir, recursive = T)}

# Strain info
if(file.exists(file.path(p$strains))){
  strains<-fread(file.path(p$strains), header = F)$V1
}else{
  strains<-strsplit(p$strains, split = ",", fixed = T)[[1]]
}
if(file.exists(file.path(p$isos))){
  isos<-fread(file.path(p$isos), header = F)$V1
}else{
  isos<-strsplit(p$isos, split = ",", fixed = T)[[1]]
}

# Coverage data format
samps<-isos
gcov<-fread(p$rawcov)
# Normalize to median gene coverage
gcov.mednorm<-data.table(gcov[, names(gcov)[!names(gcov)%in%samps], with = F], 
                         gcov[,lapply(samps, function(x) get(x)/median(get(x)))])
setnames(gcov.mednorm, c(names(gcov)[!names(gcov)%in%samps], samps))
# ## Organized data.table info so can loop/lapply through plots etc. largely copied form exploregenecoverage.R!
# dat.info<-data.table(datname = c("gcov", "gcov.mednorm"),
#                      fdescrip = c("rawcoverage", "medgenenormcoverage"),
#                      ldescrip = c("Raw coverage", "Coverage normalized to median gene's"))
# dat.list<-list(gcov, gcov.mednorm)
# names(dat.list)<-dat.info$datname

#### Just coverage ####
# Histograms just to 25% mean coverage etc
pdf(file.path(p$outdir, paste0(p$baseoutname, "_lowdnacovhists.pdf")), 6, 6)
covhistogram(gcov, samps) + xlim(c(0, 10)) + ylim(c(0, 50)) + # arbitrary/one off! thresholds
  ggtitle("Raw coverage only to 10x coverage")
covhistogram(gcov.mednorm, samps) + xlim(c(0, 0.25)) + ylim(c(0, 40)) + # arbitrary/one off! thresholds
  ggtitle("Med-norm coverage only to 25% median coverage")
invisible(dev.off())

# Number genes with actual 0 raw coverage; < 1% median; < 25% median 
## Convert to long from data. should've just made it longform data to begin with...
gcov.l<-data.table()
for(s in samps){
  onedata<-gcov[, c(names(gcov)[!names(gcov)%in%samps], s), with = F]
  setnames(onedata, s, "rawcov")
  onedata[,samp:=s]
  gcov.l<-rbind(gcov.l, onedata)
}
gcov.l[,samp:=factor(samp, levels = samps)]
setkey(gcov.l, samp)

gcovmed.l<-data.table()
for(s in samps){
  onedata<-gcov.mednorm[, c(names(gcov.mednorm)[!names(gcov.mednorm)%in%samps], s), with = F]
  setnames(onedata, s, "mednormcov")
  onedata[,samp:=s]
  gcovmed.l<-rbind(gcovmed.l, onedata)
}
gcovmed.l[,samp:=factor(samp, levels = samps)]
setkey(gcovmed.l, samp)

# Get #s
numsumm<-gcovmed.l[, .(.N, sum(mednormcov<0.25)), by = samp][gcov.l[ ,sum(rawcov==0), by = samp]]
setnames(numsumm, c("strain", "ngenes", "nUnder25pMedDNACov", "n0RawCov"))
numsumm[,pAllThatAre0Cov:=n0RawCov/ngenes]
numsumm[,pLowThatAre0Cov:=n0RawCov/nUnder25pMedDNACov]
write.table(numsumm, file.path(p$outdir, paste0(p$baseoutname, "_number0coveragesummary.txt")),
            sep = "\t", quote = F, row.names = F)

# Bar plot # low coverage with 0 coverage number highlighted
bpdat<-rbind(numsumm[,.(strain, group = "0x raw coverage", num = n0RawCov)],
             numsumm[,.(strain, group = ">0x raw coverage,\n<25% med. coverage", num = nUnder25pMedDNACov - n0RawCov)])
setnames(bpdat, "strain", "Strain")

pdf(file.path(p$outdir, paste0(p$baseoutname, "_barplot0covlowcov.pdf")), 5, 5)
ggplot(bpdat, aes(Strain, num)) + geom_bar(aes(fill = group), stat = "identity") +
  scale_fill_manual(values = c("gray30", "gray"), breaks =  c("0x raw coverage", ">0x raw coverage,\n<25% med. coverage")) +
  ylab("# genes with low DNA coverage\nin strain (< 25% of median)") + labs(fill = "") +
  theme_bw() + theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14, angle = 45, vjust = 0.55), 
                     axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
                     strip.text.y = element_text(size = 14), legend.position = "bottom",
                     legend.text = element_text(size = 12), legend.title = element_text(size = 13))
invisible(dev.off())

#### Coverage and off/DE genes ####
offgs<-fread(p$offgenes)
# Make long format, add in actual coverage values [backsolving to data I used to have...]
setkey(gcov.l, gene_id) 
setkey(gcovmed.l, gene_id)

offg.l<-rbindlist(lapply(1:length(strains), function(strnind){
  strn<-strains[strnind]
  iso<-isos[strnind]
  
  out<-offgs[grepl(strn, whichStrainsZeroExp), .(gene_id, get(paste0(strn, ".nocov")))]
  setnames(out, "V2", "nocovflagged")
  out[,Strain:=strn]
  setcolorder(out, "Strain")
  setkey(out, gene_id)
  
  # Add coverages
  out<-gcov.l[samp==iso, .(gene_id, locus, sequence_name, rawcov)][gcovmed.l[samp==iso, .(gene_id, mednormcov)]][out]
  setcolorder(out, c("Strain", "gene_id", "sequence_name", "nocovflagged", "mednormcov", "rawcov"))
  
  return(out)
}))
offg.l[, Strain := factor(Strain, levels = strains)] # each gene is off in the strain in Strain col! [and optionally other strains too]

# Summary numbers: off per strain and coverage category
offgcovnsumm<-offg.l[, .(nOffGenes = .N, nUnder25pMedDNACov = sum(nocovflagged == T), 
                         n0RawCov = sum(rawcov==0),
                         nGreater0Under25pCov = sum(rawcov > 0 & mednormcov < 0.25),
                         nOver25pMedDNACov = sum(nocovflagged == F)),
                     by = Strain]
offgcovnsumm[, `:=`(poff_0RawCov = n0RawCov/nOffGenes, 
                    poff_Greater0Under25pCov = nGreater0Under25pCov/nOffGenes,
                    p_offOver25pMedDNACov = nOver25pMedDNACov/nOffGenes)]

write.table(offgcovnsumm, file.path(p$outdir, paste0(p$baseoutname, "_offgenesX0coveragesummary.txt")),
            sep = "\t", quote = F, row.names = F)
# Barplot
bpoffgdat<-rbind(offgcovnsumm[,.(Strain, group = "zero", num = n0RawCov)],
                 offgcovnsumm[, .(Strain, group = "low", num = nGreater0Under25pCov)],
                 offgcovnsumm[, .(Strain, group = "normal", num = nOver25pMedDNACov)])
bpoffgdat[, group:=factor(group, levels = c("zero", "low", "normal"))]

pdf(file.path(p$outdir, paste0(p$baseoutname, "_barplotOffGenes0covlowcov.pdf")), 5, 5)
ggplot(bpoffgdat, aes(Strain, num)) + geom_bar(aes(fill = group), stat = "identity") +
  scale_fill_manual(values = c("gray30", "gray60", "gray80"), breaks =  c("zero", "low", "normal")) +
  ylab("# genes turned off in strain") + labs(fill = "DNA seq. coverage") +
  theme_bw() + theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 14, angle = 45, vjust = 0.55), 
                     axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
                     strip.text.y = element_text(size = 14), legend.position = "bottom",
                     legend.text = element_text(size = 12), legend.title = element_text(size = 13))
invisible(dev.off())

#### Save session info ####
cat("....exploregenecoverage_fromexons_lowend.R complete! Session information:....\n")
sessionInfo()
