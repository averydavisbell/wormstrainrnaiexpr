#!/usr/bin/env nextflow
/*
By Avery Davis Bell
Runs mosdepth for *merged exons* on BAMs - creates mered bed file too.
*/

/*
#### Set up input parameters & defaults ####
*/
// General parameters
params.gtf = "" // path to GTF containing genes for which to determine coverage from all BAMs. Columns as ws276 GTF from Wormbase.
params.sampinfo = "" // Path to sample information file. Columns SampleID (name of sample for output), bam (path to BAM file to process for this sample), bai (path to BAM .bai index file for this sample)
params.outdir = "out" // Path to output directory
params.outname = "out" // Prefix for output files that contain all samples' mosdepth information
// GTFtools parameters
params.refname = "ref" // Prefix for output file containing merged exons - i.e. reference genome name
params.chrs = "I,II,III,IV,V,X,MtDNA" // Chromosomes to process exons/genes for - need to match the GTF. Default is for C. elegans ws276
params.gtftoolsdir = "GTFtools_0.8.5" // Directory containing GTF tools python script gtftools.py
// mosdepth parameters
params.flag = "1796" // --flag (SAM flag bits to exclude) argument for mosdepth
params.mapq = "0" // -Q, mapq threshold argument for mosdepth, threshold below which read will be excluded
// Summary R script parameters
params.rscriptdir = "../.." // Directory containing exploregenecoverage_fromexons.R
params.gff = "" // Path to *genes only* gff3 file containing info on all genes from the GTF
params.gsubset = '' // OPTIONAL Path to no-header list of genes to run summary R script for - SUBSET of all genes. It will also be run for all genes.
params.gsubsetname = "" // OPTIONAL name of gene subset for output filenaming. Provide if provide --gsubset

// Create output directories
outdir = file(params.outdir)
outdir.mkdirs()
outplotdir = file(params.outdir + "/plots")
outplotdir.mkdirs()

/*
#### Channel management ####
*/

//  Channel with all needed sample info as tuple.
formatChannel = Channel.fromPath(params.sampleinfo)
  .splitCsv(header: true, sep: '\t', strip: true)
  .map{row->
    return [row.SampleID, row.bam, row.bai]
  }
  .set{sampleInfo}

/*
#### Processes ####
*/

process gtf2mergedexonbed{
  // Get merged exons bed file from GTF using GTFtools
  publishDir outdir, mode: 'copy', pattern: "*.bed.gz"

  output:
  path("*.bed") into bedfile
  path("*.bed.gz") into gzout

  """
  python ${params.gtftoolsdir}/gtftools.py -m ${params.refname}.mergedexons.bed \
  -c ${params.chrs} ${params.gtf}

  cp *.bed tmp.bed
  gzip tmp.bed
  mv tmp.bed.gz ${params.refname}.mergedexons.bed.gz
  """
}

process mosdepth{
  // Run mosdepth for genes. Also unzips bed output for downstream ease.

  input:
  path(bed) from bedfile
  tuple val(mysamp), path(bam), path(bai) from sampleInfo

  output:
  path("*.regions.bed") into bedcovs
  path("*withsampcol.mosdepth.summary.txt") into covsumms

  """
  mosdepth -t 4 \
  -b ${bed} \
  -n \
  -F ${params.flag} \
  -Q ${params.mapq} \
  ${mysamp} ${bam}

  # Output management for downstream
  gunzip *.bed.gz
  awk '{print \"${mysamp}\\t\" \$0}' ${mysamp}.mosdepth.summary.txt > ${mysamp}withsampcol.mosdepth.summary.txt
  """
}

process combinedpbeds{
  // Combines *.regions.bed.gz mosdepth outputs into one file with all samples

  publishDir outdir, mode: 'copy', overwrite: true

  input:
  file bed from bedcovs.collect()

  output:
  path("*.bed.gz") into allbed
  path("*.bed.gz") into allbed2

  """
  # Organize sample file list & get sampleIDs
  echo ${bed} | sed 's/\\s\\+/\\n/g' | sort > tmp.txt
  awk -F\".regions.bed\" '{print \$1}' tmp.txt | paste -s > samps.txt
  echo "chr\tstart\tend\tname" | paste - samps.txt > header.txt

  # Make file with gene info, one column per sample
  cat tmp.txt | xargs paste > tmpall.txt
  awk -F\"\\t\" '{for(i=5;i<=NF;i+=5)printf \"%s%s\", \$i, (i+3>NF?\"\\n\":FS)}' tmpall.txt > tmpnums.txt
  cut -f1-4 tmpall.txt | paste - tmpnums.txt > bed.tmp

  # Add header
  cat header.txt bed.tmp > ${params.outname}.mergedexons.bed

  # gzip
  gzip ${params.outname}.mergedexons.bed
  """
}

process combinedpsumms{
  // Combines *.mosdepth.summary.txt mosdepth outputs into one file with all samples
  publishDir outdir, mode: 'copy', overwrite: true

  input:
  file summ from covsumms.collect()

  output:
  path("*.mosdepth.summary.txt") into allsumm
  path("*.mosdepth.summary.txt") into allsumm2

  """
  # Combine summaries without their headers
  echo ${summ} | sed 's/\\s\\+/\\n/g' | sort > tmp.txt
  cat tmp.txt | xargs tail -n+2 -q > allsumm.tmp

  # Make new header and add in
  echo "SampleID" > sampid.tmp
  head -1 tmp.txt | xargs head -1 | cut -f2- | paste sampid.tmp - > head.tmp
  cat head.tmp allsumm.tmp > ${params.outname}.mosdepth.summary.txt
  """
}

process comboexonsexplore{
  // Runs exploregenecoverage_fromexons.R for all genes

  publishDir outdir, mode: 'copy', overwrite: true, pattern: "*.txt*"
  publishDir outplotdir, mode: 'copy', overwrite: true, pattern: "*.pdf"

  input:
  file exbed from allbed
  file summ from allsumm

  output:
  path("*") into routs

  """
  Rscript ${params.rscriptdir}/exploregenecoverage_fromexons.R \
  --exoncov ${exbed} \
  --covsumm ${summ} \
  --genegff ${params.gff} \
  --outstem ${params.outname}
  """
}

process comboexonsexploresubset{
  // Runs exploregenecoverage_fromexons.R for all genes

  publishDir outdir, mode: 'copy', overwrite: true, pattern: "*.txt*"
  publishDir outplotdir, mode: 'copy', overwrite: true, pattern: "*.pdf"

  input:
  file exbed from allbed2
  file summ from allsumm2

  output:
  path("*") into routs2

  when:
  params.gsubset != ''

  """
  Rscript ${params.rscriptdir}/exploregenecoverage_fromexons.R \
  --exoncov ${exbed} \
  --covsumm ${summ} \
  --genegff ${params.gff} \
  --outstem ${params.outname}_${params.gsubsetname} \
  --genelist ${params.gsubset}
  """
}
