#!/usr/bin/env nextflow
/*
By Avery Davis Bell
Creation of strain-specific transcriptomes and indices
*/

/*
#### Set up input parameters & defaults ####
*/
// Inputs: input/output related
params.strainlist = "" // Path to one-column file with strains to process (one per line). These must be column headers in VCFs and will be used for output labelling as well.
params.outputdir = "" // parent output directory
params.snpvcf = "/storage/coda1/p-apaaby3/0/shared/datasets/cendr-data/vcfs-20210121/WI.20210121.hard-filter.isotype.snpsonly.vcf.gz" // VCF containing all SNPs (only) for all strains of interest vs. reference genome of interest
params.indelvcf = "/storage/coda1/p-apaaby3/0/shared/datasets/cendr-data/vcfs-20210121/WI.20210121.hard-filter.isotype.indelsonly.vcf.gz" // VCF containing all INDELs (only) for all strains of interest vs. reference genome of interest
params.reffasta = "/storage/coda1/p-apaaby3/0/shared/referenceresources/celegans_n2_ws276/c_elegans.PRJNA13758.WS276.genomic.fa" // FASTA reference genome - for strain-specific genome creation; what variants were called against. Should be unzipped
params.refgtf = "/storage/coda1/p-apaaby3/0/shared/labmembers/abell65/rnaomnibus/genecountquant/emase/ws276_20210121cendr/genomeinfo/c_elegans.PRJNA13758.WS276.canonical_geneset_nounderscoretranscriptids.gtf" // GTF file for reference genome.  Transcripts must NOT have underscores in their name if salmon is to be run through EMASE. Modify transcript names to exclude them if necessary.
params.idx = 'salmon' // Which index to build? Possible values: 'bowtie2', 'salmon', 'all' (builds both), 'none' (builds none)
params.salmdecoy = 'no' // 'yes' or 'no' - add REFERENCE genome sequence as decoy sequence for salmon index building?

// Inputs: organizational
params.g2gtoolsconda = '/storage/home/hcoda1/2/abell65/.conda/envs/g2gtools' // path to python 2 conda environment where g2gtools properly set up
params.procgtfscriptdir = '/storage/coda1/p-apaaby3/0/shared/labmembers/abell65/scripts/rnaseq/ase/processdiptranscriptomes' // Directory containing Python GTF worker scripts gtfnonchroverlaps.py and gtftoemasegenemapping.py
params.salmonenv = '/storage/home/hcoda1/2/abell65/.conda/envs/salmon' // path to conda environment where salmon is installed

// Housekeeping:  create output directories
outdir = file(params.outputdir)
outdir.mkdirs()

/*
#### Channel management ####
*/
straininfo = Channel.fromList(file(params.strainlist).readLines())

/*
#### Processes ####
*/

process g2gvcf2chain{
  // Chain indels onto reference

  conda params.g2gtoolsconda

  input:
  val(mystrain) from straininfo

  output:
  tuple val(mystrain), path("*.chain") into indelchains

  """
  g2gtools vcf2chain -f ${params.reffasta} \
  -i ${params.indelvcf} \
  -s ${mystrain} \
  -o REF-to-${mystrain}.chain
  """
}

process g2gpatch{
  // Patch SNPs onto reference genome

  conda params.g2gtoolsconda

  input:
  tuple val(mystrain), path(chain) from indelchains

  output:
  tuple val(mystrain), path(chain), path("*.fa") into snppatchedfas

  """
  g2gtools patch -i ${params.reffasta} \
   -s ${mystrain} \
   -v ${params.snpvcf} \
   -o ${mystrain}.patched.fa
  """
}

process g2gtransform{
  // chain indels onto patched genome

  conda params.g2gtoolsconda

  input:
  tuple val(mystrain), path(chain), path(fasta) from snppatchedfas

  output:
  tuple val(mystrain), path(chain), path("*_final.fa") into strainfa

  """
  g2gtools transform -i ${fasta} \
  -c ${chain} \
  -o ${mystrain}_final.fa
  """
}

process g2gconvert{
  // update GTF file based on new genome

  conda params.g2gtoolsconda

  input:
  tuple val(mystrain), path(chain), path(strainfasta) from strainfa

  output:
  tuple val(mystrain), path(strainfasta), path("*.gtf"), path("*.gtf.unmapped") into gtfplus

  """
  g2gtools convert -c ${chain} \
  -i ${params.refgtf} \
  -f gtf \
  -o ${mystrain}.gtf

  if [ ! -e ${mystrain}.gtf.unmapped ]; then
    touch ${mystrain}.gtf.unmapped
  fi
  """
}

process getexclseqs{
  // For alt strain, splits GTF into sequences that should be excluded vs. included from making transcriptome (those that end after chromosome ends are excluded)
  // Runs gtfnonchroverlaps.py

  input:
  tuple val(mystrain), path(strainfasta), path(straingtf), path(strainunmappedgtf) from gtfplus

  output:
  tuple val(mystrain), path(strainfasta), path(strainunmappedgtf), path("*keep.gtf"), path("*excl.gtf") into usegtfs

  """
  bioawk -c fastx '{print \$name "\\t" length(\$seq)}' ${strainfasta} > ${mystrain}_genome_lengths.txt

  python ${params.procgtfscriptdir}/gtfnonchroverlaps.py -gtf ${straingtf} \
  -lens ${mystrain}_genome_lengths.txt \
  -okeep ${mystrain}_keep.gtf \
  -oexcl ${mystrain}_excl.gtf
  """
}

process straintranscriptome{
  // Generate strain-specific transcriptome
  // A couple processes together: gffread to extract, then some formatting

  publishDir "${params.outputdir}/${mystrain}", mode: 'copy', path: "*namesorted.fa", overwrite: true

  input:
  tuple val(mystrain), path(strainfasta), path(strainunmappedgtf), path(keepgtf), path(exclgtf) from usegtfs

  output:
  tuple val(mystrain), path(strainunmappedgtf), path(exclgtf), path("*namesorted.fa") into strtrfa

  """
  gffread -v -w ${mystrain}_transcriptome.fa -g ${strainfasta} ${keepgtf}

  bioawk -c fastx '{print ">" \$name ORS \$seq}' ${mystrain}_transcriptome.fa | \
  seqkit sort -n -2 > ${mystrain}_transcriptome_namesorted.fa
  """
}

process straintrnslengths{
  // Get lengths of each transcript for one strain, adding in any that were excluded

  publishDir "${params.outputdir}/${mystrain}", mode: 'copy', path: "${mystrain}_transcript_lengths.txt", overwrite: true

  input:
  tuple val(mystrain), path(strainunmappedgtf), path(exclgtf), path(straintrfa) from strtrfa

  output:
  tuple val(mystrain), path(straintrfa), path("*lengths.txt") into strtrfalns
  tuple val(mystrain), path(straintrfa), path("*lengths.txt") into strtrfalns2
  tuple val(mystrain), path(straintrfa), path("*lengths.txt") into strtrfalns3

  """
  bioawk -c fastx '{print \$name "\\t" length(\$seq)}' ${straintrfa} \
  > ${mystrain}_transcript_lengths.tmp

  nadd=\$(cat $strainunmappedgtf $exclgtf | awk '\$3=="transcript"'| wc -l)
  if [[ \$nadd -gt 0 ]]; then
    cat ${strainunmappedgtf} ${exclgtf} | awk '\$3=="transcript"' | cut -f9 | \
    cut -d\\; -f2 | awk '{print \$2}' |  tr -d '"' | awk '{print \$0}' | \
    cat ${mystrain}_transcript_lengths.tmp - | sort -k1 > ${mystrain}_transcript_lengths.txt
  else
    mv ${mystrain}_transcript_lengths.tmp ${mystrain}_transcript_lengths.txt
  fi
  """
}

process bowtie2idx{
  // Builds bowtie2 single end index for diploid transcriptome
  // Only runs if --idx is 'bowtie2' or 'all'

  // save outputs
  publishDir "${params.outputdir}/${mystrain}", mode: 'copy', path: "*"

  input:
  tuple val(mystrain), path(straintrfa), path(ln) from strtrfalns

  output:
  path("*") into bowtieidxs

  when:
  params.idx == 'bowtie2' || params.idx == 'all'

  """
  bowtie2-build -f --offrate 4 --threads 4 ${straintrfa} ${mystrain}
  """
}

process salmondecoyprep{
  // Creates fasta of diploid transcriptome + reference genome (as decoy), + chromosome names for use as salmon index's decoy
  // Only runs if --idx is 'salmon' or 'all' AND --salmdecoy is 'yes'

  input:
  tuple val(mystrain), path(straintrfa), path(ln) from strtrfalns2

  output:
  tuple val(mystrain), path("diptranscrrefgen.fa"), path("decoys.txt") into salmonin

  when:
  params.salmdecoy == 'yes' && (params.idx == 'salmon' || params.idx == 'all')

  """
  grep "^>" <(gunzip -c ${params.reffasta}) | cut -d " " -f 1 > decoys.txt
  sed -i.bak -e 's/>//g' decoys.txt
  zcat ${params.reffasta} | cat ${straintrfa} - > diptranscrrefgen.fa
  """
}

process salmonidxdecoy{
  // builds salmon index (using ref genome as decoy)
  // Only runs if --idx is 'salmon' or 'all' AND --salmdecoy is 'yes'

  conda params.salmonenv

  // save outputs
  publishDir "${params.outputdir}/${mystrain}", mode: 'copy', path: "*", overwrite: true

  input:
  tuple val(mystrain), path(trgfa), path(decoys) from salmonin

  output:
  path("*") into salmonidxwithdec

  when:
  params.salmdecoy == 'yes' && (params.idx == 'salmon' || params.idx == 'all')

  """
  salmon index \
  -t ${trgfa} \
  -i salmon_idx \
  --decoys ${decoys} \
  -k 31 --keepDuplicates \
  -p 4
  """
}

process salmonidxnodecoy{
  // builds salmon index with NO DECOY
  // Only runs if --idx is 'salmon' or 'all' AND --salmdecoy is 'no'

  conda params.salmonenv

  // save outputs
  publishDir "${params.outputdir}/${mystrain}", mode: 'copy', path: "*", overwrite: true

  input:
  tuple val(mystrain), path(straintrfa), path(ln) from strtrfalns3

  output:
  path("*") into salmonidxnodec

  when:
  params.salmdecoy == 'no' && (params.idx == 'salmon' || params.idx == 'all')

  """
  salmon index \
  -t ${straintrfa} \
  -i salmon_idx \
  -k 31 --keepDuplicates \
  -p 4
  """
}
