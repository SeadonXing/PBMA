#!/bin/bash

export PATH=/Software/anaconda3/envs/Numbat/bin:$PATH

datdir=/Project/BrM/CellRanger_count
refdir=/Project/BrM/Numbat/Ref_Software

samples=(BrM02 BrM03 BrM04 BrM05 BrM06 BrM07 BrM08 BrM09 BrM10) # Change to all SampleID

for pat in ${samples[@]}
do
	cp ${datdir}/${pat}_count/outs/possorted_genome_bam.bam ./${pat}_possorted_genome_bam.bam
	cp ${datdir}/${pat}_count/outs/possorted_genome_bam.bam.bai ./${pat}_possorted_genome_bam.bam.bai
	cp ${datdir}/${pat}_count/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ./${pat}_barcodes.tsv.gz

Rscript pileup_and_phase.R \
	--label ${pat} \
	--samples ${pat} \
	--bams ${pat}_possorted_genome_bam.bam \
	--barcodes ${pat}_barcodes.tsv.gz \
	--gmap ${refdir}/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
	--eagle ${refdir}/Eagle_v2.4.1/eagle \
	--snpvcf ${refdir}/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
	--paneldir ${refdir}/1000G_hg38 \
	--outdir ./${pat} \
	--ncores 32

  rm ${pat}_possorted_genome_bam.*
  rm ${pat}_barcodes.tsv.gz

done
