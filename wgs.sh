#!/bin/bash

soft/bowtie2-build databases/all_genome.fasta databases/all_genome.fasta

mkdir out
mkdir out/sam_data
mkdir out/bam_data
mkdir out/gene_data

soft/bowtie2 --fast --end-to-end -x databases/all_genome.fasta -p 6 -1 fastq/EchE_R1.fastq.gz -2 fastq/EchE_R2.fastq.gz -S out/sam_data/EchE.sam

samtools view -b --threads 6 -1 out/sam_data/EchE.sam -o out/bam_data/EchE.bam

samtools sort --threads 6 out/bam_data/EchE.bam -o out/bam_data/EchE.sorted.bam

samtools index out/bam_data/EchE.sorted.bam

samtools idxstats --threads 6 out/bam_data/EchE.sorted.bam > out/bam_data/EchE.idxstats 

soft/megahit --mem-flag 0 --k-list 21 -1 fastq/EchE_R1.fastq.gz -2 fastq/EchE_R2.fastq.gz -o out/megahit_data/

prodigal -i out/megahit_data/final.contigs.fa -d out/gene_data/pred_genes.fasta

sed "s:>:*\n>:g" out/gene_data/pred_genes.fasta | sed -n "/partial=00/,/*/p"|grep -v "*" > out/gene_data/genes_full.fna

blastn -db databases/resfinder.fna -query out/gene_data/genes_full.fna -perc_identity 80 -evalue 0.001 -qcov_hsp_perc 80 -outfmt '6 qseqid sseqid pident qcovs evalue' -best_hit_score_edge 0.001
