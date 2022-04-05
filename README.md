# Longshort-analysis
## Preprocessing with long read data:
### Convert existing bam file to fastq file
bedtools bamtofastq -i m54189_170915_112634.subreads.bam -fq m1.subread.fastq
### Align m1 long reads to full genome (GRCh38)
minimap2 -ax splice -uf -C5 GRCh38.p13.genome.fa m1.subread.fastq | samtools sort -o m1.genome.align.bam
### Use stringtie to get gtf file for longread-genome alignment with v39 annotation
stringtie -L -o m1.v39.genome.gtf -G gencode.v39.annotation.gtf m1.genome.align.bam
### use gffread to get corresponding transcript fasta file for long reads, which is going to be used as long read library
gffread m1.v39.genome.gtf  -g GRCh38.p13.genome.fa -w m1.v39.0322.fasta
