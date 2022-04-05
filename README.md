# Longshort-analysis pipeline
## Preprocessing with long read data
### Convert existing bam file to fastq file
bedtools bamtofastq -i m54189_170915_112634.subreads.bam -fq m1.subread.fastq

### Align m1 long reads to full genome (GRCh38)
minimap2 -ax splice -uf -C5 GRCh38.p13.genome.fa m1.subread.fastq | samtools sort -o m1.genome.align.bam

### Use stringtie to get gtf file for longread-genome alignment with v39 annotation
stringtie -L -o m1.v39.genome.gtf -G gencode.v39.annotation.gtf m1.genome.align.bam

### use gffread to get corresponding transcript fasta file for long reads, which is going to be used as long read library
gffread m1.v39.genome.gtf  -g GRCh38.p13.genome.fa -w m1.v39.0322.fasta

## Kallisto quantification
### generate index for full transcriptome
kallisto index -i v39fullindex gencode.v39.transcript.fa

### full transcriptome-short reads analysis
kallisto quant -i v39fullindex -o /storage/yhhuang/users/yhsz/2022_longread/analysis0321/onestep --pseudobam GX069-AO_1.fastq.gz GX069-AO_2.fastq.gz

### generate index for long read
kallisto index -i m1v39index m1.v39.0322.fasta

### long-short reads analysis
kallisto quant -i m1v39index -o /storage/yhhuang/users/yhsz/2022_longread/analysis0321/longshort --pseudobam GX069-AO_1.fastq.gz GX069-AO_2.fastq.gz

## Merging GTF files
### Trail one: Gffcompare
gffcompare -r gencode.v39.annotation.gtf -R m1_subreads_stgt2.gtf

### Trail two: StringTie merging mode, guided
Use Stringtie to merge predicted transcripts from all libraries into a unified transcriptome.

stringtie --rf --merge -p 4 -o stringtie_merged.gtf -G gencode.v39.annotation.gtf m1_subreads_stgt2.gtf

Compare reference guided transcripts to the known annotations

gffcompare -r m1_subreads_stgt2.gtf -o gffcompare stringtie_merged.gtf

### Trail three: StringTie merging mode, de novo
stringtie --rf --merge -p 4 -o denovo_full_m1_merge.gtf gencode.v39.annotation.gtf m1_subreads_stgt2.gtf 

gffcompare -r m1_subreads_stgt2.gtf -o gffcompare denovo_full_m1_merge.gtf 

## Using merged GTF as guidance to assign estimated/raw counts for each transcript

estimated counts are obtained in abundance file generated with Kallisto

raw count can be obtained with following command lines:

samtools idxstats sortedonestep.bam > idxstat.txt

awk '{print $1, $3}' longshortidxstat.txt > idxstatmapped.txt  (third column is mapped raw count)

These data are processed to csv files, which can be easily read in pandas

## The jupyter lab book for processing data to get scatter plot can be found in attached file
