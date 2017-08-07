#!/bin/sh
#$ -S /bin/bash
#$ -pe smp 4
#$ -cwd
#$ -N hla_typing
# number of threads
#NT=4


SAMPLE=$1
mkdir -p "hla_optitype"
cd "hla_optitype"

exec 1>log_hla.optitype.$SAMPLE.out
exec 2>log_hla.optitype.$SAMPLE.err

mkfifo fq1
mkfifo fq2
ln -s ../normal.bam normal.bam

picard SortSam I="normal.bam" O="/dev/stdout" SORT_ORDER=queryname | \
bedtools bamtofastq -i "/dev/stdin" -fq "fq1" -fq2 "fq2" &
#picard SamToFastq I="/dev/stdin" F="fq1" F2="fq2" &
bwa mem -M -t 2 "$HLA_DNA_REFERENCE" "fq1" | samtools view -bSho "out1.bam" - &
bwa mem -M -t 2 "$HLA_DNA_REFERENCE" "fq2" | samtools view -bSho "out2.bam" - &
wait;

rm normal.bam
rm *.bai

samtools bam2fq out1.bam > out1.fastq
samtools bam2fq out2.bam > out2.fastq

rm out1.bam
rm out2.bam

cp "$HLA_CONFIG" .
optitype --input out1.fastq out2.fastq --dna -c config.ini --outdir hla

rm *.fastq

mv hla/*/*.tsv $SAMPLE.HLA.tsv
mv hla/*/*.pdf $SAMPLE.HLA.pdf

rm -r hla

cd ..
