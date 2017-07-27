#!/bin/sh
#$ -S /bin/bash
#$ -pe smp 2
#$ -l mem=10G
#$ -cwd
#$ -o log/log.o
#$ -e log/log.e

. config.sh
. /etc/profile.d/sge.sh

#exec 1>log_dna.out

sample=$(awk NR==$SGE_TASK_ID filelist.txt)
read IF_GERMLINE IF_SOMATIC SAMPLE <<< $sample

mkdir -p "$SAMPLE"
cd "$SAMPLE"

exec 2>log_dna_realignment.$SAMPLE.err
exec 1>log_dna_realignment.$SAMPLE.out

cp "$BAM_ROOT/$IF_GERMLINE" normal.bam
INDEX="$IF_GERMLINE"."bai"
cp "$BAM_ROOT/$INDEX" normal.bai
cp "$BAM_ROOT/$IF_SOMATIC" tumor.bam
INDEX="$IF_SOMATIC"."bai"
cp "$BAM_ROOT/$INDEX" tumor.bai

#### process .bam file
for i in *.bam
do
id=$(echo $i | awk -F'[.]' '{print $1}')
  picard RevertSam I=$id.bam \
                    O=$id.tmp.ubam \
                    SANITIZE=true \
                    ATTRIBUTE_TO_CLEAR=XT \
                    ATTRIBUTE_TO_CLEAR=XN \
                    ATTRIBUTE_TO_CLEAR=AS \
                    ATTRIBUTE_TO_CLEAR=OC \
                    ATTRIBUTE_TO_CLEAR=OP \
                    SORT_ORDER=queryname \
                    RESTORE_ORIGINAL_QUALITIES=true \
                    REMOVE_DUPLICATE_INFORMATION=true \
                    REMOVE_ALIGNMENT_INFORMATION=true \
                    TMP_DIR=tmp
 picard AddOrReplaceReadGroups \
 I=$id.tmp.ubam \
 O=$id.ubam \
 RGID=$RG \
 RGLB=$LN \
 RGPL=$P \
 RGPU=$PU \
 RGSM=$id \
 TMP_DIR=tmp

rm $id.tmp.ubam
rm $id.bam
rm $id.bai
done

# Mark adapter sequences using MarkIlluminaAdapters 
for i in *.ubam
do
        id=$(echo $i | awk -F'[.]' '{print $1}')
        picard MarkIlluminaAdapters I=$id.ubam \
                                    O=$id.markadapt.ubam \
                                    M=$id.markadapt.metrics.txt \
                                    TMP_DIR=tmp
rm $id.ubam
done
## Alignment by BWA MEM
set -o pipefail
for i in *.markadapt.ubam
do
	id=$(echo $i | awk -F'[.]' '{print $1}')
picard SamToFastq I="$id.markadapt.ubam" FASTQ="/dev/stdout" CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true TMP_DIR=tmp | bwa mem -M -t 2 -p "$REF_FASTA" "/dev/stdin" | \
samtools view -bSho $id.out.bam - ;
picard MergeBamAlignment ALIGNED_BAM="$id.out.bam" UNMAPPED_BAM="$id.markadapt.ubam" OUTPUT="$id.bam" R="$REF_FASTA" \
	CREATE_INDEX=true \
	ADD_MATE_CIGAR=true \
	CLIP_ADAPTERS=false \
	CLIP_OVERLAPPING_READS=true \
	INCLUDE_SECONDARY_ALIGNMENTS=true \
	MAX_INSERTIONS_OR_DELETIONS=-1 \
	PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
	ATTRIBUTES_TO_RETAIN=XS \
	TMP_DIR=tmp
rm $id.markadapt.ubam
rm $id.out.bam
done

########################### Mark duplicates #################################################
for i in *.bam
do
                id=$(echo $i | awk -F'[.]' '{print $1}')
                picard MarkDuplicatesWithMateCigar \
                        INPUT=$id.bam \
                        OUTPUT=$id.markdupl.bam \
                        METRICS_FILE=$id.markdupl.metrics.txt \
                        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
                        CREATE_INDEX=true \
                        TMP_DIR=tmp
rm $id.bam
rm $id.bai
done

######################## local realignment around indels #####################################
for i in *.markdupl.bam
do
        id=$(echo $i | awk -F'[.]' '{print $1}')
        gatk \
                -T RealignerTargetCreator \
                -R $REF_FASTA \
                -known $GOLD_INDEL \
                -known $KNOWN_INDEL \
                -I $id.markdupl.bam \
                -o $id.realignertargetcreator.intervals
done
#### Realign reads using IndelRealigner
for i in *.markdupl.bam
do
        id=$(echo $i | awk -F'[.]' '{print $1}')
        gatk \
                -T IndelRealigner \
                -R $REF_FASTA \
                -targetIntervals $id.realignertargetcreator.intervals \
                -known $KNOWN_INDEL \
                -known $GOLD_INDEL \
                -I $id.markdupl.bam \
                -o $id.indelrealign.bam
rm $id.markdupl.bam
rm $id.markdupl.bai
done

##  Recalibrate base quality scores = run BQSR    
## Analyze patterns of covariation in the sequence dataset
for i in *.indelrealign.bam
do
        id=$(echo $i | awk -F'[.]' '{print $1}')
        gatk \
                -T BaseRecalibrator \
                -R $REF_FASTA \
                -I $id.indelrealign.bam \
                -knownSites $KNOWN_INDEL \
                -knownSites $GOLD_INDEL \
                -knownSites $PHASE1_SNP \
		-knownSites $HAPMAP \
                -knownSites $DBSNP \
		-knownSites $OMNI \
                -o $id.recal_data.table
done

#### Do a second pass to analyze covariation remaining after recalibration
for i in *.indelrealign.bam
do
        id=$(echo $i | awk -F'[.]' '{print $1}')
        gatk \
                -T BaseRecalibrator \
                -R $REF_FASTA \
                -I $id.indelrealign.bam \
                -knownSites $KNOWN_INDEL \
                -knownSites $GOLD_INDEL \
                -knownSites $PHASE1_SNP \
                -knownSites $HAPMAP \
                -knownSites $DBSNP \
                -knownSites $OMNI \
		-BQSR $id.recal_data.table \
                -o $id.post_recal_data.table
done

for i in *.indelrealign.bam
do
        id=$(echo $i | awk -F'[.]' '{print $1}')
        gatk \
                -T AnalyzeCovariates \
                -R $REF_FASTA \
                -before $id.recal_data.table \
                -after $id.post_recal_data.table \
                -plots $id.recalibration_plots.pdf
done
####  Apply the recalibration to the sequence data
for i in *.indelrealign.bam
do
        id=$(echo $i | awk -F'[.]' '{print $1}')
        gatk \
                -T PrintReads \
                -R $REF_FASTA \
                -I $id.indelrealign.bam \
                -BQSR $id.recal_data.table \
                -o $id.bam
rm $id.indelrealign.bam
rm $id.indelrealign.bai
done

for i in *.bam
do
id=$(echo $i | awk -F'[.]' '{print $1}')

mv $id.bam $BAM_OUT/$SAMPLE.$id.bam
mv $id.bai $BAM_OUT/$SAMPLE.$id.bai

done

