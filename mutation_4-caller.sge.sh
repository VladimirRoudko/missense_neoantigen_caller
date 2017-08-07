#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 2
#$ -cwd
#$ -N mutation_caller
#$ -e log/
#$ -o log/
NT=2

. /etc/profile.d/sge.sh
. ~/.bashrc
. config.sh

sample=$(awk NR==$SGE_TASK_ID samplelist)
read IF_GERMLINE IF_SOMATIC SAMPLE <<< $sample
#cd "$SAMPLE"

mkdir -p "$SAMPLE"
cd "$SAMPLE"

exec 1>log_mutation.$SAMPLE.out
exec 2>log_mutation.$SAMPLE.err


if [[ $BAM_ROOT eq $BAM_OUT ]]; then
cp "$BAM_ROOT/$IF_GERMLINE" normal.bam
cp "$BAM_ROOT/${IF_GERMLINE%.bam}.bai" normal.bam.bai
#samtools index normal.bam
cp "$BAM_ROOT/$IF_SOMATIC" tumor.bam
cp "$BAM_ROOT/${IF_SOMATIC%.bam}.bai" tumor.bam.bai
#samtools index tumor.bam
else 
cp "$BAM_OUT/$SAMPLE.normal.bam" normal.bam
cp "$BAM_OUT/$SAMPLE.normal.bai" normal.bam.bai
#samtools index normal.bam
cp "$BAM_OUT/$SAMPLE.tumor.bam" tumor.bam
cp "$BAM_OUT/$SAMPLE.tumor.bai" tumor.bam.bai
#samtools index tumor.bam
fi

# Run HLA typing by optitype:
sh ../$HLA_SCRIPT $SAMPLE $IF_SOMATIC $IF_GERMLINE

# Run mutation callers
mkdir -p {mutect,varscan,somatic-sniper,vcf,strelka,samtools}



# === varscan ===
cd varscan
ln -s ../normal.bam
ln -s ../normal.bam.bai
ln -s ../tumor.bam
ln -s ../tumor.bam.bai
samtools mpileup -B -q 1 -f "$REF_FASTA" normal.bam tumor.bam \
	> nt.mpileup

varscan somatic nt.mpileup vs --mpileup 1 \
	--min-coverage 5 --min-var-freq 0.01 --somatic-p-value 0.05 \
	--output-vcf 1

for x in vs*.vcf
do
  varscan processSomatic "$x" &
done
wait

for x in vs*.Somatic.hc.vcf
do
  varscan somaticFilter "$x" --indel-file vs.indel.vcf \
    --output-file "${x}.filter" &
done
wait

: <<'END'
for x in *.filter
do
  awk '{print $1 "\t" $2 "\t" $2}' "${x}" | tail -n +2 > "${x}.mod"
  bam-readcount -q 1 -b 20 -f "$REF_FASTA" -l "${x}.mod" \
    tumor.bam > "${x}.readcount" &
done
wait

for x in *.filter
do
  perl /work/software/varscan/fpfilter.pl "$x" "${x}.readcount" \
    --output-basename "${x}.fpfilter" &
done
wait
END

rm nt.mpileup
cd ..

mv varscan/vs.indel.Somatic.hc.vcf.filter vcf/varscan.indel.vcf
mv varscan/vs.snp.Somatic.hc.vcf.filter vcf/varscan.snvs.vcf

for i in vcf/varscan.*.vcf
do
snpeff -c "$SNPEFF_CONFIG" -v "$SNPEFF_DB" $i > $i.snpeff.vcf & wait;
sed 's/PASS/PASS\tVS/' $i.snpeff.vcf > $i.snpeff.vcf.tmp
mv $i.snpeff.vcf.tmp $i.snpeff.vcf    
done

# === somatic-sniper ===
cd somatic-sniper

ln -s ../normal.bam
ln -s ../normal.bam.bai
ln -s ../tumor.bam
ln -s ../tumor.bam.bai

bam-somaticsniper -q 1 -Q 25 -G -L -F vcf -f "$REF_FASTA" \
  tumor.bam normal.bam ss.snv.vcf

cd ..

mv somatic-sniper/ss.snv.vcf vcf/ssniper.snvs.vcf

for i in vcf/ssniper.*.vcf
do
while read line
do
NORMAL=$(echo $line | awk '{print $(NF-1)}')
TUMOR=$(echo $line | awk '{print $NF}')
NREFCOV1=$(echo $NORMAL | awk -F"[,:]" '{print $4}')
NREFCOV2=$(echo $NORMAL | awk -F"[,:]" '{print $5}')
NALTCOV1=$(echo $NORMAL | awk -F"[,:]" '{print $6}')
NALTCOV2=$(echo $NORMAL | awk -F"[,:]" '{print $7}')
NCOV=$(echo $NORMAL | awk -F"[,:]" '{print $3}')
NALTCOV=$(expr "$NALTCOV1" + "$NALTCOV2")
TCOV=$(echo $TUMOR | awk -F":" '{print $3}')
TALTCOV1=$(echo $TUMOR | awk -F"[,:]" '{print $6}')
TALTCOV2=$(echo $TUMOR | awk -F"[,:]" '{print $6}')
TALTCOV=$(expr "$TALTCOV1" + "$TALTCOV2")

if [ $(echo "$TCOV > 10" | bc) -ne 0 ] && [ $(echo "$TALTCOV > 3" | bc) -ne 0 ] && [ $(echo "$NCOV > 7" | bc) -ne 0 ] && [ $(echo "$NALTCOV < 2" | bc) -ne 0 ]; then 
echo $line >> $i.qc
fi
done < $i
tr ' ' '\t' < $i.qc > $i.qc.tmp
mv $i.qc.tmp $i.qc
snpeff -c "$SNPEFF_CONFIG" -v "$SNPEFF_DB" $i.qc > $i.qc.snpeff.vcf & wait;
sed 's/\.\t\./\.\t\.\tSS/' $i.qc.snpeff.vcf > $i.qc.snpeff.vcf.tmp
mv $i.qc.snpeff.vcf.tmp $i.qc.snpeff.vcf
done

# === samtools ===
cd samtools

ln -s ../normal.bam
ln -s ../normal.bam.bai
ln -s ../tumor.bam
ln -s ../tumor.bam.bai

samtools mpileup -D -S -A -C50 -uf "$REF_FASTA" normal.bam tumor.bam | \
bcftools view -bvcgT pair - > var.raw.bcf

bcftools view var.raw.bcf | \
vcfutils.pl varFilter -d30 > var.flt.vcf

grep "0/0.*0/1" var.flt.vcf >> var.somatic.vcf
grep "0/0.*1/1" var.flt.vcf >> var.somatic.vcf 

rm var.raw.bcf
rm var.flt.vcf

cd ..

mv samtools/var.somatic.vcf vcf/samtools.snvs.vcf

for i in vcf/samtools.*.vcf
do
while read line
do
NORMAL=$(echo $line | awk '{print $(NF-1)}')
TUMOR=$(echo $line | awk '{print $NF}')
STAT=$(echo $line | awk '{print $8}')
DP4=$(echo $STAT | awk -F';' '{print $7}')
NCOV=$(echo $NORMAL | awk -F":" '{print $3}')
TCOV=$(echo $TUMOR | awk -F":" '{print $3}')

if [ $(echo "$TCOV > 10" | bc) -ne 0 ] && [ $(echo "$NCOV > 7" | bc) -ne 0 ] && ! [[ $line =~ "INDEL" ]]; then
echo $line >> $i.qc
fi
done < $i

tr ' ' '\t' < $i.qc > $i.qc.tmp
mv $i.qc.tmp $i.qc

snpeff -c "$SNPEFF_CONFIG" -v "$SNPEFF_DB" $i.qc > $i.qc.snpeff.vcf & wait;
sed 's/\.\tD/\.\tSA\tD/' $i.qc.snpeff.vcf > $i.qc.snpeff.vcf.tmp
mv $i.qc.snpeff.vcf.tmp $i.qc.snpeff.vcf

done


# === mutect ===
cd mutect

ln -s ../normal.bam
ln -s ../normal.bam.bai
ln -s ../tumor.bam
ln -s ../tumor.bam.bai

gatk4 Mutect2 \
  -R "$REF_FASTA" -normal normal -tumor tumor \
   -O mutect2.vcf --dbsnp "$DBSNP" --cosmic "$COSMIC" \
   -I tumor.bam -I normal.bam
gatk4 FilterMutectCalls --variant mutect2.vcf --output mutect2.filter.vcf
grep "PASS" mutect2.filter.vcf >> mutect2.called.vcf
cd ..

mv mutect/mutect2.called.vcf vcf/mutect.vcf

for i in vcf/mutect.vcf
do
snpeff -c "$SNPEFF_CONFIG" -v "$SNPEFF_DB" $i > $i.snpeff.vcf & wait;
sed 's/PASS/PASS\tMU/' $i.snpeff.vcf > $i.snpeff.vcf.tmp
mv $i.snpeff.vcf.tmp $i.snpeff.vcf    
done


### merge outputs:
for i in vcf/*.snpeff.vcf
do
cat $i >> vcf/$SAMPLE.merged.vcf
done
sort -V -k1,1 -k2,2 < vcf/$SAMPLE.merged.vcf > vcf/$SAMPLE.merged.sorted.vcf

rm *.bam
rm *.bai
