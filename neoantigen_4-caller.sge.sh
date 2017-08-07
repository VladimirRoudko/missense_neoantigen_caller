#!/bin/sh
#$ -S /bin/bash
#$ -pe smp 1
#$ -cwd
#$ -N neoantigen_caller
#$ -e log/
#$ -o log/

# qsub -t 1-N neoantigen_caller.sge.sh ...

# number of threads
NT=1

. /etc/profile.d/sge.sh

. config.sh

sample=$(awk NR==$SGE_TASK_ID samplelist)
read IF_GERMLINE IF_SOMATIC SAMPLE <<< $sample
#cd "$SAMPLE"

mkdir -p "$SAMPLE"
cd "$SAMPLE"

exec 1>log_neoantigen.$SAMPLE.out
exec 2>log_neoantigen.$SAMPLE.err

cp vcf/$SAMPLE.merged.sorted.vcf mutation.vcf
cp hla_optitype/$SAMPLE.HLA.tsv normal.hla.txt

### extract information from *.vcf file.  DATA extracted:
### 1 - LINE, names of callers, which called mutation
### 3 - LINE, ENSG_ID
### 4 - LINE, ENST_ID
### 5 - NUMBER, position of mutation
### 6 - WT aa residue
### 7 - MUT aa residue
### data is space delimited. stored in 

grep "missense" mutation.vcf |
  awk 'BEGIN {OFS="\t"; FS="\t"} {print $1"_"$2,$8,$1"_"$2"_"$10"_"$11"_"$12"_"$9}' |
awk 'END {for (k in c) print k, c[k], d[k] } { k = $1; c[k]=c[k]$2; d[k] = $3 }' |
sed "s/missense_variant/ /g" |
while read line; do out=$(echo $line | awk -F'[ ]' '{print $2,$3,$4}'); echo $out; done | while read line; do out1=$(echo $line | awk -F'[ ]' '{print $1,$2}'); out2=$(echo $line | awk -F'[ ]' '{print $3}'); out3=$(echo $out2 | awk -F'[ |]' '{print $4,$6,$10}'); out4=$(echo $out3 | awk -F'[ .]' '{print $1,$2,$NF}'); echo -e "$out1 $out4"; done |
while read line 
do 
NUMBER=$(echo $line | awk '{print $NF}' | sed 's/[^0-9]*//g');
MUTANT=$(echo $line | awk '{print $NF}' | sed 's/[0-9]/ /g' | awk '{print $NF}') ;
res=$(Rscript $CONVERT $MUTANT)
M=$(echo $res | awk -F'"' '{print $2}')
WT=$(echo $line | awk '{print $NF}' | sed 's/[0-9]/ /g' | awk '{print $1}') ;
res=$(Rscript $CONVERT $WT)
W=$(echo $res | awk -F'"' '{print $2}')
LINE=$(echo $line | awk '{print $1,$2,$3,$4}');
echo -e "$LINE $NUMBER $W $M" >> mutation_coordinates.txt; 
done

##### extract protein sequences containing missense mutation: 
Rscript $BIOMART mutation_coordinates.txt protein.fasta
tr -d '*' < protein.fasta > tmp.txt; 
mv tmp.txt protein.fasta;
sed -i '/^$/d' protein.fasta     

##### extract peptides.  Extract WT and MUT peptides, 
##### for 9-mers: extract -8/+8 around mutant aa. -> 17 aa peptide
##### add also the line with quality stats of the mutation at the end
while read line;
do
ID_SEQ=$(echo $line | awk -F' ' '{print $4}')
ID_GENE=$(echo $line | awk -F' ' '{print $3}')
ANN=$(echo $line | awk -F' ' '{print $2}')
CALLER=$(echo $line | awk -F' ' '{print $1}')
if grep -q $ID_SEQ "protein.fasta"; then
SEQ=$(grep -A1 $ID_SEQ "protein.fasta" | tail -n 1)
MUT_POS=$(echo $line | awk -F' ' '{print $5}')
WT_AA=$(echo $line | awk -F' ' '{print $6}')
MUT_AA=$(echo $line | awk -F' ' '{print $7}')
BORDER="8"
START=$(expr "$MUT_POS" - "$BORDER")
END=$(expr "$MUT_POS" + "$BORDER")
WT_PEPTIDE=$(samtools faidx "protein.fasta" $ID_SEQ:$START-$END)
WT_PEPTIDE=$(echo $WT_PEPTIDE | awk '{print $2}')
MUT_PEPTIDE=$(echo $WT_PEPTIDE | sed -r "s/^(.{$BORDER})$WT_AA(.*)$/\1$MUT_AA\2/")
echo -e "$CALLER\t$ID_GENE\t$ID_SEQ\t$WT_AA\t$MUT_AA\t$WT_PEPTIDE\t$MUT_PEPTIDE\t$ANN" >> peptide.9-17.stats.txt
fi
done < mutation_coordinates.txt

rm "protein.fasta"
rm "protein.fasta.fai"

## make fasta file from peptide sequences
ID2="1"
ONE="1"
while read line;
do
CALLER=$(echo $line | awk '{print $1}')
ID=$(echo $line | awk '{print $3}' | sed 's/[^0-9]*//g')
WT=$(echo $line | awk '{print $6}')
MUT=$(echo $line | awk '{print $7}')
ANN=$(echo $line | awk '{print $8}')
echo ">W"$CALLER"_$ID2" >> peptide.9-17.fasta
echo "$WT" >> peptide.9-17.fasta
echo ">M"$CALLER"_$ID2" >> peptide.9-17.fasta
echo "$MUT" >> peptide.9-17.fasta
ID2=$(($ID2+$ONE))
done < peptide.9-17.stats.txt

##### extract HLA information from HLA predictions
A_1=$(tail -1 normal.hla.txt | awk '{print $2}')
A_2=$(tail -1 normal.hla.txt | awk '{print $3}')

B_1=$(tail -1 normal.hla.txt | awk '{print $4}')
B_2=$(tail -1 normal.hla.txt | awk '{print $5}')

C_1=$(tail -1 normal.hla.txt | awk '{print $6}')
C_2=$(tail -1 normal.hla.txt | awk '{print $7}')

HLA=$(echo "$A_1 $A_2 $B_1 $B_2 $C_1 $C_2" | sed 's/\*//g' | sed 's/\://g' | sed "s/'//g")
HLApan=$(echo "$A_1 $A_2 $B_1 $B_2 $C_1 $C_2" | sed 's/\*//g' | sed "s/'//g")
list_MHC=$(netMHC -listMHC)
list_MHCpan=$(netMHCpan -listMHC)

############################### RUN netMHC with predicted peptides
for type in $HLA; 
do 
if [[ "$list_MHC" =~ "$type" ]]; 
then 
netMHC -a "HLA-$type" -f "peptide.9-17.fasta" -l 9 -s 1 -xls 1 -xlsfile "$type.peptide.9-17.xls"
fi
done

for output_file in *.xls
do
hla_type=$(head -1 $output_file | awk '{print $1}')
while read line
do
N_binders=$(echo $line | awk '{print $NF}')
if [[ $N_binders -eq "1" ]]; then
echo $line >> "$hla_type.peptide.binders.txt"
fi
done < $output_file
while read line
do
##### add ANN from peptide.stat file to the final output
ID=$(echo $line | awk '{print $3}' | awk -F'_' '{print $2}')
line_ANN=$(awk NR==$ID peptide.9-17.stats.txt)
ANN=$(echo $line_ANN | awk '{print $4,$5,$NF}')
echo -e "$hla_type\t$line\t$ANN" >> "$SAMPLE.MHCI.tumor.binders.9mer.netHMC.txt"
done < "$hla_type.peptide.binders.txt"
rm "$hla_type.peptide.binders.txt"
done
mkdir "netMHC"
mv *.xls "netMHC"

############################### RUN netMHCpan with predicted peptides
for type in $HLApan; 
do 
if [[ "$list_MHCpan" =~ "$type" ]]; 
then 
netMHCpan -a "HLA-$type" -f "peptide.9-17.fasta" -l 9 -s 1 -xls 1 -xlsfile "$type.peptide.9-17.xls"
fi
done

for output_file in *.xls
do
hla_type=$(head -1 $output_file | awk '{print $1}')
while read line
do
N_binders=$(echo $line | awk '{print $NF}')
if [[ $N_binders -eq "1" ]]; then
echo $line >> "$hla_type.peptide.binders.txt"
fi
done < $output_file
while read line
do
##### add ANN from peptide.stat file to the final output
ID=$(echo $line | awk '{print $3}' | awk -F'_' '{print $2}')
line_ANN=$(awk NR==$ID peptide.9-17.stats.txt)
ANN=$(echo $line_ANN | awk '{print $4,$5,$NF}')
echo -e "$hla_type\t$line\t$ANN" >> "$SAMPLE.MHCI.tumor.binders.9mer.netHMCpan.txt"
done < "$hla_type.peptide.binders.txt"
rm "$hla_type.peptide.binders.txt"
done
mkdir "netMHCpan"
mv *.xls "netMHCpan"

mkdir "neoantigens"
mv $SAMPLE.MHCI.* "neoantigens"

cd neoantigens

#### merge all outputs in one merged file:
for i in *.txt
do
while read line;
do
flag=$(echo $line | awk '{print $5}')
re='^[0-9]+'
if [[ $flag =~ $re ]]; then
text=$(echo $line | awk '{print $1,$2,$3,$4,$5,$(NF-2),$(NF-1),$NF}')
echo $text >> "$SAMPLE.merged.txt"
else
text=$(echo $line | awk '{print $1,$2,$3,$4,$7,$(NF-2),$(NF-1),$NF}')
echo $text >> "$SAMPLE.merged.txt"  
fi
done < $i
done

### output format for marta:
#echo -e "HLA\tmutAA_position\tpeptide\tnormal_AA\tmutant_AA\tID\tKD\tchromosome\tposition\tn_ref\tn_alt\tt_ref\tt_alt\tannotation" >> $SAMPLE.output.txt
while read line
do
hla=$(echo $line | awk '{print $1}')
mutation=$(echo $line | awk '{print $2}')
middle="9"
mutation=$(($middle - $mutation))
peptide=$(echo $line | awk '{print $3}')
ID=$(echo $line | awk '{print $4}')
KD=$(echo $line | awk '{print $5}')
if ! [[ $peptide == *"X"* ]]; then
WT_AA=$(echo $line | awk '{print $(NF-2)}')
MUT_AA=$(echo $line | awk '{print $(NF-1)}')
annotation=$(echo $line | awk '{print $NF}')
chr=$(echo $annotation | awk -F'_' '{print $1}')
position=$(echo $annotation | awk -F'_' '{print $2}')

if [[ $annotation == *"FOXOG"* ]]; then
n_ref=$(echo $annotation | awk -F'/' '{print $2}' | awk -F'[:,]' '{print $2}')
n_alt=$(echo $annotation | awk -F'/' '{print $2}' | awk -F'[:,]' '{print $3}')
t_ref=$(echo $annotation | awk -F'/' '{print $3}' | awk -F'[:,]' '{print $2}')
t_alt=$(echo $annotation | awk -F'/' '{print $3}' | awk -F'[:,]' '{print $3}')
echo "$hla $mutation $peptide $WT_AA $MUT_AA $ID $KD $chr $position $n_ref $n_alt $t_ref $t_alt $annotation" >> $SAMPLE.output.txt
fi

if [[ $annotation == *"BCOUNT"* ]]; then
n_ref1=$(echo $annotation | awk -F'/' '{print $3}' | awk -F'[:,]' '{print $3}')
n_ref2=$(echo $annotation | awk -F'/' '{print $3}' | awk -F'[:,]' '{print $4}')
n_alt1=$(echo $annotation | awk -F'/' '{print $3}' | awk -F'[:,]' '{print $5}')
n_alt2=$(echo $annotation | awk -F'/' '{print $3}' | awk -F'[:,]' '{print $6}')
t_ref1=$(echo $annotation | awk -F'/' '{print $5}' | awk -F'[:,]' '{print $3}')
t_ref2=$(echo $annotation | awk -F'/' '{print $5}' | awk -F'[:,]' '{print $4}')
t_alt1=$(echo $annotation | awk -F'/' '{print $5}' | awk -F'[:,]' '{print $5}')
t_alt2=$(echo $annotation | awk -F'/' '{print $5}' | awk -F'[:,]' '{print $6}')
n_ref=$(($n_ref1 + $n_ref2))
n_alt=$(($n_alt1 + $n_alt2))
t_ref=$(($t_ref1 + $t_ref2))
t_alt=$(($t_alt1 + $t_alt2))
echo "$hla $mutation $peptide $WT_AA $MUT_AA $ID $KD $chr $position $n_ref $n_alt $t_ref $t_alt $annotation" >> $SAMPLE.output.txt
fi

if [[ $annotation == *"FREQ"* ]]; then
n_ref=$(echo $annotation | awk -F'/' '{print $2}' | awk -F'[:]' '{print $4}')
n_alt=$(echo $annotation | awk -F'/' '{print $2}' | awk -F'[:]' '{print $5}')
t_ref=$(echo $annotation | awk -F'/' '{print $3}' | awk -F'[:]' '{print $4}')
t_alt=$(echo $annotation | awk -F'/' '{print $3}' | awk -F'[:]' '{print $5}')
echo "$hla $mutation $peptide $WT_AA $MUT_AA $ID $KD $chr $position $n_ref $n_alt $t_ref $t_alt $annotation" >> $SAMPLE.output.txt
fi

if [[ $annotation == *"RPB="* ]]; then
n_ref=$(echo $annotation | awk -F'[_]' '{print $4}' | awk -F':' '{print $3}')
t_dp=$(echo $annotation | awk -F'[_]' '{print $5}' | awk -F':' '{print $3}')
t_alt1=$(echo $annotation | awk -F'[_]' '{print $6}' | awk -F'[;]' '{print $7}' | awk -F',' '{print $3}')
t_alt2=$(echo $annotation | awk -F'[_]' '{print $6}' | awk -F'[;]' '{print $7}' | awk -F',' '{print $4}')
t_alt=$(expr "$t_alt1" + "$t_alt2")
if [ $(echo "$t_dp > $t_alt" | bc) -ne 0 ]; then
t_ref=$(expr "$t_dp" - "$t_alt")
else
t_ref="NA"
fi
n_alt="0"
echo "$hla $mutation $peptide $WT_AA $MUT_AA $ID $KD $chr $position $n_ref $n_alt $t_ref $t_alt $annotation" >> $SAMPLE.output.txt
fi

fi
done < $SAMPLE.merged.txt

#### Filter output
while read line
do
hla=$(echo $line | awk '{print $1}')
mutAA_position=$(echo $line | awk '{print $2}')
ONE=1
BORDER=$(expr $mutAA_position - $ONE)
MUT_peptide=$(echo $line | awk '{print $3}')
ID=$(echo $line | awk '{print $6}')
MUT_KD=$(echo $line | awk '{print $7}')
WT_AA=$(echo $line | awk '{print $4}')
MUT_AA=$(echo $line | awk '{print $5}')
chr=$(echo $line | awk '{print $8}')
position=$(echo $line | awk '{print $9}')
n_ref=$(echo $line | awk '{print $(10)}')
n_alt=$(echo $line | awk '{print $(11)}')
t_ref=$(echo $line | awk '{print $(12)}')
t_alt=$(echo $line | awk '{print $(13)}')
annotation=$(echo $line | awk '{print $(14)}')
THRESHOLD="500"
if [[ $ID == "M"* ]]; then
if [[ $(bc <<< "$MUT_KD <= $THRESHOLD") -eq 1 ]]; then
WT_PEPTIDE=$(echo $MUT_peptide | sed -r "s/^(.{$BORDER})$MUT_AA(.*)$/\1$WT_AA\2/")
search=$(echo "$mutAA_position $WT_PEPTIDE $WT_AA $MUT_AA")
NORMAL=$(grep "$search" $SAMPLE.output.txt)
if [[ $NORMAL == "" ]]; then 
WT_peptide="NA"
WT_KD="NA"
else
WT_peptide=$(echo $NORMAL | awk '{print $3}')
WT_KD=$(echo $NORMAL | awk '{print $7}')
fi
echo "$SAMPLE $hla $mutAA_position $WT_peptide $MUT_peptide $WT_AA $MUT_AA $WT_KD $MUT_KD $ID $chr $position $n_ref $n_alt $t_ref $t_alt $annotation" >> $SAMPLE.neoantigens.txt
fi
fi
done < $SAMPLE.output.txt
mv $SAMPLE.output.txt $SAMPLE.temp2.txt
mv $SAMPLE.merged.txt $SAMPLE.temp1.txt

cd ..

mkdir "tmp"
mv *.txt tmp
mv *.fasta tmp
mv log_* tmp
mv *.vcf tmp
mv config.ini tmp








