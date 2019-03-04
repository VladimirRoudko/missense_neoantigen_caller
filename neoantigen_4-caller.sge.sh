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

sample=$(awk NR==$SGE_TASK_ID sample_list)
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
  awk 'BEGIN {OFS="\t"; FS="\t"} {print $1"_"$2"_"$4"_"$5,$8,$1"_"$2"_"$4"_"$5"_"$10"_"$11"_"$12"_"$9}' |
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
echo -e "$CALLER\t$ID_GENE\t$ID_SEQ\t$WT_AA\t$MUT_AA\t$WT_PEPTIDE\t$MUT_PEPTIDE\t$ANN" >> peptide.stats.txt
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
if ! [[ -z "$MUT" ]]; then
echo ">W"$CALLER"_$ID2" >> peptide.fasta
echo "$WT" >> peptide.fasta
echo ">M"$CALLER"_$ID2" >> peptide.fasta
echo "$MUT" >> peptide.fasta
ID2=$(($ID2+$ONE))
else
ID2=$(($ID2+$ONE))
fi
done < peptide.stats.txt

##### extract HLA information from HLA predictions
A_1=$(tail -1 normal.hla.txt | awk '{print $2}')
A_2=$(tail -1 normal.hla.txt | awk '{print $3}')

B_1=$(tail -1 normal.hla.txt | awk '{print $4}')
B_2=$(tail -1 normal.hla.txt | awk '{print $5}')

C_1=$(tail -1 normal.hla.txt | awk '{print $6}')
C_2=$(tail -1 normal.hla.txt | awk '{print $7}')

HLAv4=$(echo "$A_1 $A_2 $B_1 $B_2 $C_1 $C_2" | sed 's/\*//g' | sed 's/\://g' | sed "s/'//g")
HLAv3=$(echo "$A_1 $A_2 $B_1 $B_2 $C_1 $C_2" | sed 's/\*//g' | sed "s/'//g")
HLAPANv3=$(echo "$A_1 $A_2 $B_1 $B_2 $C_1 $C_2" | sed 's/\*//g' | sed "s/'//g")
HLAPANv4=$(echo "$A_1 $A_2 $B_1 $B_2 $C_1 $C_2" | sed 's/\*//g' | sed "s/'//g")

list_MHCv3=$(netMHCv3 -A)
list_MHCv4=$(netMHCv4 -listMHC)
list_MHCPANv3=$(netMHCPANv3 -listMHC)
list_MHCPANv4=$(netMHCPANv4 -listMHC)

############################ RUN netMHC-3.4 with predicted peptides
for type in $HLAv3
do 
if [[ "$list_MHCv3" =~ "$type" ]];
then
netMHCv3 -a "HLA-$type" -l 9 --xls="$type.peptide.xls" "peptide.fasta" > $type.MHCv3.4.peptide.txt
fi
done
rm *.xls

for file in *.peptide.txt
do
HLA_TYPE=$(echo $file | awk -F'.' '{print $1}')
HLA_TYPE=$(echo $HLA_TYPE | sed 's/\://g')
while read line
do
MPOS=$(echo $line | awk '{print $1}')
MKD=$(echo $line | awk '{print $4}')
MPEPTIDE=$(echo $line | awk '{print $2}')
MID=$(echo $line | awk '{print $(NF-1)}')
THRESHOLD=500
if [[ $(bc <<< "$MKD <= $THRESHOLD") -eq 1  ]]; then
if [[ $MID =~ ^M ]]; then
echo -e "$MID\t$MPOS\t$MKD\t$MPEPTIDE" >> "$HLA_TYPE.netMHCv3.txt"
WID=$(echo $MID | sed -r "s/^M/W/")
grep "\<$WID\>" $file |
while read text
do
WPOS=$(echo $text | awk '{print $1}')
WPEPTIDE=$(echo $text | awk '{print $2}')
if [[ $WPOS == $MPOS ]]; then
WKD=$(echo $text | awk '{print $4}')
echo -e "$WID\t$WPOS\t$WKD\t$WPEPTIDE" >> "$HLA_TYPE.netMHCv3.txt"
fi
done
fi
fi
done < $file

##### add ANN from peptide.stat file to the final output
while read line
do
ID=$(echo $line | awk '{print $1}' | awk -F'_' '{print $2}')
LINE_ANN=$(awk NR==$ID peptide.stats.txt)
ANN=$(echo $LINE_ANN | awk '{print $4,$5,$NF}')
echo -e "$HLA_TYPE\t"netMHCv3.4"\t$line\t$ANN" >> "$SAMPLE.netMHCv3.txt"
done < "$HLA_TYPE.netMHCv3.txt"
rm "$HLA_TYPE.netMHCv3.txt"
done
mkdir -p "netMHCv3"
mv *.peptide.txt "netMHCv3"

############################### RUN netMHC-4.0 with predicted peptides
for type in $HLAv4; 
do 
if [[ "$list_MHCv4" =~ "$type" ]]; 
then 
netMHCv4 -a "HLA-$type" -f "peptide.fasta" -l 9 -s 1 -xls 1 -xlsfile "$type.MHCv4.peptide.xls"
fi
done

for file in *.xls
do
HLA_TYPE=$(head -1 $file | awk '{print $1}')
HLA_TYPE=$(echo $HLA_TYPE | sed 's/\://g')
while read line
do
MPOS=$(echo $line | awk '{print $1}')
MKD=$(echo $line | awk '{print $4}')
MID=$(echo $line | awk '{print $3}')
MPEPTIDE=$(echo $line | awk '{print $2}')
THRESHOLD=500
if [[ $(bc <<< "$MKD <= $THRESHOLD") -eq 1  ]]; then
if [[ $MID =~ ^M ]]; then
echo -e "$MID\t$MPOS\t$MKD\t$MPEPTIDE" >> "$HLA_TYPE.netMHCv4.txt"
WID=$(echo $MID | sed -r "s/^M/W/")
grep "\<$WID\>" $file |
while read text
do
WPOS=$(echo $text | awk '{print $1}')
if [[ $WPOS == $MPOS ]]; then
WKD=$(echo $text | awk '{print $4}')
WPEPTIDE=$(echo $text | awk '{print $2}')
echo -e "$WID\t$WPOS\t$WKD\t$WPEPTIDE" >> "$HLA_TYPE.netMHCv4.txt"
fi
done
fi
fi
done < $file

##### add ANN from peptide.stat file to the final output
while read line
do
ID=$(echo $line | awk '{print $1}' | awk -F'_' '{print $2}')
line_ANN=$(awk NR==$ID peptide.stats.txt)
ANN=$(echo $line_ANN | awk '{print $4,$5,$NF}')
echo -e "$HLA_TYPE\t"netMHCv4.0"\t$line\t$ANN" >> "$SAMPLE.netMHCv4.txt"
done < "$HLA_TYPE.netMHCv4.txt"
rm "$HLA_TYPE.netMHCv4.txt"
done
mkdir -p "netMHCv4"
mv *.xls "netMHCv4"

############################### RUN netMHCpan-3.0 with predicted peptides
for type in $HLAPANv3; 
do 
if [[ "$list_MHCPANv3" =~ "$type" ]]; 
then 
netMHCPANv3 -a "HLA-$type" -f "peptide.fasta" -l 9 -s 1 -xls 1 -xlsfile "$type.MHCPANv3.peptide.xls"
fi
done

for file in *.xls
do
HLA_TYPE=$(head -1 $file | awk '{print $1}')
HLA_TYPE=$(echo $HLA_TYPE | sed 's/\://g')
while read line
do
MPOS=$(echo $line | awk '{print $1}')
MKD=$(echo $line | awk '{print $6}')
MID=$(echo $line | awk '{print $3}')
MPEPTIDE=$(echo $line | awk '{print $2}')
THRESHOLD=500
if [[ $(bc <<< "$MKD <= $THRESHOLD") -eq 1  ]]; then
if [[ $MID =~ ^M ]]; then
echo -e "$MID\t$MPOS\t$MKD\t$MPEPTIDE" >> "$HLA_TYPE.netMHCPANv3.txt"
WID=$(echo $MID | sed -r "s/^M/W/")
grep "\<$WID\>" $file |
while read text
do
WPOS=$(echo $text | awk '{print $1}')
if [[ $WPOS == $MPOS ]]; then
WKD=$(echo $text | awk '{print $6}')
WPEPTIDE=$(echo $text | awk '{print $2}')
echo -e "$WID\t$WPOS\t$WKD\t$WPEPTIDE" >> "$HLA_TYPE.netMHCPANv3.txt"
fi
done
fi
fi
done < $file

##### add ANN from peptide.stat file to the final output
while read line
do
ID=$(echo $line | awk '{print $1}' | awk -F'_' '{print $2}')
line_ANN=$(awk NR==$ID peptide.stats.txt)
ANN=$(echo $line_ANN | awk '{print $4,$5,$NF}')
echo -e "$HLA_TYPE\t"netMHCPANv3.0"\t$line\t$ANN" >> "$SAMPLE.netMHCPANv3.txt"
done < "$HLA_TYPE.netMHCPANv3.txt"
rm "$HLA_TYPE.netMHCPANv3.txt"
done
mkdir -p "netMHCPANv3"
mv *.xls "netMHCPANv3"


############################### RUN netMHCpan-4.0 with predicted peptides
for type in $HLAPANv4;
do
if [[ "$list_MHCPANv4" =~ "$type" ]];
then
netMHCPANv4 -a "HLA-$type" -f "peptide.fasta" -l 9 -s 1 -BA -xls 1 -xlsfile "$type.MHCPANv4.peptide.xls"
fi
done
for type in $HLAv4;
do
if [[ "$list_MHCPANv4" =~ "$type" ]];
then
netMHCPANv4 -a "HLA-$type" -f "peptide.fasta" -l 9 -s 1 -BA -xls 1 -xlsfile "$type.MHCPANv4.peptide.xls"
fi
done

for file in *.xls
do
HLA_TYPE=$(head -1 $file | awk '{print $1}')
HLA_TYPE=$(echo $HLA_TYPE | sed 's/\://g')
while read line
do
MPOS=$(echo $line | awk '{print $1}')
MSCORE=$(echo $line | awk '{print $7}')
MID=$(echo $line | awk '{print $3}')
MPEPTIDE=$(echo $line | awk '{print $2}')
THRESHOLD=500
if [[ $(bc <<< "$MSCORE <= $THRESHOLD") -eq 1  ]]; then
if [[ $MID =~ ^M ]]; then
WID=$(echo $MID | sed -r "s/^M/W/")
echo -e "$MID\t$MPOS\t$MSCORE\t$MPEPTIDE" >> "$HLA_TYPE.netMHCPANv4.txt"
grep "\<$WID\>" $file |
while read text
do
WPOS=$(echo $text | awk '{print $1}')
if [[ $WPOS == $MPOS ]]; then
WSCORE=$(echo $text | awk '{print $7}')
WPEPTIDE=$(echo $text | awk '{print $2}')
echo -e "$WID\t$WPOS\t$WSCORE\t$WPEPTIDE" >> "$HLA_TYPE.netMHCPANv4.txt"
fi
done
fi
fi
done < $file

##### add ANN from peptide.stat file to the final output
while read line
do
ID=$(echo $line | awk '{print $1}' | awk -F'_' '{print $2}')
line_ANN=$(awk NR==$ID peptide.stats.txt)
ANN=$(echo $line_ANN | awk '{print $4,$5,$NF}')
echo -e "$HLA_TYPE\t"netMHCPANv4.0"\t$line\t$ANN" >> "$SAMPLE.netMHCPANv4.txt"
done < "$HLA_TYPE.netMHCPANv4.txt"
rm "$HLA_TYPE.netMHCPANv4.txt"
done
mkdir -p "netMHCPANv4"
mv *.xls "netMHCPANv4"

mkdir "neoantigens"
mv $SAMPLE.netMHC* "neoantigens"
cd "neoantigens"

#### merge all outputs in one merged file:
for i in *.txt
do
flag=$(echo $i | awk -F'.' '{print $1}')
while read line
do
echo -e "$flag\t$line" >> "$SAMPLE.merged.txt"
done < $i
done

### output format for marta:
#echo -e "SAMPLE\tTOOL\tHLA\tmutAA_position\tpeptide\tnormal_AA\tmutant_AA\tID\tKD\tchromosome\tposition\tn_ref\tn_alt\tt_ref\tt_alt\tannotation" >> $SAMPLE.output.txt
while read line
do
sample=$(echo $line | awk '{print $1}')
hla=$(echo $line | awk '{print $2}')
program=$(echo $line | awk '{print $3}')
mutation=$(echo $line | awk '{print $5}')
middle="9"
mutation=$(($middle - $mutation))
peptide=$(echo $line | awk '{print $7}')
ID=$(echo $line | awk '{print $4}')
KD=$(echo $line | awk '{print $6}')
if ! [[ $peptide == *"X"* ]]; then
WT_AA=$(echo $line | awk '{print $(NF-2)}')
MUT_AA=$(echo $line | awk '{print $(NF-1)}')
annotation=$(echo $line | awk '{print $NF}')
chr=$(echo $annotation | awk -F'_' '{print $1}')
position=$(echo $annotation | awk -F'_' '{print $2}')
ref_allele=$(echo $annotation | awk -F'_' '{print $3}')
alt_allele=$(echo $annotation | awk -F'_' '{print $4}')

if [[ $annotation == *"FOXOG"* ]]; then
n_ref=$(echo $annotation | awk -F'/' '{print $2}' | awk -F'[:,]' '{print $2}')
n_alt=$(echo $annotation | awk -F'/' '{print $2}' | awk -F'[:,]' '{print $3}')
t_ref=$(echo $annotation | awk -F'/' '{print $3}' | awk -F'[:,]' '{print $2}')
t_alt=$(echo $annotation | awk -F'/' '{print $3}' | awk -F'[:,]' '{print $3}')
echo "$sample $program $hla $mutation $peptide $WT_AA $MUT_AA $ID $KD $chr $position $ref_allele $alt_allele $n_ref $n_alt $t_ref $t_alt $annotation" >> $SAMPLE.output.txt
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
echo "$sample $program $hla $mutation $peptide $WT_AA $MUT_AA $ID $KD $chr $position $ref_allele $alt_allele $n_ref $n_alt $t_ref $t_alt $annotation" >> $SAMPLE.output.txt
fi

if [[ $annotation == *"FREQ"* ]]; then
n_ref=$(echo $annotation | awk -F'/' '{print $2}' | awk -F'[:]' '{print $4}')
n_alt=$(echo $annotation | awk -F'/' '{print $2}' | awk -F'[:]' '{print $5}')
t_ref=$(echo $annotation | awk -F'/' '{print $3}' | awk -F'[:]' '{print $4}')
t_alt=$(echo $annotation | awk -F'/' '{print $3}' | awk -F'[:]' '{print $5}')
echo "$sample $program $hla $mutation $peptide $WT_AA $MUT_AA $ID $KD $chr $position $ref_allele $alt_allele $n_ref $n_alt $t_ref $t_alt $annotation" >> $SAMPLE.output.txt
fi

if [[ $annotation == *"RPB="* ]]; then
n_ref=$(echo $annotation | awk -F'[_]' '{print $6}' | awk -F':' '{print $3}')
t_dp=$(echo $annotation | awk -F'[_]' '{print $7}' | awk -F':' '{print $3}')
t_alt1=$(echo $annotation | awk -F'[_]' '{print $8}' | awk -F'[;]' '{print $7}' | awk -F',' '{print $3}')
t_alt2=$(echo $annotation | awk -F'[_]' '{print $8}' | awk -F'[;]' '{print $7}' | awk -F',' '{print $4}')
t_alt=$(expr "$t_alt1" + "$t_alt2")
if [ $(echo "$t_dp > $t_alt" | bc) -ne 0 ]; then
t_ref=$(expr "$t_dp" - "$t_alt")
else
t_ref="NA"
fi
n_alt="0"
echo "$sample $program $hla $mutation $peptide $WT_AA $MUT_AA $ID $KD $chr $position $ref_allele $alt_allele $n_ref $n_alt $t_ref $t_alt $annotation" >> $SAMPLE.output.txt
fi

fi
done < $SAMPLE.merged.txt

#### Filter output
echo -e "SAMPLE\tTOOL\tHLA_allele\tmutAA_position\tWT_peptide\tMUT_peptide\tWT_AA\tMUT_AA\tWT_KD\tMUT_KD\tID\tchr\tposition\tref_allele\talt_allele\tNr\tNa\tTr\tTa\tannotation" >> $SAMPLE.neoantigens.txt
while read line
do
SAMPLE=$(echo $line | awk '{print $1}')
Mtool=$(echo $line | awk '{print $2}')
Mhla=$(echo $line | awk '{print $3}')
mutAA_position=$(echo $line | awk '{print $4}')
Mpeptide=$(echo $line | awk '{print $5}')
MID=$(echo $line | awk '{print $8}')
MKD=$(echo $line | awk '{print $9}')
WT_AA=$(echo $line | awk '{print $6}')
MUT_AA=$(echo $line | awk '{print $7}')
chr=$(echo $line | awk '{print $10}')
position=$(echo $line | awk '{print $(11)}')
ref_allele=$(echo $line | awk '{print $(12)}')
alt_allele=$(echo $line | awk '{print $(13)}')
n_ref=$(echo $line | awk '{print $(14)}')
n_alt=$(echo $line | awk '{print $(15)}')
t_ref=$(echo $line | awk '{print $(16)}')
t_alt=$(echo $line | awk '{print $(17)}')
annotation=$(echo $line | awk '{print $(18)}')
if [[ $MID =~ ^M ]]; then
WID=$(echo $MID | sed -r "s/^M/W/")
while read line2
do
ID=$(echo $line2 | awk '{print $8}')
if [[ $ID == $WID ]]; then
Wtool=$(echo $line2 | awk '{print $2}')
if [[ $Wtool == $Mtool ]]; then
AA_position=$(echo $line2 | awk '{print $4}')
if [[ $AA_position == $mutAA_position ]]; then
WKD=$(echo $line2 | awk '{print $9}')
Wpeptide=$(echo $line2 | awk '{print $5}')
echo -e "$SAMPLE\t$Mtool\t$Mhla\t$mutAA_position\t$Wpeptide\t$Mpeptide\t$WT_AA\t$MUT_AA\t$WKD\t$MKD\t$MID\t$chr\t$position\t$ref_allele\t$alt_allele\t$n_ref\t$n_alt\t$t_ref\t$t_alt\t$annotation" >> $SAMPLE.neoantigens.txt
break;
fi
fi
fi
done < $SAMPLE.output.txt
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








