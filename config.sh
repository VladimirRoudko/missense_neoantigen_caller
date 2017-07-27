STRELKA_INSTALL_DIR=/work/software/strelka

#REF_FASTA="/work/scratch/vladimir/annotations/b37/human_g1k_v37_decoy.fasta"
#REF_FASTA="/work/scratch/mutation/DATA/flu_ref.fa"
REF_FASTA="/work/scratch/TCGA-REF/GRCh38.d1.vd1.fa"


GOLD_INDEL="/work/scratch/TCGA-REF/Mills_and_1000G_gold_standard.indels.hg38.vcf"
KNOWN_INDEL="/work/scratch/TCGA-REF/Homo_sapiens_assembly38.known_indels.vcf"
DBSNP="/work/scratch/TCGA-REF/Homo_sapiens_assembly38.dbsnp138.vcf"
PHASE1_SNP="/work/scratch/TCGA-REF/1000G_phase1.snps.high_confidence.hg38.vcf"
HAPMAP="/work/scratch/TCGA-REF/hapmap_3.3.hg38.vcf"
OMNI="/work/scratch/TCGA-REF/1000G_omni2.5.hg38.vcf"
COSMIC="/work/scratch/TCGA-REF/CosmicALL_chr_M_sorted.vcf"

SNPEFF_CONFIG="/work/software/snpeff/snpEff.config"
SNPEFF_DB="GRCh38.82"

BAM_ROOT="/work/DATA2/hammer/lung/realigned_bam"
#HLA_REF_FASTA="/work/scratch/HLA-REF/IMGTHLA/ATHLATES_DB/fasta/hla_cds_nuc.fasta"
HLA_REF_FASTA="/work/software/athlates/db/ref/hla.clean.fasta"
#BED_FILE_PATH="/work/scratch/HLA-REF/IMGTHLA/ATHLATES_DB/bed"
BED_FILE_PATH="/work/software/athlates/db/bed"
#MSA_FILE_PATH="/work/scratch/HLA-REF/IMGTHLA/ATHLATES_DB/msa"
MSA_FILE_PATH="/work/software/athlates/db/msa"

HLA_SCRIPT="hla_typing.sge.sh"
HLA_DNA_REFERENCE="/work/software/optitype/data/hla_reference_dna.fasta"
HLA_CONFIG="/work/software/optitype/config.ini"



CONVERT="../convert.three-one-aa.R"
BIOMART="../biomart.enst.protein.R"

mkdir -p log

RG="RG"
LN="LN"
P="illumina"
PU="120"
SN="SN"
SC="mountsinai"
n="DNAseq"
t="bam"
m="PE"


#### only for DNA realignment:
#mkdir -p /work/DATA2/hammer/lung/realigned_bam
#BAM_OUT="/work/DATA2/hammer/lung/realigned_bam"








