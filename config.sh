STRELKA_INSTALL_DIR=/work/software/strelka

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

HLA_SCRIPT="hla_4-caller.sge.sh"
HLA_DNA_REFERENCE="/work/software/optitype/data/hla_reference_dna.fasta"
HLA_CONFIG="/work/software/optitype/config.ini"

CONVERT="../convert.three-one-aa.R"
BIOMART="../biomart.enst.protein.R"

RG="RG"
LN="LN"
P="illumina"
PU="120"
SN="SN"
SC="mountsinai"
n="DNAseq"
t="bam"
m="PE"


mkdir -p log
#### only for DNA realignment:
mkdir -p /work/DATA2/neoantigen_pipeline/bladder/realigned_bam
BAM_OUT="/work/DATA2/neoantigen_pipeline/bladder/realigned_bam"
BAM_ROOT="/work/TCGA/msk-bladder/bams"







