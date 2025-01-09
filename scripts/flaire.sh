#!/bin/bash
#$ -pe smp 16
#$ -N vir_FLAIR
#$ -cwd

# Abort on error
set -e
# ILLUMINA DATA =  /home/pthorpe001/ggs_simpson/demeth_mutants_temperature_variable/data/processed_illumina/reads

# Script to run flair on a whole load of bam files. Reads already mapped to the genome
# Long reads!

# ##################################################
# INSTALL instructions ....
# conda create -n flair -c conda-forge -c bioconda flair cd-hit
# https://flair.readthedocs.io/en/latest/requirements.html


#conda activate flair


##################################################
# .... FILL THESE IN .....
# Define the path to the directory and other resources
##################################################

GENOME_FA="../genomic_data/TAIR10_chr_all.fas"
ANNOTATION_GTF="../genomic_data/atRTD3_TS_21Feb22_transfix.gtf"
THREADS=16
EXP_NAME="vir_temp"
bams=*sorted.bam

# Directory where the BAM files are stored, assumed to be the current directory
BAM_DIR="."

# Directory to store the outputs, here also in the current directory
OUTPUT_DIR="./${EXP_NAME}_atRTD3_flair_outputs"
mkdir -p $OUTPUT_DIR

# ... NO MORE TO FILE IN FROM HERE ... 
##################################################


# Loop through each BAM file in the directory
for BAM_FILE in ${bams}
do
    # Extract the base name of the file without the path and extension
    BASE_NAME=$(basename $BAM_FILE .bam)

    ##################################################
    # extract reads
    cmd="samtools fasta -F 4 ${BAM_FILE} > ${BASE_NAME}.fq"
    output_file="${BASE_NAME}.fq"

    # Check if the output file exists and is not empty
    if [ ! -s "$output_file" ]; then
        echo "Running command: $cmd"
        eval $cmd
    else
        echo "Output file $output_file already exists and is not empty."
    fi

    ##################################################
    # Correct step

    bam2Bed12 -i ${BAM_FILE} > ${BAM_FILE}.bed12
    output_file=$OUTPUT_DIR/${BASE_NAME}*_all_corrected.bed
    cmd="flair correct -g $GENOME_FA  
        -q ${BAM_FILE}.bed12 
        --threads ${THREADS} 
        -o $OUTPUT_DIR/${BASE_NAME}_corrected
        --gtf ${ANNOTATION_GTF}"  

    # Check if the output file exists and is not empty
    if [ ! -s "$output_file" ]; then
        echo "Running command: $cmd"
        eval $cmd
    else
        echo "Output file $output_file already exists and is not empty."
    fi  

    ##################################################
    # Collapse step with GTF option
    output_file=$OUTPUT_DIR/${BASE_NAME}*.isoforms.fa

    cmd="flair collapse --reads ${BASE_NAME}.fq 
        -g $GENOME_FA 
        -q $OUTPUT_DIR/${BASE_NAME}_corrected_all_corrected.bed 
        -o $OUTPUT_DIR/${BASE_NAME}_collapsed 
        --gtf ${ANNOTATION_GTF}  
        --threads ${THREADS}
        --trust_ends"

        # Check if the output file exists and is not empty
    if [ ! -s "$output_file" ]; then
        echo $cmd
        eval $cmd
    else
        echo "Output file $output_file already exists and is not empty."
    fi 

done

#################################################
# THE ISOFORM need to be cat in to one file that cdhit at 1.0 to remove reducnacy here. 
output_file="${OUTPUT_DIR}/${EXP_NAME}.flair.fasta"
# make one isofome file
echo "concatanating your fa files to one transcriptome"
cat ${OUTPUT_DIR}/*.fa > ${OUTPUT_DIR}/${EXP_NAME}.FULL.fasta

# remove redundancy
echo "removing reduncdancy at 100% identity"
cmd="cd-hit-est 
    -i ${OUTPUT_DIR}/${EXP_NAME}.FULL.fasta
     -c 1.0 
     -o ${output_file} -M 0 -T ${THREADS}"

    # Check if the output file exists and is not empty
if [ ! -s "$output_file" ]; then
    echo "Running command: $cmd"
    eval $cmd
else
    echo "Output file $output_file already exists and is not empty."
fi
# remove the tmp file
rm ${OUTPUT_DIR}/${EXP_NAME}.FULL.fasta

# Quantify step -  you will have to prepare this manually. 
# manifest file : Tab delimited file containing sample id, condition, batch, reads.fq
# FROM THE DOCS: Note: Do not use underscores in the first three fields, see below for details.


#################################################
# make a final bed file:
output_file=${OUTPUT_DIR}/sorted_combined.bed
# Check if the output file exists and is not empty
if [ ! -s "$output_file" ]; then
    echo "combining bed file "
    cat ${OUTPUT_DIR}/*.isoforms.bed > ${OUTPUT_DIR}/combined.bed
    sort -k1,1 -k2,2n -k3,3n ${OUTPUT_DIR}/combined.bed > ${OUTPUT_DIR}/sorted_combined.bed
else
    echo "Output file $output_file already exists and is not empty."
fi


cat *.gtf > combined.gtf
sort -k1,1 -k4,4n combined.gtf > sorted.gtf
awk '!seen[$0]++' sorted.gtf > unique.gtf



#################################################
# Quantify
output_file=${OUTPUT_DIR}/${EXP_NAME}.tpm.tsv
cmd="flair quantify  
    -r manifest
    -o $OUTPUT_DIR/${EXP_NAME}_quantify 
    --threads ${THREADS} 
    --isoforms ${OUTPUT_DIR}/${EXP_NAME}.flair.fasta
    --tpm 
    --trust_ends
    --check_splice
    --isoform_bed ${OUTPUT_DIR}/sorted_combined.bed"

    # names dont match in the bed file. --isoform_bed ${OUTPUT_DIR}/sorted_combined.bed"

# Check if the output file exists and is not empty
if [ ! -s "$output_file" ]; then
    echo "Running command: " 
    eval $cmd
else
    echo "Output file $output_file already exists and is not empty."
fi



#################################################
# diffsplice

rm -rf ${OUTPUT_DIR}/${EXP_NAME}_flair_diffSplice
cmd="flair diffSplice 
    --threads  ${THREADS}
    --isoforms ${OUTPUT_DIR}/sorted_combined.bed
    --counts_matrix $OUTPUT_DIR/${EXP_NAME}_quantify.counts.tsv
    --test
    -o ${OUTPUT_DIR}/${EXP_NAME}_flair_diffSplice"

echo ${cmd}
eval ${cmd}


# filter for conditions of interest
awk 'BEGIN {FS=OFS="\t"} 
    NR==1 {
        for (i=1; i<=NF; i++) {
            if ($i == "ID") id_index = i
            else if ($i ~ /col0low|vir1low/) col[i]=1
        }
        printf("%s%s", $id_index, OFS);  # Print the ID header first
        for (i in col) printf("%s%s", $i, (i==NF) ? "\n" : OFS);  # Then print other headers
        next;  # Skip to the next record after printing the header
    }
    {
        printf("%s%s", $id_index, OFS);  # Print the ID data first
        for (i in col) {
            printf("%s%s", $i, (i==NF) ? "\n" : OFS);  # Then print other data
        }
    }'  $OUTPUT_DIR/${EXP_NAME}_quantify.counts.tsv > $OUTPUT_DIR/${EXP_NAME}vir_low_col0_low_quantify.counts.tsv



# expression is too low, so reduce thresholds --drim1 we dont have 6 samples anyway!

cmd="flair diffSplice 
    --threads  ${THREADS}
    --isoforms ${OUTPUT_DIR}/sorted_combined.bed
    --counts_matrix $OUTPUT_DIR/vir_low_col0_low.quantify.counts.tsv
    --test
    --drim1 2
    --drim3 10 
    --drim4 6 
    -o ${OUTPUT_DIR}/${EXP_NAME}_flair_diffSplice_vir_low_vs_col0_low"

echo ${cmd}
eval ${cmd}

echo " ... finished ... "




cmd="flair diffSplice 
    --threads  ${THREADS}
    --isoforms ${OUTPUT_DIR}/sorted_combined.bed
    --counts_matrix $OUTPUT_DIR/vir_low_col0_all.quantify.counts.tsv
    --test
    --drim1 2
    --drim3 10 
    --drim4 6 
    -o ${OUTPUT_DIR}/${EXP_NAME}_flair_diffSplice_vir_low_vs_col0_ALL"

echo ${cmd}
eval ${cmd}


$OUTPUT_DIR/${EXP_NAME}_quantify.counts.tsv > $OUTPUT_DIR/${EXP_NAME}reordered_quantify.counts.tsv


# YOU HAVE TO RENAME THEM as 2_Condition


cmd="flair diffSplice 
    --threads  ${THREADS}
    --isoforms ${OUTPUT_DIR}/sorted_combined.bed
    --counts_matrix alter_names.quantify.counts.tsv
    --test
    -o ${OUTPUT_DIR}/${EXP_NAME}_flair_diffSplice_SAMPLE_CONDITION"

echo ${cmd}
eval ${cmd}

