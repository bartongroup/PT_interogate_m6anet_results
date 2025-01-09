#!/bin/bash
#$ -pe smp 16
#$ -N virSQANTI3
#$ -cwd

# Abort on error
#set -e



# Script to run flair on a whole load of bam files. Reads already mapped to the genome
# Long reads!

# ##################################################
# INSTALL instructions ....
# https://github.com/ConesaLab/SQANTI3
#  git clone https://github.com/ConesaLab/SQANTI3.git
# cd SQANTI3
# conda env create -f SQANTI3.conda_env.yml
#conda activate SQANTI3.env

#mamba install -c gffread genometools

FIX_cupckae_import_error="
in the cupcake directory:
~/apps/SQANTI3/utilities/cupcake
add this import 
from utilities.cupcake.io import BioReaders  

to err_correct_w_genome.py "

# echo " you may need to do this: " ${FIX_cupckae_import_error}

##################################################
##################################################
# .... FILL THESE IN .....
# Define the path to the directory and other resources
##################################################

GENOME_FA="../genomic_data/TAIR10_chr_all.fas"
ANNOTATION_GTF="../genomic_data/atRTD3_TS_21Feb22_transfix.gtf"
THREADS=16
EXP_NAME="vir_temp"
SQANTI_PATH=~/apps/SQANTI3
DATA_DIR="${EXP_NAME}_atRTD3_flair_outputs"
ISOFORMs=${DATA_DIR}/${EXP_NAME}.flair.fasta
OUTPUT_DIR=${EXP_NAME}_SQANTI3

# ... NOTHING TO FILL IN FROM HERE
##################################################




mkdir -p $OUTPUT_DIR



# convert gff3 to gtf
#gffread transcript.fasta.sorted.gff3 -T -o transcript.fasta.sorted.gtf



# rmoved for now     --expression MATRIX.mtraix      --chunks 8
#  --expression ./vir_temp_flair_outputs/_quantify.counts.tsv
 # try the mapping phase?
cmd="python ${SQANTI_PATH}/sqanti3_qc.py 
    --dir ${OUTPUT_DIR}_1_skipORF_w_fasta
    --output ${EXP_NAME}
    --force_id_ignore
    --chunks 1
    --ratio_TSS_metric median
    --report both
    --aligner_choice minimap2
    --cpus ${THREADS}
     --expression ./vir_temp_flair_outputs/${EXP_NAME}_quantify.counts.tsv
    --skipORF
    --fasta ${ISOFORMs}
    ${ANNOTATION_GTF}
    ${GENOME_FA}"

echo ${cmd}
eval ${cmd}

info="
gt gff3 -sort  -sortlines -checkids -tidy -fixregionboundaries -addintrons -retainids ${DATA_DIR}/unique.gtf -o  ${DATA_DIR}sorted.gtf

/cluster/gjb_lab/pthorpe001/conda/envs/genometools/bin/bin/gt gff3 -sort \
 -sortlines -checkids -tidy -fixregionboundaries -addintrons -retainids \
 -o   ../genomic_data/atRTD3_TS_21Feb22_transfix.gt_fixed.gtf \
  ../genomic_data/atRTD3_TS_21Feb22_transfix.gtf



../genomic_data/atRTD3_TS_21Feb22_transfix.gtf

this is with all the gtf files 

cat vir1*.gtf > vir1combined.gtf
sort -k1,1 -k4,4n vir1combined.gtf > vir1sorted.gtf
awk '!seen[$0]++' vir1sorted.gtf > vir1unique.gtf

cat col0*.gtf > col0combined.gtf
sort -k1,1 -k4,4n col0combined.gtf > col0sorted.gtf
awk '!seen[$0]++' col0sorted.gtf > col0unique.gtf"

 # try the mapping phase?
cmd="python ${SQANTI_PATH}/sqanti3_qc.py 
    --dir ${OUTPUT_DIR}_2_All
    --output ${EXP_NAME}
    --force_id_ignore
    --chunks 1
    --ratio_TSS_metric median
    --report both
    --aligner_choice minimap2
    --cpus ${THREADS}
    --expression ./vir_temp_flair_outputs/${EXP_NAME}_quantify.counts.tsv
    ${DATA_DIR}/unique.gtf
    ${ANNOTATION_GTF}
    ${GENOME_FA}"

echo ${cmd}
eval ${cmd}

 # try the mapping phase?
cmd="python ${SQANTI_PATH}/sqanti3_qc.py 
    --dir ${OUTPUT_DIR}_2_ALL_skip_ORF
    --output ${EXP_NAME}
     --skipORF
    --force_id_ignore
    --chunks 1
    --ratio_TSS_metric median
    --report both
    --aligner_choice minimap2
    --cpus ${THREADS}
    --expression ./vir_temp_flair_outputs/${EXP_NAME}_quantify.counts.tsv
    ${DATA_DIR}/unique.gtf
    ${ANNOTATION_GTF}
    ${GENOME_FA}"

echo ${cmd}
eval ${cmd}

# vir
 # try the mapping phase?
cmd="python ${SQANTI_PATH}/sqanti3_qc.py 
    --dir ${OUTPUT_DIR}_vir
    --output ${EXP_NAME}
    --force_id_ignore
    --chunks 1
    --ratio_TSS_metric median
    --report both
    --aligner_choice minimap2
    --cpus ${THREADS}
    --expression ./vir_temp_flair_outputs/${EXP_NAME}_quantify.counts.tsv
    ${DATA_DIR}/vir1unique.gtf
    ${ANNOTATION_GTF}
    ${GENOME_FA}"

echo ${cmd}
eval ${cmd}

 # try the mapping phase?
cmd="python ${SQANTI_PATH}/sqanti3_qc.py 
    --dir ${OUTPUT_DIR}_vir_skip_orf
    --output ${EXP_NAME}
    --force_id_ignore
    --skipORF
    --chunks 1
    --ratio_TSS_metric median
    --report both
    --aligner_choice minimap2
    --cpus ${THREADS}
    --expression ./vir_temp_flair_outputs/${EXP_NAME}_quantify.counts.tsv
    ${DATA_DIR}/vir1unique.gtf
    ${ANNOTATION_GTF}
    ${GENOME_FA}"

echo ${cmd}
eval ${cmd}


 # try the mapping phase?
cmd="python ${SQANTI_PATH}/sqanti3_qc.py 
    --dir ${OUTPUT_DIR}_WT_skip_ORF
    --output ${EXP_NAME}
    --force_id_ignore
    --skipORF
    --chunks 1
    --ratio_TSS_metric median
    --report both
    --aligner_choice minimap2
    --cpus ${THREADS}
    --expression ./vir_temp_flair_outputs/${EXP_NAME}_quantify.counts.tsv
    ${DATA_DIR}/col0unique.gtf
    ${ANNOTATION_GTF}
    ${GENOME_FA}"

echo ${cmd}
eval ${cmd}


 # try the mapping phase?
cmd="python ${SQANTI_PATH}/sqanti3_qc.py 
    --dir ${OUTPUT_DIR}_WT_skip_ORF
    --output ${EXP_NAME}
    --force_id_ignore
    --chunks 1
    --ratio_TSS_metric median
    --report both
    --aligner_choice minimap2
    --cpus ${THREADS}
    --expression ./vir_temp_flair_outputs/${EXP_NAME}_quantify.counts.tsv
    ${DATA_DIR}/col0unique.gtf
    ${ANNOTATION_GTF}
    ${GENOME_FA}"

echo ${cmd}
eval ${cmd}




cmd="python ${SQANTI_PATH}/sqanti3_qc.py 
    --dir ${OUTPUT_DIR}_NO_EXPRESSION_1_skipORF_w_fasta
    --output ${EXP_NAME}
    --force_id_ignore
    --chunks 1
    --ratio_TSS_metric median
    --report both
    --aligner_choice minimap2
    --cpus ${THREADS}
    --skipORF
    --fasta ${ISOFORMs}
    ${ANNOTATION_GTF}
    ${GENOME_FA}"

echo ${cmd}
eval ${cmd}

info="
gt gff3 -sortlines -tidy -retainids ${DATA_DIR}/unique.gtf >  ${DATA_DIR}sorted.gtf

this is with all the gtf files 

cat vir1*.gtf > vir1combined.gtf
sort -k1,1 -k4,4n vir1combined.gtf > vir1sorted.gtf
awk '!seen[$0]++' vir1sorted.gtf > vir1unique.gtf

cat col0*.gtf > col0combined.gtf
sort -k1,1 -k4,4n col0combined.gtf > col0sorted.gtf
awk '!seen[$0]++' col0sorted.gtf > col0unique.gtf"

 # try the mapping phase?
cmd="python ${SQANTI_PATH}/sqanti3_qc.py 
    --dir ${OUTPUT_DIR}_NO_EXPRESSION_2_All
    --output ${EXP_NAME}
    --force_id_ignore
    --chunks 1
    --ratio_TSS_metric median
    --report both
    --aligner_choice minimap2
    --cpus ${THREADS}
    
    ${DATA_DIR}/unique.gtf
    ${ANNOTATION_GTF}
    ${GENOME_FA}"

echo ${cmd}
eval ${cmd}

 # try the mapping phase?
cmd="python ${SQANTI_PATH}/sqanti3_qc.py 
    --dir ${OUTPUT_DIR}_NO_EXPRESSION_2_ALL_skip_ORF
    --output ${EXP_NAME}
     --skipORF
    --force_id_ignore
    --chunks 1
    --ratio_TSS_metric median
    --report both
    --aligner_choice minimap2
    --cpus ${THREADS}
    
    ${DATA_DIR}/unique.gtf
    ${ANNOTATION_GTF}
    ${GENOME_FA}"

echo ${cmd}
eval ${cmd}

# vir
 # try the mapping phase?
cmd="python ${SQANTI_PATH}/sqanti3_qc.py 
    --dir ${OUTPUT_DIR}_NO_EXPRESSION_vir
    --output ${EXP_NAME}
    --force_id_ignore
    --chunks 1
    --ratio_TSS_metric median
    --report both
    --aligner_choice minimap2
    --cpus ${THREADS}
    
    ${DATA_DIR}/vir1unique.gtf
    ${ANNOTATION_GTF}
    ${GENOME_FA}"

echo ${cmd}
eval ${cmd}

 # try the mapping phase?
cmd="python ${SQANTI_PATH}/sqanti3_qc.py 
    --dir ${OUTPUT_DIR}_NO_EXPRESSION_vir_skip_orf
    --output ${EXP_NAME}
    --force_id_ignore
    --skipORF
    --chunks 1
    --ratio_TSS_metric median
    --report both
    --aligner_choice minimap2
    --cpus ${THREADS}
    
    ${DATA_DIR}/vir1unique.gtf
    ${ANNOTATION_GTF}
    ${GENOME_FA}"

echo ${cmd}
eval ${cmd}


 # try the mapping phase?
cmd="python ${SQANTI_PATH}/sqanti3_qc.py 
    --dir ${OUTPUT_DIR}_NO_EXPRESSION_WT_skip_ORF
    --output ${EXP_NAME}
    --force_id_ignore
    --skipORF
    --chunks 1
    --ratio_TSS_metric median
    --report both
    --aligner_choice minimap2
    --cpus ${THREADS}
    
    ${DATA_DIR}/col0unique.gtf
    ${ANNOTATION_GTF}
    ${GENOME_FA}"

echo ${cmd}
eval ${cmd}


 # try the mapping phase?
cmd="python ${SQANTI_PATH}/sqanti3_qc.py 
    --dir ${OUTPUT_DIR}_NO_EXPRESSION_WT_skip_ORF
    --output ${EXP_NAME}
    --force_id_ignore
    --chunks 1
    --ratio_TSS_metric median
    --report both
    --aligner_choice minimap2
    --cpus ${THREADS}
    
    ${DATA_DIR}/col0unique.gtf
    ${ANNOTATION_GTF}
    ${GENOME_FA}"

echo ${cmd}
eval ${cmd}
