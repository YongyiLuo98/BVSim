#!/bin/bash
#SBATCH -J hg38_bvsim_chr1-21_benchmark
#SBATCH -N 1 -c 32

#SBATCH --output=~/benchmark_summary.log

source ~/conda.sh
conda activate BVSim

BASE_OUTPUT_DIR="~/hg38"
BASE_LOG_DIR="~/hg38"
REGION_BED_DIR="~/BVSim/empirical/hg38"
SUMMARY_FILE="$BASE_LOG_DIR/benchmark_summary.csv"

echo "chromosome,seq_index,start_time,end_time,duration_sec,max_rss_mb,max_vm_mb" > "$SUMMARY_FILE"

for CHR in {1..22}; do
    CHR_NAME="chr$CHR"
    SEQ_INDEX=$((CHR - 1))
    
    # 动态路径
    OUTPUT_DIR="$BASE_OUTPUT_DIR/$CHR_NAME"
    LOG_OUT="$BASE_LOG_DIR/hg38_${CHR_NAME}_test.out"
    LOG_ERR="$BASE_LOG_DIR/hg38_${CHR_NAME}_test.err"
    mkdir -p "$OUTPUT_DIR"

    # 记录开始时间
    START_TIME=$(date +%s)
    echo "[$(date)] Starting $CHR_NAME with seq_index=$SEQ_INDEX..."

    /usr/bin/time -v bvsim \
        -mimic \
        -hg38 "$CHR_NAME" \
        -ref "/storage01/users/s1155146014/TGS/data/original_data/hg38/hg38_main_chroms.fa" \
        -seq_index "$SEQ_INDEX" \
        -save "$OUTPUT_DIR/" \
        -seed 0 \
        -rep "$SEQ_INDEX" \
        -cores 32 \
        -len_bins 500000 \
        -mode empirical \
        -snp 0.0001 \
        -snv_del 0.00001 \
        -snv_ins 0.00001 \
        -p_del_region 0.810 \
        -p_ins_region 0.828 \
        -region_bed_url "$REGION_BED_DIR/${CHR_NAME}_TR_unique.bed" \
        > "$LOG_OUT" 2> "$LOG_ERR"

    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))

    MAX_RSS=$(grep "Maximum resident set size" "$LOG_ERR" | awk '{print $6/1024}')  
    MAX_VM=$(grep "Virtual memory (kbytes)" "$LOG_ERR" | awk '{print $6/1024}')    

    
    echo "$CHR_NAME,$SEQ_INDEX,$START_TIME,$END_TIME,$DURATION,$MAX_RSS,$MAX_VM" >> "$SUMMARY_FILE"
    echo "[$(date)] Finished $CHR_NAME. Duration: $DURATION sec, Max RSS: $MAX_RSS MB, Max VM: $MAX_VM MB"
done

conda deactivate