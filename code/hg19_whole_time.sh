#!/bin/bash
#SBATCH -J bvsim_chr1-21_benchmark
#SBATCH -p chpc
#SBATCH -N 1 -c 32
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --exclusive
#SBATCH --nodelist=chpc-cn026
#SBATCH --output=/storage01/users/s1155146014/TGS/code/task03/benchmark_summary.log

# 加载环境
source /users/s1155146014/miniconda3/etc/profile.d/conda.sh
conda activate BVSim

# 基础路径
BASE_OUTPUT_DIR="/storage01/users/s1155146014/TGS/data/test_data/task03"
BASE_LOG_DIR="/storage01/users/s1155146014/TGS/code/task03"
REGION_BED_DIR="/storage01/users/s1155146014/BVSim/empirical/hg19"
SUMMARY_FILE="$BASE_LOG_DIR/benchmark_summary.csv"

# 初始化汇总文件（CSV格式）
echo "chromosome,seq_index,start_time,end_time,duration_sec,max_rss_mb,max_vm_mb" > "$SUMMARY_FILE"

# 循环处理 chr1 到 chr21
for CHR in {1..21}; do
    CHR_NAME="chr$CHR"
    SEQ_INDEX=$((CHR - 1))
    
    # 动态路径
    OUTPUT_DIR="$BASE_OUTPUT_DIR/$CHR_NAME"
    LOG_OUT="$BASE_LOG_DIR/hg19_${CHR_NAME}_test.out"
    LOG_ERR="$BASE_LOG_DIR/hg19_${CHR_NAME}_test.err"
    mkdir -p "$OUTPUT_DIR"

    # 记录开始时间
    START_TIME=$(date +%s)
    echo "[$(date)] Starting $CHR_NAME with seq_index=$SEQ_INDEX..."

    # 运行 BVSim 并通过 /usr/bin/time 监控资源
    /usr/bin/time -v bvsim \
        -ref "/storage01/users/s1155146014/TGS/data/original_data/hg19/hs37d5.fasta" \
        -seq_index "$SEQ_INDEX" \
        -save "$OUTPUT_DIR/" \
        -seed 0 \
        -rep "$SEQ_INDEX" \
        -cores 32 \
        -len_bins 500000 \
        -hg19 "$CHR_NAME" \
        -mode empirical \
        -snp 0.0001 \
        -snv_del 0.00001 \
        -snv_ins 0.00001 \
        -p_del_region 0.810 \
        -p_ins_region 0.828 \
        -region_bed_url "$REGION_BED_DIR/${CHR_NAME}_TR_unique.bed" \
        > "$LOG_OUT" 2> "$LOG_ERR"

    # 记录结束时间和资源使用
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    
    # 从 time 输出中提取内存指标（Max RSS 和 Max VM）
    MAX_RSS=$(grep "Maximum resident set size" "$LOG_ERR" | awk '{print $6/1024}')  # 转换为MB
    MAX_VM=$(grep "Virtual memory (kbytes)" "$LOG_ERR" | awk '{print $6/1024}')    # 转换为MB

    # 写入汇总文件
    echo "$CHR_NAME,$SEQ_INDEX,$START_TIME,$END_TIME,$DURATION,$MAX_RSS,$MAX_VM" >> "$SUMMARY_FILE"
    echo "[$(date)] Finished $CHR_NAME. Duration: $DURATION sec, Max RSS: $MAX_RSS MB, Max VM: $MAX_VM MB"
done

conda deactivate