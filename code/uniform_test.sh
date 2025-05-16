#!/bin/bash
#SBATCH -J uniform_benchmark
#SBATCH -N 1 -c 1
#SBATCH --output=uniform_benchmark_%j.out
#SBATCH --error=uniform_benchmark_%j.err

# 加载conda环境
source ~/conda.sh
conda activate BVSim

# 创建结果目录
RESULT_BASE="~/benchmark_results"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RESULT_DIR="${RESULT_BASE}/uniform_${TIMESTAMP}"
mkdir -p ${RESULT_DIR}

for i in {1..10}; do
    echo "Running replicate $i"
    
    RUN_DIR="${RESULT_DIR}/run_${i}"
    mkdir -p ${RUN_DIR}
    
    /usr/bin/time -f "elapsed_time %e\nmax_rss %M\nswap %W" -o ${RUN_DIR}/time.txt \
    bvsim \
        -uniform \
        -ref /storage01/users/s1155146014/BVSim/empirical/sub_hg19_chr1.fasta \
        -seq_index 0 \
        -save ${RUN_DIR} \
        -seed ${i} \
        -rep ${i} \
        -snp 100 \
        -snv_del 10 \
        -snv_ins 10 \
        2>&1 | tee ${RUN_DIR}/log.txt
done


echo "Generating summary report..."
echo "Replicate,Time(sec),MaxMemory(kb),Swap(kb)" > ${RESULT_DIR}/summary.csv
for i in {1..10}; do
    RUN_DIR="${RESULT_DIR}/run_${i}"
    time_sec=$(grep "elapsed_time" ${RUN_DIR}/time.txt | awk '{print $2}')
    max_mem=$(grep "max_rss" ${RUN_DIR}/time.txt | awk '{print $2}')
    swap=$(grep "swap" ${RUN_DIR}/time.txt | awk '{print $2}')
    echo "$i,$time_sec,$max_mem,$swap" >> ${RESULT_DIR}/summary.csv
done

awk -F',' '
NR>1 {
    time[NR]=$2; sum_t+=$2; sum_t2+=$2*$2;
    mem[NR]=$3; sum_m+=$3; sum_m2+=$3*$3;
    swap[NR]=$4; sum_s+=$4; sum_s2+=$4*$4;
} 
END {
    mean_t=sum_t/10; mean_m=sum_m/10; mean_s=sum_s/10;
    sd_t=sqrt(sum_t2/10 - mean_t^2);
    sd_m=sqrt(sum_m2/10 - mean_m^2);
    sd_s=sqrt(sum_s2/10 - mean_s^2);
    print "Average,"mean_t","mean_m","mean_s;
    print "SD,"sd_t","sd_m","sd_s;
}' ${RESULT_DIR}/summary.csv >> ${RESULT_DIR}/summary.csv

echo "Benchmark completed. Results saved to ${RESULT_DIR}"
conda deactivate