#!/bin/bash

# 定义输入文件和输出路径
VCF="/disk18T3/project18/home_project18/data/original_data/TGS/HGSVC/freeze3.sv.alt.vcf.gz"
OUTPUT_DIR="/disk18T3/project18/data/test_data/TGS/HGSVC/"

# 创建输出目录（如果不存在）
mkdir -p "$OUTPUT_DIR"

# 获取样本列表
SAMPLES=$(bcftools query -l "$VCF")

# 为每个样本创建BED文件
for SAMPLE in $SAMPLES; do
    OUTPUT_FILE="${OUTPUT_DIR}${SAMPLE}.bed"
    
    # 方法1：使用--min-ac=1确保样本中存在变异
    bcftools view -s "$SAMPLE" --min-ac=1 "$VCF" | \
    bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\t%SVLEN\n' -i 'SVTYPE=="DEL" || SVTYPE=="INS"' | \
    awk '$5 >= 50 || $5 <= -50 {len=($5<0?-$5:$5); print $1"\t"$2"\t"$3"\t"$4"\t"len}' | \
    sort -k1,1 -k2,2n > "$OUTPUT_FILE"
    
    # 方法2：或者使用-c 1只保留样本为变异的记录（更严格）
    # bcftools view -s "$SAMPLE" -c 1 "$VCF" | \
    # bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\t%SVLEN\n' -i 'SVTYPE=="DEL" || SVTYPE=="INS"' | \
    # awk '$5 >= 50 || $5 <= -50 {len=($5<0?-$5:$5); print $1"\t"$2"\t"$3"\t"$4"\t"len}' | \
    # sort -k1,1 -k2,2n > "$OUTPUT_FILE"
    
    echo "Created BED file for $SAMPLE with $(wc -l < "$OUTPUT_FILE") records"
done

echo "All BED files have been generated in $OUTPUT_DIR"