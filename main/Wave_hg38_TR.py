#!/usr/bin/env python
# coding: utf-8
import time
start_time0 = time.time()
# process_time0 = time.process_time()
#packages
import numpy as np
import pandas as pd
import random
import pysam
import sys
import re
from random import sample, choice
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import copy
from copy import deepcopy
import os
from datetime import datetime
#from multiprocessing import Pool
import multiprocessing.pool as pool

from scipy.stats import pareto
import time
from datetime import timedelta
# Import necessary libraries
from pandas import concat
import seaborn as sns
from scipy.stats import entropy

import psutil
import time
import os
import multiprocessing
import pickle
import shutil
from multiprocessing import Pool, Manager

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import gc
import multiprocessing.shared_memory as shm
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)
from typing import List, Tuple, Union
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
print('hg38 mode')


def load_reference_sequence(fasta_path: str, seq_index: int) -> Tuple[str, str]:
    """直接读取解压后的FASTA文件（不再处理.gz）"""
    try:
        with pysam.FastaFile(fasta_path) as fasta_file:
            if seq_index >= len(fasta_file.references):
                raise ValueError(f"Sequence index {seq_index} out of range")
            seqname = fasta_file.references[seq_index]
            return seqname, fasta_file.fetch(seqname)
    except Exception as e:
        raise RuntimeError(f"Failed to load reference sequence: {str(e)}")


def plot_sv_distribution_in_tr_region(args, df_sv, sv_type, chr_id, chr_length, len_bins, df_TR, output_dir=None, rep=None):
    """
    绘制SV类型在TR区域内和区域外的分布 (模拟数据)
    
    参数:
    - df_sv: 包含SV信息的DataFrame (模拟数据)
    - sv_type: 要绘制的SV类型 ('Deletion'或'Insertion')
    - chr_id: 染色体ID
    - chr_length: 染色体长度
    - len_bins: 分箱大小
    - df_TR: 包含TR区域的DataFrame，必须包含'start'和'end'列
    - output_dir: 输出目录 (可选)
    - rep: 重复次数标识 (可选)
    """
    # 筛选指定类型的SV
    df = df_sv[df_sv['SV_type'] == sv_type].copy()
    if df.empty:
        print(f"No {sv_type} data for {chr_id}")
        return
    
    # 根据输入的分箱方式重新构建bins
    devided_length = int(chr_length // len_bins)
    remainder_length = chr_length - devided_length * len_bins
    bins = np.arange(0, chr_length - remainder_length, len_bins).tolist()
    bins.append(chr_length)
    
    # 使用set.union()函数合并所有TR区域的点，得到一个包含所有点的集合
    df_TR['points'] = df_TR.apply(get_points, axis=1)
    TR_set = set.union(*df_TR['points'])

    # 将模拟数据分为TR区域内和区域外
    df['is_in_TR'] = df['Original_start'].apply(lambda x: x in TR_set)
    sv_counts_in_TR = df[df['is_in_TR']].groupby(pd.cut(df[df['is_in_TR']]['Original_start'], bins)).size()
    sv_counts_out_TR = df[~df['is_in_TR']].groupby(pd.cut(df[~df['is_in_TR']]['Original_start'], bins)).size()

    # 确保计数长度与bins一致
    sv_counts_in_TR = sv_counts_in_TR.reindex(pd.IntervalIndex.from_breaks(bins), fill_value=0)
    sv_counts_out_TR = sv_counts_out_TR.reindex(pd.IntervalIndex.from_breaks(bins), fill_value=0)

    # 计算TR区域内的SV比例
    total_svs = len(df)
    svs_in_TR = df['is_in_TR'].sum()
    p_in_TR = svs_in_TR / total_svs if total_svs > 0 else 0

    # 创建图形
    plt.figure(figsize=(10, 6))
    
    # 计算X轴位置和范围
    x_max = chr_length / 1e6  # 转换为Mbp单位
    x = np.array([(bins[i] + bins[i + 1]) / 2 for i in range(len(bins) - 1)]) / 1e6  # 将点放在bin中心

    # 绘制TR区域内和区域外的分布
    plt.plot(x, sv_counts_in_TR.values, 'o-', color='green', label=f'{sv_type} in TR Region')
    plt.plot(x, sv_counts_out_TR.values, 'o-', color='orange', label=f'{sv_type} outside TR Region')
    
    # 设置图形属性
    plt.title(f'Simulated {sv_type} Distribution in TR Region vs Non-TR Region on Sequence {chr_id}')
    plt.xlabel(f'Genomic Position (Mbp, bins = {len_bins/1e3:.0f}kb)')
    plt.ylabel('Number of SVs')
    plt.xlim(0, x_max)
    xticks = np.arange(0, x_max + 1, 5)  # 每5Mbp一个刻度
    plt.xticks(xticks, rotation=45, ha='right')
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # 添加图例，包含比例信息
    green_line = plt.Line2D([0], [0], color='green', marker='o', lw=2, 
                          label=f'{sv_type} in TR Region ({p_in_TR:.1%})')
    orange_line = plt.Line2D([0], [0], color='orange', marker='o', lw=2, 
                           label=f'{sv_type} outside TR Region ({1-p_in_TR:.1%})')
    plt.legend(handles=[green_line, orange_line], loc='upper right')
    
    # 调整布局防止标签被截断
    plt.tight_layout()
    
    # 创建输出目录(如果不存在)并保存图形
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        # 修改输出路径格式
        output_path = os.path.join(output_dir, 'BV_' + str(args.rep) + '_seq_'+ f"{chr_id}_{sv_type}_distribution_comparison_TR.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved {sv_type} TR region distribution comparison plot to {output_path}")
    plt.close()

def write_template_fasta_con(args, seqname, consensus_):
    # 生成输出文件路径
    output_path = args.save + 'BV_' + str(args.rep) + "_seq_"+str(seqname) +".fasta"
    # Prepare the new sequence
    sequences = [consensus_]
    new_sequences = []
    for sequence in sequences:
        record = SeqRecord(Seq(re.sub('[^GATCN-]', "", str(sequence).upper())), id=seqname, name=seqname, description="<custom description>")
        new_sequences.append(record)

    # Write the new sequence to a file
    # with open(args.save + 'BV_' + str(args.rep) + "_seq_"+str(seqname) +".fasta", "w") as output_handle:
    #     SeqIO.write(new_sequences, output_handle, "fasta")
    # Write the new sequence to a file
    with open(output_path, "w") as output_handle:
        SeqIO.write(new_sequences, output_handle, "fasta")
    
    # 打印输出文件路径
    print(f"Saved consensus sequence to {output_path}")

def write_vcf(args, df, seqname, start_base, end_base):
    # 生成输出文件路径
    output_path = args.save + 'BV_' + str(args.rep) + '_seq_' + str(seqname) + ".vcf"
    
    # Get the current date
    current_date = datetime.now().strftime('%Y%m%d')
    
    # Write the DataFrame to a VCF file
    with open(output_path, 'w') as f:
    # with open(args.save + 'BV_' + str(args.rep) + '_seq_' + str(seqname) + ".vcf", 'w') as f:
        # VCF header
        f.write('##fileformat=VCFv4.3\n')  # 更新为 VCF 4.3
        f.write('##fileDate=' + current_date + '\n')
        f.write('##source=Wave_hg38_TR.py\n')
        f.write('##reference=' + args.ref + ':' + str(start_base) + '-' + str(end_base) + '\n')
        f.write('##contig=<ID=' + str(seqname) + ',length=' + str(end_base - start_base + 1) + '>\n')
        
        # INFO fields
        f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        f.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the variant">\n')
        f.write('##INFO=<ID=CSV_TYPE,Number=1,Type=String,Description="Type of CSV">\n')
        f.write('##INFO=<ID=CSV_INDEX,Number=1,Type=Integer,Description="Index of CSV">\n')
        f.write('##INFO=<ID=SUB,Number=1,Type=String,Description="Substitution">\n')
        f.write('##INFO=<ID=smallINS,Number=1,Type=String,Description="Small insertion">\n')
        f.write('##INFO=<ID=LEN,Number=1,Type=Integer,Description="Length of the variant">\n')
        f.write('##INFO=<ID=smallDEL,Number=1,Type=String,Description="Small deletion">\n')
        f.write('##INFO=<ID=DEL,Number=1,Type=String,Description="Deletion">\n')
        f.write('##INFO=<ID=INS,Number=1,Type=String,Description="Insertion">\n')
        f.write('##INFO=<ID=INV,Number=1,Type=String,Description="Inversion">\n')
        
        # FORMAT field
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        
        # Column headers
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_ID\n')
        
        # Write data
        if df.empty:
            print("Warning: DataFrame is empty. No data will be written to VCF.")
        else:
            df.to_csv(f, sep='\t', index=False, header=False)
            
    # 打印输出文件路径
    print(f"Saved VCF file to {output_path}")

def normalize_mode(mode_input: str) -> str:
    """Normalize the mode input and handle typos/errors."""
    mode_input = str(mode_input).strip().lower()
    if mode_input.startswith('prob'):
        return 'probability'
    elif mode_input.startswith('emp'):
        return 'empirical'
    else:
        print(f"Warning: Unknown mode '{mode_input}'. Defaulting to 'empirical'.")
        return 'empirical'
    

def plot_sv_distribution(mode, args, df_sv, sv_type, chr_id, chr_length, len_bins, empirical_data=None, n_bed_path=None, output_dir=None):
    """
    绘制SV类型在染色体上的分布，并比较模拟数据和实际数据
    
    参数:
    - df_sv: 包含SV信息的DataFrame (模拟数据)
    - sv_type: 要绘制的SV类型 ('Deletion'或'Insertion')
    - chr_id: 染色体ID
    - chr_length: 染色体长度
    - len_bins: 分箱大小 (来自args.len_bins)
    - empirical_data: 实际数据 (单个样本的计数列表或多样本的均值和标准差字典)
    - n_bed_path: N区域BED文件路径 (可选)
    - output_dir: 输出目录 (来自args.save)
    """
    # 筛选指定类型的SV
    df = df_sv[df_sv['SV_type'] == sv_type].copy()
    if df.empty:
        print(f"No {sv_type} data for {chr_id}")
        return
    
    # 根据输入的分箱方式重新构建bins
    devided_length = int(chr_length // len_bins)
    remainder_length = chr_length - devided_length * len_bins
    bins = np.arange(0, chr_length - remainder_length, len_bins).tolist()
    bins.append(chr_length)
    
    # 使用新的bins对模拟数据统计每个bin的SV数量
    df['dist_bin'] = pd.cut(df['Original_start'], bins)
    sv_counts = df.groupby('dist_bin').size()
    
    # 确保模拟数据的计数长度与bins一致
    sv_counts = sv_counts.reindex(pd.IntervalIndex.from_breaks(bins), fill_value=0)
    
    # 创建图形
    plt.figure(figsize=(10, 6))
    
    # 计算X轴位置和范围
    x_max = chr_length / 1e6  # 转换为Mbp单位
    x = np.array([(bins[i] + bins[i + 1]) / 2 for i in range(len(bins) - 1)]) / 1e6  # 将点放在bin中心
    
    # 绘制模拟数据分布
    sim_line = plt.plot(x, sv_counts.values, 'o-', color='blue', label=f'Simulated {sv_type} Count')[0]
    
    # 绘制实际数据分布
    if empirical_data is not None:
        if isinstance(empirical_data, dict):
            # 多样本情况 - 绘制均值和标准差范围
            plt.plot(x, empirical_data['mean'], 'o-', color='green', label=f'Samples {sv_type} Mean')
            # plt.fill_between(x, 
            #                  empirical_data['mean'] - empirical_data['std'], 
            #                  empirical_data['mean'] + empirical_data['std'], 
            #                  color='red', alpha=0.2, label=f'Empirical {sv_type} ±1 SD')
            # 计算95%置信区间
            confidence_interval_lower = empirical_data['mean'] - 1.96 * empirical_data['std']
            confidence_interval_upper = empirical_data['mean'] + 1.96 * empirical_data['std']
            
            # 绘制95%置信区间
            plt.fill_between(x, confidence_interval_lower, confidence_interval_upper, 
                             color='green', alpha=0.5, label='95% Confidence Interval')
        else:
            # 单样本情况 - 直接绘制计数
            plt.plot(x, empirical_data, 'o-', color='red', label=f'Sample {sv_type} Count')
    
    # 如果有N区域数据，绘制灰色背景
    if n_bed_path and os.path.exists(n_bed_path):
        df_N = pd.read_csv(n_bed_path, sep='\t', names=['chr', 'start', 'end'])
        for _, row in df_N.iterrows():
            start_norm = row['start'] / 1e6  # 转换为Mbp
            end_norm = row['end'] / 1e6     # 转换为Mbp
            plt.axvspan(start_norm, end_norm, color='gray', alpha=0.3)
    
    # 设置图形属性
    # plt.title(f'{sv_type} Distribution Comparison on {chr_id} (bin size: {len_bins/1e3:.0f}kb)')
    # plt.xlabel('Genomic Position (Mbp)')

    # 确定样本类型
    if empirical_data is not None:
        sample_type = "Multiple Samples" if isinstance(empirical_data, dict) else "Single Sample"
    else:
        sample_type = "No Empirical Data"

    # 根据模式生成模式描述
    if mode == 'probability':
        mode_desc = "Probability-based Simulation"
    elif mode == 'empirical':
        mode_desc = "Sample-matched Distribution"  # 或 "Empirical-matched Distribution"
    else:
        mode_desc = f"Invalid Mode ({mode}), showing Empirical"

    # 组合标题 - 选择你喜欢的格式
    title_v1 = f"{mode_desc} of {sv_type} SVs for {chr_id} from {sample_type}"
    title_v2 = f"{sv_type} SV distribution in {chr_id} using {mode_desc} ({sample_type})"
    title_v3 = f"{mode_desc} of {sv_type} for {chr_id} based on {sample_type}"

    plt.title(title_v2)  # 选择最合适的版本
    # plt.title(f'Simulated {sv_type} Distribution on Sequence {chr_id}')
    plt.xlabel(f'Genomic Position (Mbp, bins = {len_bins/1e3:.0f}kb)')
    plt.ylabel('Number of SVs')
    
    plt.xlim(0, x_max)
    xticks = np.arange(0, x_max + 1, 5)  # 每5Mbp一个刻度
    plt.xticks(xticks, rotation=90, ha='right')
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # 添加图例
    handles = [sim_line]
    if empirical_data is not None:
        if isinstance(empirical_data, dict):
            handles.append(plt.Line2D([], [], color='green', marker='o', label=f'Samples {sv_type} Mean'))
            handles.append(mpatches.Patch(color='green', alpha=0.5, label='95% Confidence Interval'))
        else:
            handles.append(plt.Line2D([], [], color='red', marker='o', label=f'Sample {sv_type} Count'))
    
    if n_bed_path and os.path.exists(n_bed_path):
        gray_patch = mpatches.Patch(color='gray', alpha=0.3, label='N Region')
        handles.append(gray_patch)
    
    plt.legend(handles=handles, loc='upper right')
    
    # 调整布局防止标签被截断
    plt.tight_layout()
    
    # 创建输出目录(如果不存在)并保存图形
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'BV_' + str(args.rep) + '_seq_'+ f"{chr_id}_{sv_type}_distribution_comparison.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved {sv_type} distribution comparison plot to {output_path}")

# 将N positions转换为BED格式并保存
def create_N_bed(n_positions, chr_id, output_dir):
    # 将连续的N positions合并为区间
    n_regions = []
    if len(n_positions) > 0:
        n_positions = sorted(n_positions)
        start = n_positions[0]
        prev = start
        
        for pos in n_positions[1:]:
            if pos != prev + 1:
                n_regions.append((start, prev))
                start = pos
            prev = pos
        n_regions.append((start, prev))
    
    # 创建DataFrame并保存为BED文件
    bed_data = [(chr_id, start, end) for start, end in n_regions]
    df_N = pd.DataFrame(bed_data, columns=['chr', 'start', 'end'])
    n_bed_path = os.path.join(output_dir, f"{chr_id}_N.bed")
    df_N.to_csv(n_bed_path, sep='\t', header=False, index=False)
    return n_bed_path


def validate_empirical_data(empirical_data: List[float], 
                          data_type: str, 
                          expected_length: int) -> bool:
    """Validate the empirical data input."""
    if len(empirical_data) != expected_length:
        print(f'Warning: empirical {data_type} info has wrong dimension '
              f'(expected {expected_length}, got {len(empirical_data)}).')
        return False
    
    if any(p < 0 for p in empirical_data):
        print(f'Warning: empirical {data_type} contains negative values.')
        return False
    
    return True
def calculate_probabilities(empirical_data: List[float], 
                          data_type: str) -> List[float]:
    """Calculate probabilities from empirical data with safety checks."""
    total = sum(empirical_data)
    if total <= 0:
        print(f'Warning: Sum of {data_type} empirical data is {total}. '
              'Using uniform probabilities.')
        return [1.0/len(empirical_data)] * len(empirical_data)
    return [p/total for p in empirical_data]

def process_sv_data(empirical_data: List[float],
                   sv_total_arg: int,
                   data_type: str,
                   number_seg: int,
                   args_sum: bool,
                   mode: str) -> Tuple[List[int], int]:
    """
    Process structural variant data with comprehensive error handling.
    Returns integer counts and total.
    """
    if not validate_empirical_data(empirical_data, data_type, number_seg):
        return [0] * number_seg, 0
    
    try:
        if mode == 'probability':
            probabilities = calculate_probabilities(empirical_data, data_type)
            
            total_sv = max(sum(empirical_data), sv_total_arg) if args_sum else sv_total_arg
            
            if total_sv <= 0:
                print(f'Warning: Total {data_type.upper()} SV is {total_sv}. '
                      'Using empirical values instead.')
                # Convert empirical data to integers
                empirical_ints = [int(round(x)) for x in empirical_data]
                total_sv = sum(empirical_ints)
                return empirical_ints, total_sv
            
            print(f'Total {data_type.upper()}: {total_sv}')
            sv_per_segment = np.random.multinomial(total_sv, probabilities)
            # return list(sv_per_segment), total_sv
            return [int(x) for x in sv_per_segment], int(total_sv)
        elif mode == 'empirical':
            # Round empirical data to integers
            empirical_ints = [int(round(x)) for x in empirical_data]
            total_sv = sum(empirical_ints)
            print(f'Total {data_type.upper()}: {total_sv}')
            # return empirical_ints, total_sv
            return [int(x) for x in empirical_ints], total_sv
        
        else:
            print(f'Error: Invalid mode {mode} for {data_type}. Using empirical values.')
            empirical_ints = [int(round(x)) for x in empirical_data]
            total_sv = sum(empirical_ints)
            # return empirical_ints, total_sv
            return [int(x) for x in empirical_ints], total_sv
    except Exception as e:
        print(f'Error processing {data_type} data: {str(e)}. Using empirical values.')
        empirical_ints = [int(round(x)) for x in empirical_data]
        total_sv = sum(empirical_ints)
        # return empirical_ints, total_sv
        return [int(x) for x in empirical_ints], total_sv

# 定义 SV_type 到缩写的映射规则
def get_svtype_abbreviation(sv_type):
    if sv_type == 'Small_Del':
        return 'smallDEL'
    elif sv_type == 'Small_Ins':
        return 'smallINS'
    else:
        return sv_type[:3].upper()  # 其他类型取前三个字母并转为大写

# 生成 INFO 列
def generate_info_column(sv_types, sv_lens):
    info_columns = []
    for sv_type, sv_len in zip(sv_types, sv_lens):
        if sv_type == 'Substitution':
            info_columns.append('SUB')
        elif sv_type == 'Small_Del':
            info_columns.append(f'smallDEL;LEN={sv_len}')
        elif sv_type == 'Small_Ins':
            info_columns.append(f'smallINS;LEN={sv_len}')
        elif sv_type == 'Deletion':
            info_columns.append(f'SVTYPE=DEL;SVLEN={sv_len}')
        elif sv_type == 'Insertion':
            info_columns.append(f'SVTYPE=INS;SVLEN={sv_len}')
        elif sv_type == 'Inversion':
            info_columns.append(f'SVTYPE=INV;SVLEN={sv_len}')  # 明确添加INV类型
        else:
            info_columns.append('UNKNOWN')
    return np.array(info_columns)

def check_start_end(row):
    try:
        start = int(row['start'])
        end = int(row['end'])
    except ValueError:
        print('Error: "start" or "end" value is not an integer.')
        return False
    
    if start >= end:
        print('Warning: The "start" value is greater than or equal to the "end" value.')
        return False
    return True

# Define a function that receives a row of data and returns all points in the 'start' and 'end' regions of the row. End is not included.
def get_points(row):
    return set(range(row['start'], row['end']))
   
    
def DNA_complement(sequence):
    trantab = str.maketrans('ATCGatcg','TAGCtagc')
    string = sequence.translate(trantab)
    return string
 
def error_handler(e):
    print('Error occurred:', e)

# 定义 SV_type 到缩写的映射规则
def SV_write_relative(SV_table_merged, ll_c, tem_ins_dic):
    tem_SV_table_merged = SV_table_merged[SV_table_merged.iloc[:, 1] == ll_c]
    # original start: start of A
    list_start1 = list(tem_SV_table_merged.iloc[:, 3])
    # start of B (dulplication or balanced trans)
    list_start2 = list(tem_SV_table_merged.iloc[:, 6])
    whole_start_abs = list_start1 + list_start2
    # order SV from left
    # set:merge repeated sites (e.g. 5 mismatch 5.5 ins)
    # ! 筛选出大于-1的
    whole_start_abs = [item for item in whole_start_abs if item > 0]
    whole_start_abs_set = sorted(list(set(whole_start_abs)))
    present_len = 0
    last_bone = 0
    # inserted term

    for ll_var_index in range(len(whole_start_abs_set)):

        # ! time 找到对应的行
        tem_SV_table_merged2 = tem_SV_table_merged[
            (tem_SV_table_merged['Original_start'] == whole_start_abs_set[ll_var_index]) | \
            (tem_SV_table_merged['New_start'] == whole_start_abs_set[ll_var_index])
        ]

        for xun_nei_row in range(len(tem_SV_table_merged2)):
            tem_row = tem_SV_table_merged2.iloc[xun_nei_row, :]  # 获取行数据
            stand_line = int(tem_row.iloc[0])  # 使用.iloc获取位置索引
            # A
            bone1s = tem_row.iloc[3]
            bone1e = tem_row.iloc[4]
            # B
            bone2s = tem_row.iloc[6]
            bone2e = tem_row.iloc[7]
            if whole_start_abs_set[ll_var_index] in list_start1:
                # 处理不同的SV类型
                sv_type = tem_row.iloc[2]  # 使用.iloc获取SV类型
                if sv_type in ['Substitution', 'Small_Ins', 'Small_Del', 'Deletion', 'Insertion', 'Inversion']:
                    if sv_type in ['Deletion', 'Small_Del']:
                        inster_number_bone = bone1s - last_bone - 1
                        present_len += inster_number_bone
                        last_bone = bone1e
                        SV_table_merged.iloc[stand_line, 10] = -1
                        SV_table_merged.iloc[stand_line, 11] = -1
                        SV_table_merged.iloc[stand_line, 12] = -1
                        SV_table_merged.iloc[stand_line, 13] = -1
                    elif sv_type == 'Substitution':
                        inster_number_bone = bone1s - last_bone
                        present_len += inster_number_bone
                        last_bone = bone1e
                        SV_table_merged.iloc[stand_line, 10] = present_len
                        SV_table_merged.iloc[stand_line, 11] = present_len
                        SV_table_merged.iloc[stand_line, 12] = -1
                        SV_table_merged.iloc[stand_line, 13] = -1
                    elif sv_type in ['Small_Ins', 'Insertion']:
                        inster_number_bone = bone1s - last_bone
                        Ins_len_present = len(tem_ins_dic[bone1s])
                        SV_table_merged.iloc[stand_line, 10] = present_len + inster_number_bone + 1
                        SV_table_merged.iloc[stand_line, 12] = -1
                        present_len += inster_number_bone + Ins_len_present
                        last_bone = bone1e
                        SV_table_merged.iloc[stand_line, 11] = present_len
                        SV_table_merged.iloc[stand_line, 13] = -1
                    else:  # Inversion
                        inster_number_bone = bone1s - last_bone
                        SV_table_merged.iloc[stand_line, 10] = present_len + inster_number_bone
                        SV_table_merged.iloc[stand_line, 12] = -1
                        present_len += bone1e - last_bone
                        SV_table_merged.iloc[stand_line, 11] = present_len
                        SV_table_merged.iloc[stand_line, 13] = -1
                        last_bone = bone1e
                elif sv_type == 'Duplication':
                    inster_number_bone = bone1s - last_bone
                    tem_plate_len = tem_row.iloc[5]
                    SV_table_merged.iloc[stand_line, 10] = present_len + inster_number_bone
                    present_len += inster_number_bone + tem_plate_len - 1
                    SV_table_merged.iloc[stand_line, 11] = present_len
                    last_bone = bone1e
                elif sv_type == 'Translocation':
                    if tem_row.iloc[9] == 1:
                        inster_number_bone = bone1s - last_bone - 1
                        SV_table_merged.iloc[stand_line, 10] = present_len + inster_number_bone + 1
                        Ins_len_present = len(tem_ins_dic[bone1s - 1])
                        present_len += inster_number_bone + Ins_len_present
                        SV_table_merged.iloc[stand_line, 11] = present_len
                        last_bone = bone1e
                    else:
                        inster_number_bone = bone1s - last_bone - 1
                        present_len += inster_number_bone
                        SV_table_merged.iloc[stand_line, 10] = present_len + 1
                        SV_table_merged.iloc[stand_line, 11] = present_len + 1
                        last_bone = bone1e

            else:  # 处理B部分
                sv_type = tem_row.iloc[2]
                if sv_type == 'Duplication':
                    inster_number_bone = bone2s - last_bone
                    SV_table_merged.iloc[stand_line, 12] = present_len + inster_number_bone + 1
                    Ins_len_present = len(tem_ins_dic[bone2s])
                    present_len += inster_number_bone + Ins_len_present
                    SV_table_merged.iloc[stand_line, 13] = present_len
                    last_bone = bone2e
                elif sv_type == 'Translocation':
                    if tem_row.iloc[9] == 1:
                        inster_number_bone = bone2s - last_bone - 1
                        SV_table_merged.iloc[stand_line, 12] = present_len + inster_number_bone + 1
                        Ins_len_present = len(tem_ins_dic[bone2s - 1])
                        present_len += inster_number_bone + Ins_len_present
                        SV_table_merged.iloc[stand_line, 13] = present_len
                        last_bone = bone2e
                    else:
                        inster_number_bone = bone2s - last_bone
                        SV_table_merged.iloc[stand_line, 12] = present_len + inster_number_bone + 1
                        Ins_len_present = len(tem_ins_dic[bone2s])
                        present_len += inster_number_bone + Ins_len_present
                        SV_table_merged.iloc[stand_line, 13] = present_len
                        SV_table_merged.iloc[stand_line, 10] = -1
                        SV_table_merged.iloc[stand_line, 11] = -1
                        last_bone = bone2e
    return SV_table_merged
  
# def SV_write_relative(SV_table_merged,ll_c,tem_ins_dic):
#     tem_SV_table_merged = SV_table_merged[SV_table_merged.iloc[:,1]==ll_c]
#     #original start: start of A
#     list_start1 = list(tem_SV_table_merged.iloc[:,3])
#     #start of B (dulplication or balanced trans)
#     list_start2 = list(tem_SV_table_merged.iloc[:,6])
#     whole_start_abs = list_start1+list_start2
#     #order SV from left
#     #set:merge repeated sites (e.g. 5 mismatch 5.5 ins)
#     #! 筛选出大于-1的
#     whole_start_abs = [item for item in whole_start_abs if item > 0]
#     whole_start_abs_set = sorted(list(set(whole_start_abs)))
#     present_len = 0
#     last_bone = 0
#     #inserted term
    
#     for ll_var_index in range(len(whole_start_abs_set)):
        
#         #! time 找到对应的行
#         tem_SV_table_merged2 = tem_SV_table_merged[(tem_SV_table_merged['Original_start']==whole_start_abs_set[ll_var_index]) |\
#                                         (tem_SV_table_merged['New_start']==whole_start_abs_set[ll_var_index])]

#         for xun_nei_row in range(len(tem_SV_table_merged2)):
#             tem_row = tem_SV_table_merged2.iloc[xun_nei_row,:]
#             stand_line = int(tem_row[0])
#             #A
#             bone1s = tem_row[3]
#             bone1e = tem_row[4]
#             #B
#             bone2s = tem_row[6]
#             bone2e = tem_row[7]
#             if whole_start_abs_set[ll_var_index] in list_start1:
#             #ls_satrt1_index_df = int(list_start1.index(whole_start_abs_set[ll_var_index]))
#             #stand_line = int(SV_table_merged.iloc[ls_satrt1_index_df,0])
#             #tem_row = SV_table_merged.iloc[ls_satrt1_index_df,:]
#                 #class of SV
#                 if tem_row[2] in ['Substitution','Small_Ins','Small_Del','Deletion','Insertion','Inversion']:
#                     if tem_row[2] in ['Deletion','Small_Del']:
#                         inster_number_bone = bone1s-last_bone-1
#                         #index for consensus before start of current variation
#                         present_len = present_len + inster_number_bone
#                         #update last_bone as end of current variation
#                         last_bone = bone1e
#                         #deleted base has no new axis on consensus
#                         SV_table_merged.iloc[stand_line,10] = -1
#                         SV_table_merged.iloc[stand_line,11] = -1
#                         SV_table_merged.iloc[stand_line,12] = -1
#                         SV_table_merged.iloc[stand_line,13] = -1
#                     elif tem_row[2] in ['Substitution']:
#                         inster_number_bone = bone1s-last_bone
#                         #one to one map
#                         present_len = present_len + inster_number_bone
#                         #bone1s=bone1e=5
#                         last_bone = bone1e
#                         SV_table_merged.iloc[stand_line,10] = present_len
#                         SV_table_merged.iloc[stand_line,11] = present_len
#                         SV_table_merged.iloc[stand_line,12] = -1
#                         SV_table_merged.iloc[stand_line,13] = -1
#                     elif tem_row[2] in ['Small_Ins','Insertion']:
#                         inster_number_bone = bone1s-last_bone
#                         Ins_len_present = len(tem_ins_dic[bone1s])
#                         #inserted position on consensus: one pos:+1, inserted after current base
#                         SV_table_merged.iloc[stand_line,10] = present_len + inster_number_bone+1
#                         SV_table_merged.iloc[stand_line,12] = -1
#                         #on consensus: end of previous SV+ number of normal base+ inserted length
#                         present_len = present_len + inster_number_bone+Ins_len_present
#                         #end of current SV
#                         last_bone = bone1e
#                         SV_table_merged.iloc[stand_line,11] = present_len
#                         SV_table_merged.iloc[stand_line,13] = -1
#                     else:## this is the inversion
#                         inster_number_bone = bone1s-last_bone
#                         SV_table_merged.iloc[stand_line,10] = present_len + inster_number_bone
#                         SV_table_merged.iloc[stand_line,12] = -1
#                         #no loss from last_bone to bone1e
#                         #????
#                         #present_len = present_len + bone1e - last_bone
#                         present_len = present_len + bone1e - last_bone
#                         SV_table_merged.iloc[stand_line,11] = present_len
#                         SV_table_merged.iloc[stand_line,13] = -1 
#                         last_bone = bone1e
                        
#                 elif tem_row[2] in ['Duplication']:
#                         #copy A to B (A no change)
#                         #5-0=5
#                         inster_number_bone = bone1s-last_bone
#                         #Ins_len_present = len(tem_ins_dic[bone2s])
#                         #length of the copied: A
#                         #=6
#                         tem_plate_len = SV_table_merged.iloc[stand_line,5]
#                         #0+5=5
#                         SV_table_merged.iloc[stand_line,10] = present_len + inster_number_bone
#                         #SV_table_merged.iloc[stand_line,12] = present_len + inster_number_bone+1
#                         present_len = present_len + inster_number_bone + tem_plate_len-1
#                         #0+5+6-1=10
#                         SV_table_merged.iloc[stand_line,11] = present_len 
#                         #SV_table_merged.iloc[stand_line,13] = present_len
#                         last_bone = bone1e
                        
#                 elif tem_row[2] in ['Translocation']:
#                     #balanced translocation
#                     #A:5-10, B:12-18
#                     if tem_row[9] == 1:
#                         #ins B to A's pos:5-0-1=4
#                         inster_number_bone = bone1s-last_bone-1
#                         #0+4+1=5,the start of copied base is 5
#                         SV_table_merged.iloc[stand_line,10] = present_len + inster_number_bone+1
#                         #length of B: 18-12+1=7
#                         Ins_len_present = len(tem_ins_dic[bone1s-1])
#                         #0+4+7=11
#                         #end of A:current SV end=11
#                         present_len = present_len + inster_number_bone + Ins_len_present
#                         SV_table_merged.iloc[stand_line,11] = present_len
#                         last_bone = bone1e
#                     #!unbalanced trans:
#                     else:
#                         inster_number_bone = bone1s-last_bone-1
#                         #index for consensus before start of current variation
#                         present_len = present_len + inster_number_bone
                        
#                         #deleted base has no new axis on consensus
#                         SV_table_merged.iloc[stand_line,10] = present_len+1
#                         SV_table_merged.iloc[stand_line,11] = present_len+1
                        
#                         #update last_bone as end of current variation
#                         last_bone = bone1e
        
        
#             else:### in the list2: pos of B (only duplication and trans)
                
#                 if tem_row[2] in ['Duplication']:
#                     #if SV_table_merged.iloc[stand_line,10]==0:
#                         #bone2s:B_start
#                         #same as ins
#                         inster_number_bone = bone2s-last_bone
#                         #SV_table_merged.iloc[stand_line,10] = present_len + inster_number_bone+1
#                         #SV_table_merged.iloc[stand_line,12] = present_len + inster_number_bone+1
#                         SV_table_merged.iloc[stand_line,12] = present_len+inster_number_bone+1
#                         Ins_len_present = len(tem_ins_dic[bone2s])
#                         present_len = present_len + inster_number_bone+Ins_len_present
                        
#                         SV_table_merged.iloc[stand_line,13] = present_len
#                         last_bone = bone2e
#                 elif tem_row[2] in ['Translocation']:
#                     #balanced: similar to A
#                     if  tem_row[9] == 1:
#                         inster_number_bone = bone2s-last_bone-1
#                         SV_table_merged.iloc[stand_line,12] = present_len + inster_number_bone+1
#                         #inserted A's length
#                         Ins_len_present = len(tem_ins_dic[bone2s-1])
#                         present_len = present_len + inster_number_bone + Ins_len_present
#                         SV_table_merged.iloc[stand_line,13] = present_len
#                         last_bone = bone2e
#                     #unbalanced
#                     else:
#                         inster_number_bone = bone2s-last_bone-1
#                         inster_number_bone = bone2s-last_bone
#                         #A is a del
#                         SV_table_merged.iloc[stand_line,10] = -1
#                         SV_table_merged.iloc[stand_line,12] = present_len + inster_number_bone+1
#                         #length of A
#                         #Ins_len_present = len(tem_ins_dic[bone2s-1])
#                         #Ins_dic_sv_seg[ins_trans_loc] = copy.deepcopy(''.join(tem_seq_post[r_s:(r_s+l_s)]))
#                         #similar to insertion
#                         Ins_len_present = len(tem_ins_dic[bone2s])
#                         present_len = present_len + inster_number_bone + Ins_len_present
#                         #A is a del
#                         SV_table_merged.iloc[stand_line,11] = -1
#                         SV_table_merged.iloc[stand_line,13] = present_len
#                         last_bone = bone2e
#     return SV_table_merged
    
def gen_consensus(xun, starts_seg, ends_seg, tem_seq_post, del_SV_per_segment, ins_SV_per_segment, del_snv_per_segment,ins_snv_per_segment,snp_per_segment, pai_pro_tem_, args, condition_dist_sv_del, len_SV_del,condition_dist_sv_ins,len_SV_ins,ins_selection,base_list,substitution_matrix,mis_selection,times,tmp_dir, TR_set, p_del_TR, p_ins_TR):   
# 设置进程特定的种子
    process_seed = args.seed * 10000 + xun
    random.seed(process_seed)
    np.random.seed(process_seed)

    loop_index_seg = int(0)
    
    #copy the tem_seq_post with already the deleted part from Trans, Inv and DUP
    #!
    start_seg = starts_seg[xun]
    end_seg = ends_seg[xun]
    tem_seq_post_seg = copy.deepcopy(tem_seq_post[start_seg:end_seg])
    print('length of no.'+str(xun) + ' seq: '+str(len(tem_seq_post_seg)))
    #dictionary, restore ins sites and corresponding content
    Ins_dic_sv_seg = {}

    # 从临时文件加载 unblock_region_vec[xun]
    
    unblock_region_vec_xun = np.load(os.path.join(tmp_dir, 'unblock_region_vec_{}.npy'.format(xun)))

    unblock_region_seg = list(np.copy(unblock_region_vec_xun))
    #! long del
    #print(str(xun)+'long del')
    ## long del part, we need get rid of the condition that overlapped the region of inserted trans
    #p_rand_del_sv = np.random.rand(1)

    True_del_number = del_SV_per_segment[xun]
    
    True_ins_number = ins_SV_per_segment[xun]
    

    len_unblock_region = len(unblock_region_seg)

    l_s_vec_ini = np.random.multinomial(n=len_unblock_region, pvals=pai_pro_tem_)

    if args.snv_del is None:
        del_snv_number = int(l_s_vec_ini[0])
    else:
        True_del_snv_number = del_snv_per_segment[xun]
        del_snv_number = min(True_del_snv_number,len_unblock_region)
        
    if args.snv_ins is None:
        ins_snv_number = (np.random.binomial(len(unblock_region_seg),diff_ins_prob_ins_real,1))[0]
    else:
        True_ins_snv_number = ins_snv_per_segment[xun]
        ins_snv_number = min(True_ins_snv_number, len_unblock_region)
        
    if args.snp is None:
        snp=int(l_s_vec_ini[1])
    else:
        True_snp_number = snp_per_segment[xun]
        snp = min(True_snp_number, len_unblock_region)
    
    max_length_numpy = True_del_number + True_ins_number + del_snv_number + ins_snv_number + snp
    
    # 初始化 numpy 数组
    # SV_table_seg = np.empty((int(max_length_numpy*2), 14), dtype=object)
    # VCF_table_seg = np.empty((int(max_length_numpy*2), 8), dtype=object)
    Unified_table_seg = np.empty((max_length_numpy, 10), dtype=object)

    
    
    circular_count_del_break = 0 
    
    
    #!
    if not unblock_region_seg:# 如果segment是空集，输出警告
        print("Warning: no available positions for variations in no."+str(xun)+" segment")
    else:
        #! 如果不是空集才会计算
        # 把unblock_region_seg转换成集合
        unblock_region_set = set(unblock_region_seg)

        # 计算unblock_region_set和TR_set的交集，保存到unblock_region_TR
        unblock_region_TR = list(unblock_region_set & TR_set)

        # 计算unblock_region_set和unblock_region_TR的差集，保存到unblock_region_nonTR
        unblock_region_nonTR = list(unblock_region_set - set(unblock_region_TR))

    # 对每一段进行处理
    if not unblock_region_seg:# 如果segment是空集，输出警告
        print("Warning: no available positions for long deletions in no."+str(xun)+" segment")
        circular_count_del_break = 1
    elif True_del_number == 0:
        circular_count_del_break = 1
    else:
        circular_count_del_break = 0
        len_seg_refine2 = max(unblock_region_seg)
        #print('SV DEL number:'+str(True_del_number))

    if not circular_count_del_break:
        print('SV DEL number in no.'+str(xun)+' segment: '+str(True_del_number))
        for del_index in range(True_del_number):
            # 根据概率在unblock_region_TR和unblock_region_nonTR中选择起点
            if np.random.rand() < p_del_TR:
                r_s = np.random.choice(unblock_region_TR)
            else:
                r_s = np.random.choice(unblock_region_nonTR)

            l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_del)
            l_s_index = list(l_s_vec).index(1)
            l_s = len_SV_del[int(l_s_index)]
            circular_count_del = 0
            circular_count_del_break = 0
            while (r_s+l_s>=len_seg_refine2) or (r_s+l_s-1 not in unblock_region_seg) or (unblock_region_seg.index(r_s+l_s-1)-unblock_region_seg.index(r_s)<l_s-1):
                if np.random.rand() < p_del_TR:
                    r_s = np.random.choice(unblock_region_TR)
                else:
                    r_s = np.random.choice(unblock_region_nonTR)

                l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_del)
                l_s_index = list(l_s_vec).index(1)
                l_s = len_SV_del[int(l_s_index)]
                circular_count_del = circular_count_del + 1
                if circular_count_del > args.times:
                    circular_count_del_break = 1
                    print("Warning: No."+str(del_index)+ "  long deletion sampling times exceeds "+str(times)+' times. Try: reduce number of variations, increase times or reference length.')
                    break
                else:
                    circular_count_del_break = 0
            if not circular_count_del_break:
                circular_count_del = 0
                #block the part deleted in the first long del
                unblock_region_seg = list(set(unblock_region_seg) - set(range(r_s-1,r_s+l_s+1)))
                unblock_region_TR = list(set(unblock_region_TR) - set(range(r_s-1,r_s+l_s+1)))
                unblock_region_nonTR = list(set(unblock_region_nonTR) - set(range(r_s-1,r_s+l_s+1)))

                # SV_table_seg[SV_loop_seg] = [SV_loop_seg,ll_c,'Deletion',r_s,r_s+l_s-1,l_s,-1,-1,-1,-1,0,0,0,0]
                # SV_loop_seg = SV_loop_seg + 1 

                # VCF_table_seg[VCF_loop_seg] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop_seg), ''.join(tem_seq_post_seg[r_s-start_seg:r_s-start_seg+l_s]) + tem_seq_post_seg[r_s-start_seg+l_s], tem_seq_post_seg[r_s-start_seg+l_s], '.', 'PASS', 'SVTYPE=DEL;SVLEN='+str(l_s)]

                # VCF_loop_seg = VCF_loop_seg + 1  

                Unified_table_seg[loop_index_seg] = ['Deletion',r_s,r_s+l_s-1,l_s,-1,-1,-1,-1, ''.join(tem_seq_post_seg[r_s-start_seg:r_s-start_seg+l_s]) + tem_seq_post_seg[r_s-start_seg+l_s], tem_seq_post_seg[r_s-start_seg+l_s]]
                loop_index_seg = loop_index_seg + 1
                tem_seq_post_seg[r_s-start_seg:(r_s-start_seg+l_s)] = '-'*l_s          
            else:
                continue
        
    
    #! long insertion
    whole_insertion_term = []

    circular_count_ins_break = 0
    ins_pos_collection_seg = []
    True_ins_number = ins_SV_per_segment[xun]

    if not unblock_region_seg:
        print("Warning: no available positions for long insertions in no."+str(xun)+" segment")
        circular_count_ins_break = 1
    elif True_ins_number == 0:
        circular_count_ins_break = 1
    else:
        circular_count_ins_break = 0
        
    ins_SV_pos = []  # 默认初始化为一个空列表
    if not circular_count_ins_break:
        print('SV INS number in No.'+str(xun)+' segment: '+str(True_ins_number))
        if not unblock_region_seg:
            print("Warning: no available positions for long insertions in no."+str(xun)+" segment")
        else:
            len_seg_refine = max(unblock_region_seg)
            if len_seg_refine-1 in unblock_region_seg:
                unblock_region_seg.remove(len_seg_refine-1)
        if len(unblock_region_seg) < True_ins_number:
            ins_SV_pos = unblock_region_seg
            print('Warning: no enough positions for '+str(args.sv_ins)+' long insertions. Only '+str(len(ins_SV_pos))+' long insertions will be generated.')
        else:
            num_ins_TR = int(True_ins_number * p_ins_TR)
            num_ins_nonTR = True_ins_number - num_ins_TR
            if len(unblock_region_TR) < num_ins_TR:
                num_ins_TR = len(unblock_region_TR)
                num_ins_nonTR = True_ins_number - num_ins_TR
            if len(unblock_region_nonTR) < num_ins_nonTR:
                num_ins_nonTR = len(unblock_region_nonTR)
            # 检查unblock_region_TR和unblock_region_nonTR是否为空
            if unblock_region_TR and unblock_region_nonTR:
                ins_SV_pos_TR = random.sample(unblock_region_TR, num_ins_TR)
                ins_SV_pos_nonTR = random.sample(unblock_region_nonTR, num_ins_nonTR)
                ins_SV_pos = ins_SV_pos_TR + ins_SV_pos_nonTR
        ins_pos_collection_seg.extend(ins_SV_pos)

        for remain_index2 in ins_pos_collection_seg:
            l_i_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_ins)
            l_i_ins = int(list(l_i_vec).index(1))
            l_i = len_SV_ins[l_i_ins]
            tem_ins = ''
            for jins in range(l_i):
                bexixuan_ = choice(ins_selection)
                tem_ins = tem_ins + bexixuan_

            Ins_dic_sv_seg[remain_index2] = tem_ins

            # 检查元素是否存在于列表中
            if remain_index2 in unblock_region_seg:
                unblock_region_seg.remove(remain_index2)
            if remain_index2 in unblock_region_TR:
                unblock_region_TR.remove(remain_index2)
            if remain_index2 in unblock_region_nonTR:
                unblock_region_nonTR.remove(remain_index2)
            if remain_index2+1 in unblock_region_seg:
                unblock_region_seg.remove(remain_index2+1)
            if remain_index2+1 in unblock_region_TR:
                unblock_region_TR.remove(remain_index2+1)
            if remain_index2+1 in unblock_region_nonTR:
                unblock_region_nonTR.remove(remain_index2+1)
            whole_insertion_term.append(tem_ins)

            # SV_table_seg[SV_loop_seg] = [SV_loop_seg,ll_c,'Insertion',remain_index2,remain_index2,l_i,-1,-1,-1,-1,0,0,0,0]
            # SV_loop_seg = SV_loop_seg + 1

            # VCF_table_seg[VCF_loop_seg] = [str(chr_id), str(remain_index2), 'rs' + str(VCF_loop_seg), tem_seq_post_seg[remain_index2-start_seg], tem_seq_post_seg[remain_index2-start_seg] + tem_ins.upper(), '.', 'PASS', 'SVTYPE=INS;SVLEN='+str(l_i)]
            # VCF_loop_seg = VCF_loop_seg + 1
            Unified_table_seg[loop_index_seg] = ['Insertion',remain_index2,remain_index2,l_i,-1,-1,-1,-1, tem_seq_post_seg[remain_index2-start_seg], tem_seq_post_seg[remain_index2-start_seg] + tem_ins.upper()]
            loop_index_seg = loop_index_seg + 1
    #!

    #! small del

    ### finish the process of the duplication and insetion
    ## possible var bone
    #snp_considered_sites = copy.deepcopy(unblock_region_seg)
    # unblock_region_snv= copy.deepcopy(unblock_region_seg)#deletions
    # undel_region_snv=copy.deepcopy(undel_region_seg)#substitutions
    # left_del_region_snv=copy.deepcopy(left_del_region_seg)#insertions
    #!!
    # 对每一段进行处理
    if not unblock_region_seg:# 如果segment是空集，输出警告
        print("Warning: no available positions for small deletions in no."+str(xun)+" segment")
        circular_count_micro_del_break = 1
    else:
        circular_count_micro_del_break = 0
        
    # 计算需要选择的位点的数量
    if not circular_count_micro_del_break:
        ## small deletion part
        len_unblock_region = len(unblock_region_seg)
        l_s_vec_ini = np.random.multinomial(n=len_unblock_region, pvals=pai_pro_tem_)
        


        # Define the maximum small deletion length
        max_small_deletion_length = 5  # Maximum small deletion length

        # Calculate the space required for each small deletion
        min_gap = 1  # Minimum gap between each small deletion

        # Calculate the maximum number of small deletions
        max_micro_dels = (len_unblock_region + min_gap) // (max_small_deletion_length + min_gap)
        if args.snv_del is None:
            del_snv_number = int(l_s_vec_ini[0])
            #del_snv_number = test_number
        else:
            if True_del_snv_number > max_micro_dels:
                del_snv_number = min(True_del_snv_number, max_micro_dels)
                print("Warning: The input for small deletions is too large and has been automatically reduced. \
                    Each small deletion event requires space, including gaps between them.")
          
        len_seg_refine = max(unblock_region_seg)
        # for each deletion
        for m_del in range(0,del_snv_number):
            r_s = sample(unblock_region_seg,1)[0]
            l_s_vec = np.random.multinomial(n=1, pvals=[6/8,1/8,1/16,1/32,1/32])
            #the length of this deletion
            l_s = list(l_s_vec).index(1)+1
            #resample if the end of this deletion is out of range

            # count the number of times that we resample
            circular_count_del = 0
            # indicator: if resampling exceeds 50 times
            circular_count_del_break = 0
            while (r_s+l_s>len_seg_refine) or (r_s+l_s-1 not in unblock_region_seg) or (unblock_region_seg.index(r_s+l_s-1)-unblock_region_seg.index(r_s)<l_s-1):
                # select the possibile deleltion point
                r_s = sample(unblock_region_seg,1)[0]
                l_s_vec = np.random.multinomial(n=1, pvals=[6/8,1/8,1/16,1/32,1/32])
                #the length of this deletion
                l_s = list(l_s_vec).index(1)+1
                #l_s = choice(del_length_list)
                # count the times of resampling
                circular_count_del = circular_count_del + 1
                if circular_count_del> args.times:
                    circular_count_del_break = 1
                    print("Warning: No."+str(m_del)+ "  small deletion sampling times exceeds "+str(times)+' times. Try: reduce number of variations, increase times or reference length.')
                    break
                else:
                    circular_count_del_break = 0
            # if sample a del start and length that is not overlapped with translocation(s)
            # update needed pos collections and prepare for the second sampling of long deletions
            if not circular_count_del_break:
                #initialize the counting of resampling times
                circular_count_del = 0
                
                #block the part deleted in the first long del
                unblock_region_seg = list(set(unblock_region_seg) -set(range(r_s-1,r_s+l_s+1)))
            
                # SV_table_seg[SV_loop_seg] = [SV_loop_seg,ll_c,'Small_Del',r_s,r_s+l_s-1,l_s,-1,-1,-1,-1,0,0,0,0]
                # SV_loop_seg = SV_loop_seg + 1

                #The ‘REF’ column is set to the original segment plus the base that is left after the deletion. 
                #The ‘ALT’ column is set to the base that is left after the deletion.
                # Add a row to the VCF_table for the small deletion
                #VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), tem_seq_post[r_s:r_s+l_s] + tem_seq_post[r_s+l_s], tem_seq_post[r_s+l_s], '.', 'PASS', 'SVTYPE=microDEL']
                # VCF_table_seg[VCF_loop_seg] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop_seg), ''.join(tem_seq_post_seg[r_s-start_seg:r_s-start_seg+l_s]) + tem_seq_post_seg[r_s-start_seg+l_s], tem_seq_post_seg[r_s-start_seg+l_s], '.', 'PASS', 'microDEL;LEN='+str(l_s)]

                # VCF_loop_seg = VCF_loop_seg + 1  
                
                Unified_table_seg[loop_index_seg] = ['Small_Del',r_s,r_s+l_s-1,l_s,-1,-1,-1,-1, ''.join(tem_seq_post_seg[r_s-start_seg:r_s-start_seg+l_s]) + tem_seq_post_seg[r_s-start_seg+l_s], tem_seq_post_seg[r_s-start_seg+l_s]]
                loop_index_seg = loop_index_seg + 1
                #replace deleted bases with -
                tem_seq_post_seg[r_s-start_seg:(r_s-start_seg+l_s)] = '-'*l_s
            else:
                break
                
        
    ### small insertions
    # 对每一段进行处理
    #if not left_del_region_snv:# 如果segment是空集，输出警告
    if not unblock_region_seg:# 如果segment是空集，输出警告
        print("Warning: no available positions for small insertions in no."+str(xun)+" segment")
        circular_count_micro_ins_break = 1
    else:
        circular_count_micro_ins_break = 0
        
    # 计算需要选择的位点的数量
    if not circular_count_micro_ins_break:
        len_unblock_region = len(unblock_region_seg)
        if args.snv_ins is None:
            #ins_snv_number = (np.random.binomial(len(unblock_region_seg),diff_ins_prob_ins_real,1))[0]
            diff_ins_prob_ins_real = 0.00038
            ins_snv_number = (np.random.binomial(len_unblock_region,diff_ins_prob_ins_real,1))[0]
            #ins_snv_number = test_number
        else:
            ins_snv_number = min(True_ins_snv_number, len_unblock_region)
        #all possitions for ins, a loop
        for number_ins in range(0,ins_snv_number):
            #the positions of insertions
            #remain_index2 = sample(left_del_region_snv,1)[0]
            remain_index2 = sample(unblock_region_seg,1)[0]
            #length of insertion
            l_i_vec = np.random.multinomial(n=1, pvals=[6/8,1/8,1/16,1/32,1/32])
            l_i = list(l_i_vec).index(1)+1
            #l_i=choice(ins_length_list)
            tem_ins=''
            #a loop to choose inserted base
            for j in range(l_i):
                bexixuan_ = choice(ins_selection)
                tem_ins = tem_ins + bexixuan_
                
            #record ins position and ins segment
            Ins_dic_sv_seg[remain_index2] = tem_ins
            
            #update index set
            # if remain_index2 in left_del_region_snv:
            #     left_del_region_snv.remove(remain_index2)
            if remain_index2 in unblock_region_seg:
                unblock_region_seg.remove(remain_index2)
            if remain_index2+1 in unblock_region_seg:
                unblock_region_seg.remove(remain_index2+1)
            #! substitution cannot be on the inserted pos
            # if remain_index2 in undel_region_snv:
            #     undel_region_snv.remove(remain_index2)
            whole_insertion_term.append(tem_ins)

            #num_collection.append(l_i)
            #only one record in the table for each remain_index2
            # SV_table_seg[SV_loop_seg] = [SV_loop_seg,ll_c,'Small_Ins',remain_index2,remain_index2,l_i,-1,-1,-1,-1,0,0,0,0]
            # SV_loop_seg = SV_loop_seg + 1

            # # Add a row to the VCF_table for the small insertion
            # VCF_table_seg[VCF_loop_seg] = [str(chr_id), str(remain_index2), 'rs' + str(VCF_loop_seg), tem_seq_post_seg[remain_index2-start_seg], tem_seq_post_seg[remain_index2-start_seg] + tem_ins.upper(), '.', 'PASS', 'microINS;LEN='+str(l_i)]
            # VCF_loop_seg = VCF_loop_seg + 1  
            Unified_table_seg[loop_index_seg] = ['Small_Ins',remain_index2,remain_index2,l_i,-1,-1,-1,-1, tem_seq_post_seg[remain_index2-start_seg], tem_seq_post_seg[remain_index2-start_seg] + tem_ins.upper()]
            loop_index_seg = loop_index_seg + 1
        #snp=int(l_s_vec_ini[1])
    
    #! start SNP
    ### small insertions
    # 对每一段进行处理
    #if not undel_region_snv:# 如果segment是空集，输出警告
    if not unblock_region_seg:# 如果segment是空集，输出警告
        print("Warning: no available positions for substitutions in no."+str(xun)+" segment")
        circular_count_snv_break = 1
    else:
        circular_count_snv_break = 0
        
    # 计算需要选择的位点的数量
    if not circular_count_snv_break:
        len_unblock_region = len(unblock_region_seg)
        if args.snp is None:
            l_s_vec_ini = np.random.multinomial(n=len_unblock_region, pvals=pai_pro_tem_)
            snp=int(l_s_vec_ini[1])
            #snp=test_number
        else:
            snp = min(True_snp_number, len_unblock_region)
        ### substitution
        for number_mis in range(0,snp):
            #pos of substitution
            # ll = sample(undel_region_snv,1)[0]
            # undel_region_snv.remove(ll)
            ll = sample(unblock_region_seg,1)[0]
            unblock_region_seg.remove(ll)
            #the base that is replaced
            # if tem_seq_post_seg[ll-start_seg] == 'N':
            #     #bexixuan_ = 'N'
            #     print('Error: SNP in Gap region')
            # else:
            if tem_seq_post_seg[ll-start_seg].upper() in base_list:
                ref_id = base_list.index(tem_seq_post_seg[ll-start_seg].upper())
            
                #ref_id=base_list.index(tem_seq_post_seg[ll-start_seg].upper())
                #selection the ref_id's probability distribution
                prob_dist = substitution_matrix[ref_id]
                #sample a column_index 
                column_index = np.random.choice(4, p=prob_dist)
                #choose a position for mismatch
                bexixuan_ = mis_selection[column_index]
                #! update index set for adding INDEL
                # if ll in left_del_region_snv:
                #     left_del_region_snv.remove(ll)
                if ll in unblock_region_seg:
                    unblock_region_seg.remove(ll)

                # SV_table_seg[SV_loop_seg] = [SV_loop_seg,ll_c,'Substitution',ll,ll,1,-1,-1,-1,-1,0,0,0,0]
                # SV_loop_seg = SV_loop_seg + 1
                # # Add a row to the VCF_table for the substitution
                # VCF_table_seg[VCF_loop_seg] = [str(chr_id), str(ll), 'rs' + str(VCF_loop_seg), tem_seq_post_seg[ll-start_seg], bexixuan_.upper(), '.', 'PASS', 'SUB']
                # VCF_loop_seg = VCF_loop_seg + 1
                
                Unified_table_seg[loop_index_seg] = ['Substitution',ll,ll,1,-1,-1,-1,-1, tem_seq_post_seg[ll-start_seg], bexixuan_.upper()]
                loop_index_seg = loop_index_seg + 1

                tem_seq_post_seg[int(ll-start_seg)] = copy.deepcopy(bexixuan_)
            else:
                print("Error: Invalid base for substitution")
        
    # # 找出所有元素都是 None 或空字符串的行
    # mask1 = np.all((SV_table_seg == None) | (SV_table_seg == ''), axis=1)

    # # 删除这些行
    # SV_table_seg = SV_table_seg[~mask1]

    # # 现在你可以将 SV_table 转换为 DataFrame
    # SV_table_seg = pd.DataFrame(SV_table_seg, columns=['Index','Index_con','SV_type','Original_start',\
    #                                     'Original_end','Len_SV','New_start','New_end','New_len_SV','Balanced Trans Flag','relative start1','relative end1','relative start2','relative end2'])

    # # 对 VCF_table 重复相同的步骤
    # mask2 = np.all((VCF_table_seg == None) | (VCF_table_seg == ''), axis=1)
    # VCF_table_seg = VCF_table_seg[~mask2]
    # VCF_table_seg = pd.DataFrame(VCF_table_seg, columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
    
        
    for idx in sorted(Ins_dic_sv_seg, reverse=True):
        tem_seq_post_seg.insert(idx-start_seg, Ins_dic_sv_seg[idx])

    # 找出所有元素都是 None 或空字符串的行
    mask1 = np.all((Unified_table_seg == None) | (Unified_table_seg == ''), axis=1)
    # 删除这些行
    Unified_table_seg = Unified_table_seg[~mask1]
    # 输出处理后的表格形状
    print('Table shape:', Unified_table_seg.shape)
    print('No.'+str(xun)+' seg is saving files')
    # return SV_table_seg, VCF_table_seg, unblock_region_seg, Ins_dic_sv_seg, tem_seq_post_seg
    return Unified_table_seg, Ins_dic_sv_seg, tem_seq_post_seg

# 定义一个全局变量来存储最大的磁盘使用率
def monitor_disk_usage(max_disk_usage):
    while True:
        disk_usage = psutil.disk_usage('/').percent
        with max_disk_usage.get_lock():
            if disk_usage > max_disk_usage.value:
                max_disk_usage.value = disk_usage
        time.sleep(1)

# Monitor memory usage
def monitor_memory(threshold, max_mem_usage, max_mem_usage_gb, tmp_dir):
    """
    Monitor the memory usage, stop the program and release resources if memory usage exceeds the threshold.
    """
    mem_info = psutil.virtual_memory()  # get virtual memory information
    mem_usage = mem_info.percent  # get virtual memory usage in percentage
    max_mem_usage = max(mem_usage, max_mem_usage)  # update max memory usage
    max_mem_usage_gb = max(max_mem_usage / 100 * mem_info.total / (1024 ** 3), max_mem_usage_gb)  # calculate max memory usage in GB
    if mem_usage > threshold:
        print(f"Memory usage exceeded! Current memory usage: {mem_usage}%")
        # Stop multiprocessing jobs
        for p in multiprocessing.active_children():
            p.terminate()
        # Remove temporary directory
        #shutil.rmtree('tmp')
        # Remove temporary directory
        shutil.rmtree(tmp_dir)

        raise SystemExit("Program stopped due to excessive memory usage.")
    return max_mem_usage, max_mem_usage_gb



def parse_args():
    parser = argparse.ArgumentParser(description='BVSim')
    parser.add_argument('-ref', type=str, help='Input reference local path', default='default_ref')
    parser.add_argument('-seq_index', type=int, default=0, 
                      help='Index of sequence to use (0-based). Default: 0 (first sequence)')
    parser.add_argument('-save', type=str, help='local path for saving', default=os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'save')+ '/')
    parser.add_argument('-seed', type=int, help='Seed for random number generator', default=999)
    parser.add_argument('-times', type=int, help='Number of times', default=10)
    parser.add_argument('-rep', type=int, help='Replication ID', default=99)
    parser.add_argument('-sv_trans', type=int, help='Number of trans SV', default=5)
    parser.add_argument('-sv_inver', type=int, help='Number of inversion SV', default=5)
    parser.add_argument('-sv_dup', type=int, help='True duplication number', default=5)
    parser.add_argument('-sv_del', type=int, help='Number of deletion SV', default=5)
    parser.add_argument('-sv_ins', type=int, help='True insertion number', default=5)
    parser.add_argument('-snp', type=float, help='SNV number or probability', default=5)
    parser.add_argument('-snv_del', type=float, help='SNV deletion number or probability', default=5)
    parser.add_argument('-snv_ins', type=float, help='SNV insertion number or probability', default=5)
    parser.add_argument('-notblockN', action='store_true', help='Do not Block N positions')
    parser.add_argument('-write', action='store_true', help='Write full results')
    parser.add_argument('-block_region_bed_url', type=str, help='local path of the block region BED file', default=None)
    parser.add_argument('-cores', type=int, help='Number of kernels for parallel processing', default=1)
   
    parser.add_argument('-len_bins', type=int, help='Length of bins for parallel processing, must be >0 and <reference length', default=50000)
    parser.add_argument('-delmin', type=int, help='Minimum deletion length', default=50)
    parser.add_argument('-delmax', type=int, help='Maximum deletion length', default=60)
    parser.add_argument('-insmin', type=int, help='Minimum insertion length', default=50)
    parser.add_argument('-insmax', type=int, help='Maximum insertion length', default=450)
    parser.add_argument('-dupmin', type=int, help='Minimum duplication length', default=50)
    parser.add_argument('-dupmax', type=int, help='Maximum duplication length', default=450)
    parser.add_argument('-invmin', type=int, help='Minimum inversion length', default=50)
    parser.add_argument('-invmax', type=int, help='Maximum inversion length', default=450)
    parser.add_argument('-transmin', type=int, help='Minimum translocation length', default=50)
    parser.add_argument('-transmax', type=int, help='Maximum translocation length', default=450)
    parser.add_argument('-mode', type=str, help='Mode for calculating probabilities', default='probability')
    parser.add_argument('-sum', action='store_true', help='total indel SV equals sum of the input bed')
    #parser.add_argument('-indel_input_bed', '--indel_input_bed', type=str, help='Input BED file for indels', default='~/data/test_data/TGS/hg002/chr21_SV_Tier1.bed')
    parser.add_argument('-p_del_region', type=float, help='Probability of SV DEL in the user-defined region for deletion', default=0.5)
    parser.add_argument('-p_ins_region', type=float, help='Probability of SV INS in the user-defined region for insertion', default=0.5)
    parser.add_argument('-region_bed_url', type=str, help='local path of the BED file for the user-defined region', default=None)
    parser.add_argument('-hg38', type=str, help='Chromosome name', required=False)
    return parser.parse_args()

# test_number = 2
def main():
    start_time1 = time.time()
    args = parse_args()

    # 设置main_url_empirical为固定的路径
    main_url_empirical = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'empirical') + '/'

    max_disk_usage = multiprocessing.Value('d', 0.0)
    max_mem_usage = 0  # initialize max memory usage
    max_mem_usage_gb = 0  # initialize max memory usage in GB

    #! start monitor process
    monitor_process = multiprocessing.Process(target=monitor_disk_usage, args=(max_disk_usage,))
    monitor_process.start()

    print('start')
    # Create a temporary directory if not exists under args.save
    timestamp = time.strftime("%Y%m%d_%H%M%S")  # 格式化时间戳
    os.makedirs(args.save, exist_ok=True)
    tmp_dir = os.path.join(args.save, f'tmp_{timestamp}')  # 使用时间戳作为文件夹名称

    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # random.seed(args.seed)
    # np.random.seed(args.seed) 
    # 设置全局随机数种子
    random.seed(args.seed)
    np.random.seed(args.seed)
    os.environ['PYTHONHASHSEED'] = str(args.seed)

    # 打印种子信息
    print(f"Global Random Seed: {args.seed}")

    # fasta文件路径
    # fasta_file_path = args.ref
    # fasta_file = pysam.FastaFile(fasta_file_path)

    # # 获取最长的序列
    # seqname = fasta_file.references[0]
    # ref_seq = fasta_file.fetch(seqname)
    # 1. 加载参考序列
    seqname, ref_seq = load_reference_sequence(args.ref, args.seq_index)
    print(f"Loaded reference sequence: {seqname} (length: {len(ref_seq)})")
    #! whole chr
    chr_id=seqname
    # print('reference\'s name:'+str(seqname))
    # start_base=43531239
    # end_base=43601240
    #len:46709983
    start_base=0
    end_base=len(ref_seq)
    real_con1 = copy.deepcopy(ref_seq[start_base:end_base+1])
    chr_length = len(real_con1)
    del ref_seq
    # print('Length of ref:'+str(chr_length))

    all_positions = set(range(len(real_con1)))
    # print('N pos:'+str(len(n_positions)))
   
    # all_region_sv=list(all_positions - n_positions_set)
    if args.notblockN:
        all_region_sv= list(all_positions)
    else:
        # collection of N positions
        n_positions = [i for i, char in enumerate(real_con1) if char == 'N']
        n_positions_set = set(n_positions)
        n_bed_path = create_N_bed(n_positions, chr_id, tmp_dir)
        #all_positions = set(range(len(real_con1)))
        print('N pos:'+str(len(n_positions)))
        
        all_region_sv=list(all_positions - n_positions_set)

        del n_positions
    
    block_region_bed = args.block_region_bed_url
        
    if block_region_bed is None or block_region_bed == 'None':
        print("block_region_bed_url is None")
    else:
        print("block_region_bed_url: ")
        # 读取BED文件
        df_block = pd.read_csv(block_region_bed, sep='\t', header=None, names=['start', 'end'])

        # 使用.apply()函数应用上面定义的函数，得到每一行'start'和'end'值的检查结果
        df_block['check'] = df_block.apply(check_start_end, axis=1)

        # 检查是否所有行的'start'值都小于'end'值
        if df_block['check'].all():
            print("All rows are valid.")
        else:
            print("Error: Not all rows are valid.")
            sys.exit(1)


        # 使用.apply()函数应用上面定义的函数，得到每一行'start'和'end'区域中的所有点
        df_block['points'] = df_block.apply(get_points, axis=1)

        # 使用set.union()函数合并所有的点，得到一个包含所有点的集合
        block_set = set.union(*df_block['points'])
        all_region_sv = list(set(all_region_sv)-block_set)
        # 现在，'block_set'集合中包含了所有的点
        print('length of the blocked region:'+str(len(block_set))) 
        
    unblock_region_sv = copy.deepcopy(all_region_sv)
    
    tem_seq_post = copy.deepcopy(list(real_con1))
    del all_region_sv
    #! Cut the sequence for multiprocessing
    #! fixed bin size
    #! fixed bins length
    
    # 计算整除部分和余数部分
    devided_length = int(chr_length // args.len_bins)
    remainder_length = chr_length - devided_length * args.len_bins

    # 创建新的bins数组，最后一个bin的长度大于5e5
    bins = np.arange(0, chr_length - remainder_length, args.len_bins).tolist()
    bins.append(chr_length)
    
    # 将 all_region_sv 转换为 numpy 数组并排序
    all_region_sv_array = np.sort(np.array(copy.deepcopy(list(all_positions))))

    del all_positions
    
    if args.len_bins <= 0:
        raise ValueError(f"'len_bins' must be > 0, got {args.len_bins}")
        sys.exit(1)  # 非0状态码表示错误退出
    if args.len_bins > len(all_region_sv_array):
        print(
            f"ERROR: 'len_bins' ({args.len_bins}) > reference length ({len(all_region_sv_array)})",
            file=sys.stderr
        )
        sys.exit(1)  # 非0状态码表示错误退出
        
    # 计算能被 args.len_bins 整除的最大元素数量和余数
    number_seg, remainder = divmod(len(all_region_sv_array), int(args.len_bins))

    # 使用 numpy.array_split 将前 number_seg * args.len_bins 个元素划分为每 args.len_bins 个一组
    segments_initial = np.array_split(all_region_sv_array[:number_seg * int(args.len_bins)], number_seg)

    # 初始化存储起始和结束位置的列表
    starts_seg = np.empty(number_seg, dtype=int)
    ends_seg = np.empty(number_seg, dtype=int)

    # 将剩余的元素全部加到最后一个分组中
    if remainder > 0:
        #segments_initial.append(all_region_sv_array[number_seg * int(args.len_bins):])
        # 将剩余的元素全部加到最后一个分组中
        segments_initial[-1] = np.concatenate((segments_initial[-1], all_region_sv_array[number_seg * int(args.len_bins):]))

    # 将 numpy 数组转换回 Python 列表，并计算每个segment的起始和结束位置
    for i, segment in enumerate(segments_initial):
        starts_seg[i] = np.min(segment)
        ends_seg[i] = np.max(segment) + 1  # 加 1 是因为 Python 的切片操作是左闭右开的

    print('number of segments:'+str(number_seg))

    # 从.npy文件中读取数据
    len_SV_del = np.load(main_url_empirical+'len_SV_del_cell.npy').tolist()
    condition_dist_sv_del = np.load(main_url_empirical+'condition_dist_sv_del_cell.npy').tolist()

    #modifiable deletion length
    # delmin = args.delmin
    # delmax = args.delmax
    #modifiable deletion length
    delmin = max(50, int(args.delmin))
    delmax = int(args.delmax)
    # 检查和调整 delmax
    if delmin >= delmax:
        print(f"Warning: delmin ({delmin}) is greater than delmax ({delmax}). Adjusting delmax to {delmin + 50}.")
        delmax = delmin + 50  # 自动调整为比 delmin 大50的值
    # 选择长度在delmin和delmax之间的子集
    selected_lengths = [all_del_len for all_del_len in len_SV_del if delmin <= all_del_len <= delmax]

    # 找到这些长度在len_SV_del中的索引
    selected_indices = [len_SV_del.index(del_len) for del_len in selected_lengths]

    # 根据这些索引截取condition_dist_sv_del
    selected_probabilities = [condition_dist_sv_del[del_index] for del_index in selected_indices]

    # 归一化新的概率向量
    total = sum(selected_probabilities)
    normalized_probabilities = [prob_del_len/total for prob_del_len in selected_probabilities]

    # 更新len_SV_del和condition_dist_sv_del
    len_SV_del = selected_lengths
    condition_dist_sv_del = normalized_probabilities

    # 从.npy文件中读取数据
    len_SV_ins = np.load(main_url_empirical+'len_SV_ins_cell.npy').tolist()
    condition_dist_sv_ins = np.load(main_url_empirical+'condition_dist_sv_ins_cell.npy').tolist()

    # insmin = args.insmin
    # insmax = args.insmax
    #! 对概率进行归一化
    #! length probabilities
    insmin = max(50, int(args.insmin))
    insmax = int(args.insmax)

    # 检查和调整 insmax（同款逻辑）
    if insmin >= insmax:
        print(f"Warning: insmin ({insmin}) is greater than insmax ({insmax}). Adjusting insmax to {insmin + 50}.")
        insmax = insmin + 50  # 自动调整为比 insmin 大50的值


    # 选择长度在insmin和insmax之间的子集
    selected_lengths = [l for l in len_SV_ins if insmin <= l <= insmax]

    # 找到这些长度在len_SV_ins中的索引
    selected_indices = [len_SV_ins.index(l) for l in selected_lengths]

    # 根据这些索引截取condition_dist_sv_ins
    selected_probabilities = [condition_dist_sv_ins[i] for i in selected_indices]

    # 归一化新的概率向量
    total = sum(selected_probabilities)
    normalized_probabilities = [p/total for p in selected_probabilities]

    # 更新len_SV_ins和condition_dist_sv_ins
    len_SV_ins = selected_lengths
    condition_dist_sv_ins = normalized_probabilities


    #record the position and content of insertions, to update the final consensus
    Ins_dic_sv = {}
    #restore sites of ins for trans and dup in case we sample new ins positions there
    location_insert_trans = []
    location_insert_dup = []
    #restore sites of del for translocation, inversion
    location_del_trans = []
    location_del_inv = []
    ### long translocation part
    #p_rand_trans_sv = np.random.rand(1)
    #circular_count_trans = 0

    #! Translocation, Inversion and Duplication
    # #record original positions and translocated positions
    # SV_table = pd.DataFrame(columns=['Index','Index_con','SV_type','Original_start',\
    #                                     'Original_end','Len_SV','New_start','New_end','New_len_SV','Balanced Trans Flag','relative start1','relative end1','relative start2','relative end2'])

    # #record location, original base (REF) and alternative base (ALT) for substitution, INDEL
    # # Define a DataFrame to store the information
    # VCF_table = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
    
    max_csv_len = args.sv_trans + args.sv_dup
    max_sim_len = args.sv_inver
    # 初始化 numpy 数组
    SV_table = np.empty((int(max_csv_len+2), 8), dtype=object)  # 减少列数，去掉 'Index' 和 'Index_con'
    VCF_table = np.empty((int(max_csv_len*4), 4), dtype=object)  # 减少列数
    # 初始化统一的 numpy 数组
    Unified_table = np.empty((max_sim_len, 10), dtype=object)
    
    
    ll_c=0
    len_seg_refine = len(real_con1)
    
    #deletion probability
    # del_sv_createria = 8.9/9
    # trans_sv_createria = 8.9/9
    # inver_sv_createria = 8.9/9
    # ins_location_p = 1/90000
    # dup_location_p = 1/90000


    SV_loop = 0
    VCF_loop = 0
    CSV_loop = 0
    #! loop index
    loop_index = 0
    times = args.times
    updated_con = []
    ratio_b_trans = 0.5
    #ratio_re_dup = 0.5
    # Whole_INS_con = []#restore Ins_dic_sv for each consensus


    # num_con=1
    # mutation_rate_list = 1
    diff_ins_prob_correct_real = 0.9989
    diff_ins_prob_mis_real = 0.0139
    diff_ins_prob_del_real = 0.00038
    diff_ins_prob_ins_real = 0.00038


    substitution_matrix = np.array([
        [0.000000000000000000e+00, 1.449913454716791894e-01, 1.725324178807925157e-01, 6.824762366475283226e-01],
        [1.444905385549744570e-01, 0.000000000000000000e+00, 6.853654028089395389e-01, 1.701440586360860041e-01],
        [1.641872252310401514e-01, 6.652683754966074448e-01, 0.000000000000000000e+00, 1.705443992723524316e-01],
        [6.642860882078220897e-01, 1.643992946044358361e-01, 1.713146171877421298e-01, 0.000000000000000000e+00]
    ])

    ratio_b_trans = 0.5
    ratio_re_dup = 1

    base_list=["A","T","C","G"]
    ins_selection = ['a', 't', 'c', 'g']
    mis_selection = ['a','t','c','g']

    base_list=["A","T","C","G"]

    # pdel_SV = 0.01
    #len_SV_del = [50,100,500]
    # len_SV_trans = [50,100,500]
    
    #len_SV_ins = [50,100,500]
    #condition_dist_sv_del = [4/5,1/10,1/10]
    # condition_dist_sv_trans = [4/5,1/10,1/10]
    
    # len_SV_inver = [50,100,500]
    # condition_dist_sv_inver = [4/5,1/10,1/10]
    #condition_dist_sv_ins = [4/5,1/10,1/10]
    # number_sv_del = [4/5,1/10,1/10]
    # number_sv_trans = [4/5,1/10,1/10]
    # number_sv_inver = [4/5,1/10,1/10]
    # number_sv_dup = [4/5,1/10,1/10]
    # copied_base_sv_prob = [4/5,1/10,1/10]
    # copied_base_sv_base = [50,100,500]
    # 设置参数
    dupmin = max(50, int(args.dupmin))
    dupmax = int(args.dupmax)
    # 检查和调整 dupmax
    if dupmin >= dupmax:
        print(f"Warning: dupmin ({dupmin}) is greater than dupmax ({dupmax}). Adjusting dupmax to {dupmin + 50}.")
        dupmax = dupmin + 50  # 自动调整 dupmax

    # 计算整数范围
    dup_range = np.arange(dupmin, dupmax + 1)
    
    # 处理倒位长度范围
    invmin = max(50, int(args.invmin))  # 确保最小值不低于50bp且为整数
    invmax = int(args.invmax)           # 确保为整数

    # 检查和调整 invmax（完全同款逻辑）
    if invmin >= invmax:
        print(f"Warning: invmin ({invmin}) is greater than invmax ({invmax}). "
            f"Adjusting invmax to {invmin + 50}.")
        invmax = invmin + 50  # 自动调整为比 invmin 大50bp
        
    # 计算整数范围
    inv_range = np.arange(invmin, invmax + 1)
    # 计算每个整数的概率（均匀分布）
    uniform_prob = 1 / len(inv_range)
    # 创建新的概率向量
    inv_sv_prob = [uniform_prob] * len(inv_range)
    
    len_SV_inver = inv_range.tolist()
    condition_dist_sv_inver = inv_sv_prob
    
    transmin = max(50, int(args.transmin))
    transmax = int(args.transmax)

    # 检查和调整 transmax（同款逻辑）
    if transmin >= transmax:
        print(f"Warning: transmin ({transmin}) is greater than transmax ({transmax}). Adjusting transmax to {transmin + 50}.")
        transmax = transmin + 50  # 自动调整为比 transmin 大50的值

    # 计算整数范围
    trans_range = np.arange(transmin, transmax + 1)
    # 计算每个整数的概率（均匀分布）
    uniform_prob_trans = 1 / len(trans_range)
    # 创建新的概率向量
    trans_sv_prob = [uniform_prob_trans] * len(trans_range)
    
    len_SV_trans = trans_range.tolist()
    condition_dist_sv_trans = trans_sv_prob
    
    #! 避免重复计算
    pai_pro_tem = [diff_ins_prob_del_real,diff_ins_prob_mis_real, diff_ins_prob_correct_real]
    pai_pro_tem_ = list(np.array(pai_pro_tem)/sum(pai_pro_tem))

    #! Translocation
    #args.sv_trans = 3
    print('Trans:'+str(args.sv_trans))

    for s in range(0, args.sv_trans):
        ###loop for each translocation
        ###the deletion part of long trans: cut A
        r_s = sample(unblock_region_sv, 1)[0]
        l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_trans)
        l_s_index = list(l_s_vec).index(1)
        l_s = len_SV_trans[int(l_s_index)]
        #l_s = np.random.poisson(lamnuma_del, 1)[0]

        #resample a start point and length if the long trans exceeds the terminal or previous deletions and this deletion overlaps
        circular_count_trans=0
        # indicator: if resampling exceeds 50 times
        circular_count_trans_break = 0
        while (r_s+l_s > len_seg_refine) or (r_s+l_s-1 not in unblock_region_sv) or ((unblock_region_sv.index(r_s+l_s-1)-unblock_region_sv.index(r_s)<l_s-1)):
            r_s = sample(unblock_region_sv, 1)[0]
            l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_trans)
            l_s_index = list(l_s_vec).index(1)
            l_s = len_SV_trans[int(l_s_index)]
            # count the times of resampling
            circular_count_trans = circular_count_trans + 1
            if circular_count_trans>args.times:
                circular_count_trans_break = 1
                print('Warning: No.'+str(s)+ " translocation sampling times exceeds "+str(times)+' times. Try: reduce number of variations, increase times or reference length.')
                break
            else:
                circular_count_trans_break = 0
            #l_s = np.random.poisson(lamnuma_del, 1)[0]

        if not circular_count_trans_break:
            circular_count_trans=0
            ## update the index set after choosing A (cut)
            #for next deletion, cannot be right or left neighbored
            unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s+1)))
        
            #!update for segments
        
            
            #! collect cutted pos
            location_del_trans.append(r_s)
            
            p_rand_trans_b = np.random.rand(1)
            ###rato_b_trans: prob of balanced trans
            # unbalanced trans: cut A and paste to B
            if p_rand_trans_b>ratio_b_trans:### original place is also deleted (cut A)
                ### select one possible insertion sites: ins pos of B
                #ins_trans_loc = sample(left_del_region_sv,1)[0]
                ins_trans_loc = sample(unblock_region_sv,1)[0]
                #inserted segements are the cut A region sampled before
                #inserted pos is what we sampled for unbalanced trans
                #paste region A to pos B (ins_trans_loc)
                Ins_dic_sv[ins_trans_loc] = copy.deepcopy(''.join(tem_seq_post[r_s:(r_s+l_s)]))

                #deletion and inserted segments cannot be left, right neighbor and overlapped.
                #remove positions around region B
                unblock_region_sv = list(set(unblock_region_sv) -set(range(ins_trans_loc,ins_trans_loc+2)))
                #inserted pos cannot be repeated: remove ins pos of B
                #left_del_region_sv = list(set(left_del_region_sv)-{ins_trans_loc})
                #! substitution cannot be on the inserted pos
                # if ins_trans_loc in undel_region_sv:
                #     undel_region_sv.remove(ins_trans_loc)
                if ins_trans_loc in unblock_region_sv:
                    unblock_region_sv.remove(ins_trans_loc)
                
                #collect inserted pos
                location_insert_trans.append(ins_trans_loc)
                #write SV table
                # SV_table.loc[SV_loop] = [SV_loop,ll_c,'Translocation',r_s,r_s+l_s-1,l_s,ins_trans_loc,ins_trans_loc,l_s,0,0,0,0,0]### cut and paste
                SV_table[SV_loop] = ['Translocation',r_s,r_s+l_s-1,l_s,ins_trans_loc,ins_trans_loc,l_s,0]### cut and paste
                SV_loop = SV_loop + 1
                
                # VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s:r_s+l_s]) + tem_seq_post[r_s+l_s], tem_seq_post[r_s+l_s], '.', 'PASS', 'SVTYPE=DEL;SVLEN='+str(l_s)+';CSV_TYPE=unbalanedTrans;CSV_INDEX='+str(CSV_loop)]
                VCF_table[VCF_loop] = [r_s, ''.join(tem_seq_post[r_s:r_s+l_s]) + tem_seq_post[r_s+l_s], tem_seq_post[r_s+l_s], 'SVTYPE=DEL;SVLEN='+str(l_s)+';CSV_TYPE=unbalanedTrans;CSV_INDEX='+str(CSV_loop)]
                VCF_loop = VCF_loop+ 1

                # VCF_table.loc[VCF_loop] = [str(chr_id), str(ins_trans_loc), 'rs' + str(VCF_loop), tem_seq_post[ins_trans_loc], tem_seq_post[ins_trans_loc] + ''.join(tem_seq_post[r_s:(r_s+l_s)]), '.', 'PASS', 'SVTYPE=INS;SVLEN='+str(l_s)+';CSV_TYPE=unbalanedTrans;CSV_INDEX='+str(CSV_loop)]
                VCF_table[VCF_loop] = [ins_trans_loc, tem_seq_post[ins_trans_loc], tem_seq_post[ins_trans_loc] + ''.join(tem_seq_post[r_s:(r_s+l_s)]), 'SVTYPE=INS;SVLEN='+str(l_s)+';CSV_TYPE=unbalanedTrans;CSV_INDEX='+str(CSV_loop)]
                VCF_loop = VCF_loop + 1
                
                CSV_loop = CSV_loop + 1
                
                tem_seq_post[r_s:(r_s+l_s)] = '-'*l_s
            else:
                ### balanced translocation: cut A and paste B, cut B and paste A
                ### select deletion start of B
                r_s2 = sample(unblock_region_sv,1)[0]
                # sample length of cut B
                l_s_vec2 = np.random.multinomial(n=1, pvals=condition_dist_sv_trans)
                l_s_index2 = list(l_s_vec2).index(1)
                #length of tanslocation region B
                l_s2 = len_SV_trans[int(l_s_index2)]
                #resample a start point and length if the long trans exceeds the terminal or previous deletions and B overlaps
                circular_count_trans=0
                while (r_s2+l_s2>len_seg_refine) or (r_s2+l_s2-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s2+l_s2-1)-unblock_region_sv.index(r_s2)<l_s2-1):
                    r_s2 = sample(unblock_region_sv,1)[0]
                    l_s_vec2 = np.random.multinomial(
                        n=1, pvals=condition_dist_sv_trans)
                    l_s_index2 = list(l_s_vec2).index(1)
                    l_s2 = len_SV_trans[int(l_s_index2)]
                    # count the times of resampling
                    circular_count_trans = circular_count_trans + 1
                    if circular_count_trans>args.times:
                        circular_count_trans_break = 1
                        print('Warning: No.'+str(s)+ " balanced translocation sampling times exceeds "+str(times)+' times. Try: reduce number of variations, increase times or reference length.')
                        break
                    else:
                        circular_count_trans_break = 0

                if not circular_count_trans_break:
                    circular_count_trans = 0
                    #inserted pos is one smaller than deleted base for balanced trans
                    #paste region B to pos A (r_s-1)
                    Ins_dic_sv[r_s-1] = copy.deepcopy(''.join(tem_seq_post[r_s2:(r_s2+l_s2)]))
                    #paste region A to pos B (r_s2-1)
                    Ins_dic_sv[r_s2-1] = copy.deepcopy(''.join(tem_seq_post[r_s:(r_s+l_s)]))
                    #later del no left pos, interval and right pos
                    unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s2-1,r_s2+l_s2+1)))
                    # #later SNP no in del part
                    # undel_region_sv = list(set(undel_region_sv) -set(range(r_s2,r_s2+l_s2)))
                    # #later ins no left neighbor
                    # left_del_region_sv = list(set(left_del_region_sv)-set(range(r_s2-1,r_s2+l_s2)) )
                    #!update for segments
                    #del_sites2=list(set(range(r_s2-1,r_s2+l_s2+1)))
                    # for index_seg2 in range(number_seg):
                    #     segments[index_seg2]=list(set(segments[index_seg2])-set(del_sites2))

                    #! collect cutted pos
                    location_del_trans.append(r_s2)

                    # collect inserted locations
                    location_insert_trans.append(r_s-1)
                    location_insert_trans.append(r_s2-1)
                    
                    ins_trans_loc1 = r_s2-1
                    inserted_string1 = copy.deepcopy(''.join(tem_seq_post[r_s:(r_s+l_s)]))
                    
                    ins_trans_loc2 = r_s-1
                    inserted_string2 = copy.deepcopy(''.join(tem_seq_post[r_s2:(r_s2+l_s2)]))

                    # write SV table
                    #SV_table.loc[SV_loop] = [SV_loop,ll_c,'Translocation',r_s,r_s+l_s-1,l_s,r_s2,r_s2+l_s2-1,l_s2,1,0,0,0,0]### copy and paste
                    SV_table[SV_loop] = ['Translocation',r_s,r_s+l_s-1,l_s,r_s2,r_s2+l_s2-1,l_s2,1]### copy and paste
                    SV_loop = SV_loop + 1
                    
                    #VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s:r_s+l_s]) + tem_seq_post[r_s+l_s], tem_seq_post[r_s+l_s], '.', 'PASS', 'SVTYPE=DEL;SVLEN='+str(l_s)+';CSV_TYPE=balanedTrans;CSV_INDEX='+str(CSV_loop)]
                    VCF_table[VCF_loop] = [r_s, ''.join(tem_seq_post[r_s:r_s+l_s]) + tem_seq_post[r_s+l_s], tem_seq_post[r_s+l_s], 'SVTYPE=DEL;SVLEN='+str(l_s)+';CSV_TYPE=balanedTrans;CSV_INDEX='+str(CSV_loop)]
                    VCF_loop = VCF_loop+ 1
                    
                    #VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s2:r_s2+l_s2]) + tem_seq_post[r_s2+l_s2], tem_seq_post[r_s2+l_s2], '.', 'PASS', 'SVTYPE=DEL;SVLEN='+str(l_s2)+';CSV_TYPE=balanedTrans;CSV_INDEX='+str(CSV_loop)]
                    VCF_table[VCF_loop] = [r_s, ''.join(tem_seq_post[r_s2:r_s2+l_s2]) + tem_seq_post[r_s2+l_s2], tem_seq_post[r_s2+l_s2], 'SVTYPE=DEL;SVLEN='+str(l_s2)+';CSV_TYPE=balanedTrans;CSV_INDEX='+str(CSV_loop)]
                    VCF_loop = VCF_loop+ 1

                    #VCF_table.loc[VCF_loop] = [str(chr_id), str(ins_trans_loc1), 'rs' + str(VCF_loop), tem_seq_post[ins_trans_loc1], tem_seq_post[ins_trans_loc1] + inserted_string1, '.', 'PASS', 'SVTYPE=INS;SVLEN='+str(l_s)+';CSV_TYPE=balanedTrans;CSV_INDEX='+str(CSV_loop)]
                    VCF_table[VCF_loop] = [ins_trans_loc1, tem_seq_post[ins_trans_loc1], tem_seq_post[ins_trans_loc1] + inserted_string1, 'SVTYPE=INS;SVLEN='+str(l_s)+';CSV_TYPE=balanedTrans;CSV_INDEX='+str(CSV_loop)]
                    VCF_loop = VCF_loop + 1
                    
                    #VCF_table.loc[VCF_loop] = [str(chr_id), str(ins_trans_loc2), 'rs' + str(VCF_loop), tem_seq_post[ins_trans_loc2], tem_seq_post[ins_trans_loc2] + inserted_string2, '.', 'PASS', 'SVTYPE=INS;SVLEN='+str(l_s2)+';CSV_TYPE=balanedTrans;CSV_INDEX='+str(CSV_loop)]
                    VCF_table[VCF_loop] = [ins_trans_loc2, tem_seq_post[ins_trans_loc2], tem_seq_post[ins_trans_loc2] + inserted_string2, 'SVTYPE=INS;SVLEN='+str(l_s2)+';CSV_TYPE=balanedTrans;CSV_INDEX='+str(CSV_loop)]
                    VCF_loop = VCF_loop + 1
                    
                    CSV_loop = CSV_loop + 1
                    
                    
                    tem_seq_post[r_s2:(r_s2+l_s2)] = '-'*l_s2

                    tem_seq_post[r_s:(r_s+l_s)] = '-'*l_s
                else:
                    break
        else:
            break

    #! Inversion
    ## long inversion we need get rid of the condition that overlapped the region of inserted trans
    #p_rand_inver_sv = np.random.rand(1)
    #args.sv_inver =  1
    print('Inversion:'+str(args.sv_inver))

    for inver_num in range (0, args.sv_inver):
        #sample the start and length of inversion region
        r_s = sample(unblock_region_sv, 1)[0]
        l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_inver)
        l_s_index = list(l_s_vec).index(1)
        l_s = len_SV_inver[int(l_s_index)]
        #l_s = np.random.poisson(lamnuma_del, 1)[0]

        # count the number of times that we resample
        circular_count_inver = 0
        # indicator: if resampling exceeds 50 times
        circular_count_inver_break = 0
        ### situation when resample is needed
        ### situation when resampling is needed
        ### situation 1: if inversion exceeds the end of real_con1 
        ### situation 2: if end of inversion is not in undel_region_sv
        ### situation 3: overlap of previous deletions, insertions (translocations) and inversions

        while (r_s+l_s>len_seg_refine) or (r_s+l_s-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s+l_s-1)-unblock_region_sv.index(r_s)<l_s-1):
            # select the possibile deleltion point
            r_s = sample(unblock_region_sv, 1)[0]
            l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_inver)
            l_s_index = list(l_s_vec).index(1)
            l_s = len_SV_inver[int(l_s_index)]
            # count the times of resampling
            circular_count_inver = circular_count_inver + 1
            if circular_count_inver>args.times:
                circular_count_inver_break = 1
                print("Warning: No."+str(inver_num)+ " inversion sampling times exceeds "+str(times)+' times. Try: reduce number of variations, increase times or reference length.')
                break
            else:
                circular_count_inver_break = 0

        if not circular_count_inver_break:
            circular_count_inver = 0
            unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s+1)))
            # undel_region_sv = list(set(undel_region_sv) -set(range(r_s,r_s+l_s)))
            # left_del_region_sv = list(set(left_del_region_sv)-set(range(r_s-1,r_s+l_s)) )

            location_del_inv.append(r_s)
            
            original_string = copy.deepcopy(real_con1[r_s:r_s+l_s])
            #original_string_reverse = original_string[::-1]
            original_string_reverse = DNA_complement(original_string[::-1])
            for ll_rever in range(r_s,r_s+l_s):
                tem_seq_post[ll_rever] = copy.deepcopy(original_string_reverse[int(ll_rever-r_s)])

            # SV_table.loc[SV_loop] = [SV_loop,ll_c,'Inversion',r_s,r_s+l_s-1,l_s,r_s,r_s+l_s-1,l_s,-1,0,0,0,0]
            # SV_loop = SV_loop + 1
            
            # VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), ''.join(original_string), ''.join(original_string_reverse), '.', 'PASS', 'SVTYPE=INV;SVLEN='+str(l_s)]
            # VCF_loop= VCF_loop + 1
            
            Unified_table[loop_index] = ['Inversion',r_s,r_s+l_s-1,l_s,r_s,r_s+l_s-1,l_s,-1,''.join(original_string), ''.join(original_string_reverse)]
            loop_index = loop_index + 1
        else:
            break
        
    #! Duplication
    #args.sv_dup= 2
    print('DUP:'+str(args.sv_dup))
    # 打印所有被选择的位点的信息
    inserted_site_collection = location_insert_trans
    all_selected_dup_SV= sample(unblock_region_sv,args.sv_dup)
    # 打印所有被选择的位点的信息

    dup_site_collection=all_selected_dup_SV
    #print('Dup sites'+str(dup_site_collection))     

    # other_sites = list(set(left_del_region_sv)-set(inserted_site_collection)-set(dup_site_collection))
    other_sites = list(set(unblock_region_sv)-set(inserted_site_collection)-set(dup_site_collection))
    #duplication: copy and paste
    for remain_index2 in dup_site_collection:##others### questions

        tem_copied_base = np.random.choice(dup_range)

        circular_count_dup = 0
        circular_count_dup_break = 0
        # Ensure the copied region is entirely within unblock_region_sv
        while len(set(range(remain_index2-tem_copied_base+1,remain_index2+1))&set(unblock_region_sv)) < tem_copied_base:
            # select the possibile copied length 
            #sites that do note have general insertions

            remain_index2 = sample(other_sites,1)[0]

            tem_copied_base = np.random.choice(dup_range)
        
            circular_count_dup = circular_count_dup + 1
            if circular_count_dup>args.times:
                circular_count_dup_break = 1
                print("Warning: duplication sampling times exceeds "+str(times)+' times. Try: reduce number of variations, increase times or reference length.')
                break
            else:
                circular_count_dup_break = 0

        #to decide type of duplication
        p_re_dup_random_num = np.random.rand(1)

        #ratio_re_dup:neighbor duplication (right after the copied area)
        if not circular_count_dup_break:
            circular_count_dup=0
            #! Remove the copied region from unblock_region_sv (关键修改点)
            unblock_region_sv = list(set(unblock_region_sv) - set(range(remain_index2 - tem_copied_base + 1, remain_index2 + 1)))
            #neighbor duplication
            if p_re_dup_random_num<ratio_re_dup:
                #later ins cannot be sampled from this point in other_sites
                if remain_index2 in other_sites:
                    other_sites.remove(remain_index2)
                # ins and del cannot neighbor
                if remain_index2 in unblock_region_sv:
                    unblock_region_sv.remove(remain_index2)
                if remain_index2+1 in unblock_region_sv:
                    unblock_region_sv.remove(remain_index2+1)
                #! substitution cannot be on the inserted pos
                # if remain_index2 in undel_region_sv:
                #     undel_region_sv.remove(remain_index2)
                #! 0412
                copied_region = range(remain_index2 - tem_copied_base + 1, remain_index2 + 1)
                unblock_region_sv = [x for x in unblock_region_sv if x not in copied_region]
                    
                location_insert_dup.append(remain_index2)

                Ins_dic_sv[remain_index2] = ''.join(tem_seq_post[remain_index2-tem_copied_base+1:remain_index2+1])
                ins_len_dup = len(Ins_dic_sv[remain_index2])
                
                #SV_table.loc[SV_loop] = [SV_loop,ll_c,'Duplication',remain_index2-tem_copied_base+1,remain_index2,tem_copied_base,remain_index2,remain_index2,ins_len_dup,-1,0,0,0,0]
                SV_table[SV_loop] = ['Duplication',remain_index2-tem_copied_base+1,remain_index2,tem_copied_base,remain_index2,remain_index2,ins_len_dup,-1]
                SV_loop = SV_loop + 1
                
                #VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index2), 'rs' + str(VCF_loop), tem_seq_post[remain_index2], tem_seq_post[remain_index2]+''.join(tem_seq_post[remain_index2-tem_copied_base+1:remain_index2+1]), '.', 'PASS', 'SVTYPE=DUP;SVLEN='+str(tem_copied_base)]
                VCF_table[VCF_loop] = [remain_index2, tem_seq_post[remain_index2], tem_seq_post[remain_index2]+''.join(tem_seq_post[remain_index2-tem_copied_base+1:remain_index2+1]), 'SVTYPE=DUP;SVLEN='+str(tem_copied_base)]
                VCF_loop= VCF_loop + 1
            #remote duplication
            else:
                #sample a pos for pasted region
                remain_index22 = sample(other_sites,1)[0]
                
                location_insert_dup.append(remain_index22)
                
                other_sites.remove(remain_index22)
                # if remain_index22 in left_del_region_sv:
                #     left_del_region_sv.remove(remain_index22)
                if remain_index22 in unblock_region_sv:
                    unblock_region_sv.remove(remain_index22)
                if remain_index22+1 in unblock_region_sv:
                    unblock_region_sv.remove(remain_index22+1)
                #! substitution cannot be on the inserted pos
                # if remain_index22 in undel_region_sv:
                #     undel_region_sv.remove(remain_index22)
                Ins_dic_sv[remain_index22] = ''.join(tem_seq_post[remain_index2-tem_copied_base+1:remain_index2+1])
                ins_len_dup = len(Ins_dic_sv[remain_index22])

                #SV_table.loc[SV_loop] = [SV_loop,ll_c,'Duplication',remain_index2-tem_copied_base+1,remain_index2,tem_copied_base,remain_index22,remain_index22,ins_len_dup,-1,0,0,0,0]
                SV_table[SV_loop] = ['Duplication',remain_index2-tem_copied_base+1,remain_index2,tem_copied_base,remain_index22,remain_index22,ins_len_dup,-1]
                SV_loop = SV_loop + 1
                
                #VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index22), 'rs' + str(VCF_loop), tem_seq_post[remain_index22], tem_seq_post[remain_index22]+''.join(tem_seq_post[remain_index2-tem_copied_base+1:remain_index2+1]), '.', 'PASS', 'SVTYPE=DUP;SVLEN=CSV_TYPE=DisDup;CSV_INDEX='+str(CSV_loop)]
                VCF_table[VCF_loop] = [remain_index22,tem_seq_post[remain_index22], tem_seq_post[remain_index22]+''.join(tem_seq_post[remain_index2-tem_copied_base+1:remain_index2+1]), 'SVTYPE=DUP;SVLEN=CSV_TYPE=DisDup;CSV_INDEX='+str(CSV_loop)]
                VCF_loop= VCF_loop + 1
                
                CSV_loop = CSV_loop+1  

    # Save results to temporary files
    # np.save(f'tmp/SV_table.npy', SV_table)
    # np.save(f'tmp/VCF_table.npy', VCF_table)
    del real_con1
    #! end CSV
    
    #! 把三个区域分段，最后取并集
    #! 存储方式改为numpy数组

    end_time0 = time.time()

    # process_end_time0 = time.process_time()

    elapsed_time0 = end_time0 - start_time0
    formatted_time0 = str(timedelta(seconds=elapsed_time0))

    # process_formatted_time0 = process_end_time0 - process_time0
    # process_0_time = str(timedelta(seconds=process_formatted_time0))

    print(f"Max disk usage during writing parameters was: {max_disk_usage.value}%")
    max_disk_usage.value = 0.0  # 重置最大磁盘使用率
    print(f"Initialization used：{formatted_time0}")
    # print(f"Initialization, Trans执行实际时间：{process_0_time}")

    start_time1 = time.time()
    # process_time1 = time.process_time()

    # 初始化一个空的numpy数组来存储交集
    unblock_region_vec = np.empty(number_seg, dtype=object)
    # 将 Python 列表转换为 numpy 数组
    unblock_region_sv = np.array(unblock_region_sv)


    # 确保 numpy 数组已排序
    unblock_region_sv.sort()

    unblock_vec_lengths = []

    for i, segment in enumerate(segments_initial):
        unblock_region_vec[i] = segment[np.isin(segment, unblock_region_sv)]

        # 记录 unblock_region_vec[i] 的长度
        unblock_vec_lengths.append(len(unblock_region_vec[i]))

        # 将 unblock_region_vec[i] 存储到临时文件
        np.save(os.path.join(tmp_dir, 'unblock_region_vec_{}.npy'.format(i)), unblock_region_vec[i])

    seg_probabilities = [p_length/sum(unblock_vec_lengths) for p_length in unblock_vec_lengths]
    
    # mode = args.mode # or probability weights
    mode = normalize_mode(args.mode)
    #!
    # # 读取BED文件
    # #df_hg002 = pd.read_csv('~/data/test_data/TGS/hg002/SV_Tier1_tel_dist_update.bed', sep='\t', header=None, names=['chr', 'pos', 'SV type','chr length','TR Group'])
    # indel_url = args.indel_input_bed
    # df_chr = pd.read_csv(indel_url, sep='\t', header=None, names=['pos', 'SV type'])
    # # 对'pos'列进行排序
    # df_chr = df_chr.sort_values('pos')

    # # 检查'pos'的最大值是否超过'chr_length'
    # max_pos = df_chr['pos'].max()
    # if max_pos > chr_length:
    #     print("Warning: The maximum position in the input file exceeds the chromosome length.")
    # #df_chr = df_hg002[df_hg002['chr'] == 'chr'+str(chr_id)]
    # #chr_length = len(real_con1)

    # # 计算'pos'间隔为5e5bp的Bins的SV综合
    # INS_df = df_chr[df_chr['SV type'] == 'INS'].copy()  # 创建副本
    # #INS_df['dist_bin'] = pd.cut(INS_df['pos'], bins=bins)
    # INS_df.loc[:, 'dist_bin'] = pd.cut(INS_df['pos'], bins=bins)
    # INS_df_grouped = INS_df.groupby('dist_bin').size()
    # ins_empirical = list(INS_df_grouped.values)

    # del INS_df, INS_df_grouped

    # # 计算'pos'间隔为5e5bp的Bins的SV综合
    # DEL_df = df_chr[df_chr['SV type'] == 'DEL'].copy()
    # #DEL_df['dist_bin'] = pd.cut(DEL_df['pos'], bins=bins)
    # DEL_df.loc[:, 'dist_bin'] = pd.cut(DEL_df['pos'], bins=bins)
    # DEL_df_grouped = DEL_df.groupby('dist_bin').size()

    # del_empirical = list(DEL_df_grouped.values)

    # del DEL_df, DEL_df_grouped
    
    chr_name = args.hg38
    #! cell population mean per segment
    # 样本名列表
    file_list = ['CHM1','CHM13','HG00268','HG00514','HG00733','NA19240','HG02818','HG01352','HG02059','NA12878','HG04217','HG02106','HG00268','AK1','HX1']

    # 为其他样本创建一个列名列表
    column_names_others = ['chr', 'pos', 'pos_end', 'SV type','SV length']

    # 读取索引文件,chr_length
    index_file=pd.read_csv(main_url_empirical+'hg38_index.bed', sep='\t', header=None, names=['chr', 'length', 'offset','linebases','linewidth'])

    # 初始化一个空列表来保存所有样本的结果
    INS_all_samples_tel_dist = []
    DEL_all_samples_tel_dist = []

    # 循环处理每个文件
    for file_name in file_list:
        # # 根据样本名选择列名
        # if file_name == 'HG002_INS_hg38':
        #     column_names = column_names_HG002
        #     df = pd.read_csv(f'/home/adduser/data/test_data/TGS/cell_16samples/{file_name}.bed', sep='\t', header=None, names=column_names)
        # else:
        column_names = column_names_others
        df = pd.read_csv(main_url_empirical+f'{file_name}.bed', sep='\t', header=None, names=column_names)
        
        # 只处理chr1
        df = df[df['chr'] == chr_name].copy()
        chr_length = index_file[index_file['chr'] == chr_name]['length'].iloc[0]

        df_INS = df[df['SV type'] == 'INS']
        df_DEL = df[df['SV type'] == 'DEL']
        # 将结果添加到列表中
        INS_all_samples_tel_dist.append(df_INS['pos'].tolist())
        DEL_all_samples_tel_dist.append(df_DEL['pos'].tolist())

        # 保存结果到.npy文件
        file_path_INS = os.path.join(tmp_dir, f'INS_pos_{file_name}_tel_dist.npy')
        np.save(file_path_INS, df_INS['pos'].tolist())
        
        # 保存结果到.npy文件Y
        file_path_DEL = os.path.join(tmp_dir, f'DEL_pos_{file_name}_tel_dist.npy')
        np.save(file_path_DEL, df_DEL['pos'].tolist())
        # 将 unblock_region_vec[i] 存储到临时文件
        
    # 创建一个空的列表来保存所有样本的数据
    all_data_INS = []
    all_data_DEL = []
    
    # 循环处理每个文件
    for file_name in file_list:
        # 构造文件路径
        file_path_INS = os.path.join(tmp_dir, f'INS_pos_{file_name}_tel_dist.npy')
        # 读取.npy文件
        data_INS = np.load(file_path_INS)
        
        # 计算每个bins中的点的数量
        counts, _ = np.histogram(data_INS, bins=bins)
        
        # 将数据保存到一个与file_name相关的变量中
        vars()[f'{file_name}_INS_counts'] = counts.tolist()

        # 获取每个bins中的点的数量
        data = vars()[f'{file_name}_INS_counts']
        
        # 将数据添加到all_data列表中
        all_data_INS.append(data)
        
        # 构造文件路径
        file_path_DEL = os.path.join(tmp_dir, f'DEL_pos_{file_name}_tel_dist.npy')
        # 读取.npy文件
        data_DEL = np.load(file_path_DEL)
        
        # 计算每个bins中的点的数量
        counts, _ = np.histogram(data_DEL, bins=bins)
        
        # 将数据保存到一个与file_name相关的变量中
        vars()[f'{file_name}_DEL_counts'] = counts.tolist()

        # 获取每个bins中的点的数量
        data = vars()[f'{file_name}_DEL_counts']
        
        # 将数据添加到all_data列表中
        all_data_DEL.append(data)

    # 转换为numpy数组
    all_data_INS = np.array(all_data_INS)

    # 计算每个bins的均值和标准差
    means_INS = np.mean(all_data_INS, axis=0)
    std_devs_INS = np.std(all_data_INS, axis=0)
    # 创建一个空的列表来保存抽取的值
    INS_sampled_values = []

    # 循环处理每个bins
    for mean, std_dev in zip(means_INS, std_devs_INS):
        # 从正态分布中抽取一个值
        value = np.random.normal(mean, std_dev)
        
        # 整数化并确保最小值为0
        value = max(0, int(value))
        
        # 将值添加到sampled_values列表中
        INS_sampled_values.append(value)

    
    # ins_empirical = INS_sampled_values

    # 转换为numpy数组
    all_data_DEL = np.array(all_data_DEL)

    # 计算每个bins的均值和标准差
    means_DEL = np.mean(all_data_DEL, axis=0)
    std_devs_DEL = np.std(all_data_DEL, axis=0)
    # 创建一个空的列表来保存抽取的值
    DEL_sampled_values = []

    # 循环处理每个bins
    for mean, std_dev in zip(means_DEL, std_devs_DEL):
        # 从正态分布中抽取一个值
        value = np.random.normal(mean, std_dev)
        
        # 整数化并确保最小值为0
        value = max(0, int(value))
        
        # 将值添加到sampled_values列表中
        DEL_sampled_values.append(value)

    #! 0409 区分mode
    mode = normalize_mode(args.mode)

    if mode == 'empirical':
        # 使用总体均值
        ins_empirical = means_INS
        del_empirical = means_DEL
    elif mode == 'probability':
        # 从正态分布中抽样
        INS_sampled_values = []
        for mean, std_dev in zip(means_INS, std_devs_INS):
            value = max(0, int(np.random.normal(mean, std_dev)))
            INS_sampled_values.append(value)
        ins_empirical = INS_sampled_values
        
        DEL_sampled_values = []
        for mean, std_dev in zip(means_DEL, std_devs_DEL):
            value = max(0, int(np.random.normal(mean, std_dev)))
            DEL_sampled_values.append(value)
        del_empirical = DEL_sampled_values
    else:
        # 默认使用总体均值
        print(f"Warning: Invalid mode '{mode}'. Using empirical means.")
        ins_empirical = means_INS
        del_empirical = means_DEL
    
    # del_empirical = DEL_sampled_values
    
    # 对概率进行归一化
    #args.sum
    # if len(del_empirical) == number_seg:
    #     if mode == 'probability':
    #         del_probabilities = [p_del/sum(del_empirical) for p_del in del_empirical]
    #         if args.sum:
    #             total_del_sv = sum(del_empirical)
    #         else:
    #             total_del_sv = args.sv_del
    #         print('Total del:'+str(total_del_sv))
    #         del_SV_per_segment = np.random.multinomial(total_del_sv, del_probabilities)

    #     elif mode == 'empirical':
    #         total_del_sv = sum(del_empirical)
    #         del_SV_per_segment = del_empirical
    #         print('Total del:'+str(total_del_sv))
    #     else:
    #         print('MODE TYPE ERROR')
    # else:
    #     print('warning: empirical del info wrong dimension')

    # if len(ins_empirical) == number_seg:
    #     if mode == 'probability':
    #         ins_probabilities = [p_ins/sum(ins_empirical) for p_ins in ins_empirical]
    #         if args.sum:
    #             total_ins_sv = sum(ins_empirical)
    #         else:
    #             total_ins_sv = args.sv_ins
    #         print('Total ins:'+str(total_ins_sv))
    #         ins_SV_per_segment = np.random.multinomial(total_ins_sv, ins_probabilities)
    #     elif mode == 'empirical':
    #         total_ins_sv = sum(ins_empirical)
    #         ins_SV_per_segment = ins_empirical
    #         print('Total ins:'+str(total_ins_sv))
    #     else:
    #         print('MODE TYPE ERROR')
    # else:
    #     print('warning: empirical ins info wrong dimension')

    # 对概率进行归一化
    # if len(del_empirical) == number_seg:
    #     if mode == 'probability':
    #         del_probabilities = [p_del/sum(del_empirical) for p_del in del_empirical]
    #         if args.sum:
    #             total_del_sv = sum(del_empirical)
    #         else:
    #             total_del_sv = args.sv_del
    #         print('Total del:'+str(total_del_sv))
    #         # 生成整数结果并四舍五入
    #         del_SV_per_segment = np.round(np.random.multinomial(total_del_sv, del_probabilities)).astype(int)
            
    #     elif mode == 'empirical':
    #         # 将empirical数据四舍五入为整数
    #         del_SV_per_segment = np.round(del_empirical).astype(int)
    #         total_del_sv = sum(del_SV_per_segment)
    #         print('Total del:'+str(total_del_sv))
    #     else:
    #         print('MODE TYPE ERROR')
    # else:
    #     print('warning: empirical del info wrong dimension')

    # if len(ins_empirical) == number_seg:
    #     if mode == 'probability':
    #         ins_probabilities = [p_ins/sum(ins_empirical) for p_ins in ins_empirical]
    #         if args.sum:
    #             total_ins_sv = sum(ins_empirical)
    #         else:
    #             total_ins_sv = args.sv_ins
    #         print('Total ins:'+str(total_ins_sv))
    #         # 生成整数结果并四舍五入
    #         ins_SV_per_segment = np.round(np.random.multinomial(total_ins_sv, ins_probabilities)).astype(int)
            
    #     elif mode == 'empirical':
    #         # 将empirical数据四舍五入为整数
    #         ins_SV_per_segment = np.round(ins_empirical).astype(int)
    #         total_ins_sv = sum(ins_SV_per_segment)
    #         print('Total ins:'+str(total_ins_sv))
    #     else:
    #         print('MODE TYPE ERROR')
    # else:
    #     print('warning: empirical ins info wrong dimension')

    # del ins_empirical, del_empirical
    
    mode = normalize_mode(args.mode)
    # Process deletions
    del_SV_per_segment, total_del_sv = process_sv_data(
        empirical_data=del_empirical,
        sv_total_arg=args.sv_del,
        data_type='del',
        number_seg=number_seg,
        args_sum=args.sum,
        mode=mode
    )

    # Process insertions
    ins_SV_per_segment, total_ins_sv = process_sv_data(
        empirical_data=ins_empirical,
        sv_total_arg=args.sv_ins,
        data_type='ins',
        number_seg=number_seg,
        args_sum=args.sum,
        mode=mode
    )

    # Handle cases where processing failed
    if del_SV_per_segment is None:
        print("Error: Invalid deletion data. Please check your input.")
        del_SV_per_segment = [0] * number_seg  # Default to zero counts
        total_del_sv = 0

    if ins_SV_per_segment is None:
        print("Error: Invalid insertion data. Please check your input.")
        ins_SV_per_segment = [0] * number_seg  # Default to zero counts
        total_ins_sv = 0
        
    # del_SV_per_segment = np.random.multinomial(args.sv_del, seg_probabilities)
    # ins_SV_per_segment = np.random.multinomial(args.sv_ins, seg_probabilities)
    print('DEL counts per segment:'+str(del_SV_per_segment))
    print('INS counts per segment:'+str(ins_SV_per_segment))

    len_unblock_region = len(unblock_region_sv)
 
    if 0 <= args.snv_del < 1:
        del_snv_number = int(len_unblock_region * args.snv_del)
    else:
        del_snv_number = min(int(args.snv_del), len_unblock_region)

    if 0 <= args.snv_ins < 1:
        ins_snv_number = int(len_unblock_region * args.snv_ins)
    else:
        ins_snv_number = min(int(args.snv_ins), len_unblock_region)

    if 0 <= args.snp < 1:
        snp = int(len_unblock_region * args.snp)
    else:
        snp = min(int(args.snp), len_unblock_region)

    print('Small del:'+str(int(del_snv_number)))
    print('Small ins:'+str(int(ins_snv_number)))
    print('SNP:'+str(int(snp)))
    del_snv_per_segment = np.random.multinomial(del_snv_number, seg_probabilities)
    ins_snv_per_segment = np.random.multinomial(ins_snv_number, seg_probabilities)     
    snp_per_segment = np.random.multinomial(snp, seg_probabilities)  
    
    # 定义概率
    #! TR region
    p_del_TR = args.p_del_region
    p_ins_TR = args.p_ins_region
    
    #TR_region_bed = '~/data/test_data/TGS/hg002/chr21_TR_unique.bed'
    TR_region_bed = args.region_bed_url
    print('Assign DEL prob. of ' + str(p_del_TR) + '; INS prob of ' + str(p_ins_TR) + ' for the input region')
    # 读取BED文件
    df_TR = pd.read_csv(TR_region_bed, sep='\t', header=None, names=['start', 'end'])

    # 使用.apply()函数应用上面定义的函数，得到每一行'start'和'end'值的检查结果
    df_TR['check'] = df_TR.apply(check_start_end, axis=1)

    # 检查是否所有行的'start'值都小于'end'值
    if df_TR['check'].all():
        print("All rows of the region BED file are valid.")
    else:
        print("Error: Not all rows of the region BED file are valid.")
        sys.exit(1)


    # 使用.apply()函数应用上面定义的函数，得到每一行'start'和'end'区域中的所有点
    df_TR['points'] = df_TR.apply(get_points, axis=1)

    # 使用set.union()函数合并所有的点，得到一个包含所有点的集合
    TR_set = set.union(*df_TR['points'])

    # 现在，'TR_set'集合中包含了所有的点
    print('length of the TR region: '+str(len(TR_set)))
            
    end_time1 = time.time()

    # process_end_time1 = time.process_time()

    elapsed_time1 = end_time1 - start_time1
    formatted_time1 = str(timedelta(seconds=elapsed_time1))

    # process_formatted_time1 = process_end_time1 - process_time1
    # process_first_time = str(timedelta(seconds=process_formatted_time1))

    print(f"Max disk usage during writing parameters was: {max_disk_usage.value}%")
    max_disk_usage.value = 0.0  # 重置最大磁盘使用率
    print(f"Cutting the sequence used：{formatted_time1}")

    # process_time2 = time.process_time()
    start_time2 = time.time()

    #! delete variables to release 内存
    # del ref_seq, n_positions, n_positions_set, all_positions, all_region_sv
    # del segments_initial, all_region_sv_array
    # del unblock_region_sv, undel_region_sv, left_del_region_sv
    #! define the function for each segment

    # Use the function in your code
    max_mem_usage, max_mem_usage_gb = monitor_memory(90, max_mem_usage, max_mem_usage_gb, tmp_dir)  # stop the program if memory usage exceeds 90%

    print(f"Max memory usage during the first section was: {max_mem_usage}% ({max_mem_usage_gb} GB)")

    
    


    #! end of the function
    
    max_mem_usage, max_mem_usage_gb = monitor_memory(90, max_mem_usage, max_mem_usage_gb, tmp_dir)  # stop the program if memory usage exceeds 90%

    
    #! 计算每个进程的复杂度
    #complexities = [del_SV + ins_SV for del_SV, ins_SV in zip(del_SV_per_segment, ins_SV_per_segment)]
    # # 计算每个进程的复杂度
    # complexities = [del_SV + ins_SV + length for del_SV, ins_SV, length in zip(del_SV_per_segment, ins_SV_per_segment, unblock_vec_lengths)]
    
    # # 创建一个包含进程索引和复杂度的列表
    # indexed_complexities = list(enumerate(complexities))

    # # 根据复杂度对进程进行排序
    # sorted_complexities = sorted(indexed_complexities, key=lambda x: x[1], reverse=True)
    # #number of kernel
    # pool = Pool(args.cores)
    # #xun_list=list(range(number_seg))
    # # 创建一个包含进程索引的列表，按照复杂度从大到小排序
    # xun_list = [index for index, complexity in sorted_complexities]
    # print('Optimal order of running due to computing complexity：'+str(xun_list)) 
    
    # #xun_list = list(range(number_seg))
    # # args.cores = 5  # set the number of processes equal to the number of cores

    # pool = Pool(args.cores)
    # results = []

    # # for xun in xun_list:
    # #     result = pool.apply_async(gen_consensus, args=(xun, process_dict), callback=lambda _: monitor_memory(90, max_mem_usage, max_mem_usage_gb), error_callback=error_handler)
    # #     results.append((xun, result))
    # for xun in xun_list:
    #     # 将所有需要的变量作为参数传递给gen_consensus函数
    #     result = pool.apply_async(gen_consensus, args=(xun, chr_id, ll_c, starts_seg, ends_seg, tem_seq_post, del_SV_per_segment, ins_SV_per_segment, del_snv_per_segment, ins_snv_per_segment, snp_per_segment, pai_pro_tem_, args, condition_dist_sv_del, len_SV_del, condition_dist_sv_ins, len_SV_ins, ins_selection,base_list,substitution_matrix,mis_selection,times, TR_set, p_del_TR, p_ins_TR,tmp_dir), callback=lambda _: monitor_memory(90, max_mem_usage, max_mem_usage_gb,tmp_dir), error_callback=error_handler)
    #     results.append((xun, result))

    # pool.close()
    # pool.join()


    # print(f"Max memory usage during the program run was: {max_mem_usage}% ({max_mem_usage_gb} GB)")

#! 并行运算
    # 计算每个进程的复杂度
    complexities = [del_SV + ins_SV + length for del_SV, ins_SV, length in zip(del_SV_per_segment, ins_SV_per_segment, unblock_vec_lengths)]

    # 创建一个包含进程索引和复杂度的列表
    indexed_complexities = list(enumerate(complexities))

    # 根据复杂度排序任务
    sorted_complexities = sorted(indexed_complexities, key=lambda x: x[1], reverse=True)

    # 按复杂度高低生成任务列表
    xun_list = [index for index, complexity in sorted_complexities]
    print('Optimal order of running due to computing complexity:', xun_list)

    # 记录程序开始时间
    start_time2 = time.time()

    # 使用 ProcessPoolExecutor 动态分配任务
    results = {}  # 用于存储结果，键为 xun，值为任务返回的结果
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = {
            executor.submit(
                gen_consensus,
                xun, starts_seg, ends_seg, tem_seq_post,
                del_SV_per_segment, ins_SV_per_segment, del_snv_per_segment,
                ins_snv_per_segment, snp_per_segment, pai_pro_tem_, args,
                condition_dist_sv_del, len_SV_del, condition_dist_sv_ins,
                len_SV_ins, ins_selection, base_list, substitution_matrix,
                mis_selection, times, tmp_dir,  TR_set, p_del_TR, p_ins_TR
            ): xun for xun in range(len(starts_seg))
        }
    # 动态处理任务结果
        for future in as_completed(futures):
            xun = futures[future]
            try:
                # 修改返回值解包方式
                unified_table_seg, ins_dic_sv_seg, tem_seq_post_seg = future.result()
                if unified_table_seg is None or unified_table_seg.size == 0:
                    print(f"Task {xun} returned empty or None unified_table_seg.")
                results[xun] = (unified_table_seg, ins_dic_sv_seg, tem_seq_post_seg)
                print(f"Task {xun} completed successfully.")
            except Exception as e:
                print(f"Task {xun} failed with exception: {e}")
    end_time2 = time.time()
    # process_end_time2 = time.process_time()

    elapsed_time2 = end_time2 - start_time2
    formatted_time2 = str(timedelta(seconds=elapsed_time2))

    # process_formatted_time2 = process_end_time2 - process_time2
    # process_second_time = str(timedelta(seconds=process_formatted_time2))

    print(f"Parallel computing used：{formatted_time2}")
    # print(f"并行处理 CPU执行实际时间：{process_second_time}")

    # process_time3 = time.process_time()
    start_time3 = time.time()
    
    import gc
    gc.collect()


    # # 获取结果并按照 xun 的值进行排序
    # results = sorted([(xun, result.get()) for xun, result in results], key=lambda x: x[0])

    # # 提取结果
    # SV_table_segs = [result[0] for xun, result in results]
    # VCF_table_segs = [result[1] for xun, result in results]
    # # unblock_region_segs = [result[2] for xun, result in results]
    # Ins_dic_sv_segs = [result[3] for xun, result in results]
    # tem_seq_post_segs = [result[4] for xun, result in results]


    # # 优化拼接数据框的代码
    # SV_table_combined = pd.concat(SV_table_segs, ignore_index=True)
    # SV_table_merged = pd.concat([SV_table, SV_table_combined], ignore_index=True)
    # VCF_table_combined = pd.concat(VCF_table_segs, ignore_index=True)
    # VCF_table_merged = pd.concat([VCF_table, VCF_table_combined], ignore_index=True)
    # 按照 xun 的顺序提取结果
    sorted_results = [results[xun] for xun in sorted(results)]  # 确保结果按 xun 顺序排列

    # 解构子任务返回的结果
    # SV_table_segs = [data[0] for data in sorted_results]  # 每个进程的 SV 表格
    # VCF_table_segs = [data[1] for data in sorted_results]  # 每个进程的 VCF 表格
    unified_table_segs = [data[0] for data in sorted_results]  # 每个进程的 Unified_table_seg（numpy 数组）
    Ins_dic_sv_segs = [data[1] for data in sorted_results]  # 插入字典
    tem_seq_post_segs = [data[2] for data in sorted_results]  # 处理后的序列

    end_time3 = time.time()
    # process_end_time3 = time.process_time()
    start_time4 = time.time()
    
    mask_uni = np.all((Unified_table == None) | (Unified_table == ''), axis=1)
    Unified_table = Unified_table[~mask_uni]

    # 调试输出
    print("Number of segments:", len(unified_table_segs))
    for i, arr in enumerate(unified_table_segs):
        print(f"Segment {i} shape: {arr.shape}")

    # 合并所有子任务的结果
    # 1. 合并 numpy 数组
    if unified_table_segs:
        # 过滤掉空数组
        unified_table_segs = [arr for arr in unified_table_segs if arr.size > 0]
        if unified_table_segs:
            # 检查列数是否一致
            cols = Unified_table.shape[1]
            for arr in unified_table_segs:
                if arr.shape[1] != cols:
                    raise ValueError(f"Column mismatch: Expected {cols} columns, but got {arr.shape[1]} columns.")
            # 将原始的 Unified_table 和所有子任务的 unified_table_segs 拼接
            unified_table_combined = np.vstack([Unified_table] + unified_table_segs)
            print(f"Combined table shape: {unified_table_combined.shape}")
        else:
            unified_table_combined = Unified_table
            print("All segments are empty. Using original Unified_table.")
    else:
        unified_table_combined = Unified_table
        print("No new segments to combine. Using original Unified_table.")

    # 合并字典
    # 使用字典解包操作快速合并多个字典
    Ins_dic_sv_combined = {**Ins_dic_sv, **{k: v for seg in Ins_dic_sv_segs for k, v in seg.items()}}

    # 拼接 DNA 序列
    # 使用 itertools.chain 快速拼接多个列表
    from itertools import chain
    tem_seq_post_update = list(chain.from_iterable(tem_seq_post_segs))

    # 找出所有元素都是 None 或空字符串的行
    mask = np.all((SV_table == None) | (SV_table == ''), axis=1)

    # 删除这些行
    SV_table = SV_table[~mask]

    # 现在你可以将 SV_table 转换为 DataFrame
    # SV_table_merged = pd.DataFrame(SV_table, columns=['Index','Index_con','SV_type','Original_start',\
                                        # 'Original_end','Len_SV','New_start','New_end','New_len_SV','Balanced Trans Flag','relative start1','relative end1','relative start2','relative end2'])
    
    # 拆分 Unified_table
    # 提取前 8 列作为 SV_sub
    SV_sub = unified_table_combined[:, :8]  # 前 8 列

    # 合并 SV_table 和 SV_sub
    # 假设 SV_table 已经初始化并填充了数据
    SV_table_merged = np.vstack((SV_table, SV_sub))  # 垂直合并

    # 将 SV_table 转换为 DataFrame
    SV_table_merged = pd.DataFrame(SV_table_merged, columns=['SV_type', 'Original_start', 'Original_end', 
                                                 'Len_SV', 'New_start', 'New_end', 
                                                 'New_len_SV', 'Balanced Trans Flag'])
    
    # 添加 'Index_con' 列并填充 ll_c 的值
    SV_table_merged['Index_con'] = ll_c  # ll_c 是你定义的变量

    # 添加 'relative start1', 'relative end1', 'relative start2', 'relative end2' 列并填充为 0
    SV_table_merged['relative start1'] = 0
    SV_table_merged['relative end1'] = 0
    SV_table_merged['relative start2'] = 0
    SV_table_merged['relative end2'] = 0

    # 定义 SV 类型的排序
    sv_type_order = ['Translocation', 'Inversion', 'Duplication', 'Deletion', 'Insertion', 'Small_Del', 'Small_Ins', 'Substitution']

    # 将 'SV type' 转换为有序的分类变量
    SV_table_merged['SV_type'] = pd.Categorical(SV_table_merged['SV_type'], categories=sv_type_order, ordered=True)

    # 按照 'SV type' 和 'Original_start' 排序
    SV_table_merged.sort_values(by=['SV_type', 'Original_start'], inplace=True)

    # 重置索引，并丢弃原来的索引
    SV_table_merged.reset_index(drop=True, inplace=True)

    # 再次重置索引，将新的索引添加为一个列，然后将这个新的列的名字改为 'Index'
    SV_table_merged.reset_index(inplace=True)
    SV_table_merged.rename(columns={'index': 'Index'}, inplace=True)

    # 调整列顺序为指定的顺序
    final_columns = ['Index', 'Index_con', 'SV_type', 'Original_start', 'Original_end', 
                    'Len_SV', 'New_start', 'New_end', 'New_len_SV', 'Balanced Trans Flag', 
                    'relative start1', 'relative end1', 'relative start2', 'relative end2']
    SV_table_merged = SV_table_merged[final_columns]

    # 提取 SV_type 和 Len_SV 列
    sv_types = unified_table_combined[:, 0]  # SV_type 列
    sv_lens = unified_table_combined[:, 3]   # Len_SV 列

    # 生成 INFO 列
    info_column = generate_info_column(sv_types, sv_lens)
    # 提取 VCF_sub
    # VCF_sub 需要包含 'POS' 和 'INFO' 列，并根据规则生成
    VCF_sub = np.empty((unified_table_combined.shape[0], 4), dtype=object)  # 初始化 VCF_sub
    VCF_sub[:, 0] = unified_table_combined[:, 1]  # 'POS' 列 = 'Original_start'
    VCF_sub[:, 1] = unified_table_combined[:, 8]  # 'REF' 列
    VCF_sub[:, 2] = unified_table_combined[:, 9]  # 'ALT' 列
    VCF_sub[:, 3] = info_column  # 'INFO' 列
    # # 根据 SV_type 和 Len_SV 生成 INFO 列
    # for i in range(Unified_table.shape[0]):
    #     sv_type = Unified_table[i, 0]  # SV_type
    #     sv_len = Unified_table[i, 3]   # Len_SV
    #     VCF_sub[i, 3] = get_info_from_sv_type(sv_type, sv_len)  # INFO 列


    # 对 VCF_table 重复相同的步骤
    mask = np.all((VCF_table == None) | (VCF_table == ''), axis=1)
    VCF_table = VCF_table[~mask]

    # 合并 VCF_table 和 VCF_sub
    # 假设 VCF_table 已经初始化并填充了数据
    VCF_table_merged = np.vstack((VCF_table, VCF_sub))  # 垂直合并


    # VCF_table_merged = pd.DataFrame(VCF_table, columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','FORMAT','SAMPLE_ID'])
    # 将 VCF_table 转换为 DataFrame
    VCF_table_merged = pd.DataFrame(VCF_table_merged, columns=['POS', 'REF', 'ALT', 'INFO'])

    # 添加 'CHROM' 列并填充为 str(chr_id) 的值
    VCF_table_merged['CHROM'] = str(chr_id)  # chr_id 是你定义的变量

    # 将 'POS' 列转换为整数类型
    VCF_table_merged['POS'] = VCF_table_merged['POS'].astype(int)

    # 按照 POS 升序排序
    VCF_table_merged = VCF_table_merged.sort_values(by='POS')

    # 添加 'ID' 列并填充从 0 开始的编号
    VCF_table_merged['ID'] = range(len(VCF_table_merged))

    # 添加 'QUAL', 'FILTER', 'FORMAT', 'SAMPLE_ID' 列并填充相应的值
    VCF_table_merged['QUAL'] = '.'
    VCF_table_merged['FILTER'] = 'PASS'
    VCF_table_merged['FORMAT'] = 'GT'
    VCF_table_merged['SAMPLE_ID'] = '1/1'

    # 重新排列列顺序（如果需要）
    VCF_table_merged = VCF_table_merged[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE_ID']]
    # 后续代码保持不变
    print("All tasks completed and results integrated successfully.")
    

    #!
    #! check the simulated count
    # Count the number of deletions in SV_table_merged
    # 按照顺序计数每种变异类型
    # deletion_count = SV_table_merged[SV_table_merged['SV_type'] == 'Deletion'].shape[0]
    # insertion_count = SV_table_merged[SV_table_merged['SV_type'] == 'Insertion'].shape[0]
    # print('Simulated deletions:'+str(deletion_count))
    # print('Simulated insertions:'+str(insertion_count))

    counts = {}
    total = 0
    
    for sv_type in sv_type_order:
        count = SV_table_merged[SV_table_merged['SV_type'] == sv_type].shape[0]
        counts[sv_type] = count
        total += count
        print(f'Simulated {sv_type.lower()}s:' + str(count))

    # 输出总数
    print('Total simulated Variations:' + str(total))

    #! 0409 图像
    # 根据输入类型确定empirical_data的格式
    
    # 多样本输入情况
    del_data = {'mean': means_DEL, 'std': std_devs_DEL}  # DEL均值和标准差
    ins_data = {'mean': means_INS, 'std': std_devs_INS}  # INS均值和标准差
    
    # 公共参数
    common_params = {
        'chr_id': chr_id,
        'chr_length': chr_length,
        'len_bins': args.len_bins,
        'n_bed_path': n_bed_path if not args.notblockN else None,
        'output_dir': args.save
    }

    # 绘制Deletion分布图
    plot_sv_distribution(mode, args,
        SV_table_merged,
        'Deletion',
        empirical_data=del_data,
        **common_params
    )

    # 绘制Insertion分布图
    plot_sv_distribution(mode, args,
        SV_table_merged,
        'Insertion',
        empirical_data=ins_data,
        **common_params
    )

    plot_sv_distribution_in_tr_region(args,
    df_sv=SV_table_merged,
    sv_type='Deletion',
    chr_id=seqname,
    chr_length=chr_length,
    len_bins=args.len_bins,
    df_TR=df_TR,
    output_dir=args.save
)
    plot_sv_distribution_in_tr_region(args,
    df_sv=SV_table_merged,
    sv_type='Insertion',
    chr_id=seqname,
    chr_length=chr_length,
    len_bins=args.len_bins,
    df_TR=df_TR,
    output_dir=args.save
)
    
        
   
        
    end_time4 = time.time()

    start_time5 = time.time()
    
    tem_seq_post_up = tem_seq_post_update.copy()
    # for idx in sorted(Ins_dic_sv_combined, reverse=True):
    #     tem_seq_post_up.insert(idx+1, Ins_dic_sv_combined[idx])
    #! 插入并行运算外产生的inserted segments
    for idx in sorted(Ins_dic_sv, reverse=True):
        tem_seq_post_up.insert(idx, Ins_dic_sv[idx])

    tem_seq_post_up_string = ''.join(tem_seq_post_up)

    tem_seq_post_up_string= tem_seq_post_up_string.replace('-','')

    updated_con.append(copy.deepcopy(tem_seq_post_up_string))
    
    print('Length of the simulated sequence: '+str(len(tem_seq_post_up_string)))
    
   
    tem_ins_dic = Ins_dic_sv_combined
    
    
    write_template_fasta_con(args, seqname, updated_con[0])
    write_vcf(args, VCF_table_merged, seqname, start_base, end_base)
    
    #! final table
    if args.write:
        print('finalize table')
        SV_table_merged = SV_write_relative(SV_table_merged,ll_c,tem_ins_dic)
        output_path = args.save +'BV_' + str(args.rep) + '_seq_' + str(seqname) + '_SVtable_full.csv'
        SV_table_merged.to_csv(output_path, header=True, index=False)
        print(f"Saved full SV table to {output_path}")
    else:
        output_path = args.save +'BV_' + str(args.rep) + '_seq' + str(seqname) + '_SVtable.csv'
        SV_table_merged.to_csv(output_path, header=True, index=False)
        print(f"Saved SV table to {output_path}")
        
        # Save the dictionary as a .npy file
        dic_path = args.save+'BV_'+str(args.rep) + '_seq_' + str(seqname) + '_tem_ins_dic.npy'
        np.save(dic_path, tem_ins_dic)
        print(f"Saved tem_ins_dic to {dic_path}")
    
    end_time5 = time.time()
    # Use the function in your code
    max_mem_usage, max_mem_usage_gb = monitor_memory(90, max_mem_usage, max_mem_usage_gb, tmp_dir)  # stop the program if memory usage exceeds 90%

    print(f"Max memory usage during the whole process run was: {max_mem_usage}% ({max_mem_usage_gb} GB)")
    # At the end of the program, remove the temporary directory
    shutil.rmtree(tmp_dir)

    #! 主程序结束后，停止磁盘监控进程
    monitor_process.terminate()

    print(f"Initialization used：{formatted_time0}")
    # print(f"Initialization执行实际时间：{process_0_time}")

    print(f"Cutting the sequence used：{formatted_time1}")
    # print(f"取交集 CPU执行实际时间：{process_first_time}")

    print(f"Parallel computing used：{formatted_time2}")

    elapsed_time3 = end_time3 - start_time3
    formatted_time3 = str(timedelta(seconds=elapsed_time3))

    # process_formatted_time3 = process_end_time3 - process_time3
    # process_third_time = str(timedelta(seconds=process_formatted_time3))

    print(f"Collecting the results time used：{formatted_time3}")
    # print(f"写出结果CPU实际运行时间：{process_third_time}")
    elapsed_time4 = end_time4 - start_time4
    formatted_time4 = str(timedelta(seconds=elapsed_time4))
    print(f"Merging the results used：{formatted_time4}")

    elapsed_time5 = end_time5 - start_time5
    formatted_time5 = str(timedelta(seconds=elapsed_time5))
    print(f"Writing the results used：{formatted_time5}")

    total_time = end_time5 - start_time0
    formatted_time6 = str(timedelta(seconds=total_time))
    print(f"Total time was：{formatted_time6}")

    # 记得在结束时关闭文件
    sys.stdout.close()
    
    
    
if __name__ == "__main__":
    main()



