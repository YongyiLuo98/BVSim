import argparse
import os
import sys
import pandas as pd
import pysam
from typing import List, Dict, Tuple
import copy
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

from scipy.stats import pareto
import time
from datetime import timedelta
# Import necessary libraries
from pandas import concat

print("Exact positions from input table")

def format_time(seconds):
    """Convert seconds to HH:MM:SS format"""
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    seconds = int(seconds % 60)
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"

# def check_and_filter_overlapping_deletions(del_df):
#     # Sort by length descending
#     del_df = del_df.sort_values('Len_SV', ascending=False).copy()
    
#     kept_regions = []
#     to_keep = []
    
#     for idx, row in del_df.iterrows():
#         current_start = row['Original_start']
#         current_end = row['Original_end']
#         overlap_found = False
        
#         # Check against kept regions using vectorized operations
#         if kept_regions:
#             # Create a DataFrame from kept_regions for vectorized comparison
#             kept_df = pd.DataFrame(kept_regions, columns=['start', 'end', 'length'])
#             # Find overlaps
#             overlaps = ((current_end >= kept_df['start']) & 
#                         (kept_df['end'] >= current_start))
            
#             if overlaps.any():
#                 overlap_found = True
#                 overlap_info = kept_df[overlaps].iloc[0]
#                 print(f"Warning: Deletion region {current_start}-{current_end} (len={current_end-current_start+1}) "
#                       f"overlaps with {overlap_info['start']}-{overlap_info['end']} (len={overlap_info['length']}). "
#                       f"Keeping the larger one.")
        
#         if not overlap_found:
#             kept_regions.append((current_start, current_end, current_end - current_start + 1))
#             to_keep.append(idx)
    
#     return del_df.loc[to_keep]

def DNA_complement(sequence):
    trantab = str.maketrans('ATCGatcg','TAGCtagc')
    string = sequence.translate(trantab)
    return string

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
        f.write('##source=uniform.py\n')
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
        f.write('##INFO=<ID=DUP,Number=1,Type=String,Description="Duplication">\n')
        # FILTER field (新增部分)
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        f.write('##FILTER=<ID=.,Description="No filter applied">\n')  # 可选 
        # FORMAT field
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        # FILTER field (新增部分)
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        f.write('##FILTER=<ID=.,Description="No filter applied">\n')  # 可选 
        # Column headers
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_ID\n')
        
        # Write data
        if df.empty:
            print("Warning: DataFrame is empty. No data will be written to VCF.")
        else:
            df.to_csv(f, sep='\t', index=False, header=False)
            
    # 打印输出文件路径
    print(f"Saved VCF file to {output_path}")
    
                
# 定义 SV_type 到缩写的映射规则
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

def parse_args():
    parser = argparse.ArgumentParser(description='BVSim - Variant Generation from User Input')
    
    # 基本参数
    parser.add_argument('-ref', type=str, required=True, help='Input reference fasta file path')
    parser.add_argument('-seq_index', type=int, default=0, 
                      help='Index of sequence to use (0-based). Default: 0 (first sequence)')
    parser.add_argument('-variant_table', type=str, required=True,
                      help='Path to variant table CSV/TSV file')
    parser.add_argument('-rep', type=int, help='Replication ID', default=99)
    parser.add_argument('-save', type=str, 
                      default=os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'save')+'/',
                      help='Output directory path')
    parser.add_argument('-seed', type=int, help='Seed for random number generator', default=999)
    parser.add_argument('-snp', type=float, help='SNV number or probability', default=5)
    parser.add_argument('-snv_del', type=float, help='SNV deletion number or probability', default=5)
    parser.add_argument('-snv_ins', type=float, help='SNV insertion number or probability', default=5)
    parser.add_argument('-notblockN', action='store_true', help='Do not Block N positions')
    parser.add_argument('-write', action='store_true', help='Write full results')
    parser.add_argument('-block_region_bed_url', '--block_region_bed_url', type=str, help='local path of the block region BED file', default=None)
    
    # 验证参数
    parser.add_argument('-validate_only', action='store_true',
                      help='Only validate input without generating variants')
    
    return parser.parse_args()

def validate_input_table(df: pd.DataFrame) -> Tuple[bool, str]:
    """验证输入表格格式是否正确"""
    required_columns = [
        'Index', 'Index_con', 'SV_type', 
        'Original_start', 'Original_end', 'Len_SV',
        'New_start', 'New_end', 'New_len_SV',  # 新增 New_len_SV
        'Balanced Trans Flag'
    ]
    
    # 检查必需列是否存在
    missing_cols = [col for col in required_columns if col not in df.columns]
    if missing_cols:
        return False, f"Missing required columns: {', '.join(missing_cols)}"
    
    # 检查 SV_type 有效性
    valid_sv_types = {
        'Substitution', 'Small_Ins', 'Small_Del',
        'Deletion', 'Insertion', 'Inversion',
        'Translocation', 'Duplication'
    }
    
    invalid_types = set(df['SV_type'].unique()) - valid_sv_types
    if invalid_types:
        return False, f"Invalid SV_type values: {', '.join(invalid_types)}"
    
    # 检查数值列是否为整数类型
    numeric_columns = [
        'Original_start', 'Original_end', 'Len_SV',
        'New_start', 'New_end', 'New_len_SV',
        'Balanced Trans Flag'
    ]
    for col in numeric_columns:
        if not pd.api.types.is_integer_dtype(df[col]):
            try:
                df[col] = df[col].astype('int64')
            except (ValueError, TypeError):
                return False, f"Column '{col}' must contain integer values"
    
    # 检查 Balanced Trans Flag 只能是 0 或 1 或 -1
    invalid_flags = set(df['Balanced Trans Flag'].unique()) - {0, 1, -1}
    if invalid_flags:
        return False, f"Balanced Trans Flag must be 0 or 1, found: {invalid_flags}"
    
    return True, ""

def read_variant_table(file_path: str) -> pd.DataFrame:
    """Read variant table file and standardize SV_type format
    
    Args:
        file_path: Input file path
        
    Returns:
        Standardized DataFrame
        
    Features:
        1. Automatically detects CSV or BED format
        2. Handles different column name cases (SV_type, sv_type, etc.)
        3. Converts shorthand SV types to standard forms
        4. Removes unsupported variant types
        5. Preserves original indices for tracking removed rows
    """
    # Define standard forms and possible shorthands
    SV_TYPE_MAPPING = {
        'Translocation': {'Translocation', 'TRANS', 'TRA'},
        'Inversion': {'Inversion', 'INV'},
        'Duplication': {'Duplication', 'DUP'},
        'Deletion': {'Deletion', 'DEL'},
        'Insertion': {'Insertion', 'INS'},
        'Small_Del': {'Small_Del', 'SMALL_DEL', 'MICRO_DEL'},
        'Small_Ins': {'Small_Ins', 'SMALL_INS', 'MICRO_INS'},
        'Substitution': {'Substitution', 'SUB', 'SNP'}
    }
    
    # Create reverse mapping: shorthand->standard
    shorthand_to_standard = {}
    for standard, shorthands in SV_TYPE_MAPPING.items():
        for shorthand in shorthands:
            shorthand_to_standard[shorthand] = standard
    
    try:
        # First try to detect if it's a BED file (no header, at least 3 columns)
        try:
            # Try reading as BED (no header, tab-separated)
            df = pd.read_csv(file_path, sep='\t', header=None, 
                           names=['chrom', 'start', 'end'], 
                           usecols=[0,1,2])
            
            # If successful, it's a BED file - create minimal variant table
            df['SV_type'] = 'Deletion'  # BED files typically represent deletions
            df['Original_start'] = df['start']
            df['Original_end'] = df['end']
            df['Len_SV'] = df['end'] - df['start']
            
            # Fill other required columns with defaults
            for col in ['New_start', 'New_end', 'New_len_SV', 'Balanced Trans Flag']:
                df[col] = -1
            
        except:
            # If BED reading fails, try as CSV/TSV with header
            try:
                # Try tab separator first (TSV)
                df = pd.read_csv(file_path, sep='\t', header=0, engine='python')
                if len(df.columns) == 1:  # If only one column, try comma
                    df = pd.read_csv(file_path, sep=',', header=0, engine='python')
            except:
                # Final attempt with whitespace
                df = pd.read_csv(file_path, delim_whitespace=True, header=0)
        
        # Standardize column names (case insensitive, space to underscore)
        df.columns = df.columns.str.strip().str.lower().str.replace(' ', '_')
        
        # Check for required SV type column (handle different cases)
        sv_type_col = None
        for col in df.columns:
            if col.lower() in ['sv_type', 'type', 'variant_type']:
                sv_type_col = col
                break
                
        if not sv_type_col:
            raise ValueError("Input file must contain a variant type column (SV_type, Type, etc.)")
        
        # Rename to standard 'SV_type'
        df = df.rename(columns={sv_type_col: 'SV_type'})
        
        # Standardize other common column names
        col_mapping = {
            'index': 'Index',
            'index_con': 'Index_con',
            'original_start': 'Original_start',
            'start': 'Original_start',
            'original_end': 'Original_end',
            'end': 'Original_end',
            'len_sv': 'Len_SV',
            'length': 'Len_SV',
            'new_start': 'New_start',
            'new_end': 'New_end',
            'new_len_sv': 'New_len_SV',
            'balanced_trans_flag': 'Balanced Trans Flag'
        }
        
        for old, new in col_mapping.items():
            if old in df.columns:
                df = df.rename(columns={old: new})
        
        # Record original info
        original_count = len(df)
        original_indices = set(df.index)
        
        # Standardize SV_type: convert to string, strip, uniform case
        df['SV_type'] = df['SV_type'].astype(str).str.strip().str.replace(' ', '_')
        
        # Create conversion column: map shorthands to standard forms
        df['SV_type_standard'] = df['SV_type'].map(shorthand_to_standard)
        
        # Separate valid and invalid rows
        valid_df = df[df['SV_type_standard'].notna()].copy()
        invalid_df = df[df['SV_type_standard'].isna()].copy()
        
        # Apply standard forms
        valid_df['SV_type'] = valid_df['SV_type_standard']
        valid_df.drop(columns=['SV_type_standard'], inplace=True)
        
        # Report filtering
        if not invalid_df.empty:
            invalid_types = invalid_df['SV_type'].unique()
            invalid_indices = invalid_df.index.tolist()
            print(f"Warning: Removed {len(invalid_df)} rows with unsupported SV types.")
            print(f"Unsupported types: {', '.join(invalid_types)}")
            print(f"Affected row indices (0-based): {invalid_indices}")
        
        print(f"Original rows: {original_count}, Kept rows: {len(valid_df)}")
        
        return valid_df
    
    except Exception as e:
        raise ValueError(f"Failed to read input file {file_path}: {str(e)}")


def load_reference_sequence(fasta_path: str, seq_index: int) -> Tuple[str, str]:
    """加载参考序列，支持选择特定索引的序列"""
    try:
        with pysam.FastaFile(fasta_path) as fasta_file:
            if seq_index >= len(fasta_file.references):
                raise ValueError(f"Sequence index {seq_index} out of range (max {len(fasta_file.references)-1})")
            
            seqname = fasta_file.references[seq_index]
            return seqname, fasta_file.fetch(seqname)
    except Exception as e:
        raise RuntimeError(f"Failed to load reference sequence: {str(e)}")

def main():
    args = parse_args()
    os.makedirs(args.save, exist_ok=True)
    # Start timer for entire process
    total_start_time = time.time()
    print("\n=== Starting Variant Simulation Process ===")
    print("Initializing data structures...")
    init_start = time.time()
    
    # 设置全局随机数种子
    random.seed(args.seed)
    np.random.seed(args.seed)
    os.environ['PYTHONHASHSEED'] = str(args.seed)
    
    # 打印种子信息
    print(f"Global Random Seed: {args.seed}")
    
    # 1. 加载参考序列
    seqname, ref_seq = load_reference_sequence(args.ref, args.seq_index)
    print(f"Loaded reference sequence: {seqname} (length: {len(ref_seq)})")
    
    # 2. 加载变异表格
    # 定义数值列，确保它们被读取为 int 类型
    dtype_dict = {
        'Original_start': 'int64',
        'Original_end': 'int64',
        'Len_SV': 'int64',
        'New_start': 'int64',
        'New_end': 'int64',
        'New_len_SV': 'int64',  # 新增这一行
        'Balanced Trans Flag': 'int64'
    }

    # variant_df = pd.read_csv(
    #     args.variant_table,
    #     sep='\t' if args.variant_table.endswith('.tsv') else ',',
    #     dtype=dtype_dict  # 强制指定数据类型
    # )
    variant_df = read_variant_table(args.variant_table)
    print("Raw data after reading:")  # 调试信息
    print(variant_df.head())
    print(variant_df.columns)
    
    is_valid, validation_msg = validate_input_table(variant_df)
    if not is_valid:
        raise ValueError(f"Invalid variant table: {validation_msg}")
        
    # 确保关键列存在
    required_cols = ['index', 'index_con', 'sv_type', 
                    'original_start', 'original_end', 'len_sv']
    for col in required_cols:
        if col not in variant_df.columns:
            variant_df[col] = -1  # 或适当的默认值

    
    print(variant_df.columns)
    # 再次检查数值列是否真的转换为 int
    for col in dtype_dict:
        if col in variant_df.columns:  # 先检查列是否存在
            if not pd.api.types.is_integer_dtype(variant_df[col]):
                variant_df[col] = variant_df[col].astype('int64')  # 强制转换


    
    # 3. 验证表格格式
    is_valid, validation_msg = validate_input_table(variant_df)
    if not is_valid:
        raise ValueError(f"Invalid variant table: {validation_msg}")
    
    print(f"Loaded {len(variant_df)} variants from input table")
    
    # 如果只是验证模式则退出
    if args.validate_only:
        print("Validation successful - exiting")
        sys.exit(0)
        
    # 4. 初始化数据结构
    
    #! whole chr
    chr_id=seqname
    start_base=0
    end_base=len(ref_seq)
    real_con1 = copy.deepcopy(ref_seq).upper()
    chr_length = len(real_con1)
    
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

    
    # 5. 按变异类型分类处理
    sv_groups = variant_df.groupby('SV_type')
    # 获取最长的序列
    # seqname = fasta_file.references[0]
    # BestRefSeq = fasta_file.fetch(seqname)

    del ref_seq
    print('Length of ref:'+str(chr_length))

    #!Block N base positions
    all_positions = set(range(len(real_con1)))
    
    if args.notblockN:
        all_region_sv= list(all_positions)
    else:
        # collection of N positions
        n_positions = [i for i, char in enumerate(real_con1) if char == 'N']
        n_positions_set = set(n_positions)
        #all_positions = set(range(len(real_con1)))
        print('N pos:'+str(len(n_positions)))
        
        all_region_sv=list(all_positions - n_positions_set)
        
        del n_positions
        
    block_region_bed = args.block_region_bed_url
        
    if block_region_bed is None or block_region_bed == 'None':
        print("block_region_bed_url is None")
    else:
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
        
        
    #! initialization
    unblock_region_sv = copy.deepcopy(all_region_sv)
    tem_seq_post = copy.deepcopy(list(real_con1))
    
    del all_region_sv

    SV_loop = 0
    VCF_loop = 0
    # 逐行写入变异数据的逻辑
    loop_index = 0
    CSV_loop = 0
    
    #Whole_INS_con = []#restore Ins_dic_sv for each consensus
    updated_con = []
    #record the position and content of insertions, to update the final consensus
    Ins_dic_sv = {}
    #restore sites of ins for trans and dup in case we sample new ins positions there
    location_insert_trans = []
    location_insert_dup = []
    #restore sites of del for translocation, inversion
    location_del_trans = []
    location_del_inv = []

    #! default input
    ins_selection = ['a', 't', 'c', 'g']
    mis_selection = ['a','t','c','g']

    base_list=["A","T","C","G"]


    # 从variant_df中计算各种变异的数量
    sv_trans = len(variant_df[variant_df['SV_type'] == 'Translocation'])
    sv_inver = len(variant_df[variant_df['SV_type'] == 'Inversion'])
    sv_dup = len(variant_df[variant_df['SV_type'] == 'Duplication'])
    sv_del = len(variant_df[variant_df['SV_type'] == 'Deletion'])
    sv_ins = len(variant_df[variant_df['SV_type'] == 'Insertion'])
    del_snv_number = len(variant_df[variant_df['SV_type'] == 'Small_Del'])
    ins_snv_number = len(variant_df[variant_df['SV_type'] == 'Small_Ins'])
    snp = len(variant_df[variant_df['SV_type'] == 'Substitution'])

    # 计算最大长度
    max_csv_len = sv_trans + sv_dup
    max_sim_len = sv_inver + sv_del + sv_ins + del_snv_number + ins_snv_number + snp + sv_dup + args.snp + args.snv_ins + args.snv_del

    # 初始化 numpy 数组
    SV_table = np.empty((int(max_csv_len+2), 8), dtype=object)  # 减少列数，去掉 'Index' 和 'Index_con'
    VCF_table = np.empty((int(max_csv_len*4), 4), dtype=object)  # 减少列数
    Unified_table = np.empty((int(max_sim_len*2), 10), dtype=object)
    print(f"Debug - Unified_table shape: {Unified_table.shape}")
    print(f"Debug - Current loop_index: {loop_index}")
    
    # --- 后续处理逻辑将在这里添加 ---
    # 将按变异类型分别处理：
    # Process each SV type separately

    init_end = time.time()
    print(f"Initialization complete. Reference length: {chr_length}")
    print(f"Time taken: {format_time(init_end - init_start)}")

    # 4. 处理每种变异类型
    print("\n[4/6] Processing variant types...")
    process_start = time.time()
    
    # ==================== TRANSLOCATION PROCESSING ====================
    trans_start = time.time()
    # ==================== TRANSLOCATION PROCESSING ====================
    trans_df = variant_df[variant_df['SV_type'] == 'Translocation'].copy()

    for _, row in trans_df.iterrows():
        original_start = int(row['Original_start'])
        original_end = int(row['Original_end'])
        len_sv = int(row['Len_SV'])
        new_start = int(row['New_start'])
        new_end = int(row['New_end'])
        is_balanced = int(row['Balanced Trans Flag'])
        
        # Check if regions are available
        if any(pos not in unblock_region_sv for pos in range(original_start, original_end+1)):
            print(f"Warning: Skipping unavailable original range in translocation: {original_start}-{original_end}")
            continue
        
        if is_balanced and any(pos not in unblock_region_sv for pos in range(new_start, new_end+1)):
            print(f"Warning: Skipping unavailable new range in balanced translocation: {new_start}-{new_end}")
            continue
        
        # For translocation, we always delete the original region
        deleted_sequence = ''.join(tem_seq_post[original_start:original_end+1])
        #! move to later
        # tem_seq_post[original_start:original_end+1] = ['-'] * len_sv
        location_del_trans.append(original_start)
        
        # Update blocked regions
        unblock_region_sv = list(set(unblock_region_sv) - set(range(original_start-1, original_end+2)))
        
        if is_balanced:
            # Balanced translocation - swap two regions
            # Delete the new region
            swapped_sequence = ''.join(tem_seq_post[new_start:new_end+1])
            # tem_seq_post[new_start:new_end+1] = ['-'] * (new_end - new_start + 1)
            location_del_trans.append(new_start)
            
            # Insert the original sequence at new location
            ins_loc1 = new_start - 1
            Ins_dic_sv[ins_loc1] = deleted_sequence
            location_insert_trans.append(ins_loc1)
            
            # Insert the swapped sequence at original location
            ins_loc2 = original_start - 1
            Ins_dic_sv[ins_loc2] = swapped_sequence
            location_insert_trans.append(ins_loc2)
            
            # Update blocked regions for new location
            unblock_region_sv = list(set(unblock_region_sv) - set(range(new_start-1, new_end+2)))
            
            # Write to SV table
            SV_table[SV_loop] = [
                'Translocation',
                original_start,
                original_end,
                len_sv,
                new_start,
                new_end,
                new_end - new_start + 1,
                1  # Balanced flag
            ]
            SV_loop += 1
            
            # Write DEL records to VCF table (both deletions)
            VCF_table[VCF_loop] = [
                str(original_start),
                deleted_sequence + tem_seq_post[original_end+1],
                tem_seq_post[original_end+1],
                f"SVTYPE=DEL;SVLEN={len_sv};CSV_TYPE=balancedTrans;CSV_INDEX={CSV_loop}"
            ]
            VCF_loop += 1
            
            VCF_table[VCF_loop] = [
                str(new_start),
                swapped_sequence + tem_seq_post[new_end+1],
                tem_seq_post[new_end+1],
                f"SVTYPE=DEL;SVLEN={new_end-new_start+1};CSV_TYPE=balancedTrans;CSV_INDEX={CSV_loop}"
            ]
            VCF_loop += 1
            
            # Write INS records to VCF table (both insertions)
            VCF_table[VCF_loop] = [
                str(ins_loc1),
                tem_seq_post[ins_loc1],
                tem_seq_post[ins_loc1] + deleted_sequence,
                f"SVTYPE=INS;SVLEN={len_sv};CSV_TYPE=balancedTrans;CSV_INDEX={CSV_loop}"
            ]
            VCF_loop += 1
            
            VCF_table[VCF_loop] = [
                str(ins_loc2),
                tem_seq_post[ins_loc2],
                tem_seq_post[ins_loc2] + swapped_sequence,
                f"SVTYPE=INS;SVLEN={new_end-new_start+1};CSV_TYPE=balancedTrans;CSV_INDEX={CSV_loop}"
            ]
            VCF_loop += 1
            
            tem_seq_post[original_start:original_end+1] = ['-'] * len_sv
            tem_seq_post[new_start:new_end+1] = ['-'] * (new_end - new_start + 1)
            
        else:
            #! Unbalanced translocation - just insert at new location
            ins_loc = new_start  # 关键修改：不再减1
            Ins_dic_sv[ins_loc] = deleted_sequence
            location_insert_trans.append(ins_loc)
            
            # Update blocked regions around insertion point
            unblock_region_sv = list(set(unblock_region_sv) - set(range(ins_loc, ins_loc+2)))
            
            # Write to SV table
            SV_table[SV_loop] = [
                'Translocation',
                original_start,
                original_end,
                len_sv,
                new_start,
                new_start,  # For unbalanced, insertion is at single point
                len_sv,
                0  # Not balanced
            ]
            SV_loop += 1
            
            # Write DEL record to VCF table
            VCF_table[VCF_loop] = [
                str(original_start),
                deleted_sequence + tem_seq_post[original_end+1],
                tem_seq_post[original_end+1],
                f"SVTYPE=DEL;SVLEN={len_sv};CSV_TYPE=unbalancedTrans;CSV_INDEX={CSV_loop}"
            ]
            VCF_loop += 1
            
            # Write INS record to VCF table
            VCF_table[VCF_loop] = [
                str(ins_loc),
                tem_seq_post[ins_loc],
                tem_seq_post[ins_loc] + deleted_sequence,
                f"SVTYPE=INS;SVLEN={len_sv};CSV_TYPE=unbalancedTrans;CSV_INDEX={CSV_loop}"
            ]
            VCF_loop += 1
            
            tem_seq_post[original_start:original_end+1] = ['-'] * len_sv
        
        CSV_loop += 1
        
    trans_end = time.time()
    print(f"Processed {len(trans_df)} translocations - {format_time(trans_end - trans_start)}")

    # ==================== INVERSION PROCESSING ====================
    inv_start = time.time()
    # ==================== INVERSION PROCESSING ====================
    inv_df = variant_df[variant_df['SV_type'] == 'Inversion'].copy()

    for _, row in inv_df.iterrows():
        # 从表格中获取反转区域信息
        inv_start = int(row['Original_start'])
        inv_end = int(row['Original_end'])
        inv_len = inv_end - inv_start + 1
        
        # 检查区域是否可用（避免与已有变异冲突）
        if any(pos not in unblock_region_sv for pos in range(inv_start, inv_end+1)):
            print(f"Warning: Skipping unavailable inversion region: {inv_start}-{inv_end}")
            continue
        
        # 获取原始序列并生成反向互补序列
        original_sequence = ''.join(tem_seq_post[inv_start:inv_end+1])
        inverted_sequence = DNA_complement(original_sequence[::-1])  # 反向互补
        
        # 更新阻塞区域（避免后续变异重叠）
        # 更新全局阻塞区域（最大范围剔除：原区域±1bp）
        blocked_range = range(inv_start - 1, inv_end + 2)  # 扩展1bp边界
        unblock_region_sv = list(set(unblock_region_sv) - set(blocked_range))
        
        # 记录到 Unified_table
        Unified_table[loop_index] = [
            'Inversion',        # 变异类型
            inv_start,          # 原始起始位置
            inv_end,            # 原始结束位置
            inv_len,            # 变异长度
            inv_start,          # 新起始位置（反转后相同）
            inv_end,            # 新结束位置（反转后相同）
            inv_len,            # 新长度（不变）
            -1,                 # 特殊标记（-1表示反转）
            original_sequence,  # 原始序列
            inverted_sequence   # 反向互补序列
        ]
        loop_index += 1
        
        # 实际修改参考序列
        tem_seq_post[inv_start:inv_end+1] = list(inverted_sequence)
        
    
    inv_end = time.time()
    print(f"Processed {len(inv_df)} inversions - {format_time(inv_end - inv_start)}")

    # ==================== DUPLICATION PROCESSING ====================
    dup_start = time.time()   
    # ==================== DUPLICATION PROCESSING ====================
    dup_df = variant_df[variant_df['SV_type'] == 'Duplication'].copy()

    for _, row in dup_df.iterrows():
        # 从表格读取坐标
        dup_source_start = int(row['Original_start'])  # 被复制区域起始
        dup_source_end = int(row['Original_end'])      # 被复制区域结束
        dup_target_pos = int(row['New_start'])         # 插入位置（前一个碱基的位置）
        dup_len = dup_source_end - dup_source_start + 1
        
        # 检查被复制区域是否全部可用
        if not all(pos in unblock_region_sv for pos in range(dup_source_start, dup_source_end + 1)):
            print(f"Warning: Duplication source region {dup_source_start}-{dup_source_end} overlaps. Skipping.")
            continue
        
        # 检查插入位置是否可用
        if dup_target_pos not in unblock_region_sv:
            print(f"Warning: Duplication target position {dup_target_pos} is blocked. Skipping.")
            continue

        # 获取被复制的序列
        duplicated_seq = ''.join(tem_seq_post[dup_source_start:dup_source_end + 1])
        
        unblock_region_sv = list(set(unblock_region_sv) - set(range(dup_source_start - 1, dup_source_end + 2)))  # ±1bp
        unblock_region_sv = list(set(unblock_region_sv) - set(range(dup_target_pos, dup_target_pos + 2)))  # 插入点及+1

        # 判断是邻近复制还是远程复制
        if abs(dup_target_pos - dup_source_end) <= 1:  # 邻近复制
            # 更新阻塞区域
            
            
            # 记录插入点和序列
            location_insert_dup.append(dup_target_pos)
            Ins_dic_sv[dup_target_pos] = duplicated_seq
            ins_len_dup = len(duplicated_seq)
            
            # 写入SV_table和VCF_table
            SV_table[SV_loop] = ['Duplication', dup_source_start, dup_source_end, dup_len, dup_target_pos, dup_target_pos, ins_len_dup, -1]
            SV_loop += 1
            
            VCF_table[VCF_loop] = [
                dup_target_pos,
                tem_seq_post[dup_target_pos],
                tem_seq_post[dup_target_pos] + duplicated_seq,
                f'SVTYPE=DUP;SVLEN={dup_len}'
            ]
            VCF_loop += 1
            
        else:  # 远程复制
            remain_index22 = dup_target_pos
            location_insert_dup.append(remain_index22)
            
            # 更新阻塞区域
            # unblock_region_sv = list(set(unblock_region_sv) - set(range(dup_source_start - 1, dup_source_end + 2)))  # ±1bp
            # unblock_region_sv = list(set(unblock_region_sv) - set(range(remain_index22, remain_index22 + 2)))  # 插入点及+1
            
            # 记录插入点和序列
            Ins_dic_sv[remain_index22] = duplicated_seq
            ins_len_dup = len(duplicated_seq)
            
            # 写入SV_table和VCF_table
            SV_table[SV_loop] = ['Duplication', dup_source_start, dup_source_end, dup_len, remain_index22, remain_index22, ins_len_dup, -1]
            SV_loop += 1
            
            VCF_table[VCF_loop] = [
                remain_index22,
                tem_seq_post[remain_index22],
                tem_seq_post[remain_index22] + duplicated_seq,
                f'SVTYPE=DUP;SVLEN={dup_len};CSV_TYPE=DisDup;CSV_INDEX={CSV_loop}'
            ]
            VCF_loop += 1
            CSV_loop += 1

    
    dup_end = time.time()
    print(f"Processed {len(dup_df)} duplications - {format_time(dup_end - dup_start)}")

    del_start = time.time()
    # ==================== LONG DELETION PROCESSING ====================
    del_df = variant_df[variant_df['SV_type'] == 'Deletion'].copy()
    print(del_df)
    # 应用重叠检查
    original_count = len(del_df)
    # del_df = check_and_filter_overlapping_deletions(del_df)
    # if len(del_df) < original_count:
    #     print(f"Filtered {original_count - len(del_df)} overlapping deletions. Kept {len(del_df)} deletions.")
    print(f"Unblock regions count: {len(unblock_region_sv)}")
    print(f"First 10 unblock regions: {unblock_region_sv[:10]}")
    # 继续原有处理逻辑
    for _, row in del_df.iterrows():
        # 从表格中获取删除区域坐标
        del_start = int(row['Original_start'])
        del_end = int(row['Original_end'])
        del_len = del_end - del_start + 1
        
        # 检查删除区域是否全部可用 (使用unblock_region_sv)
        if not all(pos in unblock_region_sv for pos in range(del_start, del_end + 1)):
            print(f"Warning: Deletion region {del_start}-{del_end} overlaps with blocked region. Skipping.")
            continue
        
        # 获取删除前后的上下文序列（用于VCF格式）
        deleted_sequence = ''.join(tem_seq_post[del_start:del_end + 1])
        remaining_base = tem_seq_post[del_end + 1] if del_end + 1 < len(tem_seq_post) else 'N'
        
        # 更新阻塞区域（关键修改：仅使用unblock_region_sv）
        blocked_range = range(del_start - 1, del_end + 2)  # 扩展1bp边界
        unblock_region_sv = list(set(unblock_region_sv) - set(blocked_range))
        
        # 写入Unified_table
        Unified_table[loop_index] = [
            'Deletion',         # 变异类型
            del_start,          # 删除起始位置
            del_end,            # 删除结束位置
            del_len,            # 删除长度
            -1, -1, -1, -1,     # 无插入相关字段
            deleted_sequence + remaining_base,  # REF = 被删序列+右侧第一个碱基
            remaining_base      # ALT = 右侧第一个碱基
        ]
        loop_index += 1
        
        # 标记删除区域（用'-'填充）
        tem_seq_post[del_start:del_end + 1] = ['-'] * del_len
    
    del_end = time.time()
    print(f"Processed {len(del_df)} deletions (filtered {original_count - len(del_df)} overlaps) - {format_time(del_end - del_start)}")

    # ==================== INSERTION PROCESSING ====================
    ins_start = time.time()
    # ==================== INSERTION PROCESSING ====================
    ins_df = variant_df[variant_df['SV_type'] == 'Insertion'].copy()

    for _, row in ins_df.iterrows():
        # 从表格中获取插入位置和长度
        ins_pos = int(row['Original_start'])  # 插入位置（VCF规范：插入点前一个碱基的位置）
        ins_len = int(row['Len_SV'])          # 插入序列长度
        
        # 检查插入位置是否可用（需检查插入点及其右邻）
        if ins_pos not in unblock_region_sv or (ins_pos + 1) not in unblock_region_sv:
            print(f"Warning: Insertion position {ins_pos} is blocked. Skipping.")
            continue
        
        # 生成随机插入序列（或从表格读取预设序列）
        if 'Inserted_Sequence' in row and pd.notna(row['Inserted_Sequence']):
            inserted_seq = row['Inserted_Sequence'].upper()
        else:
            inserted_seq = ''.join(np.random.choice(['A', 'T', 'C', 'G'], size=ins_len))
        
        # 更新阻塞区域（关键修改：仅使用unblock_region_sv）
        # 阻塞插入点及其右邻（防止其他变异重叠）
        unblock_region_sv = list(set(unblock_region_sv) - {ins_pos, ins_pos + 1})
        
        # 记录插入序列（统一通过Ins_dic_sv管理）
        Ins_dic_sv[ins_pos] = inserted_seq
        
        # 写入Unified_table
        Unified_table[loop_index] = [
            'Insertion',            # 变异类型
            ins_pos,                # 插入位置（VCF规范：前一个碱基）
            ins_pos,                # 结束位置（单点插入）
            ins_len,                # 插入序列长度
            -1, -1, -1, -1,         # 无复制/删除相关字段
            tem_seq_post[ins_pos],   # REF：插入点原碱基
            tem_seq_post[ins_pos] + inserted_seq  # ALT：REF+插入序列
        ]
        loop_index += 1
        
    ins_end = time.time()
    print(f"Processed {len(ins_df)} insertions - {format_time(ins_end - ins_start)}")

    # ==================== SMALL VARIANT PROCESSING ====================
    small_start = time.time()
    small_del_df = variant_df[variant_df['SV_type'] == 'Small_Del'].copy()

    if not small_del_df.empty:
        # 处理用户提供的small deletion
        for _, row in small_del_df.iterrows():
            # 从表格中获取删除位置和长度
            del_start = int(row['Original_start'])
            del_len = int(row['Len_SV'])
            del_end = del_start + del_len - 1
            
            # 检查删除区域是否全部可用
            if not all(pos in unblock_region_sv for pos in range(del_start, del_end + 1)):
                print(f"Warning: Small deletion region {del_start}-{del_end} overlaps. Skipping.")
                continue
            
            # 获取删除前后的上下文序列
            deleted_sequence = ''.join(tem_seq_post[del_start:del_end + 1])
            remaining_base = tem_seq_post[del_end + 1] if del_end + 1 < len(tem_seq_post) else 'N'
            
            # 更新阻塞区域
            blocked_range = range(del_start - 1, del_end + 2)
            unblock_region_sv = list(set(unblock_region_sv) - set(blocked_range))
            
            # 写入Unified_table
            Unified_table[loop_index] = [
                'Small_Del', del_start, del_end, del_len,
                -1, -1, -1, -1,
                deleted_sequence + remaining_base,
                remaining_base
            ]
            loop_index += 1
            
            # 标记删除区域
            tem_seq_post[del_start:del_end + 1] = ['-'] * del_len
    else:
        # 如果没有用户提供的small deletion，则随机生成
        if not unblock_region_sv:
            print("Warning: no available positions for small deletions.")
        else:
            # 计算参数
            len_unblock_region = len(unblock_region_sv)
            max_small_deletion_length = 5
            min_gap = 1
            DIVISOR = max_small_deletion_length + min_gap
            max_micro_dels = (len_unblock_region + min_gap) // DIVISOR
            
            # 确定要生成的deletion数量
            if args.snv_del is None:
                pai_pro_tem = [diff_ins_prob_del_real, diff_ins_prob_mis_real, diff_ins_prob_correct_real]
                pai_pro_tem_ = list(np.array(pai_pro_tem) / sum(pai_pro_tem))
                l_s_vec_ini = np.random.multinomial(n=len_unblock_region, pvals=pai_pro_tem_)
                del_snv_number = int(l_s_vec_ini[0])
                print(f'Info: Generating {del_snv_number} small deletions based on default probability.')
            elif args.snv_del < 1:
                del_snv_number = int((len_unblock_region // DIVISOR) * args.snv_del)
                print(f'Info: Generating {del_snv_number} small deletions ({args.snv_del*100}% of available capacity).')
            else:
                input_del = int(args.snv_del)
                del_snv_number = min(input_del, max_micro_dels)
                
                if input_del > max_micro_dels:
                    print(f"Info: Requested {input_del} small deletions but only {max_micro_dels} can be accommodated.")
                    print(f"Info: Adjusted small deletion count from {input_del} to {del_snv_number} based on available space.")
                else:
                    print(f"Info: Generating {del_snv_number} small deletions as requested.")
            
            print(f'Random small dels: {del_snv_number}')
            
            # 生成随机deletion
            len_seg_refine = max(unblock_region_sv)
            for m_del in range(del_snv_number):
                # 检查是否还有可用位置
                if len(unblock_region_sv) < 2:  # 至少需要2个位置(删除位置+右邻)
                    print("Warning: no more available positions for small deletions.")
                    break
                    
                # 随机选择删除起始位置
                r_s = sample(unblock_region_sv, 1)[0]
                
                # 随机选择删除长度
                l_s_vec = np.random.multinomial(n=1, pvals=[6/8,1/8,1/16,1/32,1/32])
                l_s = list(l_s_vec).index(1) + 1
                
                # 检查删除区域是否有效
                circular_count_del = 0
                circular_count_del_break = 0
                while (r_s + l_s > len_seg_refine) or \
                    (r_s + l_s - 1 not in unblock_region_sv) or \
                    (unblock_region_sv.index(r_s + l_s - 1) - unblock_region_sv.index(r_s) < l_s - 1):
                    r_s = sample(unblock_region_sv, 1)[0]
                    l_s_vec = np.random.multinomial(n=1, pvals=[6/8,1/8,1/16,1/32,1/32])
                    l_s = list(l_s_vec).index(1) + 1
                    
                    circular_count_del += 1
                    if circular_count_del > args.times:
                        circular_count_del_break = 1
                        print(f"Warning: No.{m_del} small deletion sampling exceeds {args.times} times.")
                        break
                
                if circular_count_del_break:
                    continue
                    
                # 计算删除结束位置
                del_start = r_s
                del_end = r_s + l_s - 1
                del_len = l_s
                
                # 获取序列信息
                deleted_sequence = ''.join(tem_seq_post[del_start:del_end + 1])
                remaining_base = tem_seq_post[del_end + 1] if del_end + 1 < len(tem_seq_post) else 'N'
                
                # 更新阻塞区域
                blocked_range = range(del_start - 1, del_end + 2)
                unblock_region_sv = list(set(unblock_region_sv) - set(blocked_range))
                
                # 写入Unified_table
                Unified_table[loop_index] = [
                    'Small_Del', del_start, del_end, del_len,
                    -1, -1, -1, -1,
                    deleted_sequence + remaining_base,
                    remaining_base
                ]
                loop_index += 1
                
                # 标记删除区域
                tem_seq_post[del_start:del_end + 1] = ['-'] * del_len

    # ==================== SMALL INSERTION PROCESSING ====================
    small_ins_df = variant_df[variant_df['SV_type'] == 'Small_Ins'].copy()

    if not small_ins_df.empty:
        # 处理用户提供的small insertion
        for _, row in small_ins_df.iterrows():
            # 从表格中获取插入位置和长度
            ins_pos = int(row['Original_start'])
            ins_len = int(row['Len_SV'])
            
            # 检查插入位置是否可用
            if ins_pos not in unblock_region_sv or (ins_pos + 1) not in unblock_region_sv:
                print(f"Warning: Small insertion position {ins_pos} is blocked. Skipping.")
                continue
            
            # 生成随机插入序列（或从表格读取预设序列）
            if 'Inserted_Sequence' in row and pd.notna(row['Inserted_Sequence']):
                inserted_seq = row['Inserted_Sequence'].upper()
            else:
                inserted_seq = ''.join(np.random.choice(['A', 'T', 'C', 'G'], size=ins_len))
            
            # 更新阻塞区域
            unblock_region_sv = list(set(unblock_region_sv) - {ins_pos, ins_pos + 1})
            
            # 记录插入序列
            Ins_dic_sv[ins_pos] = inserted_seq
            
            # 写入Unified_table
            Unified_table[loop_index] = [
                'Small_Ins',            # 变异类型
                ins_pos,                # 插入位置
                ins_pos,                # 结束位置
                ins_len,                # 插入序列长度
                -1, -1, -1, -1,        # 无相关字段
                tem_seq_post[ins_pos],  # REF
                tem_seq_post[ins_pos] + inserted_seq  # ALT
            ]
            loop_index += 1
    else:
        # # 如果没有用户提供的small insertion，则随机生成
        # if not unblock_region_sv:
        #     print("Warning: no available positions for small insertions.")
        # else:
        #     # 计算需要生成的insertion数量
        #     len_unblock_region = len(unblock_region_sv)
            
        #     if args.snv_ins is None:
        #         # 如果没有指定数量，使用默认概率计算
        #         ins_snv_number = (np.random.binomial(len_unblock_region, diff_ins_prob_ins_real, 1))[0]
        #     elif 0 <= args.snv_ins < 1:
        #         # 如果指定的是比例
        #         ins_snv_number = int(len_unblock_region * args.snv_ins)
        #     else:
        #         # 如果指定的是绝对数量
        #         ins_snv_number = min(int(args.snv_ins), len_unblock_region)
            
        #     print(f'Random small ins: {ins_snv_number}')
            
        # 如果没有用户提供的small insertion，则随机生成
        if not unblock_region_sv:
            print("Warning: no available positions for small insertions.")
        else:
            # 计算需要生成的insertion数量
            len_unblock_region = len(unblock_region_sv)
            
            if args.snv_ins is None:
                # 如果没有指定数量，使用默认概率计算
                ins_snv_number = (np.random.binomial(len_unblock_region, diff_ins_prob_ins_real, 1))[0]
                print(f'Info: Generating {ins_snv_number} small insertions based on default probability.')
            elif 0 <= args.snv_ins < 1:
                # 如果指定的是比例
                ins_snv_number = int(len_unblock_region * args.snv_ins)
                print(f'Info: Generating {ins_snv_number} small insertions ({args.snv_ins*100}% of available positions).')
            else:
                # 如果指定的是绝对数量
                original_request = int(args.snv_ins)
                ins_snv_number = min(original_request, len_unblock_region)
                
                if original_request > len_unblock_region:
                    print(f"Info: Requested {original_request} small insertions but only {len_unblock_region} positions available.")
                    print(f"Info: Adjusted small insertion count from {original_request} to {ins_snv_number} based on available positions.")
                else:
                    print(f"Info: Generating {ins_snv_number} small insertions as requested.")
            
            print(f'Random small ins: {ins_snv_number}')
            
            # 生成随机insertion
            for _ in range(ins_snv_number):
                # 检查是否还有可用位置
                if not unblock_region_sv:
                    print("Warning: no more available positions for small insertions.")
                    break
                    
                # 随机选择插入位置
                ins_pos = sample(unblock_region_sv, 1)[0]
                
                # 随机选择插入长度（按照您原来的概率分布）
                l_i_vec = np.random.multinomial(n=1, pvals=[6/8,1/8,1/16,1/32,1/32])
                ins_len = list(l_i_vec).index(1) + 1
                
                # 生成随机插入序列
                inserted_seq = ''.join(np.random.choice(['A', 'T', 'C', 'G'], size=ins_len))
                
                # 检查插入位置是否可用（包括右邻位置）
                if ins_pos not in unblock_region_sv or (ins_pos + 1) not in unblock_region_sv:
                    continue
                
                # 更新阻塞区域
                unblock_region_sv = list(set(unblock_region_sv) - {ins_pos, ins_pos + 1})
                
                # 记录插入序列
                Ins_dic_sv[ins_pos] = inserted_seq
                
                # 写入Unified_table
                Unified_table[loop_index] = [
                    'Small_Ins',            # 变异类型
                    ins_pos,                # 插入位置
                    ins_pos,                # 结束位置
                    ins_len,                # 插入序列长度
                    -1, -1, -1, -1,        # 无相关字段
                    tem_seq_post[ins_pos],  # REF
                    tem_seq_post[ins_pos] + inserted_seq  # ALT
                ]
                loop_index += 1

    # ==================== SUBSTITUTION (SNP) PROCESSING ====================
    snp_df = variant_df[variant_df['SV_type'] == 'Substitution'].copy()
    snp_df = snp_df.drop_duplicates(subset=['Original_start'])  # 去重
    
    # print(snp_df['Original_start'])

    # 确保有可用的位置
    if not unblock_region_sv:
        print("Warning: no available positions for substitutions.")
    else:
        # 处理每个SNP
        for _, row in snp_df.iterrows():
            pos = int(row['Original_start'])  # 假设位置信息在Original_start列
            original_base = tem_seq_post[pos].upper()
            
            # 检查位置是否可用
            if pos not in unblock_region_sv:
                print(f"Warning: position {pos} not available for substitution")
                continue
            
            # 检查是否是有效碱基
            if original_base == 'N':
                print(f"Error: SNP in Gap region at position {pos}")
                continue
            elif original_base not in base_list:
                print(f"Error: Invalid base at position {pos}")
                continue
            
            # 从变异表格中获取替代碱基，如果没有则随机选择
            if 'Substituted_Base' in row and pd.notna(row['Substituted_Base']):
                substituted_base = row['Substituted_Base'].upper()
                # 验证替代碱基是否有效
                if substituted_base not in base_list:
                    print(f"Error: Invalid substituted base {substituted_base} at position {pos}")
                    continue
            else:
                # 随机选择替代碱基
                ref_id = base_list.index(original_base)
                prob_dist = substitution_matrix[ref_id]
                column_index = np.random.choice(4, p=prob_dist)
                substituted_base = mis_selection[column_index]
            
            # 更新序列
            tem_seq_post[pos] = substituted_base
            
            # 从可用位置中移除
            if pos in unblock_region_sv:
                unblock_region_sv.remove(pos)
            
            # 写入Unified表格
            Unified_table[loop_index] = [
                'Substitution',
                pos,
                pos,
                1,  # Length
                -1, -1, -1, -1,  # Placeholders
                original_base,
                substituted_base.upper()
            ]
            loop_index += 1
            
    
    small_end = time.time()
    print(f"Processed small variants - {format_time(small_end - small_start)}")

    process_end = time.time()
    print(f"\nTotal variant processing time: {format_time(process_end - process_start)}")

    # 5. 生成最终输出
    print(" Generating final outputs...")
    output_start = time.time()
    
    
    tem_seq_post_up = copy.deepcopy(tem_seq_post)
    for idx in sorted(Ins_dic_sv, reverse=True):
        #idx = 4981
        tem_seq_post_up.insert(
            idx+1, Ins_dic_sv[idx])


    tem_seq_post_up_string = ''.join(tem_seq_post_up)
    tem_seq_post_up_string= tem_seq_post_up_string.replace('-','')
    updated_con.append(copy.deepcopy(tem_seq_post_up_string))
    print('Length of the simulated sequence: '+str(len(tem_seq_post_up_string)))
    # 找出所有元素都是 None 或空字符串的行
    mask = np.all((SV_table == None) | (SV_table == ''), axis=1)

    # 删除这些行
    SV_table = SV_table[~mask]

    mask_uni = np.all((pd.isna(Unified_table)) | (Unified_table == ''), axis=1)
    Unified_table = Unified_table[~mask_uni]
    # print(Unified_table)
    # 现在你可以将 SV_table 转换为 DataFrame
    # SV_table_merged = pd.DataFrame(SV_table, columns=['Index','Index_con','SV_type','Original_start',\
                                        # 'Original_end','Len_SV','New_start','New_end','New_len_SV','Balanced Trans Flag','relative start1','relative end1','relative start2','relative end2'])
    
    # 拆分 Unified_table
    # 提取前 8 列作为 SV_sub
    SV_sub = Unified_table[:, :8]  # 前 8 列
    # print(SV_sub)

    # 合并 SV_table 和 SV_sub
    # 假设 SV_table 已经初始化并填充了数据
    SV_table_merged = np.vstack((SV_table, SV_sub))  # 垂直合并

    # 将 SV_table 转换为 DataFrame
    SV_table_merged = pd.DataFrame(SV_table_merged, columns=['SV_type', 'Original_start', 'Original_end', 
                                                'Len_SV', 'New_start', 'New_end', 
                                                'New_len_SV', 'Balanced Trans Flag'])
    # print(len(SV_table))
    # print(len(SV_table_merged))
    # print(len(SV_sub))
    
    # 添加 'Index_con' 列并填充 ll_c 的值
    ll_c = args.seq_index
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
    sv_types = Unified_table[:, 0]  # SV_type 列
    sv_lens = Unified_table[:, 3]   # Len_SV 列

    # 生成 INFO 列
    info_column = generate_info_column(sv_types, sv_lens)
    # 提取 VCF_sub
    # VCF_sub 需要包含 'POS' 和 'INFO' 列，并根据规则生成
    VCF_sub = np.empty((Unified_table.shape[0], 4), dtype=object)  # 初始化 VCF_sub
    VCF_sub[:, 0] = Unified_table[:, 1]  # 'POS' 列 = 'Original_start'
    VCF_sub[:, 1] = Unified_table[:, 8]  # 'REF' 列
    VCF_sub[:, 2] = Unified_table[:, 9]  # 'ALT' 列
    VCF_sub[:, 3] = info_column  # 'INFO' 列
    # # 根据 SV_type 和 Len_SV 生成 INFO 列
    # for i in range(Unified_table.shape[0]):
    #     sv_type = Unified_table[i, 0]  # SV_type
    #     sv_len = Unified_table[i, 3]   # Len_SV
    #     VCF_sub[i, 3] = get_info_from_sv_type(sv_type, sv_len)  # INFO 列


    # 对 VCF_table 重复相同的步骤
    mask = np.all((pd.isna(VCF_table)) | (VCF_table == ''), axis=1)
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

    
    tem_ins_dic = Ins_dic_sv

    #! visualize simulated pattern
    counts = {}
    total = 0
    
    for sv_type in sv_type_order:
        count = SV_table_merged[SV_table_merged['SV_type'] == sv_type].shape[0]
        counts[sv_type] = count
        total += count
        print(f'Simulated {sv_type.lower()}s:' + str(count))

    # 输出总数
    print('Total simulated Variations:' + str(total))
    
    #! 0409

    
    # 按照 'pos' 列排序
    #VCF_table_merged = copy.deepcopy(VCF_table)
    VCF_table_merged.sort_values(by='POS', inplace=True)

    # 重置索引并将旧索引添加为 'Index' 列
    VCF_table_merged.reset_index(inplace=True, drop=True)

    # 更新 'ID' 列
    VCF_table_merged['ID'] = 'rs' + VCF_table_merged.index.astype(str)


    write_template_fasta_con(args, seqname, updated_con[0])
    write_vcf(args, VCF_table_merged, seqname, start_base, end_base)
    
    #! final table
    # 强制转换数值列
    numeric_columns = [
        'Original_start', 'Original_end', 'Len_SV',
        'New_start', 'New_end', 'New_len_SV',
        'Balanced Trans Flag',
        'relative start1', 'relative end1',
        'relative start2', 'relative end2'
    ]

    for col in numeric_columns:
        try:
            SV_table_merged[col] = pd.to_numeric(SV_table_merged[col], errors='raise').astype('int64')
        except Exception as e:
            raise ValueError(f"无法将列 '{col}' 转换为整数: {str(e)}")
        
    # 预期列顺序
    EXPECTED_COLUMNS = [
        'Index', 'Index_con', 'SV_type', 'Original_start', 'Original_end',
        'Len_SV', 'New_start', 'New_end', 'New_len_SV', 'Balanced Trans Flag',
        'relative start1', 'relative end1', 'relative start2', 'relative end2'
    ]

    # 验证列顺序和数据类型
    if not SV_table_merged.columns.tolist() == EXPECTED_COLUMNS:
        raise ValueError(
            f"列顺序不匹配!\n"
            f"当前顺序: {SV_table_merged.columns.tolist()}\n"
            f"预期顺序: {EXPECTED_COLUMNS}"
        )

    # 验证关键列数据类型
    if not pd.api.types.is_numeric_dtype(SV_table_merged['Original_start']):
        raise TypeError("Original_start 列必须为数值类型")
    
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
    
    output_end = time.time()
    print(f"Output generation complete - {format_time(output_end - output_start)}")
    total_end_time = time.time()
    print(f"\n=== Simulation Completed Successfully ===")
    print(f"Total execution time: {format_time(total_end_time - total_start_time)}")
    # print(f"Results saved to: {args.save}")
            
    
    
    
    
    
    
    
    # 6. 最终输出
    
    print(f"Successfully generated modified sequence at {args.save}")
        
        
        
        
    

if __name__ == "__main__":
    main()