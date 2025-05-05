#!/usr/bin/env python3
import argparse
import pandas as pd
import gzip
import re
import os
import time
from typing import Optional, Tuple, List, Dict, Union

def parse_args():
    parser = argparse.ArgumentParser(description="Filter SVs from VCF and convert to CSV format")
    parser.add_argument("-vcf", required=True, help="Input VCF/VCF.GZ file path")
    parser.add_argument("-chr", required=True, help="Target chromosome (e.g., 21)")
    parser.add_argument("-select", help="Selection criteria (e.g., 'AF>0.001', 'SVLEN>=100')")
    parser.add_argument("-min_len", type=int, default=50, help="Minimum SV length (bp)")
    parser.add_argument("-sv_types", nargs="+", default=["DEL", "INS", "DUP", "INV"], 
                      help="SV types to include")
    parser.add_argument("-save", type=str, help="Output directory", 
                      default=os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "save") + "/")
    return parser.parse_args()

def format_time(seconds: float) -> str:
    """Convert seconds to HH:MM:SS format"""
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{int(hours):02d}:{int(minutes):02d}:{int(seconds):02d}"

def parse_selection_condition(select_str: Optional[str]) -> Optional[Tuple[str, str, Union[float, int]]]:
    """Parse selection condition into (field, operator, value)"""
    if not select_str:
        return None
    
    match = re.match(r"^([A-Za-z_]+)([<>=!]+)([\d\.]+)$", select_str)
    if not match:
        raise ValueError(f"Invalid selection condition format: {select_str}. Expected like 'AF>0.001'")
    
    field, operator, value = match.groups()
    return field, operator, float(value) if '.' in value else int(value)

def parse_info_field(info_str: str) -> Dict[str, Union[str, bool]]:
    """Parse VCF INFO field into dictionary"""
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, val = item.split('=', 1)
            info_dict[key] = val
        else:
            info_dict[item] = True
    return info_dict

def passes_filters(
    chrom: str,
    svtype: str,
    svlen: int,
    info_dict: Dict[str, Union[str, bool]],
    chrom_filter: str,
    min_len: int,
    sv_types: List[str],
    select_cond: Optional[Tuple[str, str, Union[float, int]]]
) -> bool:
    """Check if variant passes all filters"""
    # Chromosome filter
    if chrom != chrom_filter:
        return False
    
    # SV type filter
    if svtype not in sv_types:
        return False
    
    # SV length filter
    if abs(svlen) < min_len:
        return False
    
    # Selection condition filter
    if select_cond:
        field, op, val = select_cond
        if field in info_dict:
            try:
                field_val = float(info_dict[field]) if isinstance(info_dict[field], str) and '.' in info_dict[field] else int(info_dict[field])
                if op == '>' and not (field_val > val):
                    return False
                elif op == '>=' and not (field_val >= val):
                    return False
                elif op == '<' and not (field_val < val):
                    return False
                elif op == '<=' and not (field_val <= val):
                    return False
                elif op == '==' and not (field_val == val):
                    return False
                elif op == '!=' and not (field_val != val):
                    return False
            except (ValueError, TypeError):
                return False
        else:
            return False
    
    return True

def parse_vcf_line(line: str) -> Optional[Dict[str, Union[str, int]]]:
    """Parse single VCF line into variant dict"""
    if line.startswith('#'):
        return None
    
    fields = line.strip().split('\t')
    if len(fields) < 8:
        return None
    
    info_dict = parse_info_field(fields[7])
    
    # Required SV fields
    svtype = info_dict.get('SVTYPE', '')
    if not svtype:
        return None
    
    try:
        svlen = int(info_dict.get('SVLEN', 0))
        end = int(info_dict.get('END', int(fields[1]) + abs(svlen)))
    except (ValueError, TypeError):
        return None
    
    return {
        'chr': fields[0],
        'Original_start': int(fields[1]),
        'Original_end': end,
        'SV_type': svtype,
        'Len_SV': svlen,
        'ALT': fields[4],
        'REF': fields[3],
        'INFO': fields[7]
    }

def read_vcf_file(vcf_path: str) -> List[Dict[str, Union[str, int]]]:
    """Read and parse VCF file"""
    variants = []
    opener = gzip.open if vcf_path.endswith('.gz') else open
    
    start_time = time.time()
    with opener(vcf_path, 'rt') as f:
        for line in f:
            variant = parse_vcf_line(line)
            if not variant:
                continue
            variants.append(variant)
    
    end_time = time.time()
    elapsed = format_time(end_time - start_time)
    print(f"VCF loading completed in {elapsed}")
    return variants

def filter_variants(
    variants: List[Dict[str, Union[str, int]]],
    chrom_filter: str,
    min_len: int,
    sv_types: List[str],
    select_cond: Optional[Tuple[str, str, Union[float, int]]]
) -> List[Dict[str, Union[str, int]]]:
    """Filter variants based on criteria"""
    start_time = time.time()
    filtered = [
        var for var in variants 
        if passes_filters(
            var['chr'],
            var['SV_type'],
            var['Len_SV'],
            parse_info_field(var['INFO']),
            chrom_filter,
            min_len,
            sv_types,
            select_cond
        )
    ]
    end_time = time.time()
    elapsed = format_time(end_time - start_time)
    print(f"Variant filtering completed in {elapsed}")
    return filtered

def format_output_csv(raw_df: pd.DataFrame, chrom: str) -> pd.DataFrame:
    """Convert filtered variants to exact.py input format (CSV version)"""
    start_time = time.time()
    df = raw_df.copy()
    
    # Handle INS type (end = start for consistency)
    ins_mask = df["SV_type"] == "INS"
    df.loc[ins_mask, "Original_end"] = df.loc[ins_mask, "Original_start"]
    
    # Create exact.py format DataFrame (keeping 1-based coordinates for CSV)
    csv_df = pd.DataFrame({
        'Index': range(len(df)),
        'Index_con': chrom,  # Use chromosome as Index_con
        'SV_type': df['SV_type'],
        'Original_start': df['Original_start'],  # Keep as 1-based for CSV
        'Original_end': df['Original_end'],      # Keep as 1-based for CSV
        'Len_SV': df['Len_SV'],
        'New_start': -1,
        'New_end': -1,
        'New_len_SV': -1,
        'Balanced Trans Flag': -1,
        'relative start1': -1,
        'relative end1': -1,
        'relative start2': -1,
        'relative end2': -1
    })
    
    end_time = time.time()
    elapsed = format_time(end_time - start_time)
    print(f"CSV formatting completed in {elapsed}")
    return csv_df

def generate_output_filename(vcf_path: str, chrom: str, select: Optional[str], min_len: int, sv_types: List[str]) -> str:
    """Generate output filename with parameters (now CSV instead of BED)"""
    base_name = os.path.splitext(os.path.basename(vcf_path))[0]
    if vcf_path.endswith('.gz'):
        base_name = os.path.splitext(base_name)[0]  # Remove .gz if present
    
    filename_parts = [base_name, f"{chrom}"]
    
    # Add selection criteria if present
    if select:
        filename_parts.append(re.sub(r'[^\w\-_.]', '_', select))
    
    # Add other parameters
    filename_parts.append(f"min{min_len}")
    filename_parts.append("_".join(sorted(sv_types)))
    
    return f"{'_'.join(filename_parts)}.csv"  # Changed from .bed to .csv

def main():
    args = parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.save, exist_ok=True)
    print(f"Output will be saved to: {args.save}")
    
    try:
        select_cond = parse_selection_condition(args.select) if args.select else None
    except ValueError as e:
        print(f"Error: {e}")
        return
    
    total_start = time.time()
    try:
        print(f"Reading VCF file: {args.vcf}")
        variants = read_vcf_file(args.vcf)
        print(f"Loaded {len(variants)} variants before filtering")
        
        filtered = filter_variants(variants, args.chr, args.min_len, args.sv_types, select_cond)
        print(f"Found {len(filtered)} variants after filtering")
        
        if not filtered:
            print("No variants passed filters")
            return
            
        raw_df = pd.DataFrame(filtered)
        csv_df = format_output_csv(raw_df, args.chr)  # Changed from format_output_bed
        
        # Generate output filename
        out_name = generate_output_filename(
            args.vcf, 
            args.chr, 
            args.select, 
            args.min_len, 
            args.sv_types
        )
        out_path = os.path.join(args.save, out_name)
        
        # Save as CSV with header
        csv_df.to_csv(out_path, index=False, header=True)  # No need for sep="\t", uses comma by default
        print(f"Successfully saved {len(csv_df)} variants to:\n{out_path}")
        
    except Exception as e:
        print(f"Error: {str(e)}")
        if hasattr(e, '__traceback__'):
            import traceback
            traceback.print_exc()
    finally:
        total_end = time.time()
        elapsed = format_time(total_end - total_start)
        print(f"Total execution time: {elapsed}")

if __name__ == "__main__":
    main()