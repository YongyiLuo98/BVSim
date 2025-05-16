import argparse
import subprocess
import os
import yaml
import sys

def parse_args():
    # 第一阶段：解析 -config 参数
    initial_parser = argparse.ArgumentParser(add_help=False)
    initial_parser.add_argument(
        '-config', 
        type=str,
        default=os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'code') + '/bvsim_config.yaml', help='YAML config file path')
    config_args, remaining_argv = initial_parser.parse_known_args()

    # 主参数解析器
    parser = argparse.ArgumentParser(
        description='BVSim version 1.0.0',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # ================== 互斥操作模式 ==================
    action_group = parser.add_mutually_exclusive_group(required=True)
    action_group.add_argument('-vcf', action='store_true', help='Run VCF processing')
    action_group.add_argument('-exact', action='store_true', help='Exact variant generation')
    action_group.add_argument('-csv', action='store_true', help='Run CSV processing')
    action_group.add_argument('-wave', action='store_true', help='Run Wave model')
    action_group.add_argument('-wave_region', action='store_true', help='Run Wave region model')
    action_group.add_argument('-mimic', action='store_true', help='Realistic genome simulation (requires -hg38/-hg19)')
    action_group.add_argument('-uniform', action='store_true', help='Uniform distribution mode')

       # ================== 模式专用参数 ==================
    # Mimic模式参数
    parser.add_argument('-hg38', type=str, help='Chromosome name (e.g. chr1-chr22) for hg38 genome')
    parser.add_argument('-hg19', type=str, help='Chromosome name (e.g. chr1-chr22) for hg19 genome')
    parser.add_argument('-cell', action='store_true', help='Use CELL dataset list')
    parser.add_argument('-hgsvc', action='store_true', help='Use HGSVC dataset list')

    # Add config file argument
    parser.add_argument('-config', type=str, help='Path to YAML configuration file', default=os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'code')+ '/bvsim_config.yaml')
    parser.add_argument('-ref', type=str, help='Input reference local path', default=os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'empirical')+ '/sub_hg19_chr1.fasta')
    # parser.add_argument('-save', type=str, help='local path for saving', default=os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'save')+ '/')
    # 修改 -save 参数的定义
    save_default = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 
        '..', 
        'save'
    ) + '/'  # 默认值强制添加斜杠
    print(save_default)
    parser.add_argument(
        '-save', 
        type=str, 
        help='local path for saving', 
        default=save_default
    )
    parser.add_argument('-seed', type=int, help='Global seed for random number generator (non-negative integer)', default=999)
    parser.add_argument('-times', type=int, help='Maximum sampling times (positive integer)', default=10)
    parser.add_argument('-rep', type=int, help='Replication ID (non-negative integer for naming the files)', default=99)
    
    parser.add_argument('-seq_index', type=int, default=0, help='Index of sequence to use (0-based). Default: 0 (first sequence)')
    
    
    # VCF-specific arguments
    
    parser.add_argument('-vcf_file', type=str, help='Input VCF file path (required for VCF mode)')
    parser.add_argument('-chr', type=str, help='Target chromosome name to filter from VCF file (e.g., chr21) (required for VC mode)')
    parser.add_argument('-select', type=str, help='Selection criteria (e.g., "AF>0.001", "SVLEN>=100") (required for VC mode)')
    parser.add_argument('-min_len', type=int, default=50, help='Minimum SV length (bp) (positice integer, required for VC mode)')
    parser.add_argument('-sv_types', nargs='+', default=["DEL", "INS", "DUP", "INV"], 
                       help='SV types to include (required for VC mode)')
    
    # 变异生成模式选择
    
   
    parser.add_argument('-variant_table', type=str, help='Path to variant table CSV file (required for exact mode)')
    parser.add_argument('-validate_only', action='store_true',
                      help='Only validate input without generating variants')
    
    parser.add_argument('-sv_trans', type=int, help='Number of trans SV (non-negative integer)', default=5)
    parser.add_argument('-sv_inver', type=int, help='Number of inversion SV (non-negative integer)', default=5)
    parser.add_argument('-sv_dup', type=int, help='True duplication number (non-negative integer)', default=5)
    parser.add_argument('-sv_del', type=int, help='Number of deletion SV (non-negative integer)', default=5)
    parser.add_argument('-sv_ins', type=int, help='True insertion number (non-negative integer)', default=5)
    parser.add_argument('-snp', type=float, help='SNV number (non-negative integer) or probability (between 0 and 1)', default=5)
    parser.add_argument('-snv_del', type=float, help='SNV deletion number (non-negative integer) or probability (between 0 and 1)', default=5)
    parser.add_argument('-snv_ins', type=float, help='SNV insertion number (non-negative integer) or probability (between 0 and 1)', default=5)
    parser.add_argument('-notblockN', action='store_true', help='Do not Block N positions')
    parser.add_argument('-write', action='store_true', help='Write relative positions')
    parser.add_argument('-block_region_bed_url', '--block_region_bed_url', type=str, help='local path of the block region BED file', default=None)
    parser.add_argument('-cores', type=int, help='Number of kernels for parallel processing (positive integer, required for uniform-parallel/wave/wave-region mode to set up parallel computing)', default=1)
    
    parser.add_argument('-len_bins', type=int, help='Length of bins for parallel processing, must be positive integer and smaller than reference length (required for uniform-parallel/wave/wave-region mode to set up parallel computing)', default=50000)
    parser.add_argument('-delmin', type=int, help='Minimum deletion length (integer, not smaller than 50)', default=50)
    parser.add_argument('-delmax', type=int, help='Maximum deletion length (integer, larger than delmin)', default=60)
    parser.add_argument('-insmin', type=int, help='Minimum insertion length (integer, not smaller than 50)', default=50)
    parser.add_argument('-insmax', type=int, help='Maximum insertion length (integer, larger than insmin)', default=450)
    parser.add_argument('-dupmin', type=int, help='Minimum duplication length (integer, not smaller than 50)', default=50)
    parser.add_argument('-dupmax', type=int, help='Maximum duplication length (integer, larger than dupmin)', default=450)
    parser.add_argument('-invmin', type=int, help='Minimum inversion length (integer, not smaller than 50)', default=50)
    parser.add_argument('-invmax', type=int, help='Maximum inversion length (integer, larger than invmin)', default=450)
    parser.add_argument('-transmin', type=int, help='Minimum translocation length (integer, not smaller than 50)', default=50)
    parser.add_argument('-transmax', type=int, help='Maximum translocation length (integer, larger than transmin)', default=450) 
    
    #CSV
    parser.add_argument('-csv_num', type=int, help='Number for each type of CSV, superior to -csv_total_num (non-negative integer)', default=0)
    parser.add_argument('-csv_total_num', type=int, help='Total number for CSV, assign number of each type by empirical weights (non-negative integer)', default=0)
    parser.add_argument('-num_ID1_csv', type=int, help='Number of ID1 (TanInvDup) (non-negative integer)', default=5)
    parser.add_argument('-num_ID2_csv', type=int, help='Number of ID2 (DisInvDup) (non-negative integer)', default=5)
    parser.add_argument('-num_ID3_csv', type=int, help='Number of ID3 (dispersed duplications) (non-negative integer)', default=5)
    parser.add_argument('-num_ID4_csv', type=int, help='Number of ID4 (DelInv+InvDel) (non-negative integer)', default=5)
    parser.add_argument('-num_ID5_csv', type=int, help='Number of ID5 (DEL+ DisInvDup) (non-negative integer)', default=5)
    parser.add_argument('-num_ID6_csv', type=int, help='Number of ID6 (DEL+ DisDup) (non-negative integer)', default=5)
    parser.add_argument('-num_ID7_csv', type=int, help='Number of ID7 (TanDup+DEL) (non-negative integer)', default=5)
    parser.add_argument('-num_ID8_csv', type=int, help='Number of ID8 (TanInvDup+DEL) (non-negative integer)', default=5)
    parser.add_argument('-num_ID9_csv', type=int, help='Number of ID9 (TanDup + DEL + INV) (non-negative integer)', default=5)
    parser.add_argument('-num_ID10_csv', type=int, help='Number of ID10 (TanInvDup + DEL + INV) (non-negative integer)', default=5)
    #ID11-18
    parser.add_argument('-num_ID11_csv', type=int, help='Number of ID11: paired-Deletion Inversion (delInvdel) (non-negative integer)', default=5)
    parser.add_argument('-num_ID12_csv', type=int, help='Number of ID12: Inversion with 5 Flanking Duplication (dupInv) (non-negative integer)', default=5)
    parser.add_argument('-num_ID13_csv', type=int, help='Number of ID13: Inversion with 3 Flanking Duplication (Invdup) (non-negative integer)', default=5)
    parser.add_argument('-num_ID14_csv', type=int, help='Number of ID14: Paired-duplication inversion (dupInvdup) (non-negative integer)', default=5)
    parser.add_argument('-num_ID15_csv', type=int, help='Number of ID15: Inversion with 5 Flanking Duplication and 3 Flanking Deletion (dupInvdel) (non-negative integer)', default=5)
    parser.add_argument('-num_ID16_csv', type=int, help='Number of ID16: Inversion with 5 Flanking Deletion and 3 Flanking Duplication (delInvdup) (non-negative integer)', default=5)
    parser.add_argument('-num_ID17_csv', type=int, help='Number of ID17: Inverted Duplication with Flanking Triplication (dupTRIPdup-INV) (non-negative integer)', default=5)
    parser.add_argument('-num_ID18_csv', type=int, help='Number of ID18: Insertion with Deletion (INSdel) (non-negative integer)', default=5)
    #define length
    parser.add_argument('-mu_ID1', type=int, help='Mean of length for CSV ID1 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID1', type=int, help='Sigma of length for CSV ID1 (non-negative integer)', default=100)

    parser.add_argument('-mu_ID2', type=int, help='Mean of length for CSV ID2 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID2', type=int, help='Sigma of length for CSV ID2 (non-negative integer)', default=100)

    parser.add_argument('-mu_ID3', type=int, help='Mean of length for CSV ID3 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID3', type=int, help='Sigma of length for CSV ID3 (non-negative integer)', default=100)

    parser.add_argument('-mu_ID4', type=int, help='Mean of length for CSV ID4 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID4', type=int, help='Sigma of length for CSV ID4 (non-negative integer)', default=100)

    parser.add_argument('-mu_ID5', type=int, help='Mean of length for CSV ID5 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID5', type=int, help='Sigma of length for CSV ID5 (non-negative integer)', default=100)

    parser.add_argument('-mu_ID6', type=int, help='Mean of length for CSV ID6 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID6', type=int, help='Sigma of length for CSV ID6 (non-negative integer)', default=100)

    parser.add_argument('-mu_ID7', type=int, help='Mean of length for CSV ID7 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID7', type=int, help='Sigma of length for CSV ID7 (non-negative integer)', default=100)

    parser.add_argument('-mu_ID8', type=int, help='Mean of length for CSV ID8 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID8', type=int, help='Sigma of length for CSV ID8 (non-negative integer)', default=100)

    parser.add_argument('-mu_ID9', type=int, help='Mean of length for CSV ID9 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID9', type=int, help='Sigma of length for CSV ID9 (non-negative integer)', default=100)

    parser.add_argument('-mu_ID10', type=int, help='Mean of length for CSV ID10 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID10', type=int, help='Sigma of length for CSV ID10 (non-negative integer)', default=100)

    parser.add_argument('-mu_ID11', type=int, help='Mean of length for CSV ID11 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID11', type=int, help='Sigma of length for CSV ID11 (non-negative integer)', default=100)
    
    parser.add_argument('-mu_ID12', type=int, help='Mean of length for CSV ID12 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID12', type=int, help='Sigma of length for CSV ID12 (non-negative integer)', default=100)
    
    parser.add_argument('-mu_ID13', type=int, help='Mean of length for CSV ID13 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID13', type=int, help='Sigma of length for CSV ID13 (non-negative integer)', default=100)
    
    parser.add_argument('-mu_ID14', type=int, help='Mean of length for CSV ID14 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID14', type=int, help='Sigma of length for CSV ID14 (non-negative integer)', default=100)
    
    parser.add_argument('-mu_ID15', type=int, help='Mean of length for CSV ID15 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID15', type=int, help='Sigma of length for CSV ID15 (non-negative integer)', default=100)
    
    parser.add_argument('-mu_ID16', type=int, help='Mean of length for CSV ID16 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID16', type=int, help='Sigma of length for CSV ID16 (non-negative integer)', default=100)
    
    parser.add_argument('-mu_ID17', type=int, help='Mean of length for CSV ID17 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID17', type=int, help='Sigma of length for CSV ID17 (non-negative integer)', default=100)
    
    parser.add_argument('-mu_ID18', type=int, help='Mean of length for CSV ID18 (integer, larger than 100)', default=1000)
    parser.add_argument('-sigma_ID18', type=int, help='Sigma of length for CSV ID18 (non-negative integer)', default=100)
    
    
    parser.add_argument('-mode', type=str, help='Mode for calculating number of SVs per bin (empirical/probability)', default='probability')
    parser.add_argument('-sum', action='store_true', help='total indel SV equals sum of the input (single sample) or mean of the input (multiple samples)')
    parser.add_argument('-indel_input_bed', type=str, help='Input BED file for indels (required if input single sample for wave or wave-region mode)',default=None)
    parser.add_argument('-file_list', type=str, nargs='+', default=['NA19240_chr21', 'HG02818_chr21', 'NA19434_chr21'],
                        help='List of sample files (e.g. NA19240_chr21.bed, HG02818_chr21.bed, NA19434_chr21.bed in empirical folder) (required if multiple samples for wave or wave-region mode)')

    
    parser.add_argument('-p_del_region', type=float, help='Probability of SV DEL (between 0 and 1) in the user-defined region for deletion (required for wave-region mode)', default=0.5)
    parser.add_argument('-p_ins_region',  type=float, help='Probability of SV INS (between 0 and 1) in the user-defined region for insertion (required for wave-region mode)', default=0.5)
    parser.add_argument('-region_bed_url', type=str, help='local path of the BED file for the user-defined region (required for wave-region mode)', default=None)
   
    
    
    # 第三阶段：加载 YAML 配置
    final_args = {}
    if config_args.config and os.path.exists(config_args.config):
        try:
            with open(config_args.config, 'r') as f:
                yaml_config = yaml.safe_load(f)
                # print("Loaded YAML config:", yaml_config)
                
                # ================== 展平嵌套结构 ==================
                flat_config = {}
                for section, values in yaml_config.items():
                    if isinstance(values, dict):
                        for k, v in values.items():
                            flat_key = f"{section}_{k}" if section not in ['general', 'variants'] else k
                            flat_config[flat_key] = v
                    else:
                        flat_config[section] = values
                        
                # ================== 新增代码：路径转换（操作 flat_config） ==================
                # 获取项目根目录
                main_dir = os.path.dirname(os.path.realpath(__file__))  # main.py 所在目录：BVSim/main/
                project_root = os.path.abspath(os.path.join(main_dir, ".."))  # 项目根目录
                
                def resolve_path(rel_path):
                    return os.path.abspath(os.path.join(project_root, rel_path))
                
                # 处理路径转换
                path_keys = ['ref', 'save', 'block_region_bed_url']
                for key in path_keys:
                    if key in flat_config and flat_config[key]:
                        flat_config[key] = resolve_path(flat_config[key])
                
                # ================== 应用配置 ==================
                for key, value in flat_config.items():
                    if value is not None:
                        final_args[key] = value

        except Exception as e:
            print(f"Error loading config: {e}")
            sys.exit(1)

    # 合并参数（命令行参数优先）
    args = parser.parse_args(args=remaining_argv, namespace=argparse.Namespace(**final_args))
    from pathlib import Path
    
    # 处理 -save 参数
    save_path = Path(args.save).resolve()  # 解析绝对路径并标准化
    args.save = str(save_path).rstrip('/') + '/'  # 确保以 / 结尾
    return args

def validate_arguments(args):
    """参数校验逻辑"""
    # Chromosome格式校验
    chroms = [f'chr{i}' for i in range(1,23)]
    if args.hg38 and (not args.hg38.startswith('chr') or args.hg38 not in chroms):
        raise argparse.ArgumentError('-hg38', f'Must be one of {chroms}')
    if args.hg19 and (not args.hg19.startswith('chr') or args.hg19 not in chroms):
        raise argparse.ArgumentError('-hg19', f'Must be one of {chroms}')

    # 互斥校验
    if args.cell and args.hgsvc:
        raise argparse.ArgumentError(None, '-cell and -hgsvc are mutually exclusive')
    if args.mimic and not (args.hg38 or args.hg19):
        raise argparse.ArgumentError('-mimic', 'Requires -hg38 or -hg19')
    if args.hg38 and args.hg19:
        raise argparse.ArgumentError(None, '-hg38 and -hg19 are mutually exclusive')

def main():
    args = parse_args()
    try:
        validate_arguments(args)
    except argparse.ArgumentError as e:
        print(f"参数错误: {e}")
        sys.exit(1)

    # 获取脚本目录
    main_dir = os.path.dirname(os.path.realpath(__file__))
    
    # 构建基础命令
    if args.vcf:  # 现在args.vcf直接包含文件路径
        script_name = "VCF.py"
        cmd = [
            "python", 
            os.path.join(main_dir, script_name),
            "-vcf_file", os.path.abspath(args.vcf_file),  # 直接使用args.vcf作为文件路径
            "-save", args.save,
            "-min_len", str(args.min_len),
            "-sv_types"
        ] + args.sv_types
        
        if hasattr(args, 'chr') and args.chr:
            cmd.extend(["-chr", args.chr])
        if hasattr(args, 'select') and args.select:
            cmd.extend(["-select", args.select])
            
        # 执行VCF.py并退出
        subprocess.run(cmd)
        return  # 确保执行后直接退出
    
    elif args.exact:
        # print("Entering exact mode...")  # 确认进入exact模式
        script_name = "exact.py"
        # 构建命令行参数列表
        cmd = [
            "python", 
            os.path.join(main_dir, script_name),
            "-ref", args.ref,
            "-seed", str(args.seed), 
             "-rep", str(args.rep), 
            "-variant_table", args.variant_table,
            "-save", args.save,
            "-seq_index", str(args.seq_index)
        ]
        
        # 添加可选参数
        if args.seed is not None:
            cmd.extend(["-seed", str(args.seed)])
        if args.block_region_bed_url is not None:
            cmd.extend(["-block_region_bed_url", args.block_region_bed_url])
        if args.validate_only:
            cmd.append("-validate_only")
            
        # # 执行exact.py并退出
        # subprocess.run(cmd)
        # return  # 确保执行后直接退出
        
    
    elif args.csv:
        script_name = "CSV.py"
        cmd = [
    "python", os.path.join(main_dir, script_name), 
    "-ref", args.ref, 
    "-save", args.save, 
    "-seed", str(args.seed), 
    "-times", str(args.times), 
    "-rep", str(args.rep), 
    "-sv_trans", str(args.sv_trans), 
    "-sv_inver", str(args.sv_inver), 
    "-sv_dup", str(args.sv_dup), 
    "-sv_del", str(args.sv_del), 
    "-sv_ins", str(args.sv_ins), 
    "-delmin", str(args.delmin), 
    "-delmax", str(args.delmax), 
    "-insmin", str(args.insmin), 
    "-insmax", str(args.insmax),
    "-dupmin", str(args.dupmin), 
    "-dupmax", str(args.dupmax), 
    "-invmin", str(args.invmin), 
    "-invmax", str(args.invmax),
    "-transmin", str(args.transmin), 
    "-transmax", str(args.transmax),
    "-csv_num", str(args.csv_num), 
    "-csv_total_num", str(args.csv_total_num),
    "-num_ID1_csv", str(args.num_ID1_csv), "-mu_ID1", str(args.mu_ID1), "-sigma_ID1", str(args.sigma_ID1),
    "-num_ID2_csv", str(args.num_ID2_csv), "-mu_ID2", str(args.mu_ID2), "-sigma_ID2", str(args.sigma_ID2),
    "-num_ID3_csv", str(args.num_ID3_csv), "-mu_ID3", str(args.mu_ID3), "-sigma_ID3", str(args.sigma_ID3),
    "-num_ID4_csv", str(args.num_ID4_csv), "-mu_ID4", str(args.mu_ID4), "-sigma_ID4", str(args.sigma_ID4),
    "-num_ID5_csv", str(args.num_ID5_csv), "-mu_ID5", str(args.mu_ID5), "-sigma_ID5", str(args.sigma_ID5),
    "-num_ID6_csv", str(args.num_ID6_csv), "-mu_ID6", str(args.mu_ID6), "-sigma_ID6", str(args.sigma_ID6),
    "-num_ID7_csv", str(args.num_ID7_csv), "-mu_ID7", str(args.mu_ID7), "-sigma_ID7", str(args.sigma_ID7),
    "-num_ID8_csv", str(args.num_ID8_csv), "-mu_ID8", str(args.mu_ID8), "-sigma_ID8", str(args.sigma_ID8),
    "-num_ID9_csv", str(args.num_ID9_csv), "-mu_ID9", str(args.mu_ID9), "-sigma_ID9", str(args.sigma_ID9),
    "-num_ID10_csv", str(args.num_ID10_csv), "-mu_ID10", str(args.mu_ID10), "-sigma_ID10", str(args.sigma_ID10),
    "-num_ID11_csv", str(args.num_ID11_csv), "-mu_ID11", str(args.mu_ID11), "-sigma_ID11", str(args.sigma_ID11),
    "-num_ID12_csv", str(args.num_ID12_csv), "-mu_ID12", str(args.mu_ID12), "-sigma_ID12", str(args.sigma_ID12),
    "-num_ID13_csv", str(args.num_ID13_csv), "-mu_ID13", str(args.mu_ID13), "-sigma_ID13", str(args.sigma_ID13),
    "-num_ID14_csv", str(args.num_ID14_csv), "-mu_ID14", str(args.mu_ID14), "-sigma_ID14", str(args.sigma_ID14),
    "-num_ID15_csv", str(args.num_ID15_csv), "-mu_ID15", str(args.mu_ID15), "-sigma_ID15", str(args.sigma_ID15),
    "-num_ID16_csv", str(args.num_ID16_csv), "-mu_ID16", str(args.mu_ID16), "-sigma_ID16", str(args.sigma_ID16),
    "-num_ID17_csv", str(args.num_ID17_csv), "-mu_ID17", str(args.mu_ID17), "-sigma_ID17", str(args.sigma_ID17),
    "-num_ID18_csv", str(args.num_ID18_csv), "-mu_ID18", str(args.mu_ID18), "-sigma_ID18", str(args.sigma_ID18),
    "-seq_index", str(args.seq_index)]
        
        # # 执行exact.py并退出
        # subprocess.run(cmd)
        # return  # 确保执行后直接退出

    elif args.wave:
        script_name = "Wave.py"
        cmd = ["python", os.path.join(main_dir, script_name), "-ref", args.ref, "-save", args.save, "-seed", str(args.seed), 
               "-times", str(args.times), "-rep", str(args.rep), "-sv_trans", str(args.sv_trans), "-sv_inver", str(args.sv_inver), "-sv_dup", 
               str(args.sv_dup), "-sv_del", str(args.sv_del), "-sv_ins", str(args.sv_ins), "-cores", str(args.cores), "-len_bins", 
               str(args.len_bins), 
               "-delmin", str(args.delmin), 
               "-delmax", str(args.delmax), 
               "-insmin", str(args.insmin), 
               "-insmax", str(args.insmax),
               "-dupmin", str(args.dupmin), 
               "-dupmax", str(args.dupmax), 
               "-invmin", str(args.invmin), 
                "-invmax", str(args.invmax),
                "-transmin", str(args.transmin), 
                "-transmax", str(args.transmax),
               "-mode", args.mode,
                "-seq_index", str(args.seq_index)]
        # 添加 indel_input_bed 参数
        if args.block_region_bed_url:
            cmd.extend(["-block_region_bed_url", str(args.block_region_bed_url)])
        if args.indel_input_bed:
            cmd.extend(["-indel_input_bed", args.indel_input_bed])
        # 添加 file_list 参数
        cmd.extend(["-file_list"] + args.file_list)  # 将 file_list 展开为多个参数
        if args.sum:
            cmd.append("-sum")
    elif args.wave_region:
        script_name = "Wave_TR.py"
        cmd = ["python", os.path.join(main_dir, script_name), "-ref", args.ref, "-save", args.save, "-seed", str(args.seed), 
               "-times", str(args.times), "-rep", str(args.rep), "-sv_trans", str(args.sv_trans), "-sv_inver", str(args.sv_inver), "-sv_dup", 
               str(args.sv_dup), "-sv_del", str(args.sv_del), "-sv_ins", str(args.sv_ins), "-cores", str(args.cores), "-len_bins", 
               str(args.len_bins), 
               "-delmin", str(args.delmin), 
               "-delmax", str(args.delmax), 
               "-insmin", str(args.insmin), 
               "-insmax", str(args.insmax),
                "-dupmin", str(args.dupmin), 
                "-dupmax", str(args.dupmax), 
                "-invmin", str(args.invmin), 
                "-invmax", str(args.invmax),
                "-transmin", str(args.transmin), 
                "-transmax", str(args.transmax),
               "-mode", args.mode, 
               "-p_del_region", str(args.p_del_region), 
               "-p_ins_region", str(args.p_ins_region),
                "-seq_index", str(args.seq_index)]
        # 添加 indel_input_bed 参数
        if args.block_region_bed_url:
            cmd.extend(["-block_region_bed_url", str(args.block_region_bed_url)])
        if args.region_bed_url:
            cmd.extend(["-region_bed_url", args.region_bed_url])
        if args.indel_input_bed:
            cmd.extend(["-indel_input_bed", args.indel_input_bed])
        # 添加 file_list 参数
        cmd.extend(["-file_list"] + args.file_list)  # 将 file_list 展开为多个参数
        if args.sum:
            cmd.append("-sum")
    elif args.mimic:
        if args.hg38:
            script_name = "Wave_hg38_TR.py"
            cmd = ["python", os.path.join(main_dir, script_name), "-ref", args.ref, "-save", args.save, "-seed", str(args.seed), 
                "-times", str(args.times), "-rep", str(args.rep), "-sv_trans", str(args.sv_trans), "-sv_inver", str(args.sv_inver), "-sv_dup", 
                str(args.sv_dup), "-sv_del", str(args.sv_del), "-sv_ins", str(args.sv_ins), "-cores", str(args.cores), "-len_bins", 
                str(args.len_bins), 
                "-delmin", str(args.delmin), 
                "-delmax", str(args.delmax), 
                "-insmin", str(args.insmin), 
                "-insmax", str(args.insmax),
                "-dupmin", str(args.dupmin), 
                    "-dupmax", str(args.dupmax), 
                    "-invmin", str(args.invmin), 
                    "-invmax", str(args.invmax),
                    "-transmin", str(args.transmin), 
                    "-transmax", str(args.transmax),
                "-mode", args.mode, "-p_del_region", str(args.p_del_region), "-p_ins_region", str(args.p_ins_region),
                "-hg38", str(args.hg38),
                    "-seq_index", str(args.seq_index)]
            # 添加 indel_input_bed 参数
            if args.block_region_bed_url:
                cmd.extend(["-block_region_bed_url", str(args.block_region_bed_url)])
            if args.region_bed_url:
                cmd.extend(["-region_bed_url", args.region_bed_url])
            if args.indel_input_bed:
                cmd.extend(["-indel_input_bed", args.indel_input_bed])
            if args.sum:
                cmd.append("-sum")
            if args.cell:
                cmd.append("-cell")
            elif args.hgsvc:
                cmd.append("-hgsvc")
        else:
            script_name = "Wave_hg19_TR.py"
            cmd = ["python", os.path.join(main_dir, script_name), "-ref", args.ref, "-save", args.save, "-seed", str(args.seed), 
                "-times", str(args.times), "-rep", str(args.rep), "-sv_trans", str(args.sv_trans), "-sv_inver", str(args.sv_inver), "-sv_dup", 
                str(args.sv_dup), "-sv_del", str(args.sv_del), "-sv_ins", str(args.sv_ins), "-cores", str(args.cores), "-len_bins", 
                str(args.len_bins), 
                "-delmin", str(args.delmin), 
                "-delmax", str(args.delmax), 
                "-insmin", str(args.insmin), 
                "-insmax", str(args.insmax),
                "-dupmin", str(args.dupmin), 
                    "-dupmax", str(args.dupmax), 
                    "-invmin", str(args.invmin), 
                    "-invmax", str(args.invmax),
                    "-transmin", str(args.transmin), 
                    "-transmax", str(args.transmax),
                "-mode", args.mode, "-p_del_region", str(args.p_del_region), "-p_ins_region", str(args.p_ins_region),
                "-hg19", str(args.hg19),
                    "-seq_index", str(args.seq_index)]
            # 添加 indel_input_bed 参数
            if args.block_region_bed_url:
                cmd.extend(["-block_region_bed_url", str(args.block_region_bed_url)])
            if args.region_bed_url:
                cmd.extend(["-region_bed_url", args.region_bed_url])
            if args.indel_input_bed:
                cmd.extend(["-indel_input_bed", args.indel_input_bed])
            if args.sum:
                cmd.append("-sum")

    elif args.uniform:
        if  args.cores is not None and args.cores > 1:
            script_name = "uniform_parallel.py"
            cmd = ["python", os.path.join(main_dir, script_name), "-ref", args.ref, "-save", args.save, "-seed", str(args.seed), 
                "-times", str(args.times), "-rep", str(args.rep), "-sv_trans", str(args.sv_trans), "-sv_inver", str(args.sv_inver), "-sv_dup", 
                str(args.sv_dup), "-sv_del", str(args.sv_del), "-sv_ins", str(args.sv_ins), "-cores", str(args.cores), "-len_bins", 
                str(args.len_bins), 
                "-delmin", str(args.delmin), 
                "-delmax", str(args.delmax), 
                "-insmin", str(args.insmin), 
                "-insmax", str(args.insmax),
                "-dupmin", str(args.dupmin), 
                    "-dupmax", str(args.dupmax), 
                    "-invmin", str(args.invmin), 
                    "-invmax", str(args.invmax),
                    "-transmin", str(args.transmin), 
                    "-transmax", str(args.transmax),
                "-block_region_bed_url", str(args.block_region_bed_url),
                    "-seq_index", str(args.seq_index)]
        else:
            script_name = "uniform.py"
            cmd = ["python", os.path.join(main_dir, script_name), "-ref", args.ref, "-save", args.save, 
                "-seed", str(args.seed), "-times", str(args.times), "-rep", str(args.rep), "-sv_trans", str(args.sv_trans), 
                "-sv_inver", str(args.sv_inver), "-sv_dup", str(args.sv_dup), "-sv_del", str(args.sv_del), "-sv_ins", str(args.sv_ins), 
                "-snp", str(args.snp), "-snv_del", str(args.snv_del), "-snv_ins", str(args.snv_ins), 
                "-delmin", str(args.delmin), 
                    "-delmax", str(args.delmax), 
                    "-insmin", str(args.insmin), 
                    "-insmax", str(args.insmax),
                    "-dupmin", str(args.dupmin), 
                    "-dupmax", str(args.dupmax), 
                    "-invmin", str(args.invmin), 
                    "-invmax", str(args.invmax),
                    "-transmin", str(args.transmin), 
                    "-transmax", str(args.transmax),
                "-block_region_bed_url", str(args.block_region_bed_url),
                    "-seq_index", str(args.seq_index)]

    # ================== 添加公共参数 ==================
    common_params = [
        "-snp", str(args.snp),
        "-snv_del", str(args.snv_del),
        "-snv_ins", str(args.snv_ins)
    ]
    if args.notblockN:
        common_params.append("-notblockN")
    if args.write:
        common_params.append("-write")
    
    cmd += common_params

    # 执行命令
    print("Run:", ' '.join(cmd))
    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()