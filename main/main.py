import argparse
import subprocess
import os

def parse_args():
    parser = argparse.ArgumentParser(description='GenoWave')
    parser.add_argument('-ref', type=str, help='Input reference local path', default=os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'empirical')+ '/sub_hg19_chr1.fasta')
    parser.add_argument('-save', type=str, help='local path for saving', default=os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'save')+ '/')
    parser.add_argument('-seed', type=int, help='Seed for random number generator', default=999)
    parser.add_argument('-times', type=int, help='Number of times', default=10)
    parser.add_argument('-rep', type=int, help='Replication ID', default=5)
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
    parser.add_argument('-block_region_bed_url', '--block_region_bed_url', type=str, help='local path of the block region BED file', default=None)
    parser.add_argument('-cores', type=int, help='Number of kernels for parallel processing', default=1)
    parser.add_argument('-len_bins', type=int, help='Length of bins for parallel processing', default=50000)
    parser.add_argument('-delmin', type=int, help='Minimum deletion length', default=50)
    parser.add_argument('-delmax', type=int, help='Maximum deletion length', default=60)
    parser.add_argument('-insmin', type=int, help='Minimum insertion length', default=50)
    parser.add_argument('-insmax', type=int, help='Maximum insertion length', default=450)
    
    parser.add_argument('-csv', action='store_true', help='Run csv.py script')
    #CSV
    parser.add_argument('-csv_num', type=int, help='Number for each type of CSV, superior to -csv_total_num', default=0)
    parser.add_argument('-csv_total_num', type=int, help='Total number for CSV, assign number of each type by empirical weights', default=0)
    parser.add_argument('-num_ID1_csv', type=int, help='Number of ID1 (TanInvDup)', default=5)
    parser.add_argument('-num_ID2_csv', type=int, help='Number of ID2 (DisInvDup)', default=5)
    parser.add_argument('-num_ID3_csv', type=int, help='Number of ID3 (dispersed duplications)', default=5)
    parser.add_argument('-num_ID4_csv', type=int, help='Number of ID4 (DelInv+InvDel)', default=5)
    parser.add_argument('-num_ID5_csv', type=int, help='Number of ID5 (DEL+ DisInvDup)', default=5)
    parser.add_argument('-num_ID6_csv', type=int, help='Number of ID6 (DEL+ DisDup)', default=5)
    parser.add_argument('-num_ID7_csv', type=int, help='Number of ID7 (TanDup+DEL)', default=5)
    parser.add_argument('-num_ID8_csv', type=int, help='Number of ID8 (TanInvDup+DEL)', default=5)
    parser.add_argument('-num_ID9_csv', type=int, help='Number of ID9 (TanDup + DEL + INV)', default=5)
    parser.add_argument('-num_ID10_csv', type=int, help='Number of ID10 (TanInvDup + DEL + INV)', default=5)
    #ID11-18
    parser.add_argument('-num_ID11_csv', type=int, help='Number of ID11: paired-Deletion Inversion (delInvdel)', default=5)
    parser.add_argument('-num_ID12_csv', type=int, help='Number of ID12: Inversion with 5 Flanking Duplication (dupInv)', default=5)
    parser.add_argument('-num_ID13_csv', type=int, help='Number of ID13: Inversion with 3 Flanking Duplication (Invdup)', default=5)
    parser.add_argument('-num_ID14_csv', type=int, help='Number of ID14: Paired-duplication inversion (dupInvdup)', default=5)
    parser.add_argument('-num_ID15_csv', type=int, help='Number of ID15: Inversion with 5 Flanking Duplication and 3 Flanking Deletion (dupInvdel)', default=5)
    parser.add_argument('-num_ID16_csv', type=int, help='Number of ID16: Inversion with 5 Flanking Deletion and 3 Flanking Duplication (delInvdup)', default=5)
    parser.add_argument('-num_ID17_csv', type=int, help='Number of ID17: Inverted Duplication with Flanking Triplication (dupTRIPdup-INV) ', default=5)
    parser.add_argument('-num_ID18_csv', type=int, help='Number of ID18: Insertion with Deletion (INSdel)', default=5)
    #define length
    parser.add_argument('-mu_ID1', type=int, help='Mean of length for CSV ID1', default=1000)
    parser.add_argument('-sigma_ID1', type=int, help='Sigma of length for CSV ID1', default=100)

    parser.add_argument('-mu_ID2', type=int, help='Mean of length for CSV ID2', default=1000)
    parser.add_argument('-sigma_ID2', type=int, help='Sigma of length for CSV ID2', default=100)

    parser.add_argument('-mu_ID3', type=int, help='Mean of length for CSV ID3', default=1000)
    parser.add_argument('-sigma_ID3', type=int, help='Sigma of length for CSV ID3', default=100)

    parser.add_argument('-mu_ID4', type=int, help='Mean of length for CSV ID4', default=1000)
    parser.add_argument('-sigma_ID4', type=int, help='Sigma of length for CSV ID4', default=100)

    parser.add_argument('-mu_ID5', type=int, help='Mean of length for CSV ID5', default=1000)
    parser.add_argument('-sigma_ID5', type=int, help='Sigma of length for CSV ID5', default=100)

    parser.add_argument('-mu_ID6', type=int, help='Mean of length for CSV ID6', default=1000)
    parser.add_argument('-sigma_ID6', type=int, help='Sigma of length for CSV ID6', default=100)

    parser.add_argument('-mu_ID7', type=int, help='Mean of length for CSV ID7', default=1000)
    parser.add_argument('-sigma_ID7', type=int, help='Sigma of length for CSV ID7', default=100)

    parser.add_argument('-mu_ID8', type=int, help='Mean of length for CSV ID8', default=1000)
    parser.add_argument('-sigma_ID8', type=int, help='Sigma of length for CSV ID8', default=100)

    parser.add_argument('-mu_ID9', type=int, help='Mean of length for CSV ID9', default=1000)
    parser.add_argument('-sigma_ID9', type=int, help='Sigma of length for CSV ID9', default=100)

    parser.add_argument('-mu_ID10', type=int, help='Mean of length for CSV ID10', default=1000)
    parser.add_argument('-sigma_ID10', type=int, help='Sigma of length for CSV ID10', default=100)

    parser.add_argument('-mu_ID11', type=int, help='Mean of length for CSV ID11', default=1000)
    parser.add_argument('-sigma_ID11', type=int, help='Sigma of length for CSV ID11', default=100)
    
    parser.add_argument('-mu_ID12', type=int, help='Mean of length for CSV ID12', default=1000)
    parser.add_argument('-sigma_ID12', type=int, help='Sigma of length for CSV ID12', default=100)
    
    parser.add_argument('-mu_ID13', type=int, help='Mean of length for CSV ID13', default=1000)
    parser.add_argument('-sigma_ID13', type=int, help='Sigma of length for CSV ID13', default=100)
    
    parser.add_argument('-mu_ID14', type=int, help='Mean of length for CSV ID14', default=1000)
    parser.add_argument('-sigma_ID14', type=int, help='Sigma of length for CSV ID14', default=100)
    
    parser.add_argument('-mu_ID15', type=int, help='Mean of length for CSV ID15', default=1000)
    parser.add_argument('-sigma_ID15', type=int, help='Sigma of length for CSV ID15', default=100)
    
    parser.add_argument('-mu_ID16', type=int, help='Mean of length for CSV ID16', default=1000)
    parser.add_argument('-sigma_ID16', type=int, help='Sigma of length for CSV ID16', default=100)
    
    parser.add_argument('-mu_ID17', type=int, help='Mean of length for CSV ID17', default=1000)
    parser.add_argument('-sigma_ID17', type=int, help='Sigma of length for CSV ID17', default=100)
    
    parser.add_argument('-mu_ID18', type=int, help='Mean of length for CSV ID18', default=1000)
    parser.add_argument('-sigma_ID18', type=int, help='Sigma of length for CSV ID18', default=100)
    
    parser.add_argument('-wave', action='store_true', help='Run Wave.py script')
    parser.add_argument('-mode', type=str, help='Mode for calculating probabilities', default='probability')
    parser.add_argument('-sum', action='store_true', help='total indel SV equals sum of the input bed')
    parser.add_argument('-indel_input_bed', type=str, help='Input BED file for indels', default='~/data/test_data/TGS/hg002/chr21_SV_Tier1.bed')

    parser.add_argument('-wave_region', action='store_true', help='Run Wave_TR.py script')
    parser.add_argument('-p_del_region', type=float, help='Probability of SV DEL in the user-defined region for deletion', default=0.5)
    parser.add_argument('-p_ins_region',  type=float, help='Probability of SV INS in the user-defined region for insertion', default=0.5)
    parser.add_argument('-region_bed_url', type=str, help='local path of the BED file for the user-defined region', default='~/data/test_data/TGS/hg002/chr21_TR_unique.bed')
    parser.add_argument('-hg38', type=str, help='Chromosome name', required=False)
    return parser.parse_args()
    

def main():
    args = parse_args()
    
     # 获取main.py文件所在的路径
    main_dir = os.path.dirname(os.path.realpath(__file__))
    
    if args.csv:
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
    "-num_ID18_csv", str(args.num_ID18_csv), "-mu_ID18", str(args.mu_ID18), "-sigma_ID18", str(args.sigma_ID18)]

    elif args.wave:
        script_name = "Wave.py"
        cmd = ["python", os.path.join(main_dir, script_name), "-ref", args.ref, "-save", args.save, "-seed", str(args.seed), 
               "-times", str(args.times), "-rep", str(args.rep), "-sv_trans", str(args.sv_trans), "-sv_inver", str(args.sv_inver), "-sv_dup", 
               str(args.sv_dup), "-sv_del", str(args.sv_del), "-sv_ins", str(args.sv_ins), "-cores", str(args.cores), "-len_bins", 
               str(args.len_bins), "-delmin", str(args.delmin), "-delmax", str(args.delmax), "-insmin", str(args.insmin), "-insmax", str(args.insmax),
               "-mode", args.mode, "-indel_input_bed", args.indel_input_bed, "-block_region_bed_url", str(args.block_region_bed_url)]
        if args.sum:
            cmd.append("-sum")
    elif args.wave_region:
        script_name = "Wave_TR.py"
        cmd = ["python", os.path.join(main_dir, script_name), "-ref", args.ref, "-save", args.save, "-seed", str(args.seed), 
               "-times", str(args.times), "-rep", str(args.rep), "-sv_trans", str(args.sv_trans), "-sv_inver", str(args.sv_inver), "-sv_dup", 
               str(args.sv_dup), "-sv_del", str(args.sv_del), "-sv_ins", str(args.sv_ins), "-cores", str(args.cores), "-len_bins", 
               str(args.len_bins), "-delmin", str(args.delmin), "-delmax", str(args.delmax), "-insmin", str(args.insmin), "-insmax", str(args.insmax),
               "-mode", args.mode, "-indel_input_bed", args.indel_input_bed, "-p_del_region", str(args.p_del_region), "-p_ins_region", str(args.p_ins_region),
               "-region_bed_url", args.region_bed_url, "-block_region_bed_url", str(args.block_region_bed_url)]
        if args.sum:
            cmd.append("-sum")
            
    elif args.hg38:
        script_name = "Wave_hg38_TR.py"
        cmd = ["python", os.path.join(main_dir, script_name), "-ref", args.ref, "-save", args.save, "-seed", str(args.seed), 
               "-times", str(args.times), "-rep", str(args.rep), "-sv_trans", str(args.sv_trans), "-sv_inver", str(args.sv_inver), "-sv_dup", 
               str(args.sv_dup), "-sv_del", str(args.sv_del), "-sv_ins", str(args.sv_ins), "-cores", str(args.cores), "-len_bins", 
               str(args.len_bins), "-delmin", str(args.delmin), "-delmax", str(args.delmax), "-insmin", str(args.insmin), "-insmax", str(args.insmax),
               "-mode", args.mode, "-p_del_region", str(args.p_del_region), "-p_ins_region", str(args.p_ins_region),
               "-region_bed_url", args.region_bed_url,"-block_region_bed_url", str(args.block_region_bed_url), "-hg38", str(args.hg38)]
        if args.sum:
            cmd.append("-sum")
            
    elif args.cores is not None and args.cores > 1:
        script_name = "uniform_parallel.py"
        cmd = ["python", os.path.join(main_dir, script_name), "-ref", args.ref, "-save", args.save, "-seed", str(args.seed), 
               "-times", str(args.times), "-rep", str(args.rep), "-sv_trans", str(args.sv_trans), "-sv_inver", str(args.sv_inver), "-sv_dup", 
               str(args.sv_dup), "-sv_del", str(args.sv_del), "-sv_ins", str(args.sv_ins), "-cores", str(args.cores), "-len_bins", 
               str(args.len_bins), "-delmin", str(args.delmin), "-delmax", str(args.delmax), "-insmin", str(args.insmin), "-insmax", str(args.insmax), 
               "-block_region_bed_url", str(args.block_region_bed_url)]
    else:
        script_name = "uniform.py"
        cmd = ["python", os.path.join(main_dir, script_name), "-ref", args.ref, "-save", args.save, 
               "-seed", str(args.seed), "-times", str(args.times), "-rep", str(args.rep), "-sv_trans", str(args.sv_trans), 
               "-sv_inver", str(args.sv_inver), "-sv_dup", str(args.sv_dup), "-sv_del", str(args.sv_del), "-sv_ins", str(args.sv_ins), 
               "-snp", str(args.snp), "-snv_del", str(args.snv_del), "-snv_ins", str(args.snv_ins), "-block_region_bed_url", str(args.block_region_bed_url)]
            
    if args.snp is not None:
        cmd.extend(["-snp", str(args.snp)])
    if args.snv_del is not None:
        cmd.extend(["-snv_del", str(args.snv_del)])
    if args.snv_ins is not None:
        cmd.extend(["-snv_ins", str(args.snv_ins)])
    if args.notblockN:
        cmd.append("-notblockN")
    if args.write:
        cmd.append("-write")
    subprocess.run(cmd)

if __name__ == "__main__":
    main()

