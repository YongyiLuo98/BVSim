# any input reference
# complex SV
# dataframe SV_table
#!/usr/bin/env python
# coding: utf-8
# Complex SV
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

from scipy.stats import pareto
import time
from datetime import timedelta
# Import necessary libraries
from pandas import concat

import argparse
import math

start_time1 = time.time()
print('Complex Structure Variation Mode')

def check_start_end(row):
    if row['start'] >= row['end']:
        print('Warning: The "start" value of the .bed file is greater than or equal to the "end" value.')
        return False
    return True

# 定义一个函数，该函数接收一行数据，返回该行'start'和'end'区域中的所有点.不包含end
def get_points(row):
    return set(range(row['start'], row['end']))

def DNA_complement(sequence):
    trantab = str.maketrans('ATCGatcg','TAGCtagc')
    string = sequence.translate(trantab)
    return string

base_list=["A","T","C","G"]
#! default input
ins_selection = ['a', 't', 'c', 'g']
mis_selection = ['a','t','c','g']

base_list=["A","T","C","G"]

pdel_SV = 0.01
#len_SV_del = [50,100,500]
# len_SV_trans = [50,100,500]
# len_SV_inver = [50,100,500]
#len_SV_ins = [50,100,500]
#condition_dist_sv_del = [4/5,1/10,1/10]
# condition_dist_sv_trans = [4/5,1/10,1/10]
# condition_dist_sv_inver = [4/5,1/10,1/10]
#condition_dist_sv_ins = [4/5,1/10,1/10]
number_sv_del = [4/5,1/10,1/10]
number_sv_trans = [4/5,1/10,1/10]
number_sv_inver = [4/5,1/10,1/10]
number_sv_dup = [4/5,1/10,1/10]
copied_base_sv_prob = [4/5,1/10,1/10]
copied_base_sv_base = [50,100,500]

#deletion probability
del_sv_createria = 8.9/9
trans_sv_createria = 8.9/9
inver_sv_createria = 8.9/9
ins_location_p = 1/90000
dup_location_p = 1/90000
# len_seg_refine = len(real_con1)


ratio_b_trans = 0.5
ratio_re_dup = 0.5
ratio_re_InvDup = 0.5

#num_con=1
#mutation_rate_list = 1
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
ratio_re_dup = 0.6
#! Translocation
def translocation(SV_table, VCF_table, unblock_region_sv, SV_loop, CSV_loop, VCF_loop, Ins_dic_sv, tem_seq_post, sv_trans, condition_dist_sv_trans, len_SV_trans, times, ll_c, len_seg_refine, ratio_b_trans, chr_id):
    print('Trans:'+str(sv_trans))
    for s in range(0, sv_trans):
        ### 循环每次转位
        ### 长转位的删除部分：切割 A
        r_s = sample(unblock_region_sv, 1)[0]
        l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_trans)
        l_s_index = list(l_s_vec).index(1)
        l_s = len_SV_trans[int(l_s_index)]
        # 如果长转位超过终端或先前的删除，并且此删除重叠，则重新采样起点和长度
        circular_count_trans=0
        # 指示器：如果重新采样超过50次
        circular_count_trans_break = 0
        while (r_s+l_s > len_seg_refine) or (r_s+l_s-1 not in unblock_region_sv) or ((unblock_region_sv.index(r_s+l_s-1)-unblock_region_sv.index(r_s)<l_s-1)):
            r_s = sample(unblock_region_sv, 1)[0]
            l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_trans)
            l_s_index = list(l_s_vec).index(1)
            l_s = len_SV_trans[int(l_s_index)]
            # 计算重新采样的次数
            circular_count_trans = circular_count_trans + 1
            if circular_count_trans>times:
                print("Warning: No."+str(s)+ " translocation sampling times exceeds "+str(times)+' times. Try: reduce number of variations, increase times or reference length.')
                circular_count_trans_break = 1
                break
            else:
                circular_count_trans_break = 0
        if not circular_count_trans_break:
            circular_count_trans=0
            ## 选择 A (切割)后更新索引集
            # 对于下一个删除，不能是右或左邻居
            unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s+1)))
            #! 收集切割位置
            #location_del_trans.append(r_s)
            p_rand_trans_b = np.random.rand(1)
            ### rato_b_trans：平衡转位的概率
            # 不平衡转位：切割 A 并粘贴到 B

            if p_rand_trans_b>ratio_b_trans:### 原始位置也被删除（切割 A）
                ### 选择一个可能的插入位点：B 的插入位点
                ins_trans_loc = sample(unblock_region_sv,1)[0]
                # 插入的片段是之前采样的切割 A 区域
                # 插入的位置是我们为不平衡转位采样的位置
                # 将区域 A 粘贴到位置 B (ins_trans_loc)
                Ins_dic_sv[ins_trans_loc] = copy.deepcopy(''.join(tem_seq_post[r_s:(r_s+l_s)]))
                # 删除和插入的片段不能是左，右邻居和重叠。
                # 删除区域 B 周围的位置
                unblock_region_sv = list(set(unblock_region_sv) -set(range(ins_trans_loc,ins_trans_loc+2)))
                if ins_trans_loc in unblock_region_sv:
                    unblock_region_sv.remove(ins_trans_loc)
                
                # 收集插入位置
                #location_insert_trans.append(ins_trans_loc)
                # if ins_trans_loc in other_sites:
                #     other_sites.remove(ins_trans_loc)
                # 写 SV 表
                SV_table.loc[SV_loop] = [SV_loop,ll_c,'Translocation',r_s,r_s+l_s-1,l_s,ins_trans_loc,ins_trans_loc,l_s,0,0,0,0,0,'.']### 切割和粘贴
                SV_loop = SV_loop + 1
                
                VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s:r_s+l_s]) + tem_seq_post[r_s+l_s], tem_seq_post[r_s+l_s], '.', 'PASS', 'SVTYPE=DEL;SVLEN='+str(l_s)+';CSV_TYPE=unbalanedTrans;CSV_INDEX='+str(CSV_loop)]
                VCF_loop = VCF_loop+ 1

                VCF_table.loc[VCF_loop] = [str(chr_id), str(ins_trans_loc), 'rs' + str(VCF_loop), tem_seq_post[ins_trans_loc], tem_seq_post[ins_trans_loc] + ''.join(tem_seq_post[r_s:(r_s+l_s)]), '.', 'PASS', 'SVTYPE=INS;SVLEN='+str(l_s)+';CSV_TYPE=unbalanedTrans;CSV_INDEX='+str(CSV_loop)]
                VCF_loop = VCF_loop + 1
                
                CSV_loop = CSV_loop + 1
                # 复制 [r_s,r_s+l_s] 的片段并替换 [ins_trans_loc,ins_trans_loc+l_s] 的片段
                tem_seq_post[r_s:(r_s+l_s)] = '-'*l_s
            else:
                ### 平衡转位：切割 A 并粘贴 B，切割 B 并粘贴 A
                ### 选择 B 的删除起点
                r_s2 = sample(unblock_region_sv,1)[0]
                # 采样切割 B 的长度
                l_s_vec2 = np.random.multinomial(n=1, pvals=condition_dist_sv_trans)
                l_s_index2 = list(l_s_vec2).index(1)
                # 转位区域 B 的长度
                l_s2 = len_SV_trans[int(l_s_index2)]
                # 如果长转位超过终端或先前的删除，并且 B 重叠，则重新采样起点和长度
                circular_count_trans=0
                while (r_s2+l_s2>len_seg_refine) or (r_s2+l_s2-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s2+l_s2-1)-unblock_region_sv.index(r_s2)<l_s2-1):
                    r_s2 = sample(unblock_region_sv,1)[0]
                    l_s_vec2 = np.random.multinomial(
                        n=1, pvals=condition_dist_sv_trans)
                    l_s_index2 = list(l_s_vec2).index(1)
                    l_s2 = len_SV_trans[int(l_s_index2)]
                    # 计算重新采样的次数
                    circular_count_trans = circular_count_trans + 1
                    if circular_count_trans>times:
                        circular_count_trans_break = 1
                        print('Warning: No.'+str(s)+ " balanced translocation sampling times exceeds "+str(times)+' times. Try: reduce number of variations, increase times or reference length.')
                        break
                    else:
                        circular_count_trans_break = 0
                if not circular_count_trans_break:
                    circular_count_trans = 0
                    # 插入的位置比删除的基础小一，用于平衡转位
                    # 将区域 B 粘贴到位置 A (r_s-1)
                    Ins_dic_sv[r_s-1] = copy.deepcopy(''.join(tem_seq_post[r_s2:(r_s2+l_s2)]))
                    # 将区域 A 粘贴到位置 B (r_s2-1)
                    Ins_dic_sv[r_s2-1] = copy.deepcopy(''.join(tem_seq_post[r_s:(r_s+l_s)]))
                    # 后面的删除没有左位置，间隔和右位置
                    unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s2-1,r_s2+l_s2+1)))
                    #! 更新片段
                    #del_sites2=list(set(range(r_s2-1,r_s2+l_s2+1)))
                    # for index_seg2 in range(number_seg):
                    #     segments[index_seg2]=list(set(segments[index_seg2])-set(del_sites2))
                    
                    # 收集插入位置
                    #location_insert_trans.append(r_s-1)
                    #location_insert_trans.append(r_s2-1)
                    ins_trans_loc1 = r_s2-1
                    inserted_string1 = copy.deepcopy(''.join(tem_seq_post[r_s:(r_s+l_s)]))
                    
                    ins_trans_loc2 = r_s-1
                    inserted_string2 = copy.deepcopy(''.join(tem_seq_post[r_s2:(r_s2+l_s2)]))
                    # 写 SV 表
                    SV_table.loc[SV_loop] = [SV_loop,ll_c,'Translocation',r_s,r_s+l_s-1,l_s,r_s2,r_s2+l_s2-1,l_s2,1,0,0,0,0,'.']### 复制和粘贴
                    SV_loop = SV_loop + 1
                    
                    VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s:r_s+l_s]) + tem_seq_post[r_s+l_s], tem_seq_post[r_s+l_s], '.', 'PASS', 'SVTYPE=DEL;SVLEN='+str(l_s)+';CSV_TYPE=balanedTrans;CSV_INDEX='+str(CSV_loop)]
                    VCF_loop = VCF_loop+ 1
                    
                    VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s2:r_s2+l_s2]) + tem_seq_post[r_s2+l_s2], tem_seq_post[r_s2+l_s2], '.', 'PASS', 'SVTYPE=DEL;SVLEN='+str(l_s2)+';CSV_TYPE=balanedTrans;CSV_INDEX='+str(CSV_loop)]
                    VCF_loop = VCF_loop+ 1

                    VCF_table.loc[VCF_loop] = [str(chr_id), str(ins_trans_loc1), 'rs' + str(VCF_loop), tem_seq_post[ins_trans_loc1], tem_seq_post[ins_trans_loc1] + inserted_string1, '.', 'PASS', 'SVTYPE=INS;SVLEN='+str(l_s)+';CSV_TYPE=balanedTrans;CSV_INDEX='+str(CSV_loop)]
                    VCF_loop = VCF_loop + 1
                    
                    VCF_table.loc[VCF_loop] = [str(chr_id), str(ins_trans_loc2), 'rs' + str(VCF_loop), tem_seq_post[ins_trans_loc2], tem_seq_post[ins_trans_loc2] + inserted_string2, '.', 'PASS', 'SVTYPE=INS;SVLEN='+str(l_s2)+';CSV_TYPE=balanedTrans;CSV_INDEX='+str(CSV_loop)]
                    VCF_loop = VCF_loop + 1
                    
                    CSV_loop = CSV_loop + 1
                    
                    ##location_del_trans.append(r_s2)
                    # 复制 [r_s,r_s+l_s] 的片段并交换 [r_s2,r_s2+l_s2] 的片段
                    tem_seq_post[r_s2:(r_s2+l_s2)] = '-'*l_s2
                    tem_seq_post[r_s:(r_s+l_s)] = '-'*l_s
                    
                else:
                    break
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post

#! Tandem Duplication

# 均值和方差
# mu_Tandup, sigma_Tandup = 1000, 100
# # 生成高斯分布的随机数
# gaussian_Tandup = np.random.normal(mu_Tandup, sigma_Tandup, 1000)
# # 离散化和分bin
# gaussian_Tandup = np.round(gaussian_Tandup / 100) * 100
# # 取绝对值并向上取整
# gaussian_Tandup = [math.ceil(abs(x)) for x in gaussian_Tandup]

# def tandem_duplication(SV_table, VCF_table,True_dup_number, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post, gaussian_Tandup, times, ll_c, chr_id):
def tandem_duplication(SV_table, VCF_table,True_dup_number, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post, dup_range, times, ll_c, chr_id):
    print('DUP:'+str(True_dup_number))
    all_selected_dup_SV= sample(unblock_region_sv,True_dup_number)
    #set for resample
    other_sites = list(set(unblock_region_sv)-set(all_selected_dup_SV))
    for remain_index2 in all_selected_dup_SV:
        
        #length of copied range
        # tem_copied_base = np.random.choice(gaussian_Tandup)
        tem_copied_base = np.random.choice(dup_range)
        circular_count_dup = 0
        circular_count_dup_break = 0
        # resample if there is any other variation in this range
        # if all bases that are copied are in unblock_region_sv, then no resampling is needed
        # number of the common elements of two sets
        while len(set(range(remain_index2-tem_copied_base+1,remain_index2+1))&set(unblock_region_sv)) < tem_copied_base:
            #sites that do note have general insertions
            remain_index2 = sample(other_sites,1)[0]
            #length of copied range
            # tem_copied_base = np.random.choice(gaussian_Tandup)
            tem_copied_base = np.random.choice(dup_range)
            circular_count_dup = circular_count_dup + 1
            if circular_count_dup>times:
                circular_count_dup_break = 1
                print(" warning: Tandem duplication sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_dup_break = 0
        if not circular_count_dup_break:
            circular_count_dup=0
            if remain_index2 in other_sites:
                other_sites.remove(remain_index2)
            if remain_index2 in unblock_region_sv:
                unblock_region_sv.remove(remain_index2)
            if remain_index2+1 in unblock_region_sv:
                unblock_region_sv.remove(remain_index2+1)
            
            Ins_dic_sv[remain_index2] = ''.join(tem_seq_post[remain_index2-tem_copied_base+1:remain_index2+1])
            ins_len_dup = len(Ins_dic_sv[remain_index2])
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Duplication',remain_index2-tem_copied_base+1,remain_index2,tem_copied_base,remain_index2,remain_index2,ins_len_dup,-1,0,0,0,0,'.']
            SV_loop = SV_loop + 1
            VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index2), 'rs' + str(VCF_loop), tem_seq_post[remain_index2], tem_seq_post[remain_index2]+''.join(tem_seq_post[remain_index2-tem_copied_base+1:remain_index2+1]), '.', 'PASS', 'SVTYPE=DUP;SVLEN='+str(tem_copied_base)]
            VCF_loop= VCF_loop + 1
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post

#! Canonical Inversion
#! Inversion

# 均值和方差
# mu_inv, sigma_inv = 1000, 100
# # 生成高斯分布的随机数
# gaussian_inv = np.random.normal(mu_inv, sigma_inv, 1000)
# # 离散化和分bin
# gaussian_inv = np.round(gaussian_inv / 100) * 100
# # 取绝对值并向上取整
# gaussian_inv = [math.ceil(abs(x)) for x in gaussian_inv]

#def inversion_process(sv_inver, SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post, chr_id, len_seg_refine,ll_c,times,real_con1, gaussian_inv):
def inversion_process(sv_inver, SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post, chr_id, len_seg_refine,ll_c,times,real_con1, condition_dist_sv_inver, len_SV_inver):
    print('Canonical Inversion:'+str(sv_inver))
    for inver_num in range (0, sv_inver):
        #sample the start and length of inversion region
        r_s = sample(unblock_region_sv, 1)[0]
        #l_s = np.random.choice(gaussian_inv)
        #l_s = np.random.poisson(lamnuma_del, 1)[0]
        l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_inver)
        l_s_index = list(l_s_vec).index(1)
        l_s = len_SV_inver[int(l_s_index)]
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
            #l_s = np.random.choice(gaussian_inv)
            l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_inver)
            l_s_index = list(l_s_vec).index(1)
            l_s = len_SV_inver[int(l_s_index)]
            # count the times of resampling
            circular_count_inver = circular_count_inver + 1
            if circular_count_inver>times:
                circular_count_inver_break = 1
                print(" warning: No."+str(inver_num)+ " inversion sampling times exceeds "+str(times)+' times. Try: reduce number of variations, increase times or reference length.')
                
                break
            else:
                circular_count_inver_break = 0

        if not circular_count_inver_break:
            circular_count_inver = 0
            unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s+1)))

            original_string = copy.deepcopy(real_con1[r_s:r_s+l_s])
            #original_string_reverse = original_string[::-1]
            original_string_reverse = DNA_complement(original_string[::-1])


            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Inversion',r_s,r_s+l_s-1,l_s,r_s,r_s+l_s-1,l_s,-1,0,0,0,0,'.']
            SV_loop = SV_loop + 1

            VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), ''.join(original_string), ''.join(original_string_reverse), '.', 'PASS', 'SVTYPE=INV;SVLEN='+str(l_s)]
            VCF_loop= VCF_loop + 1

            for ll_rever in range(r_s,r_s+l_s):
                tem_seq_post[ll_rever] = copy.deepcopy(original_string_reverse[int(ll_rever-r_s)])
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post

#! ID1: TanInvDup (Tandem Inverted Dup), 11 (73.8kb)

def ID1_TanInvDup_process(unblock_region_sv, True_TanInvDup_number, times, real_con1, SV_table, VCF_table, ll_c, chr_id, tem_seq_post, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv,gaussian_ID1):
    print("ID1 (Tandem Inverted Dup):"+str(True_TanInvDup_number))
    #inserted_site_collection = sample(unblock_region_sv,full_number_base)
    all_selected_InvDup_SV= sample(unblock_region_sv,True_TanInvDup_number)
    # 打印所有被选择的位点的信息
    other_sites = list(set(unblock_region_sv)-set(all_selected_InvDup_SV))
    #InvDuplication: copy and paste
    for remain_index2 in all_selected_InvDup_SV:##others### questions
        #print(remain_index2)
        #number of copies
        #InvDup_num_index = np.random.multinomial(n=1, pvals=number_sv_InvDup)
        #InvDup_num = int(list(InvDup_num_index).index(1))+1
        
        # length of InvDuplication
        #tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
        #tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
        tem_copied_base = np.random.choice(gaussian_ID1)
        ##if all bases that are copied are in unblock_region_sv, then no resampling is needed
        ## number of the common elements of two sets
        # if len(set(range(remain_index2-tem_copied_base+1,remain_index2+1))&set(unblock_region_sv)) < tem_copied_base:
        #     repeat_falg_tran_ins4 = 1

        circular_count_InvDup = 0
        circular_count_InvDup_break = 0
        while len(set(range(remain_index2-tem_copied_base+1,remain_index2+1))&set(unblock_region_sv)) < tem_copied_base:
            # select the possibile copied length 
            #sites that do note have general insertions

            remain_index2 = sample(other_sites,1)[0]
            #tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
            #tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
            tem_copied_base = np.random.choice(gaussian_ID1)
            circular_count_InvDup = circular_count_InvDup + 1
            if circular_count_InvDup>times:
                circular_count_InvDup_break = 1
                print("Warning: ID1 (tandem inverted duplication) sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_InvDup_break = 0

        #ratio_re_InvDup:neighbor InvDuplication (right after the copied area)
        if not circular_count_InvDup_break:
            circular_count_InvDup=0
            original_string_dup = copy.deepcopy(real_con1[remain_index2-tem_copied_base+1:remain_index2+1])
            #! original_string_dup_reverse = original_string_dup[::-1]
            original_string_dup_reverse = DNA_complement(original_string_dup[::-1])
            
            # Tandem InvDuplication
            
            #later ins cannot be sampled from this point in other_sites
            if remain_index2 in other_sites:
                other_sites.remove(remain_index2)
            # ins and del cannot neighbor
            if remain_index2 in unblock_region_sv:
                unblock_region_sv.remove(remain_index2)
            if remain_index2+1 in unblock_region_sv:
                unblock_region_sv.remove(remain_index2+1)
                
            #location_insert_dup.append(remain_index2)

            #remain_index2 is a new ins position, later cannot ins here
            if remain_index2 in unblock_region_sv:
                unblock_region_sv.remove(remain_index2)
            Ins_dic_sv[remain_index2] = ''.join(copy.deepcopy(original_string_dup_reverse))
            #Ins_dic_sv[remain_index2] = ''.join(tem_seq_post[remain_index2-tem_copied_base+1:remain_index2+1])*InvDup_num
            ins_len_InvDup = len(Ins_dic_sv[remain_index2])
            
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'TanInvDup',remain_index2-tem_copied_base+1,remain_index2,tem_copied_base,remain_index2,remain_index2,ins_len_InvDup,-1,0,0,0,0,'CSV_TYPE=ID1;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            
            #ins
            VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index2), 'rs' + str(VCF_loop), tem_seq_post[remain_index2], tem_seq_post[remain_index2]+''.join(copy.deepcopy(original_string_dup_reverse)), '.', 'PASS', 'SVTYPE=TanInvDup;SVLEN='+str(ins_len_InvDup)+';CSV_TYPE=ID1;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1
            
            CSV_loop = CSV_loop + 1
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post

#! ID2: DisInvDup (Dispersed Inverted Dup), 11 (73.8kb)

def ID2_DisInvDup_process(unblock_region_sv, True_DisInvDup_number, times, real_con1, SV_table, VCF_table, ll_c, chr_id, tem_seq_post, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv,gaussian_ID2):
    print("ID2 (Dispersed Inverted Dup):"+str(True_DisInvDup_number))
    #inserted_site_collection = sample(unblock_region_sv,full_number_base)
    all_selected_InvDup_SV= sample(unblock_region_sv,True_DisInvDup_number)
    # 打印所有被选择的位点的信息
    other_sites = list(set(unblock_region_sv)-set(all_selected_InvDup_SV))
    #InvDuplication: copy and paste
    for remain_index2 in all_selected_InvDup_SV:##others### questions
        #print(remain_index2)
        #number of copies
        #InvDup_num_index = np.random.multinomial(n=1, pvals=number_sv_InvDup)
        #InvDup_num = int(list(InvDup_num_index).index(1))+1
        
        # length of InvDuplication
        # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
        # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
        tem_copied_base = np.random.choice(gaussian_ID2)

        ##if all bases that are copied are in unblock_region_sv, then no resampling is needed
        ## number of the common elements of two sets
        # if len(set(range(remain_index2-tem_copied_base+1,remain_index2+1))&set(unblock_region_sv)) < tem_copied_base:
        #     repeat_falg_tran_ins4 = 1

        circular_count_InvDup = 0
        circular_count_InvDup_break = 0
        while len(set(range(remain_index2-tem_copied_base+1,remain_index2+1))&set(unblock_region_sv)) < tem_copied_base:
            # select the possibile copied length 
            #sites that do note have general insertions

            remain_index2 = sample(other_sites,1)[0]
            # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
            # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
            tem_copied_base = np.random.choice(gaussian_ID2)
            circular_count_InvDup = circular_count_InvDup + 1
            if circular_count_InvDup>times:
                circular_count_InvDup_break = 1
                print("Warning: ID2 (dispersed duplication) sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_InvDup_break = 0


        #ratio_re_InvDup:neighbor InvDuplication (right after the copied area)
        if not circular_count_InvDup_break:
            circular_count_InvDup=0
            original_string_dup = copy.deepcopy(real_con1[remain_index2-tem_copied_base+1:remain_index2+1])
            #! original_string_dup_reverse = original_string_dup[::-1]
            original_string_dup_reverse = DNA_complement(original_string_dup[::-1])
    
            #remote InvDuplication
            #sample a pos for pasted region
            remain_index22 = sample(other_sites,1)[0]
            
            #location_insert_dup.append(remain_index22)
            
            other_sites.remove(remain_index22)
            if remain_index22 in unblock_region_sv:
                unblock_region_sv.remove(remain_index22)
            
            if remain_index22+1 in unblock_region_sv:
                unblock_region_sv.remove(remain_index22+1)
        
            Ins_dic_sv[remain_index22] = ''.join(copy.deepcopy(original_string_dup_reverse))
            
            #Ins_dic_sv[remain_index22] = ''.join(tem_seq_post[remain_index2-tem_copied_base+1:remain_index2+1])*InvDup_num
            ins_len_InvDup = len(Ins_dic_sv[remain_index22])

            SV_table.loc[SV_loop] = [SV_loop,ll_c,'DisInvDup',remain_index2-tem_copied_base+1,remain_index2,tem_copied_base,remain_index22,remain_index22,ins_len_InvDup,-1,0,0,0,0,'CSV_TYPE=ID2;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            
            #ins
            VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index22), 'rs' + str(VCF_loop), tem_seq_post[remain_index22], tem_seq_post[remain_index22]+''.join(copy.deepcopy(original_string_dup_reverse)), '.', 'PASS', 'SVTYPE=DisInvDup;CSV_TYPE=ID2;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1
            
            CSV_loop = CSV_loop + 1
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post

#! ID3: DisDup

def ID3_disdup_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, True_DisDup_number, times, ll_c, chr_id, gaussian_ID3):
    print('ID3 (Dispersed Dup):'+str(True_DisDup_number))
    # 遍历所有的重复位点
    all_selected_disdup_SV= sample(unblock_region_sv,True_DisDup_number)
    
    other_sites = list(set(unblock_region_sv)-set(all_selected_disdup_SV))
    for remain_index2 in all_selected_disdup_SV:
    #for remain_index2 in range(True_DisDup_number):
        # 选择复制的数量
        # dup_num_index = np.random.multinomial(n=1, pvals=number_sv_dup)
        # dup_num = int(list(dup_num_index).index(1))+1
        # 选择复制的长度
        # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
        # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
        tem_copied_base = np.random.choice(gaussian_ID3)
        circular_count_dup = 0
        circular_count_dup_break = 0
        # 如果所有被复制的基因位点都在 unblock_region_sv 中，则不需要重新采样
        while len(set(range(remain_index2-tem_copied_base+1,remain_index2+1))&set(unblock_region_sv)) < tem_copied_base:
            # 选择可能被复制的长度
            remain_index2 = sample(other_sites,1)[0]
            # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
            # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
            tem_copied_base = np.random.choice(gaussian_ID3)
            circular_count_dup = circular_count_dup + 1
            if circular_count_dup>times:
                circular_count_dup_break = 1
                print("Warning: ID3 (dispersed duplication) sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_dup_break = 0

        if not circular_count_dup_break:
            circular_count_dup=0
            remain_index22 = sample(other_sites,1)[0]
            
            other_sites.remove(remain_index22)
            if remain_index22 in unblock_region_sv:
                unblock_region_sv.remove(remain_index22)
            if remain_index22+1 in unblock_region_sv:
                unblock_region_sv.remove(remain_index22+1)
            Ins_dic_sv[remain_index22] = ''.join(tem_seq_post[remain_index2-tem_copied_base+1:remain_index2+1])
            ins_len_dup = len(Ins_dic_sv[remain_index22])
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'DisDup',remain_index2-tem_copied_base+1,remain_index2,tem_copied_base,remain_index22,remain_index22,ins_len_dup,-1,0,0,0,0,'CSV_TYPE=ID3;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            
            VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index22), 'rs' + str(VCF_loop), tem_seq_post[remain_index22], tem_seq_post[remain_index22]+''.join(tem_seq_post[remain_index2-tem_copied_base+1:remain_index2+1]), '.', 'PASS', 'SVTYPE=DisDup;CSV_TYPE=ID3;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1
            CSV_loop = CSV_loop + 1
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv

#! ID4: DEL+INV
## long Delinvsion we need get rid of the condition that overlapped the region of inserted trans
#p_rand_Delinv_sv = np.random.rand(1)

len_SV_Delinv = [100,500,600]
condition_dist_sv_Delinv = [4/5,1/10,1/10]

def ID4_delinv_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_Delinv_sv, gaussian_ID4, times, ll_c, chr_id, len_seg_refine, real_con1):
    print('ID4 (DelInv+InvDel):'+str(num_Delinv_sv))
    for Delinv_num in range (0, num_Delinv_sv):
        # 采样 Delinvsion 区域的起点和长度
        r_s = sample(unblock_region_sv, 1)[0]
        # l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_Delinv)
        # l_s_index = list(l_s_vec).index(1)
        # l_s = len_SV_Delinv[int(l_s_index)]
        l_s = np.random.choice(gaussian_ID4)
        # 计算我们重新采样的次数
        circular_count_Delinv = 0
        # 指示器：如果重新采样超过50次
        circular_count_Delinv_break = 0
        ### 需要重新采样的情况
        ### 需要重新采样的情况
        ### 情况 1：如果 Delinvsion 超过 real_con1 的末端
        ### 情况 2：如果 Delinvsion 的末端不在 unblock_region_sv 中
        ### 情况 3：先前的删除、插入（转位）和 Delinvsions 的重叠
        while (r_s+l_s>len_seg_refine) or (r_s+l_s-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s+l_s-1)-unblock_region_sv.index(r_s)<l_s-1):
            # 选择可能的删除点
            r_s = sample(unblock_region_sv, 1)[0]
            # l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_Delinv)
            # l_s_index = list(l_s_vec).index(1)
            # l_s = len_SV_Delinv[int(l_s_index)]
            l_s = np.random.choice(gaussian_ID4)
            # 计算重新采样的次数
            circular_count_Delinv = circular_count_Delinv + 1
            if circular_count_Delinv>times:
                circular_count_Delinv_break = 1
                print("Warning: ID4 (DEL+INV) sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_Delinv_break = 0
        if not circular_count_Delinv_break:
            circular_count_Delinv = 0
            unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s+1)))
            cut_point = sample(range(int(l_s/4),int(l_s*3/4)),1)[0]
            
            # 决定 InvDuplication 的类型
            p_delinv_random_num = np.random.rand(1)
            
            if not circular_count_Delinv_break:
                circular_count_Delinv=0
                if p_delinv_random_num< 0.5:
                    original_string = copy.deepcopy(real_con1[r_s+cut_point:r_s+l_s])
                    original_string_reverse = DNA_complement(original_string[::-1])
                    r_s_del=r_s
                    l_s_del=cut_point
                    r_s_inv=r_s+cut_point
                    l_s_inv=l_s-cut_point
                else:
                    # INV + DEL
                    original_string = copy.deepcopy(real_con1[r_s:r_s+cut_point])
                    original_string_reverse = DNA_complement(original_string[::-1])
                    r_s_del=r_s+cut_point
                    l_s_del=l_s-cut_point
                    r_s_inv=r_s
                    l_s_inv=cut_point
                SV_table.loc[SV_loop] = [SV_loop,ll_c,'Deletion',r_s_del,r_s_del+l_s_del-1,l_s_del,r_s_del,r_s_del+l_s_del-1,l_s_del,-1,0,0,0,0,'CSV_TYPE=ID4;CSV_INDEX='+str(CSV_loop)]
                SV_loop = SV_loop + 1
                SV_table.loc[SV_loop] = [SV_loop,ll_c,'Inversion',r_s_inv,r_s_inv+l_s_inv-1,l_s_inv,r_s_inv,r_s_inv+l_s_inv-1,l_s_inv,-1,0,0,0,0,'CSV_TYPE=ID4;CSV_INDEX='+str(CSV_loop)]
                SV_loop = SV_loop + 1
                VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s_inv), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s_inv:r_s_inv+l_s_inv]), ''.join(original_string_reverse), '.', 'PASS', 'SVTYPE=INV;CSV_TYPE=ID4;CSV_INDEX='+str(CSV_loop)]
                VCF_loop= VCF_loop + 1
                VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s_del), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s_del:r_s_del+l_s_del]) + tem_seq_post[r_s_del+l_s_del], tem_seq_post[r_s_del+l_s_del], '.', 'PASS', 'SVTYPE=DEL;CSV_TYPE=ID4;CSV_INDEX='+str(CSV_loop)]
                VCF_loop = VCF_loop+ 1 
                CSV_loop = CSV_loop + 1
                for ll_rever in range(r_s_inv,r_s_inv+l_s_inv):
                    tem_seq_post[ll_rever] = copy.deepcopy(original_string_reverse[int(ll_rever-r_s_inv)])
                tem_seq_post[r_s_del:(r_s_del+l_s_del)] = '-'*l_s_del
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post


#! ID5: DEL+ DisInvDup 

def ID5_disinvdupdel_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_ID5_csv, times, ll_c, chr_id, real_con1, condition_dist_sv_del, len_SV_del, len_seg_refine, gaussian_ID5):
    print("ID5 (DEL+ DisInvDup or DisInvDup + DEL):"+str(num_ID5_csv))
    all_selected_DisDup_SV= sample(unblock_region_sv,num_ID5_csv)
    other_sites = list(set(unblock_region_sv)-set(all_selected_DisDup_SV))
    # DisDuplication: copy and paste
    for remain_index2 in all_selected_DisDup_SV:##others### questions
        # copy A, site and length of DisDup
        # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
        # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
        tem_copied_base = np.random.choice(gaussian_ID5)
        circular_count_DisDup = 0
        circular_count_DisDup_break = 0
        while len(set(range(remain_index2-tem_copied_base+1,remain_index2+1))&set(unblock_region_sv)) < tem_copied_base:
            remain_index2 = sample(other_sites,1)[0]
            # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
            # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
            tem_copied_base = np.random.choice(gaussian_ID5)
            circular_count_DisDup = circular_count_DisDup + 1
            if circular_count_DisDup>times:
                circular_count_DisDup_break = 1
                print("Warning: ID5 (DEL+ DisInvDup or DisInvDup + DEL) sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_DisDup_break = 0
        if not circular_count_DisDup_break:
            circular_count_DisDup = 0
            original_string_dup = copy.deepcopy(real_con1[remain_index2-tem_copied_base+1:remain_index2+1])
            original_string_dup_reverse = DNA_complement(original_string_dup[::-1])
            # sample the deletion site
            r_s = sample(unblock_region_sv, 1)[0]
            l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_del)
            l_s_index = list(l_s_vec).index(1)
            l_s = len_SV_del[int(l_s_index)]
            circular_count_Del = 0
            circular_count_Del_break = 0
            while (r_s+l_s>len_seg_refine) or (r_s+l_s-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s+l_s-1)-unblock_region_sv.index(r_s)<l_s-1):
                r_s = sample(unblock_region_sv, 1)[0]
                l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_del)
                l_s_index = list(l_s_vec).index(1)
                l_s = len_SV_del[int(l_s_index)]
                circular_count_Del = circular_count_Del + 1
                if circular_count_Del>times:
                    circular_count_Del_break = 1
                    break
                else:
                    circular_count_Del_break = 0
            if not circular_count_Del_break:
                circular_count_Del = 0
                unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s+1)))
                
                p_left_right = np.random.rand(1)
                if p_left_right<0.5:#left
                    remain_index22 = r_s-1
                else:#right
                    remain_index22 = r_s+l_s-1
                other_sites.remove(remain_index22)
                if remain_index22 in unblock_region_sv:
                    unblock_region_sv.remove(remain_index22)
                if remain_index22 in unblock_region_sv:
                    unblock_region_sv.remove(remain_index22)
                if remain_index22+1 in unblock_region_sv:
                    unblock_region_sv.remove(remain_index22+1)
                Ins_dic_sv[remain_index22] = ''.join(copy.deepcopy(original_string_dup_reverse))
                ins_len_DisDup = len(Ins_dic_sv[remain_index22])
                SV_table.loc[SV_loop] = [SV_loop,ll_c,'Deletion',r_s,r_s+l_s-1,l_s,r_s,r_s+l_s-1,l_s,-1,0,0,0,0,'CSV_TYPE=ID5;CSV_INDEX='+str(CSV_loop)]
                SV_loop = SV_loop + 1
                SV_table.loc[SV_loop] = [SV_loop,ll_c,'DisInvDup',remain_index2-tem_copied_base+1,remain_index2,tem_copied_base,remain_index22,remain_index22,ins_len_DisDup,-1,0,0,0,0,'CSV_TYPE=ID5;CSV_INDEX='+str(CSV_loop)]
                SV_loop = SV_loop + 1
                VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s:r_s+l_s]) + tem_seq_post[r_s+l_s], tem_seq_post[r_s+l_s], '.', 'PASS', 'SVTYPE=DEL;CSV_TYPE=ID5;CSV_INDEX='+str(CSV_loop)]
                VCF_loop = VCF_loop+ 1 
                VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index22), 'rs' + str(VCF_loop), tem_seq_post[remain_index22], tem_seq_post[remain_index22]+''.join(copy.deepcopy(original_string_dup_reverse)), '.', 'PASS', 'SVTYPE=DisInvDup;CSV_TYPE=ID5;CSV_INDEX='+str(CSV_loop)]
                VCF_loop= VCF_loop + 1
                CSV_loop = CSV_loop + 1
                tem_seq_post[r_s:(r_s+l_s)] = '-'*l_s
                
            else:
                break
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post

#! ID6: DEL+ DisDup

def ID6_disdupdel_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_ID6_csv, gaussian_ID6, times, ll_c, chr_id, real_con1, condition_dist_sv_del, len_SV_del, len_seg_refine):
    print("ID6 (DEL+ DisDup or DisDup + DEL):"+str(num_ID6_csv))
    
    all_selected_DisDup_SV= sample(unblock_region_sv,num_ID6_csv)
    other_sites = list(set(unblock_region_sv)-set(all_selected_DisDup_SV))
    # DisDuplication: copy and paste
    for remain_index2 in all_selected_DisDup_SV:##others### questions
        # copy A, site and length of DisDup
        # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
        # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
        tem_copied_base = np.random.choice(gaussian_ID6)
        
        circular_count_DisDup = 0
        circular_count_DisDup_break = 0
        while len(set(range(remain_index2-tem_copied_base+1,remain_index2+1))&set(unblock_region_sv)) < tem_copied_base:
            remain_index2 = sample(other_sites,1)[0]
            # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
            # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
            tem_copied_base = np.random.choice(gaussian_ID6)
            circular_count_DisDup = circular_count_DisDup + 1
            if circular_count_DisDup>times:
                circular_count_DisDup_break = 1
                print("Warning: ID6 (DEL+ DisDup or DisDup + DEL) sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_DisDup_break = 0
        if not circular_count_DisDup_break:
            circular_count_DisDup = 0
            original_string_dup = copy.deepcopy(real_con1[remain_index2-tem_copied_base+1:remain_index2+1])
            #original_string_dup_reverse = DNA_complement(original_string_dup[::-1])
            # sample the deletion site
            r_s = sample(unblock_region_sv, 1)[0]
            l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_del)
            l_s_index = list(l_s_vec).index(1)
            l_s = len_SV_del[int(l_s_index)]
            circular_count_Del = 0
            circular_count_Del_break = 0
            while (r_s+l_s>len_seg_refine) or (r_s+l_s-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s+l_s-1)-unblock_region_sv.index(r_s)<l_s-1):
                r_s = sample(unblock_region_sv, 1)[0]
                l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_del)
                l_s_index = list(l_s_vec).index(1)
                l_s = len_SV_del[int(l_s_index)]
                circular_count_Del = circular_count_Del + 1
                if circular_count_Del>times:
                    circular_count_Del_break = 1
                    break
                else:
                    circular_count_Del_break = 0
            if not circular_count_Del_break:
                circular_count_Del = 0
                unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s+1)))
                p_left_right = np.random.rand(1)
                if p_left_right<0.5:#left
                    remain_index22 = r_s-1
                else:#right
                    remain_index22 = r_s+l_s-1
                other_sites.remove(remain_index22)
                if remain_index22 in unblock_region_sv:
                    unblock_region_sv.remove(remain_index22)
                if remain_index22 in unblock_region_sv:
                    unblock_region_sv.remove(remain_index22)
                if remain_index22+1 in unblock_region_sv:
                    unblock_region_sv.remove(remain_index22+1)
                Ins_dic_sv[remain_index22] = ''.join(copy.deepcopy(original_string_dup))
                ins_len_DisDup = len(Ins_dic_sv[remain_index22])
                SV_table.loc[SV_loop] = [SV_loop,ll_c,'Deletion',r_s,r_s+l_s-1,l_s,r_s,r_s+l_s-1,l_s,-1,0,0,0,0,'CSV_TYPE=ID6;CSV_INDEX='+str(CSV_loop)]
                SV_loop = SV_loop + 1
                SV_table.loc[SV_loop] = [SV_loop,ll_c,'DisDup',remain_index2-tem_copied_base+1,remain_index2,tem_copied_base,remain_index22,remain_index22,ins_len_DisDup,-1,0,0,0,0,'CSV_TYPE=ID6;CSV_INDEX='+str(CSV_loop)]
                SV_loop = SV_loop + 1
                VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s:r_s+l_s]) + tem_seq_post[r_s+l_s], tem_seq_post[r_s+l_s], '.', 'PASS', 'SVTYPE=DEL;CSV_TYPE=ID6;CSV_INDEX='+str(CSV_loop)]
                VCF_loop = VCF_loop+ 1 
                VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index22), 'rs' + str(VCF_loop), tem_seq_post[remain_index22], tem_seq_post[remain_index22]+original_string_dup, '.', 'PASS', 'SVTYPE=DisDup;CSV_TYPE=ID6;CSV_INDEX='+str(CSV_loop)]
                VCF_loop= VCF_loop + 1
                CSV_loop = CSV_loop + 1
                tem_seq_post[r_s:(r_s+l_s)] = '-'*l_s
            else:
                break
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post


#! ID7: TanDup+DEL

def ID7_tandupDEL_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_ID7_csv, gaussian_ID7, times, ll_c, chr_id, real_con1, condition_dist_sv_del, len_SV_del,len_seg_refine):
    print('ID7 (TanDup+DEL):'+str(num_ID7_csv))
    # TanDuplication: copy and paste
    for circule_dup in range(0,num_ID7_csv):##others### questions
        # copy A, site and length of TanDup
        remain_index2 = sample(unblock_region_sv,1)[0]
        # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
        # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
        tem_copied_base = np.random.choice(gaussian_ID7)

        circular_count_TanDup = 0
        circular_count_TanDup_break = 0
        # if all bases that are copied are in unblock_region_sv, then no resampling is needed
        while len(set(range(remain_index2-tem_copied_base+1,remain_index2+1))&set(unblock_region_sv)) < tem_copied_base:
            remain_index2 = sample(unblock_region_sv,1)[0]
            # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
            # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
            tem_copied_base = np.random.choice(gaussian_ID7)
            circular_count_TanDup = circular_count_TanDup + 1
            if circular_count_TanDup>times:
                circular_count_TanDup_break = 1
                print("Warning: ID7 (TanDup+DEL) sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_TanDup_break = 0
        if not circular_count_TanDup_break:
            circular_count_TanDup = 0
            original_string_dup = copy.deepcopy(real_con1[remain_index2-tem_copied_base+1:remain_index2+1])
            #original_string_dup_reverse = DNA_complement(original_string_dup[::-1])
            # sample the deletion site
            r_s = remain_index2+1
            l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_del)
            l_s_index = list(l_s_vec).index(1)
            l_s = len_SV_del[int(l_s_index)]
            circular_count_Del = 0
            circular_count_Del_break = 0
            # if deletion exceeds the end of real_con1 or end of deletion is not in unblock_region_sv or overlap of previous deletions, insertions (translocations) and deletions
            while (r_s+l_s>len_seg_refine) or (r_s+l_s-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s+l_s-1)-unblock_region_sv.index(r_s)<l_s-1):
                l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_del)
                l_s_index = list(l_s_vec).index(1)
                l_s = len_SV_del[int(l_s_index)]
                circular_count_Del = circular_count_Del + 1
                if circular_count_Del > times:
                    circular_count_Del_break = 1
                    print("Warning: ID7 (TanDup+DEL) sampling times exceeds "+str(times)+' times.')
                    break
                else:
                    circular_count_Del_break = 0
            if not circular_count_Del_break:
                circular_count_Del = 0
                unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s+1)))
                if remain_index2 in unblock_region_sv:
                    unblock_region_sv.remove(remain_index2)
                if remain_index2+1 in unblock_region_sv:
                    unblock_region_sv.remove(remain_index2+1)

                # TanDup + DEL
            
                Ins_dic_sv[remain_index2] = ''.join(copy.deepcopy(original_string_dup))
                ins_len_TanDup = len(Ins_dic_sv[remain_index2])
                SV_table.loc[SV_loop] = [SV_loop,ll_c,'TanDup',remain_index2-tem_copied_base+1,remain_index2,tem_copied_base,remain_index2,remain_index2,ins_len_TanDup,-1,0,0,0,0,'CSV_TYPE=ID7;CSV_INDEX='+str(CSV_loop)]
                SV_loop = SV_loop + 1
                SV_table.loc[SV_loop] = [SV_loop,ll_c,'Deletion',r_s,r_s+l_s-1,l_s,r_s,r_s+l_s-1,l_s,-1,0,0,0,0,'CSV_TYPE=ID7;CSV_INDEX='+str(CSV_loop)]
                SV_loop = SV_loop + 1
                VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index2), 'rs' + str(VCF_loop), tem_seq_post[remain_index2], tem_seq_post[remain_index2]+''.join(copy.deepcopy(original_string_dup)), '.', 'PASS', 'SVTYPE=TanDup;CSV_TYPE=ID7;CSV_INDEX='+str(CSV_loop)]
                VCF_loop= VCF_loop + 1
                VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s:r_s+l_s]) + tem_seq_post[r_s+l_s], tem_seq_post[r_s+l_s], '.', 'PASS', 'SVTYPE=DEL;CSV_TYPE=ID7;CSV_INDEX='+str(CSV_loop)]
                VCF_loop = VCF_loop+ 1 
                CSV_loop = CSV_loop + 1
                tem_seq_post[r_s:(r_s+l_s)] = '-'*l_s
            else:
                break
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post



#! ID8: TanInvDup+DEL

def ID8_taninvdupDEL_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_ID8_csv, gaussian_ID8, times, ll_c, chr_id, real_con1, condition_dist_sv_del, len_SV_del, len_seg_refine):
    print('ID8 (TanInvDup+DEL):'+str(num_ID8_csv))
    # TanDuplication: copy and paste
    for circule_dup in range(0,num_ID8_csv):##others### questions
        # copy A, site and length of TanDup
        remain_index2 = sample(unblock_region_sv,1)[0]
        # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
        # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
        tem_copied_base = np.random.choice(gaussian_ID8)

        circular_count_TanDup = 0
        circular_count_TanDup_break = 0
        # if all bases that are copied are in unblock_region_sv, then no resampling is needed
        while len(set(range(remain_index2-tem_copied_base+1,remain_index2+1))&set(unblock_region_sv)) < tem_copied_base:
            remain_index2 = sample(unblock_region_sv,1)[0]
            # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
            # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
            tem_copied_base = np.random.choice(gaussian_ID8)
            circular_count_TanDup = circular_count_TanDup + 1
            if circular_count_TanDup>times:
                circular_count_TanDup_break = 1
                print("Warning: ID8 (TanInvDup+DEL) sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_TanDup_break = 0
        if not circular_count_TanDup_break:
            circular_count_TanDup = 0
            original_string_dup = copy.deepcopy(real_con1[remain_index2-tem_copied_base+1:remain_index2+1])
            original_string_dup_reverse = DNA_complement(original_string_dup[::-1])
            # sample the deletion site
            r_s = remain_index2+1
            l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_del)
            l_s_index = list(l_s_vec).index(1)
            l_s = len_SV_del[int(l_s_index)]
            circular_count_Del = 0
            circular_count_Del_break = 0
            # if deletion exceeds the end of real_con1 or end of deletion is not in unblock_region_sv or overlap of previous deletions, insertions (translocations) and deletions
            while (r_s+l_s>len_seg_refine) or (r_s+l_s-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s+l_s-1)-unblock_region_sv.index(r_s)<l_s-1):
                l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_del)
                l_s_index = list(l_s_vec).index(1)
                l_s = len_SV_del[int(l_s_index)]
                circular_count_Del = circular_count_Del + 1
                if circular_count_Del > times:
                    circular_count_Del_break = 1
                    print("Warning: ID8 (TanInvDup+DEL) sampling times exceeds "+str(times)+' times.')
                    break
                else:
                    circular_count_Del_break = 0
            if not circular_count_Del_break:
                circular_count_Del = 0
                unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s+1)))
                if remain_index2 in unblock_region_sv:
                    unblock_region_sv.remove(remain_index2)
                if remain_index2+1 in unblock_region_sv:
                    unblock_region_sv.remove(remain_index2+1)
               
                # TanInvDup + DEL
            
                Ins_dic_sv[remain_index2] = ''.join(copy.deepcopy(original_string_dup_reverse))
                ins_len_InvDup = len(Ins_dic_sv[remain_index2])
                SV_table.loc[SV_loop] = [SV_loop,ll_c,'TanInvDup',remain_index2-tem_copied_base+1,remain_index2,tem_copied_base,remain_index2,remain_index2,ins_len_InvDup,-1,0,0,0,0,'CSV_TYPE=ID8;CSV_INDEX='+str(CSV_loop)]
                SV_loop = SV_loop + 1
                SV_table.loc[SV_loop] = [SV_loop,ll_c,'Deletion',r_s,r_s+l_s-1,l_s,r_s,r_s+l_s-1,l_s,-1,0,0,0,0,'CSV_TYPE=ID8;CSV_INDEX='+str(CSV_loop)]
                SV_loop = SV_loop + 1
                VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index2), 'rs' + str(VCF_loop), tem_seq_post[remain_index2], tem_seq_post[remain_index2]+''.join(copy.deepcopy(original_string_dup_reverse)), '.', 'PASS', 'SVTYPE=TanInvDup;CSV_TYPE=ID8;CSV_INDEX='+str(CSV_loop)]
                VCF_loop= VCF_loop + 1
                VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s:r_s+l_s]) + tem_seq_post[r_s+l_s], tem_seq_post[r_s+l_s], '.', 'PASS', 'SVTYPE=DEL;CSV_TYPE=ID8;CSV_INDEX='+str(CSV_loop)]
                VCF_loop = VCF_loop+ 1 
                CSV_loop = CSV_loop + 1
                tem_seq_post[r_s:(r_s+l_s)] = '-'*l_s
                
            else:
                break
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post

#! ID9: TanDup + DEL + INV

def ID9_tandup_delinv_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_ID9_csv, gaussian_ID9, times, ll_c, chr_id, real_con1, condition_dist_sv_del, len_SV_del, len_seg_refine):
    # TanDuplication: copy and paste
    print('ID9 (TanDup + DEL + INV):'+str(num_ID9_csv))
    for circule_dup in range(0,num_ID9_csv):##others### questions
        # copy A, site and length of TanDup
        remain_index2 = sample(unblock_region_sv,1)[0]
        # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
        # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
        tem_copied_base = np.random.choice(gaussian_ID9)
        
        circular_count_TanDup = 0
        circular_count_TanDup_break = 0
        while len(set(range(remain_index2-tem_copied_base+1,remain_index2+1))&set(unblock_region_sv)) < tem_copied_base:
            remain_index2 = sample(unblock_region_sv,1)[0]
            # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
            # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
            tem_copied_base = np.random.choice(gaussian_ID9)
            circular_count_TanDup = circular_count_TanDup + 1
            if circular_count_TanDup>times:
                circular_count_TanDup_break = 1
                print("Warning: ID9 (TanDup + DEL + INV) sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_TanDup_break = 0
        if not circular_count_TanDup_break:
            circular_count_TanDup = 0
            original_string_dup = copy.deepcopy(real_con1[remain_index2-tem_copied_base+1:remain_index2+1])
            #! original_string_dup_reverse = original_string_dup[::-1]
            #original_string_dup_reverse = DNA_complement(original_string_dup[::-1])
            #sample del site
            r_s = remain_index2+1
            l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_del)
            l_s_index = list(l_s_vec).index(1)
            l_s = len_SV_del[int(l_s_index)]
            circular_count_Del = 0
            circular_count_Del_break = 0
            while (r_s+l_s>len_seg_refine) or (r_s+l_s-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s+l_s-1)-unblock_region_sv.index(r_s)<l_s-1):

                l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_del)
                l_s_index = list(l_s_vec).index(1)
                l_s = len_SV_del[int(l_s_index)]
                # count the times of resampling
                circular_count_Del = circular_count_Del + 1
                if circular_count_Del > times:
                    circular_count_Del_break = 1
                    print("Warning: ID9 (TanDup + DEL + INV) sampling times exceeds "+str(times)+' times.')
                    break
                else:
                    circular_count_Del_break = 0
            if not circular_count_Del_break:
                circular_count_Del = 0
                unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s+1)))
            
                cut_point = sample(range(int(l_s/4),int(l_s*3/4)),1)[0]
                
                #to decide type of InvDuplication
                p_delinv_random_num = np.random.rand(1)
                if p_delinv_random_num< 0.5:

                    
                    original_string = copy.deepcopy(real_con1[r_s+cut_point:r_s+l_s])
                    #! original_string_reverse = original_string[::-1]
                    original_string_reverse = DNA_complement(original_string[::-1])
                    
                    r_s_del=r_s
                    l_s_del=cut_point
                    r_s_inv=r_s+cut_point
                    l_s_inv=l_s-cut_point

                else:
                    #INV + DEL
                    original_string = copy.deepcopy(real_con1[r_s:r_s+cut_point])
                    #! original_string_reverse = original_string[::-1]
                    original_string_reverse = DNA_complement(original_string[::-1])

                    
                    r_s_del=r_s+cut_point
                    l_s_del=l_s-cut_point
                    r_s_inv=r_s
                    l_s_inv=cut_point

            else:
                break
                
            if remain_index2 in unblock_region_sv:
                unblock_region_sv.remove(remain_index2)
            if remain_index2 in unblock_region_sv:
                unblock_region_sv.remove(remain_index2)
            if remain_index2+1 in unblock_region_sv:
                unblock_region_sv.remove(remain_index2+1)
            
            # TanDup + DEL
            
            Ins_dic_sv[remain_index2] = ''.join(copy.deepcopy(original_string_dup))
            ins_len_TanDup = len(Ins_dic_sv[remain_index2])

            SV_table.loc[SV_loop] = [SV_loop,ll_c,'TanDup',remain_index2-tem_copied_base+1,remain_index2,tem_copied_base,remain_index2,remain_index2,ins_len_TanDup,-1,0,0,0,0,'CSV_TYPE=ID9;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Deletion',r_s_del,r_s_del+l_s_del-1,l_s_del,r_s_del,r_s_del+l_s_del-1,l_s_del,-1,0,0,0,0,'CSV_TYPE=ID9;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Inversion',r_s_inv,r_s_inv+l_s_inv-1,l_s_inv,r_s_inv,r_s_inv+l_s_inv-1,l_s_inv,-1,0,0,0,0,'CSV_TYPE=ID9;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            
            VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index2), 'rs' + str(VCF_loop), tem_seq_post[remain_index2], tem_seq_post[remain_index2]+''.join(copy.deepcopy(original_string_dup)), '.', 'PASS', 'SVTYPE=TanDup;CSV_TYPE=ID9;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1
            
            VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s_inv), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s_inv:r_s_inv+l_s_inv]), ''.join(original_string_reverse), '.', 'PASS', 'SVTYPE=INV;CSV_TYPE=ID9;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1

            VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s_del), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s_del:r_s_del+l_s_del]) + tem_seq_post[r_s_del+l_s_del], tem_seq_post[r_s_del+l_s_del], '.', 'PASS', 'SVTYPE=DEL;CSV_TYPE=ID9;CSV_INDEX='+str(CSV_loop)]
            VCF_loop = VCF_loop+ 1 
            
            CSV_loop = CSV_loop + 1
            
            for ll_rever in range(r_s_inv,r_s_inv+l_s_inv):
                    tem_seq_post[ll_rever] = copy.deepcopy(original_string_reverse[int(ll_rever-r_s_inv)])
                
            tem_seq_post[r_s_del:(r_s_del+l_s_del)] = '-'*l_s_del
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post



#! ID10: TanInvDup + DEL + INV
    
def ID10_taninvdup_delinv_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_ID10_csv, gaussian_ID10, times, ll_c, chr_id, real_con1, condition_dist_sv_del, len_SV_del, len_seg_refine):
    # TanDuplication: copy and paste
    print('ID10 (TanInvDup + DEL + INV):'+str(num_ID10_csv))
    for circule_dup in range(0,num_ID10_csv):##others### questions
        # copy A, site and length of TanDup
        remain_index2 = sample(unblock_region_sv,1)[0]
        # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
        # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
        tem_copied_base = np.random.choice(gaussian_ID10)
        circular_count_TanDup = 0
        circular_count_TanDup_break = 0
        while len(set(range(remain_index2-tem_copied_base+1,remain_index2+1))&set(unblock_region_sv)) < tem_copied_base:
            remain_index2 = sample(unblock_region_sv,1)[0]
            # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[0]
            # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
            tem_copied_base = np.random.choice(gaussian_ID10)
            circular_count_TanDup = circular_count_TanDup + 1
            if circular_count_TanDup>times:
                circular_count_TanDup_break = 1
                print("Warning: ID10 (TanInvDup + DEL + INV) sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_TanDup_break = 0
        if not circular_count_TanDup_break:
            circular_count_TanDup = 0
            original_string_dup = copy.deepcopy(real_con1[remain_index2-tem_copied_base+1:remain_index2+1])
            #! original_string_dup_reverse = original_string_dup[::-1]
            original_string_dup_reverse = DNA_complement(original_string_dup[::-1])
            #sample del site
            r_s = remain_index2+1
            l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_del)
            l_s_index = list(l_s_vec).index(1)
            l_s = len_SV_del[int(l_s_index)]
            circular_count_Del = 0
            circular_count_Del_break = 0
            while (r_s+l_s>len_seg_refine) or (r_s+l_s-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s+l_s-1)-unblock_region_sv.index(r_s)<l_s-1):

                l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_del)
                l_s_index = list(l_s_vec).index(1)
                l_s = len_SV_del[int(l_s_index)]
                # count the times of resampling
                circular_count_Del = circular_count_Del + 1
                if circular_count_Del > times:
                    circular_count_Del_break = 1
                    print("Warning: ID10 (TanInvDup + DEL + INV) sampling times exceeds "+str(times)+' times.')
                    break
                else:
                    circular_count_Del_break = 0
            if not circular_count_Del_break:
                circular_count_Del = 0
                unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s+1)))
            
                cut_point = sample(range(int(l_s/4),int(l_s*3/4)),1)[0]
                
                #to decide type of InvDuplication
                p_delinv_random_num = np.random.rand(1)
                if p_delinv_random_num< 0.5:

                    original_string = copy.deepcopy(real_con1[r_s+cut_point:r_s+l_s])
                    #! original_string_reverse = original_string[::-1]
                    original_string_reverse = DNA_complement(original_string[::-1])
                    
                    r_s_del=r_s
                    l_s_del=cut_point
                    r_s_inv=r_s+cut_point
                    l_s_inv=l_s-cut_point

                else:
                    #INV + DEL
                    original_string = copy.deepcopy(real_con1[r_s:r_s+cut_point])
                    #! original_string_reverse = original_string[::-1]
                    original_string_reverse = DNA_complement(original_string[::-1])

                    
                    r_s_del=r_s+cut_point
                    l_s_del=l_s-cut_point
                    r_s_inv=r_s
                    l_s_inv=cut_point

            else:
                break
                
            if remain_index2 in unblock_region_sv:
                unblock_region_sv.remove(remain_index2)
            if remain_index2 in unblock_region_sv:
                unblock_region_sv.remove(remain_index2)
            if remain_index2+1 in unblock_region_sv:
                unblock_region_sv.remove(remain_index2+1)

            # TanInvDup + DEL
            
            Ins_dic_sv[remain_index2] = ''.join(copy.deepcopy(original_string_dup_reverse))
            
            ins_len_InvDup = len(Ins_dic_sv[remain_index2])
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'TanInvDup',remain_index2-tem_copied_base+1,remain_index2,tem_copied_base,remain_index2,remain_index2,ins_len_InvDup,-1,0,0,0,0,'CSV_TYPE=ID10;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Deletion',r_s_del,r_s_del+l_s_del-1,l_s_del,r_s_del,r_s_del+l_s_del-1,l_s_del,-1,0,0,0,0,'CSV_TYPE=ID10;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Inversion',r_s_inv,r_s_inv+l_s_inv-1,l_s_inv,r_s_inv,r_s_inv+l_s_inv-1,l_s_inv,-1,0,0,0,0,'CSV_TYPE=ID10;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            
            #ins
            VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index2), 'rs' + str(VCF_loop), tem_seq_post[remain_index2], tem_seq_post[remain_index2]+''.join(copy.deepcopy(original_string_dup_reverse)), '.', 'PASS', 'SVTYPE=TanInvDup;CSV_TYPE=ID10;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1
            
            VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s_inv), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s_inv:r_s_inv+l_s_inv]), ''.join(original_string_reverse), '.', 'PASS', 'SVTYPE=INV;CSV_TYPE=ID10;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1

            VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s_del), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s_del:r_s_del+l_s_del]) + tem_seq_post[r_s_del+l_s_del], tem_seq_post[r_s_del+l_s_del], '.', 'PASS', 'SVTYPE=DEL;CSV_TYPE=ID10;CSV_INDEX='+str(CSV_loop)]
            VCF_loop = VCF_loop+ 1 
            
            CSV_loop = CSV_loop + 1
            
            for ll_rever in range(r_s_inv,r_s_inv+l_s_inv):
                    tem_seq_post[ll_rever] = copy.deepcopy(original_string_reverse[int(ll_rever-r_s_inv)])
                
            tem_seq_post[r_s_del:(r_s_del+l_s_del)] = '-'*l_s_del

        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post

#! ID11: paired-Deletion Inversion (delInvdel)


def ID11_delinvdel_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_delinvdel_csv,gaussian_ID11, times, ll_c, chr_id, len_seg_refine):
    print('ID11 (DelInvDel):'+str(num_delinvdel_csv))
    for Delinv_num in range (0, num_delinvdel_csv):
        # 采样 Delinvsion 区域的起点和长度
        r_s = sample(unblock_region_sv, 1)[0]
        l_s = np.random.choice(gaussian_ID11)
        
        # 计算我们重新采样的次数
        circular_count_delinvdel = 0
        # 指示器：如果重新采样超过50次
        circular_count_delinvdel_break = 0
        ### 需要重新采样的情况
        ### 情况 1：如果 Delinvsion 超过 real_con1 的末端
        ### 情况 2：如果 Delinvsion 的末端不在 unblock_region_sv 中
        ### 情况 3：先前的删除、插入（转位）和 Delinvsions 的重叠
        while (r_s+l_s>len_seg_refine) or (r_s+l_s-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s+l_s-1)-unblock_region_sv.index(r_s)<l_s-1):
            # 选择可能的删除点
            r_s = sample(unblock_region_sv, 1)[0]
            l_s = np.random.choice(gaussian_ID11)
            # 计算重新采样的次数
            circular_count_delinvdel = circular_count_delinvdel + 1
            if circular_count_delinvdel>times:
                circular_count_delinvdel_break = 1
                print("Warning: ID11 (delInvdel), sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_delinvdel_break = 0
        if not circular_count_delinvdel_break:
            circular_count_delinvdel = 0
            unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s+1)))
            cut_point1 = sample(range(int(l_s/4),int(l_s*5/12)),1)[0]
            cut_point2 = sample(range(int(l_s*7/12),int(l_s*3/4)),1)[0]
            # print('r_s:r_s+l_s'+str(r_s)+','+str(l_s))
            # print(cut_point1, cut_point2)
            
            if not circular_count_delinvdel_break:
                circular_count_delinvdel=0
                unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s+1)))
                # DEL+INV + DEL
                
                r_s_del1=r_s
                
                l_s_del1=cut_point1
                
                r_s_del2=r_s + cut_point2
                
                l_s_del2=l_s-cut_point2
                
                r_s_inv = r_s + cut_point1
                
                l_s_inv=cut_point2-cut_point1
                
                # print('r_s_del1'+str(r_s_del1))
                # print('r_s_del2'+str(r_s_del2))
                # print('r_s_inv'+str(r_s_inv))

                original_string = copy.deepcopy(''.join(tem_seq_post[r_s_inv:r_s_inv+l_s_inv]))
                original_string_reverse = DNA_complement(original_string[::-1])
                
                SV_table.loc[SV_loop] = [SV_loop,ll_c,'Deletion',r_s_del1,r_s_del1+l_s_del1-1,l_s_del1,r_s_del1,r_s_del1+l_s_del1-1,l_s,-1,0,0,0,0,'CSV_TYPE=ID11;CSV_INDEX='+str(CSV_loop)]
                SV_loop = SV_loop + 1
                SV_table.loc[SV_loop] = [SV_loop,ll_c,'Inversion',r_s_inv,r_s_inv+l_s_inv-1,l_s_inv,r_s_inv,r_s_inv+l_s_inv-1,l_s,-1,0,0,0,0,'CSV_TYPE=ID11;CSV_INDEX='+str(CSV_loop)]
                SV_loop = SV_loop + 1
                SV_table.loc[SV_loop] = [SV_loop,ll_c,'Deletion',r_s_del2,r_s_del2+l_s_del2-1,l_s_del2,r_s_del2,r_s_del2+l_s_del2-1,l_s,-1,0,0,0,0,'CSV_TYPE=ID11;CSV_INDEX='+str(CSV_loop)]
                SV_loop = SV_loop + 1
                
                VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s_inv), 'rs' + str(VCF_loop), ''.join(original_string), ''.join(original_string_reverse), '.', 'PASS', 'SVTYPE=INV;CSV_TYPE=ID11;CSV_INDEX='+str(CSV_loop)]
                VCF_loop= VCF_loop + 1
                VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s_del1), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s_del1:r_s_del1+l_s_del1]) + tem_seq_post[r_s_del1+l_s_del1], tem_seq_post[r_s_del1+l_s_del1], '.', 'PASS', 'SVTYPE=DEL;CSV_TYPE=ID11;CSV_INDEX='+str(CSV_loop)]
                VCF_loop = VCF_loop+ 1 
                VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s_del2), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s_del2:r_s_del2+l_s_del2]) + tem_seq_post[r_s_del2+l_s_del2], tem_seq_post[r_s_del2+l_s_del2], '.', 'PASS', 'SVTYPE=DEL;CSV_TYPE=ID11;CSV_INDEX='+str(CSV_loop)]
                VCF_loop = VCF_loop+ 1 
                CSV_loop = CSV_loop + 1
                for ll_rever in range(r_s_inv,r_s_inv+l_s_inv):
                    tem_seq_post[ll_rever] = copy.deepcopy(original_string_reverse[int(ll_rever-r_s_inv)])
                tem_seq_post[r_s_del1:(r_s_del1+l_s_del1)] = '-'*l_s_del1
                tem_seq_post[r_s_del2:(r_s_del2+l_s_del2)] = '-'*l_s_del2
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post

#! ID12: Inversion with 5' Flanking Duplication (dupInv)

def ID12_dupinv_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_dupInv_csv,gaussian_ID12, times, ll_c, chr_id, len_seg_refine):
    print('ID12 (DupInv):'+str(num_dupInv_csv))
    for dupinv_num in range (0, num_dupInv_csv):
        # 采样 dupinversion 区域的起点和长度
        r_s = sample(unblock_region_sv, 1)[0]
        # sample the whole length for this csv
        l_s= np.random.choice(gaussian_ID12)
        # length of inversion is a subset
        l_s_inv = sample(range(int(l_s/2),int(l_s*11/12)),1)[0]
        #print('l_s, l_s_inv:'+str(l_s)+','+str(l_s_inv))
        # 计算我们重新采样的次数
        circular_count_dupInv = 0
        # 指示器：如果重新采样超过50次
        circular_count_dupInv_break = 0
        ### 需要重新采样的情况
        ### 情况 1：如果 Delinvsion 超过 real_con1 的末端
        ### 情况 2：如果 Delinvsion 的末端不在 unblock_region_sv 中
        ### 情况 3：先前的删除、插入（转位）和 Delinvsions 的重叠
        while (r_s+l_s_inv>len_seg_refine) or (r_s+l_s_inv-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s+l_s_inv-1)-unblock_region_sv.index(r_s)<l_s_inv-1):
            # 选择可能的删除点
            r_s = sample(unblock_region_sv, 1)[0]
            l_s = np.random.choice(gaussian_ID12)
            # length of inversion is a subset, a wider range
            l_s_inv = sample(range(int(l_s/2),int(l_s*11/12)),1)[0]
            # 计算重新采样的次数
            circular_count_dupInv = circular_count_dupInv + 1
            if circular_count_dupInv>times:
                circular_count_dupInv_break = 1
                print("Warning: ID12 (dupInv), sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_dupInv_break = 0
        if not circular_count_dupInv_break:
            circular_count_dupInv = 0
            
            unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s_inv+1)))
            
            # 5' Flanking Duplication
            l_s_dup = l_s - l_s_inv
            # start of the copied region
            remain_index2 = r_s
            # inserted pos is before r_s
            remain_index22 = r_s-1
            copied_string = tem_seq_post[remain_index2:remain_index2+l_s_dup]
            
            Ins_dic_sv[remain_index22] = ''.join(copied_string)
            
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Duplication',remain_index2,remain_index2+l_s_dup-1,l_s_dup,remain_index22,remain_index22, l_s,-1,0,0,0,0,'CSV_TYPE=ID12;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index22), 'rs' + str(VCF_loop), tem_seq_post[remain_index22], tem_seq_post[remain_index22]+''.join(copied_string), '.', 'PASS', 'SVTYPE=DUP;CSV_TYPE=ID12;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1
            
            # Inversion
            r_s_inv = r_s
            original_string = copy.deepcopy(''.join(tem_seq_post[r_s_inv:r_s_inv+l_s_inv]))
            original_string_reverse = DNA_complement(original_string[::-1])
            
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Inversion',r_s_inv,r_s_inv+l_s_inv-1,l_s_inv,r_s_inv,r_s_inv+l_s_inv-1,l_s,-1,0,0,0,0,'CSV_TYPE=ID12;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            
            VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s_inv), 'rs' + str(VCF_loop), ''.join(original_string), ''.join(original_string_reverse), '.', 'PASS', 'SVTYPE=INV;CSV_TYPE=ID12;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1

            CSV_loop = CSV_loop + 1
            for ll_rever in range(r_s_inv,r_s_inv+l_s_inv):
                tem_seq_post[ll_rever] = copy.deepcopy(original_string_reverse[int(ll_rever-r_s_inv)])
            
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post

#! ID13: Inversion with 3' Flanking Duplication (Invdup)

def ID13_invdup_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_Invdup_csv,gaussian_ID13, times, ll_c, chr_id, len_seg_refine):
    print('ID13 (InvDup):'+str(num_Invdup_csv))
    for invdup_num in range (0, num_Invdup_csv):
        # 采样 inversion+dup区域的起点和长度
        r_s = sample(unblock_region_sv, 1)[0]
        # sample the whole length for this csv
        l_s= np.random.choice(gaussian_ID13)
        # length of inversion is a subset
        l_s_inv = sample(range(int(l_s/2),int(l_s*11/12)),1)[0]
        #print('l_s, l_s_inv:'+str(l_s)+','+str(l_s_inv))
        # 计算我们重新采样的次数
        circular_count_dupInv = 0
        # 指示器：如果重新采样超过50次
        circular_count_dupInv_break = 0
        ### 需要重新采样的情况
        ### 情况 1：如果 Delinvsion 超过 real_con1 的末端
        ### 情况 2：如果 Delinvsion 的末端不在 unblock_region_sv 中
        ### 情况 3：先前的删除、插入（转位）和 Delinvsions 的重叠
        while (r_s+l_s_inv>len_seg_refine) or (r_s+l_s_inv-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s+l_s_inv-1)-unblock_region_sv.index(r_s)<l_s_inv-1):
            # 选择可能的删除点
            r_s = sample(unblock_region_sv, 1)[0]
            l_s = np.random.choice(gaussian_ID13)
            # length of inversion is a subset, a wider range
            l_s_inv = sample(range(int(l_s/2),int(l_s*11/12)),1)[0]
            # 计算重新采样的次数
            circular_count_dupInv = circular_count_dupInv + 1
            if circular_count_dupInv>times:
                circular_count_dupInv_break = 1
                print("Warning: ID13 (Invdup), sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_dupInv_break = 0
        if not circular_count_dupInv_break:
            circular_count_dupInv = 0
            
            unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s_inv+1)))
            
            # 5' Flanking Duplication
            l_s_dup = l_s - l_s_inv
            # start of the copied region
            #remain_index2-tem_copied_base+1
            remain_index2 =  r_s + l_s_inv - l_s_dup
            # inserted pos is before r_s
            remain_index22 = r_s + l_s_inv
            copied_string = tem_seq_post[remain_index2:remain_index2+l_s_dup]
            
            Ins_dic_sv[remain_index22] = ''.join(copied_string)
            
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Duplication',remain_index2,remain_index2+l_s_dup-1,l_s_dup,remain_index22,remain_index22, l_s,-1,0,0,0,0,'CSV_TYPE=ID13;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index22), 'rs' + str(VCF_loop), tem_seq_post[remain_index22], tem_seq_post[remain_index22]+''.join(copied_string), '.', 'PASS', 'SVTYPE=DUP;CSV_TYPE=ID13;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1
            
            # Inversion
            r_s_inv = r_s
            original_string = copy.deepcopy(''.join(tem_seq_post[r_s_inv:r_s_inv+l_s_inv]))
            original_string_reverse = DNA_complement(original_string[::-1])
            
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Inversion',r_s_inv,r_s_inv+l_s_inv-1,l_s_inv,r_s_inv,r_s_inv+l_s_inv-1,l_s,-1,0,0,0,0,'CSV_TYPE=ID13;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            
            VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s_inv), 'rs' + str(VCF_loop), ''.join(original_string), ''.join(original_string_reverse), '.', 'PASS', 'SVTYPE=INV;CSV_TYPE=ID13;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1

           
            for ll_rever in range(r_s_inv,r_s_inv+l_s_inv):
                tem_seq_post[ll_rever] = copy.deepcopy(original_string_reverse[int(ll_rever-r_s_inv)])
            
            CSV_loop = CSV_loop + 1
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post

#! ID14: Paired-duplication inversion (dupInvdup)

def ID14_dupinvdupprocess(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_dupinvdup_csv,gaussian_ID14, times, ll_c, chr_id, len_seg_refine):
    print('ID14 (dupInvdup):'+str(num_dupinvdup_csv))
    for dupinvdupnum in range (0, num_dupinvdup_csv):
        # 采样 inversion+dup区域的起点和长度
        r_s = sample(unblock_region_sv, 1)[0]
        # sample the whole length for this csv
        l_s= np.random.choice(gaussian_ID14)
        # length of inversion is a subset
        l_s_inv = sample(range(int(l_s/2),int(l_s*11/12)),1)[0]
        #print('l_s, l_s_inv:'+str(l_s)+','+str(l_s_inv))
        # 计算我们重新采样的次数
        circular_count_dupinvdup = 0
        # 指示器：如果重新采样超过50次
        circular_count_dupinvdup_break = 0
        ### 需要重新采样的情况
        ### 情况 1：如果 Delinvsion 超过 real_con1 的末端
        ### 情况 2：如果 Delinvsion 的末端不在 unblock_region_sv 中
        ### 情况 3：先前的删除、插入（转位）和 Delinvsions 的重叠
        while (r_s+l_s_inv>len_seg_refine) or (r_s+l_s_inv-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s+l_s_inv-1)-unblock_region_sv.index(r_s)<l_s_inv-1):
            # 选择可能的删除点
            r_s = sample(unblock_region_sv, 1)[0]
            l_s = np.random.choice(gaussian_ID14)
            # length of inversion is a subset, a wider range
            l_s_inv = sample(range(int(l_s/2),int(l_s*11/12)),1)[0]
            # 计算重新采样的次数
            circular_count_dupinvdup = circular_count_dupinvdup + 1
            if circular_count_dupinvdup>times:
                circular_count_dupinvdup_break = 1
                print("Warning: ID14 (dupInvdup), sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_dupinvdup_break = 0
        if not circular_count_dupinvdup_break:
            circular_count_dupinvdup = 0
            
            unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s_inv+1)))
            
            # 5' Flanking Duplication
            l_s_dup = l_s - l_s_inv
            #length of the 1st dup
            l_s_dup1 = sample(range(int(l_s_dup/4),int(l_s_dup*3/4)),1)[0]
            #length of the 2nd dup
            l_s_dup2 = l_s_dup - l_s_dup1
            # start of the copied region
            #remain_index2-tem_copied_base+1
            # start of the copied region
            #! 5' dup
            #start of the copied region
            remain_index2_dup1 = r_s
     
            remain_index22_dup1 = r_s - 1
            copied_string1 = tem_seq_post[remain_index2_dup1:remain_index2_dup1+l_s_dup1]
            Ins_dic_sv[remain_index22_dup1] = ''.join(copied_string1)
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Duplication',remain_index2_dup1,remain_index2_dup1+l_s_dup1-1,l_s_dup1,remain_index22_dup1,remain_index22_dup1, l_s,-1,0,0,0,0,'CSV_TYPE=ID14;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index22_dup1), 'rs' + str(VCF_loop), tem_seq_post[remain_index22_dup1], tem_seq_post[remain_index22_dup1]+''.join(copied_string1), '.', 'PASS', 'SVTYPE=DUP;CSV_TYPE=ID14;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1


            remain_index2_dup2 = r_s + l_s_inv - l_s_dup2 
            remain_index22_dup2 = r_s + l_s_inv
            
            
            copied_string2 = tem_seq_post[remain_index2_dup2:remain_index2_dup2+l_s_dup2]
            Ins_dic_sv[remain_index22_dup2] = ''.join(copied_string2)
            
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Duplication',remain_index2_dup2,remain_index2_dup2+l_s_dup2-1,l_s_dup2,remain_index22_dup2,remain_index22_dup2, l_s,-1,0,0,0,0,'CSV_TYPE=ID14;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index22_dup2), 'rs' + str(VCF_loop), tem_seq_post[remain_index22_dup2], tem_seq_post[remain_index22_dup2]+''.join(copied_string2), '.', 'PASS', 'SVTYPE=DUP;CSV_TYPE=ID14;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1
            
            # Inversion
            r_s_inv = r_s
            original_string = copy.deepcopy(''.join(tem_seq_post[r_s_inv:r_s_inv+l_s_inv]))
            original_string_reverse = DNA_complement(original_string[::-1])
            
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Inversion',r_s_inv,r_s_inv+l_s_inv-1,l_s_inv,r_s_inv,r_s_inv+l_s_inv-1,l_s,-1,0,0,0,0,'CSV_TYPE=ID14;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            
            VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s_inv), 'rs' + str(VCF_loop), ''.join(original_string), ''.join(original_string_reverse), '.', 'PASS', 'SVTYPE=INV;CSV_TYPE=ID14;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1

            for ll_rever in range(r_s_inv,r_s_inv+l_s_inv):
                tem_seq_post[ll_rever] = copy.deepcopy(original_string_reverse[int(ll_rever-r_s_inv)])
            
            CSV_loop = CSV_loop + 1
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post

#! ID15: Inversion with 5' Flanking Duplication and 3' Flanking Deletion (dupInvdel)


def ID15_dupinv_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_dupInvdel_csv,gaussian_ID15, times, ll_c, chr_id, len_seg_refine):
    print('ID15 (dupInvdel):'+str(num_dupInvdel_csv))
    for dupinv_num in range (0, num_dupInvdel_csv):
        # 采样 dupinversion 区域的起点和长度
        r_s = sample(unblock_region_sv, 1)[0]
        # sample the whole length for this csv
        l_s= np.random.choice(gaussian_ID15)
        # length of inversion is a subset
        l_s_invdel = sample(range(int(l_s/2),int(l_s*11/12)),1)[0]
        #print('l_s, l_s_invdel:'+str(l_s)+','+str(l_s_invdel))
        # 计算我们重新采样的次数
        circular_count_dupInvdel = 0
        # 指示器：如果重新采样超过50次
        circular_count_dupInvdel_break = 0
        ### 需要重新采样的情况
        ### 情况 1：如果 Delinvsion 超过 real_con1 的末端
        ### 情况 2：如果 Delinvsion 的末端不在 unblock_region_sv 中
        ### 情况 3：先前的删除、插入（转位）和 Delinvsions 的重叠
        while (r_s+l_s_invdel>len_seg_refine) or (r_s+l_s_invdel-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s+l_s_invdel-1)-unblock_region_sv.index(r_s)<l_s_invdel-1):
            # 选择可能的删除点
            r_s = sample(unblock_region_sv, 1)[0]
            l_s = np.random.choice(gaussian_ID15)
            # length of inversion is a subset, a wider range
            l_s_invdel = sample(range(int(l_s/2),int(l_s*11/12)),1)[0]
            # 计算重新采样的次数
            circular_count_dupInvdel = circular_count_dupInvdel + 1
            if circular_count_dupInvdel>times:
                circular_count_dupInvdel_break = 1
                print("Warning: ID15 (dupInvdel), sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_dupInvdel_break = 0
        if not circular_count_dupInvdel_break:
            circular_count_dupInvdel = 0
            
            unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s_invdel+1)))
            cut_point = sample(range(int(l_s_invdel/4),int(l_s_invdel*3/4)),1)[0]

            r_s_del=r_s+cut_point
            l_s_del=l_s_invdel-cut_point

            r_s_inv=r_s
            l_s_inv=cut_point

            # 5' Flanking Duplication
            l_s_dup = l_s - l_s_invdel
            # start of the copied region
            remain_index2 = r_s
            # inserted pos is before r_s
            remain_index22 = r_s-1
            copied_string = tem_seq_post[remain_index2:remain_index2+l_s_dup]
            
            Ins_dic_sv[remain_index22] = ''.join(copied_string)
            
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Duplication',remain_index2,remain_index2+l_s_dup-1,l_s_dup,remain_index22,remain_index22, l_s,-1,0,0,0,0,'CSV_TYPE=ID15;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index22), 'rs' + str(VCF_loop), tem_seq_post[remain_index22], tem_seq_post[remain_index22]+''.join(copied_string), '.', 'PASS', 'SVTYPE=DUP;CSV_TYPE=ID15;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1
            
            # Inversion
            r_s_inv = r_s
            original_string = copy.deepcopy(''.join(tem_seq_post[r_s_inv:r_s_inv+l_s_inv]))
            original_string_reverse = DNA_complement(original_string[::-1])
            
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Inversion',r_s_inv,r_s_inv+l_s_inv-1,l_s_inv,r_s_inv,r_s_inv+l_s_inv-1,l_s,-1,0,0,0,0,'CSV_TYPE=ID15;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            
            VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s_inv), 'rs' + str(VCF_loop), ''.join(original_string), ''.join(original_string_reverse), '.', 'PASS', 'SVTYPE=INV;CSV_TYPE=ID15;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1

            for ll_rever in range(r_s_inv,r_s_inv+l_s_inv):
                tem_seq_post[ll_rever] = copy.deepcopy(original_string_reverse[int(ll_rever-r_s_inv)])

            # Deletion
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Deletion',r_s_del,r_s_del+l_s_del-1,l_s_del,r_s_del,r_s_del+l_s_del-1,l_s,-1,0,0,0,0,'CSV_TYPE=ID15;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s_del), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s_del:r_s_del+l_s_del]) + tem_seq_post[r_s_del+l_s_del], tem_seq_post[r_s_del+l_s_del], '.', 'PASS', 'SVTYPE=DEL;CSV_TYPE=ID15;CSV_INDEX='+str(CSV_loop)]
            VCF_loop = VCF_loop+ 1 
            tem_seq_post[r_s_del:(r_s_del+l_s_del)] = '-'*l_s_del
            
            CSV_loop = CSV_loop + 1
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post

#! ID16: Inversion with 5' Flanking Deletion and 3' Flanking Duplication (delInvdup)

def ID16_dupinv_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_delInvdup_csv,gaussian_ID16, times, ll_c, chr_id, len_seg_refine):
    print('ID16 (DupInvdup):'+str(num_delInvdup_csv))
    for dupinv_num in range (0, num_delInvdup_csv):
        # 采样 dupinversion 区域的起点和长度
        r_s = sample(unblock_region_sv, 1)[0]
        # sample the whole length for this csv
        l_s= np.random.choice(gaussian_ID16)
        # length of inversion is a subset
        l_s_invdel = sample(range(int(l_s/2),int(l_s*11/12)),1)[0]
        #print('l_s, l_s_invdel:'+str(l_s)+','+str(l_s_invdel))
        # 计算我们重新采样的次数
        circular_count_delInvdup = 0
        # 指示器：如果重新采样超过50次
        circular_count_delInvdup_break = 0
        ### 需要重新采样的情况
        ### 情况 1：如果 Delinvsion 超过 real_con1 的末端
        ### 情况 2：如果 Delinvsion 的末端不在 unblock_region_sv 中
        ### 情况 3：先前的删除、插入（转位）和 Delinvsions 的重叠
        while (r_s+l_s_invdel>len_seg_refine) or (r_s+l_s_invdel-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s+l_s_invdel-1)-unblock_region_sv.index(r_s)<l_s_invdel-1):
            # 选择可能的删除点
            r_s = sample(unblock_region_sv, 1)[0]
            l_s = np.random.choice(gaussian_ID16)
            # length of inversion is a subset, a wider range
            l_s_invdel = sample(range(int(l_s/2),int(l_s*11/12)),1)[0]
            
            # 计算重新采样的次数
            circular_count_delInvdup = circular_count_delInvdup + 1
            if circular_count_delInvdup>times:
                circular_count_delInvdup_break = 1
                print("Warning: ID16 (DupInvdup), sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_delInvdup_break = 0
        if not circular_count_delInvdup_break:
            circular_count_delInvdup = 0
            
            unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s_invdel+1)))
            cut_point = sample(range(int(l_s_invdel/4),int(l_s_invdel*3/4)),1)[0]
            # print('ID16 r_s,l_s:'+str(r_s)+','+str(l_s))
            # print('ID16 l_s_invdel:'+str(l_s_invdel))
            # print('ID16 cut point:'+str(cut_point))
            r_s_del=r_s
            l_s_del=cut_point
            r_s_inv=r_s+cut_point
            l_s_inv=l_s_invdel-cut_point

            # 3' Flanking Duplication
            l_s_dup = l_s - l_s_invdel
            # start of the copied region
            #remain_index2-tem_copied_base+1
            remain_index2 =  r_s + l_s_invdel - l_s_dup
            # inserted pos is before r_s
            remain_index22 = r_s + l_s_invdel
            copied_string = tem_seq_post[remain_index2:remain_index2+l_s_dup]
            
            Ins_dic_sv[remain_index22] = ''.join(copied_string)
            
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Duplication',remain_index2,remain_index2+l_s_dup-1,l_s_dup,remain_index22,remain_index22, l_s,-1,0,0,0,0,'CSV_TYPE=ID16;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index22), 'rs' + str(VCF_loop), tem_seq_post[remain_index22], tem_seq_post[remain_index22]+''.join(copied_string), '.', 'PASS', 'SVTYPE=DUP;CSV_TYPE=ID16;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1
            
            # Inversion
            
            original_string = copy.deepcopy(''.join(tem_seq_post[r_s_inv:r_s_inv+l_s_inv]))
            original_string_reverse = DNA_complement(original_string[::-1])
            
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Inversion',r_s_inv,r_s_inv+l_s_inv-1,l_s_inv,r_s_inv,r_s_inv+l_s_inv-1,l_s,-1,0,0,0,0,'CSV_TYPE=ID16;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            
            VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s_inv), 'rs' + str(VCF_loop), ''.join(original_string), ''.join(original_string_reverse), '.', 'PASS', 'SVTYPE=INV;CSV_TYPE=ID16;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1

            for ll_rever in range(r_s_inv,r_s_inv+l_s_inv):
                tem_seq_post[ll_rever] = copy.deepcopy(original_string_reverse[int(ll_rever-r_s_inv)])
            
            # Deletion
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Deletion',r_s_del,r_s_del+l_s_del-1,l_s_del,r_s_del,r_s_del+l_s_del-1,l_s,-1,0,0,0,0,'CSV_TYPE=ID16;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s_del), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s_del:r_s_del+l_s_del]) + tem_seq_post[r_s_del+l_s_del], tem_seq_post[r_s_del+l_s_del], '.', 'PASS', 'SVTYPE=DEL;CSV_TYPE=ID16;CSV_INDEX='+str(CSV_loop)]
            VCF_loop = VCF_loop+ 1 
            tem_seq_post[r_s_del:(r_s_del+l_s_del)] = '-'*l_s_del

            CSV_loop = CSV_loop + 1
        else:
            break

    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post

#! ID17: Inverted Duplication with Flanking Triplication (dupTRIPdup-INV) 
def ID17_dupTRIPdup_INV_process(unblock_region_sv, True_ID17_number, times,  real_con1, SV_table, VCF_table, ll_c, chr_id, tem_seq_post, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, gaussian_ID17):
    print("ID17 (dupTRIPdup-INV):"+str(True_ID17_number))
    #inserted_site_collection = sample(unblock_region_sv,full_number_base)
    all_selected_ID17_SV= sample(unblock_region_sv,True_ID17_number)
    # 打印所有被选择的位点的信息
    other_sites = list(set(unblock_region_sv)-set(all_selected_ID17_SV))
    #InvDuplication: copy and paste
    for remain_index2 in all_selected_ID17_SV:##others### questions
        # sample the whole length for this csv
        l_s= np.random.choice(gaussian_ID17)
        circular_count_ID17 = 0
        circular_count_ID17_break = 0
        while len(set(range(remain_index2-l_s+1,remain_index2+1))&set(unblock_region_sv)) < l_s:
            # select the possibile copied length 
            #sites that do note have general insertions

            remain_index2 = sample(other_sites,1)[0]
            # sample the whole length for this csv
            l_s= np.random.choice(gaussian_ID17)
            
            circular_count_ID17 = circular_count_ID17 + 1
            if circular_count_ID17>times:
                circular_count_ID17_break = 1
                print("Warning: ID17 (dupTRIPdup-INV), sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_ID17_break = 0

        #ratio_re_ID17:neighbor InvDuplication (right after the copied area)
        if not circular_count_ID17_break:
            circular_count_ID17=0
            #!
            # unblock_region_sv = list(set(unblock_region_sv) -set(range(remain_index2-l_s+1,remain_index2)))
            original_string_dup = copy.deepcopy(real_con1[remain_index2-l_s+1:remain_index2+1])
            #! original_string_dup_reverse = original_string_dup[::-1]
            original_string_dup_reverse = DNA_complement(original_string_dup[::-1])
            
            # Tandem InvDuplication
            
            #later ins cannot be sampled from this point in other_sites
            if remain_index2 in other_sites:
                other_sites.remove(remain_index2)
            # ins and del cannot neighbor
            if remain_index2 in unblock_region_sv:
                unblock_region_sv.remove(remain_index2)
            if remain_index2+1 in unblock_region_sv:
                unblock_region_sv.remove(remain_index2+1)
            
            
            #! select part of the sequence
            # select part of the sequence
            start = sample(range(len(original_string_dup)), 1)[0]
            end = sample(range(start, len(original_string_dup)), 1)[0]
            add_string = original_string_dup[start:end+1]
        
            inserted_string = ''.join(copy.deepcopy(original_string_dup_reverse))+ ''.join(copy.deepcopy(add_string))
            Ins_dic_sv[remain_index2] = inserted_string
            
            ins_len_ID17 = len(Ins_dic_sv[remain_index2])
            # print("ID17 len:"+str(ins_len_ID17))
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'ID17',remain_index2-l_s+1,remain_index2,l_s,remain_index2,remain_index2,ins_len_ID17,-1,0,0,0,0,'CSV_TYPE=ID17;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1
            
            #ins
            VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index2), 'rs' + str(VCF_loop), tem_seq_post[remain_index2], tem_seq_post[remain_index2]+ inserted_string, '.', 'PASS', 'SVTYPE=TanInvDup;CSV_TYPE=ID17;CSV_INDEX='+str(CSV_loop)]
            VCF_loop= VCF_loop + 1
            
            CSV_loop = CSV_loop + 1
            
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post


#! ID18: Insertion with Deletion (INSdel)

#! unbalancedTranslocation + 5' Flanking deletion
def ID18_INSdel_process(SV_table, VCF_table, unblock_region_sv, SV_loop,CSV_loop,gaussian_ID18, VCF_loop, Ins_dic_sv, tem_seq_post, num_insdel_csv, times, ll_c, len_seg_refine,chr_id):
    
    print('ID18 (INSdel):'+str(num_insdel_csv))
    for s in range(0, num_insdel_csv):
        ### 循环每次转位
        ### 长转位的删除部分：切割 A
        r_s = sample(unblock_region_sv, 1)[0]
        # sample the whole length for this csv
        l_s= int(np.random.choice(gaussian_ID18)/2)
        # 如果长转位超过终端或先前的删除，并且此删除重叠，则重新采样起点和长度
        circular_count_trans=0
        # 指示器：如果重新采样超过50次
        circular_count_trans_break = 0
        while (r_s+l_s > len_seg_refine) or (r_s+l_s-1 not in unblock_region_sv) or ((unblock_region_sv.index(r_s+l_s-1)-unblock_region_sv.index(r_s)<l_s-1)):
            r_s = sample(unblock_region_sv, 1)[0]
            l_s= int(np.random.choice(gaussian_ID18)/2)
            # 计算重新采样的次数
            circular_count_trans = circular_count_trans + 1
            if circular_count_trans>times:
                circular_count_trans_break = 1
                print("Warning: ID18 (INSdel), sampling times exceeds "+str(times)+' times.')
                break
            else:
                circular_count_trans_break = 0
        if not circular_count_trans_break:
            circular_count_trans=0
            ## 选择 A (切割)后更新索引集
            # 对于下一个删除，不能是右或左邻居
            unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s+1)))
            
            ### deletion region
            ### 选择 B 的删除起点
            r_s2 = sample(unblock_region_sv,1)[0]
            # 采样del的长度
            l_s2 = int(np.random.choice(gaussian_ID18)/2)
            # 如果长转位超过终端或先前的删除，并且 B 重叠，则重新采样起点和长度
            circular_count_trans=0
            while (r_s2+l_s2>len_seg_refine) or (r_s2+l_s2-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s2+l_s2-1)-unblock_region_sv.index(r_s2)<l_s2-1):
                r_s2 = sample(unblock_region_sv,1)[0]
                # 采样del的长度
                l_s2 = int(np.random.choice(gaussian_ID18)/2)
                # 计算重新采样的次数
                circular_count_trans = circular_count_trans + 1
                if circular_count_trans>times:
                    circular_count_trans_break = 1
                    break
                else:
                    circular_count_trans_break = 0

           
            # 不平衡转位：切割 A 并粘贴到 B
            ## 原始位置也被删除（切割 A）
            ### 插入位点：紧邻删除的最后一个Base
            ins_trans_loc = r_s2+l_s2-1
           
            # 插入的片段是之前采样的切割 A 区域
            # 插入的位置是deletion的最后区域
            # 将区域 A 粘贴到位置 B (ins_trans_loc)
            Ins_dic_sv[ins_trans_loc] = copy.deepcopy(''.join(tem_seq_post[r_s:(r_s+l_s)]))
            # 删除和插入的片段不能是左，右邻居和重叠。
            # 删除区域 B 周围的位置
            unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s2-1,r_s2+l_s2+1)))
            
            # 写 SV 表
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Translocation',r_s,r_s+l_s-1,l_s,ins_trans_loc,ins_trans_loc,l_s,0,0,0,0,0,'CSV_TYPE=ID18;CSV_INDEX='+str(CSV_loop)]### 切割和粘贴
            SV_loop = SV_loop + 1
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Deletion',r_s2,r_s2+l_s2-1,l_s2,r_s2,r_s2+l_s2-1,l_s,-1,0,0,0,0,'CSV_TYPE=ID18;CSV_INDEX='+str(CSV_loop)]
            SV_loop = SV_loop + 1

            VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s:r_s+l_s]) + tem_seq_post[r_s+l_s], tem_seq_post[r_s+l_s], '.', 'PASS', 'SVTYPE=DEL;CSV_TYPE=ID18;CSV_INDEX='+str(CSV_loop)]
            VCF_loop = VCF_loop+ 1
            VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s2), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s2:r_s2+l_s2]) + tem_seq_post[r_s2+l_s2], tem_seq_post[r_s2+l_s2], '.', 'PASS', 'SVTYPE=DEL;CSV_TYPE=ID18;CSV_INDEX='+str(CSV_loop)]
            VCF_loop = VCF_loop+ 1
            VCF_table.loc[VCF_loop] = [str(chr_id), str(ins_trans_loc), 'rs' + str(VCF_loop), tem_seq_post[ins_trans_loc], tem_seq_post[ins_trans_loc] + ''.join(tem_seq_post[r_s:(r_s+l_s)]), '.', 'PASS', 'SVTYPE=INS;CSV_TYPE=ID18;CSV_INDEX='+str(CSV_loop)]
            VCF_loop = VCF_loop + 1

            CSV_loop = CSV_loop + 1

            #删除被复制的片段
            tem_seq_post[r_s:(r_s+l_s)] = '-'*l_s
            #删除deletion的片段
            tem_seq_post[r_s2:(r_s2+l_s2)] = '-'*l_s2
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post


#! Long DEL

def long_del_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop,tem_seq_post, sv_del, condition_dist_sv_del, len_SV_del, len_seg_refine, times, ll_c, chr_id):
    #print('Deletion:'+str(sv_del))
    for del_num in range(0, sv_del):
        ### Sample the first long deletion
        # sample a start point in the unblock region for deletion
        r_s = sample(unblock_region_sv, 1)[0]
        # sample the index of the first deletion's length
        l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_del)
        # .index(a): find the index of the first occurrence of a given element, i.e. a, in a list.
        # find the index of the event==1
        l_s_index = list(l_s_vec).index(1)
        # length of the deletion, chosen by sampled index on the len_SV_del=[50, 100, 500]
        l_s = len_SV_del[int(l_s_index)]
        # l_s = np.random.poisson(lamnuma_del, 1)[0]

        ### situation when resampling is needed
        ### situation 1: if deletion exceeds the end of real_con1 
        ### situation 2: if end of deletion is not in unblock_region_sv
        ### situation 3: overlap of deletion and translocation
        # count the number of times that we resample
        circular_count_del = 0
        # indicator: if resampling exceeds 50 times
        circular_count_del_break = 0
        while (r_s+l_s>len_seg_refine) or (r_s+l_s-1 not in unblock_region_sv) or (unblock_region_sv.index(r_s+l_s-1)-unblock_region_sv.index(r_s)<l_s-1):
            # select the possibile deleltion point
            r_s = sample(unblock_region_sv, 1)[0]
            l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_del)
            l_s_index = list(l_s_vec).index(1)
            l_s = len_SV_del[int(l_s_index)]
            # count the times of resampling
            circular_count_del = circular_count_del + 1
            if circular_count_del>times:
                circular_count_del_break = 1
                break
            else:
                circular_count_del_break = 0
        # if sample a del start and length that is not overlapped with translocation(s)
        # update needed pos collections and prepare for the second sampling of long deletions
        if not circular_count_del_break:
            print('SV DEL:' + str(sv_del))
            # initialize the counting of resampling times
            circular_count_del = 0
            # block the part deleted in the first long del
            unblock_region_sv = list(set(unblock_region_sv) -set(range(r_s-1,r_s+l_s+1)))

            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Deletion',r_s,r_s+l_s-1,l_s,-1,-1,-1,-1,0,0,0,0,'.']
            SV_loop = SV_loop + 1      
            
            #The ‘REF’ column is set to the original segment plus the base that is left after the deletion. 
            #The ‘ALT’ column is set to the base that is left after the deletion.
            # Add a row to the VCF_table for the micro deletion
            #VCF_table[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), tem_seq_post[r_s:r_s+l_s] + tem_seq_post[r_s+l_s], tem_seq_post[r_s+l_s], '.', 'PASS', 'SVTYPE=DEL']
            VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s:r_s+l_s]) + tem_seq_post[r_s+l_s], tem_seq_post[r_s+l_s], '.', 'PASS', 'SVTYPE=DEL;SVLEN='+str(l_s)]
            VCF_loop = VCF_loop + 1
            
            #! update the sequence after record the deleted part
            #replace deleted bases with -
            tem_seq_post[r_s:(r_s+l_s)] = '-'*l_s
            
        else:
            break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, tem_seq_post
                                 
#! Long INS
def long_ins_process(unblock_region_sv,sv_ins, condition_dist_sv_ins, len_SV_ins, ins_selection, SV_table, VCF_table, ll_c, chr_id, tem_seq_post, SV_loop, VCF_loop, Ins_dic_sv):
    #print('Insertion:'+str(sv_ins))
    # long insertion part
    # indicator: if resampling exceeds 50 times
    circular_count_ins_break = 0
    ins_pos_collection = []

    # 对每一段进行处理
    if not unblock_region_sv:  # 如果segment是空集，输出警告
        print("Warning: empty unblock_region_sv")
        circular_count_ins_break = 1
    elif sv_ins == 0:
        circular_count_ins_break = 1
    else:
        circular_count_ins_break = 0

    # 计算需要选择的位点的数量
    if not circular_count_ins_break:
        print('SV INS:' + str(sv_ins))
        if not unblock_region_sv:
            print("warning: empty unblock_region_sv")
        else:
            len_seg_refine = max(unblock_region_sv)
            if len_seg_refine - 1 in unblock_region_sv:
                unblock_region_sv.remove(len_seg_refine - 1)
        # 如果当前段的长度小于需要选择的位点的数量，就选择所有的位点
        if len(unblock_region_sv) < sv_ins:
            ins_SV_pos = unblock_region_sv
        else:
            # 否则，随机选择位点
            ins_SV_pos = random.sample(unblock_region_sv, sv_ins)
        # 将被选择的位点添加到列表中
        ins_pos_collection.extend(ins_SV_pos)

        for remain_index2 in ins_pos_collection:
            # sample length of true insertion
            l_i_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_ins)
            l_i_ins = int(list(l_i_vec).index(1))
            l_i = len_SV_ins[l_i_ins]
            # initialize: inserted terms
            tem_ins = ''
            # a loop to choose inserted base
            for j in range(l_i):
                bexixuan_ = choice(ins_selection)
                tem_ins = tem_ins + bexixuan_

            Ins_dic_sv[remain_index2] = tem_ins

            if remain_index2 in unblock_region_sv:
                unblock_region_sv.remove(remain_index2)
            if remain_index2 + 1 in unblock_region_sv:
                unblock_region_sv.remove(remain_index2 + 1)
            

            # only one record in the table for each remain_index2
            SV_table.loc[SV_loop] = [SV_loop, ll_c, 'Insertion', remain_index2, remain_index2, l_i, -1, -1, -1, -1, 0, 0, 0, 0,'.']
            SV_loop = SV_loop + 1

            # Add a row to the VCF_table for the micro insertion
            VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index2), 'rs' + str(VCF_loop), tem_seq_post[remain_index2], tem_seq_post[remain_index2] + tem_ins.upper(), '.', 'PASS', 'SVTYPE=INS;SVLEN='+str(l_i)]
            VCF_loop = VCF_loop + 1

    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post

#! micro del

pai_pro_tem = [diff_ins_prob_del_real, diff_ins_prob_mis_real, diff_ins_prob_correct_real]
pai_pro_tem_ = list(np.array(pai_pro_tem) / sum(pai_pro_tem))

def micro_del_process(snv_del, unblock_region_sv, pai_pro_tem_, times, SV_table, VCF_table, ll_c, chr_id, tem_seq_post, SV_loop, VCF_loop, Ins_dic_sv):
    print('Micro del:'+str(snv_del))
    # micro_del#对每一段进行处理
    if not unblock_region_sv:  # 如果segment是空集，输出警告
        print("Warning: empty unblock_region_sv")
        circular_count_micro_del_break = 1
    else:
        circular_count_micro_del_break = 0

    # 计算需要选择的位点的数量
    if not circular_count_micro_del_break:
        ## micro deletion part
        len_unblock_region = len(unblock_region_sv)
        l_s_vec_ini = np.random.multinomial(n=len_unblock_region, pvals=pai_pro_tem_)

        max_del = len_unblock_region // 2

        if snv_del is None:
            del_snv_number = int(l_s_vec_ini[0])
        else:
            if snv_del < 1:
                input_del = int(len_unblock_region * snv_del)
            else:
                input_del = int(snv_del)

            del_snv_number = min(input_del, max_del)
            if input_del > max_del:
                print("Warning: The input for -snv_del is too large and has been automatically reduced. \
                    This is because each micro deletion event requires at least one position, and there must be at least one position between two micro deletion events. \
                        Therefore, the maximum number of micro deletions that can occur is half of the total number of positions available for events.")

        len_seg_refine = max(unblock_region_sv)
        # for each deletion
        for s in range(0, del_snv_number):
            r_s = sample(unblock_region_sv, 1)[0]
            # del_sites_indexes.append(r_s)
            # l_s = choice(del_length_list)
            # sample the length of deletion by a multinomial distribution
            l_s_vec = np.random.multinomial(n=1, pvals=[6 / 8, 1 / 8, 1 / 16, 1 / 32, 1 / 32])
            # the length of this deletion
            l_s = list(l_s_vec).index(1) + 1
            # resample if the end of this deletion is out of range

            # count the number of times that we resample
            circular_count_del = 0
            # indicator: if resampling exceeds 50 times
            circular_count_del_break = 0
            while (r_s + l_s > len_seg_refine) or (r_s + l_s - 1 not in unblock_region_sv) or (
                    unblock_region_sv.index(r_s + l_s - 1) - unblock_region_sv.index(r_s) < l_s - 1):
                # select the possibile deleltion point
                r_s = sample(unblock_region_sv, 1)[0]
                l_s_vec = np.random.multinomial(n=1, pvals=[6 / 8, 1 / 8, 1 / 16, 1 / 32, 1 / 32])
                # the length of this deletion
                l_s = list(l_s_vec).index(1) + 1
                # l_s = choice(del_length_list)
                # count the times of resampling
                circular_count_del = circular_count_del + 1
                if circular_count_del > times:
                    circular_count_del_break = 1
                    print('break')
                    break
                else:
                    circular_count_del_break = 0
            # if sample a del start and length that is not overlapped with translocation(s)
            # update needed pos collections and prepare for the second sampling of long deletions
            if not circular_count_del_break:
                # initialize the counting of resampling times
                circular_count_del = 0

                # block the part deleted in the first long del
                unblock_region_sv = list(set(unblock_region_sv) - set(range(r_s - 1, r_s + l_s + 1)))

            # for ll in del_sites_indexes:
            #     #tem_seq_post[int(ll)] = '-'
            #     SV_table[SV_loop] = [SV_loop,ll_c,'Micro_Del',ll,ll,1,-1,-1,-1,-1,0,0,0,0]
            #     SV_loop = SV_loop + 1
                SV_table.loc[SV_loop] = [SV_loop, ll_c, 'Micro_Del', r_s, r_s + l_s - 1, l_s, -1, -1, -1, -1, 0, 0, 0, 0, '.']
                SV_loop = SV_loop + 1

                # The ‘REF’ column is set to the original segment plus the base that is left after the deletion.
                # The ‘ALT’ column is set to the base that is left after the deletion.
                # Add a row to the VCF_table for the micro deletion
                # VCF_table[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), tem_seq_post[r_s:r_s+l_s] + tem_seq_post[r_s+l_s], tem_seq_post[r_s+l_s], '.', 'PASS', 'SVTYPE=microDEL']
                VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop),
                                           ''.join(tem_seq_post[r_s:r_s + l_s]) + tem_seq_post[r_s + l_s],
                                           tem_seq_post[r_s + l_s], '.', 'PASS', 'microDEL;LEN='+str(l_s)]

                VCF_loop = VCF_loop + 1

                # replace deleted bases with -
                tem_seq_post[r_s:(r_s + l_s)] = '-' * l_s
            else:
                break
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post

#! micro ins
### micro insertions
# 对每一段进行处理

def micro_ins_process(snv_ins, unblock_region_sv, diff_ins_prob_ins_real, ins_selection, SV_table, VCF_table, ll_c, chr_id, tem_seq_post, SV_loop, VCF_loop, Ins_dic_sv):
    print('Micro ins:'+str(snv_ins))
    # 如果segment是空集，输出警告
    if not unblock_region_sv:
        print("Warning: empty unblock_region_sv")
        circular_count_micro_ins_break = 1
    else:
        circular_count_micro_ins_break = 0

    # 计算需要选择的位点的数量
    if not circular_count_micro_ins_break:
        len_unblock_region = len(unblock_region_sv)
        if snv_ins is None:
            ins_snv_number = (np.random.binomial(len_unblock_region,diff_ins_prob_ins_real,1))[0]
        elif 0 <= snv_ins < 1:
            ins_snv_number = int(len_unblock_region * snv_ins)
        else:
            ins_snv_number = min(int(snv_ins), len_unblock_region)

        #all possitions for ins, a loop
        for number_ins in range(0,ins_snv_number):
            #the positions of insertions
            remain_index2 = sample(unblock_region_sv,1)[0]
            #length of insertion
            l_i_vec = np.random.multinomial(n=1, pvals=[6/8,1/8,1/16,1/32,1/32])
            l_i = list(l_i_vec).index(1)+1
            tem_ins=''
            #a loop to choose inserted base
            for j in range(l_i):
                bexixuan_ = choice(ins_selection)
                tem_ins = tem_ins + bexixuan_

            #record ins position and ins segment
            Ins_dic_sv[remain_index2] = tem_ins

            #update index set
            if remain_index2 in unblock_region_sv:
                unblock_region_sv.remove(remain_index2)
            if remain_index2+1 in unblock_region_sv:
                unblock_region_sv.remove(remain_index2+1)

            #only one record in the table for each remain_index2
            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Micro_Ins',remain_index2,remain_index2,l_i,-1,-1,-1,-1,0,0,0,0,'.']
            SV_loop = SV_loop + 1

            # Add a row to the VCF_table for the micro insertion
            VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index2), 'rs' + str(VCF_loop), tem_seq_post[remain_index2], tem_seq_post[remain_index2] + tem_ins.upper(), '.', 'PASS', 'microINS;LEN='+str(l_i)]
            VCF_loop = VCF_loop + 1
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post

#! start SNP

def snp_process(mis_snv_number, unblock_region_sv, base_list, substitution_matrix, mis_selection, SV_table, VCF_table, ll_c, chr_id, tem_seq_post, SV_loop, VCF_loop, Ins_dic_sv):
    print('SNP:'+str(mis_snv_number))
    # 如果segment是空集，输出警告
    if not unblock_region_sv:
        print("Warning: empty unblock_region_sv.")
        circular_count_snv_break = 1
    else:
        circular_count_snv_break = 0

    # 计算需要选择的位点的数量
    if not circular_count_snv_break:
        max_snp = len(unblock_region_sv)
        mis_snv_number = min(max_snp,mis_snv_number)
        print('maximum snp:'+str(max_snp))
        print('snp number:'+str(mis_snv_number))
        for number_mis in range(0, mis_snv_number):
            # pos of substitution
            ll = sample(unblock_region_sv, 1)[0]
            unblock_region_sv.remove(ll)
            # the base that is replaced
            if tem_seq_post[ll] == 'N':
                print('Error: SNP in Gap region')
            elif tem_seq_post[ll].upper() in base_list:
                ref_id = base_list.index(tem_seq_post[ll].upper())
                # selection the ref_id's probability distribution
                prob_dist = substitution_matrix[ref_id]
                # sample a column_index 
                column_index = np.random.choice(4, p=prob_dist)
                # choose a position for mismatch
                bexixuan_ = mis_selection[column_index]

                SV_table.loc[SV_loop] = [SV_loop, ll_c, 'Substitution', ll, ll, 1, -1, -1, -1, -1, 0, 0, 0, 0, '.']
                SV_loop = SV_loop + 1
                # Add a row to the VCF_table for the substitution
                VCF_table.loc[VCF_loop] = [str(chr_id), str(ll), 'rs' + str(VCF_loop), tem_seq_post[ll], bexixuan_.upper(), '.', 'PASS', 'SUB']
                VCF_loop = VCF_loop + 1

                tem_seq_post[int(ll)] = copy.deepcopy(bexixuan_)
            else:
                print("Error: Invalid base")
    return SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post

def CSV_finalize_table(SV_table_merged,ll_c, tem_ins_dic):
    tem_SV_table_merged = SV_table_merged[SV_table_merged.iloc[:,1]==ll_c]
    #tem_SV_table_merged = SV_table_merged
    #original start: start of A
    list_start1 = list(tem_SV_table_merged.iloc[:,3])
    #start of B (dulplication or balanced trans)
    list_start2 = list(tem_SV_table_merged.iloc[:,6])
    whole_start_abs = list_start1+list_start2
    #order SV from left
    #set:merge repeated sites (e.g. 5 mismatch 5.5 ins)
    whole_start_abs_set = sorted(list(set(whole_start_abs)))
    present_len = 0
    last_bone = 0
    #inserted term
    #tem_ins_dic = Whole_INS_con[ll_c]
    #as there is only one consensus for each cut
        
    for ll_var_index in range(len(whole_start_abs_set)):
        
        #! time 找到对应的行
        tem_SV_table_merged2 = tem_SV_table_merged[(tem_SV_table_merged['Original_start']==whole_start_abs_set[ll_var_index]) |\
                                        (tem_SV_table_merged['New_start']==whole_start_abs_set[ll_var_index])]

        for xun_nei_row in range(len(tem_SV_table_merged2)):
            tem_row = tem_SV_table_merged2.iloc[xun_nei_row,:]
            stand_line = int(tem_row[0])
            #A
            bone1s = tem_row[3]
            bone1e = tem_row[4]
            #B
            bone2s = tem_row[6]
            bone2e = tem_row[7]
            if whole_start_abs_set[ll_var_index] in list_start1:
            #ls_satrt1_index_df = int(list_start1.index(whole_start_abs_set[ll_var_index]))
            #stand_line = int(SV_table_merged.iloc[ls_satrt1_index_df,0])
            #tem_row = SV_table_merged.iloc[ls_satrt1_index_df,:]
                #class of SV
                if tem_row[2] in ['Substitution','Micro_Ins','Micro_Del','Deletion','Insertion','Inversion']:
                    if tem_row[2] in ['Deletion','Micro_Del']:
                        inster_number_bone = bone1s-last_bone-1
                        #index for consensus before start of current variation
                        present_len = present_len + inster_number_bone
                        #update last_bone as end of current variation
                        last_bone = bone1e
                        #deleted base has no new axis on consensus
                        SV_table_merged.iloc[stand_line,10] = -1
                        SV_table_merged.iloc[stand_line,11] = -1
                        SV_table_merged.iloc[stand_line,12] = -1
                        SV_table_merged.iloc[stand_line,13] = -1
                    elif tem_row[2] in ['Substitution']:
                        inster_number_bone = bone1s-last_bone
                        #one to one map
                        present_len = present_len + inster_number_bone
                        #bone1s=bone1e=5
                        last_bone = bone1e
                        SV_table_merged.iloc[stand_line,10] = present_len
                        SV_table_merged.iloc[stand_line,11] = present_len
                        SV_table_merged.iloc[stand_line,12] = -1
                        SV_table_merged.iloc[stand_line,13] = -1
                    elif tem_row[2] in ['Micro_Ins','Insertion']:
                        inster_number_bone = bone1s-last_bone
                        Ins_len_present = len(tem_ins_dic[bone1s])
                        #inserted position on consensus: one pos:+1, inserted after current base
                        SV_table_merged.iloc[stand_line,10] = present_len + inster_number_bone+1
                        SV_table_merged.iloc[stand_line,12] = -1
                        #on consensus: end of previous SV+ number of normal base+ inserted length
                        present_len = present_len + inster_number_bone+Ins_len_present
                        #end of current SV
                        last_bone = bone1e
                        SV_table_merged.iloc[stand_line,11] = present_len
                        SV_table_merged.iloc[stand_line,13] = -1
                    else:## this is the inversion
                        inster_number_bone = bone1s-last_bone
                        SV_table_merged.iloc[stand_line,10] = present_len + inster_number_bone
                        SV_table_merged.iloc[stand_line,12] = -1
                        #no loss from last_bone to bone1e
                        #????
                        #present_len = present_len + bone1e - last_bone
                        present_len = present_len + bone1e - last_bone
                        SV_table_merged.iloc[stand_line,11] = present_len
                        SV_table_merged.iloc[stand_line,13] = -1 
                        last_bone = bone1e
                        
                elif tem_row[2] in ['Duplication','TanInvDup','DisInvDup','DisDup','TanDup']:
                        #copy A to B (A no change)
                        #5-0=5
                        inster_number_bone = bone1s-last_bone
                        #Ins_len_present = len(tem_ins_dic[bone2s])
                        #length of the copied: A
                        #=6
                        tem_plate_len = SV_table_merged.iloc[stand_line,5]
                        #0+5=5
                        SV_table_merged.iloc[stand_line,10] = present_len + inster_number_bone
                        #SV_table_merged.iloc[stand_line,12] = present_len + inster_number_bone+1
                        present_len = present_len + inster_number_bone + tem_plate_len-1
                        #0+5+6-1=10
                        SV_table_merged.iloc[stand_line,11] = present_len 
                        #SV_table_merged.iloc[stand_line,13] = present_len
                        last_bone = bone1e
                elif tem_row[2] in ['ID17']:
                        #copy A to B (A no change)
                        #5-0=5
                        inster_number_bone = bone1s-last_bone
                        #Ins_len_present = len(tem_ins_dic[bone2s])
                        #length of the copied: A
                        #=6
                        tem_plate_len = SV_table_merged.iloc[stand_line,5]
                        #0+5=5
                        SV_table_merged.iloc[stand_line,10] = present_len + inster_number_bone
                        #SV_table_merged.iloc[stand_line,12] = present_len + inster_number_bone+1
                        present_len = present_len + inster_number_bone + tem_plate_len-1
                        #0+5+6-1=10
                        SV_table_merged.iloc[stand_line,11] = present_len 
                        #SV_table_merged.iloc[stand_line,13] = present_len
                        #!
                        last_bone = bone1e
                        
                elif tem_row[2] in ['Translocation']:
                    #balanced translocation
                    #A:5-10, B:12-18
                    if tem_row[9] == 1:
                        #ins B to A's pos:5-0-1=4
                        inster_number_bone = bone1s-last_bone-1
                        #0+4+1=5,the start of copied base is 5
                        SV_table_merged.iloc[stand_line,10] = present_len + inster_number_bone+1
                        #length of B: 18-12+1=7
                        Ins_len_present = len(tem_ins_dic[bone1s-1])
                        #0+4+7=11
                        #end of A:current SV end=11
                        present_len = present_len + inster_number_bone + Ins_len_present
                        SV_table_merged.iloc[stand_line,11] = present_len
                        last_bone = bone1e
                    #!unbalanced trans:
                    else:
                        inster_number_bone = bone1s-last_bone-1
                        #index for consensus before start of current variation
                        present_len = present_len + inster_number_bone
                        
                        #deleted base has no new axis on consensus
                        SV_table_merged.iloc[stand_line,10] = present_len+1
                        SV_table_merged.iloc[stand_line,11] = present_len+1
                        
                        #update last_bone as end of current variation
                        last_bone = bone1e
        
        
            else:### in the list2: pos of B (only duplication and trans)
                
                if tem_row[2] in ['Duplication','TanInvDup','DisInvDup','DisDup','TanDup']:
                    #if SV_table_merged.iloc[stand_line,10]==0:
                        #bone2s:B_start
                        #same as ins
                        inster_number_bone = bone2s-last_bone
                        #SV_table_merged.iloc[stand_line,10] = present_len + inster_number_bone+1
                        #SV_table_merged.iloc[stand_line,12] = present_len + inster_number_bone+1
                        SV_table_merged.iloc[stand_line,12] = present_len+inster_number_bone+1
                        Ins_len_present = len(tem_ins_dic[bone2s])
                        present_len = present_len + inster_number_bone+Ins_len_present
                        
                        SV_table_merged.iloc[stand_line,13] = present_len
                        last_bone = bone2e
                elif tem_row[2] in ['ID17']:
                    #if SV_table_merged.iloc[stand_line,10]==0:
                        #bone2s:B_start
                        #same as ins
                        inster_number_bone = bone2s-last_bone
                        #SV_table_merged.iloc[stand_line,10] = present_len + inster_number_bone+1
                        #SV_table_merged.iloc[stand_line,12] = present_len + inster_number_bone+1
                        SV_table_merged.iloc[stand_line,12] = present_len+inster_number_bone+1
                        Ins_len_present = len(tem_ins_dic[bone2s])
                        present_len = present_len + inster_number_bone+Ins_len_present
                        
                        SV_table_merged.iloc[stand_line,13] = present_len
                        last_bone = bone2e
                elif tem_row[2] in ['Translocation']:
                    #balanced: similar to A
                    if  tem_row[9] == 1:
                        inster_number_bone = bone2s-last_bone-1
                        SV_table_merged.iloc[stand_line,12] = present_len + inster_number_bone+1
                        #inserted A's length
                        Ins_len_present = len(tem_ins_dic[bone2s-1])
                        present_len = present_len + inster_number_bone + Ins_len_present
                        SV_table_merged.iloc[stand_line,13] = present_len
                        last_bone = bone2e
                    #unbalanced
                    else:
                        inster_number_bone = bone2s-last_bone-1
                        inster_number_bone = bone2s-last_bone
                        #A is a del
                        SV_table_merged.iloc[stand_line,10] = -1
                        SV_table_merged.iloc[stand_line,12] = present_len + inster_number_bone+1
                        #length of A
                        #Ins_len_present = len(tem_ins_dic[bone2s-1])
                        #Ins_dic_sv_seg[ins_trans_loc] = copy.deepcopy(''.join(tem_seq_post[r_s:(r_s+l_s)]))
                        #similar to insertion
                        Ins_len_present = len(tem_ins_dic[bone2s])
                        present_len = present_len + inster_number_bone + Ins_len_present
                        #A is a del
                        SV_table_merged.iloc[stand_line,11] = -1
                        SV_table_merged.iloc[stand_line,13] = present_len
                        last_bone = bone2e
    return SV_table_merged

def parse_args():
    parser = argparse.ArgumentParser(description='GenoWave')
    parser.add_argument('-ref', type=str, help='Input reference local path', default='default_ref')
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
    parser.add_argument('-block_region_bed_url',type=str, help='local path of the block region BED file', default=None)
    #CSV
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
    return parser.parse_args()


   
def main():
    start_time1 = time.time()
    args = parse_args()
    times = args.times
    
    # 设置main_url_empirical为固定的路径
    main_url_empirical = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'empirical') + '/'
    # save = args.save
    # rep = args.rep
    ref = args.ref
    ll_c = 0
    CSV_loop = 0
    SV_loop = 0
    VCF_loop = 0

    updated_con = []
    #record the position and content of insertions, to update the final consensus
    Ins_dic_sv = {}

    random.seed(args.seed)
    np.random.seed(args.seed) 

    # fasta文件路径
    fasta_file_path = args.ref
    fasta_file = pysam.FastaFile(fasta_file_path)

    # 获取最长的序列
    seqname = fasta_file.references[0]
    BestRefSeq = fasta_file.fetch(seqname)

    #! whole chr
    chr_id=seqname
    start_base=0
    end_base=len(BestRefSeq)
    real_con1 = copy.deepcopy(BestRefSeq[start_base:end_base+1]).upper()
    chr_length = len(real_con1)
    del BestRefSeq
    print('Length of ref:'+str(chr_length))

    #!Block N base positions

    # #collection of N positions
    # n_positions = [i for i, char in enumerate(real_con1) if char == 'N']
    # n_positions_set = set(n_positions)
    all_positions = set(range(len(real_con1)))
    # print('N pos:'+str(len(n_positions)))
   
    # all_region_sv=list(all_positions - n_positions_set)
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
        #block_region_bed = '~/data/test_data/TGS/hg002/chr21_block_unique.bed'
        #block_region_bed = args.block_region_bed_url
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
        
    #! initialization
    unblock_region_sv = copy.deepcopy(all_region_sv)
    tem_seq_post = copy.deepcopy(list(real_con1))
    
    del all_region_sv

    # 从.npy文件中读取数据
    len_SV_del = np.load(main_url_empirical+'len_SV_del.npy').tolist()
    condition_dist_sv_del = np.load(main_url_empirical+'condition_dist_sv_del.npy').tolist()

    #modifiable deletion length
    delmin = args.delmin
    delmax = args.delmax

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
    len_SV_ins = np.load(main_url_empirical+'len_SV_ins.npy').tolist()
    condition_dist_sv_ins = np.load(main_url_empirical+'condition_dist_sv_ins.npy').tolist()

    insmin = args.insmin
    insmax = args.insmax

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
    
    dupmin = args.dupmin
    dupmax = args.dupmax
    # 计算整数范围
    dup_range = np.arange(dupmin, dupmax + 1)
    # 计算每个整数的概率（均匀分布）
    # uniform_prob = 1 / len(dup_range)
    # # 创建新的概率向量
    # dup_sv_prob = [uniform_prob] * len(dup_range)
    
    # copied_base_sv_base = dup_range.tolist()
    # copied_base_sv_prob = dup_sv_prob
    
    # 设置参数
    # invmin = 84
    # invmax = 207683
    invmin = args.invmin
    invmax = args.invmax

    # 计算整数范围
    inv_range = np.arange(invmin, invmax + 1)
    # 计算每个整数的概率（均匀分布）
    uniform_prob = 1 / len(inv_range)
    # 创建新的概率向量
    inv_sv_prob = [uniform_prob] * len(inv_range)
    
    len_SV_inver = inv_range.tolist()
    condition_dist_sv_inver = inv_sv_prob
    
    transmin = args.transmin
    transmax = args.transmax

    # 计算整数范围
    trans_range = np.arange(transmin, transmax + 1)
    # 计算每个整数的概率（均匀分布）
    uniform_prob_trans = 1 / len(trans_range)
    # 创建新的概率向量
    trans_sv_prob = [uniform_prob_trans] * len(trans_range)
    
    len_SV_trans = trans_range.tolist()
    condition_dist_sv_trans = trans_sv_prob
    
    len_seg_refine = len(real_con1)
     #! start variations
    ll_c = 0
    
    pai_pro_tem = [diff_ins_prob_del_real,diff_ins_prob_mis_real, diff_ins_prob_correct_real]
    pai_pro_tem_ = list(np.array(pai_pro_tem)/sum(pai_pro_tem))
    
     # Assign the arguments to the variables used in the functions
    sv_trans = args.sv_trans
    sv_inver = args.sv_inver
    True_dup_number = args.sv_dup
    sv_del = args.sv_del
    sv_ins = args.sv_ins
    mis_snv_number = int(args.snp)
    snv_del = int(args.snv_del)
    snv_ins = int(args.snv_ins)
    
    csv_params = ['num_ID1_csv', 'num_ID2_csv', 'num_ID3_csv', 'num_ID4_csv', 'num_ID5_csv', 'num_ID6_csv', 'num_ID7_csv', 'num_ID8_csv', 'num_ID9_csv', 'num_ID10_csv', 'num_ID11_csv', 'num_ID12_csv', 'num_ID13_csv', 'num_ID14_csv', 'num_ID15_csv', 'num_ID16_csv', 'num_ID17_csv', 'num_ID18_csv']

    if args.csv_num != 0:
        for param in csv_params:
            setattr(args, param, args.csv_num)
    elif args.csv_total_num !=0:
        # 内置向量
        vector = np.array([11,11,10,78,10,11,11,11,2,2,58,7,3,45,6,10,5,4])
        # 归一化向量
        normalized_vector = vector / vector.sum()
        # 计算每个csv参数的值
        csv_values = np.round(normalized_vector * args.csv_total_num).astype(int)
        # 计算差值
        diff = args.csv_total_num - csv_values.sum()
        # 将差值均匀分配到csv_values的各个元素上
        csv_values[:diff] += 1
        for i, param in enumerate(csv_params):
            setattr(args, param, csv_values[i])
            
    True_TanInvDup_number = args.num_ID1_csv
    True_DisInvDup_number = args.num_ID2_csv
    True_DisDup_number = args.num_ID3_csv
    num_Delinv_sv = args.num_ID4_csv
    num_ID5_csv = args.num_ID5_csv
    num_ID6_csv = args.num_ID6_csv
    num_ID7_csv = args.num_ID7_csv
    num_ID8_csv = args.num_ID8_csv
    num_ID9_csv = args.num_ID9_csv
    num_ID10_csv = args.num_ID10_csv
    num_delinvdel_csv = args.num_ID11_csv
    num_dupInv_csv = args.num_ID12_csv
    num_Invdup_csv = args.num_ID13_csv
    num_dupinvdup_csv = args.num_ID14_csv
    num_dupInvdel_csv = args.num_ID15_csv 
    num_delInvdup_csv = args.num_ID16_csv
    True_ID17_number = args.num_ID17_csv
    num_insdel_csv = args.num_ID18_csv
    
    def generate_gaussian(mu, sigma, num_samples=1000, scale_factor=100):
        gaussian = np.random.normal(mu, sigma, num_samples)
        gaussian = np.round(gaussian / scale_factor) * scale_factor
        gaussian = [math.ceil(abs(x)) for x in gaussian]
        return gaussian

    mu_sigma_pairs = [
        (args.mu_ID1, args.sigma_ID1), 
        (args.mu_ID2, args.sigma_ID2), 
        (args.mu_ID3, args.sigma_ID3), 
        (args.mu_ID4, args.sigma_ID4),
        (args.mu_ID5, args.sigma_ID5), 
        (args.mu_ID6, args.sigma_ID6),
        (args.mu_ID7, args.sigma_ID7), 
        (args.mu_ID8, args.sigma_ID8),
        (args.mu_ID9, args.sigma_ID9), 
        (args.mu_ID10, args.sigma_ID10),
        (args.mu_ID11, args.sigma_ID11), 
        (args.mu_ID12, args.sigma_ID12), 
        (args.mu_ID13, args.sigma_ID13), 
        (args.mu_ID14, args.sigma_ID14),
        (args.mu_ID15, args.sigma_ID15), 
        (args.mu_ID16, args.sigma_ID16),
        (args.mu_ID17, args.sigma_ID17), 
        (args.mu_ID18, args.sigma_ID18)
    ]

    gaussian_IDs = [generate_gaussian(mu, sigma) for mu, sigma in mu_sigma_pairs]

    gaussian_ID1, gaussian_ID2, gaussian_ID3, gaussian_ID4, gaussian_ID5, gaussian_ID6, gaussian_ID7, gaussian_ID8, gaussian_ID9, gaussian_ID10, gaussian_ID11, gaussian_ID12, gaussian_ID13, gaussian_ID14, gaussian_ID15, gaussian_ID16, gaussian_ID17, gaussian_ID18 = gaussian_IDs

    #record original positions and translocated positions
    SV_table = pd.DataFrame(columns=['Index','Index_con','SV_type','Original_start',\
                                        'Original_end','Len_SV','New_start','New_end','New_len_SV','Balanced Trans Flag','relative start1','relative end1','relative start2','relative end2','INFO'])

    #record location, original base (REF) and alternative base (ALT) for substitution, INDEL
    # Define a DataFrame to store the information
    VCF_table = pd.DataFrame(columns=['CHROM', 'POS', 'CSV_TYPE=ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])

    #! start variations
    # Call the functions with the assigned variables
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post = translocation(SV_table, VCF_table, unblock_region_sv, SV_loop, CSV_loop, VCF_loop, Ins_dic_sv, tem_seq_post, sv_trans, condition_dist_sv_trans, len_SV_trans, times, ll_c, len_seg_refine, ratio_b_trans, chr_id)
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post = tandem_duplication(SV_table, VCF_table,True_dup_number, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post, dup_range, times, ll_c, chr_id)
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post = inversion_process(sv_inver, SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post, chr_id, len_seg_refine,ll_c,times,real_con1, condition_dist_sv_inver, len_SV_inver)
    #! ID1: TanInvDup (Tandem Inverted Dup), 11 (73.8kb)
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post = ID1_TanInvDup_process(unblock_region_sv, True_TanInvDup_number, times, real_con1, SV_table, VCF_table, ll_c, chr_id, tem_seq_post, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv,gaussian_ID1)
    #! ID2: DisInvDup (Dispersed Inverted Dup), 11 (73.8kb)
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post = ID2_DisInvDup_process(unblock_region_sv, True_DisInvDup_number, times, real_con1, SV_table, VCF_table, ll_c, chr_id, tem_seq_post, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv,gaussian_ID2)
    #! ID3: DisDup
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv = ID3_disdup_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, True_DisDup_number, times, ll_c, chr_id, gaussian_ID3)
    #! ID4: 
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post = ID4_delinv_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_Delinv_sv, gaussian_ID4, times, ll_c, chr_id, len_seg_refine, real_con1)
    #! ID5: DEL+ DisInvDup
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post = ID5_disinvdupdel_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_ID5_csv, times, ll_c, chr_id, real_con1, condition_dist_sv_del, len_SV_del, len_seg_refine, gaussian_ID5)
    #! ID6: DEL+ DisDup
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post = ID6_disdupdel_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_ID6_csv, gaussian_ID6, times, ll_c, chr_id, real_con1, condition_dist_sv_del, len_SV_del, len_seg_refine)
    #! ID7: TanDup+DEL
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post = ID7_tandupDEL_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_ID7_csv, gaussian_ID7, times, ll_c, chr_id, real_con1, condition_dist_sv_del, len_SV_del,len_seg_refine)
    #! ID8: TanInvDup+DEL
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post = ID8_taninvdupDEL_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_ID8_csv, gaussian_ID8, times, ll_c, chr_id, real_con1, condition_dist_sv_del, len_SV_del, len_seg_refine)
    #! ID9: TanDup + DEL + INV
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post = ID9_tandup_delinv_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_ID9_csv, gaussian_ID9, times, ll_c, chr_id, real_con1, condition_dist_sv_del, len_SV_del, len_seg_refine)
    #! ID10: TanInvDup + DEL + INV
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post = ID10_taninvdup_delinv_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_ID10_csv, gaussian_ID10, times, ll_c, chr_id, real_con1, condition_dist_sv_del, len_SV_del, len_seg_refine)
    #! ID11: paired-Deletion Inversion (delInvdel)
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post = ID11_delinvdel_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_delinvdel_csv,gaussian_ID11, times, ll_c, chr_id, len_seg_refine)
    #! ID12: Inversion with 5' Flanking Duplication (dupInv)
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post = ID12_dupinv_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_dupInv_csv,gaussian_ID12, times, ll_c, chr_id, len_seg_refine)
    #! ID13: Inversion with 3' Flanking Duplication (Invdup)
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post= ID13_invdup_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_Invdup_csv,gaussian_ID13, times, ll_c, chr_id, len_seg_refine)
    #! ID14: Paired-duplication inversion (dupInvdup)    
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post = ID14_dupinvdupprocess(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_dupinvdup_csv,gaussian_ID14, times, ll_c, chr_id, len_seg_refine)
    #! ID15: Inversion with 5' Flanking Duplication and 3' Flanking Deletion (dupInvdel)
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post = ID15_dupinv_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_dupInvdel_csv,gaussian_ID15, times, ll_c, chr_id, len_seg_refine)
    #! ID16: Inversion with 5' Flanking Deletion and 3' Flanking Duplication (delInvdup)
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post = ID16_dupinv_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post, num_delInvdup_csv,gaussian_ID16, times, ll_c, chr_id, len_seg_refine)
    #! ID17: Inverted Duplication with Flanking Triplication (dupTRIPdup-INV)  
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post = ID17_dupTRIPdup_INV_process(unblock_region_sv, True_ID17_number, times,  real_con1, SV_table, VCF_table, ll_c, chr_id, tem_seq_post, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, gaussian_ID17)

    #! ID18: Insertion with Deletion (INSdel)
    #! unbalancedTranslocation + 5' Flanking deletion   
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, CSV_loop, Ins_dic_sv, tem_seq_post = ID18_INSdel_process(SV_table, VCF_table, unblock_region_sv, SV_loop,CSV_loop,gaussian_ID18, VCF_loop, Ins_dic_sv, tem_seq_post, num_insdel_csv, times, ll_c, len_seg_refine,chr_id)

            
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, tem_seq_post = long_del_process(SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop,tem_seq_post, sv_del, condition_dist_sv_del, len_SV_del, len_seg_refine, times, ll_c, chr_id)
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post = long_ins_process(unblock_region_sv,sv_ins, condition_dist_sv_ins, len_SV_ins, ins_selection, SV_table, VCF_table, ll_c, chr_id, tem_seq_post, SV_loop, VCF_loop, Ins_dic_sv)
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post = micro_del_process(snv_del, unblock_region_sv, pai_pro_tem_, times, SV_table, VCF_table, ll_c, chr_id, tem_seq_post, SV_loop, VCF_loop, Ins_dic_sv)
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post = micro_ins_process(snv_ins, unblock_region_sv, diff_ins_prob_ins_real, ins_selection, SV_table, VCF_table, ll_c, chr_id, tem_seq_post, SV_loop, VCF_loop, Ins_dic_sv)
    SV_table, VCF_table, unblock_region_sv, SV_loop, VCF_loop, Ins_dic_sv, tem_seq_post = snp_process(mis_snv_number, unblock_region_sv, base_list, substitution_matrix, mis_selection, SV_table, VCF_table, ll_c, chr_id, tem_seq_post, SV_loop, VCF_loop, Ins_dic_sv)
    
    end_time1 = time.time()

    start_time2 = time.time()
    tem_seq_post_up = tem_seq_post.copy()
    #print('tem seq post up:'+str(len(tem_seq_post_up)))
    for idx in sorted(Ins_dic_sv, reverse=True):
        #idx = 4981
        tem_seq_post_up.insert(idx+1, Ins_dic_sv[idx])
    #print(Ins_dic_sv[111211]) 
 
    tem_seq_post_up_string = ''.join(tem_seq_post_up)
    tem_seq_post_up_string= tem_seq_post_up_string.replace('-','')
    updated_con.append(copy.deepcopy(tem_seq_post_up_string))
    print('Length of the simulated sequence: '+str(len(tem_seq_post_up_string)))
    # 删除第一列
    SV_table_merged = copy.deepcopy(SV_table)
    tem_ins_dic = Ins_dic_sv

    # # Save the dictionary as a .npy file
    # np.save(save+str(rep)+'_tem_ins_dic.npy', tem_ins_dic)
    #! final table
    if args.write:
        print('finalize table')
        SV_table_merged = CSV_finalize_table(SV_table_merged,ll_c, tem_ins_dic)
        SV_table_merged.to_csv(args.save +'BV_' + str(args.rep) + '_seq_' + str(seqname) + '_SVtable_full.csv', header=True, index=False)
   
    else:
        SV_table_merged.to_csv(args.save +'BV_' + str(args.rep) + '_seq_' + str(seqname) + '_SVtable.csv', header=True, index=False)
        # Save the dictionary as a .npy file
        np.save(args.save+'BV_'+str(args.rep)+ '_seq_'+str(seqname)+'_tem_ins_dic.npy', tem_ins_dic)
    
    # 按照 'pos' 列排序
    VCF_table_merged = copy.deepcopy(VCF_table)
    VCF_table_merged.sort_values(by='POS', inplace=True)

    # 重置索引并将旧索引添加为 'Index' 列
    VCF_table_merged.reset_index(inplace=True, drop=True)

    # 更新 'CSV_TYPE=ID' 列
    VCF_table_merged['CSV_TYPE=ID'] = 'rs' + VCF_table_merged.index.astype(str)
            
    def write_template_fasta_con(args, seqname, consensus_):
        # Prepare the new sequence
        sequences = [consensus_]
        new_sequences = []
        for sequence in sequences:
            record = SeqRecord(Seq(re.sub('[^GATCN-]', "", str(sequence).upper())), id=seqname, name=seqname, description="<custom description>")
            new_sequences.append(record)

        # Write the new sequence to a file
        with open(args.save + 'BV_' + str(args.rep) + "_seq_"+str(seqname) +".fasta", "w") as output_handle:
            SeqIO.write(new_sequences, output_handle, "fasta")

    ############vcf
        
    # def write_vcf(df, save, rep, con_id, ref, seqname, start_base, end_base):
    def write_vcf(args, df, seqname, start_base, end_base):
        # Get the current date
        current_date = datetime.now().strftime('%Y%m%d')
        # Add the additional columns to the DataFrame
        df['FORMAT'] = 'GT'
        df['SAMPLE_ID'] = '1/1'
        # Write the DataFrame to a VCF file
        with open(args.save +'BV_' + str(args.rep) + '_seq_' + str(seqname) +".vcf", 'w') as f:
            f.write('##fileformat=VCFv4.2\n')
            f.write('##fileDate=' + current_date + '\n')
            f.write('##source=uniform.py\n')
            f.write('##reference=' + args.ref + ':' + str(start_base) + '-' + str(end_base) + '\n')
            f.write('##contig=<ID='+str(seqname)+',length=' + str(end_base - start_base + 1) + '>\n')
            f.write('##source=CSV.py\n')
            f.write('##reference=' + ref + ':' + str(start_base) + '-' + str(end_base) + '\n')
            f.write('##contig=<ID='+str(seqname)+',length=' + str(end_base - start_base + 1) + '>\n')
            f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
            f.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the variant">\n')
            f.write('##INFO=<ID=CSV_TYPE,Number=1,Type=String,Description="Complex SV type">\n')
            f.write('##INFO=<ID=CSV_INDEX,Number=1,Type=Integer,Description="Index of the generated complex SV">\n')
            f.write('##INFO=<ID=ID1,Number=1,Type=String,Description="TanInvDup (Tandem Inverted Dup)">\n')
            f.write('##INFO=<ID=ID2,Number=1,Type=String,Description="DisInvDup (Dispersed Inverted Dup)">\n')
            f.write('##INFO=<ID=ID3,Number=1,Type=String,Description="DisDup">\n')
            f.write('##INFO=<ID=ID4,Number=1,Type=String,Description="DEL + INV">\n')
            f.write('##INFO=<ID=ID5,Number=1,Type=String,Description="DEL + DisInvDup">\n')
            f.write('##INFO=<ID=ID6,Number=1,Type=String,Description="DEL + DisDup">\n')
            f.write('##INFO=<ID=ID7,Number=1,Type=String,Description="TanDup + DEL">\n')
            f.write('##INFO=<ID=ID8,Number=1,Type=String,Description="TanInvDup + DEL">\n')
            f.write('##INFO=<ID=ID9,Number=1,Type=String,Description="TanDup + DEL + INV">\n')
            f.write('##INFO=<ID=ID10,Number=1,Type=String,Description="TanInvDup + DEL + INV">\n')
            f.write('##INFO=<ID=ID11,Number=1,Type=String,Description="DEL + DEL + INV">\n')
            f.write('##INFO=<ID=ID12,Number=1,Type=String,Description="DUP + INV">\n')
            f.write('##INFO=<ID=ID13,Number=1,Type=String,Description="INV + DUP">\n')
            f.write('##INFO=<ID=ID14,Number=1,Type=String,Description="DUP + INV + DUP">\n')
            f.write('##INFO=<ID=ID15,Number=1,Type=String,Description="DUP + INV + DEL">\n')
            f.write('##INFO=<ID=ID16,Number=1,Type=String,Description="DEL + INV + DUP">\n')
            f.write('##INFO=<ID=ID17,Number=1,Type=String,Description="DUPTRIPDUP + INV">\n')
            f.write('##INFO=<ID=ID18,Number=1,Type=String,Description="DEL + unbalancedTrans">\n')
            
            f.write('##INFO=<ID=SUB,Number=1,Type=String,Description="Substitution">\n')
            f.write('##INFO=<ID=microINS,Number=1,Type=String,Description="Micro insertion">\n')
            f.write('##INFO=<ID=LEN,Number=1,Type=Integer,Description="Length of the variant">\n')
            f.write('##INFO=<ID=microDEL,Number=1,Type=String,Description="Micro deletion">\n')
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_ID\n')
            df.to_csv(f, sep='\t', index=False, header=False)
        
            

            
            
    #def write_template_fasta_con(save, seqname, consensus_, rep, con_id):

    write_template_fasta_con(args, seqname, updated_con[0])
    
    
    #def write_vcf(df, save, rep, con_id, ref, start_base, end_base):
    #write_vcf(VCF_table_merged, save,rep,ll_c,ref,seqname,start_base, end_base)
    write_vcf(args, VCF_table_merged, seqname, start_base, end_base)
    end_time2 = time.time()

    elapsed_time1 = end_time1 - start_time1
    #formatted_time1 = str(timedelta(seconds=elapsed_time1))

    #print(f"Trans,INV,DUP运行时间：{formatted_time1}")

    elapsed_time2 = end_time2 - start_time2
    #formatted_time2 = str(timedelta(seconds=elapsed_time2))

    #print(f"写出结果运行时间：{formatted_time2}")
    
    elapsed_time3 = end_time2 - start_time1
    formatted_time3 = str(timedelta(seconds=elapsed_time3))

    print(f"total time：{formatted_time3}")
    # 记得在结束时关闭文件


    sys.stdout.close()
    
if __name__ == "__main__":
    main()
