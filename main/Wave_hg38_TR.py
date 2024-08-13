#!/usr/bin/env python
# coding: utf-8
import time
start_time0 = time.time()
process_time0 = time.process_time()
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

print('hg38 mode')
# 定义一个函数，该函数接收一行数据，检查'start'值是否小于'end'值
def check_start_end(row):
    if row['start'] >= row['end']:
        print('Warning: The "start" value of the .bed file is greater than or equal to the "end" value.')
        return False
    return True

# 定义一个函数，该函数接收一行数据，返回该行'start'和'end'区域中的所有点
def get_points(row):
    return set(range(row['start'], row['end']))
    
def error_handler(e):
    print('Error occurred:', e)
    
def SV_write_relative(SV_table_merged,ll_c,tem_ins_dic):
    tem_SV_table_merged = SV_table_merged[SV_table_merged.iloc[:,1]==ll_c]
    #original start: start of A
    list_start1 = list(tem_SV_table_merged.iloc[:,3])
    #start of B (dulplication or balanced trans)
    list_start2 = list(tem_SV_table_merged.iloc[:,6])
    whole_start_abs = list_start1+list_start2
    #order SV from left
    #set:merge repeated sites (e.g. 5 mismatch 5.5 ins)
    #! 筛选出大于-1的
    whole_start_abs = [item for item in whole_start_abs if item > 0]
    whole_start_abs_set = sorted(list(set(whole_start_abs)))
    present_len = 0
    last_bone = 0
    #inserted term
    
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
                        
                elif tem_row[2] in ['Duplication']:
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
                
                if tem_row[2] in ['Duplication']:
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
    
#def gen_consensus(xun, process_dict):
def gen_consensus(xun, chr_id,ll_c, starts_seg, ends_seg, tem_seq_post, del_SV_per_segment, ins_SV_per_segment, del_snv_per_segment,ins_snv_per_segment,snp_per_segment, pai_pro_tem_, args, condition_dist_sv_del, len_SV_del,condition_dist_sv_ins,len_SV_ins,ins_selection,base_list,substitution_matrix,mis_selection,times, TR_set, p_del_TR, p_ins_TR,tmp_dir):   

    SV_loop_seg = 0
    VCF_loop_seg = 0
    
    #copy the tem_seq_post with already the deleted part from Trans, Inv and DUP
    #!
    start_seg = starts_seg[xun]
    end_seg = ends_seg[xun]
    tem_seq_post_seg = copy.deepcopy(tem_seq_post[start_seg:end_seg])
    print('length of seq'+str(len(tem_seq_post_seg)))
    #dictionary, restore ins sites and corresponding content
    Ins_dic_sv_seg = {}

    # 从临时文件加载 unblock_region_vec[xun]
    unblock_region_vec_xun = np.load(os.path.join(tmp_dir, 'unblock_region_vec_{}.npy'.format(xun)))


    unblock_region_seg = list(np.copy(unblock_region_vec_xun))
    
    
    #! long del
    print(str(xun)+'long del')
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
    SV_table_seg = np.empty((int(max_length_numpy*2), 14), dtype=object)
    VCF_table_seg = np.empty((int(max_length_numpy*2), 8), dtype=object)

    
    circular_count_del_break = 0 
    
    # 设置随机数种子使得结果可以复现
    # 为每个进程设置一个基于输入的种子
    seed_seg = args.seed + xun
    random.seed(seed_seg)
    np.random.seed(seed_seg)
    
    #!
    if not unblock_region_seg:# 如果segment是空集，输出警告
        print("Warning: no available positions for variations in no."+str(xun)+"segment")
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
        print("Warning: no available positions for long deletions in no."+str(xun)+"segment")
        circular_count_del_break = 1
    elif True_del_number == 0:
        circular_count_del_break = 1
    else:
        circular_count_del_break = 0
        len_seg_refine2 = max(unblock_region_seg)
        #print('SV DEL number:'+str(True_del_number))

    if not circular_count_del_break:
        print('SV DEL number:'+str(True_del_number))
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

                SV_table_seg[SV_loop_seg] = [SV_loop_seg,ll_c,'Deletion',r_s,r_s+l_s-1,l_s,-1,-1,-1,-1,0,0,0,0]
                SV_loop_seg = SV_loop_seg + 1 

                VCF_table_seg[VCF_loop_seg] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop_seg), ''.join(tem_seq_post_seg[r_s-start_seg:r_s-start_seg+l_s]) + tem_seq_post_seg[r_s-start_seg+l_s], tem_seq_post_seg[r_s-start_seg+l_s], '.', 'PASS', 'SVTYPE=DEL;SVLEN='+str(l_s)]

                VCF_loop_seg = VCF_loop_seg + 1  

                tem_seq_post_seg[r_s-start_seg:(r_s-start_seg+l_s)] = '-'*l_s          
            else:
                continue
        
    
    #! long insertion
    print('long ins')

    
    whole_insertion_term = []

    circular_count_ins_break = 0
    ins_pos_collection_seg = []
    True_ins_number = ins_SV_per_segment[xun]

    if not unblock_region_seg:
        print("Warning: no available positions for long insertions in no."+str(xun)+"segment")
        circular_count_ins_break = 1
    elif True_ins_number == 0:
        circular_count_ins_break = 1
    else:
        circular_count_ins_break = 0

    if not circular_count_ins_break:
        print('SV INS:'+str(True_ins_number))
        if not unblock_region_seg:
            print("Warning: no available positions for long insertions in no."+str(xun)+"segment")
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

            SV_table_seg[SV_loop_seg] = [SV_loop_seg,ll_c,'Insertion',remain_index2,remain_index2,l_i,-1,-1,-1,-1,0,0,0,0]
            SV_loop_seg = SV_loop_seg + 1

            VCF_table_seg[VCF_loop_seg] = [str(chr_id), str(remain_index2), 'rs' + str(VCF_loop_seg), tem_seq_post_seg[remain_index2-start_seg], tem_seq_post_seg[remain_index2-start_seg] + tem_ins.upper(), '.', 'PASS', 'SVTYPE=INS;SVLEN='+str(l_i)]
            VCF_loop_seg = VCF_loop_seg + 1

    
    #!

    #! micro del

    ### finish the process of the duplication and insetion
    ## possible var bone
    #snp_considered_sites = copy.deepcopy(unblock_region_seg)
    # unblock_region_snv= copy.deepcopy(unblock_region_seg)#deletions
    # undel_region_snv=copy.deepcopy(undel_region_seg)#substitutions
    # left_del_region_snv=copy.deepcopy(left_del_region_seg)#insertions
    #!!
    # 对每一段进行处理
    if not unblock_region_seg:# 如果segment是空集，输出警告
        print("Warning: no available positions for micro deletions in no."+str(xun)+"segment")
        circular_count_micro_del_break = 1
    else:
        circular_count_micro_del_break = 0
        
    # 计算需要选择的位点的数量
    if not circular_count_micro_del_break:
        ## micro deletion part
        len_unblock_region = len(unblock_region_seg)
        l_s_vec_ini = np.random.multinomial(n=len_unblock_region, pvals=pai_pro_tem_)
        
        max_del = len_unblock_region // 2

        if args.snv_del is None:
            del_snv_number = int(l_s_vec_ini[0])
            #del_snv_number = test_number
        else:
            del_snv_number = min(True_del_snv_number, max_del)
            if True_del_snv_number > max_del:
                print("Warning: The input for -snv_del is too large and has been automatically reduced. \
                    This is because each micro deletion event requires at least one position, and there must be at least one position between two micro deletion events. \
                        Therefore, the maximum number of micro deletions that can occur is half of the total number of positions available for events.")
    
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
                if circular_count_del> args.args.times:
                    circular_count_del_break = 1
                    print("Warning: No."+str(m_del)+ "  micro deletion sampling times exceeds "+str(times)+' times. Try: reduce number of variations, increase times or reference length.')
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
            
                SV_table_seg[SV_loop_seg] = [SV_loop_seg,ll_c,'Micro_Del',r_s,r_s+l_s-1,l_s,-1,-1,-1,-1,0,0,0,0]
                SV_loop_seg = SV_loop_seg + 1

                #The ‘REF’ column is set to the original segment plus the base that is left after the deletion. 
                #The ‘ALT’ column is set to the base that is left after the deletion.
                # Add a row to the VCF_table for the micro deletion
                #VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), tem_seq_post[r_s:r_s+l_s] + tem_seq_post[r_s+l_s], tem_seq_post[r_s+l_s], '.', 'PASS', 'SVTYPE=microDEL']
                VCF_table_seg[VCF_loop_seg] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop_seg), ''.join(tem_seq_post_seg[r_s-start_seg:r_s-start_seg+l_s]) + tem_seq_post_seg[r_s-start_seg+l_s], tem_seq_post_seg[r_s-start_seg+l_s], '.', 'PASS', 'microDEL;LEN='+str(l_s)]

                VCF_loop_seg = VCF_loop_seg + 1  

                #replace deleted bases with -
                tem_seq_post_seg[r_s-start_seg:(r_s-start_seg+l_s)] = '-'*l_s
            else:
                break
                
        
    ### micro insertions
    # 对每一段进行处理
    #if not left_del_region_snv:# 如果segment是空集，输出警告
    if not unblock_region_seg:# 如果segment是空集，输出警告
        print("Warning: no available positions for micro insertions in no."+str(xun)+"segment")
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
            SV_table_seg[SV_loop_seg] = [SV_loop_seg,ll_c,'Micro_Ins',remain_index2,remain_index2,l_i,-1,-1,-1,-1,0,0,0,0]
            SV_loop_seg = SV_loop_seg + 1

            # Add a row to the VCF_table for the micro insertion
            VCF_table_seg[VCF_loop_seg] = [str(chr_id), str(remain_index2), 'rs' + str(VCF_loop_seg), tem_seq_post_seg[remain_index2-start_seg], tem_seq_post_seg[remain_index2-start_seg] + tem_ins.upper(), '.', 'PASS', 'microINS;LEN='+str(l_i)]
            VCF_loop_seg = VCF_loop_seg + 1  

        #snp=int(l_s_vec_ini[1])
    
    #! start SNP
    ### micro insertions
    # 对每一段进行处理
    #if not undel_region_snv:# 如果segment是空集，输出警告
    if not unblock_region_seg:# 如果segment是空集，输出警告
        print("Warning: no available positions for substitutions in no."+str(xun)+"segment")
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

                SV_table_seg[SV_loop_seg] = [SV_loop_seg,ll_c,'Substitution',ll,ll,1,-1,-1,-1,-1,0,0,0,0]
                SV_loop_seg = SV_loop_seg + 1
                # Add a row to the VCF_table for the substitution
                VCF_table_seg[VCF_loop_seg] = [str(chr_id), str(ll), 'rs' + str(VCF_loop_seg), tem_seq_post_seg[ll-start_seg], bexixuan_.upper(), '.', 'PASS', 'SUB']
                VCF_loop_seg = VCF_loop_seg + 1

                tem_seq_post_seg[int(ll-start_seg)] = copy.deepcopy(bexixuan_)
            else:
                print("Error: Invalid base")
        
    # 找出所有元素都是 None 或空字符串的行
    mask1 = np.all((SV_table_seg == None) | (SV_table_seg == ''), axis=1)

    # 删除这些行
    SV_table_seg = SV_table_seg[~mask1]

    # 现在你可以将 SV_table 转换为 DataFrame
    SV_table_seg = pd.DataFrame(SV_table_seg, columns=['Index','Index_con','SV_type','Original_start',\
                                        'Original_end','Len_SV','New_start','New_end','New_len_SV','Balanced Trans Flag','relative start1','relative end1','relative start2','relative end2'])

    # 对 VCF_table 重复相同的步骤
    mask2 = np.all((VCF_table_seg == None) | (VCF_table_seg == ''), axis=1)
    VCF_table_seg = VCF_table_seg[~mask2]
    VCF_table_seg = pd.DataFrame(VCF_table_seg, columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
    
    print(str(xun)+':save files')
    return SV_table_seg, VCF_table_seg, unblock_region_seg, Ins_dic_sv_seg, tem_seq_post_seg

def DNA_complement(sequence):
    trantab = str.maketrans('ATCGatcg','TAGCtagc')
    string = sequence.translate(trantab)
    return string

# 定义一个全局变量来存储最大的磁盘使用率
def monitor_disk_usage(max_disk_usage):
    while True:
        disk_usage = psutil.disk_usage('/').percent
        with max_disk_usage.get_lock():
            if disk_usage > max_disk_usage.value:
                max_disk_usage.value = disk_usage
        time.sleep(1)

# Monitor memory usage
def monitor_memory(threshold, max_mem_usage, max_mem_usage_gb,tmp_dir):
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
    parser = argparse.ArgumentParser(description='GenoWave')
    parser.add_argument('-ref', type=str, help='Input reference local path', default='default_ref')
    parser.add_argument('-save', type=str, help='local path for saving', default=os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'save')+ '/')
    parser.add_argument('-seed', type=int, help='Seed for random number generator', default=999)
    parser.add_argument('-times', type=int, help='Number of times', default=10)
    parser.add_argument('-rep', type=int, help='Replication ID', default=1)
    parser.add_argument('-sv_trans', type=int, help='Number of trans SV', default=2)
    parser.add_argument('-sv_inver', type=int, help='Number of inversion SV', default=3)
    parser.add_argument('-sv_dup', type=int, help='True duplication number', default=4)
    parser.add_argument('-sv_del', type=int, help='Number of deletion SV', default=8)
    parser.add_argument('-sv_ins', type=int, help='True insertion number', default=8)
    parser.add_argument('-snp', type=float, help='SNV number or probability', default=None)
    parser.add_argument('-snv_del', type=float, help='SNV deletion number or probability', default=None)
    parser.add_argument('-snv_ins', type=float, help='SNV insertion number or probability', default=None)
    parser.add_argument('-notblockN', action='store_true', help='Do not Block N positions')
    parser.add_argument('-write', action='store_true', help='Write full results')
    parser.add_argument('-block_region_bed_url', type=str, help='local path of the block region BED file', default=None)
    parser.add_argument('-cores', type=int, help='Number of kernels for parallel processing', default=1)
    parser.add_argument('-len_bins', type=int, help='Length of bins for parallel processing', default=50000)
    parser.add_argument('-delmin', type=int, help='Minimum deletion length', default=50)
    parser.add_argument('-delmax', type=int, help='Maximum deletion length', default=60)
    parser.add_argument('-insmin', type=int, help='Minimum insertion length', default=50)
    parser.add_argument('-insmax', type=int, help='Maximum insertion length', default=450)
    parser.add_argument('-mode', type=str, help='Mode for calculating probabilities', default='probability')
    parser.add_argument('-sum', action='store_true', help='total indel SV equals sum of the input bed')
    #parser.add_argument('-indel_input_bed', '--indel_input_bed', type=str, help='Input BED file for indels', default='~/data/test_data/TGS/hg002/chr21_SV_Tier1.bed')
    parser.add_argument('-p_del_region', type=float, help='Probability of SV DEL in the user-defined region for deletion', default=0.5)
    parser.add_argument('-p_ins_region', type=float, help='Probability of SV INS in the user-defined region for insertion', default=0.5)
    parser.add_argument('-region_bed_url', type=str, help='local path of the BED file for the user-defined region', default='~/data/test_data/TGS/hg002/chr21_TR_unique.bed')
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
    # Create a temporary directory if not exists
    # if not os.path.exists('tmp'):
    #     os.makedirs('tmp')

    # Create a temporary directory if not exists under args.save
    tmp_dir = os.path.join(args.save, 'tmp')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    print('Trans:'+str(args.sv_trans))
    print('Inversion:'+str(args.sv_inver))
    print('DUP:'+str(args.sv_dup))
    print('DEL:'+str(args.sv_del))
    print('INS:'+str(args.sv_ins))

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
    print('reference\'s name:'+str(seqname))
    # start_base=43531239
    # end_base=43601240
    #len:46709983
    start_base=0
    end_base=len(BestRefSeq)
    real_con1 = copy.deepcopy(BestRefSeq[start_base:end_base+1])
    chr_length = len(real_con1)
    del BestRefSeq
    print('Length of ref:'+str(chr_length))

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

    print('number_seg:'+str(number_seg))

    # 从.npy文件中读取数据
    len_SV_del = np.load(main_url_empirical+'len_SV_del_cell.npy').tolist()
    condition_dist_sv_del = np.load(main_url_empirical+'condition_dist_sv_del_cell.npy').tolist()

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
    len_SV_ins = np.load(main_url_empirical+'len_SV_ins_cell.npy').tolist()
    condition_dist_sv_ins = np.load(main_url_empirical+'condition_dist_sv_ins_cell.npy').tolist()

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
    #record original positions and translocated positions
    SV_table = pd.DataFrame(columns=['Index','Index_con','SV_type','Original_start',\
                                        'Original_end','Len_SV','New_start','New_end','New_len_SV','Balanced Trans Flag','relative start1','relative end1','relative start2','relative end2'])

    #record location, original base (REF) and alternative base (ALT) for substitution, INDEL
    # Define a DataFrame to store the information
    VCF_table = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
    
    
    ll_c=0
    len_seg_refine = len(real_con1)
    
    #deletion probability
    del_sv_createria = 8.9/9
    trans_sv_createria = 8.9/9
    inver_sv_createria = 8.9/9
    ins_location_p = 1/90000
    dup_location_p = 1/90000


    SV_loop = 0
    VCF_loop = 0
    CSV_loop = 0
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

    pdel_SV = 0.01
    #len_SV_del = [50,100,500]
    len_SV_trans = [50,100,500]
    
    #len_SV_ins = [50,100,500]
    #condition_dist_sv_del = [4/5,1/10,1/10]
    condition_dist_sv_trans = [4/5,1/10,1/10]
    
    len_SV_inver = [50,100,500]
    condition_dist_sv_inver = [4/5,1/10,1/10]
    #condition_dist_sv_ins = [4/5,1/10,1/10]
    number_sv_del = [4/5,1/10,1/10]
    number_sv_trans = [4/5,1/10,1/10]
    number_sv_inver = [4/5,1/10,1/10]
    number_sv_dup = [4/5,1/10,1/10]
    # copied_base_sv_prob = [4/5,1/10,1/10]
    # copied_base_sv_base = [50,100,500]
    # 设置参数
    dupmin = 50
    dupmax = 450

    # 计算整数范围
    dup_range = np.arange(dupmin, dupmax + 1)
    # 计算每个整数的概率（均匀分布）
    # uniform_prob = 1 / len(dup_range)
    # # 创建新的概率向量
    # dup_sv_prob = [uniform_prob] * len(dup_range)
    
    # copied_base_sv_base = dup_range.tolist()
    # copied_base_sv_prob = dup_sv_prob
    
    # 设置参数
    invmin = 50
    invmax = 450

    # 计算整数范围
    inv_range = np.arange(invmin, invmax + 1)
    # 计算每个整数的概率（均匀分布）
    uniform_prob = 1 / len(inv_range)
    # 创建新的概率向量
    inv_sv_prob = [uniform_prob] * len(inv_range)
    
    len_SV_inver = inv_range.tolist()
    condition_dist_sv_inver = inv_sv_prob
    
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
                SV_table.loc[SV_loop] = [SV_loop,ll_c,'Translocation',r_s,r_s+l_s-1,l_s,ins_trans_loc,ins_trans_loc,l_s,0,0,0,0,0]### cut and paste
                SV_loop = SV_loop + 1
                
                VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), ''.join(tem_seq_post[r_s:r_s+l_s]) + tem_seq_post[r_s+l_s], tem_seq_post[r_s+l_s], '.', 'PASS', 'SVTYPE=DEL;SVLEN='+str(l_s)+';CSV_TYPE=unbalanedTrans;CSV_INDEX='+str(CSV_loop)]
                VCF_loop = VCF_loop+ 1

                VCF_table.loc[VCF_loop] = [str(chr_id), str(ins_trans_loc), 'rs' + str(VCF_loop), tem_seq_post[ins_trans_loc], tem_seq_post[ins_trans_loc] + ''.join(tem_seq_post[r_s:(r_s+l_s)]), '.', 'PASS', 'SVTYPE=INS;SVLEN='+str(l_s)+';CSV_TYPE=unbalanedTrans;CSV_INDEX='+str(CSV_loop)]
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
                    SV_table.loc[SV_loop] = [SV_loop,ll_c,'Translocation',r_s,r_s+l_s-1,l_s,r_s2,r_s2+l_s2-1,l_s2,1,0,0,0,0]### copy and paste
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

            SV_table.loc[SV_loop] = [SV_loop,ll_c,'Inversion',r_s,r_s+l_s-1,l_s,r_s,r_s+l_s-1,l_s,-1,0,0,0,0]
            SV_loop = SV_loop + 1
            
            VCF_table.loc[VCF_loop] = [str(chr_id), str(r_s), 'rs' + str(VCF_loop), ''.join(original_string), ''.join(original_string_reverse), '.', 'PASS', 'SVTYPE=INV;SVLEN='+str(l_s)]
            VCF_loop= VCF_loop + 1
        else:
            break
        
    #! Duplication
    #args.sv_dup= 2
    print('DUP:'+str(args.sv_dup))
    # 打印所有被选择的位点的信息
    inserted_site_collection = location_insert_trans

    #inserted_site_collection = sample(left_del_region_sv,full_number_base)
    #all_selected_dup_SV= sample(left_del_region_sv,args.sv_dup)
    all_selected_dup_SV= sample(unblock_region_sv,args.sv_dup)
    # 打印所有被选择的位点的信息

    dup_site_collection=all_selected_dup_SV
    #print('Dup sites'+str(dup_site_collection))     

    # other_sites = list(set(left_del_region_sv)-set(inserted_site_collection)-set(dup_site_collection))
    other_sites = list(set(unblock_region_sv)-set(inserted_site_collection)-set(dup_site_collection))
    #duplication: copy and paste
    for remain_index2 in dup_site_collection:##others### questions
        #print(remain_index2)
        #number of copies
        # dup_num_index = np.random.multinomial(n=1, pvals=number_sv_dup)
        # dup_num = int(list(dup_num_index).index(1))+1
        # # length of duplication
        # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[1]
        # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
        tem_copied_base = np.random.choice(dup_range)
        ##if all bases that are copied are in left_del_region_sv, then no resampling is needed
        ## number of the common elements of two sets
        # if len(set(range(remain_index2-tem_copied_base+1,remain_index2+1))&set(left_del_region_sv)) < tem_copied_base:
        #     repeat_falg_tran_ins4 = 1

        circular_count_dup = 0
        circular_count_dup_break = 0
        #while len(set(range(remain_index2-tem_copied_base+1,remain_index2+1))&set(left_del_region_sv)) < tem_copied_base:
        while len(set(range(remain_index2-tem_copied_base+1,remain_index2+1))&set(unblock_region_sv)) < tem_copied_base:
            # select the possibile copied length 
            #sites that do note have general insertions

            remain_index2 = sample(other_sites,1)[0]
            # tem_copied_base_index = (np.random.multinomial(n=1, pvals=copied_base_sv_prob))[1]
            # tem_copied_base = int(copied_base_sv_base[int(tem_copied_base_index)])
            tem_copied_base = np.random.choice(dup_range)
            # l_s_vec = np.random.multinomial(n=1, pvals=condition_dist_sv_inver)
            # l_s_index = list(l_s_vec).index(1)
            # l_s = len_SV_inver[int(l_s_index)]
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
                    
                location_insert_dup.append(remain_index2)

                #remain_index2 is a new ins position, later cannot ins here
                # if remain_index2 in left_del_region_sv:
                #     left_del_region_sv.remove(remain_index2)
                Ins_dic_sv[remain_index2] = ''.join(tem_seq_post[remain_index2-tem_copied_base+1:remain_index2+1])
                ins_len_dup = len(Ins_dic_sv[remain_index2])
                SV_table.loc[SV_loop] = [SV_loop,ll_c,'Duplication',remain_index2-tem_copied_base+1,remain_index2,tem_copied_base,remain_index2,remain_index2,ins_len_dup,-1,0,0,0,0]
                SV_loop = SV_loop + 1
                
                VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index2), 'rs' + str(VCF_loop), tem_seq_post[remain_index2], tem_seq_post[remain_index2]+''.join(tem_seq_post[remain_index2-tem_copied_base+1:remain_index2+1]), '.', 'PASS', 'SVTYPE=DUP;SVLEN='+str(tem_copied_base)]
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

                SV_table.loc[SV_loop] = [SV_loop,ll_c,'Duplication',remain_index2-tem_copied_base+1,remain_index2,tem_copied_base,remain_index22,remain_index22,ins_len_dup,-1,0,0,0,0]
                SV_loop = SV_loop + 1
                
                VCF_table.loc[VCF_loop] = [str(chr_id), str(remain_index22), 'rs' + str(VCF_loop), tem_seq_post[remain_index22], tem_seq_post[remain_index22]+''.join(tem_seq_post[remain_index2-tem_copied_base+1:remain_index2+1]), '.', 'PASS', 'SVTYPE=DUP;SVLEN=CSV_TYPE=DisDup;CSV_INDEX='+str(CSV_loop)]
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

    process_end_time0 = time.process_time()

    elapsed_time0 = end_time0 - start_time0
    formatted_time0 = str(timedelta(seconds=elapsed_time0))

    process_formatted_time0 = process_end_time0 - process_time0
    process_0_time = str(timedelta(seconds=process_formatted_time0))

    print(f"Max disk usage during writing parameters was: {max_disk_usage.value}%")
    max_disk_usage.value = 0.0  # 重置最大磁盘使用率
    print(f"Initialization,Trans等待时间：{formatted_time0}")
    print(f"Initialization, Trans执行实际时间：{process_0_time}")

    start_time1 = time.time()
    process_time1 = time.process_time()

    # 初始化一个空的numpy数组来存储交集
    unblock_region_vec = np.empty(number_seg, dtype=object)
    # 将 Python 列表转换为 numpy 数组
    unblock_region_sv = np.array(unblock_region_sv)


    # 确保 numpy 数组已排序
    unblock_region_sv.sort()

    unblock_vec_lengths = []

    for i, segment in enumerate(segments_initial):
        unblock_region_vec[i] = segment[np.in1d(segment, unblock_region_sv)]

        # 记录 unblock_region_vec[i] 的长度
        unblock_vec_lengths.append(len(unblock_region_vec[i]))

        # 将 unblock_region_vec[i] 存储到临时文件
        np.save(os.path.join(tmp_dir, 'unblock_region_vec_{}.npy'.format(i)), unblock_region_vec[i])

    seg_probabilities = [p_length/sum(unblock_vec_lengths) for p_length in unblock_vec_lengths]
    
    mode = args.mode # or probability weights
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

    
    ins_empirical = INS_sampled_values

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

    
    del_empirical = DEL_sampled_values
    
    # 对概率进行归一化
    #args.sum
    if len(del_empirical) == number_seg:
        if mode == 'probability':
            del_probabilities = [p_del/sum(del_empirical) for p_del in del_empirical]
            if args.sum:
                total_del_sv = sum(del_empirical)
            else:
                total_del_sv = args.sv_del
            print('Total del:'+str(total_del_sv))
            del_SV_per_segment = np.random.multinomial(total_del_sv, del_probabilities)

        elif mode == 'empirical':
            total_del_sv = sum(del_empirical)
            del_SV_per_segment = del_empirical
            print('Total del:'+str(total_del_sv))
        else:
            print('MODE TYPE ERROR')
    else:
        print('warning: empirical del info wrong dimension')

    if len(ins_empirical) == number_seg:
        if mode == 'probability':
            ins_probabilities = [p_ins/sum(ins_empirical) for p_ins in ins_empirical]
            if args.sum:
                total_ins_sv = sum(ins_empirical)
            else:
                total_ins_sv = args.sv_ins
            print('Total ins:'+str(total_ins_sv))
            ins_SV_per_segment = np.random.multinomial(total_ins_sv, ins_probabilities)
        elif mode == 'empirical':
            total_ins_sv = sum(ins_empirical)
            ins_SV_per_segment = ins_empirical
        else:
            print('MODE TYPE ERROR')
    else:
        print('warning: empirical ins info wrong dimension')


    del ins_empirical, del_empirical
        
    # del_SV_per_segment = np.random.multinomial(args.sv_del, seg_probabilities)
    # ins_SV_per_segment = np.random.multinomial(args.sv_ins, seg_probabilities)
    print('DEL counts:'+str(del_SV_per_segment))
    print('INS counts:'+str(ins_SV_per_segment))

    if args.snv_del is not None:
        del_snv_per_segment = np.random.multinomial(args.snv_del, seg_probabilities)
    
    if args.snv_ins is not None:
        ins_snv_per_segment = np.random.multinomial(args.snv_ins, seg_probabilities)     
        
    if args.snp is not None:
        snp_per_segment = np.random.multinomial(args.snp, seg_probabilities)  
        
    # 定义概率
    #! TR region
    p_del_TR = args.p_del_region
    p_ins_TR = args.p_ins_region
    
    #TR_region_bed = '~/data/test_data/TGS/hg002/chr21_TR_unique.bed'
    TR_region_bed = args.region_bed_url
    # 读取BED文件
    df_TR = pd.read_csv(TR_region_bed, sep='\t', header=None, names=['start', 'end'])

    

    # 使用.apply()函数应用上面定义的函数，得到每一行'start'和'end'值的检查结果
    df_TR['check'] = df_TR.apply(check_start_end, axis=1)

    # 检查是否所有行的'start'值都小于'end'值
    if df_TR['check'].all():
        print("All rows are valid.")
    else:
        print("Error: Not all rows are valid.")
        sys.exit(1)


    # 使用.apply()函数应用上面定义的函数，得到每一行'start'和'end'区域中的所有点
    df_TR['points'] = df_TR.apply(get_points, axis=1)

    # 使用set.union()函数合并所有的点，得到一个包含所有点的集合
    TR_set = set.union(*df_TR['points'])

    # 现在，'TR_set'集合中包含了所有的点
    print('len of TR region'+str(len(TR_set)))
            
    end_time1 = time.time()

    process_end_time1 = time.process_time()

    elapsed_time1 = end_time1 - start_time1
    formatted_time1 = str(timedelta(seconds=elapsed_time1))

    process_formatted_time1 = process_end_time1 - process_time1
    process_first_time = str(timedelta(seconds=process_formatted_time1))

    print(f"Max disk usage during writing parameters was: {max_disk_usage.value}%")
    max_disk_usage.value = 0.0  # 重置最大磁盘使用率
    print(f"取交集 等待时间：{formatted_time1}")
    print(f"取交集 CPU执行实际时间：{process_first_time}")

    process_time2 = time.process_time()
    start_time2 = time.time()

    #! delete variables to release 内存
    # del BestRefSeq, n_positions, n_positions_set, all_positions, all_region_sv
    # del segments_initial, all_region_sv_array
    # del unblock_region_sv, undel_region_sv, left_del_region_sv
    #! define the function for each segment

    # Use the function in your code
    max_mem_usage, max_mem_usage_gb = monitor_memory(90, max_mem_usage, max_mem_usage_gb,tmp_dir)  # stop the program if memory usage exceeds 90%

    print(f"第一部分Max memory usage during the program run was: {max_mem_usage}% ({max_mem_usage_gb} GB)")

    
    


    #! end of the function
    
    max_mem_usage, max_mem_usage_gb = monitor_memory(90, max_mem_usage, max_mem_usage_gb,tmp_dir)  # stop the program if memory usage exceeds 90%

    
    #! 计算每个进程的复杂度
    #complexities = [del_SV + ins_SV for del_SV, ins_SV in zip(del_SV_per_segment, ins_SV_per_segment)]
    # 计算每个进程的复杂度
    complexities = [del_SV + ins_SV + length for del_SV, ins_SV, length in zip(del_SV_per_segment, ins_SV_per_segment, unblock_vec_lengths)]
    
    # 创建一个包含进程索引和复杂度的列表
    indexed_complexities = list(enumerate(complexities))

    # 根据复杂度对进程进行排序
    sorted_complexities = sorted(indexed_complexities, key=lambda x: x[1], reverse=True)
    #number of kernel
    pool = Pool(args.cores)
    #xun_list=list(range(number_seg))
    # 创建一个包含进程索引的列表，按照复杂度从大到小排序
    xun_list = [index for index, complexity in sorted_complexities]
    print('运行顺序：'+str(xun_list))
    
    #xun_list = list(range(number_seg))
    args.cores = 5  # set the number of processes equal to the number of cores

    pool = Pool(args.cores)
    results = []

    # for xun in xun_list:
    #     result = pool.apply_async(gen_consensus, args=(xun, process_dict), callback=lambda _: monitor_memory(90, max_mem_usage, max_mem_usage_gb), error_callback=error_handler)
    #     results.append((xun, result))
    for xun in xun_list:
        # 将所有需要的变量作为参数传递给gen_consensus函数
        result = pool.apply_async(gen_consensus, args=(xun, chr_id, ll_c, starts_seg, ends_seg, tem_seq_post, del_SV_per_segment, ins_SV_per_segment, del_snv_per_segment, ins_snv_per_segment, snp_per_segment, pai_pro_tem_, args, condition_dist_sv_del, len_SV_del, condition_dist_sv_ins, len_SV_ins, ins_selection,base_list,substitution_matrix,mis_selection,times, TR_set, p_del_TR, p_ins_TR,tmp_dir), callback=lambda _: monitor_memory(90, max_mem_usage, max_mem_usage_gb,tmp_dir), error_callback=error_handler)
        results.append((xun, result))

    pool.close()
    pool.join()


    print(f"Max memory usage during the program run was: {max_mem_usage}% ({max_mem_usage_gb} GB)")

    end_time2 = time.time()
    process_end_time2 = time.process_time()

    elapsed_time2 = end_time2 - start_time2
    formatted_time2 = str(timedelta(seconds=elapsed_time2))

    process_formatted_time2 = process_end_time2 - process_time2
    process_second_time = str(timedelta(seconds=process_formatted_time2))

    print(f"并行处理运行时间：{formatted_time2}")
    print(f"并行处理 CPU执行实际时间：{process_second_time}")

    process_time3 = time.process_time()
    start_time3 = time.time()
    
    import gc
    gc.collect()


    # 获取结果并按照 xun 的值进行排序
    results = sorted([(xun, result.get()) for xun, result in results], key=lambda x: x[0])

    # 提取结果
    SV_table_segs = [result[0] for xun, result in results]
    VCF_table_segs = [result[1] for xun, result in results]
    unblock_region_segs = [result[2] for xun, result in results]
    Ins_dic_sv_segs = [result[3] for xun, result in results]
    tem_seq_post_segs = [result[4] for xun, result in results]


    # 优化拼接数据框的代码
    SV_table_combined = pd.concat(SV_table_segs, ignore_index=True)
    SV_table_merged = pd.concat([SV_table, SV_table_combined], ignore_index=True)
    VCF_table_combined = pd.concat(VCF_table_segs, ignore_index=True)
    VCF_table_merged = pd.concat([VCF_table, VCF_table_combined], ignore_index=True)

    end_time3 = time.time()
    process_end_time3 = time.process_time()
    start_time4 = time.time()
    # Combine all elements in unblock_region_segs
    unblock_region = list(set().union(*unblock_region_segs))
    # 初始化一个空字典来存储合并的结果
    Ins_dic_sv_combined = Ins_dic_sv
    for seg in Ins_dic_sv_segs:
        Ins_dic_sv_combined.update(seg)
        
    #print('ins number'+str(len(Ins_dic_sv_combined)))
    tem_seq_post_update = []

    for seg in tem_seq_post_segs:
        tem_seq_post_update.extend(seg)

    print('拼好的序列长度'+str(len(tem_seq_post_update)))
    #!
    #! check the simulated count
    # Count the number of deletions in SV_table_merged
    deletion_count = SV_table_merged[SV_table_merged['SV_type'] == 'Deletion'].shape[0]
    insertion_count = SV_table_merged[SV_table_merged['SV_type'] == 'Insertion'].shape[0]
    print('del_count:'+str(deletion_count))
    print('ins_count:'+str(insertion_count))
    end_time4 = time.time()

    start_time5 = time.time()

    tem_seq_post_up = tem_seq_post_update.copy()
    for idx in sorted(Ins_dic_sv_combined, reverse=True):
        tem_seq_post_up.insert(idx+1, Ins_dic_sv_combined[idx])

    tem_seq_post_up_string = ''.join(tem_seq_post_up)

    tem_seq_post_up_string= tem_seq_post_up_string.replace('-','')

    updated_con.append(copy.deepcopy(tem_seq_post_up_string))

    # 删除第一列
    SV_table_merged = SV_table_merged.iloc[:, 1:]

    # 定义 SV 类型的排序
    sv_type_order = ['Translocation', 'Inversion', 'Duplication', 'Deletion', 'Insertion', 'Micro_Del', 'Micro_Ins', 'Substitution']

    # 将 'SV type' 转换为有序的分类变量
    SV_table_merged['SV_type'] = pd.Categorical(SV_table_merged['SV_type'], categories=sv_type_order, ordered=True)

    # 按照 'SV type' 和 'Original_start' 排序
    SV_table_merged.sort_values(by=['SV_type', 'Original_start'], inplace=True)

    # 重置索引，并丢弃原来的索引
    SV_table_merged.reset_index(drop=True, inplace=True)

    # 再次重置索引，将新的索引添加为一个列，然后将这个新的列的名字改为 'Index'
    SV_table_merged.reset_index(inplace=True)
    SV_table_merged.rename(columns={'index': 'Index'}, inplace=True)

    # 按照 'pos' 列排序
    VCF_table_merged.sort_values(by='POS', inplace=True)

    # 重置索引并将旧索引添加为 'Index' 列
    VCF_table_merged.reset_index(inplace=True, drop=True)

    # 更新 'ID' 列
    VCF_table_merged['ID'] = 'rs' + VCF_table_merged.index.astype(str)

    tem_ins_dic = Ins_dic_sv_combined
    
    def write_template_fasta_con(args, seqname, consensus_, con_id):
        # Prepare the new sequence
        sequences = [consensus_]
        new_sequences = []
        for sequence in sequences:
            record = SeqRecord(Seq(re.sub('[^GATCN-]', "", str(sequence).upper())), id=seqname, name=seqname, description="<custom description>")
            new_sequences.append(record)

        # Write the new sequence to a file
        with open(args.save + 'BV_' + str(args.rep) + "_con" + str(con_id) + "_chr"+str(seqname) +"_ref.fasta", "w") as output_handle:
            SeqIO.write(new_sequences, output_handle, "fasta")
    
    def write_vcf(args, df, seqname, start_base, end_base,con_id):
        # Get the current date
        current_date = datetime.now().strftime('%Y%m%d')
        # Add the additional columns to the DataFrame
        df['FORMAT'] = 'GT'
        df['SAMPLE_ID'] = '1/1'
        # Write the DataFrame to a VCF file
        with open(args.save +'BV_' + str(args.rep) +"_con" + str(con_id) +  "_" +'chr'+str(seqname) +"_INDEL.vcf", 'w') as f:
            f.write('##fileformat=VCFv4.2\n')
            f.write('##fileDate=' + current_date + '\n')
            f.write('##source=Wave_region.py\n')
            f.write('##reference=' + args.ref + ':' + str(start_base) + '-' + str(end_base) + '\n')
            f.write('##contig=<ID='+str(seqname)+',length=' + str(end_base - start_base + 1) + '>\n')
            f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
            f.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the variant">\n')
            f.write('##INFO=<ID=CSV_TYPE,Number=1,Type=String,Description="Type of CSV">\n')
            f.write('##INFO=<ID=CSV_INDEX,Number=1,Type=Integer,Description="Index of CSV">\n')
            f.write('##INFO=<ID=SUB,Number=1,Type=String,Description="Substitution">\n')
            f.write('##INFO=<ID=microINS,Number=1,Type=String,Description="Micro insertion">\n')
            f.write('##INFO=<ID=LEN,Number=1,Type=Integer,Description="Length of the variant">\n')
            f.write('##INFO=<ID=microDEL,Number=1,Type=String,Description="Micro deletion">\n')
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_ID\n')
            df.to_csv(f, sep='\t', index=False, header=False)

    write_template_fasta_con(args, seqname, updated_con[0],ll_c)
    write_vcf(args, VCF_table_merged, seqname, start_base, end_base,ll_c)
    
    #! final table
    if args.write:
        print('finalize table')
        #!
        SV_table_merged = SV_write_relative(SV_table_merged,ll_c,tem_ins_dic)

     # Call the functions
    
        #SV_table_merged.to_csv(args.save + str(args.rep) + '_SVtable_full.txt', index=0, sep='\t')
        SV_table_merged.to_csv(args.save +'BV_' + str(args.rep) +'_con'+str(ll_c)+'_chr'+str(seqname)+ '_SVtable_full.csv', header=True, index=False)
    else:
        SV_table_merged.to_csv(args.save +'BV_' + str(args.rep) +'_con'+str(ll_c)+'_chr'+str(seqname)+ '_SVtable.csv', header=True, index=False)
        # Save the dictionary as a .npy file
        np.save(args.save+'BV_'+str(args.rep)+'_con'+str(ll_c)+'_chr'+str(seqname)+'_tem_ins_dic.npy', tem_ins_dic)
    
    
    end_time5 = time.time()
    # Use the function in your code
    max_mem_usage, max_mem_usage_gb = monitor_memory(90, max_mem_usage, max_mem_usage_gb,tmp_dir)  # stop the program if memory usage exceeds 90%

    print(f"整个过程Max memory usage during the program run was: {max_mem_usage}% ({max_mem_usage_gb} GB)")
    # At the end of the program, remove the temporary directory
    shutil.rmtree(tmp_dir)

    #! 主程序结束后，停止磁盘监控进程
    monitor_process.terminate()

    print(f"Initialization等待时间：{formatted_time0}")
    print(f"Initialization执行实际时间：{process_0_time}")

    print(f"取交集 等待时间：{formatted_time1}")
    print(f"取交集 CPU执行实际时间：{process_first_time}")

    print(f"并行处理运行时间：{formatted_time2}")
    print(f"并行处理 CPU执行实际时间：{process_second_time}")

    elapsed_time3 = end_time3 - start_time3
    formatted_time3 = str(timedelta(seconds=elapsed_time3))

    process_formatted_time3 = process_end_time3 - process_time3
    process_third_time = str(timedelta(seconds=process_formatted_time3))

    print(f"整理结果等待时间：{formatted_time3}")
    # print(f"写出结果CPU实际运行时间：{process_third_time}")
    elapsed_time4 = end_time4 - start_time4
    formatted_time4 = str(timedelta(seconds=elapsed_time4))
    print(f"拼接结果时间：{formatted_time4}")

    elapsed_time5 = end_time5 - start_time5
    formatted_time5 = str(timedelta(seconds=elapsed_time5))
    print(f"写出结果时间：{formatted_time5}")

    total_time = end_time5 - start_time0
    formatted_time6 = str(timedelta(seconds=total_time))
    print(f"总时间：{formatted_time6}")

    # 记得在结束时关闭文件
    sys.stdout.close()
    
    
    
if __name__ == "__main__":
    main()



