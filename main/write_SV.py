import pandas as pd
import numpy as np
import sys
import time
from datetime import timedelta

start_time1 = time.time()
process_time1 = time.process_time()

main_url_SV_path = sys.argv[1]
main_url_SV_name = sys.argv[2]
main_url_save = main_url_SV_path
main_url_npy = sys.argv[3]
# 加载.npy文件
#tem_ins_dic = np.load(main_url_save+'PSDV_4_seed0_chr22_tem_ins_dic.npy', allow_pickle=True).item()
tem_ins_dic = np.load(main_url_SV_path+ main_url_npy, allow_pickle=True).item()
# 现在tem_ins_dic是一个字典，你可以像操作普通字典一样操作它

# 读取数据
SV_table_merged = pd.read_csv(main_url_SV_path+main_url_SV_name)
print(main_url_SV_path+main_url_SV_name)
ll_c = 0
#! final table
print('finalize table')
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

# for ll_var_index in range(len(whole_start_abs_set)):
#     #! time 找到对应的行
#     print(whole_start_abs_set[ll_var_index])
#     tem_SV_table_merged2 = tem_SV_table_merged[(tem_SV_table_merged['Original_start']==whole_start_abs_set[ll_var_index]) |\
#                                     (tem_SV_table_merged['New_start']==whole_start_abs_set[ll_var_index])]
    

#     for xun_nei_row in range(len(tem_SV_table_merged2)):
#         tem_row = tem_SV_table_merged2.iloc[xun_nei_row,:]
#         stand_line = int(tem_row[0])
#         #A
#         bone1s = tem_row[3]
#         bone1e = tem_row[4]
#         #B
#         bone2s = tem_row[6]
#         bone2e = tem_row[7]
#         if whole_start_abs_set[ll_var_index] in list_start1:
#         #ls_satrt1_index_df = int(list_start1.index(whole_start_abs_set[ll_var_index]))
#         #stand_line = int(SV_table_merged.iloc[ls_satrt1_index_df,0])
#         #tem_row = SV_table_merged.iloc[ls_satrt1_index_df,:]
#             #class of SV
#             if tem_row[2] in ['Substitution','Micro_Ins','Micro_Del','Deletion','Insertion','Inversion']:
#                 if tem_row[2] in ['Deletion','Micro_Del']:
#                     inster_number_bone = bone1s-last_bone-1
#                     #index for consensus before start of current variation
#                     present_len = present_len + inster_number_bone
#                     #update last_bone as end of current variation
#                     last_bone = bone1e
#                     #deleted base has no new axis on consensus
#                     SV_table_merged.iloc[stand_line,10] = present_len
#                     SV_table_merged.iloc[stand_line,11] = present_len
#                     SV_table_merged.iloc[stand_line,12] = -1
#                     SV_table_merged.iloc[stand_line,13] = -1
#                 elif tem_row[2] in ['Substitution']:
#                     inster_number_bone = bone1s-last_bone
#                     #one to one map
#                     present_len = present_len + inster_number_bone
#                     #bone1s=bone1e=5
#                     last_bone = bone1e
#                     SV_table_merged.iloc[stand_line,10] = present_len
#                     SV_table_merged.iloc[stand_line,11] = present_len
#                     SV_table_merged.iloc[stand_line,12] = -1
#                     SV_table_merged.iloc[stand_line,13] = -1
#                 elif tem_row[2] in ['Micro_Ins','Insertion']:
#                     inster_number_bone = bone1s-last_bone
#                     Ins_len_present = len(tem_ins_dic[bone1s])
#                     #inserted position on consensus: one pos:+1, inserted after current base
#                     SV_table_merged.iloc[stand_line,10] = present_len + inster_number_bone+1
#                     SV_table_merged.iloc[stand_line,12] = -1
#                     #on consensus: end of previous SV+ number of normal base+ inserted length
#                     present_len = present_len + inster_number_bone+Ins_len_present
#                     #end of current SV
#                     last_bone = bone1e
#                     SV_table_merged.iloc[stand_line,11] = present_len
#                     SV_table_merged.iloc[stand_line,13] = -1
#                 else:## this is the inversion
#                     inster_number_bone = bone1s-last_bone
#                     SV_table_merged.iloc[stand_line,10] = present_len + inster_number_bone
#                     SV_table_merged.iloc[stand_line,12] = -1
#                     #no loss from last_bone to bone1e
#                     #????
#                     #present_len = present_len + bone1e - last_bone
#                     present_len = present_len + bone1e - last_bone
#                     SV_table_merged.iloc[stand_line,11] = present_len
#                     SV_table_merged.iloc[stand_line,13] = -1 
#                     last_bone = bone1e
                    
#             elif tem_row[2] in ['Duplication']:
#                     #copy A to B (A no change)
#                     #5-0=5
#                     inster_number_bone = bone1s-last_bone
#                     #Ins_len_present = len(tem_ins_dic[bone2s])
#                     #length of the copied: A
#                     #=6
#                     tem_plate_len = SV_table_merged.iloc[stand_line,5]
#                     #0+5=5
#                     SV_table_merged.iloc[stand_line,10] = present_len + inster_number_bone
#                     #SV_table_merged.iloc[stand_line,12] = present_len + inster_number_bone+1
#                     present_len = present_len + inster_number_bone + tem_plate_len-1
#                     #0+5+6-1=10
#                     SV_table_merged.iloc[stand_line,11] = present_len 
#                     #SV_table_merged.iloc[stand_line,13] = present_len
#                     last_bone = bone1e
                    
#             elif tem_row[2] in ['Translocation']:
#                 #balanced translocation
#                 #A:5-10, B:12-18
#                 if tem_row[9] == 1:
#                     #ins B to A's pos:5-0-1=4
#                     inster_number_bone = bone1s-last_bone-1
#                     #0+4+1=5,the start of copied base is 5
#                     SV_table_merged.iloc[stand_line,10] = present_len + inster_number_bone+1
#                     #length of B: 18-12+1=7
#                     Ins_len_present = len(tem_ins_dic[bone1s-1])
#                     #0+4+7=11
#                     #end of A:current SV end=11
#                     present_len = present_len + inster_number_bone + Ins_len_present
#                     SV_table_merged.iloc[stand_line,11] = present_len
#                     last_bone = bone1e
#                 #!unbalanced trans:
#                 else:
#                     inster_number_bone = bone1s-last_bone-1
#                     #index for consensus before start of current variation
#                     present_len = present_len + inster_number_bone
                    
#                     #deleted base has no new axis on consensus
#                     SV_table_merged.iloc[stand_line,10] = present_len+1
#                     SV_table_merged.iloc[stand_line,11] = present_len+1
                    
#                     #update last_bone as end of current variation
#                     last_bone = bone1e
    
    
#         else:### in the list2: pos of B (only duplication and trans)
            
#             if tem_row[2] in ['Duplication']:
#                 #if SV_table_merged.iloc[stand_line,10]==0:
#                     #bone2s:B_start
#                     #same as ins
#                     inster_number_bone = bone2s-last_bone
#                     #SV_table_merged.iloc[stand_line,10] = present_len + inster_number_bone+1
#                     #SV_table_merged.iloc[stand_line,12] = present_len + inster_number_bone+1
#                     SV_table_merged.iloc[stand_line,12] = present_len+inster_number_bone+1
#                     Ins_len_present = len(tem_ins_dic[bone2s])
#                     present_len = present_len + inster_number_bone+Ins_len_present
                    
#                     SV_table_merged.iloc[stand_line,13] = present_len
#                     last_bone = bone2e
#             elif tem_row[2] in ['Translocation']:
#                 #balanced: similar to A
#                 if  tem_row[9] == 1:
#                     inster_number_bone = bone2s-last_bone-1
#                     SV_table_merged.iloc[stand_line,12] = present_len + inster_number_bone+1
#                     #inserted A's length
#                     Ins_len_present = len(tem_ins_dic[bone2s-1])
#                     present_len = present_len + inster_number_bone + Ins_len_present
#                     SV_table_merged.iloc[stand_line,13] = present_len
#                     last_bone = bone2e
#                 #unbalanced
#                 else:
#                     inster_number_bone = bone2s-last_bone-1
#                     inster_number_bone = bone2s-last_bone
#                     #A is a del
#                     SV_table_merged.iloc[stand_line,10] = -1
#                     SV_table_merged.iloc[stand_line,12] = present_len + inster_number_bone+1
#                     #length of A
#                     #Ins_len_present = len(tem_ins_dic[bone2s-1])
#                     #Ins_dic_sv_seg[ins_trans_loc] = copy.deepcopy(''.join(tem_seq_post[r_s:(r_s+l_s)]))
#                     #similar to insertion
#                     Ins_len_present = len(tem_ins_dic[bone2s])
#                     present_len = present_len + inster_number_bone + Ins_len_present
#                     #A is a del
#                     SV_table_merged.iloc[stand_line,11] = -1
#                     SV_table_merged.iloc[stand_line,13] = present_len
#                     last_bone = bone2e

# 将 Numpy 数组转回 DataFrame
#! time
# SV_table_merged_merged = pd.DataFrame(SV_table_merged_merged_np, columns=SV_table_merged.columns)

SV_table_merged.to_csv(main_url_save + 'full_'+ main_url_SV_name , header=True, index=False)

end_time1 = time.time()

elapsed_time1 = end_time1 - start_time1
formatted_time1 = str(timedelta(seconds=elapsed_time1))

print(f" relative pos writing out cost：{formatted_time1}")
