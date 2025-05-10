import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# 样本名列表
file_list = ['CHM1','CHM13','HG00268','HG00514','HG00733','NA19240','HG02818','HG01352','HG02059','NA12878','HG04217','HG02106','HG00268','AK1','HX1']

# 染色体列表
chromosomes = ['chr'+str(i) for i in range(1, 23)]

# 循环处理每条染色体
for chr_num in chromosomes:
    plt.figure(figsize=(10, 6))

    # 循环处理每个样本
    for sample in file_list:
        # 读取BED文件
        df = pd.read_csv(f'/disk18T3/project18/data/test_data/TGS/hg38/cell_15samples/{sample}_tel_dist_centro.bed', sep='\t', header=None, names=['chr', 'pos', 'SV type','chr length','TR Group','centro start','centro end'])

        # 计算 'tel dist' 列的值
        df['tel dist'] = np.where(df['pos'] < df['centro start'], df['pos'],
                                  np.where(df['pos'] >= df['centro end'], df['chr length'] - df['pos'],
                                           'CENTRO'))

        # 创建一个新的数据框，仅包含“tel dist”不等于“CENTRO”的行
        filtered_df = df[df['tel dist'] != 'CENTRO']
        centro_df = df[df['tel dist'] == 'CENTRO']

        # 转换'tel dist'列为数值类型
        df['tel dist'] = pd.to_numeric(df['tel dist'], errors='coerce')

        # 排除Y染色体
        df = df[df['chr'] != 'chrY']

        # 选出当前染色体的列
        df_chr = df[df['chr'] == chr_num]

        # 分别处理'DEL'和'INS'
        for sv_type in ['DEL', 'INS']:
            df_sv = df_chr[df_chr['SV type'] == sv_type]

            max_dist = df_sv['tel dist'].max()

            # 计算'tel dist'间隔为500,000的Bins的SV综合
            df_sv['dist_bin'] = pd.cut(df_sv['tel dist'], bins=np.arange(0, max_dist, 5e5))
            df_sv_grouped = df_sv.groupby('dist_bin').size()

            number_seg = len(df_sv_grouped)/2

            # 按照'TR Group'进行分组
            groups = df_sv.groupby('TR Group')

            # 对每个组分别绘制一条线
            for name, group in groups:
                df_sv_grouped = group.groupby('dist_bin').size()
                x= np.arange(0, number_seg, 0.5)
                y = df_sv_grouped.values
                plt.plot(x, y, marker='o', label=f'{name} - {sample}')

    # 设置x轴的标注
    plt.xticks(np.arange(0, number_seg, 5), rotation=90)

    # 添加标题和标签
    plt.title(f'Number of {sv_type} SV ({chr_num})')
    plt.xlabel('Telomere Distance (Mbp)')
    plt.ylabel('Number of SVs')

    # 添加图例
    plt.legend()

    # 保存图形到指定目录
    plt.savefig(f'/disk18T3/project18/data/test_data/TGS/hg38/plots/{sv_type}_{chr_num}.png')

    # 显示图形
    plt.show()
