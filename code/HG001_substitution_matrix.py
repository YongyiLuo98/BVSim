import pysam
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

main_url_save='~/evol_simulator/'
main_url_read='~/hg001/'
#vcf_in = VariantFile(main_url_read+"0_mapping.vcf") # auto-detect input format
#######GRCh38
vcf_in = pysam.VariantFile(main_url_read+"HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz") # auto-detect input format

base_list=["A","T","C","G"]
#contig_id_list= list(range(1, 23)) + ['X', 'Y']
contig_id_list = ['chr' + str(i) for i in list(range(1, 23))]
#contig_id_list = ['chr' + str(i) for i in list(range(1, 23))]
#overall substitution matrix
substitution_whole=np.zeros((4,4))

# Create a 4x6 grid of subplots
fig, axes = plt.subplots(nrows=11, ncols=2, figsize=(12, 8))


for chrom_index, ax in enumerate(axes.flatten()):
    # Generate random data for the heatmap
    chrom_ID=contig_id_list[chrom_index]
    chrom_ID=str(chrom_ID)
    
    ref_base=[]
    alt_base=[]
    #extract substitution from empirical
    substitution_vcf=np.zeros((4, 4))
    for record_test in vcf_in.fetch(chrom_ID):
        if record_test.alts is not None:
            # Get the REF and ALT fields
            ref_content = record_test.ref
            alt_content = record_test.alts[0]
            #only select one spot 
            if len(ref_content) ==1 and len(alt_content) ==1:
                ref_base.append(record_test.ref)
                #seq = ''.join(['C', 'C', 'A', 'C', 'G', 'A'])
                alt_base.append(record_test.alts[0])
                
    ref_seq = ''.join(ref_base)
    alt_seq=''.join(alt_base)
    if len(ref_seq)==len(alt_seq):
        len_seq=len(ref_seq)
    else:
        len_seq=0
    print(len_seq)
    for index_vcf_test in range(len_seq):
        ref_base=ref_seq[index_vcf_test]
        alt_base=alt_seq[index_vcf_test]
        if ref_base in base_list and alt_base in base_list:
            ref_id=base_list.index(ref_base)
            alt_id=base_list.index(alt_base)
            substitution_vcf[ref_id,alt_id]=substitution_vcf[ref_id,alt_id]+1
        
            # add corresponding elements
            substitution_whole=np.add(substitution_whole, substitution_vcf)
    
    #plot a heatmap for each chromosome's matrix
    matrix_plot = substitution_vcf
    row_sums = matrix_plot.sum(axis=1)
    substitution_normalized = matrix_plot / row_sums[:, np.newaxis]
    sub_matrix_vcf_DBSNP=substitution_normalized.astype(float)
    matrix_sub_plot=sub_matrix_vcf_DBSNP
    print(matrix_sub_plot) 
    
    
    data=matrix_sub_plot
    #data = np.random.rand(10, 10)

    # Plot the heatmap
    im = ax.imshow(data, cmap='YlGnBu')
    ax.set_xticks(range(len(data)))
    ax.set_yticks(range(len(data)))
    ax.set_xticklabels(['A', 'T', 'C', 'G'])
    ax.set_yticklabels(['A', 'T', 'C', 'G'])

    # Set the title of the subplot
    #ax.set_title(f'Chrom {chrom_index+1}')
    ax.set_title(f'Chrom {chrom_ID}')
    # Remove the x and y ticks
    # ax.set_xticks([])
    # ax.set_yticks([])

    # Add a colorbar to the subplot
    cbar = ax.figure.colorbar(im, ax=ax)

# Adjust the spacing between the subplots
plt.tight_layout()
plt.savefig(main_url_save+'heatmap_24chroms.jpg')
# Show the figure
plt.show()

# save the overall matrix to a CSV file
#normalize
matrix = substitution_whole
row_sums = matrix.sum(axis=1)
substitution_normalized = matrix / row_sums[:, np.newaxis]
substitution_float=substitution_normalized.astype(float)
np.savetxt(main_url_save+'substitution_HG001.csv', substitution_float, delimiter=',', header='', comments='')
plt.figure()
#matrix = np.random.rand(4, 4)
sns.heatmap(substitution_float, annot=True, cmap='YlGnBu', fmt='.2f')
# Set column and row names
plt.xticks([0.5, 1.5, 2.5, 3.5], ['A', 'T', 'C', 'G'])
plt.yticks([0.5, 1.5, 2.5, 3.5], ['A', 'T', 'C', 'G'])
# Show all numbers
plt.rcParams['font.size'] = 8
plt.savefig(main_url_save+'heatmap_sub_HG001.jpg')
plt.show()
