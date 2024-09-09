# BVSim
BVSim: A Benchmarking Variation Simulator Mimicking Human Variation Spectrum

# BVSim: A Benchmarking Variation Simulator Mimicking Human Variation Spectrum

## Table of Contents

- [Getting Started](#getting-started)

## <a name="getting-started"></a>Getting Started

To get started with BVSim, follow these steps to install and run the simulator:

```sh
# Installzation
## Clone the repository
cd your_home_path
git clone https://github.com/YongyiLuo98/BVSim.git
## Navigate to the main directory and install the package
cd your_home_path/BVSim/main/
pip install .
## To use BVSim, it is necessary to install dependencies by using conda or
conda activate (your_env)
pip install -r requirements.txt
## Alternatively, create an environment with the provided YAML file
conda env create -f environment.yml

# Verify the installation
cd your_home_path
python -m BVSim --help
python -m BVSim -h

## Run a toy example with a specified reference in the cloned folder
conda activate (your_env)
python -m BVSim -ref 'your_home_path/BVSim/empirical/sub_hg19_chr1.fasta' -seed 0 -rep 0 -write -snp 2000
## If you prefer using the default reference, simply execute
cd your_home_path
python -m BVSim


# Generate variations with specific parameters
cd your_home_path
python -m BVSim -seed 1 -rep 1 -snp 2000

# To write out the relative positions, use the following command
python your_home_path/BVSim/main/write_SV.py your_home_path/BVSim/save/ BV_1_con0_chr1_SVtable.csv BV_1_con0_chr1_tem_ins_dic.npy

# Create a block intervals BED file
cd your_home_path
echo -e "0\t1000\n3000\t4000" > block_intervals.bed

# Run the simulator with block regions
cd your_home_path
python -m BVSim -seed 1 -rep 1 -write -snp 2000 -block_region_bed_url block_intervals.bed
```

![Workflow of BVSim](flow_chart_pipline.png)

![Definitions of SVs](Fig1a_BVSim.png)

## Installation
```bash
cd your_home_path
git clone https://github.com/YongyiLuo98/BVSim.git
cd your_home_path/BVSim/main/
pip install .
```
To use the package smoothly, you can install the corresponding dependencies by either

```bash
conda activate (your_env)
pip install -r requirements.txt
```
or
```bash
conda env create -f environment.yml
```

### Verification
To verify if you can use BVSim, try to type:
```bash
cd your_home_path
python -m BVSim --help
python -m BVSim -h
```

#### Toy Example:
```bash
conda activate (your_env)
python -m BVSim -ref 'your_home_path/BVSim/empirical/sub_hg19_chr1.fasta' -seed 0 -rep 0 -write -snp 2000
```
or you can use the default reference to test the installation by type the following in your home path. If you do not give a saving path, the outputs will go to "your_home_path\BVSim\save\".

```bash
cd your_home_path
python -m BVSim 
```
## Functions and Parameters

Five modes: uniform, uniform parallel, csv, wave, wave_region

### Shared Parameters
The BVSim package provides several functions (modes) and parameters for simulating genetic variations. Here is a table that introduces all the functions and different parameters:

| Parameter | Type | Description | Default |
| --- | --- | --- | --- |
| `-ref` | str | Input reference URL | 'default_ref' |
| `-save` | str | Save URL | your_home_path/BVSim/save/ |
| `-seed` | int | Seed for random number generator | 999 |
| `-times` | int | Number of times | 10 |
| `-rep` | int | Replication ID | 5 |
| `-sv_trans` | int | Number of trans SV | 5 |
| `-sv_inver` | int | Number of inversion SV | 5 |
| `-sv_dup` | int | Number of tandem duplication | 5 |
| `-sv_del` | int | Number of SV deletion | 5 |
| `-sv_ins` | int | Number of SV insertion | 5 |
| `-snp` | float | SNV number or probability | 5 |
| `-snv_del` | float | SNV deletion number or probability | 5 |
| `-snv_ins` | float | SNV insertion number or probability | 5 |
| `-notblockN` | bool | Do not Block N positions | False |
| `-write` | bool | Write full results | False |
| `-delmin` | int | Minimum deletion length | 50 |
| `-delmax` | int | Maximum deletion length | 60 |
| `-insmin` | int | Minimum insertion length | 50 |
| `-insmax` | int | Maximum insertion length | 450 |
| `-dupmin` | int | Minimum duplication length | 50 |
| `-dupmax` | int | Maximum duplication length | 450 |
| `-invmin` | int | Minimum inversion length | 50 |
| `-invmax` | int | Maximum inversion length | 450 |
| `-block_region_bed_url` | str | local path of the block region BED file | None |

If '-write' is present, in the '....SV_table_full.csv', the relative positions of all variations with respect to the consensus will be in the columns containing 'relative'. If this is not necessary, you can drop this parameter in your command as it extends the total time if there are lots of variations (1 min/ 10000 variations). However, it is still possible to update the relative positions after the simulation. We will save the intermediate documents for this, see the example below.

#### Toy Example:
```bash
cd your_home_path
python -m BVSim -seed 1 -rep 1 -snp 2000
```
In this case you generated default number of elementary SVs and micro indels, as well as 20000 SNPs saved in the default directory. However, you did not write out the relative positions.
If you type the following you will get a file called: 'full_BV_1_con0_chr1_SVtable.csv' in the same directory.
```bash
python your_home_path/BVSim/main/write_SV.py your_home_path/BVSim/save/ BV_1_con0_chr1_SVtable.csv BV_1_con0_chr1_tem_ins_dic.npy
```
The input of the '-block_region_bed_url' should be two columns of positions(start;end) without headers seperated by '\t'. To create a bed file, you can refer to the following example. In this case, positions from 0 to 999, from 3000 to 3999 cannot have any variation, so called blocked.

#### Toy Example:
```bash
cd your_home_path
echo -e "0\t1000\n3000\t4000" > block_intervals.bed
# uniform.py
cd your_home_path
python -m BVSim -seed 1 -rep 1 -write -snp 2000 -block_region_bed_url block_intervals.bed
```

### Human genome
For human reference genome GRCh37 or GRCh38, users are recommended to call -hg19 or -hg38 in the command line for utilizing the HG002 and Cell dataset.

### Uniform mode
If you do not call any of the following parameters (-csv, -cores, -len_bins, -wave), the simulation will be generated one by one uniformly.

### Complex SV mode
Add -csv to your command, 18 types of Complex Structure Variations can be generated.

* ID1: Tandem Inverted Duplication (TanInvDup)
* ID2: Dispersed Inverted Duplication (DisInvDup)
* ID3: Dispersed Duplication (DisDup)
* ID4: Inversion with 5’ or 3’ Flanking Deletion (DEL+INV/INV+DEL)
* ID5: 5’ Deletion and Dispersed Inverted Duplication (DEL+DisInvDup)
* ID6: 5’ Deletion and Dispersed Duplication (DEL+DisDup)
* ID7: Tandem Duplication and 3’ Deletion (TanDup+DEL)
* ID8: Tandem Inverted Duplication and 3’ Deletion (TanInvDup+DEL)
* ID9: Tandem Duplication, Deletion and Inversion (TanDup+DEL+INV)
* ID10: Tandem Inverted Duplication, Deletion and Inversion (TanInvDup+DEL+INV)
* ID11: Paired-Deletion Inversion (DEL+INV+DEL)
* ID12: Inversion with 5’ Flanking Duplication (DUP+INV)
* ID13: Inversion with 3’ Flanking Duplication (INV+DUP)
* ID14: Paired-Duplication Inversion (DUP+INV+DUP)
* ID15: Inversion with 5’ Flanking Duplication and 3’ Flanking Deletion (DUP+INV+DEL)
* ID16: Inversion with 5’ Flanking Deletion and 3’ Flanking Duplication (DEL+INV+DUP)
* ID17: Inverted Duplication with Flanking Triplication (DupTripDup-INV)
* ID18: Insertion with Deletion (INSdel)
#### Toy Example:
```bash
cd your_home_path
python -m BVSim -ref 'your_home_path/BVSim/empirical/sub_hg19_chr1.fasta' -save your_saving_url -seed 1 -rep 1 -csv -write -snp 2000
```

The lengths of the CSVs follow different Gaussian distributions with modifiable means (-mu) and standard deviations (-sigma).
| Parameter | Type | Description | Default |
| --- | --- | --- | --- |
| `-csv_num` | int | Number for each type of CSV, superior to -csv_total_num | 0 |
| `-csv_total_num` | int | Total number for CSV, assign number of each type by empirical weights | 0 |
| `-num_ID1_csv to -num_ID18_csv` | int | Number of respective CSV types | 5 |
| `-mu_ID1 to -mu_ID18` | int | Mean of Gaussian distribution of CSV length | 1000 |
| `-sigma_ID1 to -sigma_ID18` | int | Standard deviation of Gaussian distribution of CSV length | 100 |

### Uniform parallel mode
Add -cores, -len_bins to your command, and write a .job file (task01.job) as follows (-c 5 means 5 cores, should be the same as -cores 5), parallel simulation will be allowed.

#### Toy Example: task01.job
```bash
#!/bin/bash
#SBATCH -J uniform_parallel
#SBATCH -N 1 -c 5
#SBATCH --output=output.txt
#SBATCH --error=err.txt

source /opt/share/etc/miniconda3-py39.sh
conda activate (your_env)
cd your_home_path
python -m BVSim -ref your_home_path/hg19/hg19_chr21.fasta -save your_home_path/test_data/BVSim/task03/ -cores 5 -len_bins 500000 -rep 3 -snp 200 -snv_del 200 -snv_ins 200 -write
conda deactivate
```
Submit the job file by:
```bash
sbatch task01.job
```

### Wave mode

If you provide an .bed file generated from an empirical .vcf file, for example from HG002, we can generate non-uniform SV INDELs with different options. The requirement of the BED file is that the first column is the breakpoint, the second is the DEL/INS label, seperated by '\t' without headers.
To generate the input BED file from a given .vcf file, try the following command in your terminal.
```bash
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz 
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz.tbi
conda activate bcftools
bcftools view -H -r 21 -i 'SVTYPE="INS" || SVTYPE="DEL"' your_home_path/hg002/HG002_SVs_Tier1_v0.6.vcf.gz | awk -v OFS='\t' '{split($8,a,";"); for (i in a) {if (a[i] ~ /^SVTYPE/) {split(a[i],b,"="); print $2, b[2]}}}' > your_home_path/hg002/chr21_SV_Tier1.bed
```
#### Toy Example: task02.job.
```bash
#!/bin/bash
#SBATCH -J full_chr21_parallel
#SBATCH -N 1 -c 5
#SBATCH --output=output_chr21_wave.txt
#SBATCH --error=err_chr21_wave.txt

source /opt/share/etc/miniconda3-py39.sh
conda activate your_env
cd your_home_path
python -m BVSim -ref your_home_path/hg19/hg19_chr21.fasta -save your_home_path/test_data/BVSim -seed 0 -rep 2 -cores 5 -len_bins 500000 -wave -indel_input_bed your_home_path/hg002/chr21_SV_Tier1.bed -mode empirical -snp 2000 -snv_del 1000 -snv_ins 100 -write
conda deactivate
```
Submit the job file by:
```bash
sbatch task02.job
```

| Parameter | Type | Description | Default |
| --- | --- | --- | --- |
| `-cores` | int | Number of kernels for parallel processing | 1 |
| `-len_bins` | int | Length of bins for parallel processing | 50000 |
| `-wave` | bool | Run Wave.py script | False |
| `-mode` | str | Mode for calculating probabilities | 'probability' |
| `-sum` | bool | Total indel SV equals sum of the input bed | False |
| `-indel_input_bed` | str | Input BED file for indels | 'your_home_path/hg002/chr21_SV_Tier1.bed' |

### Wave region mode
In this mode, we allow different INDEL probabilities in region_bed_url. For example, we want to set the SV INDELs probabilities to be higher in hg19's TR region. We can do the following things.
First, extract the TR regions' positions from the UCSC and make it two columns (start;end) in a BED file seperated by '\t'.

#### Toy Example: generation of the BED file
```bash
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/resources/hg19.simpleRepeat.bed.gz 
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/resources/hg19.simpleRepeat.bed.gz.tbi

zcat hg19.simpleRepeat.bed.gz | awk 'BEGIN{OFS="\t"} {print $1, $2, $3}' > your_home_path/hg002/windows_TR.bed
#merge overlapping intervals and remove duplicates
bedtools sort -i your_home_path/hg002/windows_TR.bed | bedtools merge -i stdin | awk 'BEGIN{OFS="\t"} {$4="TR"; print}' | uniq > your_home_path/hg002/windows_TR_unique.bed

##########your_home_path/hg002/windows_TR_unique.bed
awk '$1 == "chr21"' your_home_path/hg002/windows_TR_unique.bed > your_home_path/hg002/windows_TR_unique_chr21.bed

awk '{print $2 "\t" $3}' your_home_path/hg002/windows_TR_unique_chr21.bed > your_home_path/hg002/chr21_TR_unique.bed
```
Then, run the following job file.

#### Toy Example: task03.job.
```bash
#!/bin/bash
#SBATCH -J full_chr21_parallel
#SBATCH -N 1 -c 5
#SBATCH --output=output_chr21_wave_region.txt
#SBATCH --error=err_chr21_wave_region.txt

source /opt/share/etc/miniconda3-py39.sh
conda activate your_env
cd your_home_path
python -m BVSim -ref your_home_path/hg19/hg19_chr21.fasta -save your_home_path/test_data/BVSim -seed 0 -rep 4 -cores 5 -len_bins 500000 -wave_region -indel_input_bed your_home_path/hg002/chr21_SV_Tier1.bed -mode empirical -snp 10000 -snv_del 100 -snv_ins 100 -write -p_del_region 0.6 -p_ins_region 0.6 -region_bed_url your_home_path/hg002/chr21_TR_unique.bed
conda deactivate
```
Submit the job file by:
```bash
sbatch task03.job
```

| Parameter | Type | Description | Default |
| --- | --- | --- | --- |
| `-cores` | int | Number of kernels for parallel processing | 1 |
| `-len_bins` | int | Length of bins for parallel processing | 50000 |
| `-wave` | bool | Run Wave.py script | False |
| `-mode` | str | Mode for calculating probabilities | 'probability' |
| `-sum` | bool | Total indel SV equals sum of the input bed | False |
| `-indel_input_bed` | str | Input BED file for indels | 'your_home_path/hg002/chr21_SV_Tier1.bed' |
| `-wave_region` | bool | Run Wave_TR.py script | False |
| `-p_del_region` | float | Probability of SV DEL in the user-defined region for deletion | 0.5 |
| `-p_ins_region` | float | Probability of SV INS in the user-defined region for insertion | 0.5 |
| `-region_bed_url` | str | URL of the BED file for the user-defined region | 'your_home_path/hg002/chr21_TR_unique.bed' |

## Uninstallation

```bash
pip uninstall BVSim
```