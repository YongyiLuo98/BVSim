# BVSim: A Benchmarking Variation Simulator Mimicking Human Variation Spectrum

[![Profile views](https://komarev.com/ghpvc/?username=YongyiLuo98&repo=BVSim&label=Profile%20views&color=0e75b6&style=flat)](https://github.com/YongyiLuo98/BVSim)
## Table of Contents

- [Getting Started](#getting-started)
- [Installation](#installation)
- [Configuration Files and CodesInstallation](#configuration-files-codes)
- [General Functions and Parameters](#parameters)
  - [Shared Parameters for Simulation Modes](#shared-parameters)
    - [Output Naming Conventions](#output)
    - [Write the Relative Positions of Simulated Variations](#write)
    - [User-defined Block Regions with No Variations](#block)
  - [Mimic Mode](#human-genome)
    - [Parameters for Mimic Mode](#parameters-for-mimic-mode)
  - [Wave Mode](#wave-mode)
    - [User-defined Sample(s) and Input BED File Requirements](#requirements-for-the-bed-file)
    - [Generate a BED File for a Single Sample](#generating-a-bed-file-for-a-single-sample-in-wave-mode)
    - [Job Submission for Single Sample (BED Format)](#job-submission-for-wave-mode-single-sample)
    - [Generating BED Files for Multiple Samples](#generating-bed-files-for-multiple-samples-in-wave-mode)
    - [Job Submission for Multiple Samples (BED Format)](#job-submission-for-wave-mode-multiple-samples)
    - [Important Note on File Placement](#important-note-on-file-placement)
    - [Parameters for Wave Mode](#parameters-for-wave-mode)
  - [Wave Region Mode](#wave-region-mode)
    - [Extract User-defined Regions (e.g. TR region) and Generate the BED File](#step-1-extract-tr-regions)
    - [Job Submission for Single Sample (BED Format)](#job-submission-for-wave-region-mode-single-sample)
    - [Parameters for Wave Region Mode](#parameters-for-wave-region-mode)
  - [Complex SV Mode](#complex-sv-mode)
    - [Parameters for CSV Mode](#parameters-for-csv-mode)
  - [Uniform Mode](#uniform-mode)
  - [Uniform Parallel Mode](#uniform-parallel-mode)
    - [Parameters for Uniform parallel Mode](#parameters-for-uniform-parallel-mode)
  - [Exact Mode](#exact-mode)
   - [Parameters for Exact Mode](#parameters-for-exact-mode)
- [VCF Mode for Preprocessing](#vcf-mode)
  - [Parameters for VCF Preprocessing](#parameters-for-VCF-mode)
- [Uninstallation for Updates](#uninstallation)
- [Workflow of BVSim](#workflow)
- [Definitions of SVs Simulated by BVSim](#definitions)

## <a name="getting-started"></a>Getting Started

To get started with BVSim, follow these steps to install and run the simulator:

```sh
# Create an envrionment called BVSim and install the dependencies from the provided environment.yml file
conda env create -f environment.yml
conda activate BVSim
# Installzation
## Clone the repository in your home path
git clone https://github.com/YongyiLuo98/BVSim.git
## Navigate to the ~/BVSim/ directory and install the package
cd ~/BVSim/
conda activate BVSim
pip install -e .

# Verify the installation
conda activate BVSim
which bvsim
bvsim --help
bvsim -h

## Run a toy example with a the modifiable YAML parameters in the cloned folder
conda activate BVSim
bvsim -config ~/BVSim/code/custom_config.yaml -sv_dup 20 \

## If you prefer using the default reference, simply execute (run with default configuration file: ~/BVSim/code/bvsim_config.yaml)
bvsim

# Generate variations with specific parameters (input parameters will overwrite the default settings)
bvsim -seed 1 -rep 1 -snp 2000

# To write out the relative positions, use the following command
python ~/BVSim/main/write_SV.py ~/BVSim/save/ BV_1_con0_chr1_SVtable.csv BV_1_con0_chr1_tem_ins_dic.npy

# Create a block intervals BED file
cd your_home_path
echo -e "0\t1000\n3000\t4000" > block_intervals.bed

# Run the simulator with block regions
bvsim -seed 1 -rep 1 -write -snp 2000 -block_region_bed_url block_intervals.bed
```

## <a name="Installation"></a>Installation
### Create an envrionment called BVSim and install the dependencies
To start with, you need to install the dependent packages in an environment, for example called BVSim.
```bash
# Create an envrionment called BVSim and install the dependencies from the provided environment.yml file
conda env create -f environment.yml
conda activate BVSim
```
### Clone the Repository
Next, you need to clone the BVSim repository to your local machine. Execute the following command in your home directory:
```bash
cd ~
git clone https://github.com/YongyiLuo98/BVSim.git
```
### Navigate to the Main Directory and Install the Package
Next, navigate to the .../BVSim/main/ directory to install the package:
```bash
cd ~/BVSim/
conda activate BVSim
pip install -e .
conda deactivate
```
### Verify the Installation
After installation, you can verify it from your BVSim conda environment. Execute the following commands:
```bash
conda activate BVSim
which bvsim
bvsim --help
bvsim -h
```
Note: You can only call BVSim in the cloned repository directory, while the installation must take place in the BVSim/main/ directory.
#### Toy Example (Uniform mode):
```bash
conda activate BVSim
bvsim -ref 'your_home_path/BVSim/empirical/sub_hg19_chr1.fasta' -seed 0 -rep 0 -write -snp 2000
```
or you can use the default reference to test the installation by type the following in your home path. If you do not give a saving path, the outputs will go to "your_home_path\BVSim\save\".

```bash
bvsim 
```
## <a name="configuration-files-codes"></a>Configuration Files and Codes

### Code Structure Overview

#### Shell Scripts (Runtime Statistics & Workflows)
| File | Description | 
|------|-------------|
| `CSV.sh` | Generate CSV reports with runtime metrics (CSV mode)|
| `exact.sh` | Precision analysis pipeline (time/memory profiled) (exact mode)|
| `hg19_whole_time.sh` | Genome-wide benchmark script for hg19 (mimic mode) |
| `hg38_whole_time.sh` | Genome-wide benchmark script for hg38 (mimic mode)|
| `uniform_test.sh` | Single-node consistency tests (uniform mode)|
| `uniform_parallel_test.sh` | Multi-node parallel consistency tests (uniform-parallel mode)|
| `wave.sh` | Main waveform simulation pipeline (wave mode)|
| `wave_region.sh` | Specialized TR (Tandem Repeat) waveform analysis  (wave-region mode)|

#### TR Processing Utilities
| File | Description |
|------|-------------|
| `extract_TR_region_hg19.sh` | hg19 TR region extraction |
| `extract_TR_region_hg38.sh` | hg38 TR region extraction |
| `extract_sv_to_bed.sh` | Convert SV calls to BED format for visualization |
| `terminal_commands` | Linux command files for data processing |
#### Python Modules
| File | Description |
|------|-------------|
| `HG001_substitution_matrix.py` | Calculate nucleotide substitution matrices from HG001 data |

#### Configuration Files
| File | Description |
|------|-------------|
| `bvsim_config.yaml` | Default parameters stored in YAML configuration file |
| `custom_config.yaml` | User-configurable parameters stored in YAML configuration file |

## <a name="parameters"></a>Functions and Parameters

Seven sequence simulation modes: mimic, wave, wave_region, csv, uniform, uniform parallel, exact

### <a name="shared-parameters"></a>Shared Parameters for Simulation Modes
The BVSim package provides several functions (modes) and parameters for simulating genetic variations. Here is a table that introduces all the functions and different parameters:

| Parameter | Type | Description | Default | 
| --- | --- | --- | --- |
| `-config` | str | Path to YAML configuration file | '~/BVSim/code/bvsim_config.yaml' |
| `-ref` | str | Input reference file | '~/BVSim/empirical/sub_hg19_chr1.fasta' |
| `-seq_index` | int | Index of sequence to use (0-based), must be an integer within the range of provided FASTA file. Default: 0 (first sequence) | 0 |
| `-save` | str | Saving path | '~/BVSim/save/' |
| `-seed` | int | Global seed for random number generator (non-negative integer) | 999 |
| `-times` | int | Maximum sampling times (positive integer) | 10 |
| `-rep` | int | Replication ID (non-negative integer for naming the files) | 99 |
| `-sv_trans` | int | Number of trans SV (non-negative integer) | 5 |
| `-sv_inver` | int | Number of inversion SV (non-negative integer) | 5 |
| `-sv_dup` | int | Number of tandem duplication (non-negative integer) | 5 |
| `-sv_del` | int | Number of SV deletion (non-negative integer) | 5 |
| `-sv_ins` | int | Number of SV insertion (non-negative integer) | 5 |
| `-snp` | float | SNV number (non-negative integer) or probability (between 0 and 1) | 5 |
| `-snv_del` | float | SNV deletion number (non-negative integer) or probability (between 0 and 1) | 5 |
| `-snv_ins` | float | SNV insertion number (non-negative integer) or probability (between 0 and 1) | 5 |
| `-notblockN` | bool | Do not Block N positions | False |
| `-write` | bool | Write full results | False |
| `-delmin` | int | Minimum deletion length (integer, not smaller than 50) | 50 |
| `-delmax` | int | Maximum deletion length (integer, larger than delmin) | 60 |
| `-insmin` | int | Minimum insertion length (integer, not smaller than 50) | 50 |
| `-insmax` | int | Maximum insertion length (integer, larger than insmin) | 450 |
| `-dupmin` | int | Minimum duplication length (integer, not smaller than 50) | 50 |
| `-dupmax` | int | Maximum duplication length (integer, larger than dupmin) | 450 |
| `-invmin` | int | Minimum inversion length (integer, not smaller than 50) | 50 |
| `-invmax` | int | Maximum inversion length (integer, larger than invmin) | 450 |
| `-transmin` | int | Minimum translocation lengthMinimum translocation length (integer, not smaller than 50) | 50 |
| `-transmax` | int | Maximum translocation lengthMaximum translocation length (integer, larger than transmin) | 450 |
| `-block_region_bed_url` | str | local path of the block region BED file | None |

#### <a name="output"></a>Output Naming Conventions
When you run the simulation tool, the output files are named based on the sequence name you input or the parameter `rep` you set (repetition number). Below is a summary of the output files you can expect:

1. **FASTA File**:  
   The output FASTA file will be named as follows:
```
BV_<rep>_seq_<seqname>.fasta
```
This file contains the simulated sequence.

2. **VCF File**:  
The VCF file will be named:
```
BV_<rep>_seq_<seqname>.vcf
```
This file stores the simulated variations.

3. **SV Table**:  
The SV table will have different naming conventions depending on whether you choose to include relative positions:
- If you include relative positions (by using the `-write` flag):
  ```
  BV_<rep>_seq_<seqname>_SVtable_full.csv
  ```
- If you do not include relative positions:
  ```
  BV_<rep>_seq_<seqname>_SVtable.csv
  ```

4. **Numpy File**:  
The numpy file that records all inserted segments we need to update the relative positions will be named:
```
BV_<rep>_seq_<seqname>_tem_ins_dic.npy
```
#### <a name="write"></a>Write the Relative Positions of Simulated Variations
If you choose not to generate the relative positions during the initial simulation run (i.e., you do not include the `-write` flag), the columns for relative positions in the SV table will be empty. However, you can still update these relative positions later using the saved intermediate files.
##### Steps to Write Relative Positions After Simulation
1. **Run the Initial Simulation**:  
For example, you can execute:
```bash
bvsim -seed 1 -rep 1 -snp 2000
```
In this case you generated default number of elementary SVs and micro indels, as well as 20000 SNPs saved in the default directory with `BV_1_seq_chr1_SVtable.csv`, `BV_1_seq_chr1_tem_ins_dic.npy`.

2. **Update Relative Positions**:
You can then run the following command to generate a table with the relative positions:
```bash
python ~/BVSim/main/write_SV.py ~/BVSim/save/ BV_1_seq_chr1_SVtable.csv BV_1_seq_chr1_tem_ins_dic.npy
```
This command will create a file called called `full_BV_1_seq_chr1_SVtable.csv` in the same directory, which will contain the relative positions for all variations with respect to the consensus sequence.
By following this naming convention and steps, you can easily manage and update your output files as needed.
#### <a name="block"></a>User-defined Block Regions with No Variations
The input of the '-block_region_bed_url' should be two columns of positions(start;end) without headers seperated by '\t'. To create a bed file, you can refer to the following example. In this case, positions from 0 to 999, from 3000 to 3999 cannot have any variation, so called blocked.

#### Toy Example:
```bash
cd ~
echo -e "0\t1000\n3000\t4000" > block_intervals.bed
# uniform.py
bvsim -seed 1 -rep 1 -write -snp 2000 -block_region_bed_url block_intervals.bed
```

### <a name="human-genome"></a>Mimic Mode
For the human genome, we derive the length distributions of SVs from HG002 and the 15 representative samples. For SNPs, we embed a learned substitution transition matrix from the dbSNP database. With a user-specified bin size, BVSim learns the distribution of SV positions per interval. It can model the SVs per interval as a multinomial distribution parameterized by the observed frequencies in HG002 (GRCh37/hg19 as reference) or sample the SV numbers per interval from a Gaussian distribution with the mean and standard deviation computed across the 15 samples (GRCh38/hg38 as reference). Calling ‘-hg19’ or ‘-hg38’ and specifying the chromosome name can activate the above procedures automatically for the human genome.

In the following example, we use 5 cores and 500,000 as length of the intervals. The reference is chromosome 21 of hg19, so we call "-hg19 chr21" in the command line to utilize the default procedure. In addition, we generated 1,000 SNPs, 99 duplications, 7 inversions, 280 deletions, and 202 insertions. The ratio of deletions/insertions in the tandem repeat regions with respect to the total number is 0.810/0.828. We also set the minimum and maximum lengths of some SVs.

#### <a name="parameters-for-Mimic-mode"></a>Parameters for Mimic Mode
The table below summarizes the parameters available for Wave region mode:
| Parameter | Type | Description | Default |
| --- | --- | --- | --- |
| `-mimic` | T/F | Realistic genome simulation (requires -hg38/-hg19) | False |
| `-hg19` | Chromosome name (e.g. chr1-chr22) for hg19 genome |  None |
| `-hg38` | Chromosome name (e.g. chr1-chr22) for hg38 genome | None |
| `-cell` | T/F | Use ONLY the Cell dataset list (hg38). If neither -cell nor -hgsvc is specified, defaults to merged CELL+HGSVC (deduplicated). | False |
| `-hgsvc` | T/F | Use ONLY the HGSVC dataset list (hg38). If neither -cell nor -hgsvc is specified, defaults to merged CELL+HGSVC (deduplicated). | False |
| `-cores` | int | Number of kernels for parallel processing (positive integer, required for uniform-parallel/wave/wave-region mode to set up parallel computing) | 1 |
| `-len_bins` | int | Length of bins for parallel processing, must be positive integer and smaller than reference length (required for uniform-parallel/wave/wave-region mode to set up parallel computing) | 50000 |
| `-mode` | str | Mode for calculating probabilities (empirical/probability)| 'probability' |
| `-sum` | bool | Total indel SV equals sum of the input (single sample) or mean of the input (multiple samples) | False |
| `-indel_input_bed` | str | Input BED file for indels (required if input single sample for wave or wave-region mode) | None |
| `-file_list` | str | List of sample files (e.g. NA19240_chr21.bed, HG02818_chr21.bed, NA19434_chr21.bed in empirical folder) (required if multiple samples for wave or wave-region mode) | None |
| `-p_del_region` | float | Probability of SV DEL (between 0 and 1) in the user-defined region for deletion (required for wave-region mode)| 0.5 |
| `-p_ins_region` | float | Probability of SV INS (between 0 and 1) in the user-defined region for insertion (required for wave-region mode)| 0.5 |
| `-region_bed_url` | str | Path of the BED file for the user-defined region (required for wave-region mode)| None |

#### Toy example (-hg19)
```bash
#!/bin/bash
#SBATCH -J 0_hg19_chr21
#SBATCH -N 1 -c 5
#SBATCH --output=output.txt
#SBATCH --error=err.txt

source /opt/share/etc/miniconda3-py39.sh
conda activate BVSim
bvsim -mimic \
-ref your_home_path/hg19/hg19_chr21.fasta -save your_home_path/test_data/BVSim/ -seed 0 -rep 0 -cores 5 -len_bins 500000 -hg19 chr21 \
-mode probability -snp 1000 -sv_trans 0 -dup 99 -sv_inver 7 -sv_del 280 -sv_ins 202 -snv_del 0 -snv_ins 0 -p_del_region 0.810 -p_ins_region 0.828 \
-region_bed_url /home/project18/data/test_data/TGS/hg002/chr21_TR_unique.bed -delmin 50 -delmax 2964912 -insmin 50 -insmax 187524
conda deactivate
```

### <a name="wave-mode"></a>Wave Mode

In Wave mode, users can provide a `.bed` file generated from an empirical `.vcf` file (for example, from HG002) or multiple BED files derived from samples of a selected population (such as the 15 Cell samples). This functionality allows you to generate non-uniform insertions and deletions with various options.

#### <a name="requirements-for-the-bed-file"></a>User-defined Sample(s) and Input BED File Requirements

The BED file must adhere to the following requirements:

- **First Column**: Location (genomic position)
- **Second Column**: DEL/INS label (indicating if the variation is a deletion or insertion)
- **Third Column**: Length (absolute value of the variation)

Each column should be separated by a tab character (`\t`) and must not include headers. Additionally, each BED file should represent variations on the same sequence.

#### <a name="generating-a-bed-file-for-a-single-sample-in-wave-mode"></a>Generate a BED File for a Single Sample

To generate a single input BED file from the HG002 `.vcf` file of chromosome 21, you can use the following commands in your terminal:

```bash
# Download the VCF file and its index
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz 
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz.tbi

# Activate the bcftools environment
conda activate bcftools

# Generate the BED file using bcftools and awk
bcftools view -H -r 21 -i 'SVTYPE="INS" || SVTYPE="DEL"' /home/adduser/data/test_data/TGS/hg002/HG002_SVs_Tier1_v0.6.vcf.gz | \
awk -v OFS='\t' '{
    split($8, a, ";");
    for (i in a) {
        if (a[i] ~ /^SVTYPE/) {
            split(a[i], b, "=");
            svtype = b[2];
        }
        else if (a[i] ~ /^SVLEN/) {
            split(a[i], c, "=");
            svlen = c[2];
            if (svlen < 0) svlen = -svlen;  # Extract the absolute value of SV length
        }
    }
    print $2, svtype, svlen;  # Print the location, SV type, and absolute SV length
}' > /home/adduser/data/test_data/TGS/hg002/chr21_SV_Tier1.bed
```
##### <a name="job-submission-for-wave-mode-single-sample"></a>Job Submission for Single Sample (BED Format)

To utilize this single BED file, users should call '-indel_input_bed' in the command. Below is the example of a SLURM job script that you can use to run the Wave mode simulation with single empirical data:

```bash
#!/bin/bash
#SBATCH -J full_chr21_parallel
#SBATCH -N 1 -c 5
#SBATCH --output=output_chr21_wave.txt
#SBATCH --error=err_chr21_wave.txt

source /opt/share/etc/miniconda3-py39.sh
conda activate BVSim
bvsim -wave \
-ref ~/hg19_chr21.fasta -save your_home_path/test_data/BVSim -seed 0 -rep 2 -cores 5 -len_bins 500000 \
-indel_input_bed ~/BVSim/empirical/chr21_SV_Tier1_2.bed -mode empirical -snp 2000 -snv_del 1000 -snv_ins 100 -write
conda deactivate
```
Submit the job file by:
```bash
sbatch task02_single.job
```
#### <a name="generating-bed-files-for-multiple-samples-in-wave-mode"></a>Generating BED Files for Multiple Samples

In this section, we will outline the steps to generate `.bed` files for multiple cell samples from the original Excel spreadsheet, using the 15 Cell samples as an example.

##### Step 1: Download the Original Excel File

First, download the Excel file containing the cell samples data:

```python
import os

# Download the Excel file
os.system('wget https://ars.els-cdn.com/content/image/1-s2.0-S0092867418316337-mmc1.xlsx')

```
##### Step 2: Load and View the Data
Next, load the Excel file into a Pandas DataFrame and view the first few rows:
```python
import pandas as pd

# Read the Excel file into a DataFrame
file_path = '1-s2.0-S0092867418316337-mmc1.xlsx'
df = pd.read_excel(file_path, sheet_name=0)  # Choose the correct sheet based on the file

# Display the first 5 rows of the DataFrame
print(df.head(5))
```
##### Step 3: Filter the Data
Extract the required columns and rename the first column:
```python
# Extract the necessary columns
columns_to_keep = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'MERGE_SAMPLES']
cell_df = df[columns_to_keep]

# List of all sample strings
samples = ['CHM1', 'CHM13', 'HG00514', 'HG00733', 'NA19240', 'HG02818', 'NA19434', 'HG01352', 'HG02059', 'NA12878', 'HG04217', 'HG02106', 'HG00268', 'AK1', 'HX1']
# selected population: the African population
AFR_samples = ['NA19240', 'HG02818', 'NA19434']

# Specify the columns to save in the BED file
columns_to_save = ['POS', 'SVTYPE', 'SVLEN']

# Extract rows where CHROM equals 'chr21'
chr21_df = cell_df[cell_df['CHROM'] == 'chr21']

# Display the first 10 rows for verification
print(chr21_df.head(10))

# Generate BED files for each sample in the AFR_samples list
for sample in AFR_samples:
    # Create a new DataFrame containing only rows where 'MERGE_SAMPLES' contains the current sample
    sample_df = chr21_df[chr21_df['MERGE_SAMPLES'].str.contains(sample)]

    # Specify the path for the new BED file
    bed_file_path = f'.../BVSim/empirical/{sample}_chr21.bed'

    # Save the specified columns to a BED file
    sample_df[columns_to_save].to_csv(bed_file_path, sep='\t', header=False, index=False)

```
#### <a name="job-submission-for-wave-mode-multiple-samples"></a>Job Submission for Multiple Samples (BED Format)

We provide an example of a Job submission script using SLURM for running the Wave mode with BVSim. This script utilizes the generated multiple sample BED files. Below is the example of a SLURM job script that you can use to run the Wave mode simulation with multiple samples:

```bash
#!/bin/bash
#SBATCH -J wave
#SBATCH -N 1 -c 5
#SBATCH --output=/home/project18/code/BVSim_code/wave2_out.txt
#SBATCH --error=/home/project18/code/BVSim_code/wave2_err.txt

source /opt/share/etc/miniconda3-py39.sh
conda activate BVSim

bvsim -ref your_home_path/hg38/chr21.fasta \
-save your_home_path/BVSim/task01/ -seed 0 -rep 1 -cores 5 \
-len_bins 500000 -wave -mode empirical -snp 2000 -snv_del 1000 -snv_ins 100 \
-write -file_list NA19240_chr21 HG02818_chr21 NA19434_chr21

conda deactivate
```
#### <a name="important-note-on-file-placement"></a>Important Note on File Placement
Ensure that both the single sample and multiple sample BED files are placed in the .../BVSim/empirical/ directory. This organization simplifies the command structure, allowing you to specify only the base names of the files (without extensions) directly in the -file_list option, as demonstrated in the script above.

#### <a name="parameters-for-wave-mode"></a>Parameters for Wave Mode

| Parameter | Type | Description | Default |
| --- | --- | --- | --- |
| `-wave` | T/F | Run Wave.py script | False |
| `-cores` | int | Number of kernels for parallel processing (positive integer, required for uniform-parallel/wave/wave-region mode to set up parallel computing) | 1 |
| `-len_bins` | int | Length of bins for parallel processing, must be positive integer and smaller than reference length (required for uniform-parallel/wave/wave-region mode to set up parallel computing) | 50000 |
| `-mode` | str | Mode for calculating probabilities (empirical/probability)| 'probability' |
| `-sum` | bool | Total indel SV equals sum of the input (single sample) or mean of the input (multiple samples) | False |
| `-indel_input_bed` | str | Input BED file for indels (required if input single sample for wave or wave-region mode) | None |
| `-file_list` | str | List of sample files (e.g. NA19240_chr21.bed, HG02818_chr21.bed, NA19434_chr21.bed in empirical folder) (required if multiple samples for wave or wave-region mode) | None |

##### Mode and Sum Parameters

The `-mode` parameter determines how the simulation calculates probabilities for insertions and deletions. It accepts two values:

- **'probability'**: In this mode, probabilities for insertions and deletions are derived from the empirical data provided in the input BED files. The total number of variations can be defined by the `-sum` parameter. If `-sum` is set to `True`, the total number of insertions or deletions will be the maximum of the calculated empirical total or the specified values in `-sv_ins` or `-sv_del`. This allows for flexibility in controlling the total number of SVs in the simulation.

- **'empirical'**: When set to this mode, the simulation directly uses the empirical values from the input data without any probability calculations. The total number of variations will be the sum of the provided empirical data.

The `-sum` parameter, when enabled, alters the total number of insertions and deletions based on the specified empirical data. If disabled, the simulation uses the fixed total values defined in `-sv_ins` and `-sv_del`, regardless of the empirical input.


### <a name="wave-region-mode"></a>Wave Region Mode

In Wave region mode, you can specify different INDEL probabilities using a BED file defined by `region_bed_url`. For example, if you want to increase the insertion and deletion probabilities in the tandem repeat (TR) regions of hg19, you can follow these steps.

#### <a name="step-1-extract-tr-regions"></a>Extract User-defined Regions (e.g. TR region) and Generate the BED File

First, extract the TR regions' positions from UCSC and create a BED file with two columns (start; end) separated by a tab character (`\t`).

You can generate the BED file using the following commands:

```bash
# Download the TR regions data
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/resources/hg19.simpleRepeat.bed.gz 
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/resources/hg19.simpleRepeat.bed.gz.tbi

# Extract the relevant columns and create the BED file
zcat hg19.simpleRepeat.bed.gz | awk 'BEGIN{OFS="\t"} {print $1, $2, $3}' > your_home_path/hg002/windows_TR.bed

# Merge overlapping intervals and remove duplicates
bedtools sort -i your_home_path/hg002/windows_TR.bed | bedtools merge -i stdin | awk 'BEGIN{OFS="\t"} {$4="TR"; print}' | uniq > your_home_path/hg002/windows_TR_unique.bed

# Filter for chromosome 21
awk '$1 == "chr21"' your_home_path/hg002/windows_TR_unique.bed > your_home_path/hg002/windows_TR_unique_chr21.bed

# Create a final BED file with start and end positions
awk '{print $2 "\t" $3}' your_home_path/hg002/windows_TR_unique_chr21.bed > your_home_path/hg002/chr21_TR_unique.bed
```
#### <a name="job-submission-for-wave-region-mode-single-sample"></a>Job Submission for Single Sample (BED Format)
In this example, we set the seed to `0` and use a replication ID of `4`. The job is configured to utilize `5` cores for parallel processing, with a bin size of `500,000`. We will generate `10,000` SNPs, along with `100` micro deletions and `100` micro insertions. The probabilities for these insertions and deletions are specified in the input BED file (`-indel_input_bed`) using the empirical mode (`-mode`). Additionally, we have set the probabilities for insertions (`-p_ins_region`) and deletions (`-p_del_region`) to approximately `0.6` for the total located in the TR region defined by `-region_bed_url`.

```bash
#!/bin/bash
#SBATCH -J full_chr21_parallel
#SBATCH -N 1 -c 5
#SBATCH --output=output_chr21_wave_region.txt
#SBATCH --error=err_chr21_wave_region.txt

source /opt/share/etc/miniconda3-py39.sh
conda activate BVSim
bvsim -ref your_home_path/hg19/hg19_chr21.fasta -save your_home_path/test_data/BVSim -seed 0 -rep 4 -cores 5 -len_bins 500000 -wave_region -indel_input_bed your_home_path/hg002/chr21_SV_Tier1.bed -mode empirical -snp 10000 -snv_del 100 -snv_ins 100 -write -p_del_region 0.6 -p_ins_region 0.6 -region_bed_url your_home_path/hg002/chr21_TR_unique.bed
conda deactivate
```

Submit the job file using the following command:
```bash
sbatch task03.job
```
#### <a name="parameters-for-wave-region-mode"></a>Parameters for Wave Region Mode
The table below summarizes the parameters available for Wave region mode:
| Parameter | Type | Description | Default |
| --- | --- | --- | --- |
| `-wave_region` | T/F | Run Wave_TR.py script | False |
| `-cores` | int | Number of kernels for parallel processing (positive integer, required for uniform-parallel/wave/wave-region mode to set up parallel computing) | 1 |
| `-len_bins` | int | Length of bins for parallel processing, must be positive integer and smaller than reference length (required for uniform-parallel/wave/wave-region mode to set up parallel computing) | 50000 |
| `-mode` | str | Mode for calculating probabilities (empirical/probability)| 'probability' |
| `-sum` | bool | Total indel SV equals sum of the input (single sample) or mean of the input (multiple samples) | False |
| `-indel_input_bed` | str | Input BED file for indels (required if input single sample for wave or wave-region mode) | None |
| `-file_list` | str | List of sample files (e.g. NA19240_chr21.bed, HG02818_chr21.bed, NA19434_chr21.bed in empirical folder) (required if multiple samples for wave or wave-region mode) | None |
| `-p_del_region` | float | Probability of SV DEL (between 0 and 1) in the user-defined region for deletion (required for wave-region mode)| 0.5 |
| `-p_ins_region` | float | Probability of SV INS (between 0 and 1) in the user-defined region for insertion (required for wave-region mode)| 0.5 |
| `-region_bed_url` | str | Path of the BED file for the user-defined region (required for wave-region mode)| None |


### <a name="exact-mode"></a>Exact Mode
Users can provide a variant table with target SVs with their length, positions and type. BVSim will simulate the variants with an non-overlapping feature. If some defined variations are overlapped, BVSim will sort them by start positions and discard the latter ones.
#### <a name="parameters-for-Exact-mode"></a>Parameters for Exact Mode
| Parameter | Type | Description | Default |
| --- | --- | --- | --- |
| `-exact` | T/F | Generate exact variants from input table | False |
| `-variant_table` | str | Path to variant table CSV file (required for exact mode) | None |
| `-validate_only` | T/F | Only validate input without generating variants| False |


### <a name="complex-sv-mode"></a>Complex SV Mode
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
#### Toy Example (CSV mode):
```bash
bvsim -ref '~/BVSim/empirical/sub_hg19_chr1.fasta' -save your_saving_url -seed 1 -rep 1 -csv -write -snp 2000
```
#### <a name="parameters-for-csv-mode"></a>Parameters for CSV Mode
The lengths of the CSVs follow different Gaussian distributions with modifiable means (-mu) and standard deviations (-sigma).
| Parameter | Type | Description | Default |
| --- | --- | --- | --- |
| `-csv` | T/F | Run csv.py script | False |
| `-csv_num` | int | Number for each type of CSV (non-negative integer), superior to -csv_total_num | 0 |
| `-csv_total_num` | int | Total number for CSV (non-negative integer), assign number of each type by empirical weights | 0 |
| `-num_ID1_csv to -num_ID18_csv` | int | Number of respective CSV types (non-negative integer)| 5 |
| `-mu_ID1 to -mu_ID18` | int | Mean of Gaussian distribution of CSV length (integer, larger than 100)| 1000 |
| `-sigma_ID1 to -sigma_ID18` | int | Standard deviation of Gaussian distribution of CSV length (non-negative integer)| 100 |

### <a name="uniform-mode"></a>Uniform Mode
If you call "-uniform", the simulation will be generated one by one uniformly.
| Parameter | Type | Description | Default |
| --- | --- | --- | --- |
| `-uniform` | T/F | Uniform distribution mode | False |
#### Toy Example (Uniform mode):
```bash
conda activate BVSim
bvsim -ref 'hg19_chr1.fasta' -uniform -seed 0 -rep 0 -write -snp 2000
```

### <a name="uniform-parallel-mode"></a>Uniform Parallel Mode
Add "-uniform" and specify "-cores" and "-len_bins" to your command, and write a .job file (task01.job) as follows (-c 5 means 5 cores, should be the same as -cores 5), parallel simulation will be allowed.

#### Toy Example (Uniform-parallel mode): task01.job
```bash
#!/bin/bash
#SBATCH -J uniform_parallel
#SBATCH -N 1 -c 5
#SBATCH --output=output.txt
#SBATCH --error=err.txt

source /opt/share/etc/miniconda3-py39.sh
conda activate BVSim
bvsim -uniform -ref ~/hg19_chr21.fasta -save ~/BVSim/task03/ \
-cores 5 -len_bins 500000 -rep 3 -snp 200 -snv_del 200 -snv_ins 200 -write
conda deactivate
```
Submit the job file by:
```bash
sbatch task01.job
```
#### <a name="parameters-for-uniform-parallel-mode"></a>Parameters for Uniform parallel Mode

| Parameter | Type | Description | Default |
| --- | --- | --- | --- |
| `-uniform` | T/F | Uniform distribution mode | False |
| `-cores` | int | Number of kernels for parallel processing (positive integer, required for uniform-parallel/wave/wave-region mode to set up parallel computing) | 1 |
| `-len_bins` | int | Length of bins for parallel processing, must be positive integer and smaller than reference length (required for uniform-parallel/wave/wave-region mode to set up parallel computing) | 50000 |


#### Toy example (exact mode)
```bash
#!/bin/bash
#SBATCH -J exact
#SBATCH -N 1 -c 1
#SBATCH --output=exact_test_out.txt
#SBATCH --error=exact_test_err.txt

source /opt/share/etc/miniconda3-py39.sh
conda activate BVSim
bvsim -exact \
-ref ~/BVSim/empirical/sub_hg19_chr1.fasta -seq_index 0 \
-seed 0 -rep 1 \
-variant_table ~/BVSim/empirical/BV_22_seq_1_SVtable_full.csv \
-write -notblockN
conda deactivate
```
We require the following format for the input '-variant_table'. A translocation is defined by two regions, A (Original_start-Original_end) and B (New_start-New_end), either balanced (exchanged, Balanced Trans Flag = 1) or unbalanced (A is lost and B is inserted to A’s start point, Balanced Trans Flag = 0) in the sequence. A duplication is defined by the copied region (Original_start-Original_end) and inserted position (New_start=New_end). An inversion is determined by one region (Original_start=New_start, Original_end = New_end). Each long or small insertion/deletion is defined by the start point and length. In the table, '-1' means unapplicable. We also ensure that the SVs and small variants will not be simulated from the regions related to other variants. Complex SVs combine these simple variants with spatial proximity.

This table's format is identical as the output of other modes (except for VCF mode). The relative positions (last four columns) are not required for 'exact mode'. So, users can utilize the randome generations' output and input some empirical SVs they want. Please pay attention, if your input has overlapped variants, the processing time will be longer. Please ensure your input does not contain any overlapped variants for smooth running.

See the complete file in the ~/BVSim/empirical/BV_22_seq_1_SVtable_full.csv.
| Index | Index_con | SV_type        | Original_start | Original_end | Len_SV | New_start | New_end | New_len_SV | Balanced Trans Flag | relative start1 | relative end1 | relative start2 | relative end2 |
|-------|-----------|----------------|----------------|--------------|--------|-----------|---------|------------|---------------------|-----------------|---------------|------------------|----------------|
| 0     | 0         | Translocation  | 12612          | 12804        | 193    | 70068     | 70139   | 72         | 1                   | 12611           | 12682         | 70179            | 70371          |
| ...   | ...       | ...            | ...            | ...          | ...    | ...       | ...     | ...        | ...                 | ...             | ...           | ...              | ...            |
| 3     | 0         | Translocation  | 128154         | 128326       | 173    | 96305     | 96305   | 173        | 0                   | 128971          | 128971        | 96989            | 97161          |
| 4     | 0         | Translocation  | 136841         | 137244       | 404    | 130197    | 130395  | 199        | 1                   | 137788          | 137986        | 130841           | 131244         |
| 5     | 0         | Inversion      | 59458          | 59511        | 54     | 59458     | 59511   | 54         | -1                  | 59619           | 59672         | -1               | -1             |
| 6     | 0         | Inversion      | 156931         | 157010       | 80     | 156931    | 157010  | 80         | -1                  | 157613          | 157692        | -1               | -1             |
| 7     | 0         | Duplication    | 38630          | 38704        | 75     | 76208     | 76208   | 75         | -1                  | 38641           | 38715         | 76441            | 76515          |
| 8     | 0         | Duplication    | 135854         | 135951       | 98     | 135951    | 135951  | 98         | -1                  | 136703          | 136800        | 136801           | 136898         |
| 9     | 0         | Deletion       | 27054          | 27105        | 52     | -1        | -1      | -1         | -1                  | -1              | -1            | -1               | -1             |
| 10    | 0         | Deletion       | 38903          | 38960        | 58     | -1        | -1      | -1         | -1                  | -1              | -1            | -1               | -1             |
| 11    | 0         | Deletion       | 68108          | 68157        | 50     | -1        | -1      | -1         | -1                  | -1              | -1            | -1               | -1             |
| 12    | 0         | Deletion       | 144007         | 144066       | 60     | -1        | -1      | -1         | -1                  | -1              | -1            | -1               | -1             |
| 13    | 0         | Deletion       | 166530         | 166583       | 54     | -1        | -1      | -1         | -1                  | -1              | -1            | -1               | -1             |
| 14    | 0         | Insertion      | 28137          | 28137        | 182    | -1        | -1      | -1         | -1                  | 27967           | 28148         | -1               | -1             |
| ...   | ...       | ...            | ...            | ...          | ...    | ...       | ...     | ...        | ...                 | ...             | ...           | ...              | ...            |
| 18    | 0         | Insertion      | 189262         | 189262       | 350    | -1        | -1      | -1         | -1                  | 189982          | 190331        | -1               | -1             |
| 19    | 0         | Small_Del      | 7033           | 7033         | 1      | -1        | -1      | -1         | -1                  | -1              | -1            | -1               | -1             |
| ...   | ...       | ...            | ...            | ...          | ...    | ...       | ...     | ...        | ...                 | ...             | ...           | ...              | ...            |
| 25    | 0         | Small_Ins      | 163605         | 163605       | 1      | -1        | -1      | -1         | -1                  | 164378          | 164378        | -1               | -1             |
| 26    | 0         | Substitution   | 8865           | 8865         | 1      | -1        | -1      | -1         | -1                  | 8864            | 8864          | -1               | -1             |



### <a name="vcf-mode"></a>VCF Mode for Pre-processing

#### <a name="parameters-for-VCF-mode"></a>Parameters for VCF Preprocessing
| Parameter | Type | Description | Default |
| --- | --- | --- | --- |
| `-vcf` | T/F | Run VCF.py script | False |
| `-vcf_file` | str | Input VCF file path (required for VCF mode) | None |
| `-chr` | str | Target chromosome name to filter from VCF file (e.g., chr21) (required for VC mode) | None |
| `-select` | str | Selection criteria (e.g., "AF>0.001", "SVLEN>=100") (required for VC mode) | None |
| `-min_len` | str | Minimum SV length (bp) (positice integer, required for VC mode) | 50 |
| `-sv_type` | str | SV types to include (required for VC mode) | ["DEL", "INS", "DUP", "INV"] |

#### Key Features
This module processes VCF/VCF.gz files to filter structural variants (SVs) and generate analysis-ready CSV files. Core capabilities:
- **Multi-format Support**: Handles both compressed (.vcf.gz) and uncompressed VCF
- **Flexible Filtering**:
  - Chromosome selection (`-chr`)
  - SV type filtering (`-sv_types`, default: DEL/INS/DUP/INV)
  - Length threshold (`-min_len`, default: 50bp)
  - Custom INFO field filters (`-select`, e.g., "AF>0.001")
- **Standardized Output**: Generates CSV with key fields which are the same as the Exact mode's input: 
  Index,Index_con,SV_type,Original_start,Original_end,Len_SV,New_start,...
  
#### Output Naming Convention
Files follow this structured naming pattern:
```
  <vcf_file>_<chr>_<select>_min<min_len>_<sv_type>.csv
```
#### Toy example (VCF mode)
```bash
#!/bin/bash
#SBATCH -J VCF
#SBATCH -N 1 -c 1
#SBATCH --output=VCF_out.txt
#SBATCH --error=VCF_err.txt

source /opt/share/etc/miniconda3-py39.sh
conda activate bio
bvsim -vcf \
-vcf_file gnomad.v4.1.sv.sites.vcf.gz \
-save ~/VCF/ \
-select "FREQ_HET_fin>0.001" -min_len 50 -sv_types DEL INS -chr chr21
conda deactivate
```
Input​​: gnomad.v4.1.sv.sites.vcf.gz with parameters: -chr 21 -select "FREQ_HET>0.001" -min_len 50 -sv_types DEL INS
​​Output​​: gnomad.v4.1_chr21_FREQ_HET_0.001_min50_DEL_INS.csv

Users need to pay attention that after filtering, there are overlapping SVs in the output. So, you may need to select your target SV and delete the overlapped ones.


## <a name="uninstallation"></a>Uninstallation for Updates
To update to the latest version of BVSim, you can uninstall and delete the cloned files. Then, try to clone from the new repository and install again.
```bash
conda activate BVSim
pip uninstall bvsim
```
## <a name="workflow"></a>Workflow of BVSim
The following figure illustrates the workflow of BVSim, encapsulated within a dashed box, and demonstrates how the output files interact with read simulators, the alignment tool Minimap2, Samtools, and evaluation tools.
![Workflow of BVSim](flow_chart_pipline.png)
## <a name="definitions"></a>Definitions of SVs Simulated by BVSim
![Definitions of SVs](Fig1a_BVSim.png)
