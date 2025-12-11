# Bulk-RNA-seq-for-Nanopore-long-reads
This pipeline is for specifical for Nanopre sequencing

- [Workflow](#Workflow)

# Workflow

# Software Installation

## 1. Install Dorado

This workflow uses **Dorado v1.0.0**. Download the Linux release from the link below:

**https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.0.0-linux-x64.tar.gz**

After downloading, extract the package and add the executable to your PATH if needed.

---

## 2. Install Other Required Software via Conda
The `appendix` directory contains `environment-ont.yaml`, mapping scripts, reference genome files, and other supplementary resources.  
You can clone this directory using `git clone` and then set `appendix_path` to the corresponding location on your own system.

```bash
git clone https://github.com/Zehui312/Bulk-RNA-seq-for-Nanopore-long-reads.git
appendix_path=/research/groups/ma1grp/home/zyu/work_2025/RNA_direct_10_Oct/bai_project/appendix
environment_file=${appendix_path}/environment-ont.yaml

conda env create -f ${environment_file}
conda activate BacDrop_Ont
```
## 3. Download the Basecalling Model
Download the basecalling model `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` and place it inside ${appendix_path}/model_files:
```
dorado download --model ${appendix_path}/model_files/dna_r10.4.1_e8.2_400bps_sup@v5.0.0
```