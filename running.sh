
conda activate BacDrop_Ont
cd /research/groups/ma1grp/home/zyu/work_2025/my_github/Bulk-RNA-seq-for-Nanopore-long-reads

output=/research/groups/ma1grp/home/zyu/work_2025/my_github/running_test/bulk_long_test
pod5_files=/research/groups/ma1grp/home/zyu/work_2025/RNA_direct_10_Oct/bai_project/JB251030/2025-10-30_ESBL2/20251030_1241_P2S-02829-A_PBG67787_fbc5abd6/pod5_skip
appendix=/research/groups/ma1grp/home/zyu/work_2025/my_github/Bulk-RNA-seq-for-Nanopore-long-reads/appendix
sample_name=JB

bsub -P jb -J jb -n 2 -R "rusage[mem=4GB]" -eo jb.err -oo jb.out "
sh Bulk_RNA_seq_long_pipeline.sh -o $output -p $pod5_files -a $appendix -s $sample_name"