appendix_path=/research/groups/ma1grp/home/zyu/work_2025/RNA_direct_10_Oct/bai_project/appendix
basecalling_module=${appendix_path}/model_files/dna_r10.4.1_e8.2_400bps_sup@v5.0.0
environment_file=${appendix_path}/environment-ont.yaml

#=================================================================
#++++++++Step 0 download module file and create environment ++++++
#=================================================================
if [ ! -d ${basecalling_module} ]; then
    echo "Downloading dorado basecalling model file..."
    cd ${appendix_path}/model_files
    dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.0.0 
else
    echo "Dorado basecalling model file already exists."
fi 

conda env create -f ${environment_file}
conda activate BacDrop_Ont