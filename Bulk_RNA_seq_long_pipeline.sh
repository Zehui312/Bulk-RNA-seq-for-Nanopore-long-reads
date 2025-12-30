#!/bin/bash


# Function to show usage
usage() {
    echo "Usage: $0 -o OUTPUT_PATH -p POD5_PATH -a APPENDIX_PATH -s SAMPLE_NAME"
    echo "  -o: Output path"
    echo "  -p: POD5 input path"
    echo "  -a: Appendix path"
    echo "  -s: Sample name"
    echo "  -h: Show this help message"
    exit 1
}

# Parse command line options
while getopts "o:p:a:s:h" opt; do
    case $opt in
        o)
            out_put_path="$OPTARG"
            ;;
        p)
            pod5_path="$OPTARG"
            ;;
        a)
            appendix_path="$OPTARG"
            ;;
        s)
            sample_name="$OPTARG"
            ;;
        h)
            usage
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done

# Check if all required parameters are provided
if [[ -z "$out_put_path" || -z "$pod5_path" || -z "$appendix_path" || -z "$sample_name" ]]; then
    echo "Error: All parameters are required."
    usage
fi



# out_put_path="/research/groups/ma1grp/home/zyu/work_2025/my_github/running_test/bulk_long_test"
# pod5_path="/research/groups/ma1grp/home/zyu/work_2025/RNA_direct_10_Oct/bai_project/JB251030/2025-10-30_ESBL2/20251030_1241_P2S-02829-A_PBG67787_fbc5abd6/pod5_skip"
# appendix_path=/research/groups/ma1grp/home/zyu/work_2025/my_github/Bulk-RNA-seq-for-Nanopore-long-reads/appendix
# sample_name=JB

## Define paths to tools and files
# Step 1 basecalling
basecalling_module=${appendix_path}/model_files/dna_r10.4.1_e8.2_400bps_sup@v5.0.0

# Step 3-1  mapping
mapping_script=${appendix_path}/script/mapping_ont_bulk.sh
reference_genome=${appendix_path}/ref

# Step 3-2 3-3 DESeq2 BAM visualisation
metadata_file=${appendix_path}/meta_data.csv
tracks_file=${appendix_path}/tracks.ini
region_meta=${appendix_path}/region_metadata.csv
contig_meta=${appendix_path}/contig_meta.csv
# check_jobs_shell=${appendix_path}/script/check_successful_task.sh
# Step 5 DESeq2 analysis
degseq2_script=${appendix_path}/script/Deseq2.R

# Other scripts
jobs_check_shell=${appendix_path}/script/Jobs_check.sh


#=================================================================
#+++++++++++++++++++++++Step 1 basecalling +++++++++++++++++++++++
#=================================================================
mkdir ${out_put_path}/1_1_dorado
cd ${out_put_path}/1_1_dorado

# Step 1 generate dorado basecalling script 
ls ${pod5_path} | grep "pod5" | cut -f 1 -d "." | awk -v basecalling_module=${basecalling_module} -v pod5_path=${pod5_path} '{print "dorado basecaller "basecalling_module" "pod5_path"/"$1".pod5 --no-trim --emit-fastq > "$1".fastq"}' > basecalling.sh

# Step 2 submit dorado basecalling jobs
count=1
while read runcode; do
    bsub -q gpu -gpu "num=1/host" -R a100_80g -R "rusage[mem=16GB]" -P ${sample_name}_basecalling_${count} -J ${sample_name}_basecalling_${count} -eo ${sample_name}_basecalling_${count}.err -oo ${sample_name}_basecalling_${count}.out $runcode
    count=$((count+1))
done < basecalling.sh

sh ${jobs_check_shell} -f basecalling.sh -l ${sample_name}_basecalling

# #=================================================================
# #+++++++++++++++++++++++Step 1_2 QC and stat++++++++++++++++++++++
# #=================================================================
mkdir -p ${out_put_path}/1_2_QC_stat
cd ${out_put_path}/1_2_QC_stat
# Step 1 merge dorado fastq files
cat ${out_put_path}/1_1_dorado/*fastq > ${sample_name}.fastq

# Step 2 filter fastq files based on length and quality
echo "NanoFilt -l 100 --maxlength 3000 -q 10 ${sample_name}.fastq > ${sample_name}.filtered.fastq" > ${sample_name}_filter.sh
bsub -P ${sample_name}_filter -J ${sample_name}_filter -n 4 -R "rusage[mem=8GB]" -eo ${sample_name}_filter.err -oo ${sample_name}_filter.out "
sh ${sample_name}_filter.sh"

sh ${jobs_check_shell} -f ${sample_name}_filter.sh -l ${sample_name}_filter

# Step 3 NanoPlot and Seqkit stat
bsub -P ${sample_name}_NanoPlot -J ${sample_name}_NanoPlot -n 4 -R "rusage[mem=32GB]" -eo ${sample_name}_NanoPlot.err -oo ${sample_name}_NanoPlot.out "
NanoPlot --fastq ${sample_name}.fastq -o ${sample_name}_NanoPlot "

bsub -P ${sample_name}_Stats -J ${sample_name}_Stats -n 2 -R "rusage[mem=4GB]" -eo ${sample_name}_Stats.err -oo ${sample_name}_Stats.out "
seqkit stat *.fastq > ${sample_name}_QC_stat.txt"


# #=================================================================
# #+++++++++++++++++++++++Step 2  Demultiplexing +++++++++++++++++++
# #=================================================================
mkdir -p ${out_put_path}/2_Demultiplexing
cd ${out_put_path}/2_Demultiplexing


dorado demux --output-dir ${out_put_path}/2_Demultiplexing --emit-fastq --kit-name SQK-RPB114-24 ${out_put_path}/1_2_QC_stat

# dorado demux --output-dir ${out_put_path}/2_mapping --emit-fastq --kit-name SQK-RPB114-24 PBG67787_skip_fbc5abd6_b3420257_0.fastq
bsub -P stat -J stat -n 2 -R "rusage[mem=8GB]" -eo stat.err -oo stat.out "
seqkit stat *.fastq > Demux_stat.txt
"

#=================================================================
#+++++++++++++++++++++++Step 3-1 mapping +++++++++++++++++++++++++
#=================================================================
mkdir -p ${out_put_path}/3_1_mapping
cd ${out_put_path}/3_1_mapping

du -sh ${out_put_path}/2_Demultiplexing/*fastq | while read line; do
    size=$(echo $line | awk '{print $1}')
    file=$(echo $line | awk '{print $2}')
    file_name=$(basename $file | sed 's/^.*_//g')
    if [[ $size == *G* ]]; then
       ln -s $file ${file_name}
    fi
done

rm *unclassified.fastq
# rename 's/^.*_//' *.fastq
ln -s ${reference_genome}/*fna* .

ls *.fna | while read ref_genome; do
    ls *.fastq | while read fqfile; do
    echo "sh ${mapping_script} -i ${fqfile} -r ${ref_genome}" 
done
done > run_mapping.sh

count=1
while read runcode; do
    bsub -P ${sample_name}_mapping_${count} -J ${sample_name}_mapping_${count} -n 4 -R "rusage[mem=8GB]" -eo ${sample_name}_mapping_${count}.err -oo ${sample_name}_mapping_${count}.out $runcode
    count=$((count + 1))
done < run_mapping.sh

sh ${jobs_check_shell} -f run_mapping.sh -l ${sample_name}_mapping

#=================================================================
#+++++++++++++++++++++++Step 3-2 bamCoverage ++++++++++++++++++++++
#=================================================================
mkdir -p ${out_put_path}/3_2_bam_coverage
cd ${out_put_path}/3_2_bam_coverage


ln -s ${out_put_path}/3_1_mapping/*sorted_reads.bam ./
ln -s ${out_put_path}/3_1_mapping/*sorted_reads.bam.bai ./


ls *bam | while read bamfile; do
    sample_name=$(echo $bamfile | sed 's/_sorted_reads.bam//g')
    echo "bamCoverage -b ${bamfile} -o ${sample_name}.bw --normalizeUsing CPM --binSize 50 --ignoreDuplicates "
done > run_bamCoverage.sh


count=1
while read runcode; do
    bsub -P bamCov_${count} -J bamCov_${count} -n 2 -R "rusage[mem=8GB]" -eo bamCov_${count}.err -oo bamCov_${count}.out $runcode
    count=$((count + 1))
done < run_bamCoverage.sh

sh ${jobs_check_shell} -f run_bamCoverage.sh -l bamCov

#=================================================================
#+++++++++++++++++++++++Step 3-3 visualisation bamfiles ++++++++++
#=================================================================
mkdir -p ${out_put_path}/3_3_bam_visualisation
cd ${out_put_path}/3_3_bam_visualisation

# make_tracks_file --trackFiles *bw  --out tracks.ini 
cat $region_meta | tail -n +2 |cut -f 1 -d ","| sort |uniq | while read ref_id;
do 
    sed 's/ref/'"$ref_id"'/g' ${tracks_file} > tracks_${ref_id}.ini
    grep ${ref_id} $region_meta |awk -F, '{print $3"\t"$4"\t"$5"\t"$2}' > region_${ref_id}.bed
done

ls  ${out_put_path}/3_2_bam_coverage/*bw | while read bwfile; do
    sample_name=$(basename $bwfile | sed 's/\.sorted_reads.*$//g')
    ln -s $bwfile ${out_put_path}/3_3_bam_visualisation/${sample_name}.bw
done

cat $contig_meta | tail -n +2 |cut -f 1 -d ","| sort |uniq | while read ref_id;
do
cat $contig_meta | grep ${ref_id} |cut -f 2 -d "," | while read contig_id;
do
length=$(cat $contig_meta | grep ${ref_id} |grep ${contig_id} |cut -f 3 -d "," |tr -d "\n" )
echo "pyGenomeTracks --tracks tracks_${ref_id}.ini --region ${contig_id}:1-${length} --dpi 300 --outFileName ${ref_id}-${contig_id}_tracks.pdf"
done
done > run_generate_tracks.sh
sh run_generate_tracks.sh

# Zoom in region 2
cat $region_meta | tail -n +2 |cut -f 1 -d ","| sort |uniq | while read ref_id;
do
    sed 's/max_value = 30/max_value = 15/g' tracks_${ref_id}.ini > zoom_in_${ref_id}.ini
    start=$(grep ${ref_id} $region_meta | grep "Region_2" | awk -F, '{print $4}')
    extent_start=$(($start - 3500))
    end=$(grep ${ref_id} $region_meta | grep "Region_2" | awk -F, '{print $5}')
    extent_end=$(($end + 60000))
    contig=$(grep ${ref_id} $region_meta | grep "Region_2" | awk -F, '{print $3}')
    echo "pyGenomeTracks --tracks zoom_in_${ref_id}.ini --region ${contig}:${extent_start}-${extent_end} --dpi 300 --outFileName zoom_in_${ref_id}_region_2.pdf"
done > run_zoom_in_region_2.sh

sh run_zoom_in_region_2.sh

rm *bw 
mkdir log_files
mv *bed *ini log_files/
#=================================================================
#+++++++++++++++++++++++Step 4 Feature Count +++++++++++++++++++++
#=================================================================
mkdir -p ${out_put_path}/4_featureCount
cd ${out_put_path}/4_featureCount
mkdir overlap_CDS nonoverlap_CDS

#overlap_CDS
cd ${out_put_path}/4_featureCount/overlap_CDS
ln -s ${out_put_path}/3_1_mapping/*sorted_reads.bam ./

ls |grep barcode|cut -f 2 -d '.' | sort | uniq | while read ref; do
echo "featureCounts -O -d 30 -D 3000 -t CDS,ncRNA,tmRNA,regulatory_region,oriT,oriC -g ID -a ${reference_genome}/${ref}.gff3 -o ${ref}_gene.count -T 4 -L " $(ls *${ref}.sorted_reads.bam | tr '\n' ' ') " "
done > run_featureCount_overlap.sh

count=1
while read runcode; do
    bsub -P fc_overlap_${count} -J fc_overlap_${count} -n 2 -R "rusage[mem=16GB]" -eo fc_overlap_${count}.err -oo fc_overlap_${count}.out $runcode
    count=$((count + 1))
done < run_featureCount_overlap.sh

sh ${jobs_check_shell} -f run_featureCount_overlap.sh -l fc_overlap


#nonoverlap_CDS
cd ${out_put_path}/4_featureCount/nonoverlap_CDS
ln -s ${out_put_path}/3_1_mapping/*sorted_reads.bam ./

ls |grep barcode|cut -f 2 -d '.' | sort | uniq | while read ref; do
echo "featureCounts -d 30 -D 3000 -t CDS,ncRNA,tmRNA,tRNA,regulatory_region,oriT,oriC -g ID -a ${reference_genome}/${ref}.gff3 -o ${ref}_gene.count -T 4 -L " $(ls *${ref}.sorted_reads.bam | tr '\n' ' ') " "
done > run_featureCount_nonoverlap.sh


count=1
while read runcode; do
    bsub -P fc_${count} -J fc_${count} -n 2 -R "rusage[mem=16GB]" -eo fc_${count}.err -oo fc_${count}.out $runcode
    count=$((count + 1))
done < run_featureCount_nonoverlap.sh

sh ${jobs_check_shell} -f run_featureCount_nonoverlap.sh -l fc
#=================================================================
#+++++++++++++++++++++++Step 5 DESeq2 analysis +++++++++++++++++++
#=================================================================
mkdir -p ${out_put_path}/5_DESeq2
cd ${out_put_path}/5_DESeq2



mkdir -p  ${out_put_path}/5_DESeq2/count_table/nonoverlap_CDS  ${out_put_path}/5_DESeq2/count_table/overlap_CDS 
cp ${out_put_path}/4_featureCount/nonoverlap_CDS/*gene.count ${out_put_path}/5_DESeq2/count_table/nonoverlap_CDS/
cp ${out_put_path}/4_featureCount/overlap_CDS/*gene.count ${out_put_path}/5_DESeq2/count_table/overlap_CDS/

ls ${out_put_path}/5_DESeq2/count_table/overlap_CDS  | sed 's/_gene.count//g' | while read sample; do
    sed -i "s/\.${sample}\.sorted_reads\.bam//g" ${out_put_path}/5_DESeq2/count_table/nonoverlap_CDS/*gene.count
    sed -i "s/\.${sample}\.sorted_reads\.bam//g" ${out_put_path}/5_DESeq2/count_table/overlap_CDS/*gene.count
done



ls ${out_put_path}/5_DESeq2/count_table/nonoverlap_CDS/ | while read line; do
    sample_name=$(echo $line | sed 's/_gene.count//g')     
    echo "Rscript ${degseq2_script} -c ${out_put_path}/5_DESeq2/count_table/nonoverlap_CDS/${line} -m ${metadata_file} -r ${sample_name} -g ${reference_genome}/${sample_name}.gff3 -o ${out_put_path}/5_DESeq2/nonoverlap_CDS_results"
done > run_deseq2.sh

ls ${out_put_path}/5_DESeq2/count_table/overlap_CDS/ | while read line; do
    sample_name=$(echo $line | sed 's/_gene.count//g')     
    echo "Rscript ${degseq2_script} -c ${out_put_path}/5_DESeq2/count_table/overlap_CDS/${line} -m ${metadata_file} -r ${sample_name} -g ${reference_genome}/${sample_name}.gff3 -o ${out_put_path}/5_DESeq2/overlap_CDS_results"
done >> run_deseq2.sh

count=1
while read runcode; do
    bsub -P deseq2_${count} -J deseq2_${count} -n 2 -R "rusage[mem=8GB]" -eo deseq2_${count}.err -oo deseq2_${count}.out $runcode
    count=$((count + 1))
done < run_deseq2.sh

sh ${jobs_check_shell} -f run_deseq2.sh -l deseq2

#=================================================================
#+++++++++++++++++++++++Step 6 Summarize analysis ++++++++++++++++
#=================================================================
mkdir -p ${out_put_path}/6_Summarize_analysis
cd ${out_put_path}/6_Summarize_analysis
cp ${out_put_path}/1_2_QC_stat/${sample_name}_NanoPlot/NanoStats.txt 1_1_NanoPlot_stats.txt
cp ${out_put_path}/1_2_QC_stat/${sample_name}_QC_stat.txt 1_2_reads_QC_stat.txt
cat ${out_put_path}/2_Demultiplexing/Demux_stat.txt | sed 's/^.*barcode/barcode/g' > 2_Demultiplexing_stat.txt

echo -e "sample_name\treference\tmapped_reads\tprimary_mapped" > 3_1_mapping_stats.txt
ls ${out_put_path}/3_1_mapping/*alignment_stats.txt | while read line; do
    file_name=$(basename $line | sed 's/.alignment_stats.txt//g')
    sample_name=$(cut -f 1 -d '.' <<< "$file_name")
    reference=$(cut -f 2 -d '.' <<< "$file_name")
    mapped_reads=$(grep "mapped (" $line|grep -v "primary mapped" | awk -F':' '{print $1}'|sed 's/^.*(//g' )
    primary_mapped=$(grep "primary mapped" $line | awk -F':' '{print $1}'|sed 's/^.*(//g' )
    
    echo -e "${sample_name}\t${reference}\t${mapped_reads}\t${primary_mapped}"  >> 3_1_mapping_stats.txt
done


mkdir -p 4_featureCount_summary/overlap_CDS 4_featureCount_summary/nonoverlap_CDS
cp ${out_put_path}/4_featureCount/overlap_CDS/*summary ${out_put_path}/6_Summarize_analysis/4_featureCount_summary/overlap_CDS/
cp ${out_put_path}/4_featureCount/nonoverlap_CDS/*summary ${out_put_path}/6_Summarize_analysis/4_featureCount_summary/nonoverlap_CDS/

sed -i 's/\.sorted_reads\.bam//g' 4_featureCount_summary/*/*summary
