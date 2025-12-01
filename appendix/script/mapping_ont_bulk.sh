#!/bin/bash

# Parse command line options using getopt
OPTS=$(getopt -o i:r:h --long input:,reference:,help -n 'mapping_ont_bulk.sh' -- "$@")

if [ $? != 0 ]; then
    echo "Failed parsing options." >&2
    exit 1
fi

eval set -- "$OPTS"

# Default values
input_fastq=""
ref_genome=""

# Function to display help
show_help() {
    echo "Usage: $0 -i INPUT_FASTQ -r REF_GENOME"
    echo "  -i, --input      Input FASTQ file"
    echo "  -r, --reference  Reference genome file"
    echo "  -h, --help       Show this help message"
    exit 0
}

# Parse options
while true; do
    case "$1" in
        -i | --input)
            input_fastq="$2"
            shift 2
            ;;
        -r | --reference)
            ref_genome="$2"
            shift 2
            ;;
        -h | --help)
            show_help
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Invalid option: $1" >&2
            exit 1
            ;;
    esac
done

# Check if required arguments are provided
if [ -z "$input_fastq" ] || [ -z "$ref_genome" ]; then
    echo "Error: Both input FASTQ file and reference genome are required."
    show_help
fi

#=================================================================
#+++++++++++++++++++++++Step 1 paramter ++++++++++++++++++++++++++
#=================================================================
fastq_name=$(echo ${input_fastq} | cut -d "." -f 1)
ref_name=$(echo ${ref_genome} | cut -d "." -f 1)
sample_name="${fastq_name}.${ref_name}"

echo "================= Parameters ================="
echo "input_fastq: ${input_fastq}"
echo "ref_genome: ${ref_genome}"
echo "sample_name: ${sample_name}"
echo "============================================="
#=================================================================
#+++++++++++++++++++++++Step 2 Mapping +++++++++++++++++++++++++++
#=================================================================


# Step 0: Sort the extracted FASTQ file
fastp -i ${input_fastq} -o ${sample_name}.sort.fastq -l 30 -j ${sample_name}.sorted.json -h ${sample_name}.sorted.html

# Step 1: Align Nanopore FASTQ reads to the reference genome
minimap2 -ax map-ont ${ref_genome} ${sample_name}.sort.fastq --secondary=no > ${sample_name}_aligned_reads.sam
#-a output in the SAM format (PAF by default)
#-x map-ont preset for mapping ONT reads
#--secondary=no do not output secondary alignments

# Step 2: Convert SAM to BAM and sort
samtools view -S -b ${sample_name}_aligned_reads.sam | samtools sort -o ${sample_name}.sorted_reads.bam

# Step 3: Index the sorted BAM file
samtools index ${sample_name}.sorted_reads.bam

# Step 4: (Optional) View alignment stats
samtools flagstat ${sample_name}.sorted_reads.bam > ${sample_name}.alignment_stats.txt