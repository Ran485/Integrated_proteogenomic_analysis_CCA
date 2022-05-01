#!/bin/bash
# use conda
# file: CCA_gene_fusion.pbs
#PBS -N Arriba-Fusion
#PBS -l nodes=1:ppn=5
#PBS -l mem=10gb
#PBS -l walltime=48:00:00
#PBS -q default
#PBS -j oe
# import all environments virables
#PBS -V

# Enter job's working directory
# Raw fastq data and reference database dir
raw_fastq_dir='/public/home/galaxy/wes_cancer/gene_fusion/biosoft/starFusion/CCA_fusion_data'
clean_fastq_dir='/public/home/galaxy/wes_cancer/gene_fusion/biosoft/arriba_fusion_output/CCA_clean_fq'
arriba_fusions_result_dir='/public/home/galaxy/wes_cancer/gene_fusion/biosoft/arriba_fusion_output/CCA_fusion_res'
# The shell script file directory
shell_dir='/public/home/galaxy/wes_cancer/gene_fusion/biosoft/arriba_fusion_output/shell_script'
# Set the CPU working cores numbers
thread_nums=5
# sample='R21026748-F18-18-176_combined'
# run trim_galore.sh function 质控和去接头,需要base环境运行
echo "Arriba-Fusion for ${sample} is start running !"
$shell_dir/trim_galore_cluster.sh -s $sample -r $raw_fastq_dir -c $clean_fastq_dir -t $thread_nums
# Activating the `gene_fusion` virtual environment
source activate
conda activate gene_fusion
# run arriba_for_cleanfq.sh  对质控过后的fq文件进行基因融合的calling
$shell_dir/arriba_for_cleanfq_cluster.sh -s $sample -a $arriba_fusions_result_dir -c $clean_fastq_dir -t $thread_nums
# Exit the virtual environment
conda deactivate
echo "Arriba-Fusion for ${sample} is succesfully done !"
