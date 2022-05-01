#!/bin/bash
# use conda
# file: example.pbs
#PBS -N STAR-Fusion
#PBS -l nodes=1:ppn=12
#PBS -l mem=50gb
#PBS -l walltime=240:00:00
#PBS -q default
#PBS -j oe
# import all environments virables
#PBS -V

# enter job's working directory
cd /public/home/galaxy/wes_cancer/gene_fusion/biosoft/starFusion/DLBC_fusion
source activate
conda activate gene_fusion
# run STAR-Fusion function
echo "STAR-Fusion for ${sample} is start running !"
# 服务器地址
refer_lib='/public/home/galaxy/wes_cancer/gene_fusion/biosoft/starFusion/db/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/'
out_dir='/public/home/galaxy/wes_cancer/gene_fusion/biosoft/starFusion/DLBC_fusion_outdir_20220303'

startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`
# echo "STAR-Fusion for ${sample} is start running !"
STAR-Fusion --genome_lib_dir $refer_lib \
            --CPU 10 \
            --left_fq ${sample}_R1.fastq.gz \
            --right_fq ${sample}_R2.fastq.gz \
            --output_dir $out_dir/${sample}_gene_fusion \
            --FusionInspector validate
echo "STAR-Fusion for ${sample} succesfully done !"
endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`
sumTime=$[ $endTime - $startTime ]
echo "$startTime ---> $endTime" "Total:$sumTime minutes."
# 换行
echo -e '\n\n\n\n\n\n' 