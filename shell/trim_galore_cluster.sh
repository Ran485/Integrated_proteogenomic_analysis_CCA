#!/bin/bash
## trim_galore_cluster.sh
# Specify the folder where the fastq raw file is located
# 获取参数 sample, raw_fastq_dir, clean_fastq_dir
while getopts ":s:r:c:t:" opt
do
    case $opt in
        s)
        sample=${OPTARG}
				echo $sample
        ;;
        r)
        raw_fastq_dir=${OPTARG}
				echo $raw_fastq_dir
        ;;
        c)
        clean_fastq_dir=${OPTARG}
				echo $clean_fastq_dir
				;;
        t)
        thread_nums=${OPTARG}
				echo $thread_nums
        ;;
        ?)
        echo "未知参数"
        exit 1;;
    esac
done
# raw_fastq_dir=$raw_fastq_dir
# clean_fastq_dir=$clean_fastq_dir
if [ ! -d $clean_fastq_dir ];then
    mkdir $clean_fastq_dir
else
    echo "${clean_fastq_dir} 文件目录已经存在!"
fi
# 判断质控后的文件是否存在，存在则跳过`质控和去接头`这个步骤
cleaned_fq1=$clean_fastq_dir/${sample}_R1_val_1.fq.gz
cleaned_fq2=$clean_fastq_dir/${sample}_R2_val_2.fq.gz
if [ ! -f $cleaned_fq1 ] || [ ! -f $cleaned_fq2 ];then
    echo "trim_galore_cluster for ${sample} is start running !"
    fq1=$raw_fastq_dir/${sample}_R1.fastq.gz
    fq2=$raw_fastq_dir/${sample}_R2.fastq.gz
    trim_galore  --paired -q 20 --phred33 --length 20 --stringency 3 -e 0.1 --gzip --cores $thread_nums -o $clean_fastq_dir/  $fq1  $fq2 >> $clean_fastq_dir/${sample}_trim.log 2>&1
    echo "trim_galore_cluster for ${sample} is succesfully done !"
else
    echo "${cleaned_fq1}; ${cleaned_fq1} 文件都已经存在!"
fi
