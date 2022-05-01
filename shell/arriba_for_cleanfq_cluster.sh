#!/bin/bash
## arriba_for_cleanfq_cluster.sh
## Specify the path to the software and reference database
STAR_index_idx='/public/home/galaxy/wes_cancer/gene_fusion/biosoft/starFusion/db/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx'
ref_annot_gtf='/public/home/galaxy/wes_cancer/gene_fusion/biosoft/starFusion/db/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf'
ref_genome_fa='/public/home/galaxy/wes_cancer/gene_fusion/biosoft/starFusion/db/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa'
arriba='/public/home/galaxy/wes_cancer/gene_fusion/biosoft/arriba_v2.2.1'
draw_fusions_shell='/public/home/galaxy/wes_cancer/gene_fusion/biosoft/arriba_fusion_output/shell_script'

## 获取参数 sample, raw_fastq_dir, clean_fastq_dir
while getopts ":s:a:c:t:" opt
do
    case $opt in
        s)
        sample=${OPTARG}
                echo $sample
        ;;
        a)
        arriba_fusions_result_dir=${OPTARG}
                echo $arriba_fusions_result_dir
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
## Specify the folder where the fastq file is located
file_dir=$clean_fastq_dir
# arriba_fusions_result_dir=$arriba_fusions_result_dir
# 检查文件目录是否存在，不存在则创建
if [ ! -d $arriba_fusions_result_dir ];then
    mkdir $arriba_fusions_result_dir
else
    echo "文件目录已经存在"
fi
## Add file suffixes
suffix='_arriba_result'
cd $arriba_fusions_result_dir
sample_outpath=$sample$suffix
# 检查sample的输出路径是否存在
if [ ! -d $sample_outpath ];then
    mkdir $sample_outpath
    echo "成功创建 ${sample_outpath} 文件夹!"
else
    echo "${sample_outpath} 文件目录已经存在!"
fi
# 切换到`sample_outpath`目录下
cd $sample_outpath
echo "当前目录: " `pwd`
# 判断之前是否已经完成融合的检测
if [ ! -f $sample_outpath/fusions.tsv ] || [ ! -f $sample_outpath/fusions.pdf ];then
    echo "arriba for ${sample} is start running !"
    startTime=`date +%Y%m%d-%H:%M`
    startTime_s=`date +%s`
    file=$( basename $sample )
    fq1=$file_dir/${file}_R1_val_1.fq.gz
    fq1=$file_dir/${file}_R2_val_2.fq.gz
    $arriba/run_arriba.sh $STAR_index_idx/ $ref_annot_gtf $ref_genome_fa $arriba/database/blacklist_hg19_hs37d5_GRCh37_v2.2.1.tsv.gz $arriba/database/known_fusions_hg19_hs37d5_GRCh37_v2.2.1.tsv.gz $arriba/database/protein_domains_hg19_hs37d5_GRCh37_v2.2.1.gff3 $thread_nums $fq1 $fq2
    # draw_fusions circos plots
    $draw_fusions_shell/draw_fusions_cluster.sh
    # 由于bam文件较大，用完之后需要及时清理
    rm -f Aligned.sortedByCoord.out.bam
    cd ../
    echo "arriba for ${sample} succesfully done !"
    endTime=`date +%Y%m%d-%H:%M`
    endTime_s=`date +%s`
    sumTime=$[ $endTime_s - $startTime_s ]
    echo "$startTime ---> $endTime" "Total:$sumTime seconds."
    # 换行
    echo -e '\n\n\n\n\n\n'
else
    echo "${sample}; 样本基因融合检测已经完成!"
    echo "命令终止,退出执行～"
fi

