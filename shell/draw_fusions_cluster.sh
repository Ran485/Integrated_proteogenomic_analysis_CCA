#!/bin/bash
# draw_fusions_cluster.sh
ref_annot_gtf='/public/home/galaxy/wes_cancer/gene_fusion/biosoft/starFusion/db/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf'
arriba='/public/home/galaxy/wes_cancer/gene_fusion/biosoft/arriba_v2.2.1'
## 激活R环境
source activate R4.1.3 

startTime=`date +%Y%m%d-%H:%M`
startTime_s=`date +%s`

$arriba/draw_fusions.R \
    --fusions=fusions.tsv \
    --alignments=Aligned.sortedByCoord.out.bam \
    --output=fusions.pdf \
    --annotation=$ref_annot_gtf \
    --cytobands=$arriba/database/cytobands_hg19_hs37d5_GRCh37_v2.2.1.tsv \
    --proteinDomains=$arriba/database/protein_domains_hg19_hs37d5_GRCh37_v2.2.1.gff3
# 退出当前环境
conda deactivate
endTime=`date +%Y%m%d-%H:%M`
endTime_s=`date +%s`
sumTime=$[ $endTime_s - $startTime_s ]
echo "$startTime ---> $endTime" "Total:$sumTime seconds."
# 换行
echo -e '\n\n\n\n\n\n'
