#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : CCA_data_mapping.py
@Time     : 2022/04/28 14:54:51
@Author   : RanPeng 
@Version  : pyhton 3.8.6
@Contact  : 2502388440@hotmail.com
@License  : (C)Copyright 2019-2020, DingLab-CHINA-SHNAGHAI
@Function : None
'''
# here put the import lib


import pandas as pd
import numpy as np
import time
import os


def merge_multi_omics_data(metadata_path=None,
                           df_mapping=None, 
                           outpath=None,
                           mutation=False,
                           CNA=False,
                           RNA=False,
                           Protein_profile=False,
                           Protein_profile_TN_ratio=False,
                           Phosphosite=False,
                           Phosphoprotein=False,
                           time_name=None,
                           table_name=None):
    if mutation:
    ## 基因突变的数据
        gene_mutation = pd.read_csv(
            metadata_path + 'Gene_Mutation_top15000_matrix.csv', header=0, index_col=0)
        gene_mutation = gene_mutation.T
        print(gene_mutation)

        merge_gene_mutation = pd.merge(
            df_mapping, gene_mutation, left_on='gene_T', right_index=True, how="left")
        merge_gene_mutation = merge_gene_mutation.T
        print(merge_gene_mutation)
        merge_gene_mutation.to_csv(
            outpath + table_name + '_' + 'merge_gene_mutation_'+time_name+'.csv')
    if CNA:
        ## CNV的数据
        CNV = pd.read_csv(metadata_path + 'CNA_All.csv', header=0, index_col=0)
        CNV = CNV.T
        print(CNV)

        merge_CNV = pd.merge(df_mapping, CNV, left_on='gene_T',
                            right_index=True, how="left")
        merge_CNV = merge_CNV.T
        print(merge_CNV)
        merge_CNV.to_csv(outpath + table_name + '_' +
                        'merge_CNV_'+time_name+'.csv')

        ## # Arm_focal  amplication deletion的数据
        ## Amplitude Threshold
        Focal_arm_events_Threshold = pd.read_csv(
            metadata_path + 'Focal_arm_events_Threshold.csv', header=0, index_col=0)
        Focal_arm_events_Threshold = Focal_arm_events_Threshold.T
        print(Focal_arm_events_Threshold)

        merge_Focal_arm_events_Threshold = pd.merge(
            df_mapping, Focal_arm_events_Threshold, left_on='gene_T', right_index=True, how="left")
        merge_Focal_arm_events_Threshold = merge_Focal_arm_events_Threshold.T
        print(merge_Focal_arm_events_Threshold)
        merge_Focal_arm_events_Threshold.to_csv(
            outpath + table_name + '_' + 'merge_Focal_arm_events_Threshold_'+time_name+'.csv')

        ## Actual Copy Change Given
        Focal_arm_events_CN_value = pd.read_csv(
            metadata_path + 'Focal_arm_events_CN_value.csv', header=0, index_col=0)
        Focal_arm_events_CN_value = Focal_arm_events_CN_value.T
        print(Focal_arm_events_CN_value)

        merge_Focal_arm_events_CN_value = pd.merge(
            df_mapping, Focal_arm_events_CN_value, left_on='gene_T', right_index=True, how="left")
        merge_Focal_arm_events_CN_value = merge_Focal_arm_events_CN_value.T
        print(merge_Focal_arm_events_CN_value)
        merge_Focal_arm_events_CN_value.to_csv(
            outpath + table_name + '_' + 'merge_Focal_arm_events_CN_value_'+time_name+'.csv')

        ## Chromosome Arm 的数据
        Chromosome_Arm = pd.read_csv(
            metadata_path + 'broad_values_by_arm_All.csv', header=0, index_col=0)
        Chromosome_Arm = Chromosome_Arm.T
        print(Chromosome_Arm)

        merge_Chromosome_Arm = pd.merge(
            df_mapping, Chromosome_Arm, left_on='gene_T', right_index=True, how="left")
        merge_Chromosome_Arm = merge_Chromosome_Arm.T
        print(merge_Chromosome_Arm)

        merge_Chromosome_Arm.to_csv(
            outpath + table_name + '_' + 'merge_Chromosome_Arm_'+time_name+'.csv')
    if RNA:
        ## RNA 的数据
        RNA = pd.read_csv(metadata_path + 'RNA_T_del50%.csv', #RNA_FPKM_codegene_deduplicate_del20%
                        header=0, index_col=0)
        # RNA = pd.read_csv(metadata_path + 'RNA_FPKM_codegene_deduplicate_del20%.csv', 
        #                 header=0, index_col=0)
        RNA = RNA.T
        print(RNA)
        # df_mapping["firT"] = df_mapping["firT"].astype("str")
        # protein["firT"] = protein["firT"].astype("int64")

        merge_RNA = pd.merge(df_mapping, RNA, left_on='RNA_T',
                            right_index=True, how="left")
        merge_RNA = merge_RNA.T
        print(merge_RNA)
        merge_RNA.to_csv(outpath + table_name + '_' +
                        'merge_RNA_Delna50%_'+time_name+'.csv')
    if Protein_profile:
        ## profiling proteome 的数据
        # protein = pd.read_csv(
        #     metadata_path + 'protein_Tumor_217patients_Delna30%.csv', header=0, index_col=0)
        protein = pd.read_csv(
            metadata_path + 'protein_Tumor_217patients_15901_genes_all.csv', header=0, index_col=0)
        protein = protein.T 
        print(protein) # protein_Tumor_217patients_15901_genes
        df_mapping["firT"] = df_mapping["firT"].astype("str")
        # protein["firT"] = protein["firT"].astype("int64")

        merge_protein = pd.merge(
            df_mapping, protein, left_on='firT', right_index=True, how="left")
        merge_protein = merge_protein.T
        print(merge_protein)
        merge_protein.to_csv(outpath + table_name + '_' +
                            'merge_protein_Delna30%_1_'+time_name+'.csv')
    if Protein_profile_TN_ratio:
        ## TN ratio 数据
        protein_log2_TN_ratio = pd.read_csv(
            metadata_path + '197_TN_log2ratio_sort.csv', header=0, index_col=0)
        protein_log2_TN_ratio = protein_log2_TN_ratio.T
        print(protein_log2_TN_ratio)
        df_mapping["firT"] = df_mapping["firT"].astype("str")
        # protein["firT"] = protein["firT"].astype("int64")

        merge_protein_log2_TN_ratio = pd.merge(
            df_mapping, protein_log2_TN_ratio, left_on='firT', right_index=True, how="left")
        merge_protein_log2_TN_ratio = merge_protein_log2_TN_ratio.T
        print(merge_protein_log2_TN_ratio)
        merge_protein_log2_TN_ratio.to_csv(
            outpath + table_name + '_' + 'protein_log2_TN_ratio_'+time_name+'.csv')
    if Phosphosite:
        ## phosphosite 的数据
        phosphosite = pd.read_csv(
            metadata_path + 'phosphosite_dropna20_FOT_20220410.csv', header=0, index_col=0)
        phosphosite = phosphosite.T
        print(phosphosite)

        merge_phosphosite = pd.merge(
            df_mapping, phosphosite, left_on='phos_T', right_index=True, how="left")
        merge_phosphosite = merge_phosphosite.T
        print(merge_phosphosite)
        merge_phosphosite.to_csv(outpath + table_name +
                                '_' + 'merge_phosphosite_'+time_name+'.csv')
    if Phosphoprotein:
        ## phosphoprotein 的数据
        phosphoprotein = pd.read_csv(
            metadata_path + 'phosphoprotein_abudance_FOT(2021-2-3).csv', header=0, index_col=0)
        phosphoprotein = phosphoprotein.T
        print(phosphoprotein)

        merge_phosphoprotein = pd.merge(
            df_mapping, phosphoprotein, left_on='phos_T', right_index=True, how="left")
        merge_phosphoprotein = merge_phosphoprotein.T
        print(merge_phosphoprotein)
        merge_phosphoprotein.to_csv(
            outpath + table_name + '_' + 'merge_phosphoprotein_'+time_name+'.csv')


def data_wash(path, outpath):
    ## 对数据进行清洗
    for root, dirs, files in os.walk(path, topdown=False):
        for f in files:
            if f.endswith('.csv'):
                df = pd.read_csv(os.path.join(root, f))
                df.drop(df.index[0:98], inplace=True)
                df.to_csv(os.path.join(outpath, f), index=False)


if __name__ == "__main__":
    ## import and output path
    ## inpath = '/Users/ranpeng/Desktop/2020-10-26/amp_genes.conf_90.txt'
    metadata_path = r'/Volumes/Samsung_T5/downloads/CCA/CCA_mapping_data/'
    outpath = r'/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/DPCR1/data/'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    else:
        print('file outpath existed')
    time_name = time.strftime('%Y%m%d', time.localtime(time.time()))
    print(time_name)
    ## -------------- clinical information的表格 -------------------
    clinical_data = pd.read_csv("/Volumes/Samsung_T5/downloads/CCA/CCA_mapping_data/clinical_merge20220406.csv", index_col=0)
    ## -------------- 需要mapping的分类的表格 -------------------
    mapping_index = pd.read_csv(
        "/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/DPCR1/DPCR1_index.csv", index_col=0)
    mapping_index = pd.merge(mapping_index, clinical_data,
                            left_index=True, right_index=True, how="left")
    print(mapping_index)
    ## ------------ mapping table name -----------------
    table_name = 'DPCR1_Mut'
    ## ------------ run the merge function -----------------
    merge_multi_omics_data(metadata_path=metadata_path,
                           df_mapping=mapping_index, 
                           outpath=outpath,
                           mutation=False,
                           CNA=False,
                           RNA=True,
                           Protein_profile=True,
                           Protein_profile_TN_ratio=False,
                           Phosphosite=True,
                           Phosphoprotein=False,
                           time_name=time_name,
                           table_name=table_name)
    ## ------------ data wash -----------------
    data_wash(path=outpath, outpath=outpath)
    print('Successfully mapped all data')
