#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : fusion_data_merge.py
@Time     : 2022/04/28 21:40:39
@Author   : RanPeng 
@Version  : pyhton 3.8.6
@Contact  : 2502388440@hotmail.com
@License  : (C)Copyright 2019-2020, DingLab-CHINA-SHNAGHAI
@Function : None
'''
# here put the import lib
import os
import pandas as pd


def get_confidence_medium(data):
    data_CI_medium = data.loc[data['confidence'] != 'low']
    return data_CI_medium


def find_all_files(input_folder):
    fuison_events = pd.DataFrame()
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.endswith('fusions.tsv'):
                filepath = os.path.join(root, file)
                sample = filepath.split('/')[-2].split('_')[0]
                print(sample)
                fusions = pd.read_table(filepath)
                # print(fusions)
                fusions = get_confidence_medium(fusions)
                fusions.insert(0,'sample',sample) ## insert sample infomation
                fuison_events = fuison_events.append(fusions)
    return fuison_events


def filter_kinase_fusions(data, kinase_file):
    kinase_list = pd.read_excel(kinase_file)
    kinase_list = kinase_list['gene'].tolist()
    data1 = data.loc[data['#gene1'].isin(kinase_list)]
    data2 = data.loc[data['gene2'].isin(kinase_list)]
    fusion_kinase = pd.concat([data1, data2])
    return fusion_kinase


if __name__ == '__main__':
    input_folder = r'/Volumes/Samsung_T5/downloads/arriba_fusion/CCA_fusion_res/'
    outpath = r'/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_data/other_drug_targrt_kinase/'
    kinase_file = r'/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/FGFR2_data/other_drug_targrt_kinase/Approved_Drug_Kinase.xlsx'
    fuison_events_merge = find_all_files(input_folder)
    fusion_kinase = filter_kinase_fusions(fuison_events_merge, kinase_file)
    fuison_events_merge.to_excel(outpath + 'fuison_events_merge.xlsx')
    fusion_kinase.to_excel(outpath + 'fusion_kinase.xlsx')