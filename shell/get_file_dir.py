#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : get_file_dir.py
@Time     : 2022/03/26 12:23:37
@Author   : RanPeng 
@Version  : python 3.8.6
@Contact  : 2502388440@hotmail.com
@License  : (C)Copyright 2021-2022, DingLab-CHINA-SHANGHAI
@Function : None
'''
# here put the import lib
import os
import pandas as pd
import numpy as np



fastq_repos_dict = {
                    'Cohort_No' : [],
                    'fq1' : [],
                    'fq2' : [],
                    'Sample' : []
                    }
# 
def find_all_files(input_folder):
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.endswith('_R1.fastq.gz') or file.endswith('_R2.fastq.gz'):
                if file.split('_')[0] not in fastq_repos_dict['Cohort_No']:
                    fastq_repos_dict['Cohort_No'].append(file.split('_')[0])
                    fastq_repos_dict['Sample'].append(file.split('_')[0])
                if file.endswith('_R1.fastq.gz') and (file.split('_')[0] in file):
                    if file not in list(map(lambda x: os.path.split(x)[-1], fastq_repos_dict['fq1'])):
                        fastq_repos_dict['fq1'].append(os.path.join(root, file))
                    else:
                        print(f'\n{file}: duplicate fq1 file')
                        print(f'The path to the duplicate fq1 file is {os.path.join(root, file)}\n')
                elif file.endswith('_R2.fastq.gz') and (file.split('_')[0] in file):
                    if file not in list(map(lambda x: os.path.split(x)[-1], fastq_repos_dict['fq2'])):
                        fastq_repos_dict['fq2'].append(os.path.join(root, file))
                    else:
                        print(f'{file}: duplicate fq2 file')
                        print(f'The path to the duplicate fq2 file is {os.path.join(root, file)}')
                else:
                    fastq_repos_dict['fq1'].append(np.nan)
                    fastq_repos_dict['fq2'].append(np.nan)
                    print(f'{file} exists an error!')
    print(fastq_repos_dict)
    fastq_repos_dict1 = pd.DataFrame(fastq_repos_dict)
    # fastq_repos_dict1 = pd.DataFrame().from_dict(fastq_repos_dict, orient='index').T
    fastq_repos_dict1.to_excel('fastq_repos_dirs.xlsx', index=False)
    print("Successfully created output file: {}".format('fastq_repos_dirs.xlsx'))


if __name__ == '__main__':
    import sys
    input_folder = sys.argv[1]
    find_all_files(input_folder)