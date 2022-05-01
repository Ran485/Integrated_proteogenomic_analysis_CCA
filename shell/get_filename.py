#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : get_filename.py
@Time     : 2022/03/18 22:48:54
@Author   : RanPeng 
@Version  : python 3.8.6
@Contact  : 2502388440@hotmail.com
@License  : (C)Copyright 2021-2022, DingLab-CHINA-SHANGHAI
@Function : None
'''
# here put the import lib
import os
import enum
import pandas as pd


# Enum for size units
class SIZE_UNIT(enum.Enum):
    BYTES = 1
    KB = 2
    MB = 3
    GB = 4


def convert_unit(size_in_bytes, unit):
    """ Convert the size from bytes to other units like KB, MB or GB"""
    if unit == SIZE_UNIT.KB:
        return size_in_bytes / 1024
    elif unit == SIZE_UNIT.MB:
        return size_in_bytes / (1024 * 1024)
    elif unit == SIZE_UNIT.GB:
        return size_in_bytes / (1024 * 1024 * 1024)
    else:
        return size_in_bytes


def get_file_size(file_name, size_type=SIZE_UNIT.BYTES):
    """ Get file in size in given unit like KB, MB or GB"""
    size = os.path.getsize(file_name)
    return convert_unit(size, size_type)


# get all the filename and the file size of a directory, save as a csv file
def get_filename(dir, format):
    '''
    dir: the directory you want to get the filename
    format: the format of the filename you want to get
    example:
    get_filename(dir='/Users/ranpeng/Desktop/', format='csv')
    :return:
    '''
    filename_list = []
    filesize_list = []
    for root, dirs, files in os.walk(dir):
        for file in files:
            if file.endswith(format):
                filename_list.append(file)
                filesize_list.append(round(get_file_size(os.path.join(root, file), SIZE_UNIT.MB), 4))
                # filesize_list = [convert_unit(filesize, 2) for filesize in filesize_list]
    pd.DataFrame(zip(filename_list, filesize_list), columns=['filename', 'filesize']).to_csv(os.path.join(dir, 'filename_list.csv'), index=True)
    print(f'filename_list.csv has been successfully saved in {dir}!')


if __name__ == '__main__':
    import sys
    dir = sys.argv[1]
    format = sys.argv[2]
    # dir = r/Users/ranpeng/Desktop/CCA/CCA_script/python_script/gene_fusion/'
    format = format.strip()
    get_filename(dir=dir, format=format)