#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : EDA.py
@Time     : 2022/03/11 18:37:45
@Author   : RanPeng 
@Version  : python 3.8.6
@Contact  : 2502388440@hotmail.com
@License  : (C)Copyright 2021-2022, DingLab-CHINA-SHANGHAI
@Function : None
'''
# here put the import lib
import pandas as pd
import numpy as np
# import lux
import matplotlib.pyplot as plt
import seaborn as sns
import os, gc, re, warnings, sys
import missingno as msno

def split_data_dtypes(data):
    """查看特征的数值类型有哪些，对象类型有哪些"""
    numerical_columns = list(data.select_dtypes(exclude=['object']).columns)
    category_columns = list(
        filter(lambda x: x not in numerical_columns, list(data.columns)))
    print(f'There are {len(numerical_columns)} numerical columns in dataset.')
    print(f'There are {len(category_columns)} category columns in dataset.')
    print('The list of numerical columns are : {}'.format(numerical_columns))
    print('The list of category columns are : {}'.format(category_columns))
    return numerical_columns, category_columns


def split_numerical_serial_fea(data, features):
    """划分数值型变量中的连续变量和离散型变量"""
    numerical_serial_fea = []
    numerical_noserial_fea = []
    for fea in features:
        temp = data[fea].nunique()
        if temp <= 10:
            numerical_noserial_fea.append(fea)
            continue
        numerical_serial_fea.append(fea)
    print(f'The numerical serial features in the dataset are: \n{numerical_serial_fea}.\n')
    print(f'The category variables features in dataset are: \n{numerical_noserial_fea}.')
    return numerical_serial_fea, numerical_noserial_fea