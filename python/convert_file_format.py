import pandas as pd
import numpy as np
import time
import os

def covert_file_format(path=None, outpath=None, format=None):
    """对数据进行清洗"""
    for root, dirs, files in os.walk(path, topdown=False):
        for f in files:
            filename = f.split('.')[0]
            if f.endswith(format):
                if format == '.csv':
                    df = pd.read_csv(os.path.join(root, f))
                elif format == '.xlsx':
                    df = pd.read_excel(os.path.join(root, f))
                else:
                    df = pd.read_table(os.path.join(root, f))
                # df.drop(df.index[0:98], inplace=True)
                df.to_excel(outpath + filename + '.xlsx', index=False)


if __name__ == '__main__':
    path = '/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/DPCR1/'
    outpath = '/Volumes/Samsung_T5/downloads/Integrated_multiomics_analysis_CCA/data/response/DPCR1/pathway_data/'
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    covert_file_format(path = path, outpath = outpath, format='.tab')