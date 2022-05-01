import pandas as pd
import numpy as np
import time


def read_file(input_path, index_col=0):
    # excel file use 'pd.read_excel', csv/txt file use 'pd.read_csv'
    if input_path.endswith(".xlsx"):
        df = pd.read_excel(input_path, index_col=index_col)
    elif input_path.endswith(".csv"):
        df = pd.read_csv(input_path, index_col=index_col)
    else:
        df = pd.read_csv(input_path, sep="\t", index_col=index_col)
    return df


def mapping_data(anatation, input_data, output_dir, file_name):
    """
    Returns a dataframe arranged by anatation information
    """
    anatation = read_file(anatation, index_col=0)
    input_data = read_file(input_data, index_col=0)
    # input_data = input_data.T

    merge = pd.merge(anatation,
                     input_data,
                     left_index=True,
                     right_index=True,
                     how="left")
    # merge = merge.T
    merge.to_csv(output_dir + file_name + time_name + '.csv', index = True)

## 合并数据
if __name__ == "__main__":

    time_name = time.strftime('%Y-%m-%d', time.localtime(time.time()))
    anatation = r"/Users/ranpeng/Desktop/CCA/Data/Fig3/TF_active/network data/merge_index.xlsx"
    input_data = r"/Users/ranpeng/Desktop/CCA/Data/Fig3/TF_active/results/RNA_TN_log2/Tumor/diff_gene_results_pvalue.csv"
    output_dir = r"/Users/ranpeng/Desktop/CCA/Data/Fig3/TF_active/network data/"
    file_name = "TF_TG_pvalue_"
    mapping_data(anatation, input_data, output_dir, file_name)