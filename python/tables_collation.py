import pandas as pd
import numpy as np

def save_multisheet_excel(append_file, input_file, filetype = 'excel', sheet_name="这是第一个sheet"):
    """
    read table and write it to the multisheet excel file
    """
    
    # excel file use 'pd.read_excel', csv/txt file use 'pd.read_csv'
    if filetype == 'excel':
        df = pd.read_excel(input_file)
    elif filetype == 'csv':
        df = pd.read_csv(input_file)
    else:
        df = pd.read_csv(input_file, sep="\t")

    # 对 df 的数据进行相应的处理, 一般操作的话将其注释掉
    # df.replace(np.nan, 0, inplace=True)

    # return df
    with pd.ExcelWriter(append_file, mode="a", engine="openpyxl") as writer:
        df.to_excel(writer, index = False, sheet_name = sheet_name)


if __name__ == '__main__':
       input_file = r"/Users/ranpeng/Desktop/CCA/CCA_mapping_data/phosphoprotein_abudance_FOT(2021-2-3).csv"
       append_file = r"/Users/ranpeng/Desktop/CCA/Supplementary_Tables/Table1.xlsx"
       save_multisheet_excel(append_file, input_file, filetype = 'csv', sheet_name="phosphoprotein_abudance")
