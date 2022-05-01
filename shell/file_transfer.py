# -*- encoding: utf-8 -*-
'''
@File     : file_transfer.py
@Time     : 2022/03/02 20:11:12
@Author   : RanPeng 
@Version  : python 3.8.6
@Contact  : 2502388440@hotmail.com
@License  : (C)Copyright 2021-2022, DingLab-CHINA-SHANGHAI
@Function : None
'''
# here put the import lib

import os
import shutil
import time
from functools import wraps

# Transfer all files ending with _R1 or _R2 in the given directory to the specified folder
# Usage: python3 gene_fusion_preprocess.py [input_folder] [output_folder]
def find_all_files(input_folder):
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.endswith('.fq.gz') or file.endswith('.fastq.gz'):
                fullname = os.path.join(root, file)
                yield fullname

# creat a timer decorator
def timer(func):
    """Return the function start and finished cost time """
    @wraps(func)
    def wrapper(*args,**kwargs):
        print(f"\n{'*'*15} {func.__name__!r} function begins running... {'*'*15}", "\n")
        start_time = time.time()
        result = func(*args,**kwargs)
        print(f"{'*'*10} {func.__name__!r} function was succesfully run in {round((time.time()-start_time), 2)}s {'*'*10}\n")
        return result
    return wrapper

@timer
def file_transfer(output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print("Successfully created output directory: {}".format(output_folder))
    else:
        print("Output directory already exists: {}".format(output_folder))

    for file in find_all_files(input_folder):
        if not (os.path.isfile(output_folder + file)):
            # shutil.move(file, output_folder) # move file to output_folder
            shutil.copyfile(file, output_folder + file.split('/')[-1])
            # shutil.copy(file, output_folder) # copy the file to the output folder,spped slowly
            print("Successfully copied file: {}".format(file))
        else:
            print("File already exists: {}".format(output_folder + file))


if __name__ == '__main__':
    import sys
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print("Successfully created output directory: {}".format(output_folder)) 
    print(f"Input folder is: {input_folder}, \nOutput folder is: {output_folder}")
    file_transfer(output_folder)