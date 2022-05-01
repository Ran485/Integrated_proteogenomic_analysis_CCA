#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : file_rename.py
@Time     : 2022/03/18 21:08:46
@Author   : RanPeng 
@Version  : python 3.8.6
@Contact  : 2502388440@hotmail.com
@License  : (C)Copyright 2021-2022, DingLab-CHINA-SHANGHAI
@Function : None
'''
# here put the import lib
import os
import sys
import datetime

# rename file according to the given excel file
# Usage: python3 file_rename.py [input_folder] [rename_configfile]
def save_log(path='./'):
    '''
    path， it is a path for save your log about fuction print
    example:
    use  make_print_to_file()   and the   all the information of funtion print , will be write in to a log file
    :return:
    '''
    class Logger(object):
        def __init__(self, filename="Default.log", path="./"):
            self.terminal = sys.stdout
            self.log = open(
                os.path.join(path, filename),
                "a",
                encoding='utf8',
            )

        def write(self, message):
            self.terminal.write(message)
            self.log.write(message)

        def flush(self):
            pass

    fileName = datetime.datetime.now().strftime('day_' + '%Y_%m_%d_%H') + '.log'
    file_modify_time = datetime.datetime.now().strftime(' Date:' + '%Y/%m/%d; Time:%H:%M:%S') + '; Run by: ' + os.getlogin() + ' '
    sys.stdout = Logger(fileName + '.log', path=path)

    # 这里输出之后的所有的输出的print 内容即将写入日志
    print(file_modify_time.center(80, '*'))



def rename_file(input_folder, rename_configfile, reverse=False):
    with open(rename_configfile, 'r', encoding="utf-8") as f:
        lines = f.readlines()[1:]
        for line in lines:
            line = line.strip()
            if line:
                if reverse:
                    original_name = line.split('\t')[1]
                    new_name = line.split('\t')[0]
                else:
                    original_name = line.split('\t')[0]
                    new_name = line.split('\t')[1]
                if os.path.isfile(input_folder + original_name):
                    if os.path.isfile(input_folder + new_name):
                        print(f"{new_name} file: already exists!\n")
                    else:
                        os.rename(input_folder + original_name, input_folder + new_name)
                        print(f"{original_name} was successfully renamed to: {new_name}\n")
                else:
                    print(f"{original_name} file: does not exist!\n")


if __name__ == '__main__':
    import sys
    input_folder = sys.argv[1]
    file = sys.argv[2]
    reverse = sys.argv[3]
    save_log()
    if reverse.lower()=='true' or reverse.lower()=='t':
        rename_file(input_folder, file, reverse=True)
    else:
        rename_file(input_folder, file)