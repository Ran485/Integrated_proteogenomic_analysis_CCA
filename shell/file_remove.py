#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : file_remove.py
@Time     : 2022/03/19 19:51:48
@Author   : RanPeng 
@Version  : python 3.8.6
@Contact  : 2502388440@hotmail.com
@License  : (C)Copyright 2021-2022, DingLab-CHINA-SHANGHAI
@Function : None
'''
# here put the import lib
import os


# remove the file in the directory if the file size is lower than the given size
def remove_file(dir,format, size):
    '''
    dir: the directory you want to get the filename
    size: the size you want to remove the file
    example:
    remove_file(dir='/Users/ranpeng/Desktop/', size=10)
    :return:
    '''
    for root, dirs, files in os.walk(dir):
        for file in files:
            if file.endswith(format):
                if os.path.getsize(os.path.join(root, file)) < size * 1024 * 1024 * 1024:
                    os.remove(os.path.join(root, file))
                    print("remove file: {}\n".format(os.path.join(root, file)))


if __name__ == '__main__':
    import sys
    dir = sys.argv[1]
    format = sys.argv[2]
    size = int(sys.argv[3])
    remove_file(dir,format, size)