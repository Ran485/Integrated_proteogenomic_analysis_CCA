#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File     : qsub.py
@Time     : 2022/03/03 16:45:25
@Author   : RanPeng 
@Version  : python 3.8.6
@Contact  : 2502388440@hotmail.com
@License  : (C)Copyright 2021-2022, DingLab-CHINA-SHANGHAI
@Function : None
'''
# here put the import lib
import csv
import pwd
import subprocess
import time
import sys
import os
import datetime


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

    fileName = datetime.datetime.now().strftime('job_samples_submission_' + '%Y_%m_%d_%H_%M')
    sys.stdout = Logger(fileName + '.log', path=path)

    # 这里输出之后的所有的输出的print 内容即将写入日志
    print(fileName.center(90, '*'))


# param_file = './job_params.csv'
# cwd = '/home/bioinfo/RNA-seq'
def qsub(param_file, cwd):
    """Batch submission of transcriptome data processing jobs

    Args:
        param_file (csv file): param_file = './job_params.csv'
        cwd (dir): Current working directory.
    """
    with open(param_file, 'r') as f:
        reader_list = csv.reader(f, delimiter=' ', quotechar='|')
        for sample in reader_list:
            sample = "".join(sample)  # 去掉列表内的中括号
            qsub_cmd = 'qsub -N {0} -d {1} -v sample={0} arriba_fusion_template_cluster.sh'.format(
                sample, cwd)
            exit_status = subprocess.call(qsub_cmd, shell=True)
            if exit_status is 1:
                print('Job "{}" failed to submit'.format(qsub_cmd))
            else:
                print('Job "{}" successfully submitted\n'.format(qsub_cmd))
            time.sleep(1)

    print("Done submitting all jobs!")


if __name__ == '__main__':
    import sys
    param_file = sys.argv[1]
    cwd = sys.argv[2]
    # save Server-assigned nodes and corresponding sample numbers information
    save_log(path=cwd) # save log
    print(f"param_file is: {param_file}, \nCurrent working directory is: {os.getcwd()}\n")
    qsub(param_file, cwd)
