'''
Author: your name
Date: 2022-02-24 23:41:56
LastEditTime: 2022-02-24 23:52:30
LastEditors: Please set LastEditors
Description: 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
FilePath: /罗俊一/script/XGBoost/XGB_script/log_save.py
'''
import sys
import os
import config_file as cfg_file
import sys
import datetime

def save_print_to_file(path='./'):
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

    fileName = datetime.datetime.now().strftime('day_' + '%Y_%m_%d')
    sys.stdout = Logger(fileName + '.log', path=path)

    # 这里输出之后的所有的输出的print 内容即将写入日志
    print(fileName.center(60, '*'))


# if __name__ == '__main__':

#     save_print_to_file(path='/Users/ranpeng/Desktop/Desktop/项目文件/对外服务/罗俊一/script/results/XGBoost/clinical_output_20220224/')
#     # 这里输出之后的所有的输出的print 内容即将写入日志
#     print("1234124")
#     print("--")
#     print(":;;;")
#     print("")
#     print("阿斯顿发11111111111111111")
#     print("zzzzz")
