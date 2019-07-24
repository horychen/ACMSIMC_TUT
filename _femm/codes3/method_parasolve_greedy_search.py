# coding:u8
# method_parasolve_greedy_search

import sys
import subprocess

dir_femm_temp = sys.argv[1]
number_of_instantces = int(sys.argv[2])
str_stack_length = sys.argv[3]

# debug
# print sys.argv
# import os
# os.system("pause")
# quit()

# manager
proc = subprocess.Popen([sys.executable, 'parasolve_greedy_search_manager.py', 
                         str(number_of_instantces), dir_femm_temp, str_stack_length], bufsize=-1)

# 等了也是白等，直接就运行返回了
proc.wait() 



