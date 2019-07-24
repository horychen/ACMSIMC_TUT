# coding:u8

import os
import sys
import femm
from time import time
# from itertools import izip
# from numpy import exp, pi # sqrt
# from numpy import savetxt, c_

def get_slipfreq_torque():
    # call this after mi_analyze
    femm.mi_loadsolution()

    # Physical Amount on the Rotor
    femm.mo_groupselectblock(100) # rotor iron
    femm.mo_groupselectblock(101) # rotor bars
    # Fx = femm.mo_blockintegral(18) #-- 18 x (or r) part of steady-state weighted stress tensor force
    # Fy = femm.mo_blockintegral(19) #--19 y (or z) part of steady-state weighted stress tensor force
    torque = femm.mo_blockintegral(22) #-- 22 = Steady-state weighted stress tensor torque
    freq = femm.mo_getprobleminfo()[1]
    femm.mo_clearblock()
    femm.mo_close()

    return freq, torque

id_solver = int(sys.argv[1])
dir_femm_temp = sys.argv[2][1:-1]
print('ParaSolve', id_solver, end=' ')


femm.openfemm(True) # bHide
# this is essential to reduce elements counts from >50000 to ~20000.
femm.callfemm_noeval('smartmesh(0)')
print('smartmesh is off')

tic = time()
fem_file_path = dir_femm_temp + 'femm_temp_%d.fem'%(id_solver)

# # debug
# print fem_file_path
# # print dir_femm_temp
# os.system('pause')
# quit()


femm.opendocument(fem_file_path)
counter_loop = 0
while True:
    counter_loop += 1
    try:
        femm.mi_analyze(1) # None for inherited. 1 for a minimized window,
        freq, torque = get_slipfreq_torque()
    except Exception as error:
        print('Error occurs when analyzing the femm file.')
        print(error.args)
        if counter_loop > 3:
            raise error
    else:
        break
femm.mi_close()
toc = time()
print('id_solver=%d:'% (id_solver), toc - tic, 's')
femm.closefemm()

# removing files is left for manager to decide
    # os.remove(fem_file_path)
    # os.remove(fem_file_path[:-4]+'.ans')

with open(dir_femm_temp + "femm_temp_%d.txt"%(id_solver), 'w') as handle_torque:
    # write results to a data file (write to partial files to avoid compete between parallel instances)
    handle_torque.write("%g\n%g\n" % ( freq, torque))
