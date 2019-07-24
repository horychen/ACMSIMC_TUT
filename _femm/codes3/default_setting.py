import os
fea_config_dict = {
    ##########################
    # Sysetm Control
    ##########################
    'designer.Show': True,
    'OnlyTableResults':False, # modified later according to pc_name # the only reason we want save mesh results are to obtain voltage profile for power factor evaluation
    'local_sensitivity_analysis':False,
    'local_sensitivity_analysis_number_of_variants': 20,
    'flag_optimization':True, # also use true for sensitivity analysis
    'Restart':False, # restart from frequency analysis is not needed, because SSATA is checked and JMAG 17103l version is used.
    'MultipleCPUs':True,
        # True for:
            # multiple cpu (SMP=2)
            # use directSolver over ICCG Solver

    ##########################
    # FEA Setting
    ##########################
    'TranRef-StepPerCycle':40, # FEMM: 5 deg   # 360 to be precise as FEMM: 0.5 deg 
    # 'FrequencyRange':range(1,6), # the first generation for PSO
    'number_of_steps_1stTTS':32, # TSS actually...
    'number_of_steps_2ndTTS':32, # use a multiples of 4! # 8*32 # steps for half period (0.5). That is, we implement two time sections, the 1st section lasts half slip period and the 2nd section lasts half fandamental period.
    'number_cycles_prolonged':1, # 150
    'FEMM_Coarse_Mesh':True,

    ##########################
    # Optimization
    ##########################
    'JMAG_Scheduler':False, # multi-cores run
    'delete_results_after_calculation': False, # check if True can we still export Terminal Voltage? 如果是True，那么就得不到Terminal Voltage和功率因数了！
    'use_weights':None,

    ##########################
    # Design Specifications
    ##########################
    'Active_Qr':None, #16, #32,
    'PoleSpecific': True,
    'use_drop_shape_rotor_bar': True,
    'DPNV': True,
    # 'DPNV_separate_winding_implementation': None, # this is obsolete feature (it is not good because the copper loss is different from the reality and only works for Qs=24, p=2 case)
    'mimic_separate_winding_with_DPNV_winding':False,
    'End_Ring_Resistance':0, # 0 for consistency with FEMM with pre-determined currents # 9.69e-6, # this is still too small for Chiba's winding
    'Steel': 'M19Gauge29', #'M15','Arnon5', 
                           # 75 deg Celsus: If you modify the temperature here, you should update the initial design (the im.DriveW_Rs should be updated and it used in JMAG FEM Coil)
    # This is actually the resistivity!
    # This is actually the resistivity!
    # This is actually the resistivity!
    'Bar_Conductivity':1/((3.76*75+873)*1e-9/55.), # @75 Degree Celsius # 1/((3.76*100+873)*1e-9/55.) for Copper, where temperature is 25 or 100 deg Celsius.
    # 'Bar_Conductivity':40e6, # 40e6 for Aluminium; 
}
def where_am_i(fea_config_dict):
    def get_pc_name():
        import platform
        import socket
        n1 = platform.node()
        n2 = socket.gethostname()
        n3 = os.environ["COMPUTERNAME"]
        if n1 == n2 == n3:
            return n1
        elif n1 == n2:
            return n1
        elif n1 == n3:
            return n1
        elif n2 == n3:
            return n2
        else:
            raise Exception("Computer names are not equal to each other.")

    dir_interpreter = os.path.abspath('')
    print(dir_interpreter)
    if 'codes3' in dir_interpreter:
        dir_parent = dir_interpreter + '/../' # '../' <- relative path is not always a good option
    else:
        if '_femm' in dir_interpreter:
            dir_parent = dir_interpreter + '/'
        else:
            raise Exception('Do not run the script from this directory: %s'%(dir_interpreter))
    dir_codes = dir_parent + 'codes3/'
    dir_femm_files = dir_parent + 'femm_files/'
    pc_name = get_pc_name()
    os.chdir(dir_codes)
    # print(dir_parent, dir_codes, pc_name, sep='\n'); quit()
    fea_config_dict['dir_interpreter']       = dir_interpreter
    fea_config_dict['dir_parent']            = dir_parent
    fea_config_dict['dir_codes']             = dir_codes
    fea_config_dict['dir_femm_files']        = dir_femm_files
    fea_config_dict['pc_name']               = pc_name

    fea_config_dict['delete_results_after_calculation'] = False # True for saving disk space (but you lose voltage profile and element data)
    if fea_config_dict['Restart'] == False:
        fea_config_dict['OnlyTableResults'] = True  # save disk space for my PC
    # However, we need field data for iron loss calculation
    fea_config_dict['OnlyTableResults'] = False 

where_am_i(fea_config_dict) # add path to fea_config_dict

from sys import path as sys_path
# for importing your package
sys_path.append(fea_config_dict['dir_codes'])

import utility
import pyrhonen_procedure_as_function
from pylab import np, plt
import logging

# run_list = [1,1,0,0,0] # use JMAG only
run_list = [0,1,0,0,0] # use FEMM to search for breakdown slip
fea_config_dict['jmag_run_list'] = run_list
def build_model_name_prefix(fea_config_dict, UID=None):
    if fea_config_dict['flag_optimization'] == True:
        fea_config_dict['model_name_prefix'] = 'OP'
    else:
        fea_config_dict['model_name_prefix'] = ''

    if fea_config_dict['PoleSpecific'] == True:
        fea_config_dict['model_name_prefix'] += '_PS'
    else:
        fea_config_dict['model_name_prefix'] += '_SC'

    fea_config_dict['model_name_prefix'] += '_Qr%d'%(fea_config_dict['Active_Qr'])

    if fea_config_dict['End_Ring_Resistance'] == 0:
        fea_config_dict['model_name_prefix'] += '_NoEndRing'

    if UID is not None:
        fea_config_dict['model_name_prefix'] += 'UID_%4d'%(UID)

    return fea_config_dict['model_name_prefix']
# print(build_model_name_prefix(fea_config_dict))

# FEMM Static FEA
fea_config_dict['femm_deg_per_step'] = None # 0.25 * (360/4) / utility.lcm(24/4., fea_config_dict['Active_Qr']/4.) # at least half period
# fea_config_dict['femm_deg_per_step'] = 1 * (360/4) / utility.lcm(24/4., fea_config_dict['Active_Qr']/4.) # at least half period
# fea_config_dict['femm_deg_per_step'] = 0.1 #0.5 # deg
# print 'femm_deg_per_step is', fea_config_dict['femm_deg_per_step'], 'deg (Qs=24, p=2)'


