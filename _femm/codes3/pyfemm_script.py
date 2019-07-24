# coding:utf-8
#execfile('D:/OneDrive - UW-Madison/c/codes/pyfemm_script.py')
#execfile(r'K:\jchen782\JMAG\c\codes/pyfemm_script.py')

''' 1. General Information & Packages Loading
'''
# execfile('D:/OneDrive - UW-Madison/c/codes/default_setting.py') # Absolute path is needed for running in JMAG
# exec(compile(open('./default_setting.py').read(), './default_setting.py', 'exec')) # Relative path is enough if run outside JMAG
filename = './default_setting.py'
exec(compile(open(filename, "rb").read(), filename, 'exec'), globals(), locals())


# Debug
# if os.path.exists('d:/femm42/PS_Qr32_NoEndRing_M19Gauge29_DPNV_1e3Hz'):
#     os.system('bash -c "rm -r /mnt/d/femm42/PS_Qr32_NoEndRing_M19Gauge29_DPNV_1e3Hz"')
# if os.path.exists('d:/OneDrive - UW-Madison/c/pop/Tran2TSS_PS_Opti.txt'):
#     os.system('bash -c "mv /mnt/d/OneDrive\ -\ UW-Madison/c/pop/Tran2TSS_PS_Opti.txt /mnt/d/OneDrive\ -\ UW-Madison/c/pop/initial_design.txt"')

logger = utility.myLogger(fea_config_dict['dir_codes'], prefix='iemdc_')

fea_config_dict['flag_optimization'] = True # we need to generate and exploit sw.init_pop
fea_config_dict['Active_Qr'] = 36

run_list = [1,1,1,1,0] # Static FEA with JMAG is too slow, so use FEMM to do that part
# run_list = [1,1,1,0,0] # TranRef is replaced by Tran2TSSProlong

run_folder = r'run#99/' # Test run for single one design
run_folder = r'run#98/' # Test run for the loop
run_folder = r'run#97/' # Test run for the loop
run_folder = r'run#96/' # Test run for the loop
run_folder = r'run#95/' # Test run for the loop
run_folder = r'run#94/' # Test run for the loop

run_folder = r'run#93/' # get TranRef back and build rotor current plot

if 'Severson' in fea_config_dict['pc_name']:
    # check number_cycles_prolonged is applied in number_of_total_steps (o)
    fea_config_dict['number_cycles_prolonged'] = 0 #150 # 1
    # app.Show() (o)
    fea_config_dict['designer.Show'] = True

    # serverson01    
    if '01' in fea_config_dict['pc_name']:
        # run_folder = r'run#299/' # some bug

        fea_config_dict['number_cycles_prolonged'] = 150
        run_folder = r'run#298/' # 2 deg, 24 steps, 32 steps

        # fea_config_dict['number_cycles_prolonged'] = 0
        # fea_config_dict['number_of_steps_2ndTTS'] = 400
        # run_folder = r'run#297/' # sensitivity of solving steps (fine step case): 0.5 deg, 48 steps, 400 steps (error occurs: refarray is 3 not 5)
        # run_folder = r'run#296/' # sensitivity of solving steps (fine step case): 0.5 deg, 48 steps, 400 steps
        # run_folder = r'run#295/'

    # serverson02 
    elif '02' in fea_config_dict['pc_name']:
        # run_folder = r'run#399/' # some bug
        # run_folder = r'run#398/' # 2 deg, 24 steps, 32 steps

        fea_config_dict['number_of_steps_2ndTTS'] = 16
        run_folder = r'run#397/' # This is not used
        run_folder = r'run#395/' # sensitivity of solving steps (coarse step case): 5 deg, 12 steps, 16 steps
    else:
        raise Exception('Where are you?')

fea_config_dict['run_folder'] = run_folder
fea_config_dict['jmag_run_list'] = run_list
fea_config_dict['DPNV_separate_winding_implementation'] = True

# rebuild the name
build_model_name_prefix(fea_config_dict)

# fea_config_dict['femm_deg_per_step'] = 0.25 * (360/4) / utility.lcm(24/4., fea_config_dict['Active_Qr']/4.) # at least half period
# fea_config_dict['femm_deg_per_step'] = 1 * (360/4) / utility.lcm(24/4., fea_config_dict['Active_Qr']/4.) # at least half period
fea_config_dict['femm_deg_per_step'] = 2 #0.5 #0.1 # deg # 5 deg will cause error when DFT getting iron loss results
logger.info('femm_deg_per_step is %g deg (Qs=24, p=2).'%(fea_config_dict['femm_deg_per_step']))

fea_config_dict['ec_rotate_divisions_per_slot_pitch'] = 24
logger.info('ec_rotate_divisions_per_slot_pitch is %g deg (Qs=24, p=2).'%(fea_config_dict['ec_rotate_divisions_per_slot_pitch']))

''' 2. Initilize Swarm and Initial Pyrhonen's Design (Run this part in JMAG)
''' # 1e-1也还是太小了（第三次报错），至少0.5mm长吧 # 1e-1 is the least geometry value. a 1e-2 will leads to：转子闭口槽极限，会导致edge过小，从而报错：small arc entity exists.png
de_config_dict = {  'bounds':     [[3,9], [0.5,4], [5e-1,3], [2.5, 6], [5e-1,3], [1,10], [5e-1,3]], 
                    'mut':        0.8,
                    'crossp':     0.7,
                    'popsize':    50, # 100
                    'iterations': 20 } # begin at 5
# get initial design as im
sw = population.swarm(fea_config_dict, de_config_dict=de_config_dict)
# sw.show(which='all')
# print sw.im.show()


# generate the initial generation
sw.generate_pop()

im_initial = sw.im
# print im_initial.l21
# print im_initial.l22

''' 3. Initialize FEMM Solver
'''
# define problem
logger.info('Running Script for FEMM with %s'%(run_folder))
solver_jmag = FEMM_Solver.FEMM_Solver(im_initial, flag_read_from_jmag=True, freq=0) # static
solver_femm = FEMM_Solver.FEMM_Solver(im_initial, flag_read_from_jmag=False, freq=2.23) # eddy+static



# # debug
# if not solver_femm.has_results(solver_femm.dir_run + 'sweeping/'):
#     solver_femm.run_frequency_sweeping(range(1,6))
#     data_femm = solver_femm.show_results(bool_plot=False) # this write rotor currents to file, which will be used later for static FEA
# else:
#     solver_femm.show_results_eddycurrent(True) # get right slip
# print solver_femm.get_iron_loss(MAX_FREQUENCY=50e3)
# print solver_femm.get_iron_loss(MAX_FREQUENCY=30e3)
# print solver_femm.get_iron_loss(MAX_FREQUENCY=10e3)
# print solver_femm.get_iron_loss(MAX_FREQUENCY=5e3)
# from pylab import show; show()
# quit()
    # solver_femm.read_current_conditions_from_FEMM()
    # solver_femm.get_copper_loss()
    # quit()

    # solver_femm.keep_collecting_static_results_for_optimization()
    # # # # quit()
    # # # # # load Arnon5 from file
    # # # # sw.print_array()
    # # # # quit()





# 50 Random Design Evaluation for IEMDC 2019
from time import time as clock_time
from time import sleep
import numpy as np
from pylab import show
from collections import OrderedDict

min_b, max_b = np.asarray(sw.de_config_dict['bounds']).T 
diff = np.fabs(min_b - max_b)
pop_denorm = min_b + sw.init_pop * diff

for ind, individual_denorm in enumerate(pop_denorm):

    im_variant = population.bearingless_induction_motor_design.local_design_variant(im_initial, \
                    0, ind, individual_denorm) # due to compatability issues: a new child class is used instead

    if not sw.has_results(im_variant):
        # break
        sw.run(im_variant, individual_index=ind, run_list=run_list)
        # im_variant.csv_previous_solve = sw.dir_csv_output_folder + im_variant.get_individual_name() + u"Freq" + '_circuit_current.csv'

    # FEMM Static Solver with pre-determined rotor currents from JMAG
    tic = clock_time()
    solver_jmag = FEMM_Solver.FEMM_Solver(im_variant, individual_index=ind, flag_read_from_jmag=True, freq=0, bool_static_fea_loss=False) # static
    if not solver_jmag.has_results():
        print('run_rotating_static_FEA')
        # utility.blockPrint()
        solver_jmag.run_rotating_static_FEA()
        solver_jmag.parallel_solve()
        # utility.enablePrint()

    # collecting parasolve with post-process
    # wait for .ans files
    # data_solver_jmag = solver_jmag.show_results_static(bool_plot=False) # this will wait as well?
    while not solver_jmag.has_results():
        print(clock_time())
        sleep(3)
    results_dict = {}
    for f in [f for f in os.listdir(solver_jmag.dir_run) if 'static_results' in f]:
        data = np.loadtxt(solver_jmag.dir_run + f, unpack=True, usecols=(0,1,2,3))
        for i in range(len(data[0])):
            results_dict[data[0][i]] = (data[1][i], data[2][i], data[3][i]) 
    keys_without_duplicates = [key for key, item in results_dict.items()]
    keys_without_duplicates.sort()
    with open(solver_jmag.dir_run + "no_duplicates.txt", 'w') as fw:
        for key in keys_without_duplicates:
            fw.writelines('%g %g %g %g\n' % (key, results_dict[key][0], results_dict[key][1], results_dict[key][2]))
    data_solver_jmag = np.array([ keys_without_duplicates, 
                                     [results_dict[key][0] for key in keys_without_duplicates], 
                                     [results_dict[key][1] for key in keys_without_duplicates], 
                                     [results_dict[key][2] for key in keys_without_duplicates]])
    # print data_solver_jmag
    toc = clock_time()
    print(ind, 'tic: %g. toc: %g. diff:%g' % (tic, toc, toc-tic))



    # for iemdc 2019 Qr36rotor_current 
    # sw.show_results(femm_solver_data=data_solver_jmag)
    # data = solver_jmag.show_results(bool_plot=False)
    # sw.show_results(femm_solver_data=data)
    sw.show_results_iemdc19(im_variant, femm_solver_data=data_solver_jmag, femm_rotor_current_function=solver_jmag.get_rotor_current_function())
    # sw.timeStepSensitivity()
    from pylab import show
    show()
    quit()


    # JMAG results (EC-Rotate and Tran2TSS and Tran2TSSProlongRef)
    data_results = utility.collect_jmag_Tran2TSSProlong_results(im_variant, sw.dir_csv_output_folder, sw.fea_config_dict, sw.axeses, femm_solver_data=data_solver_jmag)
    # show()
    # quit()
    sw.fig_main.savefig(sw.dir_run + im_variant.get_individual_name() + 'results.png', dpi=150)
    utility.pyplot_clear(sw.axeses)

    # write to file for inspection
    with open(sw.dir_run + 'iemdc_data.txt', 'a') as f:
        f.write(','.join(['%g'%(el) for el in [ind] + data_results]) + '\n')

    # break

# 绘制K线图表征最大误差和最小误差
exec(compile(open('./iemdc_random_design_plot.py').read(), './iemdc_random_design_plot.py', 'exec'))
quit()












''' 4. Show results, if not exist then produce it
'''
from numpy import arange
if True: # this generate plots for iemdc19
    data = solver_jmag.show_results(bool_plot=False)
    # sw.show_results(femm_solver_data=data)
    sw.show_results_iemdc19(femm_solver_data=data, femm_rotor_current_function=solver_jmag.get_rotor_current_function())
    # sw.timeStepSensitivity()
else:
    # FEMM only
    # FEMM Static Solver with pre-determined rotor currents from FEMM
    if not solver_femm.has_results():
        solver_femm.run_frequency_sweeping(arange(0.5,5.1,0.5))
        data_femm = solver_femm.show_results(bool_plot=False) # this write rotor currents to file, which will be used later for static FEA

        quit()
    #     solver_femm.run_rotating_static_FEA()
    #     solver_femm.parallel_solve()


    # JMAG 
    if not sw.has_results(im_initial, 'Freq'):
        sw.run(im_initial, run_list=sw.fea_config_dict['jmag_run_list'])
        if not sw.has_results(im_initial, 'Freq'):
            raise Exception('Something went south with JMAG.')

    print('Rotor current solved by FEMM is shown here:')
    solver_femm.read_current_conditions_from_FEMM()

    print('Rotor current solved by JMAG is shown here:')
    solver_jmag.read_current_from_EC_FEA()
    for key, item in solver_jmag.dict_rotor_current_from_EC_FEA.items():
        if '1' in key:
            from math import sqrt
            print(key, sqrt(item[0]**2+item[1]**2))

    # FEMM Static Solver with pre-determined rotor currents from JMAG
    if not solver_jmag.has_results():
        solver_jmag.run_rotating_static_FEA()
        solver_jmag.parallel_solve()

    # Compare JMAG and FEMM's rotor currents to see DPNV is implemented correctly in bot JMAG and FEMM    
    data_femm = solver_femm.show_results_static(bool_plot=False)
    data_jmag = solver_jmag.show_results_static(bool_plot=False)

    if data_femm is None or data_jmag is None:
        if data_femm is None:
            print('data_femm is', data_femm, 'try run again')
        if data_jmag is None:
            print('data_jmag is', data_jmag, 'try run again')
    else:
        from pylab import subplots
        fig, axes = subplots(3, 1, sharex=True)
        for data in [data_femm, data_jmag]:
            ax = axes[0]; ax.plot(data[0]*0.1, data[1], alpha=0.5, label='torque'); ax.legend(); ax.grid()
            ax = axes[1]; ax.plot(data[0]*0.1, data[2], alpha=0.5, label='Fx'); ax.legend(); ax.grid()
            ax = axes[2]; ax.plot(data[0]*0.1, data[3], alpha=0.5, label='Fy'); ax.legend(); ax.grid()



# sw.write_to_file_fea_config_dict()
from pylab import show; show()





# Run JCF from command linie instead of GUI
# if not sw.has_results(im_initial, study_type='Tran2TSS'):
#     os.system(r'set InsDir=D:\Program Files\JMAG-Designer17.1/')
#     os.system(r'set WorkDir=D:\JMAG_Files\JCF/')
#     os.system(r'cd /d "%InsDir%"')
#     for jcf_file in os.listdir(sw.dir_jcf):
#         if 'Tran2TSS' in jcf_file and 'Mesh' in jcf_file:
#             os.system('ExecSolver "%WorkDir%' + jcf_file + '"')
if False: # Test and Debug

    ''' 3. Eddy Current FEA with FEMM
    '''
    solver = FEMM_Solver.FEMM_Solver(deg_per_step, im_initial, dir_codes, dir_femm_files + run_folder)

    solver.read_current_from_EC_FEA()
    for key, item in solver.dict_rotor_current_from_EC_FEA.items():
        if '1' in key:
            from math import sqrt
            print(key, sqrt(item[0]**2+item[1]**2))

    if solver.has_results():
        solver.show_results()
    else:
        solver.run_frequency_sweeping(arange(2, 5.1, 0.25), fraction=2)
        solver.parallel_solve(6, bool_watchdog_postproc='JMAG' not in dir_interpreter)


'''
Two Shared process - Max CPU occupation 32.6%
Freq: 2:16 (13 Steps)
Tran2TSS: 3:38 (81 Steps)
Freq-FFVRC: 1:33 (10 Steps)
TranRef: 1:49:02 (3354 Steps) -> 3:18:19 if no MultiCPU
StaticJMAG: 
StaticFEMM: 15 sec one eddy current solve
            7 sec one static solve
'''

a=''' Loss Study of Transient FEA
# -*- coding: utf-8 -*-
app = designer.GetApplication()
app.GetModel(u"PS_ID32").GetStudy(u"Loss_Tran2TSS-Prolong").GetCondition(u"P1").SetValue(u"BasicFrequencyType", 1)
app.GetModel(u"PS_ID32").GetStudy(u"Loss_Tran2TSS-Prolong").GetCondition(u"P1").SetValue(u"RevolutionSpeed", 15000)
app.GetModel(u"PS_ID32").GetStudy(u"Loss_Tran2TSS-Prolong").GetCondition(u"P1").SetValue(u"StressType", 0)
app.GetModel(u"PS_ID32").GetStudy(u"Loss_Tran2TSS-Prolong").GetCondition(u"P1").ClearParts()
sel = app.GetModel(u"PS_ID32").GetStudy(u"Loss_Tran2TSS-Prolong").GetCondition(u"P1").GetSelection()
sel.SelectPart(35)
app.GetModel(u"PS_ID32").GetStudy(u"Loss_Tran2TSS-Prolong").GetCondition(u"P1").AddSelected(sel)
app.SetCurrentStudy(u"Loss_Tran2TSS-Prolong")
app.GetModel(u"PS_ID32").GetStudy(u"Loss_Tran2TSS-Prolong").CreateCondition(u"Ironloss", u"P2")
app.GetModel(u"PS_ID32").GetStudy(u"Loss_Tran2TSS-Prolong").GetCondition(u"P2").SetValue(u"BasicFrequencyType", 2)
app.GetModel(u"PS_ID32").GetStudy(u"Loss_Tran2TSS-Prolong").GetCondition(u"P2").SetValue(u"BasicFrequency", 3)
app.GetModel(u"PS_ID32").GetStudy(u"Loss_Tran2TSS-Prolong").GetCondition(u"P2").SetValue(u"StressType", 0)
app.GetModel(u"PS_ID32").GetStudy(u"Loss_Tran2TSS-Prolong").GetCondition(u"P2").ClearParts()
sel = app.GetModel(u"PS_ID32").GetStudy(u"Loss_Tran2TSS-Prolong").GetCondition(u"P2").GetSelection()
sel.SelectPart(2)
app.GetModel(u"PS_ID32").GetStudy(u"Loss_Tran2TSS-Prolong").GetCondition(u"P2").AddSelected(sel)
app.SetCurrentStudy(u"Loss_Tran2TSS-Prolong")
'''

