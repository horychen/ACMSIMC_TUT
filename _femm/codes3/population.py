# -*- coding: utf-8 -*-
# execfile(r'D:\Users\horyc\OneDrive - UW-Madison\ec_rotate.py') # , {'__name__': 'load'})

from math import cos, sin, pi
import csv
import logging
import numpy as np  # for de
import os
from pylab import plot, legend, grid, figure, subplots, array, mpl, show
import utility
from VanGogh import VanGogh
import FEMM_Solver

from time import time as clock_time

from pyrhonen_procedure_as_function import get_material_data
import winding_layout

EPS = 1e-2 # unit: mm


class swarm(object):

    def __init__(self, fea_config_dict, de_config_dict=None):
        # directories part I
        self.dir_parent             = fea_config_dict['dir_parent']
        self.initial_design_file    = self.dir_parent + 'pop/' + r'initial_design.txt'

        # load initial design using the obsolete class bearingless_induction_motor_design
        self.im_list = []
        with open(self.initial_design_file, 'r') as f: 
            print(self.initial_design_file)
            for row in self.csv_row_reader(f):
                im = bearingless_induction_motor_design([row[0]]+[float(el) for el in row[1:]], fea_config_dict, model_name_prefix=fea_config_dict['model_name_prefix'])
                self.im_list.append(im)
        for im in self.im_list:
            if im.Qr == fea_config_dict['Active_Qr']:
                self.im = im
                print('Take the first Active_Qr match as initial design')
                break
        try: 
            self.im
        except:
            print('There is no design matching Active_Qr.')
            msg = 'Please activate one initial design. Refer %s.' % (self.initial_design_file)
            logger = logging.getLogger(__name__)
            logger.warn(msg)
            raise Exception('no match for Active_Qr: %d'%(fea_config_dict['Active_Qr']))

        # directories part II
        if im.DriveW_Freq == 1000: # New design for 1000 Hz machine. some patch for my scrappy codes (lot of bugs are fixed, we need a new name).
            im.model_name_prefix += '_1e3Hz'
            self.model_name_prefix = im.model_name_prefix
            fea_config_dict['model_name_prefix'] = im.model_name_prefix
            print('[Updated] New model_name_prefix is', self.model_name_prefix)
        else:
            # design objective
            self.model_name_prefix = fea_config_dict['model_name_prefix']

        # csv output folder
        if fea_config_dict['flag_optimization'] == False:
            self.dir_csv_output_folder  = self.dir_parent + 'csv/' + self.model_name_prefix + '/'
        else:
            self.dir_csv_output_folder  = self.dir_parent + 'csv_opti/' + self.model_name_prefix + '/'
        if not os.path.exists(self.dir_csv_output_folder):
            os.makedirs(self.dir_csv_output_folder)

        if fea_config_dict['local_sensitivity_analysis'] == True:
            self.run_folder         = fea_config_dict['run_folder'][:-1] + 'lsa/'
        else:
            self.run_folder         = fea_config_dict['run_folder']

        self.dir_run                = self.dir_parent + 'pop/' + self.run_folder

        if fea_config_dict['flag_optimization'] == True:
            self.dir_project_files  = fea_config_dict['dir_project_files'] + self.run_folder
        else:
            self.dir_project_files  = fea_config_dict['dir_project_files']

        self.dir_jcf                = self.dir_project_files + 'jcf/'
        self.pc_name                = fea_config_dict['pc_name']
        self.fea_config_dict        = fea_config_dict

        # dict of optimization
        self.de_config_dict = de_config_dict

        # swarm state control
        self.bool_first_time_call_de = True  # no use
        self.bool_auto_recovered_run = False # auto-recovered run

        # post-process feature
        self.fig_main, self.axeses = subplots(2, 2, sharex=True, dpi=150, figsize=(16, 8), facecolor='w', edgecolor='k')
        utility.pyplot_clear(self.axeses)

    def init_femm_solver(self):
        # and let jmag know about it
        self.femm_solver = FEMM_Solver.FEMM_Solver(self.im, flag_read_from_jmag=False, freq=2.233) # eddy+static

    def write_to_file_fea_config_dict(self):
        with open(self.dir_run + '../FEA_CONFIG-%s.txt'%(self.fea_config_dict['model_name_prefix']), 'w') as f:
            for key, val in self.fea_config_dict.items():
                # print key, val
                f.write('%s:%s\n' % (key, str(val)) )

    def generate_pop(self, specified_initial_design_denorm=None):
        # csv_output folder is used for optimziation
        self.dir_csv_output_folder  = self.fea_config_dict['dir_parent'] + 'csv_opti/' + self.run_folder

        # check if it is a new run
        if not os.path.exists(self.dir_run):
            logger = logging.getLogger(__name__)
            logger.debug('There is no run yet. Generate the run folder for pop as %s...', self.run_folder)
            os.makedirs(self.dir_run)
        if not os.path.exists(self.dir_csv_output_folder):
            try:  
                os.makedirs(self.dir_csv_output_folder)
            except OSError as e:
                logging.getLogger(__name__).error("Creation of the directory %s failed" % self.dir_csv_output_folder, exc_info=True)
                raise e
            else:
                logging.getLogger(__name__).info("Successfully created the directory %s " % self.dir_csv_output_folder)

        # the run folder has been created. check the pop data files
        self.index_interrupt_beginning = 0 # 不出意外，就从第一个个体开始迭代。但是很多时候跑到第150个个体的时候跑断了，你总不会想要把这一代都删了，重新跑吧？
        for file in os.listdir(self.dir_run):
            if 'ongoing' in file:
                # remove gen and fit files
                # os.remove(self.dir_run + file) 
                if 'gen#' in file:
                    print(file)
                    self.ongoing_pop_denorm = []
                    with open(self.dir_run + file, 'r') as f:
                        for row in self.csv_row_reader(f):
                            if len(row)>0: # there could be empty row, since we use append mode and write down f.write('\n')
                                self.ongoing_pop_denorm.append([float(el) for el in row])

                    self.ongoing_pop_denorm = np.asarray(self.ongoing_pop_denorm)

                    logger = logging.getLogger(__name__)
                    logger.warn('Unfinished iteration is found with ongoing files in run folder.') # Make sure the project is not opened in JMAG, and we are going to remove the project files and ongoing files.')
                    # logger.warn(u'不出意外，就从第一个个体开始迭代。但是很多时候跑到第150个个体的时候跑断了，你总不会想要把这一代都删了，重新跑吧？')
                    # os.remove(u"D:/JMAG_Files/" + run_folder[:-1] + file[:-12] + ".jproj")

                    # logger.debug('\n'.join(','.join('%.16f'%(x) for x in y) for y in self.ongoing_pop_denorm))
                    print('List ongoing_pop_denorm here:', self.ongoing_pop_denorm.tolist())

                if 'fit#' in file:
                    print(file)
                    self.ongoing_fitness = []
                    with open(self.dir_run + file, 'r') as f: 
                        for row in self.csv_row_reader(f):
                            if len(row)>0: # there could be empty row, since we use append mode and write down f.write('\n')
                                self.ongoing_fitness.append(float(row[0])) # not neccessary to be an array. a list is enough
                    print('List ongoing_fitness here:', self.ongoing_fitness)
                    # ongoing_pop (yes) and ongoing_fitness will be directly used in de.

                if 'liv#' in file:
                    print(file)
                    self.ongoing_living_pop_denorm = []
                    with open(self.dir_run + file, 'r') as f: 
                        for row in self.csv_row_reader(f):
                            if len(row)>0: 
                                self.ongoing_living_pop_denorm.append([float(el) for el in row])
                    print('List ongoing_living_pop_denorm here:', self.ongoing_living_pop_denorm)

                    file = file[:3] + '_fit' + file[3:8] + '.txt' # liv_fit#xxxx.txt has always no tag '-ongoing'.
                    print(file)
                    self.ongoing_living_fitness = []
                    with open(self.dir_run + file, 'r') as f: 
                        for row in self.csv_row_reader(f):
                            if len(row)>0: 
                                self.ongoing_living_fitness.append(float(row[0]))
                    print('List ongoing_living_fitness here:', self.ongoing_living_fitness)
            
                    # decide at which j the fobj will called in de()
                    self.index_interrupt_beginning = len(self.ongoing_living_pop_denorm)
                    if len(self.ongoing_living_pop_denorm) != len(self.ongoing_pop_denorm):
                        raise Exception('It seemed that Swarm failed to write living pop after writing the generation pop. Manually delete the extra lines in your gen# and fit# files to be consistent with liv# file.')

        # search for completed generation files
        generations_complete = [file[4:8] for file in os.listdir(self.dir_run) if 'gen' in file and not 'ongoing' in file]
        generations_complete_and_onging = [file[4:8] for file in os.listdir(self.dir_run) if 'gen' in file]

        # initialize for de-normalization
        popsize = self.de_config_dict['popsize']
        # the least popsize is 4
        if popsize<=3:
            logger = logging.getLogger(__name__)
            logger.error('The popsize must be greater than 3 so the choice function can pick up three other individuals different from the individual under mutation among the population')
            raise Exception('Specify a popsize larger than 3.')
        bounds = self.de_config_dict['bounds']
        dimensions = len(bounds)
        min_b, max_b = np.asarray(bounds).T 
        diff = np.fabs(min_b - max_b)
        if self.index_interrupt_beginning != 0:
                  # pop_denorm = min_b + pop * diff =>
            self.ongoing_pop = (self.ongoing_pop_denorm - min_b) / diff


        # check for number of generations_complete
        if len(generations_complete) == 0:
            logger = logging.getLogger(__name__)
            logger.debug('There is no swarm yet. Generate the initial random swarm...')
    
            # generate the initial random swarm from the initial design
            self.init_pop = np.random.rand(popsize, dimensions) # normalized design parameters between 0 and 1

            def local_sensitivity_analysis(self, specified_initial_design_denorm):
                # 敏感性检查：以基本设计为准，检查不同的参数取极值时的电机性能变化！这是最简单有效的办法。七个设计参数，那么就有14种极值设计。
                if specified_initial_design_denorm is None:
                    initial_design_denorm = np.array( utility.Pyrhonen_design(self.im).design_parameters_denorm )
                else:
                    initial_design_denorm = specified_initial_design_denorm
                initial_design = (initial_design_denorm - min_b) / diff
                print(initial_design_denorm.tolist())
                print(initial_design.tolist())
                base_design = initial_design.tolist()
                print('base_design:', base_design, '\n-------------')
                # quit()
                number_of_variants = self.fea_config_dict['local_sensitivity_analysis_number_of_variants']
                self.init_pop = [initial_design] # include initial design!
                for i in range(len(base_design)): # 7 design parameters
                    for j in range(number_of_variants+1): # 21 variants interval
                        # copy list
                        design_variant = base_design[::]
                        design_variant[i] = j * 1./number_of_variants
                        self.init_pop.append(design_variant)
                # for ind, el in enumerate(self.init_pop):
                #     print ind, el

            if self.fea_config_dict['local_sensitivity_analysis'] == True:
                local_sensitivity_analysis(self, specified_initial_design_denorm)

            self.init_pop_denorm = min_b + self.init_pop * diff
            self.init_fitness = None

            # and save to file named gen#0000.txt
            self.number_current_generation = 0
            self.write_population_data(self.init_pop_denorm)
            logger = logging.getLogger(__name__)
            logger.debug('Initial pop (de-normalized) is saved as %s', self.dir_run + 'gen#0000.txt')

        else:
            # number_current_generation begins at 0.
            self.number_current_generation = max([int(el) for el in generations_complete])

            # 第零代是特殊的，可能跑完了，也可能没跑完，所以要把ongoing也数进来判断。
            if len(generations_complete_and_onging) == 1:
                self.init_pop_denorm = self.pop_reader(self.get_gen_file(self.number_current_generation)) # gen#0000
                self.init_fitness = None
                logger = logging.getLogger(__name__)
                logger.debug('The initial pop (i.e., gen#%s) is found in run folder. Use it.', generations_complete[0])
            else:
                # get the latest generation of swarm data and 
                # restore the living pop from file liv#xxxx.txt
                # then combine if the size of living pop is smaller than that of last-gen
                self.init_pop_denorm, self.init_fitness = self.read_completed_living_pop(self.number_current_generation) # this is completed
                                      # also read fitness of last completed living pop as init_fitness 

                logger = logging.getLogger(__name__)
                logger.debug('The latest finished generation is gen#%d', self.number_current_generation)
                logger.debug('Read in living_pop as init_pop with a size of %d.' % (len(self.init_pop_denorm)))

            # normalize for init_pop
            self.init_pop = (np.array(self.init_pop_denorm) - min_b) / diff

            if self.fea_config_dict['local_sensitivity_analysis'] == False:
                # 检查popsize是否变大了 # TODO: 文件完整性检查
                solved_popsize = len(self.init_pop)
                if popsize > solved_popsize:
                    logger.warn('The popsize is changed from last run. From %d to %d' % (solved_popsize, popsize))
                    logger.debug('init_pop before:\n' + '\n'.join(','.join('%.16f'%(x) for x in y) for y in self.init_pop))

                    # append new random designs 
                    self.init_pop = self.init_pop.tolist()
                    self.init_pop += np.random.rand(popsize-solved_popsize, dimensions).tolist() # normalized design parameters between 0 and 1
                    self.init_pop_denorm = min_b + np.array(self.init_pop) * diff

                    logger.debug('init_pop after:\n' + '\n'.join(','.join('%.16f'%(x) for x in y) for y in self.init_pop))

                    # and save to file named gen#xxxx.txt as well as liv#xxxx.txt
                    # because those just born are in the meantime living!
                    # we will over-write to file later along with fitness data
                    # self.append_population_data(self.init_pop_denorm[solved_popsize:])
                    # logger.debug('init_pop before:\n' + '\n'.join(','.join('%.16f'%(x) for x in y) for y in self.init_pop))
                elif popsize < solved_popsize:
                    logger = logging.getLogger(__name__)
                    logger.warn('Reducing popsize during a run---feature is not supported. If you are testing local sensitivity analysis, then delete the run folder (under pop/ and csv_opti/) and run again.')
                    raise Exception('Reducing popsize during a run---feature is not supported. If you are testing local sensitivity analysis, then delete the run folder (under pop/ and csv_opti/) and run again.')

        # # add a database using the swarm_data.txt file
        # if os.path.exists(self.dir_run+'swarm_data.txt'):
        #     self.database = utility.SwarmDataAnalyzer(dir_run=self.dir_run)
        #     self.database.the_generation_that_i_am_worried_about = self.number_current_generation+1
        # else: 
        #     self.database = None

        # pop: make sure the data type is array
        self.init_pop = np.asarray(self.init_pop)

        # app for jmag
        self.app = None
        self.jmag_control_state = False # indicating that by default, the jmag designer is already opened but the project file is not yet loaded or created.

        # set self.bool_run_in_JMAG_Script_Editor
        self.designer_init()

        logger = logging.getLogger(__name__)
        logger.info('Swarm is generated.')

    def designer_init(self):
        try:
            if self.app is None:
                # run inside JMAG Designer 
                import designer
                self.app = designer.GetApplication() 
                self.bool_run_in_JMAG_Script_Editor = True
        except:
            import win32com.client
            self.app = win32com.client.Dispatch('designer.Application.171')
            if self.fea_config_dict['designer.Show'] == True:
                self.app.Show()
            else:
                self.app.Hide()
            # self.app.Quit()
            self.bool_run_in_JMAG_Script_Editor = False

        def add_steel(self):
            print('[First run on this computer detected]', self.fea_config_dict['Steel'], 'is added to jmag material library.')

            if 'M15' in self.fea_config_dict['Steel']:
                add_M1xSteel(self.app, self.dir_parent, steel_name="M-15 Steel")
            elif 'M19' in self.fea_config_dict['Steel']:
                add_M1xSteel(self.app, self.dir_parent)
            elif 'Arnon5' == self.fea_config_dict['Steel']:
                add_Arnon5(self.app, self.dir_parent)        

        # too avoid tons of the same material in JAMG's material library
        if not os.path.exists(self.dir_parent + '.jmag_state.txt'):
            with open(self.dir_parent + '.jmag_state.txt', 'w') as f:
                f.write(self.fea_config_dict['pc_name'] + '/' + self.fea_config_dict['Steel'] + '\n')
            add_steel(self)
        else:
            with open(self.dir_parent + '.jmag_state.txt', 'r') as f:
                for line in f.readlines():
                    if self.fea_config_dict['pc_name'] + '/' + self.fea_config_dict['Steel'] not in line:
                        add_steel(self)

    def get_gen_file(self, no_generation, ongoing=False):
        if ongoing == True:
            return self.dir_run + 'gen#%04d-ongoing.txt'%(int(no_generation))
        else:
            return self.dir_run + 'gen#%04d.txt'%(int(no_generation))

    def get_fit_file(self, no_generation, ongoing=False):
        if ongoing == True:
            return self.dir_run + 'fit#%04d-ongoing.txt'%(int(no_generation))
        else:
            return self.dir_run + 'fit#%04d.txt'%(int(no_generation))

    def get_liv_file(self, no_generation, ongoing=False):
        if ongoing == True:
            return self.dir_run + 'liv#%04d-ongoing.txt'%(int(no_generation))
        else:
            return self.dir_run + 'liv#%04d.txt'%(int(no_generation))

    def rename_onging_files(self, no_generation):
        os.rename(  self.dir_run + 'fit#%04d-ongoing.txt'%(int(no_generation)),
                    self.dir_run + 'fit#%04d.txt'%(int(no_generation)))
        os.rename(  self.dir_run + 'gen#%04d-ongoing.txt'%(int(no_generation)),
                    self.dir_run + 'gen#%04d.txt'%(int(no_generation)))
        os.rename(  self.dir_run + 'liv#%04d-ongoing.txt'%(int(no_generation)),
                    self.dir_run + 'liv#%04d.txt'%(int(no_generation)))

    def whole_row_reader(self, reader):
        for row in reader:
            yield row[:]

    def csv_row_reader(self, handle):
        read_iterator = csv.reader(handle, skipinitialspace=True)
        return self.whole_row_reader(read_iterator)

    def show(self, which=0, toString=False):
        out_string = ''

        if which == 'all':
            for el in self.im_list:
                out_string += el.show(toString)
        else:
            self.im_list[which].show(toString)

        return out_string

    def get_project_name(self, individual_index=None):
        if self.fea_config_dict['flag_optimization'] == False or individual_index is None:
            return self.fea_config_dict['model_name_prefix']
        else:
            return self.run_folder[:-1]+'gen#%04dind#%04d' % (self.number_current_generation, individual_index)

    # @unblockPrinting
    def fobj(self, individual_index, individual_denorm):
        print('Call fobj with gen#%dind#%d'%(self.number_current_generation, individual_index))
        # based on the individual_denorm data, create design variant of the initial design of Pyrhonen09
        logger = logging.getLogger(__name__)
        mylatch = self.im
        self.im = None # to avoid to use this reference by mistake
        im_variant = bearingless_induction_motor_design.local_design_variant(mylatch, \
                        self.number_current_generation, individual_index, individual_denorm) # due to compatability issues: a new child class is used instead
        self.im_variant = im_variant # for write im_variant.ID to file corresponding to the living individual_denorm # also for command line access debug purpose
        im = im_variant # for Tran2TSS (应该给它弄个函数调用的)
        im_variant.individual_name = im_variant.get_individual_name() 

        self.project_name = self.get_project_name(individual_index)

        self.jmag_control_state = False

        # local scripts
        def open_jmag():
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # Initialize JMAG Designer
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # every model is a new project, so as to avoid the situation of 100 models in one project (occupy RAM and slow). add individual_index to project_name
            self.designer_init()
            app = self.app
            if self.jmag_control_state == False: # initilize JMAG Designer
                expected_project_file = self.dir_project_files + "%s.jproj"%(self.project_name)
                if not os.path.exists(expected_project_file):
                    app.NewProject("Untitled")
                    app.SaveAs(expected_project_file)
                    logger.debug('Create JMAG project file: %s'%(expected_project_file))
                else:
                    if self.number_current_generation > 0:
                        print('An existing model is found. This is not expected. In order to make sure femm and jmag share the exact model, the femm results will be deleted if found.')
                        if os.path.exists(self.femm_output_file_path):
                            os.remove(self.femm_output_file_path)
                            print('Removed:', self.femm_output_file_path)
                            if os.path.exists(self.femm_output_file_path[:-4]+'.fem'):
                                os.remove(self.femm_output_file_path[:-4]+'.fem')
                                print('Removed:', self.femm_output_file_path[:-4]+'.fem')
                    app.Load(expected_project_file)
                    logger.debug('Load JMAG project file: %s'%(expected_project_file))
                    logger.debug('Existing models of %d are found in %s', app.NumModels(), app.GetDefaultModelFolderPath())

                    model = app.GetCurrentModel()
                    app.DeleteModel(model.GetName())

                        # this `if' is obselete. it is used when a project contains 100 models.
                        # if app.NumModels() <= individual_index:
                        #     logger.warn('Some models are not plotted because of bad bounds (some lower bound is too small)! individual_index=%d, NumModels()=%d. See also the fit#%04d.txt file for 99999. There will be no .png file for these individuals either.', individual_index, app.NumModels(), self.number_current_generation)

                        # print app.NumStudies()
                        # print app.NumAnalysisGroups()
                        # app.SubmitAllModelsLocal() # we'd better do it one by one for easing the programing?

            self.jmag_control_state = True # indicating that the jmag project is already created
            return app

        def draw_jmag():
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # Draw the model in JMAG Designer
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            DRAW_SUCCESS = self.draw_jmag_model(individual_index, 
                                                im_variant,
                                                im_variant.individual_name)
            if DRAW_SUCCESS == 0:
                # TODO: skip this model and its evaluation
                cost_function = 99999 # penalty
                logging.getLogger(__name__).warn('Draw Failed for %s-%s\nCost function penalty = %g.%s', self.project_name, im_variant.individual_name, cost_function, self.im_variant.show(toString=True))
                raise Exception('Draw Failed: Are you working on the PC? Sometime you by mistake operate in the JMAG Geometry Editor, then it fails to draw.')
                return None
            elif DRAW_SUCCESS == -1:
                # The model already exists
                print('Model Already Exists')
                logging.getLogger(__name__).debug('Model Already Exists')
            # Tip: 在JMAG Designer中DEBUG的时候，删掉模型，必须要手动save一下，否则再运行脚本重新load project的话，是没有删除成功的，只是删掉了model的name，新导入进来的model name与project name一致。
        
            # JMAG
            if app.NumModels()>=1:
                model = app.GetModel(im_variant.individual_name)
            else:
                logger.error('there is no model yet!')
                raise Exception('why is there no model yet?')
            return model

        def exe_frequency():
            # Freq Sweeping for break-down Torque Slip
            if model.NumStudies() == 0:
                study = im_variant.add_study(app, model, self.dir_csv_output_folder, choose_study_type='frequency')
            else:
                # there is already a study. then get the first study.
                study = model.GetStudy(0)

            self.mesh_study(im_variant, app, model, study)
            self.run_study(im_variant, app, study, clock_time())

            # evaluation based on the csv results
            try:
                slip_freq_breakdown_torque, breakdown_torque, breakdown_force = self.check_csv_results(study.GetName())
            except IOError as e:
                msg = 'CJH: The solver did not exit with results, so reading the csv files reports an IO error. It is highly possible that some lower bound is too small.'
                logger.error(msg + self.im_variant.show(toString=True))
                print(msg)
                # raise e
                breakdown_torque = 0
                breakdown_force = 0

            # self.fitness_in_physics_data # torque density, torque ripple, force density, force magnitude error, force angle error, efficiency, material cost 
            return slip_freq_breakdown_torque, breakdown_torque, breakdown_force

        def load_transeint():
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # Load Results for Tran2TSS
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            results = utility.build_str_results(self.axeses, im_variant, self.project_name, tran2tss_study_name, self.dir_csv_output_folder, self.fea_config_dict, self.femm_solver)
            if results is not None:
                self.fig_main.savefig(self.dir_run + im_variant.individual_name + 'results.png', dpi=150)
                utility.pyplot_clear(self.axeses)
                # show()
            return results

        # this should be summoned even before initializing femm, and it will decide whether the femm results are reliable
        original_study_name = im_variant.individual_name + "Freq"
        self.dir_femm_temp = self.dir_csv_output_folder + 'femm_temp/'
        self.femm_output_file_path = self.dir_femm_temp + original_study_name + '.csv'
        app = open_jmag() # will set self.jmag_control_state to True

        ################################################################
        # Begin from where left: Frequency Study
        ################################################################
        # Freq Study: you can choose to not use JMAG to find the breakdown slip.
        slip_freq_breakdown_torque = None
        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        # Eddy Current Solver for Breakdown Torque and Slip
        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        if self.fea_config_dict['jmag_run_list'][0] == 0:
            # FEMM # In this case, you have to set im_variant.slip_freq_breakdown_torque by FEMM Solver
            # check for existing results

            if os.path.exists(self.femm_output_file_path):
                # raise Exception('This should never be reached')

                # # 本来想在这里多加一层保险，如果上一次优化是在运行完 wait_greedy_search 之后中断的，那么在 femm_temp/ 下就会有 ID16-5-8Freq.csv 和 ID16-5-8Freq.fem 存在。
                # # 但是，从第二代以后开始，所有在gen#文件里的个体都是随机生成的，也就是说，每次中断重新跑，这一代的这个个体都是新的 trial_denorm，所以旧的 ID16-5-8Freq.csv 和 ID16-5-8Freq.fem 必须被删除。
                # # 但是下面这样的代码是有问题，因为正常情况到了这里，self.check_csv_results(tran2tss_study_name, returnBoolean=True) 永远都应该返回 False。
                # # if not self.check_csv_results(tran2tss_study_name, returnBoolean=True):
                # #     os.remove(femm_output_file_path)
                # #     if os.path.exists(femm_output_file_path[-4:]+'.fem'):
                # #         os.remove(femm_output_file_path[-4:]+'.fem')
                # # 正确的做法是在刚刚初始化 pop 的时候检查！

                with open(self.femm_output_file_path, 'r') as f:
                    data = f.readlines()

                    slip_freq_breakdown_torque                  = float(data[0][:-1])
                    breakdown_torque                            = float(data[1][:-1])

                    self.femm_solver.stator_slot_area           = float(data[2][:-1])
                    self.femm_solver.rotor_slot_area            = float(data[3][:-1])

                    self.femm_solver.vals_results_rotor_current = []
                    for row in data[4:]:
                        index = row.find(',')
                        self.femm_solver.vals_results_rotor_current.append(float(row[:index]) + 1j*float(row[index+1:-1]))
                    # self.femm_solver.list_rotor_current_amp = [abs(el) for el in vals_results_rotor_current]
                    # print 'debug,'
                    # print self.femm_solver.vals_results_rotor_current
            else:

                # no direct returning of results, wait for it later when you need it.
                femm_tic = clock_time()
                self.femm_solver.__init__(im_variant, flag_read_from_jmag=False, freq=50.0)
                if im_variant.DriveW_poles == 2:
                    self.femm_solver.greedy_search_for_breakdown_slip( self.dir_femm_temp, original_study_name, 
                                                                        bool_run_in_JMAG_Script_Editor=self.bool_run_in_JMAG_Script_Editor, fraction=1) # 转子导条必须形成通路
                else:
                    self.femm_solver.greedy_search_for_breakdown_slip( self.dir_femm_temp, original_study_name, 
                                                                        bool_run_in_JMAG_Script_Editor=self.bool_run_in_JMAG_Script_Editor, fraction=2)

                # this is the only if path that no slip_freq_breakdown_torque is assigned!
                # this is the only if path that no slip_freq_breakdown_torque is assigned!
                # this is the only if path that no slip_freq_breakdown_torque is assigned!
        else:
            # check for existing results
            temp = self.check_csv_results(original_study_name)
            if temp is None:
                model = draw_jmag()
                slip_freq_breakdown_torque, breakdown_torque, breakdown_force = exe_frequency()
            else:
                slip_freq_breakdown_torque, breakdown_torque, breakdown_force = temp
                toc = clock_time()


        ################################################################
        # Begin from where left: Transient Study
        ################################################################
        count_tran2tss_loop = 0
        while True: # this is added in case of, e.g., iron loss data are not present while other data are already there.
            count_tran2tss_loop += 1 
            if count_tran2tss_loop == 4:
                raise Exception('This counter should not reach 4.')
            tran2tss_study_name = im_variant.individual_name + 'Tran2TSS'
            # bool_skip_transient = False

            # check whether or not the transient problem is already solved.
            if not self.check_csv_results(tran2tss_study_name, returnBoolean=True, file_suffix='_iron_loss_loss.csv'):

                if self.jmag_control_state:
                    if count_tran2tss_loop == 3:
                        # if this is reached, means the while clause is working, so we delete the existing model and draw again
                        # print('In normal case, jmag_control_state should be False here, but now it is True meaning the loop of While True is triggered. Something went wrong with JMAG solver. Loop count:', count_tran2tss_loop)
                        msg = 'The second try has failed. A manual check is needed with the jmag project %s. The most likely place that goes wrong is that a linear magnetic is somehow applied so that the iron loss solver will invariably fail.' % (self.project_name)
                        # utility.send_notification(text=msg)
                        raise Exception(msg)

                    elif count_tran2tss_loop > 1:
                        # model exists so delete it
                        model = app.GetCurrentModel()
                        app.DeleteModel(model.GetName())
                        # and draw it again
                        model = draw_jmag()

                    else:
                        model = draw_jmag()
                else:
                    raise Exception('Wrong jmag_control_state.')

                #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
                # TranFEAwi2TSS for ripples and iron loss
                #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
                # add or duplicate study for transient FEA denpending on jmag_run_list
                if self.fea_config_dict['jmag_run_list'][0] == 0:
                    if False: # this will wait for breakdown slip before setting up the pre-processor
                        if slip_freq_breakdown_torque is None:
                            # wait for femm to finish, and get your slip of breakdown
                            slip_freq_breakdown_torque, breakdown_torque, breakdown_force = self.femm_solver.wait_greedy_search(femm_tic)
                        else:
                            # femm already has results and has assign value to slip_freq_breakdown_torque
                            pass

                        # FEMM+JMAG
                        im_variant.update_mechanical_parameters(slip_freq_breakdown_torque)
                        study = im_variant.add_TranFEAwi2TSS_study( slip_freq_breakdown_torque, app, model, self.dir_csv_output_folder, tran2tss_study_name, logger)
                        self.mesh_study(im_variant, app, model, study)
                    else:

                        # FEMM+JMAG
                        study = im_variant.add_TranFEAwi2TSS_study( 50.0, app, model, self.dir_csv_output_folder, tran2tss_study_name, logger)
                        self.mesh_study(im_variant, app, model, study)

                        if slip_freq_breakdown_torque is None:
                            # wait for femm to finish, and get your slip of breakdown
                            slip_freq_breakdown_torque, breakdown_torque, breakdown_force = self.femm_solver.wait_greedy_search(femm_tic)
                        else:
                            # femm already has results and has assign value to slip_freq_breakdown_torque
                            pass
                        # Now we have the slip, set it up!
                        im_variant.update_mechanical_parameters(slip_freq_breakdown_torque) # do this for records only
                        if im_variant.the_slip != slip_freq_breakdown_torque / im_variant.DriveW_Freq:
                            raise Exception('Check update_mechanical_parameters().')
                        study.GetDesignTable().GetEquation("slip").SetExpression("%g"%(im_variant.the_slip))
                    self.run_study(im_variant, app, study, clock_time())
                else:
                    # JMAG+JMAG
                    # model = app.GetCurrentModel()
                    im_variant.update_mechanical_parameters(slip_freq_breakdown_torque)
                    if False: 
                        # this is not a good option if the excitation is different for transient study from frequency study
                        self.duplicate_TranFEAwi2TSS_from_frequency_study(im_variant, slip_freq_breakdown_torque, app, model, original_study_name, tran2tss_study_name, logger)
                    else:
                        study = im_variant.add_TranFEAwi2TSS_study( slip_freq_breakdown_torque, app, model, self.dir_csv_output_folder, tran2tss_study_name, logger)
                        self.mesh_study(im_variant, app, model, study)
                        self.run_study(im_variant, app, study, clock_time())

                # export Voltage if field data exists.
                if self.fea_config_dict['delete_results_after_calculation'] == False:
                    # Export Circuit Voltage
                    ref1 = app.GetDataManager().GetDataSet("Circuit Voltage")
                    app.GetDataManager().CreateGraphModel(ref1)
                    app.GetDataManager().GetGraphModel("Circuit Voltage").WriteTable(self.dir_csv_output_folder + im_variant.individual_name + "_EXPORT_CIRCUIT_VOLTAGE.csv")
            else:
                # if csv already exists, im_variant.slip_freq_breakdown_torque is still None till here
                # print 'DEBUG::::::::::::', im_variant.slip_freq_breakdown_torque
                if slip_freq_breakdown_torque is None:
                    # wait for femm to finish, and get your slip of breakdown
                    slip_freq_breakdown_torque, breakdown_torque, breakdown_force = self.femm_solver.wait_greedy_search(femm_tic)
                im_variant.update_mechanical_parameters(slip_freq_breakdown_torque)

            ################################################################
            # Load data for cost function evaluation
            ################################################################
            results_to_be_unpacked = load_transeint()
            if results_to_be_unpacked is not None:
                str_results, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle, jmag_loss_list, femm_loss_list, power_factor, total_loss, cost_function = results_to_be_unpacked

                # write design evaluation data to file
                with open(self.dir_run + 'swarm_data.txt', 'a') as f:
                    f.write(str_results)

                self.im = mylatch
                return cost_function
            else:
                continue

    def fobj_test(self, individual_index, individual_denorm):

        def fmodel(x, w):
            return w[0] + w[1]*x + w[2] * x**2 + w[3] * x**3 + w[4] * x**4 + w[5] * x**5 + w[6] * x**6
        def rmse(y, individual_denorm):
            y_pred = fmodel(x, individual_denorm)
            return np.sqrt(sum((y - y_pred)**2) / len(y))

        try:
            self.weights
        except:
            self.weights = [0.5*sum(el) for el in self.bounds]
            print(self.bounds)
            print(self.weights)

        x = np.linspace(0, 6.28, 50) 
        y = fmodel(x, w=self.weights)

        return rmse(y, individual_denorm)
        # plt.scatter(x, y)
        # plt.plot(x, np.cos(x), label='cos(x)')
        # plt.legend()
        # plt.show()

    def de(self):
        fobj = self.fobj
        # fobj = self.fobj_test

        if self.bool_first_time_call_de == True:
            self.bool_first_time_call_de = False

        # if self.fea_config_dict['flag_optimization'] == False:
        #     raise Exception('Please set flag_optimization before calling de.')

        # ''' 
        #     https://pablormier.github.io/2017/09/05/a-tutorial-on-differential-evolution-with-python/
        #     我对于DE的感觉，就是如果受体本身就是原最优个体的时候，如果试验体优于受体，不应该采用greedy selection把它抛弃了。
        #     read also https://nathanrooy.github.io/posts/2017-08-27/simple-differential-evolution-with-python/
        #     Notes: In step #3, this implementation differs from most DE algorithms in that we cycle through each member of the swarm, generate a donor vector, then perform selection. In this setup, every member of the swarm becomes a target vector at some point which means that every individual has the possibility of being replaced. In most DE implementations, the target vector is randomly chosen. In my experience using DE (both in aerodynamic shape optimization, as well as in deep learning applications), I have found that the current implementation works much better. The standard DE implementation might be slightly more stochastic in nature, but whatever… If, you’re interested in the “standard” DE implementation swap out the lines in the above code with the following:
        # '''
        bounds     = self.de_config_dict['bounds']
        mut        = self.de_config_dict['mut']
        crossp     = self.de_config_dict['crossp']
        popsize    = self.de_config_dict['popsize']
        iterations = self.de_config_dict['iterations']
        iterations -= self.number_current_generation # make this iterations the total number of iterations seen from the user
        logger = logging.getLogger(__name__)
        logger.debug('DE Configuration:\n\t' + '\n\t'.join('%.4f,%.4f'%tuple(el) for el in bounds) + '\tmut=%.4f, crossp=%.4f, popsize=%d, iterations=%d' % (mut,crossp,popsize,iterations) \
                     +'''\n\t# stator_tooth_width_b_ds\n\t# air_gap_length_delta\n\t# Width_RotorSlotOpen \n\t# rotor_tooth_width_b_dr \n\t# Length_HeadNeckRotorSlot\n\t# Angle_StatorSlotOpen\n\t# Width_StatorTeethHeadThickness''')
    
        self.bounds = np.array(bounds) # for debug purpose in fobj_test
    
        # mut \in  [0.5, 2.0]
        if mut < 0.5 or mut > 2.0:
            logger = logging.getLogger(__name__)
            logger.warn('Coefficient mut is generally between 0.5 and 2.0.')

        pop = self.init_pop # modification #1
        dimensions = len(bounds)
        min_b, max_b = np.asarray(bounds).T 
        # [[-5 -5 -5 -5]
        #  [ 5  5  5  5]]
        diff = np.fabs(min_b - max_b)
        pop_denorm = min_b + pop * diff
        # if self.number_current_generation == 0:
        #     print 'gen#0000 with initail design as the first individual:'
        #     for el in pop_denorm:
        #         print el.tolist()
        # 出现以下BUG：    
        # 2019-01-25 00:44:42,292 - root - ERROR - Optimization aborted.
        # Traceback (most recent call last):
        #   File "D:/OneDrive - UW-Madison/c/codes/opti_script.py", line 297, in <module>
        #     for result in de_generator:
        #   File "D:/OneDrive - UW-Madison/c/codes/population.py", line 885, in de
        #     print el.tolist()
        # IOError: [Errno 9] Bad file descriptor



        # 判断：如果是第一次，那就需要对现有pop进行生成fitness；如果续上一次的运行，则读入fitness。
        if self.init_fitness is None:
            if self.number_current_generation == 0:
                # there is no fitness file yet. run evaluation for the initial pop        
                self.jmag_control_state = False # demand to initialize the jamg designer
                fitness = np.asarray( [fobj(index, individual_denorm) for index, individual_denorm in enumerate(pop_denorm)] ) # modification #2

                logger = logging.getLogger(__name__)
                logger.debug('Generating fitness data for the initial population: %s', self.get_fit_file(self.number_current_generation))

                # write fit#0000 (fitness results to file for the initial pop)
                with open(self.get_fit_file(self.number_current_generation), 'w') as f:
                    f.write('\n'.join('%.16f'%(x) for x in fitness)) 
                # write liv_fit#0000 
                with open(self.dir_run+'liv_fit#%04d.txt'%(self.number_current_generation), 'w') as f:
                    f.write('\n'.join('%.16f'%(x) for x in fitness)) 
                # write liv#0000 (liv_id#0000 is trivial so not needed)
                self.write_population_data(pop_denorm, fname=self.get_liv_file(self.number_current_generation))
                logger.debug('Copy the 1st generation (gen#%4d) as living pop and its fitness to file.' % (self.number_current_generation))        
        else:
            # this is a continued run. load the last completed living pop's fitness data
            fitness = self.init_fitness

            if len(self.init_pop) > len(self.init_fitness):
                # 检查popsize是否变大了 # TODO: 文件完整性检查
                # in case the last iteration is not done or the popsize is incread by user after last run of optimization
                solved_popsize = len(self.init_fitness)
                print('popsize=', popsize)
                print('solved_popsize=', solved_popsize)
                if popsize > solved_popsize:
                    logger.debug('Popsize changed. New fitness for the newly come individuals should be generated...')
                    fitness_part2 = np.asarray( [fobj(index+solved_popsize, individual) for index, individual in enumerate(pop_denorm[solved_popsize:])] ) # modification #2
                
                    print('DEBUG fitness_part1:\n', fitness)
                    print('DEBUG fitness_part2:\n', fitness_part2.tolist())
                    fitness += fitness_part2.tolist()

                    # write gen#xxxx
                    self.write_population_data(pop_denorm)
                    # write fit#xxxx
                    with open(self.get_fit_file(self.number_current_generation), 'w') as f:
                        f.write('\n'.join('%.16f'%(x) for x in fitness)) 
                    # write liv_fit#xxxx
                    with open(self.dir_run+'liv_fit#%04d.txt'%(self.number_current_generation), 'w') as f:
                        f.write('\n'.join('%.16f'%(x) for x in fitness)) 
                    # write liv#xxxx (those just born are also living)
                    self.write_population_data(pop_denorm, fname=self.get_liv_file(self.number_current_generation))
                    logger.debug('At (gen#%4d) popsize changed to %d. Newly come data are written to liv_fit#, liv#, gen#, fit#' % (self.number_current_generation, popsize))
        # # debug
        # print '---------------------'
        # print fitness
        # print pop_denorm
        # quit()

        # auto-recovered run
        if self.bool_auto_recovered_run == True: 
            logger.debug('Auto-recovering the optimization...')
            self.bool_auto_recovered_run = False
            print('---------------------')
            print(pop_denorm) 
            print(self.pop_denorm)
            print('---------------------')
            print(fitness) 
            print(self.fitness)
            os.system('pause') 
            # pop_denorm = self.pop_denorm
            # fitness = self.fitness
            quit()

        # make sure fitness is an array
        fitness = np.asarray(fitness)
        best_idx = np.argmin(fitness)
        best = pop_denorm[best_idx]
            # return min_b + pop * diff, fitness, best_idx

        # for restarter or auto-recovered run
        self.pop_denorm = pop_denorm
        self.fitness    = fitness

        # # for easily identifying where the living pop is born
        # self.pop_id = [] # init
        # try:
        #     if self.popid == []:
        #         if self.im_variant is not None: 
        #             self.pop_id = [self.im_variant.ID[:-self.im_variant.ID[::-1].find('-')]+str(i) for i in range(popsize)]
        #             for el in self.pop_id:
        #                 print self.popid
        #                 self.write_living_individual_id(el)
        # except:
        #     pass

        if self.fea_config_dict['local_sensitivity_analysis'] == True:
            return
            raise Exception('The local sensitivity analysis is done.')

        # Begin DE with pop, fitness of a completed generation
        for i in range(iterations):

            self.number_current_generation += 1 # modification #3
            logger = logging.getLogger(__name__)
            logger.debug('Iteration i=%d for this execution of the script, but it is the %d DE iteration for this run %s.', i, self.number_current_generation, self.run_folder)
            # demand to initialize the jamg designer because number_current_generation has changed and a new jmag project is required.
            self.jmag_control_state = False

            for j in range(popsize): # j is the index of individual

                # get trial individual
                if i==0 and j < self.index_interrupt_beginning:
                    logger.debug('Skip individual gen#%04dind#%04d' % (self.number_current_generation, j)) 

                    # check for debugging purpose only
                    # 我们在调用 read_completed_living_pop() 的时候，init_pop（也就是这里的pop）的前半部分，是活着的living_pop，而后半部分则是用gen#xxx.txt文件中对应的位置去补充的。
                    # legacy ongoing is found during an interrupted run, so the the first iteration should continue from legacy results.
                    trial = self.ongoing_pop[j] 
                    trial_denorm = min_b + trial * diff
                    f = self.ongoing_fitness[j]

                    if f < self.ongoing_living_fitness[j]:
                        print(f, self.ongoing_living_fitness[j])
                        logger.debug('%g < %g' % (f, self.ongoing_living_fitness[j]) )
                        raise Exception('pop fit < liv fit?')
                    else:
                        # the living fitness must be smaller or the same to f
                        # print '--------------'
                        # print trial_denorm.tolist(), f, 'larger?'
                        # print pop_denorm[j].tolist(), self.ongoing_living_fitness[j], 'smaller?'
                        pop_denorm[j] = self.ongoing_living_pop_denorm[j] # read from liv-ongoing#xxxx.txt
                        pop[j] = (pop_denorm[j] - min_b) / diff
                        # print pop_denorm[j].tolist()
                        # These tolist will cause unexpected IO exception 
                        # 2019-01-19 00:14:59,193 - root - ERROR - Optimization aborted.
                        # Traceback (most recent call last):
                        #   File "K:\jchen782\c\codes/opti_script.py", line 221, in <module>
                        #     for result in de_generator:
                        #   File "K:/jchen782/c/codes/population.py", line 958, in de
                        #     print pop_denorm[j].tolist()
                        # IOError: [Errno 9] Bad file descriptor

                else:
                    logger.debug('Build individual #%d' % (j)) 

                    # # debug compare this individual to log file's line "last-liv(?):------"
                    # print 'Build individual #%d' % (j) 
                    # print (min_b + pop[j] * diff).tolist()
                    # quit()

                    idxs = [idx for idx in range(popsize) if idx != j]
                    # print 'idxs', idxs
                    a, b, c = pop[np.random.choice(idxs, 3, replace = False)] # we select three other vectors that are not the current best one, let’s call them a, b and c
                    mutant = np.clip(a + mut * (b - c), 0, 1)

                    cross_points = np.random.rand(dimensions) < crossp
                    if not np.any(cross_points): # 如果运气不好，全都不更新也不行，至少找一个位置赋值为True。
                        cross_points[np.random.randint(0, dimensions)] = True

                    # normal run
                    trial = np.where(cross_points, mutant, pop[j])
                    trial_denorm = min_b + trial * diff

                    # quit()
                    # get fitness value for the trial individual 
                    f = fobj(j, trial_denorm)

                    # write ongoing results (gen代表新随机生成的用于比较的这一代，而liv代表物竞天择以后活下来的这一代。)
                    self.write_individual_data(trial_denorm) # we write individual data after fitness is evaluated in order to make sure the two files are synchronized
                                                             # this means that the pop data file on disk does not necessarily correspondes to the current generation of pop.
                    self.write_individual_fitness(f)

                    if f < fitness[j]: # then update the pop
                        fitness[j] = f
                        pop[j] = trial # greedy selection (see de_implementation.py)

                        if f < fitness[best_idx]:
                            best_idx = j
                            best = trial_denorm

                        # try:
                        #     # for easily identifying where the living pop is born
                        #     if self.im_variant is not None:
                        #         self.write_living_individual_id(self.pop_id[j]+'->'+self.im_variant.ID)
                        # except:
                        #     pass

                    # for restart or auto-recovered run
                    self.pop_denorm[j] = min_b + pop[j] * diff
                    self.fitness[j]    = fitness[j]
                    self.write_living_individual_data(self.pop_denorm[j])
                    self.write_living_individual_fitness(fitness[j])
                    # quit()

            # one generation is finished
            self.rename_onging_files(self.number_current_generation)

            yield best, fitness[best_idx] # de verion 1
            # yield min_b + pop * diff, fitness, best_idx # de verion 2

        # TODO: 跑完一轮优化以后，必须把de_config_dict和当前的代数存在文件里，否则gen文件里面存的normalized data就没有物理意义了。

    def write_living_individual_data(self, pop_j_denorm):
        with open(self.get_liv_file(self.number_current_generation, ongoing=True), 'a') as f:
            f.write('\n' + ','.join('%.16f'%(y) for y in pop_j_denorm)) # convert 1d array to string

    def write_living_individual_fitness(self, fitness_j):
        with open(self.dir_run + 'liv_fit#%04d.txt'%(self.number_current_generation), 'a') as f:
            f.write('\n%.16f'%(fitness_j))

    def write_living_individual_id(self, str_id):
        with open(self.dir_run + 'liv_id#%04d.txt'%(self.number_current_generation), 'a') as f:
            f.write(str_id+'\n')

    def pop_reader(self, fname):
        pop = []
        with open(fname, 'r') as f:
            for row in self.csv_row_reader(f):
                if len(row)>0: # there could be empty row, since we use append mode and write down f.write('\n')
                    pop.append([float(el) for el in row])
        return pop

    def read_completed_living_pop(self, no_current_generation):
        logger = logging.getLogger(__name__)
        # if no_current_generation == 0:
        #     raise Exception('This is reached unexpectedly. However, this function returns right results.')
        #     return self.pop_reader(self.get_gen_file(no_current_generation))
        # else:
        if True:
            fname_last_gen = self.get_gen_file(no_current_generation)
            fname_last_liv = self.dir_run + 'liv#%04d.txt'%(no_current_generation)
            fname_last_liv_fit = self.dir_run + 'liv_fit#%04d.txt'%(no_current_generation)

            logger.debug('fname_last_gen=%s'% fname_last_gen)
            logger.debug('fname_last_liv=%s'% fname_last_liv)
            logger.debug('fname_last_liv_fit=%s'% fname_last_liv_fit)

        last_generation_pop_denorm = self.pop_reader(fname_last_gen)
        last_living_pop_denorm     = self.pop_reader(fname_last_liv)
        last_living_fitness        = []
        with open(fname_last_liv_fit, 'r') as f: 
            for row in self.csv_row_reader(f):
                if len(row)>0: 
                    last_living_fitness.append(float(row[0]))

        size_last =  len(last_generation_pop_denorm)
        size_living = len(last_living_pop_denorm)

        msg = 'Read in living pop of gen#%d\n'%(no_current_generation)
        msg += 'last-gen(%d):------\n\t'%(size_last)   + '\n\t'.join([str(el) for el in last_generation_pop_denorm]) + '\n'
        msg += 'last-liv(%d):------\n\t'%(size_living) + '\n\t'.join([str(el) for el in last_living_pop_denorm]) + '\n'

        try: 
            self.ongoing_pop_denorm
        except:
            logger.debug('The last execution is interrupted with a clean slate.')
            # self.ongoing_pop_denorm不存在。
            # 这种情况，说明优化脚本中断的点，刚好是上一代刚好跑完的时候，这个时候dir_run下的ongoing文件都被删掉了。那么我们就以上一代的living_pop作为这一代的ongoing_pop。
            self.ongoing_pop_denorm        = np.array(last_generation_pop_denorm)
            self.ongoing_living_pop_denorm = last_living_pop_denorm

        size_ongoing_last   = len(self.ongoing_pop_denorm       )
        self.size_ongoing_living = len(self.ongoing_living_pop_denorm)
        msg += 'ongoing-gen(%d):------\n\t'%(size_ongoing_last)        + '\n\t'.join([str(el) for el in self.ongoing_pop_denorm.tolist()]) + '\n'
        msg += 'ongoing-liv(%d):------\n\t'%(self.size_ongoing_living) + '\n\t'.join([str(el) for el in self.ongoing_living_pop_denorm]) + '\n'

        if size_living < size_last:
            raise Exception('size_living < size_last')
            #     msg += 'Append\n\t' + '\n\t'.join([str(el) for el in last_generation_pop_denorm[-(size_last - size_living):]]) + '\n'
            #     last_living_pop_denorm += last_generation_pop_denorm[-(size_last - size_living):]
            # msg += 'combined:------\n\t' + '\n\t'.join([str(el) for el in last_living_pop_denorm]) + '\n'
        if len(last_living_pop_denorm) != len(last_generation_pop_denorm):
            raise Exception('living and last gen do not match')

        logger.debug(msg)
        print('last_living_fitness=', last_living_fitness)
        return last_living_pop_denorm, last_living_fitness

    def write_population_data(self, pop, fname=None):
        if fname is None:
            with open(self.get_gen_file(self.number_current_generation), 'w') as f:
                f.write('\n'.join(','.join('%.16f'%(x) for x in y) for y in pop)) # convert 2d array to string
        else:
            with open(fname, 'w') as f:
                f.write('\n'.join(','.join('%.16f'%(x) for x in y) for y in pop)) # convert 2d array to string        

    def append_population_data(self, pop): # for increased popsize from last run
        with open(self.get_gen_file(self.number_current_generation), 'a') as f:
            f.write('\n')
            f.write('\n'.join(','.join('%.16f'%(x) for x in y) for y in pop)) # convert 2d array to string

    def write_individual_data(self, trial_individual):
        with open(self.get_gen_file(self.number_current_generation, ongoing=True), 'a') as f:
            f.write('\n' + ','.join('%.16f'%(y) for y in trial_individual)) # convert 1d array to string

    def write_individual_fitness(self, fitness_scalar):
        with open(self.get_fit_file(self.number_current_generation, ongoing=True), 'a') as f:
            f.write('\n%.16f'%(fitness_scalar))



    def draw_jmag_model(self, individual_index, im_variant, model_name, bool_trimDrawer_or_vanGogh=True, doNotRotateCopy=False):

        if individual_index == -1: # 后处理是-1
            print('Draw model for post-processing')
            if individual_index+1 + 1 <= self.app.NumModels():
                logger = logging.getLogger(__name__)
                logger.debug('The model already exists for individual with index=%d. Skip it.', individual_index)
                return -1 # the model is already drawn

        elif individual_index+1 <= self.app.NumModels(): # 一般是从零起步
            logger = logging.getLogger(__name__)
            logger.debug('The model already exists for individual with index=%d. Skip it.', individual_index)
            return -1 # the model is already drawn

        # open JMAG Geometry Editor
        self.app.LaunchGeometryEditor()
        geomApp = self.app.CreateGeometryEditor()
        # geomApp.Show()
        geomApp.NewDocument()
        doc = geomApp.GetDocument()
        ass = doc.GetAssembly()

        # draw parts
        try:
            if bool_trimDrawer_or_vanGogh:
                d = TrimDrawer(im_variant) # 传递的是地址哦
                d.doc, d.ass = doc, ass
                d.plot_shaft("Shaft")

                d.plot_rotorCore("Rotor Core")
                d.plot_cage("Cage")

                d.plot_statorCore("Stator Core")

                d.plot_coil("Coil")
                # d.plot_airWithinRotorSlots(u"Air Within Rotor Slots")
            else:
                d = VanGogh_JMAG(im_variant, doNotRotateCopy=doNotRotateCopy) # 传递的是地址哦
                d.doc, d.ass = doc, ass
                d.draw_model()
            self.d = d
        except Exception as e:
            print('See log file to plotting error.')
            logger = logging.getLogger(__name__)
            logger.error('The drawing is terminated. Please check whether the specified bounds are proper.', exc_info=True)

            # print 'Draw Failed'
            # if self.pc_name == 'Y730':
            #     # and send the email to hory chen
            #     raise e

            # or you can skip this model and continue the optimization!
            return False # indicating the model cannot be drawn with the script.

        # Import Model into Designer
        doc.SaveModel(True) # True=on : Project is also saved. 
        model = self.app.GetCurrentModel() # model = self.app.GetModel(u"IM_DEMO_1")
        model.SetName(model_name)
        model.SetDescription(im_variant.model_name_prefix + '\n' + im_variant.show(toString=True))

        if doNotRotateCopy:
            im_variant.pre_process_structural(self.app, self.d.listKeyPoints)
        else:
            im_variant.pre_process(self.app)

        model.CloseCadLink() # this is essential if you want to create a series of models
        return True

    def check_csv_results(self, study_name, returnBoolean=False, file_suffix='_torque.csv'): # '_iron_loss_loss.csv'
        # print self.dir_csv_output_folder + study_name + '_torque.csv'
        if not os.path.exists(self.dir_csv_output_folder + study_name + file_suffix):
            if returnBoolean == False:
                print('Nothing is found when looking into:', self.dir_csv_output_folder + study_name + file_suffix)
                return None
            else:
                return False
        else:
            if returnBoolean == True:
                return True

        try:
            # check csv results 
            l_slip_freq = []
            l_TorCon    = []
            l_ForCon_X  = []
            l_ForCon_Y  = []

            self.fitness_in_physics_data = []
            with open(self.dir_csv_output_folder + study_name + '_torque.csv', 'r') as f: 
                for ind, row in enumerate(self.csv_row_reader(f)):
                    if ind >= 5:
                        try:
                            float(row[0])
                        except:
                            continue
                        l_slip_freq.append(float(row[0]))
                        l_TorCon.append(float(row[1]))

            with open(self.dir_csv_output_folder + study_name + '_force.csv', 'r') as f: 
                for ind, row in enumerate(self.csv_row_reader(f)):
                    if ind >= 5:
                        try:
                            float(row[0])
                        except:
                            continue
                        # l_slip_freq.append(float(row[0]))
                        l_ForCon_X.append(float(row[1]))
                        l_ForCon_Y.append(float(row[2]))

            # self.fitness_in_physics_data.append(l_slip_freq)
            # self.fitness_in_physics_data.append(l_TorCon)

            breakdown_force = max(np.sqrt(np.array(l_ForCon_X)**2 + np.array(l_ForCon_Y)**2))

            index, breakdown_torque = utility.get_index_and_max(l_TorCon)
            slip_freq_breakdown_torque = l_slip_freq[index]
            return slip_freq_breakdown_torque, breakdown_torque, breakdown_force
        except NameError as e:
            logger = logging.getLogger(__name__)
            logger.error('No CSV File Found.', exc_info=True)
            raise e

    def plot_csv_results_for_all(self):
        # 想要观察所有解的滑差分布情况，合理规划滑差步长的选择以及寻找范围！
        from pylab import plot, show, legend, grid, figure
        max_breakdown_torque = 0.0
        max_breakdown_torque_file = None
        max_breakdown_force_amp = 0.0
        max_breakdown_force_amp_file = None
        # self.fitness_in_physics_data = []
        for file in os.listdir(self.dir_csv_output_folder):
            if 'lock' not in file and self.model_name_prefix in file: # check model_name_prefix in file because there is D:/Users/horyc/OneDrive - UW-Madison/csv_opti/run#1/Freq_#32-0-0-FFVRC-RCR_force.csv
                if '_torque' in file:
                    l_slip_freq, l_TorCon = self.read_csv_results(self.dir_csv_output_folder + file)
                    figure(1)
                    plot(l_slip_freq, l_TorCon, label=file)
                    temp = max(l_TorCon)
                    if temp > max_breakdown_torque:
                        max_breakdown_torque = temp
                        max_breakdown_torque_file = file
                elif '_force' in file:
                    l_slip_freq_2, l_ForCon_XY = self.read_csv_results(self.dir_csv_output_folder + file)
                    l_ForCon_X = [el[0] for el in l_ForCon_XY]
                    l_ForCon_Y = [el[1] for el in l_ForCon_XY]
                    figure(2)
                    plot(l_slip_freq_2, l_ForCon_X, label=file)
                    figure(3)
                    plot(l_slip_freq_2, l_ForCon_Y, label=file)
                    temp = max(np.sqrt(np.array(l_ForCon_X)**2 + np.array(l_ForCon_Y)**2))
                    if temp > max_breakdown_force_amp:
                        max_breakdown_force_amp = temp
                        max_breakdown_force_amp_file = file
        print(max_breakdown_torque, max_breakdown_torque_file, type(max_breakdown_torque))
        print(max_breakdown_force_amp, max_breakdown_force_amp_file, type(max_breakdown_force_amp))
        for i in range(1, 3+1):
            figure(i)
            legend()
            grid()
        show()

        ''' the slip? the individual? I guess max torque and force are produced with the smallest air gap!
        '''

    def read_csv_results(self, file_location):
        l_slip_freq = [] # or other x-axis variable such as time and angle.
        l_data    = []
        with open(file_location, 'r') as f: 
            for ind, row in enumerate(self.csv_row_reader(f)):
                if ind >= 5:
                    # print file_location
                    # print row
                    l_slip_freq.append(float(row[0]))
                    l_data.append([float(el) for el in row[1:]])
        return l_slip_freq, l_data



    def logging_1d_array_for_debug(self, a, a_name):
        # self.logging_1d_array_for_debug(trial)
        logger = logging.getLogger(__name__)
        logger.debug(a_name + ','.join('%.16f'%(y) for y in a)) # convert 1d array to string

    ''' API for FEMM Solver
    '''
    def get_breakdown_results(self, im, study_type=1):
        if study_type == 1:
            study_name = "Freq"
        else:
            raise Exception('not supported study_type.')

        return self.check_csv_results(study_name)

    def has_results(self, im, study_type='Freq'):
        # short study_name because of one model one project convension
        if   study_type == 'Freq': 
            study_name =  "Freq"
        elif study_type == 'Tran2TSS':
            study_name =  "Tran2TSS"
        elif study_type == 'Freq-FFVRC':
            study_name =  "Freq-FFVRC"
        else:
            raise Exception('not supported study_type.')

        # im.csv_previous_solve = self.dir_csv_output_folder + study_name + '_circuit_current.csv'
        im.csv_previous_solve = self.dir_csv_output_folder + im.get_individual_name()+ study_name + '_circuit_current.csv'
        bool_temp = os.path.exists(im.csv_previous_solve)

        if study_type == 'Freq':
            if bool_temp == True:
                # this is no longer needed, because the TranFEAwi2TSS will do this for ya.
                # slip_freq_breakdown_torque, _, _ = self.get_breakdown_results(im)
                slip_freq_breakdown_torque, _, _ = self.check_csv_results(im.get_individual_name()+ study_name)

                im.update_mechanical_parameters(slip_freq_breakdown_torque)
                # print 'slip_freq_breakdown_torque:', slip_freq_breakdown_torque, 'Hz'

        return bool_temp

    def run(self, im, individual_index=0, run_list=[1,1,0,0,0]):
        ''' run_list: toggle solver for Freq, Tran2TSS, Freq-FFVRC, TranRef, Static'''
        self.im = None
        logger = logging.getLogger(__name__)

        # Settings 
        self.jmag_control_state = False # new one project one model convension

        # initialize JMAG Designer
        self.designer_init()
        app = self.app
        self.project_name = self.get_project_name(individual_index=individual_index)
        # self.project_name = self.get_project_name()
        im.model_name = im.get_individual_name() 
        if self.jmag_control_state == False: # initilize JMAG Designer
            expected_project_file = self.dir_project_files + "%s.jproj"%(self.project_name)
            if not os.path.exists(expected_project_file):
                app.NewProject("Untitled")
                app.SaveAs(expected_project_file)
                logger.debug('Create JMAG project file: %s'%(expected_project_file))
            else:
                app.Load(expected_project_file)
                logger.debug('Load JMAG project file: %s'%(expected_project_file))
                logger.debug('Existing models of %d are found in %s', app.NumModels(), app.GetDefaultModelFolderPath())

                # if app.NumModels() <= individual_index:
                #     logger.warn('Some models are not plotted because of bad bounds (some lower bound is too small)! individual_index=%d, NumModels()=%d. See also the fit#%04d.txt file for 99999. There will be no .png file for these individuals either.', individual_index, app.NumModels(), self.number_current_generation)

                # print app.NumStudies()
                # print app.NumAnalysisGroups()
                # app.SubmitAllModelsLocal() # we'd better do it one by one for easing the programing?

        # draw the model in JMAG Designer
        DRAW_SUCCESS = self.draw_jmag_model( individual_index, 
                                        im,
                                        im.model_name)
        # print 'TEST VanGogh for JMAG.'
        # return
        self.jmag_control_state = True # indicating that the jmag project is already created
        if DRAW_SUCCESS == 0:
            # TODO: skip this model and its evaluation
            cost_function = 99999
            logger.warn('Draw Failed for'+'%s-%s: %g', self.project_name, im.model_name, cost_function)
            return cost_function
        elif DRAW_SUCCESS == -1:
            # The model already exists
            print('Model Already Exists')
            logger.debug('Model Already Exists')
        # Tip: 在JMAG Designer中DEBUG的时候，删掉模型，必须要手动save一下，否则再运行脚本重新load project的话，是没有删除成功的，只是删掉了model的name，新导入进来的model name与project name一致。
        if app.NumModels()>=1:
            model = app.GetModel(im.model_name)
        else:
            logger.error('there is no model yet!')
            print('why is there no model yet?')

        # Freq Sweeping for break-down Torque Slip
        # remember to export the B data using subroutine 
        # and check export table results only
        if model.NumStudies() == 0:
            study = im.add_study(app, model, self.dir_csv_output_folder, choose_study_type='frequency')
        else:
            # there is already a study. then get the first study.
            study = model.GetStudy(0)

        # Freq Study: you can choose to not use JMAG to find the breakdown slip.
        # In this case, you have to set im.slip_freq_breakdown_torque by FEMM Solver
        print('debug: run_list[0]', run_list[0])
        if run_list[0] == 0:
            # Does femm has already done the frequency sweeping for breakdown torque?
            if im.slip_freq_breakdown_torque is None:
                raise Exception('run_list[0] is 0, so you have to run FEMM solver first to get slip_freq_breakdown_torque')
            print("FEMM's slip_freq_breakdown_torque is used for Tran2TSS.")
            slip_freq_breakdown_torque = im.slip_freq_breakdown_torque
        else:
            # Use JMAG to sweeping the frequency
            # Does study has results?
            if study.AnyCaseHasResult():
                slip_freq_breakdown_torque, breakdown_torque, breakdown_force = self.check_csv_results(study.GetName())
            else:
                # mesh
                im.add_mesh(study, model)

                # Export Image
                    # for i in range(app.NumModels()):
                    #     app.SetCurrentModel(i)
                    #     model = app.GetCurrentModel()
                    #     app.ExportImage(r'D:\Users\horyc\OneDrive - UW-Madison\pop\run#10/' + model.GetName() + '.png')
                app.View().ShowAllAirRegions()
                # app.View().ShowMeshGeometry() # 2nd btn
                app.View().ShowMesh() # 3rn btn
                app.View().Zoom(3)
                app.View().Pan(-im.Radius_OuterRotor, 0)
                app.ExportImageWithSize(self.dir_run + model.GetName() + '.png', 2000, 2000)
                app.View().ShowModel() # 1st btn. close mesh view, and note that mesh data will be deleted because only ouput table results are selected.

                # run
                study.RunAllCases()
                app.Save()

                # evaluation based on the csv results
                slip_freq_breakdown_torque, breakdown_torque, breakdown_force = self.check_csv_results(study.GetName())
                # self.fitness_in_physics_data # torque density, torque ripple, force density, force magnitude error, force angle error, efficiency, material cost 

        # this will be used for other duplicated studies
        original_study_name = study.GetName()
        im.csv_previous_solve = self.dir_csv_output_folder + original_study_name + '_circuit_current.csv'
        im.update_mechanical_parameters(slip_freq_breakdown_torque, syn_freq=im.DriveW_Freq)


        # EC Rotate
        tic = clock_time()
        ecrot_study_name = original_study_name + "-FFVRC"
        if model.NumStudies()<2:
            # EC Rotate: Rotate the rotor to find the ripples in force and torque # 不关掉这些云图，跑第二个study的时候，JMAG就挂了：app.View().SetVectorView(False); app.View().SetFluxLineView(False); app.View().SetContourView(False)
            casearray = [0 for i in range(1)]
            casearray[0] = 1
            model.DuplicateStudyWithCases(original_study_name, ecrot_study_name, casearray)

            app.SetCurrentStudy(ecrot_study_name)
            study = app.GetCurrentStudy()
            divisions_per_slot_pitch = self.fea_config_dict['ec_rotate_divisions_per_slot_pitch']  # 24
            study.GetStep().SetValue("Step", divisions_per_slot_pitch) 
            study.GetStep().SetValue("StepType", 0)
            study.GetStep().SetValue("FrequencyStep", 0)
            study.GetStep().SetValue("Initialfrequency", slip_freq_breakdown_torque)

                # study.GetCondition(u"RotCon").SetValue(u"MotionGroupType", 1)
            study.GetCondition("RotCon").SetValue("Displacement", + 360.0/im.Qr/divisions_per_slot_pitch)

            # https://www2.jmag-international.com/support/en/pdf/JMAG-Designer_Ver.17.1_ENv3.pdf
            study.GetStudyProperties().SetValue("DirectSolverType", 1)

            if run_list[2] == True:
                # model.RestoreCadLink()
                study.Run()
                app.Save()
                # model.CloseCadLink()
            else:
                pass # if the jcf file already exists, it pops a msg window
                # study.WriteAllSolidJcf(self.dir_jcf, im.model_name+study.GetName()+'Solid', True) # True : Outputs cases that do not have results 
                # study.WriteAllMeshJcf(self.dir_jcf, im.model_name+study.GetName()+'Mesh', True)
        logger.debug('EC-Rotate spent %g sec.'%(clock_time() - tic))


        # Transient FEA wi 2 Time Step Section
        tic = clock_time()
        tran2tss_study_name = original_study_name[:-4] + "Tran2TSS"
        if model.NumStudies()<3:
            model.DuplicateStudyWithType(original_study_name, "Transient2D", tran2tss_study_name)
            app.SetCurrentStudy(tran2tss_study_name)
            study = app.GetCurrentStudy()

            # 上一步的铁磁材料的状态作为下一步的初值，挺好，但是如果每一个转子的位置转过很大的话，反而会减慢非线性迭代。
            # 我们的情况是：0.33 sec 分成了32步，每步的时间大概在0.01秒，0.01秒乘以0.5*497 Hz = 2.485 revolution...
            study.GetStudyProperties().SetValue("NonlinearSpeedup", 0) 

            # 2 sections of different time step
            # IEMDC
            number_cycles_prolonged = self.fea_config_dict['number_cycles_prolonged'] # 150 or 1
            number_of_steps_1stTTS = self.fea_config_dict['number_of_steps_1stTTS'] 
            number_of_steps_2ndTTS = self.fea_config_dict['number_of_steps_2ndTTS'] 
            DM = app.GetDataManager()
            DM.CreatePointArray("point_array/timevsdivision", "SectionStepTable")
            if number_cycles_prolonged == 0: # 不延长了，就是普通的Tran2TSS
                refarray = [[0 for i in range(3)] for j in range(3)]
                refarray[0][0] = 0
                refarray[0][1] =    1
                refarray[0][2] =        50
                refarray[1][0] = 0.5/slip_freq_breakdown_torque
                refarray[1][1] =    number_of_steps_1stTTS
                refarray[1][2] =        50
                refarray[2][0] = refarray[1][0] + 0.5/im.DriveW_Freq
                refarray[2][1] =    number_of_steps_2ndTTS  # also modify range_ss! # don't forget to modify below!
                refarray[2][2] =        50
                number_of_total_steps = 1 + number_of_steps_1stTTS + number_of_steps_2ndTTS # [Double Check] don't forget to modify here!
            else:
                refarray = [[0 for i in range(3)] for j in range(5)]
                refarray[0][0] = 0
                refarray[0][1] =    1
                refarray[0][2] =        50
                refarray[1][0] = 0.5/slip_freq_breakdown_torque
                refarray[1][1] =    number_of_steps_1stTTS
                refarray[1][2] =        50
                refarray[2][0] = refarray[1][0] + 0.5/im.DriveW_Freq
                refarray[2][1] =    number_of_steps_2ndTTS  # also modify range_ss! # don't forget to modify below!
                refarray[2][2] =        50
                refarray[3][0] = refarray[2][0] + number_cycles_prolonged/im.DriveW_Freq # =50*0.002 sec = 0.1 sec is needed to converge to TranRef
                refarray[3][1] =    number_cycles_prolonged*self.fea_config_dict['TranRef-StepPerCycle'] # =50*40, every 0.002 sec takes 40 steps 
                refarray[3][2] =        50
                refarray[4][0] = refarray[3][0] + 0.5/im.DriveW_Freq # 最后来一个超密的半周期400步
                refarray[4][1] =    400
                refarray[4][2] =        50
                number_of_total_steps = 1 + number_of_steps_1stTTS + number_of_steps_2ndTTS + number_cycles_prolonged*self.fea_config_dict['TranRef-StepPerCycle'] + 400 # [Double Check] don't forget to modify here!
            DM.GetDataSet("SectionStepTable").SetTable(refarray)
            study.GetStep().SetValue("Step", number_of_total_steps)
            study.GetStep().SetValue("StepType", 3)
            study.GetStep().SetTableProperty("Division", DM.GetDataSet("SectionStepTable"))

            # add equations
            study.GetDesignTable().AddEquation("freq")
            study.GetDesignTable().AddEquation("slip")
            study.GetDesignTable().AddEquation("speed")
            study.GetDesignTable().GetEquation("freq").SetType(0)
            study.GetDesignTable().GetEquation("freq").SetExpression("%g"%((im.DriveW_Freq)))
            study.GetDesignTable().GetEquation("freq").SetDescription("Excitation Frequency")
            study.GetDesignTable().GetEquation("slip").SetType(0)
            study.GetDesignTable().GetEquation("slip").SetExpression("%g"%(im.the_slip))
            study.GetDesignTable().GetEquation("slip").SetDescription("Slip [1]")
            study.GetDesignTable().GetEquation("speed").SetType(1)
            study.GetDesignTable().GetEquation("speed").SetExpression("freq * (1 - slip) * %d"%(60/(im.DriveW_poles/2))) ######################## 4 pole motor
            study.GetDesignTable().GetEquation("speed").SetDescription("mechanical speed of four pole")

            # speed, freq, slip
            study.GetCondition("RotCon").SetValue("AngularVelocity", 'speed')
            app.ShowCircuitGrid(True)
            study.GetCircuit().GetComponent("CS4").SetValue("Frequency", "freq")
            study.GetCircuit().GetComponent("CS2").SetValue("Frequency", "freq")

            # max_nonlinear_iteration = 50
            # study.GetStudyProperties().SetValue(u"NonlinearMaxIteration", max_nonlinear_iteration)
            study.GetStudyProperties().SetValue("ApproximateTransientAnalysis", 1) # psuedo steady state freq is for PWM drive to use
            study.GetStudyProperties().SetValue("SpecifySlip", 1)
            study.GetStudyProperties().SetValue("OutputSteadyResultAs1stStep", 0)
            study.GetStudyProperties().SetValue("Slip", "slip")

            # # add other excitation frequencies other than 500 Hz as cases
            # for case_no, DriveW_Freq in enumerate([50.0, slip_freq_breakdown_torque]):
            #     slip = slip_freq_breakdown_torque / DriveW_Freq
            #     study.GetDesignTable().AddCase()
            #     study.GetDesignTable().SetValue(case_no+1, 0, DriveW_Freq)
            #     study.GetDesignTable().SetValue(case_no+1, 1, slip)

            if False: # to save some time: you study is prolonged!
                # # Iron Loss Calculation Condition (Legacy Codes)
                # # Stator 
                # cond = study.CreateCondition(u"Ironloss", u"IronLossConStator")
                # cond.SetValue(u"RevolutionSpeed", u"freq*60/%d"%(0.5*(im.DriveW_poles)))
                # cond.ClearParts()
                # sel = cond.GetSelection()
                # sel.SelectPartByPosition(-im.Radius_OuterStatorYoke+EPS, 0 ,0)
                # cond.AddSelected(sel)
                # # Use FFT for hysteresis to be consistent with FEMM's results
                # cond.SetValue(u"HysteresisLossCalcType", 1)
                # cond.SetValue(u"PresetType", 3)
                # # Rotor
                # cond = study.CreateCondition(u"Ironloss", u"IronLossConRotor")
                # cond.SetValue(u"BasicFrequencyType", 2)
                # cond.SetValue(u"BasicFrequency", u"slip*freq")
                # cond.ClearParts()
                # sel = cond.GetSelection()
                # sel.SelectPartByPosition(-im.Radius_Shaft-EPS, 0 ,0)
                # cond.AddSelected(sel)
                # # Use FFT for hysteresis to be consistent with FEMM's results
                # cond.SetValue(u"HysteresisLossCalcType", 1)
                # cond.SetValue(u"PresetType", 3)

                # Iron Loss Calculation Condition (Working codes from optimization)
                # Stator 
                im_variant = im
                if True:
                    cond = study.CreateCondition("Ironloss", "IronLossConStator")
                    cond.SetValue("RevolutionSpeed", "freq*60/%d"%(0.5*(im_variant.DriveW_poles)))
                    cond.ClearParts()
                    sel = cond.GetSelection()
                    sel.SelectPartByPosition(-im_variant.Radius_OuterStatorYoke+EPS, 0 ,0)
                    cond.AddSelected(sel)
                    # Use FFT for hysteresis to be consistent with FEMM's results and to have a FFT plot
                    cond.SetValue("HysteresisLossCalcType", 1)
                    cond.SetValue("PresetType", 3) # 3:Custom
                    # Specify the reference steps yourself because you don't really know what JMAG is doing behind you
                    cond.SetValue("StartReferenceStep", number_of_total_steps+1-number_of_steps_2ndTTS*0.5) # 1/4 period <=> number_of_steps_2ndTTS*0.5
                    cond.SetValue("EndReferenceStep", number_of_total_steps)
                    cond.SetValue("UseStartReferenceStep", 1)
                    cond.SetValue("UseEndReferenceStep", 1)
                    cond.SetValue("Cyclicity", 4) # specify reference steps for 1/4 period and extend it to whole period
                    cond.SetValue("UseFrequencyOrder", 1)
                    cond.SetValue("FrequencyOrder", "1-50") # Harmonics up to 50th orders 
                # Check CSV reults for iron loss (You cannot check this for Freq study) # CSV and save space
                study.GetStudyProperties().SetValue("CsvOutputPath", self.dir_csv_output_folder) # it's folder rather than file!
                study.GetStudyProperties().SetValue("CsvResultTypes", "Torque;Force;LineCurrent;TerminalVoltage;JouleLoss;TotalDisplacementAngle;JouleLoss_IronLoss;IronLoss_IronLoss;HysteresisLoss_IronLoss")
                study.GetStudyProperties().SetValue("DeleteResultFiles", self.fea_config_dict['delete_results_after_calculation'])
                # Terminal Voltage/Circuit Voltage: Check for outputing CSV results 
                study.GetCircuit().CreateTerminalLabel("Terminal4U", 8, -13)
                study.GetCircuit().CreateTerminalLabel("Terminal4V", 8, -11)
                study.GetCircuit().CreateTerminalLabel("Terminal4W", 8, -9)
                study.GetCircuit().CreateTerminalLabel("Terminal2U", 23, -13)
                study.GetCircuit().CreateTerminalLabel("Terminal2V", 23, -11)
                study.GetCircuit().CreateTerminalLabel("Terminal2W", 23, -9)
                # Export Stator Core's field results only for iron loss calculation (the csv file of iron loss will be clean with this setting)
                    # study.GetMaterial(u"Rotor Core").SetValue(u"OutputResult", 0) # at least one part on the rotor should be output or else a warning "the jplot file does not contains displacement results when you try to calc. iron loss on the moving part." will pop up, even though I don't add iron loss condition on the rotor.
                # study.GetMeshControl().SetValue(u"AirRegionOutputResult", 0)
                study.GetMaterial("Shaft").SetValue("OutputResult", 0)
                study.GetMaterial("Cage").SetValue("OutputResult", 0)
                study.GetMaterial("Coil").SetValue("OutputResult", 0)
                # Rotor
                if True:
                    cond = study.CreateCondition("Ironloss", "IronLossConRotor")
                    cond.SetValue("BasicFrequencyType", 2)
                    cond.SetValue("BasicFrequency", "freq")
                        # cond.SetValue(u"BasicFrequency", u"slip*freq") # this require the signal length to be at least 1/4 of slip period, that's too long!
                    cond.ClearParts()
                    sel = cond.GetSelection()
                    sel.SelectPartByPosition(-im_variant.Radius_Shaft-EPS, 0 ,0)
                    cond.AddSelected(sel)
                    # Use FFT for hysteresis to be consistent with FEMM's results
                    cond.SetValue("HysteresisLossCalcType", 1)
                    cond.SetValue("PresetType", 3)
                    # Specify the reference steps yourself because you don't really know what JMAG is doing behind you
                    cond.SetValue("StartReferenceStep", number_of_total_steps+1-number_of_steps_2ndTTS*0.5) # 1/4 period <=> number_of_steps_2ndTTS*0.5
                    cond.SetValue("EndReferenceStep", number_of_total_steps)
                    cond.SetValue("UseStartReferenceStep", 1)
                    cond.SetValue("UseEndReferenceStep", 1)
                    cond.SetValue("Cyclicity", 4) # specify reference steps for 1/4 period and extend it to whole period
                    cond.SetValue("UseFrequencyOrder", 1)
                    cond.SetValue("FrequencyOrder", "1-50") # Harmonics up to 50th orders 


            # https://www2.jmag-international.com/support/en/pdf/JMAG-Designer_Ver.17.1_ENv3.pdf
            study.GetStudyProperties().SetValue("DirectSolverType", 1)

            if run_list[1] == True:
                study.RunAllCases()
                app.Save()
            else:
                pass # if the jcf file already exists, it pops a msg window
                # study.WriteAllSolidJcf(self.dir_jcf, im.model_name+study.GetName()+'Solid', True) # True : Outputs cases that do not have results 
                # study.WriteAllMeshJcf(self.dir_jcf, im.model_name+study.GetName()+'Mesh', True)
        logger.debug('Tran2TSSProlong spent %g sec.'%(clock_time() - tic))



        # These two studies are not needed after iemdc
        if True:
            # Transient Reference
            tranRef_study_name = "TranRef"
            if model.NumStudies()<4:
                model.DuplicateStudyWithType(tran2tss_study_name, "Transient2D", tranRef_study_name)
                app.SetCurrentStudy(tranRef_study_name)
                study = app.GetCurrentStudy()

                # 将一个滑差周期和十个同步周期，分成 400 * end_point / (1.0/im.DriveW_Freq) 份。
                end_point = 0.5/slip_freq_breakdown_torque + 10.0/im.DriveW_Freq
                # Pavel Ponomarev 推荐每个电周期400~600个点来捕捉槽效应。
                division = self.fea_config_dict['TranRef-StepPerCycle'] * end_point / (1.0/im.DriveW_Freq)  # int(end_point * 1e4)
                                                                        # end_point = division * 1e-4
                study.GetStep().SetValue("Step", division + 1) 
                study.GetStep().SetValue("StepType", 1) # regular inverval
                study.GetStep().SetValue("StepDivision", division)
                study.GetStep().SetValue("EndPoint", end_point)

                # https://www2.jmag-international.com/support/en/pdf/JMAG-Designer_Ver.17.1_ENv3.pdf
                study.GetStudyProperties().SetValue("DirectSolverType", 1)

                if run_list[3] == True:
                    study.RunAllCases()
                    app.Save()
                else:
                    pass # if the jcf file already exists, it pops a msg window
                    # study.WriteAllSolidJcf(self.dir_jcf, im.model_name+study.GetName()+'Solid', True) # True : Outputs cases that do not have results 
                    # study.WriteAllMeshJcf(self.dir_jcf, im.model_name+study.GetName()+'Mesh', True)

            # Rotating Static FEA (This can be done in FEMM)
            if run_list[4] == True:

                im.MODEL_ROTATE = True # this is used in Rotating Static FEA
                im.total_number_of_cases = 1 # 这个值取24的话，就可以得到24个不同位置下，电机的转矩-滑差曲线了，这里取1，真正的cases在StaticFEA中添加。

                # draw another model with MODEL_ROTATE as True
                if True:
                    DRAW_SUCCESS = self.draw_jmag_model( 1, # +1
                                                    im,
                                                    im.model_name + 'MODEL_ROTATE')
                    self.jmag_control_state = True # indicating that the jmag project is already created
                    if DRAW_SUCCESS == 0:
                        # TODO: skip this model and its evaluation
                        cost_function = 99999
                        return cost_function
                    elif DRAW_SUCCESS == -1:
                        # The model already exists
                        print('Model Already Exists')
                        logging.getLogger(__name__).debug('Model Already Exists')
                    # Tip: 在JMAG Designer中DEBUG的时候，删掉模型，必须要手动save一下，否则再运行脚本重新load project的话，是没有删除成功的，只是删掉了model的name，新导入进来的model name与project name一致。
                    if app.NumModels()>=2: # +1
                        model = app.GetModel(im.model_name + 'MODEL_ROTATE')
                    else:
                        logging.getLogger(__name__).error('there is no model yet!')
                        print('why is there no model yet?')
                        raise

                    if model.NumStudies() == 0:
                        study = im.add_study(app, model, self.dir_csv_output_folder, choose_study_type='static')
                    else:
                        # there is already a study. then get the first study.
                        study = model.GetStudy(0)

                im.theta = 6./180.0*pi # 5 deg
                total_number_of_cases = 2 #12 #36

                # add equations
                    # DriveW_Freq = im.DriveW_Freq
                    # slip = slip_freq_breakdown_torque / DriveW_Freq
                    # im.DriveW_Freq = DriveW_Freq
                    # im.the_speed = DriveW_Freq * (1 - slip) * 30
                    # im.the_slip = slip
                study.GetDesignTable().AddEquation("freq")
                study.GetDesignTable().AddEquation("slip")
                study.GetDesignTable().AddEquation("speed")
                study.GetDesignTable().GetEquation("freq").SetType(0)
                study.GetDesignTable().GetEquation("freq").SetExpression("%g"%((im.DriveW_Freq)))
                study.GetDesignTable().GetEquation("freq").SetDescription("Excitation Frequency")
                study.GetDesignTable().GetEquation("slip").SetType(0)
                study.GetDesignTable().GetEquation("slip").SetExpression("%g"%(im.the_slip))
                study.GetDesignTable().GetEquation("slip").SetDescription("Slip [1]")
                study.GetDesignTable().GetEquation("speed").SetType(1)
                study.GetDesignTable().GetEquation("speed").SetExpression("freq * (1 - slip) * %d"%(60/(im.DriveW_poles/2)))
                study.GetDesignTable().GetEquation("speed").SetDescription("mechanical speed of four pole")

                # rotate the rotor by cad parameters via Park Transformation
                # cad parameters cannot be duplicated! even you have a list of cad paramters after duplicating, but they cannot be used to create cases! So you must set total_number_of_cases to 1 in the first place if you want to do Rotating Static FEA in JMAG
                im.add_cad_parameters(study)
                im.add_cases_rotate_rotor(study, total_number_of_cases) 
                    # print study.GetDesignTable().NumParameters()

                # set rotor current conditions
                im.slip_freq_breakdown_torque = slip_freq_breakdown_torque
                im.add_rotor_current_condition(app, model, study, total_number_of_cases, 
                                              self.dir_csv_output_folder + original_study_name + '_circuit_current.csv')
                print(self.dir_csv_output_folder + original_study_name + '_circuit_current.csv')
                print(self.dir_csv_output_folder + original_study_name + '_circuit_current.csv')
                print(self.dir_csv_output_folder + original_study_name + '_circuit_current.csv')
                    # print self.dir_csv_output_folder + im.get_individual_name() + original_study_name + '_circuit_current.csv'
                    # print self.dir_csv_output_folder + im.get_individual_name() + original_study_name + '_circuit_current.csv'
                    # print self.dir_csv_output_folder + im.get_individual_name() + original_study_name + '_circuit_current.csv'
                # set stator current conditions
                im.add_stator_current_condition(app, model, study, total_number_of_cases, 
                                              self.dir_csv_output_folder + original_study_name + '_circuit_current.csv')

                # https://www2.jmag-international.com/support/en/pdf/JMAG-Designer_Ver.17.1_ENv3.pdf
                study.GetStudyProperties().SetValue("DirectSolverType", 1)

                model.RestoreCadLink()

                study.RunAllCases()
                app.Save()
                model.CloseCadLink()


        # Stand-alone Loss Study
        pass



        # compute the fitness 
        rotor_volume = pi*(im.Radius_OuterRotor*1e-3)**2 * (im.stack_length*1e-3)
        rotor_weight = rotor_volume * 8050 # steel 8,050 kg/m3. Copper/Density 8.96 g/cm³
        cost_function = 30e3 / ( breakdown_torque/rotor_volume ) \
                        + 1.0 / ( breakdown_force/rotor_weight )
        logger = logging.getLogger(__name__)
        logger.debug('%s-%s: %g', self.project_name, im.model_name, cost_function)

        return cost_function

    def get_csv(self, data_name):
        r'D:\Users\horyc\OneDrive - UW-Madison\csv_opti\run#36\BLIM_PS_ID36-0-0Freq_#36-0-0_circuit_current.csv'
        pass
        pass

    def read_current_from_EC_FEA(self):

        # read from eddy current results
        dict_circuit_current_complex = {}
    
        with open(self.im.get_csv('circuit_current'), 'r') as f:
            for row in self.csv_row_reader(f):
                try: 
                    float(row[0])
                except:
                    continue
                else:
                    if '%g'%(self.slip_freq_breakdown_torque) in row[0]:
                        beginning_column = 1 + 2*3*2 # title + drive/bearing * 3 phase * real/imag
                        for i in range(0, int(self.no_slot_per_pole)):
                            natural_i = i+1
                            current_phase_column = beginning_column + i * int(self.im.DriveW_poles) * 2
                            for j in range(int(self.im.DriveW_poles)):
                                natural_j = j+1
                                re = float(row[current_phase_column+2*j])
                                im = float(row[current_phase_column+2*j+1])
                                dict_circuit_current_complex["%s%d"%(rotor_phase_name_list[i], natural_j)] = (re, im)
        dict_circuit_current_amp_and_phase = {}
        for key, item in dict_circuit_current_complex.items():
            amp = np.sqrt(item[1]**2 + item[0]**2)
            phase = np.arctan2(item[0], -item[1]) # atan2(y, x), y=a, x=-b
            dict_circuit_current_amp_and_phase[key] = (amp, phase)

    def show_results(self, femm_solver_data=None):
        print('show results!')
        # from pylab import *
        # plot style
        # plt.style.use('ggplot') 
        # plt.style.use('grayscale') # print plt.style.available # get [u'dark_background', u'bmh', u'grayscale', u'ggplot', u'fivethirtyeight']
        mpl.rcParams['legend.fontsize'] = 15
        font = {'family' : 'Times New Roman', #'serif',
                'weight' : 'normal',
                'size'   : 15}
        mpl.rcParams['font.family'] = ['serif'] # default is sans-serif
        mpl.rcParams['font.serif'] = ['Times New Roman']

        # color and alpha
        # Freq-FFVRC
        # rotor current

        fig, axes = subplots(2, 1, sharex=True)

        ''' TranRef '''
        study_name = 'TranRef'
        dm = utility.read_csv_results_4_comparison__transient(study_name, path_prefix=self.dir_csv_output_folder)
        basic_info, time_list, TorCon_list, ForConX_list, ForConY_list, ForConAbs_list = dm.unpack()
        end_time = time_list[-1]

        ax = axes[0]; ax.plot(time_list, TorCon_list, alpha=0.7, label=study_name); ax.set_xlabel('Time [s]'); ax.set_ylabel('Torque [Nm]')
        ax = axes[1]; ax.plot(time_list, ForConX_list, alpha=0.7, label=study_name+'-X'); ax.plot(time_list, ForConY_list, alpha=0.7, label=study_name+'Y'); ax.set_xlabel('Time [s]'); ax.set_ylabel('Force [N]')
        ax.plot(time_list, ForConAbs_list, alpha=0.7, label=study_name+'-Abs')


        # Avg_ForCon_Vector, _, Max_ForCon_Err_Angle = self.get_force_error_angle(ForConX_list[-400:], ForConY_list[-400:])
        # sfv = utilily.suspension_force_vector(ForConX_list, ForConY_list, range_ss=400)

        # print '---------------\nTranRef-FEA \nForce Average Vector:', Avg_ForCon_Vector, '[N]'
        # # print ForCon_Angle_List, 'deg'
        # print 'Maximum Force Angle Error', Max_ForCon_Err_Angle, '[deg]'
        # print '\tbasic info:', basic_info




        ''' Tran2TSS '''
        study_name = 'Tran2TSS'
        dm = utility.read_csv_results_4_comparison__transient(study_name, path_prefix=self.dir_csv_output_folder)
        basic_info, time_list, TorCon_list, ForConX_list, ForConY_list, ForConAbs_list = dm.unpack()

        ax = axes[0]; ax.plot(time_list, TorCon_list, alpha=0.7, label=study_name); ax.set_xlabel('Time [s]'); ax.set_ylabel('Torque [Nm]')
        ax = axes[1]; ax.plot(time_list, ForConX_list, alpha=0.7, label=study_name+'-X'); ax.plot(time_list, ForConY_list, alpha=0.7, label=study_name+'Y'); ax.set_xlabel('Time [s]'); ax.set_ylabel('Force [N]')
        ax.plot(time_list, ForConAbs_list, alpha=0.7, label=study_name+'-Abs')


        # Avg_ForCon_Vector, _, Max_ForCon_Err_Angle = self.get_force_error_angle(ForConX_list[-48:], ForConY_list[-48:])
        # sfv = utilily.suspension_force_vector(ForConX_list, ForConY_list, range_ss=400)
        # print '---------------\nTran2TSS-FEA \nForce Average Vector:', Avg_ForCon_Vector, '[N]'
        # # print ForCon_Angle_List, 'deg'
        # print 'Maximum Force Angle Error', Max_ForCon_Err_Angle, '[deg]'
        # print '\tbasic info:', basic_info




        ''' Static FEA with FEMM '''
        rotor_position_in_deg = femm_solver_data[0]*0.1 
        time_list = rotor_position_in_deg/180.*pi / self.im.Omega
        number_of_repeat = int(end_time / time_list[-1])
        # print number_of_repeat, end_time, time_list[-1]
        femm_force_x = femm_solver_data[2].tolist()
        femm_force_y = femm_solver_data[3].tolist()    
        femm_force_abs = np.sqrt(np.array(femm_force_x)**2 + np.array(femm_force_y)**2 )

        # # Vector plot
        # ax = figure().gca()
        # ax.text(0,0,'FEMM')
        # for x, y in zip(femm_force_x, femm_force_y):
        #     ax.arrow(0,0, x,y)
        # xlim([0,220])
        # ylim([0,220])

        # # # Force Error Angle
        # # Avg_ForCon_Vector, _, Max_ForCon_Err_Angle = self.get_force_error_angle(femm_force_x, femm_force_y)
        # # print '---------------\nFEMM-Static-FEA \nForce Average Vector:', Avg_ForCon_Vector, '[N]'
        # # # print ForCon_Angle_List, 'deg'
        # # print 'Maximum Force Angle Error', Max_ForCon_Err_Angle, '[deg]'

        # femm_torque  = number_of_repeat * femm_solver_data[1].tolist()
        # time_one_step = time_list[1]
        # time_list    = [i*time_one_step for i in range(len(femm_torque))]
        # femm_force_x = number_of_repeat * femm_solver_data[2].tolist()
        # femm_force_y = number_of_repeat * femm_solver_data[3].tolist()
        # femm_force_abs = number_of_repeat * femm_force_abs.tolist()
        # # print len(time_list), len(femm_force_abs)
        # if not femm_solver_data is None:
        #     ax = axes[0]; ax.plot(time_list, femm_torque,label='FEMM', zorder=1)
        #     ax = axes[1]; 
        #     ax.plot(time_list, femm_force_x,label='FEMM-X', zorder=1)
        #     ax.plot(time_list, femm_force_y,label='FEMM-Y', zorder=1)
        #     ax.plot(time_list, femm_force_abs,label='FEMM-Abs', zorder=1)





        axes[0].grid()
        axes[0].legend()
        axes[1].grid()
        axes[1].legend()
        # return basic_info

    # IEMDC19
    # show_results_iemdc19 is removed here
    def show_results_iemdc19(self, im_variant, femm_solver_data=None, femm_rotor_current_function=None):
        print('show results!')
        mpl.rcParams['font.family'] = ['serif'] # default is sans-serif
        mpl.rcParams['font.serif'] = ['Times New Roman']
    
        Fs = 500.*400.
        def basefreqFFT(x, Fs, base_freq=500, ax=None, ax_time_domain=None): #频域横坐标除以基频，即以基频为单位
            def nextpow2(L):
                n = 0
                while 2**n < L:
                    n += 1
                return n
            L = len(x)
            Ts = 1.0/Fs
            t = [el*Ts for el in range(0,L)]
            if ax_time_domain != None:
                ax_time_domain.plot(t, x)

            # NFFT = 2**nextpow2(L) # this causes incorrect dc bin (too large)
            NFFT = L
            y = np.fft.fft(x,NFFT) # y is a COMPLEX defined in numpy
            Y = [2 * el.__abs__() / L for el in y] # /L for spectrum aplitude consistent with actual signal. 2* for single-sided. abs for amplitude of complem number.
            Y[0] *= 0.5 # DC does not need to be times 2
            if base_freq==None:
                # f = np.fft.fftfreq(NFFT, t[1]-t[0]) # for double-sided
                f = Fs/2.0*np.linspace(0,1,NFFT/2+1) # unit is Hz
            else:
                f = Fs/2.0/base_freq*np.linspace(0,1,NFFT/2+1) # unit is base_freq Hz

            if ax == None:
                fig, ax = subplots()
            # ax.bar(f,Y[0:int(NFFT/2)+1], width=1.5)
            ax.plot(f,Y[0:int(NFFT/2)+1])

            # fig.title('Single-Sided Amplitude Spectrum of x(t)')
            # ax.xlabel('Frequency divided by base_freq / base freq * Hz')
            #ylabel('|Y(f)|')
            # ax.ylabel('Amplitude / 1')

            # # 计算频谱
            # fft_parameters = np.fft.fft(y_data) / len(y_data)
            # # 计算各个频率的振幅
            # fft_data = np.clip(20*np.log10(np.abs(fft_parameters))[:self.fftsize/2+1], -120, 120)
        def FFT_another_implementation(TorCon_list, Fs):
            # 这个FFT的结果和上面那个差不多，但是会偏小一点！不知道为什么！

            # Number of samplepoints
            N = len(TorCon_list)
            # Sample spacing
            T = 1.0 / Fs
            yf = np.fft.fft(TorCon_list)
            yf[0] *= 0.5
            xf = np.linspace(0.0, 1.0/(2.0*T), N/2)
            fig, ax = subplots()
            ax.plot(xf, 2.0/N * np.abs(yf[:N//2]))
        global count_plot
        count_plot = 0
        def add_plot(axeses, title=None, label=None, zorder=None, time_list=None, sfv=None, torque=None, range_ss=None, alpha=0.7):

            # # Avg_ForCon_Vector, Avg_ForCon_Magnitude, Avg_ForCon_Angle, ForCon_Angle_List, Max_ForCon_Err_Angle = self.get_force_error_angle(force_x[-range_ss:], force_y[-range_ss:])
            # print '\n\n---------------%s' % (title)
            # print 'Average Force Mag:', sfv.ss_avg_force_magnitude, '[N]'
            # print 'Average Torque:', sum(torque[-range_ss:])/len(torque[-range_ss:]), '[Nm]'
            # print 'Normalized Force Error Mag: %g%%, (+)%g%% (-)%g%%' % (0.5*(sfv.ss_max_force_err_abs[0]-sfv.ss_max_force_err_abs[1])/sfv.ss_avg_force_magnitude*100,
            #                                                               sfv.ss_max_force_err_abs[0]/sfv.ss_avg_force_magnitude*100,
            #                                                               sfv.ss_max_force_err_abs[1]/sfv.ss_avg_force_magnitude*100)
            # print 'Maximum Force Error Angle: %g [deg], (+)%g deg (-)%g deg' % (0.5*(sfv.ss_max_force_err_ang[0]-sfv.ss_max_force_err_ang[1]),
            #                                                              sfv.ss_max_force_err_ang[0],
            #                                                              sfv.ss_max_force_err_ang[1])
            # print 'Extra Information:'
            # print '\tAverage Force Vector:', sfv.ss_avg_force_vector, '[N]'
            # print '\tTorque Ripple (Peak-to-Peak)', max(torque[-range_ss:]) - min(torque[-range_ss:]), 'Nm'
            # print '\tForce Mag Ripple (Peak-to-Peak)', sfv.ss_max_force_err_abs[0] - sfv.ss_max_force_err_abs[1], 'N'

            ax = axeses[0][0]; ax.plot(time_list, torque, alpha=alpha, label=label, zorder=zorder)
            ax = axeses[0][1]; ax.plot(time_list, sfv.force_abs, alpha=alpha, label=label, zorder=zorder)
            ax = axeses[1][0]; ax.plot(time_list, 100*sfv.force_err_abs/sfv.ss_avg_force_magnitude, label=label, alpha=alpha, zorder=zorder)
            ax = axeses[1][1]; ax.plot(time_list, np.arctan2(sfv.force_y, sfv.force_x)/pi*180. - sfv.ss_avg_force_angle, label=label, alpha=alpha, zorder=zorder)

            global count_plot
            count_plot += 1
            # This is used for table in latex
            print('''
            \\newcommand\\torqueAvg%s{%s}
            \\newcommand\\torqueRipple%s{%s}
            \\newcommand\\forceAvg%s{%s}
            \\newcommand\\forceErrMag%s{%s}
            \\newcommand\\forceErrAng%s{%s}
            ''' % (chr(64+count_plot), utility.to_precision(sum(torque[-range_ss:])/len(torque[-range_ss:])),
                   chr(64+count_plot), utility.to_precision(0.5*(max(torque[-range_ss:]) - min(torque[-range_ss:]))),
                   chr(64+count_plot), utility.to_precision(sfv.ss_avg_force_magnitude),
                   chr(64+count_plot), utility.to_precision(0.5*(sfv.ss_max_force_err_abs[0]-sfv.ss_max_force_err_abs[1])/sfv.ss_avg_force_magnitude*100),
                   chr(64+count_plot), utility.to_precision(0.5*(sfv.ss_max_force_err_ang[0]-sfv.ss_max_force_err_ang[1]))))

        def iemdc_add_plot(axes, title=None, label=None, zorder=None, time_list=None, sfv=None, torque=None, range_ss=None, alpha=0.7, length_time_list=None, color='k'):
            if length_time_list is not None:
                ax = axes[0]; ax.plot(time_list[:length_time_list], torque[:length_time_list], alpha=alpha, label=label, zorder=zorder, color=color)
                ax = axes[1]; ax.plot(time_list[:length_time_list], sfv.force_abs[:length_time_list], alpha=alpha, label=label, zorder=zorder, color=color)
            else:
                ax = axes[0]; ax.plot(time_list, torque, alpha=alpha, label=label, zorder=zorder, color=color)
                ax = axes[1]; ax.plot(time_list, sfv.force_abs, alpha=alpha, label=label, zorder=zorder, color=color)

        fig_main, axeses = subplots(2, 2, sharex=False, dpi=150, figsize=(16, 8), facecolor='w', edgecolor='k')
        ax = axeses[0][0]; ax.set_xlabel('(a)',fontsize=14.5); ax.set_ylabel('Torque [Nm]',fontsize=14.5)
        ax = axeses[0][1]; ax.set_xlabel('(b)',fontsize=14.5); ax.set_ylabel('Force Amplitude [N]',fontsize=14.5)
        ax = axeses[1][0]; ax.set_xlabel('Time [s]\n(c)',fontsize=14.5); ax.set_ylabel('Normalized Force Error Magnitude [%]',fontsize=14.5)
        ax = axeses[1][1]; ax.set_xlabel('Time [s]\n(d)',fontsize=14.5); ax.set_ylabel('Force Error Angle [deg]',fontsize=14.5)
        for axes in axeses:
            for ax in axes:
                ax.set_xlim([0,0.1])

        iemdc_fig1, iemdc_axes1 = subplots(2, 1, sharex=False, dpi=150, figsize=(16, 8), facecolor='w', edgecolor='k')
        ax = iemdc_axes1[0]; ax.set_xlabel('Time [s]',fontsize=14.5); ax.set_ylabel('Torque [Nm]',fontsize=14.5)
        ax = iemdc_axes1[1]; ax.set_xlabel('Time [s]',fontsize=14.5); ax.set_ylabel('Force Amplitude [N]',fontsize=14.5)

        iemdc_fig2, iemdc_axes2 = subplots(2, 1, sharex=False, dpi=150, figsize=(16, 8), facecolor='w', edgecolor='k')
        ax = iemdc_axes2[0]; ax.set_xlabel('Time [s]',fontsize=14.5); ax.set_ylabel('Torque [Nm]',fontsize=14.5)
        ax = iemdc_axes2[1]; ax.set_xlabel('Time [s]',fontsize=14.5); ax.set_ylabel('Force Amplitude [N]',fontsize=14.5)

        iemdc_fig3, iemdc_axes3 = subplots(2, 1, sharex=False, dpi=150, figsize=(16, 8), facecolor='w', edgecolor='k')
        ax = iemdc_axes3[0]; ax.set_xlabel('Time [s]',fontsize=14.5); ax.set_ylabel('Torque [Nm]',fontsize=14.5)
        ax = iemdc_axes3[1]; ax.set_xlabel('Time [s]',fontsize=14.5); ax.set_ylabel('Force Amplitude [N]',fontsize=14.5)

        iemdc_fig4, iemdc_axes4 = subplots(2, 1, sharex=False, dpi=150, figsize=(16, 8), facecolor='w', edgecolor='k')
        ax = iemdc_axes4[0]; ax.set_xlabel('Time [s]',fontsize=14.5); ax.set_ylabel('Torque [Nm]',fontsize=14.5)
        ax = iemdc_axes4[1]; ax.set_xlabel('Time [s]',fontsize=14.5); ax.set_ylabel('Force Amplitude [N]',fontsize=14.5)

        for ax in iemdc_axes1:
            ax.grid(True)
        for ax in iemdc_axes2:
            ax.grid(True)
        for ax in iemdc_axes3:
            ax.grid(True)
        for ax in iemdc_axes4:
            ax.grid(True)

        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        # TranRef400
        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        study_name = 'TranRef'
        # dm = utility.read_csv_results_4_comparison__transient(study_name, path_prefix=r'D:\JMAG_Files\TimeStepSensitivity/'+'PS_Qr%d_NoEndRing_M15_17303l/'%(int(self.im.Qr)))
        dm = utility.read_csv_results_4_comparison__transient(study_name, path_prefix=self.dir_csv_output_folder)
        basic_info, time_list, TorCon_list, ForConX_list, ForConY_list, ForConAbs_list = dm.unpack()
        sfv = utility.suspension_force_vector(ForConX_list, ForConY_list, range_ss=400) # samples in the tail that are in steady state
        add_plot( axeses,
                  title=study_name,
                  label='Transient FEA', #'TranFEARef', #400',
                  zorder=1,
                  time_list=time_list,
                  sfv=sfv,
                  torque=TorCon_list,
                  range_ss=sfv.range_ss) 
        print('\tbasic info:', basic_info)

        iemdc_add_plot( iemdc_axes1,
                          title=study_name,
                          label='Transient FEA', #'TranFEARef', #400',
                          zorder=1,
                          time_list=time_list,
                          sfv=sfv,
                          torque=TorCon_list,
                          range_ss=sfv.range_ss,
                          length_time_list=int(len(time_list)/2),
                          color='#1f77b4') 

        # Current profile of TranRef
        fig_cur, axes_cur = subplots(2,1)
        fig_cur2, ax_cur2 = subplots(1,1)
        ax_cur = axes_cur[0]
        # for key in dm.key_list:
        #     if 'A1' in key: # e.g., ConductorA1
        #         ax_cur.plot(dm.Current_dict['Time(s)'], 
        #                 dm.Current_dict[key], 
        #                 label=study_name, #+'40',
        #                 alpha=0.7)
        # Current of TranRef400
        for key in dm.key_list:
            # if 'Coil' in key:
            #     ax_cur.plot(dm.Current_dict['Time(s)'], 
            #             dm.Current_dict[key], 
            #             label=key,
            #             alpha=0.7)
            if 'A1' in key: # e.g., ConductorA1
                ax_cur.plot(dm.Current_dict['Time(s)'], 
                        dm.Current_dict[key], 
                        label=study_name+'400',
                        alpha=0.7,
                        color='blue') #'purple')
                ax_cur2.plot(dm.Current_dict['Time(s)'], 
                        dm.Current_dict[key], 
                        label='Reference Transient FEA',
                        alpha=0.5,
                        color='k',
                        lw='0.5') #'purple')
                basefreqFFT(dm.Current_dict[key], Fs, ax=axes_cur[1])
                print('\tcheck time step', time_list[1], '=',  1./Fs)
                break

        # basefreqFFT(np.sin(2*pi*1500*np.arange(0,0.1,1./Fs)), Fs)
        fig, axes = subplots(2,1)
        basefreqFFT(TorCon_list[int(len(TorCon_list)/2):], Fs, base_freq=500., ax=axes[1], ax_time_domain=axes[0])
        fig, axes = subplots(2,1)
        basefreqFFT(ForConAbs_list[int(len(TorCon_list)/2):], Fs, ax=axes[1], ax_time_domain=axes[0])


        # ''' TranRef '''
        # study_name = 'TranRef'
        # dm = utility.read_csv_results_4_comparison__transient(study_name)    
        # basic_info, time_list, TorCon_list, ForConX_list, ForConY_list, ForConAbs_list = dm.unpack()
        # add_plot( axeses,
        #           title=study_name,
        #           label='TranFEARef40',
        #           zorder=5,
        #           time_list=time_list,
        #           force_x=ForConX_list,
        #           force_y=ForConY_list,
        #           force_abs=ForConAbs_list,
        #           torque=TorCon_list,
        #           range_ss=480) # samples in the tail that are in steady state
        # print '\tbasic info:', basic_info


        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        # Tran2TSS
        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        study_name = im_variant.get_individual_name() + 'Tran2TSS'
        dm = utility.read_csv_results_4_comparison__transient(study_name, self.dir_csv_output_folder)
        basic_info, time_list, TorCon_list, ForConX_list, ForConY_list, ForConAbs_list = dm.unpack()
        end_time = time_list[-1]
        sfv = utility.suspension_force_vector(ForConX_list, ForConY_list, range_ss=48) # samples in the tail that are in steady state
        add_plot( axeses,
                  title=study_name,
                  label='Transient FEA w/ 2 Time Step Sections',
                  zorder=8,
                  time_list=time_list,
                  sfv=sfv,
                  torque=TorCon_list,
                  range_ss=sfv.range_ss)
        print('\tbasic info:', basic_info)

        iemdc_add_plot( iemdc_axes2,
                          title=study_name,
                          label='Transient FEA w/ 2 Time Step Sections',
                          zorder=8,
                          time_list=time_list,
                          sfv=sfv,
                          torque=TorCon_list,
                          range_ss=sfv.range_ss,
                          length_time_list=None,
                          color='#ff7f0e') 

        # Current of Tran2TSS
        for key in dm.key_list:
            if 'A1' in key: # e.g., ConductorA1
                ax_cur.plot(dm.Current_dict['Time(s)'], 
                        dm.Current_dict[key], 
                        label=study_name,
                        alpha=0.7)
                ax_cur2.plot(dm.Current_dict['Time(s)'], 
                        dm.Current_dict[key], 
                        label='Transient w/ 2 Time Sect.',
                        alpha=0.7,
                        color='r')
                break

        # Current of FEMM
        if femm_rotor_current_function!=None:
            ax_cur.plot(dm.Current_dict['Time(s)'], 
                    [femm_rotor_current_function(t) for t in dm.Current_dict['Time(s)']], 
                    label='FEMM',
                    alpha=1,
                    c='r')
            ax_cur2.plot(dm.Current_dict['Time(s)'], 
                    [femm_rotor_current_function(t) for t in dm.Current_dict['Time(s)']], 
                    label='Static FEA',
                    alpha=1,
                    c='b')

        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        # Static FEA with FEMM
        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        study_name = 'FEMM'
        rotor_position_in_deg = femm_solver_data[0]*0.1 
        time_list = rotor_position_in_deg/180.*pi / self.im.Omega
        length_original_time_list = len(time_list)
        number_of_repeat = int(end_time / time_list[-1]) + 2
        femm_force_x = femm_solver_data[2].tolist()
        femm_force_y = femm_solver_data[3].tolist()    
        femm_force_abs = np.sqrt(np.array(femm_force_x)**2 + np.array(femm_force_y)**2 )

        # 延拓
        femm_torque  = number_of_repeat * femm_solver_data[1].tolist()
        time_one_step = time_list[1]
        time_list    = [i*time_one_step for i in range(len(femm_torque))]
        femm_force_x = number_of_repeat * femm_solver_data[2].tolist()
        femm_force_y = number_of_repeat * femm_solver_data[3].tolist()
        femm_force_abs = number_of_repeat * femm_force_abs.tolist()

        sfv = utility.suspension_force_vector(femm_force_x, femm_force_y, range_ss=len(rotor_position_in_deg)) # samples in the tail that are in steady state
        add_plot( axeses,
                  title=study_name,
                  label='Static FEA', #'StaticFEAwiRR',
                  zorder=3,
                  time_list=time_list,
                  sfv=sfv,
                  torque=femm_torque,
                  range_ss=sfv.range_ss,
                  alpha=0.5) 

        iemdc_add_plot( iemdc_axes3,
                        title=study_name,
                        label='Static FEA', #'StaticFEAwiRR',
                        zorder=3,
                        time_list=time_list,
                        sfv=sfv,
                        torque=femm_torque,
                        range_ss=sfv.range_ss,
                        alpha=0.5,
                        length_time_list=length_original_time_list,
                        color='#2ca02c') 
        
        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        # EddyCurrent
        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        study_name = im_variant.get_individual_name() + 'Freq-FFVRC'
        dm = utility.read_csv_results_4_comparison_eddycurrent(study_name, self.dir_csv_output_folder)
        _, _, TorCon_list, ForConX_list, ForConY_list, ForConAbs_list = dm.unpack()

        rotor_position_in_deg = 360./self.im.Qr / len(TorCon_list) * np.arange(0, len(TorCon_list))
        # print rotor_position_in_deg
        time_list = rotor_position_in_deg/180.*pi / self.im.Omega
        length_original_time_list = len(time_list)
        number_of_repeat = int(end_time / time_list[-1])

        # 延拓
        ec_torque        = number_of_repeat*TorCon_list
        time_one_step = time_list[1]
        time_list     = [i*time_one_step for i in range(len(ec_torque))]
        ec_force_abs     = number_of_repeat*ForConAbs_list.tolist()
        ec_force_x       = number_of_repeat*ForConX_list
        ec_force_y       = number_of_repeat*ForConY_list

        sfv = utility.suspension_force_vector(ec_force_x, ec_force_y, range_ss=len(rotor_position_in_deg))
        add_plot( axeses,
                  title=study_name,
                  label='Eddy Current FEA', #'EddyCurFEAwiRR',
                  zorder=2,
                  time_list=time_list,
                  sfv=sfv,
                  torque=ec_torque,
                  range_ss=sfv.range_ss) # samples in the tail that are in steady state

        iemdc_add_plot( iemdc_axes4,
                        title=study_name,
                        label='Eddy Current FEA', #'EddyCurFEAwiRR',
                        zorder=2,
                        time_list=time_list,
                        sfv=sfv,
                        torque=ec_torque,
                        range_ss=sfv.range_ss,  # samples in the tail that are in steady state
                        length_time_list=length_original_time_list,
                        color='#d62728') 


        # Force Vector plot
        ax = figure().gca()
        # ax.text(0,0,'FEMM')
        for x, y in zip(femm_force_x, femm_force_y):
            ax.arrow(0,0, x,y)
        # ax.set_xlim([0,220])
        # ax.set_ylim([0,220])
        ax.grid()
        ax.set_xlabel('Force X (FEMM) [N]'); 
        ax.set_ylabel('Force Y (FEMM) [N]')

        for ax in [axeses[0][0],axeses[0][1],axeses[1][0],axeses[1][1]]:
            ax.grid()
            ax.legend(loc='lower center')
            # ax.set_xlim([0,0.35335])


        # axeses[0][1].set_ylim([260, 335])
        # axeses[1][0].set_ylim([-0.06, 0.06])

        ax_cur.legend()
        ax_cur.grid()
        ax_cur.set_xlabel('Time [s]'); 
        ax_cur.set_ylabel('Rotor slot current [A]')
            # plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=12)

        ax_cur2.legend(fontsize=18)
        ax_cur2.grid()
        ax_cur2.set_xlabel('Time [s]', fontsize=18)
        ax_cur2.set_ylabel('Rotor slot current [A]', fontsize=18)
        ax_cur2.set_xlim([0, 0.09]); 
        ax_cur2.tick_params(axis='both', which='major', labelsize=18)
            # ax_cur2.set_xticklabels(ax_cur2.get_xticklabels(), fontsize=16)

        fig_main.tight_layout()
        if int(self.im.Qr) == 36:
            pass
            # fig_cur2.savefig(r'D:\OneDrive\[00]GetWorking\31 Bearingless_Induction_FEA_Model\p2019_iemdc_bearingless_induction full paper\images\Qr36_rotor_current_1000A.png', dpi=150)
            pass
            # fig_main.savefig('FEA_Model_Comparisons.png', dpi=150)
            #     # fig_main.savefig(r'D:\OneDrive\[00]GetWorking\31 BlessIMDesign\p2019_iemdc_bearingless_induction full paper\images\FEA_Model_Comparisons.png', dpi=150)
            # fig_main.savefig(r'D:\OneDrive\[00]GetWorking\31 Bearingless_Induction_FEA_Model\p2019_iemdc_bearingless_induction full paper\images\New_FEA_Model_Comparisons.png', dpi=150)


    def timeStepSensitivity(self):
        from pylab import figure, show, subplots, xlim, ylim
        print('\n\n\n-----------timeStepSensitivity------------')

        fig, axes = subplots(2, 1, sharex=True)

        ''' Super TranRef '''
        path_prefix = r'D:\JMAG_Files\TimeStepSensitivity/'
        dm = utility.read_csv_results_4_comparison__transient('TranRef', path_prefix=path_prefix+self.run_folder)
        basic_info, time_list, TorCon_list, ForConX_list, ForConY_list, ForConAbs_list = dm.unpack()

        study_name = 'SuperTranRef'
        ax = axes[0]; ax.plot(time_list, TorCon_list, alpha=0.7, label=study_name); ax.set_xlabel('Time [s]'); ax.set_ylabel('Torque [Nm]')
        ax = axes[1]; ax.plot(time_list, ForConX_list, alpha=0.7, label=study_name+'-X'); ax.plot(time_list, ForConY_list, alpha=0.7, label=study_name+'Y'); ax.set_xlabel('Time [s]'); ax.set_ylabel('Force [N]')
        ax.plot(time_list, ForConAbs_list, alpha=0.7, label=study_name+'-Abs')


        Avg_ForCon_Vector, _, Max_ForCon_Err_Angle = self.get_force_error_angle(ForConX_list[-48:], ForConY_list[-48:])
        print('---------------\nTranRef-FEA \nForce Average Vector:', Avg_ForCon_Vector, '[N]')
        # print ForCon_Angle_List, 'deg'
        print('Maximum Force Angle Error', Max_ForCon_Err_Angle, '[deg]')
        print('basic info:', basic_info)


        ''' TranRef '''
        study_name = 'TranRef'
        dm = utility.read_csv_results_4_comparison__transient(study_name)
        basic_info, time_list, TorCon_list, ForConX_list, ForConY_list, ForConAbs_list = dm.unpack()
        end_time = time_list[-1]

        ax = axes[0]; ax.plot(time_list, TorCon_list, alpha=0.7, label=study_name); ax.set_xlabel('Time [s]'); ax.set_ylabel('Torque [Nm]')
        ax = axes[1]; ax.plot(time_list, ForConX_list, alpha=0.7, label=study_name+'-X'); ax.plot(time_list, ForConY_list, alpha=0.7, label=study_name+'Y'); ax.set_xlabel('Time [s]'); ax.set_ylabel('Force [N]')
        ax.plot(time_list, ForConAbs_list, alpha=0.7, label=study_name+'-Abs')


        Avg_ForCon_Vector, _, Max_ForCon_Err_Angle = self.get_force_error_angle(ForConX_list[-48:], ForConY_list[-48:])
        print('---------------\nTranRef-FEA \nForce Average Vector:', Avg_ForCon_Vector, '[N]')
        # print ForCon_Angle_List, 'deg'
        print('Maximum Force Angle Error', Max_ForCon_Err_Angle, '[deg]')
        print('basic info:', basic_info)


        axes[0].grid()
        axes[0].legend()
        axes[1].grid()
        axes[1].legend()




    # ECCE19
    def read_csv_results_4_optimization(self, study_name, path_prefix=None):
        if path_prefix == None:
            path_prefix = self.dir_csv_output_folder
        # print 'look into:', path_prefix

        # Torque
        basic_info = []
        time_list = []
        TorCon_list = []
        with open(path_prefix + study_name + '_torque.csv', 'r') as f:
            count = 0
            for row in self.csv_row_reader(f):
                count +=1
                if count<=8:
                    try:
                        float(row[1])
                    except:
                        continue
                    else:
                        basic_info.append((row[0], float(row[1])))
                else:
                    time_list.append(float(row[0]))
                    TorCon_list.append(float(row[1]))

        # Force
        basic_info = []
        # time_list = []
        ForConX_list = []
        ForConY_list = []
        with open(path_prefix + study_name + '_force.csv', 'r') as f:
            count = 0
            for row in self.csv_row_reader(f):
                count +=1
                if count<=8:
                    try:
                        float(row[1])
                    except:
                        continue
                    else:
                        basic_info.append((row[0], float(row[1])))
                else:
                    # time_list.append(float(row[0]))
                    ForConX_list.append(float(row[1]))
                    ForConY_list.append(float(row[2]))
        ForConAbs_list = np.sqrt(np.array(ForConX_list)**2 + np.array(ForConY_list)**2 )

        # Current
        key_list = []
        Current_dict = {}
        with open(path_prefix + study_name + '_circuit_current.csv', 'r') as f:
            count = 0
            for row in self.csv_row_reader(f):
                count +=1
                if count<=8:
                    if 'Time' in row[0]: # Time(s)
                        for key in row:
                            key_list.append(key)
                            Current_dict[key] = []
                    else:
                        continue
                else:
                    for ind, val in enumerate(row):
                        Current_dict[key_list[ind]].append(float(val))

        # Terminal Voltage 
        new_key_list = []
        if self.fea_config_dict['delete_results_after_calculation'] == False:
            # file name is by individual_name like ID32-2-4_EXPORT_CIRCUIT_VOLTAGE.csv rather than ID32-2-4Tran2TSS_circuit_current.csv
            with open(path_prefix + study_name[:-8] + "_EXPORT_CIRCUIT_VOLTAGE.csv", 'r') as f:
                count = 0
                for row in self.csv_row_reader(f):
                    count +=1
                    if count==1: # Time | Terminal1 | Terminal2 | ... | Termial6
                        if 'Time' in row[0]: # Time, s
                            for key in row:
                                new_key_list.append(key) # Yes, you have to use a new key list, because the ind below bgeins at 0.
                                Current_dict[key] = []
                        else:
                            raise Exception('Problem with csv file for terminal voltage.')
                    else:
                        for ind, val in enumerate(row):
                            Current_dict[new_key_list[ind]].append(float(val))
        key_list += new_key_list

        # Loss
        # Iron Loss
        with open(path_prefix + study_name + '_iron_loss_loss.csv', 'r') as f:
            count = 0
            for row in self.csv_row_reader(f):
                count +=1
                if count>8:
                    stator_iron_loss = float(row[3]) # Stator Core
                    break
        with open(path_prefix + study_name + '_joule_loss_loss.csv', 'r') as f:
            count = 0
            for row in self.csv_row_reader(f):
                count +=1
                if count>8:
                    stator_eddycurrent_loss = float(row[3]) # Stator Core
                    break
        with open(path_prefix + study_name + '_hysteresis_loss_loss.csv', 'r') as f:
            count = 0
            for row in self.csv_row_reader(f):
                count +=1
                if count>8:
                    stator_hysteresis_loss = float(row[3]) # Stator Core
                    break
        # Copper Loss
        rotor_copper_loss_list = []
        with open(path_prefix + study_name + '_joule_loss.csv', 'r') as f:
            count = 0
            for row in self.csv_row_reader(f):
                count +=1
                if count>8:
                    if count==9:
                        stator_copper_loss = float(row[8]) # Coil # it is the same over time, this value does not account for end coil

                    rotor_copper_loss_list.append(float(row[7])) # Cage
    
        # use the last 1/4 period data to compute average copper loss of Tran2TSS rather than use that of Freq study
        effective_part = rotor_copper_loss_list[:int(0.5*self.fea_config_dict['number_of_steps_2ndTTS'])] # number_of_steps_2ndTTS = steps for half peirod
        rotor_copper_loss = sum(effective_part) / len(effective_part)

        if self.fea_config_dict['jmag_run_list'][0] == 0:
            utility.blockPrint()
            try:
                # convert rotor current results (complex number) into its amplitude
                self.femm_solver.list_rotor_current_amp = [abs(el) for el in self.femm_solver.vals_results_rotor_current] # el is complex number
                # settings not necessarily be consistent with Pyrhonen09's design: , STATOR_SLOT_FILL_FACTOR=0.5, ROTOR_SLOT_FILL_FACTOR=1., TEMPERATURE_OF_COIL=75
                s, r = self.femm_solver.get_copper_loss(self.femm_solver.stator_slot_area, self.femm_solver.rotor_slot_area)
            except Exception as e:
                raise e
            utility.enablePrint()
        else:
            s, r = None, None

        dm = utility.data_manager()
        dm.basic_info     = basic_info
        dm.time_list      = time_list
        dm.TorCon_list    = TorCon_list
        dm.ForConX_list   = ForConX_list
        dm.ForConY_list   = ForConY_list
        dm.ForConAbs_list = ForConAbs_list
        dm.Current_dict   = Current_dict
        dm.key_list       = key_list
        dm.jmag_loss_list    = [stator_copper_loss, rotor_copper_loss, stator_iron_loss, stator_eddycurrent_loss, stator_hysteresis_loss]
        dm.femm_loss_list = [s, r]
        return dm



    def duplicate_TranFEAwi2TSS_from_frequency_study(self, im_variant, slip_freq_breakdown_torque, app, model, original_study_name, tran2tss_study_name, logger):
        if model.NumStudies()<2:
            model.DuplicateStudyWithType(original_study_name, "Transient2D", tran2tss_study_name)
            app.SetCurrentStudy(tran2tss_study_name)
            study = app.GetCurrentStudy()
            self.study = study

            # 上一步的铁磁材料的状态作为下一步的初值，挺好，但是如果每一个转子的位置转过很大的话，反而会减慢非线性迭代。
            # 我们的情况是：0.33 sec 分成了32步，每步的时间大概在0.01秒，0.01秒乘以0.5*497 Hz = 2.485 revolution...
            # study.GetStudyProperties().SetValue(u"NonlinearSpeedup", 0) # JMAG17.1以后默认使用。现在后面密集的步长还多一点（32步），前面16步慢一点就慢一点呗！

            # two sections of different time step
            if True: # ECCE19
                number_of_steps_2ndTTS = self.fea_config_dict['number_of_steps_2ndTTS'] 
                DM = app.GetDataManager()
                DM.CreatePointArray("point_array/timevsdivision", "SectionStepTable")
                refarray = [[0 for i in range(3)] for j in range(3)]
                refarray[0][0] = 0
                refarray[0][1] =    1
                refarray[0][2] =        50
                refarray[1][0] = 0.5/slip_freq_breakdown_torque #0.5 for 17.1.03l # 1 for 17.1.02y
                refarray[1][1] =    number_of_steps_2ndTTS                          # 16 for 17.1.03l #32 for 17.1.02y
                refarray[1][2] =        50
                refarray[2][0] = refarray[1][0] + 0.5/im_variant.DriveW_Freq #0.5 for 17.1.03l 
                refarray[2][1] =    number_of_steps_2ndTTS  # also modify range_ss! # don't forget to modify below!
                refarray[2][2] =        50
                DM.GetDataSet("SectionStepTable").SetTable(refarray)
                number_of_total_steps = 1 + 2 * number_of_steps_2ndTTS # [Double Check] don't forget to modify here!
                study.GetStep().SetValue("Step", number_of_total_steps)
                study.GetStep().SetValue("StepType", 3)
                study.GetStep().SetTableProperty("Division", DM.GetDataSet("SectionStepTable"))

            else: # IEMDC19
                number_cycles_prolonged = 1 # 50
                DM = app.GetDataManager()
                DM.CreatePointArray("point_array/timevsdivision", "SectionStepTable")
                refarray = [[0 for i in range(3)] for j in range(4)]
                refarray[0][0] = 0
                refarray[0][1] =    1
                refarray[0][2] =        50
                refarray[1][0] = 1.0/slip_freq_breakdown_torque
                refarray[1][1] =    32 
                refarray[1][2] =        50
                refarray[2][0] = refarray[1][0] + 1.0/im_variant.DriveW_Freq
                refarray[2][1] =    48 # don't forget to modify below!
                refarray[2][2] =        50
                refarray[3][0] = refarray[2][0] + number_cycles_prolonged/im_variant.DriveW_Freq # =50*0.002 sec = 0.1 sec is needed to converge to TranRef
                refarray[3][1] =    number_cycles_prolonged*self.fea_config_dict['TranRef-StepPerCycle'] # =50*40, every 0.002 sec takes 40 steps 
                refarray[3][2] =        50
                DM.GetDataSet("SectionStepTable").SetTable(refarray)
                study.GetStep().SetValue("Step", 1 + 32 + 48 + number_cycles_prolonged*self.fea_config_dict['TranRef-StepPerCycle']) # [Double Check] don't forget to modify here!
                study.GetStep().SetValue("StepType", 3)
                study.GetStep().SetTableProperty("Division", DM.GetDataSet("SectionStepTable"))

            # add equations
            study.GetDesignTable().AddEquation("freq")
            study.GetDesignTable().AddEquation("slip")
            study.GetDesignTable().AddEquation("speed")
            study.GetDesignTable().GetEquation("freq").SetType(0)
            study.GetDesignTable().GetEquation("freq").SetExpression("%g"%((im_variant.DriveW_Freq)))
            study.GetDesignTable().GetEquation("freq").SetDescription("Excitation Frequency")
            study.GetDesignTable().GetEquation("slip").SetType(0)
            study.GetDesignTable().GetEquation("slip").SetExpression("%g"%(im_variant.the_slip))
            study.GetDesignTable().GetEquation("slip").SetDescription("Slip [1]")
            study.GetDesignTable().GetEquation("speed").SetType(1)
            study.GetDesignTable().GetEquation("speed").SetExpression("freq * (1 - slip) * %d"%(60/(im_variant.DriveW_poles/2)))
            study.GetDesignTable().GetEquation("speed").SetDescription("mechanical speed of four pole")

            # speed, freq, slip
            study.GetCondition("RotCon").SetValue("AngularVelocity", 'speed')
            app.ShowCircuitGrid(True)
            study.GetCircuit().GetComponent("CS4").SetValue("Frequency", "freq")
            study.GetCircuit().GetComponent("CS2").SetValue("Frequency", "freq")

            # max_nonlinear_iteration = 50
            # study.GetStudyProperties().SetValue(u"NonlinearMaxIteration", max_nonlinear_iteration)
            study.GetStudyProperties().SetValue("ApproximateTransientAnalysis", 1) # psuedo steady state freq is for PWM drive to use
            study.GetStudyProperties().SetValue("SpecifySlip", 1)
            study.GetStudyProperties().SetValue("OutputSteadyResultAs1stStep", 0)
            study.GetStudyProperties().SetValue("Slip", "slip")

            # # add other excitation frequencies other than 500 Hz as cases
            # for case_no, DriveW_Freq in enumerate([50.0, slip_freq_breakdown_torque]):
            #     slip = slip_freq_breakdown_torque / DriveW_Freq
            #     study.GetDesignTable().AddCase()
            #     study.GetDesignTable().SetValue(case_no+1, 0, DriveW_Freq)
            #     study.GetDesignTable().SetValue(case_no+1, 1, slip)

            # 你把Tran2TSS计算周期减半！
            # 也要在计算铁耗的时候选择1/4或1/2的数据！（建议1/4）
            # 然后，手动添加end step 和 start step，这样靠谱！2019-01-09：注意设置铁耗条件（iron loss condition）的Reference Start Step和End Step。

            # Iron Loss Calculation Condition
            # Stator 
            if True:
                cond = study.CreateCondition("Ironloss", "IronLossConStator")
                cond.SetValue("RevolutionSpeed", "freq*60/%d"%(0.5*(im_variant.DriveW_poles)))
                cond.ClearParts()
                sel = cond.GetSelection()
                sel.SelectPartByPosition(-im_variant.Radius_OuterStatorYoke+EPS, 0 ,0)
                cond.AddSelected(sel)
                # Use FFT for hysteresis to be consistent with FEMM's results and to have a FFT plot
                cond.SetValue("HysteresisLossCalcType", 1)
                cond.SetValue("PresetType", 3) # 3:Custom
                # Specify the reference steps yourself because you don't really know what JMAG is doing behind you
                cond.SetValue("StartReferenceStep", number_of_total_steps+1-number_of_steps_2ndTTS*0.5) # 1/4 period <=> number_of_steps_2ndTTS*0.5
                cond.SetValue("EndReferenceStep", number_of_total_steps)
                cond.SetValue("UseStartReferenceStep", 1)
                cond.SetValue("UseEndReferenceStep", 1)
                cond.SetValue("Cyclicity", 4) # specify reference steps for 1/4 period and extend it to whole period
                cond.SetValue("UseFrequencyOrder", 1)
                cond.SetValue("FrequencyOrder", "1-50") # Harmonics up to 50th orders 
            # Check CSV reults for iron loss (You cannot check this for Freq study) # CSV and save space
            study.GetStudyProperties().SetValue("CsvOutputPath", self.dir_csv_output_folder) # it's folder rather than file!
            study.GetStudyProperties().SetValue("CsvResultTypes", "Torque;Force;LineCurrent;TerminalVoltage;JouleLoss;TotalDisplacementAngle;JouleLoss_IronLoss;IronLoss_IronLoss;HysteresisLoss_IronLoss")
            study.GetStudyProperties().SetValue("DeleteResultFiles", self.fea_config_dict['delete_results_after_calculation'])
            # Terminal Voltage/Circuit Voltage: Check for outputing CSV results 
            study.GetCircuit().CreateTerminalLabel("Terminal4U", 8, -13)
            study.GetCircuit().CreateTerminalLabel("Terminal4V", 8, -11)
            study.GetCircuit().CreateTerminalLabel("Terminal4W", 8, -9)
            study.GetCircuit().CreateTerminalLabel("Terminal2U", 23, -13)
            study.GetCircuit().CreateTerminalLabel("Terminal2V", 23, -11)
            study.GetCircuit().CreateTerminalLabel("Terminal2W", 23, -9)
            # Export Stator Core's field results only for iron loss calculation (the csv file of iron loss will be clean with this setting)
                # study.GetMaterial(u"Rotor Core").SetValue(u"OutputResult", 0) # at least one part on the rotor should be output or else a warning "the jplot file does not contains displacement results when you try to calc. iron loss on the moving part." will pop up, even though I don't add iron loss condition on the rotor.
            # study.GetMeshControl().SetValue(u"AirRegionOutputResult", 0)
            study.GetMaterial("Shaft").SetValue("OutputResult", 0)
            study.GetMaterial("Cage").SetValue("OutputResult", 0)
            study.GetMaterial("Coil").SetValue("OutputResult", 0)
            # Rotor
            if True:
                cond = study.CreateCondition("Ironloss", "IronLossConRotor")
                cond.SetValue("BasicFrequencyType", 2)
                cond.SetValue("BasicFrequency", "freq")
                    # cond.SetValue(u"BasicFrequency", u"slip*freq") # this require the signal length to be at least 1/4 of slip period, that's too long!
                cond.ClearParts()
                sel = cond.GetSelection()
                sel.SelectPartByPosition(-im_variant.Radius_Shaft-EPS, 0 ,0)
                cond.AddSelected(sel)
                # Use FFT for hysteresis to be consistent with FEMM's results
                cond.SetValue("HysteresisLossCalcType", 1)
                cond.SetValue("PresetType", 3)
                # Specify the reference steps yourself because you don't really know what JMAG is doing behind you
                cond.SetValue("StartReferenceStep", number_of_total_steps+1-number_of_steps_2ndTTS*0.5) # 1/4 period <=> number_of_steps_2ndTTS*0.5
                cond.SetValue("EndReferenceStep", number_of_total_steps)
                cond.SetValue("UseStartReferenceStep", 1)
                cond.SetValue("UseEndReferenceStep", 1)
                cond.SetValue("Cyclicity", 4) # specify reference steps for 1/4 period and extend it to whole period
                cond.SetValue("UseFrequencyOrder", 1)
                cond.SetValue("FrequencyOrder", "1-50") # Harmonics up to 50th orders 

            # https://www2.jmag-international.com/support/en/pdf/JMAG-Designer_Ver.17.1_ENv3.pdf
            study.GetStudyProperties().SetValue("DirectSolverType", 1)

            # it is duplicated study, so no need to set up mesh 
            # run
            self.run_study(im_variant, app, study, clock_time())
        else:
            # the results exist already?
            return 

    def run_study(self, im_variant, app, study, toc):
        logger = logging.getLogger(__name__)
        if self.fea_config_dict['JMAG_Scheduler'] == False:
            print('Run jam.exe...')
            # if run_list[1] == True:
            study.RunAllCases()
            msg = 'Time spent on %s is %g s.'%(study.GetName() , clock_time() - toc)
            logger.debug(msg)
            print(msg)
        else:
            print('Submit to JMAG_Scheduler...')
            job = study.CreateJob()
            job.SetValue("Title", study.GetName())
            job.SetValue("Queued", True)
            job.Submit(False) # Fallse:CurrentCase, True:AllCases
            logger.debug('Submit %s to queue (Tran2TSS).'%(im_variant.individual_name))
            # wait and check
            # study.CheckForCaseResults()
        app.Save()
        # if the jcf file already exists, it pops a msg window
        # study.WriteAllSolidJcf(self.dir_jcf, im_variant.model_name+study.GetName()+'Solid', True) # True : Outputs cases that do not have results 
        # study.WriteAllMeshJcf(self.dir_jcf, im_variant.model_name+study.GetName()+'Mesh', True)

        # # run
        # if self.fea_config_dict['JMAG_Scheduler'] == False:
        #     study.RunAllCases()
        #     app.Save()
        # else:
        #     job = study.CreateJob()
        #     job.SetValue(u"Title", study.GetName())
        #     job.SetValue(u"Queued", True)
        #     job.Submit(True)
        #     logger.debug('Submit %s to queue (Freq).'%(im_variant.individual_name))
        #     # wait and check
        #     # study.CheckForCaseResults()

    def mesh_study(self, im_variant, app, model, study):

        # this `if' judgment is effective only if JMAG-DeleteResultFiles is False 
        # if not study.AnyCaseHasResult(): 
        # mesh
        im_variant.add_mesh(study, model)

        # Export Image
        app.View().ShowAllAirRegions()
        # app.View().ShowMeshGeometry() # 2nd btn
        app.View().ShowMesh() # 3rn btn
        app.View().Zoom(3)
        app.View().Pan(-im_variant.Radius_OuterRotor, 0)
        app.ExportImageWithSize(self.dir_run + model.GetName() + '.png', 2000, 2000)
        app.View().ShowModel() # 1st btn. close mesh view, and note that mesh data will be deleted if only ouput table results are selected.

class bearingless_induction_motor_design(object):

    def __init__(self, row=None, fea_config_dict=None, model_name_prefix='PS'):

        # introspection (settings that may differ for initial design and variant designs)
        self.bool_initial_design = True
        self.fea_config_dict = fea_config_dict
        self.slip_freq_breakdown_torque = None
        self.MODEL_ROTATE = False 
    
        #01 Model Name
        self.model_name_prefix = model_name_prefix # do include 'PS' here
        self.name = 'BLIM'

        # self.TORQUE_CURRENT_RATIO = 0.975

        #02 Pyrhonen Data
        if row != None:
            self.ID = str(row[0])
            self.Qs = row[1]
            self.Qr = row[2]
            self.Angle_StatorSlotSpan = 360. / self.Qs # in deg.
            self.Angle_RotorSlotSpan  = 360. / self.Qr # in deg.

            self.Radius_OuterStatorYoke = row[3]
            self.Radius_InnerStatorYoke = row[4]
            self.Length_AirGap          = row[5]
            self.Radius_OuterRotor      = row[6]
            self.Radius_Shaft           = row[7]

            self.Length_HeadNeckRotorSlot = row[8]
            self.Radius_of_RotorSlot      = row[9]
            self.Location_RotorBarCenter  = row[10]
            self.Width_RotorSlotOpen      = row[11]

            self.Radius_of_RotorSlot2     = row[12]
            self.Location_RotorBarCenter2 = row[13]

            self.Angle_StatorSlotOpen            = row[14]
            self.Width_StatorTeethBody           = row[15]
            self.Width_StatorTeethHeadThickness  = row[16]
            self.Width_StatorTeethNeck           = row[17]

            self.DriveW_poles       = row[18]
            self.DriveW_zQ          = row[19] # per slot
            self.DriveW_Rs          = row[20]
            self.DriveW_CurrentAmp  = row[21] * self.fea_config_dict['TORQUE_CURRENT_RATIO']
            print('IM template: self.DriveW_CurrentAmp = %g, while row[21] = %g, and TORQUE_CURRENT_RATIO = %g' % (self.DriveW_CurrentAmp, row[21], self.fea_config_dict['TORQUE_CURRENT_RATIO']))
            self.DriveW_Freq        = row[22]

            self.stack_length       = row[23]


            # inferred design parameters
            # self.Radius_InnerStator = self.Length_AirGap + self.Radius_OuterRotor
            try:
                self.parameters_for_imposing_constraints_among_design_parameters = row[24:]
            except IndexError as e:
                logger.error('The initial design file you provided is not for the puporse of optimization.', exc_info=True)
        else:
            # this is called from shitty_design re-producer.
            return None # __init__ is required to return None. You cannot (or at least shouldn't) return something else.

        #03 Mechanical Parameters
        self.update_mechanical_parameters(slip_freq=50.0) #, syn_freq=500.)

        #04 Material Condutivity Properties
        if self.fea_config_dict is not None:
            self.End_Ring_Resistance = fea_config_dict['End_Ring_Resistance']
            self.Bar_Conductivity = fea_config_dict['Bar_Conductivity']
        self.Copper_Loss = self.DriveW_CurrentAmp**2 / 2 * self.DriveW_Rs * 3
        # self.Resistance_per_Turn = 0.01 # TODO


        #05 Windings & Excitation
        if self.fea_config_dict is not None:
            self.wily = winding_layout.winding_layout(self.fea_config_dict['DPNV'], self.Qs, self.DriveW_poles/2)

        if self.DriveW_poles == 2:
            self.BeariW_poles = 4
            if self.DriveW_zQ % 2 != 0:
                print('zQ=', self.DriveW_zQ)
                raise Exception('This zQ does not suit for two layer winding.')
        elif self.DriveW_poles == 4:
            self.BeariW_poles = 2;
        else:
            raise Exception('Not implemented error.')
        self.BeariW_turns      = self.DriveW_zQ
        self.BeariW_Rs         = self.DriveW_Rs * self.BeariW_turns / self.DriveW_zQ
        # self.BeariW_CurrentAmp = 0.025 * self.DriveW_CurrentAmp/0.975 # extra 2.5% as bearing current
        self.BeariW_CurrentAmp = self.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * (self.DriveW_CurrentAmp / self.fea_config_dict['TORQUE_CURRENT_RATIO'])
        self.BeariW_Freq       = self.DriveW_Freq


        # 临时这么处理吧，好困了
        self.fill_factor = 0.5
        self.Js = 4e6
        # CurrentAmp_in_the_slot = self.coils.mm2_slot_area * self.fill_factor * self.Js*1e-6 * np.sqrt(2) #/2.2*2.8
        # CurrentAmp_per_conductor = CurrentAmp_in_the_slot / self.DriveW_zQ
        # CurrentAmp_per_phase = CurrentAmp_per_conductor * self.wily.number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
        # variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor
        # self.DriveW_CurrentAmp = self.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
        # self.BeariW_CurrentAmp = self.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp
        # print('---Variant CurrentAmp_in_the_slot =', CurrentAmp_in_the_slot)
        # print('---variant_DriveW_CurrentAmp = CurrentAmp_per_phase =', variant_DriveW_CurrentAmp)
        # print('---self.DriveW_CurrentAmp =', self.DriveW_CurrentAmp)
        # print('---self.BeariW_CurrentAmp =', self.BeariW_CurrentAmp)
        # print('---TORQUE_CURRENT_RATIO:', self.fea_config_dict['TORQUE_CURRENT_RATIO'])
        # print('---SUSPENSION_CURRENT_RATIO:', self.fea_config_dict['SUSPENSION_CURRENT_RATIO'])



        if self.fea_config_dict is not None:
            self.dict_coil_connection = {41:self.wily.l41, 42:self.wily.l42, 21:self.wily.l21, 22:self.wily.l22} # 这里的2和4等价于leftlayer和rightlayer。

        #06 Meshing & Solver Properties
        self.max_nonlinear_iteration = 50 # 30 for transient solve
        self.meshSize_Rotor = 1.8 #1.2 0.6 # mm


        #07: Some Checking
        # if abs(self.Location_RotorBarCenter2 - self.Location_RotorBarCenter) < 0.1*(self.Radius_of_RotorSlot + self.Radius_of_RotorSlot2):
        #     logger = logging.getLogger(__name__)
        #     logger.debug('Warning: There is no need to use a drop shape rotor, because the centers of the inner and outer circles are too close: %g, %g.' % (self.Location_RotorBarCenter, self.Location_RotorBarCenter2))
        #     self.use_drop_shape_rotor_bar = False
        
        #     # for VanGogh (FEMM) to work properly （如何不把两者设置成一样的，那么FEMM画出来的（很短的droop shape slot）和JMAG画出来的（Round shape圆形槽）则不一样。
        #     self.Location_RotorBarCenter2_backup = self.Location_RotorBarCenter2
        #     self.Location_RotorBarCenter2 = self.Location_RotorBarCenter 
        #     # use the larger slot radius
        #     if self.Radius_of_RotorSlot > self.Radius_of_RotorSlot2:
        #         self.Radius_of_RotorSlot2 = self.Radius_of_RotorSlot
        #     if self.Radius_of_RotorSlot2 > self.Radius_of_RotorSlot:
        #         self.Radius_of_RotorSlot = self.Radius_of_RotorSlot2

        # else:
        if self.fea_config_dict is not None:
            self.use_drop_shape_rotor_bar = self.fea_config_dict['use_drop_shape_rotor_bar'] # Since 5/23/2019
        else:
            print('In population.py: use_drop_shape_rotor_bar is None or not given. Set it to True.\n'*3)
            self.use_drop_shape_rotor_bar = True

        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        # For a low Qr value, there is a chance that the Location_RotorBarCenter2 is larger than Location_RotorBarCenter,
        # This means the round bar should be used instead of drop shape bar.
        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        if self.use_drop_shape_rotor_bar == True:
            # logger.debug('Location_RotorBarCenter:1,2: %g, %g.'%(self.Location_RotorBarCenter, self.Location_RotorBarCenter2))

            if self.Location_RotorBarCenter2 >= self.Location_RotorBarCenter \
                or abs(self.Location_RotorBarCenter2 - self.Location_RotorBarCenter)<=0.5: # To avoid small entity (Line.4). See 2019-06-06 @ 工作日志
                logger = logging.getLogger(__name__)
                logger.debug('Location_RotorBarCenter (%g) and Location_RotorBarCenter2 (%g) suggests to use round bar.'%(self.Location_RotorBarCenter, self.Location_RotorBarCenter2))
                self.use_drop_shape_rotor_bar = False

                # for VanGogh (FEMM) to work properly
                self.Location_RotorBarCenter2_backup = self.Location_RotorBarCenter2
                self.Location_RotorBarCenter2 = self.Location_RotorBarCenter 

                # ~~use the larger slot radius~~ (this is problematic)
                # use the smaller slot radius
                if self.Radius_of_RotorSlot > self.Radius_of_RotorSlot2:
                    self.Radius_of_RotorSlot = self.Radius_of_RotorSlot2
                if self.Radius_of_RotorSlot2 > self.Radius_of_RotorSlot:
                    self.Radius_of_RotorSlot2 = self.Radius_of_RotorSlot

        if abs(self.Qs-self.Qr)<1:
            print('Warning: Must not use a same Qs and Qr, to avoid synchronous torques created by slot harmonics. - (7.111)')

        self.no_slot_per_pole = self.Qr/self.DriveW_poles
        if self.no_slot_per_pole.is_integer() == False:
            print('This slot-pole combination will not be applied with pole-specific rotor winding')

        self.RSH = []
        for pm in [-1, +1]:
            for v in [-1, +1]:
                    self.RSH.append( (pm * self.Qr*(1-self.the_slip)/(0.5*self.DriveW_poles) + v)*self.DriveW_Freq )
        # print self.Qr, ', '.join("%g" % (rsh/self.DriveW_Freq) for rsh in self.RSH), '\n'

    def update_mechanical_parameters(self, slip_freq=None, syn_freq=None):
        # This function is first introduced to derive the new slip for different fundamental frequencies.
        if syn_freq is None:
            syn_freq = self.DriveW_Freq
        else:
            if syn_freq != self.DriveW_Freq:
                raise Exception('I do not recommend to modify synchronous speed at instance level. Go update the initial design.')

        if syn_freq == 0.0: # lock rotor
            self.the_slip = 0. # this does not actually make sense
            if slip_freq == None:
                self.DriveW_Freq = self.slip_freq_breakdown_torque
                self.BeariW_Freq = self.slip_freq_breakdown_torque
            else:
                self.DriveW_Freq = slip_freq
                self.BeariW_Freq = slip_freq
        else:
            if slip_freq != None:
                # change slip
                self.the_slip = slip_freq / syn_freq
                self.slip_freq_breakdown_torque = slip_freq
            else:
                # change syn_freq so update the slip
                self.the_slip = self.slip_freq_breakdown_torque / syn_freq

            self.DriveW_Freq = syn_freq
            self.BeariW_Freq = syn_freq

        self.the_speed = self.DriveW_Freq*60. / (0.5*self.DriveW_poles) * (1 - self.the_slip) # rpm

        self.Omega = + self.the_speed / 60. * 2*pi
        self.omega = None # This variable name is devil! you can't tell its electrical or mechanical! #+ self.DriveW_Freq * (1-self.the_slip) * 2*pi
        # self.the_speed = + self.the_speed

        if self.fea_config_dict is not None:
            if self.fea_config_dict['flag_optimization'] == False: # or else it becomes annoying
                print('[Update ID:%s]'%(self.ID), self.slip_freq_breakdown_torque, self.the_slip, self.the_speed, self.Omega, self.DriveW_Freq, self.BeariW_Freq)

    @staticmethod
    def get_stator_yoke_diameter_Dsyi(stator_tooth_width_b_ds, area_stator_slot_Sus, stator_inner_radius_r_is, Qs, Width_StatorTeethHeadThickness, Width_StatorTeethNeck):
        stator_inner_radius_r_is_eff = stator_inner_radius_r_is + ( Width_StatorTeethHeadThickness + Width_StatorTeethNeck )*1e-3
        temp = (2*pi*stator_inner_radius_r_is_eff - Qs*stator_tooth_width_b_ds)
        stator_tooth_height_h_ds = ( np.sqrt(temp**2 + 4*pi*area_stator_slot_Sus*Qs) - temp ) / (2*pi)
        stator_yoke_diameter_Dsyi = 2*stator_inner_radius_r_is + 2*stator_tooth_height_h_ds
        return stator_yoke_diameter_Dsyi

    @staticmethod
    def get_rotor_tooth_height_h_dr(rotor_tooth_width_b_dr, area_rotor_slot_Sur, rotor_outer_radius_r_or, Qr, Length_HeadNeckRotorSlot, minimum__area_rotor_slot_Sur):
        rotor_outer_radius_r_or_eff = rotor_outer_radius_r_or - Length_HeadNeckRotorSlot*1e-3

        new__rotor_tooth_width_b_dr = rotor_tooth_width_b_dr
        new__area_rotor_slot_Sur = area_rotor_slot_Sur
        logger = logging.getLogger(__name__)
        thermal_penalty = 0.0
        while True:
            temp = (2*pi*rotor_outer_radius_r_or_eff - Qr*new__rotor_tooth_width_b_dr)
            # 注意，这里用的是numpy的sqrt函数，根号下为负号不会像math.sqrt那样raise Exception，而是返回一个nan。
            operand_in_sqrt = temp**2 - 4*pi*new__area_rotor_slot_Sur*Qr
            if operand_in_sqrt < 0: # if np.isnan(rotor_tooth_height_h_dr) == True:
                # 这里应该能自动选择一个最大的可行转子槽才行！
                # modify rotor geometry to make this work and loop back.

                new__area_rotor_slot_Sur -= area_rotor_slot_Sur*0.05 # decrease by 5% every loop
                logger.warn('Sur=%g too large. Try new value=%g.'%(area_rotor_slot_Sur, new__area_rotor_slot_Sur))
                if new__area_rotor_slot_Sur < minimum__area_rotor_slot_Sur: # minimum__area_rotor_slot_Sur corresponds to 8 MA/m^2 current density
                    new__area_rotor_slot_Sur += area_rotor_slot_Sur*0.05 # don't reduce new__area_rotor_slot_Sur any further
                    new__rotor_tooth_width_b_dr -= rotor_tooth_width_b_dr*0.05 # instead, decrease new__rotor_tooth_width_b_dr
                    logger.warn('Reach minimum__area_rotor_slot_Sur. Bad bound on b_dr.\n\tIn other words, b_dr=%g too wide. Try new value=%g.'%(rotor_tooth_width_b_dr, new__rotor_tooth_width_b_dr))
                thermal_penalty += 0.1
                # raise Exception('There is not enough space for rotor slot or the required rotor current density will not be fulfilled.')
            else:
                rotor_tooth_height_h_dr = ( -np.sqrt(operand_in_sqrt) + temp ) / (2*pi)
                break
        return rotor_tooth_height_h_dr, new__rotor_tooth_width_b_dr, thermal_penalty, new__area_rotor_slot_Sur

    @classmethod
    def local_design_variant(cls, im, number_current_generation, individual_index, x_denorm):
        # Never assign anything to im, you can build your self after calling cls and assign stuff to self

        # unpack x_denorm
        air_gap_length_delta          = x_denorm[0]*1e-3 # m
        stator_tooth_width_b_ds       = x_denorm[1]*1e-3 # m
        rotor_tooth_width_b_dr        = x_denorm[2]*1e-3 # m

        Width_RotorSlotOpen           = x_denorm[4]      # mm, rotor slot opening
        Length_HeadNeckRotorSlot      = x_denorm[6]

        # Constranint #2
        # stator_tooth_width_b_ds imposes constraint on stator slot height
        Width_StatorTeethHeadThickness = x_denorm[5]
        Width_StatorTeethNeck = 0.5 * Width_StatorTeethHeadThickness

        area_stator_slot_Sus    = im.parameters_for_imposing_constraints_among_design_parameters[0]
        stator_inner_radius_r_is = im.Radius_OuterRotor*1e-3 + air_gap_length_delta 
        stator_yoke_diameter_Dsyi = cls.get_stator_yoke_diameter_Dsyi(  stator_tooth_width_b_ds, 
                                                                    area_stator_slot_Sus, 
                                                                    stator_inner_radius_r_is,
                                                                    im.Qs,
                                                                    Width_StatorTeethHeadThickness,
                                                                    Width_StatorTeethNeck)

        if 'VariableStatorSlotDepth' in im.fea_config_dict['which_filter']:
            backup = stator_yoke_diameter_Dsyi
            print('[V] The original stator_yoke_diameter_Dsyi is %g m'%(stator_yoke_diameter_Dsyi))
            # TODO: 用定子槽深变量，覆盖掉这边计算出来的 stator_yoke_diameter_Dsyi 即可。
            # 同时们还要修改 Radius_OuterStatorYoke，从而保证轭部深度不变！
            stator_tooth_height_h_ds = x_denorm[7]*1e-3
            stator_yoke_diameter_Dsyi = 2*stator_inner_radius_r_is + 2*stator_tooth_height_h_ds
            My_Radius_OuterStatorYoke = im.Radius_OuterStatorYoke + 0.5*(stator_yoke_diameter_Dsyi - backup)*1e3
            print('[V] The new stator_yoke_diameter_Dsyi is %g mm (the original one is %g mm)'%(1e3*stator_yoke_diameter_Dsyi, 1e3*backup))
            print('[V] The new Radius_OuterStatorYoke is %g mm (the original one is %g mm)' % (My_Radius_OuterStatorYoke, im.Radius_OuterStatorYoke))
        else:
            My_Radius_OuterStatorYoke = im.Radius_OuterStatorYoke

        # Constranint #3
        # rotor_tooth_width_b_dr imposes constraint on rotor slot height
        area_rotor_slot_Sur = im.parameters_for_imposing_constraints_among_design_parameters[1]
        rotor_outer_radius_r_or = im.Radius_OuterRotor*1e-3
        # overwrite rotor_tooth_width_b_dr if there is not enough space for rotor slot
        rotor_tooth_height_h_dr, rotor_tooth_width_b_dr, thermal_penalty, new__area_rotor_slot_Sur = cls.get_rotor_tooth_height_h_dr(  rotor_tooth_width_b_dr,
                                                                area_rotor_slot_Sur,
                                                                rotor_outer_radius_r_or,
                                                                im.Qr,
                                                                Length_HeadNeckRotorSlot,
                                                                im.parameters_for_imposing_constraints_among_design_parameters[2])
        rotor_slot_height_h_sr = rotor_tooth_height_h_dr


        # radius of outer rotor slot
        Radius_of_RotorSlot = 1e3 * (2*pi*(im.Radius_OuterRotor - Length_HeadNeckRotorSlot)*1e-3 - rotor_tooth_width_b_dr*im.Qr) / (2*im.Qr+2*pi)
        Location_RotorBarCenter = im.Radius_OuterRotor - Length_HeadNeckRotorSlot - Radius_of_RotorSlot

        # Constraint #1: Rotor slot opening cannot be larger than rotor slot width.
        punishment = 0.0
        if Width_RotorSlotOpen>2*Radius_of_RotorSlot:
            logger = logging.getLogger(__name__)
            logger.warn('Constraint #1: Rotor slot opening cannot be larger than rotor slot width. Gen#%04d. Individual index=%d.', number_current_generation, individual_index)
            # we will plot a model with b1 = rotor_tooth_width_b_dr instead, and apply a punishment for this model
            Width_RotorSlotOpen = 0.95 * 2*Radius_of_RotorSlot # 确保相交
            punishment = 0.0

        # new method for radius of inner rotor slot 2
        Radius_of_RotorSlot2 = 1e3 * (2*pi*(im.Radius_OuterRotor - Length_HeadNeckRotorSlot - rotor_slot_height_h_sr*1e3)*1e-3 - rotor_tooth_width_b_dr*im.Qr) / (2*im.Qr-2*pi)
        Location_RotorBarCenter2 = im.Radius_OuterRotor - Length_HeadNeckRotorSlot - rotor_slot_height_h_sr*1e3 + Radius_of_RotorSlot2 

        # translate design_parameters into row variable
        row_translated_from_design_paramters = \
            [   im.ID + '-' + str(number_current_generation) + '-' + str(individual_index), # the ID is str
                im.Qs,
                im.Qr,
                My_Radius_OuterStatorYoke, #im.Radius_OuterStatorYoke,
                0.5*stator_yoke_diameter_Dsyi * 1e3, # 定子内轭部处的半径由需要的定子槽面积和定子齿宽来决定。
                x_denorm[0],                # [1] # Length_AirGap
                im.Radius_OuterRotor,
                im.Radius_Shaft, 
                Length_HeadNeckRotorSlot,            # [6]
                Radius_of_RotorSlot,                 # inferred from [2]
                Location_RotorBarCenter, 
                Width_RotorSlotOpen,                 # [4]
                Radius_of_RotorSlot2,
                Location_RotorBarCenter2,
                x_denorm[3],                # [3] # Angle_StatorSlotOpen
                x_denorm[1],                # [1] # Width_StatorTeethBody
                Width_StatorTeethHeadThickness,      # [5]
                Width_StatorTeethNeck,
                im.DriveW_poles, 
                im.DriveW_zQ, # turns per slot
                im.DriveW_Rs,    
                im.DriveW_CurrentAmp * (1.0 - 0.4*im.fea_config_dict['mimic_separate_winding_with_DPNV_winding']),
                im.DriveW_Freq,
                im.stack_length
            ] + im.parameters_for_imposing_constraints_among_design_parameters
                # im.parameters_for_imposing_constraints_among_design_parameters[0], # area_stator_slot_Sus
                # im.parameters_for_imposing_constraints_among_design_parameters[1], # area_rotor_slot_Sur 
                # im.parameters_for_imposing_constraints_among_design_parameters[2]  # minimum__area_rotor_slot_Sur
        # initialze with the class's __init__ method
        self = cls( row_translated_from_design_paramters, 
                    im.fea_config_dict, 
                    im.model_name_prefix)
        self.number_current_generation = number_current_generation
        self.individual_index = individual_index

        self.design_parameters = x_denorm

        self.stator_yoke_diameter_Dsyi = stator_yoke_diameter_Dsyi # this is used for copper loss calculation in FEMM_Solver.py
        self.rotor_slot_height_h_sr = rotor_slot_height_h_sr       # this is used for copper loss calculation in FEMM_Solver.py

        # introspection (settings that may differ for initial design and variant designs)
        self.bool_initial_design = False # no, this is a variant design of the initial design
        self.fea_config_dict = im.fea_config_dict # important FEA configuration data
        self.slip_freq_breakdown_torque = None # initialize this slip freq for FEMM or JMAG
        self.MODEL_ROTATE = False # during optimization, rotate at model level is not needed.
        logger = logging.getLogger(__name__) 
        logger.info('im_variant ID %s is initialized.', self.ID)

        self.use_drop_shape_rotor_bar = self.fea_config_dict['use_drop_shape_rotor_bar'] # since 5/23/2019

        if self.use_drop_shape_rotor_bar == True:
            # logger.debug('Location_RotorBarCenter:1,2: %g, %g.'%(self.Location_RotorBarCenter, self.Location_RotorBarCenter2))

            if self.Location_RotorBarCenter2 >= self.Location_RotorBarCenter \
                or abs(self.Location_RotorBarCenter2 - self.Location_RotorBarCenter)<=0.5: # To avoid small entity (Line.4). See 2019-06-06 @ 工作日志
                logger = logging.getLogger(__name__)
                logger.debug('Location_RotorBarCenter (%g) and Location_RotorBarCenter2 (%g) suggests to use round bar.'%(self.Location_RotorBarCenter, self.Location_RotorBarCenter2))
                self.use_drop_shape_rotor_bar = False

                # for VanGogh (FEMM) to work properly
                self.Location_RotorBarCenter2_backup = self.Location_RotorBarCenter2
                self.Location_RotorBarCenter2 = self.Location_RotorBarCenter 

                # ~~use the larger slot radius~~ (this is problematic)
                # use the smaller slot radius
                if self.Radius_of_RotorSlot > self.Radius_of_RotorSlot2:
                    self.Radius_of_RotorSlot = self.Radius_of_RotorSlot2
                if self.Radius_of_RotorSlot2 > self.Radius_of_RotorSlot:
                    self.Radius_of_RotorSlot2 = self.Radius_of_RotorSlot

        self.CurrentAmp_per_phase = None # will be used in copper loss calculation
        self.slot_area_utilizing_ratio = self.fea_config_dict['SUSPENSION_CURRENT_RATIO'] + self.fea_config_dict['TORQUE_CURRENT_RATIO']
        print('self.slot_area_utilizing_ratio:', self.slot_area_utilizing_ratio)

        # thermal penalty for reduced rotor slot area (copper hot) and rotor tooth width (iron hot)
        # try:
            # thermal_penalty
        # except:
            # thermal_penalty = 0.0
        if thermal_penalty != 0:
            if self.fea_config_dict['local_sensitivity_analysis'] == True:
                run_folder = self.fea_config_dict['run_folder'][:-1] + 'lsa/'
            else:
                run_folder = self.fea_config_dict['run_folder']

            print('Warn: Themal penalty (%g) is not written...'%(thermal_penalty))
            # with open(self.fea_config_dict['dir_parent'] + 'pop/' + run_folder + 'thermal_penalty_individuals.txt', 'a') as f:
            #     f.write(self.get_individual_name() + ',%g,%g,%g\n'%(thermal_penalty, rotor_tooth_width_b_dr, new__area_rotor_slot_Sur))
        self.thermal_penalty = thermal_penalty
        return self

    @classmethod
    def reproduce_the_problematic_design(cls, path_to_shitty_design_file):
        self = cls()
        with open(path_to_shitty_design_file, 'r') as f:
            buf = f.readlines()
            while True:
                title = buf.pop(0) # pop out the first item
                if 'Bearingless' in title: # Have we removed the title line?
                    break
            while True:
                # the last line does not end in ', \n', so handle it first
                last_line = buf.pop()
                if len(last_line) > 1:
                    exec('self.' + last_line[8:])
                    break

        count_of_space = buf[0].find('A')
        print('count_of_space=', count_of_space)

        for ind, line in enumerate(buf):
            index_equal_symbol = line.find('=') 
            # the member variable that is a string is problematic
            if 'ID' in line[:index_equal_symbol] or 'name' in line[:index_equal_symbol]: # 在等号之前出现了ID或者name，说明这个变量很可能是字符串！
                index_comma_symbol = line.find(',')
                line = line[:index_equal_symbol+2] + '"' + line[index_equal_symbol+2:index_comma_symbol] + '"' + line[index_comma_symbol:]

            # the leading 8 char are space, while the ending 3 char are ', \n'
            exec('self.' + line[count_of_space:-3]) 
        return self


    def get_rotor_volume(self, stack_length=None):
        if stack_length is None:
            return np.pi*(self.Radius_OuterRotor*1e-3)**2 * (self.stack_length*1e-3)
        else:
            return np.pi*(self.Radius_OuterRotor*1e-3)**2 * (stack_length*1e-3)

    def get_rotor_weight(self, gravity=9.8, stack_length=None):
        material_density_rho = get_material_data()[0]
        if stack_length is None:
            return gravity * self.get_rotor_volume() * material_density_rho # steel 7860 or 8050 kg/m^3. Copper/Density 8.96 g/cm³. gravity: 9.8 N/kg
        else:
            return gravity * self.get_rotor_volume(stack_length=stack_length) * material_density_rho # steel 7860 or 8050 kg/m^3. Copper/Density 8.96 g/cm³. gravity: 9.8 N/kg

    def whole_row_reader(self, reader):
        for row in reader:
            yield row[:]

    def show(self, toString=False):
        attrs = list(vars(self).items())
        key_list = [el[0] for el in attrs]
        val_list = [el[1] for el in attrs]
        the_dict = dict(list(zip(key_list, val_list)))
        sorted_key = sorted(key_list, key=lambda item: (int(item.partition(' ')[0]) if item[0].isdigit() else float('inf'), item)) # this is also useful for string beginning with digiterations '15 Steel'.
        tuple_list = [(key, the_dict[key]) for key in sorted_key]
        if toString==False:
            print('- Bearingless Induction Motor Individual #%s\n\t' % (self.ID), end=' ')
            print(', \n\t'.join("%s = %s" % item for item in tuple_list))
            return ''
        else:
            return '\n- Bearingless Induction Motor Individual #%s\n\t' % (self.ID) + ', \n\t'.join("%s = %s" % item for item in tuple_list)

    def pre_process_structural(self, app, listKeyPoints):
        # pre-process : you can select part by coordinate!

        model = app.GetCurrentModel()
        part_ID_list = model.GetPartIDs()

        if len(part_ID_list) != int(1 + 1 + 1):
            msg = 'Number of Parts is unexpected.\n' + self.show(toString=True)
            # utility.send_notification(text=msg)
            raise Exception(msg)

        id_shaft = part_ID_list[0]
        id_rotorCore = part_ID_list[1]
        partIDRange_Cage = part_ID_list[2]
        # print(part_ID_list)

        # ''' Add Part to Set for later references '''
        # Not Working ???
            # def part_set(name, x, y):
            #     model.GetSetList().CreatePartSet(name)
            #     model.GetSetList().GetSet(name).SetMatcherType("Selection")
            #     model.GetSetList().GetSet(name).ClearParts()
            #     sel = model.GetSetList().GetSet(name).GetSelection()
            #     sel.SelectPartByPosition(x,y,0.) # z=0 for 2D
            #     model.GetSetList().GetSet(name).AddSelected(sel)
            #     print(x,y)
            #     print('---')

        def part_set_by_id(name, part_id):
            model.GetSetList().CreatePartSet(name)
            model.GetSetList().GetSet(name).SetMatcherType("Selection")
            model.GetSetList().GetSet(name).ClearParts()
            sel = model.GetSetList().GetSet(name).GetSelection()
            sel.SelectPart(part_id) # z=0 for 2D
            model.GetSetList().GetSet(name).AddSelected(sel)

        # Create Set for Shaft
        part_set_by_id("ShaftSet", id_shaft )

        # Create Set for Rotor Iron Core
        part_set_by_id("RotorCoreSet", id_rotorCore)

        # Create Set for Bar
        part_set_by_id("BarSet", partIDRange_Cage)

        # Create Set for Motion Region
        def part_list_set(name, list_xy, prefix=None):
            model.GetSetList().CreatePartSet(name)
            model.GetSetList().GetSet(name).SetMatcherType("Selection")
            model.GetSetList().GetSet(name).ClearParts()
            sel = model.GetSetList().GetSet(name).GetSelection() 
            for xy in list_xy:
                sel.SelectPartByPosition(xy[0],xy[1],0.0) # z=0 for 2D
                model.GetSetList().GetSet(name).AddSelected(sel)
        part_list_set('Motion_Region', [[-self.Radius_Shaft+EPS, 0.0],
                                        [-self.Radius_Shaft-EPS, 0.0],
                                        [-self.Location_RotorBarCenter, 0.0] ])
        # not working???
            # # Create Set for Rotor Iron Core
            # part_list_set("RotorCoreSet", [[-self.Radius_Shaft-EPS, 0.0]])

            # # Create Set for Bar
            # part_list_set("BarSet", [[-self.Location_RotorBarCenter, 0.0]])

            # Create Set for Shaft
            # part_list_set("ShaftSet", [[-EPS, 0.0]])

        def edge_set(name,x,y):
            model.GetSetList().CreateEdgeSet(name)
            model.GetSetList().GetSet(name).SetMatcherType(u"Selection")
            model.GetSetList().GetSet(name).ClearParts()
            sel = model.GetSetList().GetSet(name).GetSelection()
            sel.SelectEdgeByPosition(x,y,0)
            model.GetSetList().GetSet(name).AddSelected(sel)
        model.SetVisibility(u"Rotor Core", 0)
        edge_set(u"RotorCoreShaft_Slave",  -self.Radius_Shaft, 0.0)
        model.SetVisibility(u"Rotor Core", 1)
        model.SetVisibility(u"Shaft", 0)
        edge_set(u"RotorCoreShaft_Master",  -self.Radius_Shaft, 0.0)
        model.SetVisibility(u"Shaft", 1)

        def edge_list_set(name, list_xy, prefix=None):
            model.GetSetList().CreateEdgeSet(name)
            model.GetSetList().GetSet(name).SetMatcherType("Selection")
            model.GetSetList().GetSet(name).ClearParts()
            sel = model.GetSetList().GetSet(name).GetSelection() 
            for xy in list_xy:
                # print (xy)
                # print('---')
                sel.SelectEdgeByPosition(xy[0],xy[1],0) # z=0 for 2D
                model.GetSetList().GetSet(name).AddSelected(sel)
        Rotor_Sector_Angle = 2*pi/self.Qr*0.5
        edge_list_set(u"FixtureConstraints", [[0.5*self.Radius_Shaft*-cos(Rotor_Sector_Angle), 0.5*self.Radius_Shaft*+sin(Rotor_Sector_Angle)],
                                              [0.5*self.Radius_Shaft*-cos(Rotor_Sector_Angle), 0.5*self.Radius_Shaft*-sin(Rotor_Sector_Angle)],
                                              [0.5*self.Radius_OuterRotor*-cos(Rotor_Sector_Angle), 0.5*self.Radius_OuterRotor*+sin(Rotor_Sector_Angle)],
                                              [0.5*self.Radius_OuterRotor*-cos(Rotor_Sector_Angle), 0.5*self.Radius_OuterRotor*-sin(Rotor_Sector_Angle)]])
        model.SetVisibility(u"Rotor Core", 0)
        edge_list_set(u"BarRotorCore_Master",  [[-self.Location_RotorBarCenter-self.Radius_of_RotorSlot, 0.0],
                                                [-self.Location_RotorBarCenter2+self.Radius_of_RotorSlot2, 0.0],
                                                [listKeyPoints[5][0]-0.01*EPS, +listKeyPoints[5][1]], # 偏太远了就选不到这条Edge了哦，0.01太大了，0.0001差不多
                                                [listKeyPoints[5][0]-0.01*EPS, -listKeyPoints[5][1]],
                                                [listKeyPoints[6][0]-0.01*EPS, +listKeyPoints[6][1]],
                                                [listKeyPoints[6][0]-0.01*EPS, -listKeyPoints[6][1]] ])
        model.SetVisibility(u"Rotor Core", 1)
        model.SetVisibility(u"Cage", 0)
        edge_list_set(u"BarRotorCore_Slave",  [ [-self.Location_RotorBarCenter-self.Radius_of_RotorSlot, 0.0],
                                                [-self.Location_RotorBarCenter2+self.Radius_of_RotorSlot2, 0.0],
                                                [listKeyPoints[5][0]-0.01*EPS, +listKeyPoints[5][1]],
                                                [listKeyPoints[5][0]-0.01*EPS, -listKeyPoints[5][1]],
                                                [listKeyPoints[6][0]-0.01*EPS, +listKeyPoints[6][1]],
                                                [listKeyPoints[6][0]-0.01*EPS, -listKeyPoints[6][1]] ])

        model.SetVisibility(u"Cage", 1)

    def pre_process(self, app):
        # pre-process : you can select part by coordinate!
        ''' Group '''
        def group(name, id_list):
            model.GetGroupList().CreateGroup(name)
            for the_id in id_list:
                model.GetGroupList().AddPartToGroup(name, the_id)

        model = app.GetCurrentModel() # model = app.GetModel(u"IM_DEMO_1")
        part_ID_list = model.GetPartIDs()
        # print part_ID_list

        # view = app.View()
        # view.ClearSelect()
        # sel = view.GetCurrentSelection()
        # sel.SelectPart(123)
        # sel.SetBlockUpdateView(False)

        global partIDRange_Cage

        # if len(part_ID_list) != int(1 + 1 + 1 + self.Qr + self.Qs*2 + self.Qr + 1): the last +1 is for the air hug rotor
        # if len(part_ID_list) != int(1 + 1 + 1 + self.Qr + self.Qs*2 + self.Qr):
        if len(part_ID_list) != int(1 + 1 + 1 + self.Qr + self.Qs*2):
            msg = 'Number of Parts is unexpected. Should be %d but only %d.\n'%(1 + 1 + 1 + self.Qr + self.Qs*2, len(part_ID_list)) + self.show(toString=True)
            # utility.send_notification(text=msg)
            # return msg
            raise utility.ExceptionBadNumberOfParts(msg)

        id_shaft = part_ID_list[0]
        id_rotorCore = part_ID_list[1]
        partIDRange_Cage = part_ID_list[2 : 2+int(self.Qr)]
        id_statorCore = part_ID_list[3+int(self.Qr)]
        partIDRange_Coil = part_ID_list[3+int(self.Qr) : 3+int(self.Qr) + int(self.Qs*2)]
        # partIDRange_AirWithinRotorSlots = part_ID_list[3+int(self.Qr) + int(self.Qs*2) : 3+int(self.Qr) + int(self.Qs*2) + int(self.Qr)]

        # print part_ID_list
        # print partIDRange_Cage
        # print partIDRange_Coil
        # print partIDRange_AirWithinRotorSlots
        group("Cage", partIDRange_Cage) # 59-44 = 15 = self.Qr - 1
        group("Coil", partIDRange_Coil) # 107-60 = 47 = 48-1 = self.Qs*2 - 1
        # group(u"AirWithinRotorSlots", partIDRange_AirWithinRotorSlots) # 123-108 = 15 = self.Qr - 1


        ''' Add Part to Set for later references '''
        def part_set(name, x, y):
            model.GetSetList().CreatePartSet(name)
            model.GetSetList().GetSet(name).SetMatcherType("Selection")
            model.GetSetList().GetSet(name).ClearParts()
            sel = model.GetSetList().GetSet(name).GetSelection()
            # print x,y
            sel.SelectPartByPosition(x,y,0) # z=0 for 2D
            model.GetSetList().GetSet(name).AddSelected(sel)

        # def edge_set(name,x,y):
        #     model.GetSetList().CreateEdgeSet(name)
        #     model.GetSetList().GetSet(name).SetMatcherType(u"Selection")
        #     model.GetSetList().GetSet(name).ClearParts()
        #     sel = model.GetSetList().GetSet(name).GetSelection()
        #     sel.SelectEdgeByPosition(x,y,0) # sel.SelectEdge(741)
        #     model.GetSetList().GetSet(name).AddSelected(sel)
        # edge_set(u"AirGapCoast", 0, self.Radius_OuterRotor+0.5*self.Length_AirGap)

        # Create Set for Shaft
        part_set("ShaftSet", 0.0, 0.0)

        # Create Set for 4 poles Winding
        R = 0.5*(self.Radius_InnerStatorYoke  +  (self.Radius_OuterRotor+self.Width_StatorTeethHeadThickness+self.Width_StatorTeethNeck)) 
            # THETA = (0.5*(self.Angle_StatorSlotSpan) -  0.05*(self.Angle_StatorSlotSpan-self.Angle_StatorSlotOpen))/180.*pi
        THETA = (0.5*(self.Angle_StatorSlotSpan) -  0.05*(self.Angle_StatorSlotSpan-self.Width_StatorTeethBody))/180.*pi
        X = R*cos(THETA)
        Y = R*sin(THETA)
        # l41=[ 'C', 'C', 'A', 'A', 'B', 'B', 'C', 'C', 'A', 'A', 'B', 'B', 'C', 'C', 'A', 'A', 'B', 'B', 'C', 'C', 'A', 'A', 'B', 'B', ]
        # l42=[ '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', ]
        count = 0
        for UVW, UpDown in zip(self.wily.l41,self.wily.l42):
            count += 1 
            part_set("Coil4%s%s %d"%(UVW,UpDown,count), X, Y)

            THETA += self.Angle_StatorSlotSpan/180.*pi
            X = R*cos(THETA)
            Y = R*sin(THETA)

        # Create Set for 2 poles Winding
            # THETA = (0.5*(self.Angle_StatorSlotSpan) +  0.05*(self.Angle_StatorSlotSpan-self.Angle_StatorSlotOpen))/180.*pi
        THETA = (0.5*(self.Angle_StatorSlotSpan) +  0.05*(self.Angle_StatorSlotSpan-self.Width_StatorTeethBody))/180.*pi
        X = R*cos(THETA)
        Y = R*sin(THETA)
        # l21=[ 'A', 'A', 'B', 'B', 'B', 'B', 'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'C', 'C', 'C', 'C', 'A', 'A', ]
        # l22=[ '-', '-', '+', '+', '+', '+', '-', '-', '-', '-', '+', '+', '+', '+', '-', '-', '-', '-', '+', '+', '+', '+', '-', '-', ]
        count = 0
        for UVW, UpDown in zip(self.wily.l21,self.wily.l22):
            count += 1 
            part_set("Coil2%s%s %d"%(UVW,UpDown,count), X, Y)

            THETA += self.Angle_StatorSlotSpan/180.*pi
            X = R*cos(THETA)
            Y = R*sin(THETA)



        # Create Set for Bars and Air within rotor slots
        R = self.Location_RotorBarCenter
                                                                                # Another BUG:对于这种槽型，part_set在将AirWithin添加为set的时候，会错误选择到转子导条上！实际上，AirWithinRotorSlots对于JMAG来说完全是没有必要的！
        # R_airR = self.Radius_OuterRotor - 0.1*self.Length_HeadNeckRotorSlot # if the airWithin is too big, minus EPS is not enough anymore
        THETA = pi # it is very important (when Qr is odd) to begin the part set assignment from the first bar you plot.
        X = R*cos(THETA)
        Y = R*sin(THETA)
        list_xy_bars = []
        # list_xy_airWithinRotorSlot = []
        for ind in range(int(self.Qr)):
            natural_ind = ind + 1
            # print THETA / pi *180
            part_set("Bar %d"%(natural_ind), X, Y)
            list_xy_bars.append([X,Y])
            # # # part_set(u"AirWithin %d"%(natural_ind), R_airR*cos(THETA),R_airR*sin(THETA))
            # list_xy_airWithinRotorSlot.append([R_airR*cos(THETA),R_airR*sin(THETA)])

            THETA += self.Angle_RotorSlotSpan/180.*pi
            X = R*cos(THETA)
            Y = R*sin(THETA)

        # Create Set for Motion Region
        def part_list_set(name, list_xy, prefix=None):
            model.GetSetList().CreatePartSet(name)
            model.GetSetList().GetSet(name).SetMatcherType("Selection")
            model.GetSetList().GetSet(name).ClearParts()
            sel = model.GetSetList().GetSet(name).GetSelection() 
            for xy in list_xy:
                sel.SelectPartByPosition(xy[0],xy[1],0) # z=0 for 2D
                model.GetSetList().GetSet(name).AddSelected(sel)
        # part_list_set(u'Motion_Region', [[0,0],[0,self.Radius_Shaft+EPS]] + list_xy_bars + list_xy_airWithinRotorSlot) 
        part_list_set('Motion_Region', [[0,0],[0,self.Radius_Shaft+EPS]] + list_xy_bars) 

        # Create Set for Cage
        model.GetSetList().CreatePartSet("CageSet")
        model.GetSetList().GetSet("CageSet").SetMatcherType("MatchNames")
        model.GetSetList().GetSet("CageSet").SetParameter("style", "prefix")
        model.GetSetList().GetSet("CageSet").SetParameter("text", "Cage")
        model.GetSetList().GetSet("CageSet").Rebuild()

        return True

    def add_study(self, app, model, dir_csv_output_folder, choose_study_type='frequency'):

        self.choose_study_type = choose_study_type # Transient solves for different Qr, skewing angles and short pitches.

        # Time step initialization
        if self.choose_study_type=='transient':
            self.minimal_time_interval = 1 / 16. / self.DriveW_Freq
            self.end_time = 0.2  #20 / self.DriveW_Freq
            self.no_divisiton = int(self.end_time/self.minimal_time_interval)
            self.no_steps = int(self.no_divisiton + 1)
            # print self.minimal_time_interval, end_time, no_divisiton, no_steps
            # quit()
        elif self.choose_study_type=='frequency': # freq analysis
            self.table_freq_division_refarray = [[0 for i in range(3)] for j in range(4)]
            self.table_freq_division_refarray[0][0] = 2
            self.table_freq_division_refarray[0][1] =   1
            self.table_freq_division_refarray[0][2] =    self.max_nonlinear_iteration
            self.table_freq_division_refarray[1][0] = 10
            self.table_freq_division_refarray[1][1] =   32 # 8
            self.table_freq_division_refarray[1][2] =    self.max_nonlinear_iteration
            self.table_freq_division_refarray[2][0] = 16
            self.table_freq_division_refarray[2][1] =   2
            self.table_freq_division_refarray[2][2] =    self.max_nonlinear_iteration
            self.table_freq_division_refarray[3][0] = 24
            self.table_freq_division_refarray[3][1] =   2
            self.table_freq_division_refarray[3][2] =    self.max_nonlinear_iteration

            # self.table_freq_division_refarray = [[0 for i in range(3)] for j in range(2)]
            # self.table_freq_division_refarray[0][0] = 2
            # self.table_freq_division_refarray[0][1] =   1
            # self.table_freq_division_refarray[0][2] =    self.max_nonlinear_iteration
            # self.table_freq_division_refarray[1][0] = 18
            # self.table_freq_division_refarray[1][1] =   4 # for testing the script # 16 
            # self.table_freq_division_refarray[1][2] =    self.max_nonlinear_iteration

            self.no_steps = sum([el[1] for el in self.table_freq_division_refarray])
        elif self.choose_study_type=='static': # static analysis
            pass

        # Study property
        if self.choose_study_type == 'transient':
            study_name = self.get_individual_name() + "Tran"
            model.CreateStudy("Transient2D", study_name)
            app.SetCurrentStudy(study_name)
            study = model.GetStudy(study_name)

            study.GetStudyProperties().SetValue("ModelThickness", self.stack_length) # Stack Length
            study.GetStudyProperties().SetValue("ConversionType", 0)
            study.GetStudyProperties().SetValue("ApproximateTransientAnalysis", 1) # psuedo steady state freq is for PWM drive to use
            study.GetStudyProperties().SetValue("SpecifySlip", 1)
            study.GetStudyProperties().SetValue("Slip", self.the_slip)
            study.GetStudyProperties().SetValue("OutputSteadyResultAs1stStep", 0)
            study.GetStudyProperties().SetValue("NonlinearMaxIteration", self.max_nonlinear_iteration)
            study.GetStudyProperties().SetValue("CsvOutputPath", dir_csv_output_folder) # it's folder rather than file!
            # study.GetStudyProperties().SetValue(u"CsvResultTypes", u"Torque;Force;FEMCoilFlux;LineCurrent;ElectricPower;TerminalVoltage;JouleLoss;TotalDisplacementAngle")
            study.GetStudyProperties().SetValue("CsvResultTypes", "Torque;Force;LineCurrent;TerminalVoltage;JouleLoss;TotalDisplacementAngle")
            study.GetStudyProperties().SetValue("TimePeriodicType", 2) # This is for TP-EEC but is not effective
            study.GetStep().SetValue("StepType", 1)
            study.GetStep().SetValue("Step", self.no_steps)
            study.GetStep().SetValue("StepDivision", self.no_divisiton)
            study.GetStep().SetValue("EndPoint", self.end_time)
            # study.GetStep().SetValue(u"Step", 501)
            # study.GetStep().SetValue(u"StepDivision", 500)
            # study.GetStep().SetValue(u"EndPoint", 0.5)
            # app.View().SetCurrentCase(1)
        elif self.choose_study_type=='frequency': # freq analysis
            study_name = self.get_individual_name() + "Freq"
            model.CreateStudy("Frequency2D", study_name)
            app.SetCurrentStudy(study_name)
            study = model.GetStudy(study_name)

            # Misc
            study.GetStudyProperties().SetValue("NonlinearMaxIteration", self.max_nonlinear_iteration)
            study.GetStudyProperties().SetValue("ModelThickness", self.stack_length) # Stack Length
            study.GetStudyProperties().SetValue("ConversionType", 0)

            # CSV & Output
            study.GetStudyProperties().SetValue("CsvOutputPath", dir_csv_output_folder) # it's folder rather than file!
                # study.GetStudyProperties().SetValue(u"CsvResultTypes", u"Torque;Force;FEMCoilFlux;LineCurrent;ElectricPower;TerminalVoltage;JouleLoss;TotalDisplacementAngle")
            study.GetStudyProperties().SetValue("CsvResultTypes", "Torque;Force;LineCurrent;JouleLoss")
            study.GetStudyProperties().SetValue("DeleteResultFiles", self.fea_config_dict['delete_results_after_calculation'])

            # Time step
            DM = app.GetDataManager()
            DM.CreatePointArray("point_array/frequency_vs_division", "table_freq_division")
            # DM.GetDataSet(u"").SetName(u"table_freq_division")
            DM.GetDataSet("table_freq_division").SetTable(self.table_freq_division_refarray)
            study.GetStep().SetValue("Step", self.no_steps)
            study.GetStep().SetValue("StepType", 3)
            study.GetStep().SetTableProperty("Division", DM.GetDataSet("table_freq_division"))

            # This is exclusive for freq analysis
            study.GetStudyProperties().SetValue("BHCorrection", 1)
            # print 'BHCorrection for nonlinear time harmonic analysis is turned ON.'
        elif self.choose_study_type=='static': # static analysis
            study_name = "Static"
            model.CreateStudy("Static2D", study_name)
            app.SetCurrentStudy(study_name)
            study = model.GetStudy(study_name)

            study.GetStudyProperties().SetValue("NonlinearMaxIteration", self.max_nonlinear_iteration)
            study.GetStudyProperties().SetValue("ModelThickness", self.stack_length) # Stack Length
            study.GetStudyProperties().SetValue("ConversionType", 0)
            study.GetStudyProperties().SetValue("CsvOutputPath", dir_csv_output_folder) # it's folder rather than file!
                # study.GetStudyProperties().SetValue(u"CsvResultTypes", u"Torque;Force;FEMCoilFlux;LineCurrent;ElectricPower;TerminalVoltage;JouleLoss;TotalDisplacementAngle")
            study.GetStudyProperties().SetValue("CsvResultTypes", "Torque;Force;LineCurrent;JouleLoss")
            study.GetStudyProperties().SetValue("DeleteResultFiles", self.fea_config_dict['delete_results_after_calculation'])

        # Material
        self.add_material(study)

        # Conditions - Motion
        if self.choose_study_type == 'transient':
            study.CreateCondition("RotationMotion", "RotCon")
            # study.GetCondition(u"RotCon").SetXYZPoint(u"", 0, 0, 1) # megbox warning
            study.GetCondition("RotCon").SetValue("AngularVelocity", int(self.the_speed))
            study.GetCondition("RotCon").ClearParts()
            study.GetCondition("RotCon").AddSet(model.GetSetList().GetSet("Motion_Region"), 0)

            study.CreateCondition("Torque", "TorCon")
            # study.GetCondition(u"TorCon").SetXYZPoint(u"", 0, 0, 0) # megbox warning
            study.GetCondition("TorCon").SetValue("TargetType", 1)
            study.GetCondition("TorCon").SetLinkWithType("LinkedMotion", "RotCon")
            study.GetCondition("TorCon").ClearParts()

            study.CreateCondition("Force", "ForCon")
            study.GetCondition("ForCon").SetValue("TargetType", 1)
            study.GetCondition("ForCon").SetLinkWithType("LinkedMotion", "RotCon")
            study.GetCondition("ForCon").ClearParts()
        elif self.choose_study_type=='frequency': # freq analysis
            study.CreateCondition("FQRotationMotion", "RotCon")
            # study.GetCondition(u"RotCon").SetXYZPoint(u"", 0, 0, 0)
            study.GetCondition("RotCon").ClearParts()
            study.GetCondition("RotCon").AddSet(model.GetSetList().GetSet("Motion_Region"), 0)

            study.CreateCondition("Torque", "TorCon")
            study.GetCondition("TorCon").SetValue("TargetType", 1)
            study.GetCondition("TorCon").SetLinkWithType("LinkedMotion", "RotCon")
            study.GetCondition("TorCon").ClearParts()

            study.CreateCondition("Force", "ForCon")
            study.GetCondition("ForCon").SetValue("TargetType", 1)
            study.GetCondition("ForCon").SetLinkWithType("LinkedMotion", "RotCon")
            study.GetCondition("ForCon").ClearParts()
        elif self.choose_study_type=='static': # static analysis
            # duplicating study can fail if the im instance is destroyed.
            # model.DuplicateStudyWithType(original_study_name, u"Static2D", "Static")
            # study = app.GetCurrentStudy()

            study.CreateCondition("Torque", "TorCon")
            study.GetCondition("TorCon").SetValue("TargetType", 1)
            study.GetCondition("TorCon").ClearParts()
            study.GetCondition("TorCon").AddSet(model.GetSetList().GetSet("Motion_Region"), 0)

            study.CreateCondition("Force", "ForCon")
            study.GetCondition("ForCon").SetValue("TargetType", 1)
            study.GetCondition("ForCon").ClearParts()
            study.GetCondition("ForCon").AddSet(model.GetSetList().GetSet("Motion_Region"), 0)

            # 静态场不需要用到电路和FEM Coil/Conductor，这里设置完直接返回了
            # no mesh results are needed
            study.GetStudyProperties().SetValue("OnlyTableResults", self.fea_config_dict['OnlyTableResults'])
            study.GetStudyProperties().SetValue("Magnetization", 0)
            study.GetStudyProperties().SetValue("PermeanceFactor", 0)
            study.GetStudyProperties().SetValue("DifferentialPermeability", 0)
            study.GetStudyProperties().SetValue("LossDensity", 0)
            study.GetStudyProperties().SetValue("SurfaceForceDensity", 0)
            study.GetStudyProperties().SetValue("LorentzForceDensity", 0)
            study.GetStudyProperties().SetValue("Stress", 0)
            study.GetStudyProperties().SetValue("HysteresisLossDensity", 0)
            study.GetStudyProperties().SetValue("RestartFile", 0)
            study.GetStudyProperties().SetValue("JCGMonitor", 0)
            study.GetStudyProperties().SetValue("CoerciveForceNormal", 0)
            study.GetStudyProperties().SetValue("Temperature", 0)
            study.GetStudyProperties().SetValue("IronLossDensity", 0) # 我们要铁耗作为后处理，而不是和磁场同时求解。（[搜Iron Loss Formulas] Use one of the following methods to calculate iron loss and iron loss density generated in magnetic materials in JMAG. • Calculating Iron Loss Using Only the Magnetic Field Analysis Solver (page 6): It is a method to run magnetic field analysis considering the effect of iron loss. In this method, iron loss condition is not used. • Calculating Iron Loss Using the Iron Loss Analysis Solver (page 8): It is a method to run iron loss analysis using the results data of magnetic field analysis. It will be one of the following procedures. • Run magnetic field analysis study with iron loss condition • Run iron loss analysis study with reference to the result file of magnetic field analysis This chapter describes these two methods.）


            # Linear Solver
            if False:
                # sometime nonlinear iteration is reported to fail and recommend to increase the accerlation rate of ICCG solver
                study.GetStudyProperties().SetValue("IccgAccel", 1.2) 
                study.GetStudyProperties().SetValue("AutoAccel", 0)
            else:
                # https://www2.jmag-international.com/support/en/pdf/JMAG-Designer_Ver.17.1_ENv3.pdf
                study.GetStudyProperties().SetValue("DirectSolverType", 1)

            # too many threads will in turn make them compete with each other and slow down the solve. 2 is good enough for eddy current solve. 6~8 is enough for transient solve.
            study.GetStudyProperties().SetValue("UseMultiCPU", True)
            study.GetStudyProperties().SetValue("MultiCPU", 2) # this is effective for Transient Solver and 2 is enough!

            self.study_name = study_name
            return study


        # Conditions - FEM Coils & Conductors (i.e. stator/rotor winding)
        if choose_study_type == 'frequency':
            if self.bool_3PhaseCurrentSource == True:
                msg = 'Cannot use Composite type function for the CurrentSource in circuit of JMAG. So it needs more work, e.g., two more CurrentSources.'
                logging.getLogger(__name__).warn(msg)
            self.add_circuit(app, model, study, bool_3PhaseCurrentSource=True)
        elif choose_study_type == 'transient':
            self.add_circuit(app, model, study, bool_3PhaseCurrentSource=self.bool_3PhaseCurrentSource)


        # True: no mesh or field results are needed
        study.GetStudyProperties().SetValue("OnlyTableResults", self.fea_config_dict['OnlyTableResults'])

        # Linear Solver
        if False:
            # sometime nonlinear iteration is reported to fail and recommend to increase the accerlation rate of ICCG solver
            study.GetStudyProperties().SetValue("IccgAccel", 1.2) 
            study.GetStudyProperties().SetValue("AutoAccel", 0)
        else:
            # this can be said to be super fast over ICCG solver.
            # https://www2.jmag-international.com/support/en/pdf/JMAG-Designer_Ver.17.1_ENv3.pdf
            study.GetStudyProperties().SetValue("DirectSolverType", 1)

        # This SMP is effective only if there are tons of elements. e.g., over 100,000.
        # too many threads will in turn make them compete with each other and slow down the solve. 2 is good enough for eddy current solve. 6~8 is enough for transient solve.
        study.GetStudyProperties().SetValue("UseMultiCPU", True)
        study.GetStudyProperties().SetValue("MultiCPU", 2) 

        # # this is for the CAD parameters to rotate the rotor. the order matters for param_no to begin at 0.
        # if self.MODEL_ROTATE:
        #     self.add_cad_parameters(study)

        self.study_name = study_name
        return study

    def add_mesh(self, study, model):
        # this is for multi slide planes, which we will not be using
        refarray = [[0 for i in range(2)] for j in range(1)]
        refarray[0][0] = 3
        refarray[0][1] = 1
        study.GetMeshControl().GetTable("SlideTable2D").SetTable(refarray) 

        study.GetMeshControl().SetValue("MeshType", 1) # make sure this has been exe'd: study.GetCondition(u"RotCon").AddSet(model.GetSetList().GetSet(u"Motion_Region"), 0)
        study.GetMeshControl().SetValue("RadialDivision", 4) # for air region near which motion occurs
        study.GetMeshControl().SetValue("CircumferentialDivision", 720) #1440) # for air region near which motion occurs 这个数足够大，sliding mesh才准确。
        study.GetMeshControl().SetValue("AirRegionScale", 1.05) # [Model Length]: Specify a value within the following area. (1.05 <= value < 1000)
        study.GetMeshControl().SetValue("MeshSize", 4) # mm
        study.GetMeshControl().SetValue("AutoAirMeshSize", 0)
        study.GetMeshControl().SetValue("AirMeshSize", 4) # mm
        study.GetMeshControl().SetValue("Adaptive", 0)

        study.GetMeshControl().CreateCondition("RotationPeriodicMeshAutomatic", "autoRotMesh") # with this you can choose to set CircumferentialDivision automatically

        study.GetMeshControl().CreateCondition("Part", "CageMeshCtrl")
        study.GetMeshControl().GetCondition("CageMeshCtrl").SetValue("Size", self.meshSize_Rotor)
        study.GetMeshControl().GetCondition("CageMeshCtrl").ClearParts()
        study.GetMeshControl().GetCondition("CageMeshCtrl").AddSet(model.GetSetList().GetSet("CageSet"), 0)

        study.GetMeshControl().CreateCondition("Part", "ShaftMeshCtrl")
        study.GetMeshControl().GetCondition("ShaftMeshCtrl").SetValue("Size", 10) # 10 mm
        study.GetMeshControl().GetCondition("ShaftMeshCtrl").ClearParts()
        study.GetMeshControl().GetCondition("ShaftMeshCtrl").AddSet(model.GetSetList().GetSet("ShaftSet"), 0)

        def mesh_all_cases(study):
            numCase = study.GetDesignTable().NumCases()
            for case in range(0, numCase):
                study.SetCurrentCase(case)
                if study.HasMesh() == False:
                    study.CreateMesh()
                # if case == 0:
                #     app.View().ShowAllAirRegions()
                #     app.View().ShowMeshGeometry()
                #     app.View().ShowMesh()
        if self.MODEL_ROTATE:
            if self.total_number_of_cases>1: # just to make sure
                model.RestoreCadLink()
                study.ApplyAllCasesCadParameters()

        mesh_all_cases(study)

    def get_individual_name(self):
        if self.fea_config_dict['flag_optimization'] == True:
            return "ID%s" % (self.ID)
        else:
            return "%s_ID%s" % (self.model_name_prefix, self.ID)


    # TranFEAwi2TSS
    def add_material(self, study):
        if 'M19' in self.fea_config_dict['Steel']:
            study.SetMaterialByName("Stator Core", "M-19 Steel Gauge-29")
            study.GetMaterial("Stator Core").SetValue("Laminated", 1)
            study.GetMaterial("Stator Core").SetValue("LaminationFactor", 95)
                # study.GetMaterial(u"Stator Core").SetValue(u"UserConductivityValue", 1900000)

            study.SetMaterialByName("Rotor Core", "M-19 Steel Gauge-29")
            study.GetMaterial("Rotor Core").SetValue("Laminated", 1)
            study.GetMaterial("Rotor Core").SetValue("LaminationFactor", 95)

        elif 'M15' in self.fea_config_dict['Steel']:
            study.SetMaterialByName("Stator Core", "M-15 Steel")
            study.GetMaterial("Stator Core").SetValue("Laminated", 1)
            study.GetMaterial("Stator Core").SetValue("LaminationFactor", 98)

            study.SetMaterialByName("Rotor Core", "M-15 Steel")
            study.GetMaterial("Rotor Core").SetValue("Laminated", 1)
            study.GetMaterial("Rotor Core").SetValue("LaminationFactor", 98)

        elif self.fea_config_dict['Steel'] == 'Arnon5':
            study.SetMaterialByName("Stator Core", "Arnon5-final")
            study.GetMaterial("Stator Core").SetValue("Laminated", 1)
            study.GetMaterial("Stator Core").SetValue("LaminationFactor", 96)

            study.SetMaterialByName("Rotor Core", "Arnon5-final")
            study.GetMaterial("Rotor Core").SetValue("Laminated", 1)
            study.GetMaterial("Rotor Core").SetValue("LaminationFactor", 96)

        else:
            msg = 'Warning: default material is used: DCMagnetic Type/50A1000.'
            print(msg)
            logging.getLogger(__name__).warn(msg)
            study.SetMaterialByName("Stator Core", "DCMagnetic Type/50A1000")
            study.GetMaterial("Stator Core").SetValue("UserConductivityType", 1)
            study.SetMaterialByName("Rotor Core", "DCMagnetic Type/50A1000")
            study.GetMaterial("Rotor Core").SetValue("UserConductivityType", 1)

        study.SetMaterialByName("Coil", "Copper")
        study.GetMaterial("Coil").SetValue("UserConductivityType", 1)

        study.SetMaterialByName("Cage", "Aluminium")
        study.GetMaterial("Cage").SetValue("EddyCurrentCalculation", 1)
        study.GetMaterial("Cage").SetValue("UserConductivityType", 1)
        study.GetMaterial("Cage").SetValue("UserConductivityValue", self.Bar_Conductivity)

    def add_circuit(self, app, model, study, bool_3PhaseCurrentSource=True):
        # Circuit - Current Source
        app.ShowCircuitGrid(True)
        study.CreateCircuit()

        # 4 pole motor Qs=24 dpnv implemented by two layer winding (6 coils). In this case, drive winding has the same slot turns as bearing winding
        def circuit(poles,turns,Rs,ampD,ampB,freq,phase=0, CommutatingSequenceD=0, CommutatingSequenceB=0, x=10,y=10, bool_3PhaseCurrentSource=True):
            study.GetCircuit().CreateSubCircuit("Star Connection", "Star Connection %d"%(poles), x, y) # è¿™äº›æ•°å­—æŒ‡çš„æ˜¯gridçš„ä¸ªæ•°ï¼Œç¬¬å‡ è¡Œç¬¬å‡ åˆ—çš„æ ¼ç‚¹å¤„
            study.GetCircuit().GetSubCircuit("Star Connection %d"%(poles)).GetComponent("Coil1").SetValue("Turn", turns)
            study.GetCircuit().GetSubCircuit("Star Connection %d"%(poles)).GetComponent("Coil1").SetValue("Resistance", Rs)
            study.GetCircuit().GetSubCircuit("Star Connection %d"%(poles)).GetComponent("Coil2").SetValue("Turn", turns)
            study.GetCircuit().GetSubCircuit("Star Connection %d"%(poles)).GetComponent("Coil2").SetValue("Resistance", Rs)
            study.GetCircuit().GetSubCircuit("Star Connection %d"%(poles)).GetComponent("Coil3").SetValue("Turn", turns)
            study.GetCircuit().GetSubCircuit("Star Connection %d"%(poles)).GetComponent("Coil3").SetValue("Resistance", Rs)
            study.GetCircuit().GetSubCircuit("Star Connection %d"%(poles)).GetComponent("Coil1").SetName("CircuitCoil%dU"%(poles))
            study.GetCircuit().GetSubCircuit("Star Connection %d"%(poles)).GetComponent("Coil2").SetName("CircuitCoil%dV"%(poles))
            study.GetCircuit().GetSubCircuit("Star Connection %d"%(poles)).GetComponent("Coil3").SetName("CircuitCoil%dW"%(poles))

            if bool_3PhaseCurrentSource == True: # must use this for frequency analysis

                study.GetCircuit().CreateComponent("3PhaseCurrentSource", "CS%d"%(poles))
                study.GetCircuit().CreateInstance("CS%d"%(poles), x-4, y+1)
                study.GetCircuit().GetComponent("CS%d"%(poles)).SetValue("Amplitude", ampD+ampB)
                study.GetCircuit().GetComponent("CS%d"%(poles)).SetValue("Frequency", "freq") # this is not needed for freq analysis # "freq" is a variable
                study.GetCircuit().GetComponent("CS%d"%(poles)).SetValue("PhaseU", phase)
                # Commutating sequence is essencial for the direction of the field to be consistent with speed: UVW rather than UWV
                study.GetCircuit().GetComponent("CS%d"%(poles)).SetValue("CommutatingSequence", CommutatingSequenceD) 
            else: 
                I1 = "CS%d-1"%(poles)
                I2 = "CS%d-2"%(poles)
                I3 = "CS%d-3"%(poles)
                study.GetCircuit().CreateComponent("CurrentSource", I1)
                study.GetCircuit().CreateInstance(                   I1, x-4, y+3)
                study.GetCircuit().CreateComponent("CurrentSource", I2)
                study.GetCircuit().CreateInstance(                   I2, x-4, y+1)
                study.GetCircuit().CreateComponent("CurrentSource", I3)
                study.GetCircuit().CreateInstance(                   I3, x-4, y-1)

                phase_shift_drive = -120 if CommutatingSequenceD == 1 else 120
                phase_shift_beari = -120 if CommutatingSequenceB == 1 else 120

                func = app.FunctionFactory().Composite()        
                f1 = app.FunctionFactory().Sin(ampD, freq, 0*phase_shift_drive) # "freq" variable cannot be used here. So pay extra attension here when you create new case of a different freq.
                f2 = app.FunctionFactory().Sin(ampB, freq, 0*phase_shift_beari)
                func.AddFunction(f1)
                func.AddFunction(f2)
                study.GetCircuit().GetComponent(I1).SetFunction(func)

                func = app.FunctionFactory().Composite()        
                f1 = app.FunctionFactory().Sin(ampD, freq, 1*phase_shift_drive)
                f2 = app.FunctionFactory().Sin(ampB, freq, 1*phase_shift_beari)
                func.AddFunction(f1)
                func.AddFunction(f2)
                study.GetCircuit().GetComponent(I2).SetFunction(func)

                func = app.FunctionFactory().Composite()
                f1 = app.FunctionFactory().Sin(ampD, freq, 2*phase_shift_drive)
                f2 = app.FunctionFactory().Sin(ampB, freq, 2*phase_shift_beari)
                func.AddFunction(f1)
                func.AddFunction(f2)
                study.GetCircuit().GetComponent(I3).SetFunction(func)

            study.GetCircuit().CreateComponent("Ground", "Ground")
            study.GetCircuit().CreateInstance("Ground", x+2, y+1)
        # 这里电流幅值中的0.5因子源自DPNV导致的等于2的平行支路数。没有考虑到这一点，是否会对initial design的有效性产生影响？
        # 仔细看DPNV的接线，对于转矩逆变器，绕组的并联支路数为2，而对于悬浮逆变器，绕组的并联支路数为1。

        npb = self.wily.number_parallel_branch
        nwl = self.wily.no_winding_layer # number of windign layers 
        # if self.fea_config_dict['DPNV_separate_winding_implementation'] == True or self.fea_config_dict['DPNV'] == False:
        if self.fea_config_dict['DPNV'] == False:
            # either a separate winding or a DPNV winding implemented as a separate winding
            ampD =  0.5 * (self.DriveW_CurrentAmp/npb + self.BeariW_CurrentAmp) # 为了代码能被四极电机和二极电机通用，代入看看就知道啦。
            ampB = -0.5 * (self.DriveW_CurrentAmp/npb - self.BeariW_CurrentAmp) # 关于符号，注意下面的DriveW对应的circuit调用时的ampB前还有个负号！
            if bool_3PhaseCurrentSource != True:
                raise Exception('Logic Error Detected.')
        else:
            '[B]: DriveW_CurrentAmp is set.'
            # case: DPNV as an actual two layer winding
            ampD = self.DriveW_CurrentAmp/npb
            ampB = self.BeariW_CurrentAmp
            if bool_3PhaseCurrentSource != False:
                raise Exception('Logic Error Detected.')

        circuit(self.DriveW_poles,  self.DriveW_zQ/nwl, bool_3PhaseCurrentSource=bool_3PhaseCurrentSource,
            Rs=self.DriveW_Rs,ampD= ampD,
                              ampB=-ampB, freq=self.DriveW_Freq, phase=0,
                              CommutatingSequenceD=self.wily.CommutatingSequenceD,
                              CommutatingSequenceB=self.wily.CommutatingSequenceB)
        circuit(self.BeariW_poles,  self.BeariW_turns/nwl, bool_3PhaseCurrentSource=bool_3PhaseCurrentSource,
            Rs=self.BeariW_Rs,ampD= ampD,
                              ampB=+ampB, freq=self.BeariW_Freq, phase=0,
                              CommutatingSequenceD=self.wily.CommutatingSequenceD,
                              CommutatingSequenceB=self.wily.CommutatingSequenceB,x=25) # CS4 corresponds to uauc (conflict with following codes but it does not matter.)

        # Link FEM Coils to Coil Set     
        # if self.fea_config_dict['DPNV_separate_winding_implementation'] == True or self.fea_config_dict['DPNV'] == False:
        if self.fea_config_dict['DPNV'] == False:
            def link_FEMCoils_2_CoilSet(poles,l1,l2):
                # link between FEM Coil Condition and Circuit FEM Coil
                for UVW in ['U','V','W']:
                    which_phase = "%d%s-Phase"%(poles,UVW)
                    study.CreateCondition("FEMCoil", which_phase)
                    condition = study.GetCondition(which_phase)
                    condition.SetLink("CircuitCoil%d%s"%(poles,UVW))
                    condition.GetSubCondition("untitled").SetName("Coil Set 1")
                    condition.GetSubCondition("Coil Set 1").SetName("delete")
                count = 0
                dict_dir = {'+':1, '-':0, 'o':None}
                # select the part to assign the FEM Coil condition
                for UVW, UpDown in zip(l1,l2):
                    count += 1 
                    if dict_dir[UpDown] is None:
                        # print 'Skip', UVW, UpDown
                        continue
                    which_phase = "%d%s-Phase"%(poles,UVW)
                    condition = study.GetCondition(which_phase)
                    condition.CreateSubCondition("FEMCoilData", "Coil Set %d"%(count))
                    subcondition = condition.GetSubCondition("Coil Set %d"%(count))
                    subcondition.ClearParts()
                    subcondition.AddSet(model.GetSetList().GetSet("Coil%d%s%s %d"%(poles,UVW,UpDown,count)), 0)
                    subcondition.SetValue("Direction2D", dict_dir[UpDown])
                # clean up
                for UVW in ['U','V','W']:
                    which_phase = "%d%s-Phase"%(poles,UVW)
                    condition = study.GetCondition(which_phase)
                    condition.RemoveSubCondition("delete")
            link_FEMCoils_2_CoilSet(self.DriveW_poles, 
                                    self.dict_coil_connection[int(self.DriveW_poles*10+1)], # 40 for 4 poles, +1 for UVW, 
                                    self.dict_coil_connection[int(self.DriveW_poles*10+2)])                 # +2 for up or down, 
            link_FEMCoils_2_CoilSet(self.BeariW_poles, 
                                    self.dict_coil_connection[int(self.BeariW_poles*10+1)], # 20 for 2 poles, +1 for UVW, .
                                    self.dict_coil_connection[int(self.BeariW_poles*10+2)])                 # +2 for up or down,  这里的2和4等价于leftlayer和rightlayer。
        else:
            # 两个改变，一个是激励大小的改变（本来是200A 和 5A，现在是205A和195A），
            # 另一个绕组分组的改变，现在的A相是上层加下层为一相，以前是用俩单层绕组等效的。

            # Link FEM Coils to Coil Set as double layer short pitched winding
            # Create FEM Coil Condition
            # here we map circuit component `Coil2A' to FEM Coil Condition 'phaseAuauc
            # here we map circuit component `Coil4A' to FEM Coil Condition 'phaseAubud
            for suffix, poles in zip(['GroupAC', 'GroupBD'], [2, 4]): # 仍然需要考虑poles，是因为为Coil设置Set那里的代码还没有更新。这里的2和4等价于leftlayer和rightlayer。
                for UVW in ['U','V','W']:
                    study.CreateCondition("FEMCoil", 'phase'+UVW+suffix)
                    # link between FEM Coil Condition and Circuit FEM Coil
                    condition = study.GetCondition('phase'+UVW+suffix)
                    condition.SetLink("CircuitCoil%d%s"%(poles,UVW))
                    condition.GetSubCondition("untitled").SetName("delete")
            count = 0 # count indicates which slot the current rightlayer is in.
            index = 0
            dict_dir = {'+':1, '-':0}
            coil_pitch = self.wily.coil_pitch #self.dict_coil_connection[0]
            # select the part (via `Set') to assign the FEM Coil condition
            for UVW, UpDown in zip(self.wily.l_rightlayer1, self.wily.l_rightlayer2):

                count += 1 
                if self.wily.grouping_AC[index] == 1:
                    suffix = 'GroupAC'
                else:
                    suffix = 'GroupBD'
                condition = study.GetCondition('phase'+UVW+suffix)

                # right layer
                # print (count, "Coil Set %d"%(count), end=' ')
                condition.CreateSubCondition("FEMCoilData", "Coil Set Right %d"%(count))
                subcondition = condition.GetSubCondition("Coil Set Right %d"%(count))
                subcondition.ClearParts()
                subcondition.AddSet(model.GetSetList().GetSet("Coil%d%s%s %d"%(4,UVW,UpDown,count)), 0) # poles=4 means right layer, rather than actual poles
                subcondition.SetValue("Direction2D", dict_dir[UpDown])

                # left layer
                if coil_pitch<=0:
                    raise Exception('把永磁电机circuit部分的代码移植过来！')
                if count+coil_pitch <= self.Qs:
                    count_leftlayer = count+coil_pitch
                    index_leftlayer = index+coil_pitch
                else:
                    count_leftlayer = int(count+coil_pitch - self.Qs)
                    index_leftlayer = int(index+coil_pitch - self.Qs)
                # 右层导体的电流方向是正，那么与其串联的一个coil_pitch之处的左层导体就是负！不需要再检查l_leftlayer2了~
                if UpDown == '+': 
                    UpDown = '-'
                else:
                    UpDown = '+'
                # print (count_leftlayer, "Coil Set %d"%(count_leftlayer))
                condition.CreateSubCondition("FEMCoilData", "Coil Set Left %d"%(count_leftlayer))
                subcondition = condition.GetSubCondition("Coil Set Left %d"%(count_leftlayer))
                subcondition.ClearParts()
                subcondition.AddSet(model.GetSetList().GetSet("Coil%d%s%s %d"%(2,UVW,UpDown,count_leftlayer)), 0) # poles=2 means left layer, rather than actual poles
                subcondition.SetValue("Direction2D", dict_dir[UpDown])
                # print 'coil_pitch=', coil_pitch
                # print l_rightlayer1[index], UVW
                # print l_leftlayer1[index_leftlayer]
                # print l_rightlayer1
                # print l_leftlayer1
                index += 1
                # double check
                if self.wily.l_leftlayer1[index_leftlayer] != UVW:
                    raise Exception('Bug in winding diagram.')
            # clean up
            for suffix in ['GroupAC', 'GroupBD']:
                for UVW in ['U','V','W']:
                    condition = study.GetCondition('phase'+UVW+suffix)
                    condition.RemoveSubCondition("delete")
            # raise Exception('Test DPNV PE.')

        # Condition - Conductor (i.e. rotor winding)
        for ind in range(int(self.Qr)):
            natural_ind = ind + 1
            study.CreateCondition("FEMConductor", "CdctCon %d"%(natural_ind))
            study.GetCondition("CdctCon %d"%(natural_ind)).GetSubCondition("untitled").SetName("Conductor Set 1")
            study.GetCondition("CdctCon %d"%(natural_ind)).GetSubCondition("Conductor Set 1").ClearParts()
            study.GetCondition("CdctCon %d"%(natural_ind)).GetSubCondition("Conductor Set 1").AddSet(model.GetSetList().GetSet("Bar %d"%(natural_ind)), 0)

        # Condition - Conductor - Grouping
        study.CreateCondition("GroupFEMConductor", "CdctCon_Group")
        for ind in range(int(self.Qr)):
            natural_ind = ind + 1
            study.GetCondition("CdctCon_Group").AddSubCondition("CdctCon %d"%(natural_ind), ind)

        # Link Conductors to Circuit
        # if 'PS' in self.model_name_prefix: # Pole-Specific Rotor Winding
        if self.fea_config_dict['PoleSpecific'] == True:
            def place_conductor(x,y,name):
                study.GetCircuit().CreateComponent("FEMConductor", name)
                study.GetCircuit().CreateInstance(name, x, y)
            def place_resistor(x,y,name,end_ring_resistance):
                study.GetCircuit().CreateComponent("Resistor", name)
                study.GetCircuit().CreateInstance(name, x, y)
                study.GetCircuit().GetComponent(name).SetValue("Resistance", end_ring_resistance)

            rotor_phase_name_list = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            X = 40; Y = 40;
            if self.DriveW_poles == 2:
                for i in range(int(self.no_slot_per_pole)):
                    Y += -12
                    place_conductor(X,   Y, "Conductor%s1"%(rotor_phase_name_list[i]))
                    # place_conductor(X, Y-3, u"Conductor%s2"%(rotor_phase_name_list[i]))
                    # place_conductor(X, Y-6, u"Conductor%s3"%(rotor_phase_name_list[i]))
                    place_conductor(X, Y-9, "Conductor%s2"%(rotor_phase_name_list[i]))

                    if self.End_Ring_Resistance == 0: # setting a small value to End_Ring_Resistance is a bad idea (slow down the solver). Instead, don't model it
                        # no end ring resistors to behave like FEMM model
                        study.GetCircuit().CreateWire(X+2,   Y, X+2, Y-9)
                        # study.GetCircuit().CreateWire(X-2, Y-3, X-2, Y-6)
                        # study.GetCircuit().CreateWire(X+2, Y-6, X+2, Y-9)
                        study.GetCircuit().CreateInstance("Ground", X-5, Y-2)
                        study.GetCircuit().CreateWire(X-2,   Y, X-5, Y)
                        study.GetCircuit().CreateWire(X-5,   Y, X-2, Y-9)
                    else:
                        raise Exception('With end ring is not implemented.')
            else: # poles = 4
                for i in range(int(self.no_slot_per_pole)):
                    Y += -12
                    place_conductor(X,   Y, "Conductor%s1"%(rotor_phase_name_list[i]))
                    place_conductor(X, Y-3, "Conductor%s2"%(rotor_phase_name_list[i]))
                    place_conductor(X, Y-6, "Conductor%s3"%(rotor_phase_name_list[i]))
                    place_conductor(X, Y-9, "Conductor%s4"%(rotor_phase_name_list[i]))

                    if self.End_Ring_Resistance == 0: # setting a small value to End_Ring_Resistance is a bad idea (slow down the solver). Instead, don't model it
                        # no end ring resistors to behave like FEMM model
                        study.GetCircuit().CreateWire(X+2,   Y, X+2, Y-3)
                        study.GetCircuit().CreateWire(X-2, Y-3, X-2, Y-6)
                        study.GetCircuit().CreateWire(X+2, Y-6, X+2, Y-9)
                        study.GetCircuit().CreateInstance("Ground", X-5, Y-2)
                        study.GetCircuit().CreateWire(X-2,   Y, X-5, Y)
                        study.GetCircuit().CreateWire(X-5,   Y, X-2, Y-9)
                    else:
                        place_resistor(X+4,   Y, "R_%s1"%(rotor_phase_name_list[i]), self.End_Ring_Resistance)
                        place_resistor(X-4, Y-3, "R_%s2"%(rotor_phase_name_list[i]), self.End_Ring_Resistance)
                        place_resistor(X+4, Y-6, "R_%s3"%(rotor_phase_name_list[i]), self.End_Ring_Resistance)
                        place_resistor(X-4, Y-9, "R_%s4"%(rotor_phase_name_list[i]), self.End_Ring_Resistance)
        
                        study.GetCircuit().CreateWire(X+6,   Y, X+2, Y-3)
                        study.GetCircuit().CreateWire(X-6, Y-3, X-2, Y-6)
                        study.GetCircuit().CreateWire(X+6, Y-6, X+2, Y-9)
                        study.GetCircuit().CreateWire(X-6, Y-9, X-7, Y-9)
                        study.GetCircuit().CreateWire(X-2, Y, X-7, Y)
                        study.GetCircuit().CreateInstance("Ground", X-7, Y-2)
                            #study.GetCircuit().GetInstance(u"Ground", ini_ground_no+i).RotateTo(90)
                        study.GetCircuit().CreateWire(X-7, Y, X-6, Y-9)

            for i in range(0, int(self.no_slot_per_pole)):
                natural_i = i+1
                if self.DriveW_poles == 2:
                    study.GetCondition("CdctCon %d"%(natural_i)                      ).SetLink("Conductor%s1"%(rotor_phase_name_list[i]))
                    study.GetCondition("CdctCon %d"%(natural_i+self.no_slot_per_pole)).SetLink("Conductor%s2"%(rotor_phase_name_list[i]))
                elif self.DriveW_poles == 4:
                    study.GetCondition("CdctCon %d"%(natural_i)                        ).SetLink("Conductor%s1"%(rotor_phase_name_list[i]))
                    study.GetCondition("CdctCon %d"%(natural_i+  self.no_slot_per_pole)).SetLink("Conductor%s2"%(rotor_phase_name_list[i]))
                    study.GetCondition("CdctCon %d"%(natural_i+2*self.no_slot_per_pole)).SetLink("Conductor%s3"%(rotor_phase_name_list[i]))
                    study.GetCondition("CdctCon %d"%(natural_i+3*self.no_slot_per_pole)).SetLink("Conductor%s4"%(rotor_phase_name_list[i]))
        else: # Cage
            dyn_circuit = study.GetCircuit().CreateDynamicCircuit("Cage")
            dyn_circuit.SetValue("AntiPeriodic", False)
            dyn_circuit.SetValue("Bars", int(self.Qr))
            dyn_circuit.SetValue("EndringResistance", self.End_Ring_Resistance)
            dyn_circuit.SetValue("GroupCondition", True)
            dyn_circuit.SetValue("GroupName", "CdctCon_Group")
            dyn_circuit.SetValue("UseInductance", False)
            dyn_circuit.Submit("Cage1", 23, 2)
            study.GetCircuit().CreateInstance("Ground", 25, 1)

    def add_TranFEAwi2TSS_study(self, slip_freq_breakdown_torque, app, model, dir_csv_output_folder, tran2tss_study_name, logger):
        im_variant = self
        # logger.debug('Slip frequency: %g = ' % (self.the_slip))
        self.the_slip = slip_freq_breakdown_torque / self.DriveW_Freq
        # logger.debug('Slip frequency:    = %g???' % (self.the_slip))
        study_name = tran2tss_study_name

        model.CreateStudy("Transient2D", study_name)
        app.SetCurrentStudy(study_name)
        study = model.GetStudy(study_name)

        # SS-ATA
        study.GetStudyProperties().SetValue("ApproximateTransientAnalysis", 1) # psuedo steady state freq is for PWM drive to use
        study.GetStudyProperties().SetValue("SpecifySlip", 1)
        study.GetStudyProperties().SetValue("Slip", self.the_slip) # this will be overwritted later with "slip"
        study.GetStudyProperties().SetValue("OutputSteadyResultAs1stStep", 0)
        # study.GetStudyProperties().SetValue(u"TimePeriodicType", 2) # This is for TP-EEC but is not effective

        # misc
        study.GetStudyProperties().SetValue("ConversionType", 0)
        study.GetStudyProperties().SetValue("NonlinearMaxIteration", self.max_nonlinear_iteration)
        study.GetStudyProperties().SetValue("ModelThickness", self.stack_length) # Stack Length

        # Material
        self.add_material(study)

        # Conditions - Motion
        study.CreateCondition("RotationMotion", "RotCon") # study.GetCondition(u"RotCon").SetXYZPoint(u"", 0, 0, 1) # megbox warning
        study.GetCondition("RotCon").SetValue("AngularVelocity", int(self.the_speed))
        study.GetCondition("RotCon").ClearParts()
        study.GetCondition("RotCon").AddSet(model.GetSetList().GetSet("Motion_Region"), 0)

        study.CreateCondition("Torque", "TorCon") # study.GetCondition(u"TorCon").SetXYZPoint(u"", 0, 0, 0) # megbox warning
        study.GetCondition("TorCon").SetValue("TargetType", 1)
        study.GetCondition("TorCon").SetLinkWithType("LinkedMotion", "RotCon")
        study.GetCondition("TorCon").ClearParts()

        study.CreateCondition("Force", "ForCon")
        study.GetCondition("ForCon").SetValue("TargetType", 1)
        study.GetCondition("ForCon").SetLinkWithType("LinkedMotion", "RotCon")
        study.GetCondition("ForCon").ClearParts()


        # Conditions - FEM Coils & Conductors (i.e. stator/rotor winding)
        self.add_circuit(app, model, study, bool_3PhaseCurrentSource=self.wily.bool_3PhaseCurrentSource)


        # True: no mesh or field results are needed
        study.GetStudyProperties().SetValue("OnlyTableResults", self.fea_config_dict['OnlyTableResults'])

        # Linear Solver
        if False:
            # sometime nonlinear iteration is reported to fail and recommend to increase the accerlation rate of ICCG solver
            study.GetStudyProperties().SetValue("IccgAccel", 1.2) 
            study.GetStudyProperties().SetValue("AutoAccel", 0)
        else:
            # this can be said to be super fast over ICCG solver.
            # https://www2.jmag-international.com/support/en/pdf/JMAG-Designer_Ver.17.1_ENv3.pdf
            study.GetStudyProperties().SetValue("DirectSolverType", 1)

        if self.fea_config_dict['MultipleCPUs'] == True:
            # This SMP(shared memory process) is effective only if there are tons of elements. e.g., over 100,000.
            # too many threads will in turn make them compete with each other and slow down the solve. 2 is good enough for eddy current solve. 6~8 is enough for transient solve.
            study.GetStudyProperties().SetValue("UseMultiCPU", True)
            study.GetStudyProperties().SetValue("MultiCPU", 2) 

        # # this is for the CAD parameters to rotate the rotor. the order matters for param_no to begin at 0.
        # if self.MODEL_ROTATE:
        #     self.add_cad_parameters(study)


        # 上一步的铁磁材料的状态作为下一步的初值，挺好，但是如果每一个转子的位置转过很大的话，反而会减慢非线性迭代。
        # 我们的情况是：0.33 sec 分成了32步，每步的时间大概在0.01秒，0.01秒乘以0.5*497 Hz = 2.485 revolution...
        # study.GetStudyProperties().SetValue(u"NonlinearSpeedup", 0) # JMAG17.1以后默认使用。现在后面密集的步长还多一点（32步），前面16步慢一点就慢一点呗！


        # two sections of different time step
        if True: # ECCE19
            number_of_steps_2ndTTS = self.fea_config_dict['number_of_steps_2ndTTS'] 
            DM = app.GetDataManager()
            DM.CreatePointArray("point_array/timevsdivision", "SectionStepTable")
            refarray = [[0 for i in range(3)] for j in range(3)]
            refarray[0][0] = 0
            refarray[0][1] =    1
            refarray[0][2] =        50
            refarray[1][0] = 0.5/slip_freq_breakdown_torque #0.5 for 17.1.03l # 1 for 17.1.02y
            refarray[1][1] =    number_of_steps_2ndTTS                          # 16 for 17.1.03l #32 for 17.1.02y
            refarray[1][2] =        50
            refarray[2][0] = refarray[1][0] + 0.5/im_variant.DriveW_Freq #0.5 for 17.1.03l 
            refarray[2][1] =    number_of_steps_2ndTTS  # also modify range_ss! # don't forget to modify below!
            refarray[2][2] =        50
            DM.GetDataSet("SectionStepTable").SetTable(refarray)
            number_of_total_steps = 1 + 2 * number_of_steps_2ndTTS # [Double Check] don't forget to modify here!
            study.GetStep().SetValue("Step", number_of_total_steps)
            study.GetStep().SetValue("StepType", 3)
            study.GetStep().SetTableProperty("Division", DM.GetDataSet("SectionStepTable"))

        else: # IEMDC19
            number_cycles_prolonged = 1 # 50
            DM = app.GetDataManager()
            DM.CreatePointArray("point_array/timevsdivision", "SectionStepTable")
            refarray = [[0 for i in range(3)] for j in range(4)]
            refarray[0][0] = 0
            refarray[0][1] =    1
            refarray[0][2] =        50
            refarray[1][0] = 1.0/slip_freq_breakdown_torque
            refarray[1][1] =    32 
            refarray[1][2] =        50
            refarray[2][0] = refarray[1][0] + 1.0/im_variant.DriveW_Freq
            refarray[2][1] =    48 # don't forget to modify below!
            refarray[2][2] =        50
            refarray[3][0] = refarray[2][0] + number_cycles_prolonged/im_variant.DriveW_Freq # =50*0.002 sec = 0.1 sec is needed to converge to TranRef
            refarray[3][1] =    number_cycles_prolonged*self.fea_config_dict['TranRef-StepPerCycle'] # =50*40, every 0.002 sec takes 40 steps 
            refarray[3][2] =        50
            DM.GetDataSet("SectionStepTable").SetTable(refarray)
            study.GetStep().SetValue("Step", 1 + 32 + 48 + number_cycles_prolonged*self.fea_config_dict['TranRef-StepPerCycle']) # [Double Check] don't forget to modify here!
            study.GetStep().SetValue("StepType", 3)
            study.GetStep().SetTableProperty("Division", DM.GetDataSet("SectionStepTable"))

        # add equations
        study.GetDesignTable().AddEquation("freq")
        study.GetDesignTable().AddEquation("slip")
        study.GetDesignTable().AddEquation("speed")
        study.GetDesignTable().GetEquation("freq").SetType(0)
        study.GetDesignTable().GetEquation("freq").SetExpression("%g"%((im_variant.DriveW_Freq)))
        study.GetDesignTable().GetEquation("freq").SetDescription("Excitation Frequency")
        study.GetDesignTable().GetEquation("slip").SetType(0)
        study.GetDesignTable().GetEquation("slip").SetExpression("%g"%(im_variant.the_slip))
        study.GetDesignTable().GetEquation("slip").SetDescription("Slip [1]")
        study.GetDesignTable().GetEquation("speed").SetType(1)
        study.GetDesignTable().GetEquation("speed").SetExpression("freq * (1 - slip) * %d"%(60/(im_variant.DriveW_poles/2)))
        study.GetDesignTable().GetEquation("speed").SetDescription("mechanical speed of four pole")

        # speed, freq, slip
        study.GetCondition("RotCon").SetValue("AngularVelocity", 'speed')
        if self.fea_config_dict['DPNV']==False:
            app.ShowCircuitGrid(True)
            study.GetCircuit().GetComponent("CS4").SetValue("Frequency", "freq")
            study.GetCircuit().GetComponent("CS2").SetValue("Frequency", "freq")

        # max_nonlinear_iteration = 50
        # study.GetStudyProperties().SetValue(u"NonlinearMaxIteration", max_nonlinear_iteration)
        study.GetStudyProperties().SetValue("ApproximateTransientAnalysis", 1) # psuedo steady state freq is for PWM drive to use
        study.GetStudyProperties().SetValue("SpecifySlip", 1)
        study.GetStudyProperties().SetValue("OutputSteadyResultAs1stStep", 0)
        study.GetStudyProperties().SetValue("Slip", "slip") # overwrite with variables

        # # add other excitation frequencies other than 500 Hz as cases
        # for case_no, DriveW_Freq in enumerate([50.0, slip_freq_breakdown_torque]):
        #     slip = slip_freq_breakdown_torque / DriveW_Freq
        #     study.GetDesignTable().AddCase()
        #     study.GetDesignTable().SetValue(case_no+1, 0, DriveW_Freq)
        #     study.GetDesignTable().SetValue(case_no+1, 1, slip)

        # 你把Tran2TSS计算周期减半！
        # 也要在计算铁耗的时候选择1/4或1/2的数据！（建议1/4）
        # 然后，手动添加end step 和 start step，这样靠谱！2019-01-09：注意设置铁耗条件（iron loss condition）的Reference Start Step和End Step。

        # Iron Loss Calculation Condition
        # Stator 
        if True:
            cond = study.CreateCondition("Ironloss", "IronLossConStator")
            cond.SetValue("RevolutionSpeed", "freq*60/%d"%(0.5*(im_variant.DriveW_poles)))
            cond.ClearParts()
            sel = cond.GetSelection()
            sel.SelectPartByPosition(-im_variant.Radius_OuterStatorYoke+EPS, 0 ,0)
            cond.AddSelected(sel)
            # Use FFT for hysteresis to be consistent with FEMM's results and to have a FFT plot
            cond.SetValue("HysteresisLossCalcType", 1)
            cond.SetValue("PresetType", 3) # 3:Custom
            # Specify the reference steps yourself because you don't really know what JMAG is doing behind you
            cond.SetValue("StartReferenceStep", number_of_total_steps+1-number_of_steps_2ndTTS*0.5) # 1/4 period <=> number_of_steps_2ndTTS*0.5
            cond.SetValue("EndReferenceStep", number_of_total_steps)
            cond.SetValue("UseStartReferenceStep", 1)
            cond.SetValue("UseEndReferenceStep", 1)
            cond.SetValue("Cyclicity", 4) # specify reference steps for 1/4 period and extend it to whole period
            cond.SetValue("UseFrequencyOrder", 1)
            cond.SetValue("FrequencyOrder", "1-50") # Harmonics up to 50th orders 
        # Check CSV reults for iron loss (You cannot check this for Freq study) # CSV and save space
        study.GetStudyProperties().SetValue("CsvOutputPath", dir_csv_output_folder) # it's folder rather than file!
        study.GetStudyProperties().SetValue("CsvResultTypes", "Torque;Force;LineCurrent;TerminalVoltage;JouleLoss;TotalDisplacementAngle;JouleLoss_IronLoss;IronLoss_IronLoss;HysteresisLoss_IronLoss")
        study.GetStudyProperties().SetValue("DeleteResultFiles", self.fea_config_dict['delete_results_after_calculation'])
        # Terminal Voltage/Circuit Voltage: Check for outputing CSV results 
        study.GetCircuit().CreateTerminalLabel("Terminal4U", 8, -13)
        study.GetCircuit().CreateTerminalLabel("Terminal4V", 8, -11)
        study.GetCircuit().CreateTerminalLabel("Terminal4W", 8, -9)
        study.GetCircuit().CreateTerminalLabel("Terminal2U", 23, -13)
        study.GetCircuit().CreateTerminalLabel("Terminal2V", 23, -11)
        study.GetCircuit().CreateTerminalLabel("Terminal2W", 23, -9)
        # Export Stator Core's field results only for iron loss calculation (the csv file of iron loss will be clean with this setting)
            # study.GetMaterial(u"Rotor Core").SetValue(u"OutputResult", 0) # at least one part on the rotor should be output or else a warning "the jplot file does not contains displacement results when you try to calc. iron loss on the moving part." will pop up, even though I don't add iron loss condition on the rotor.
        # study.GetMeshControl().SetValue(u"AirRegionOutputResult", 0)
        study.GetMaterial("Shaft").SetValue("OutputResult", 0)
        study.GetMaterial("Cage").SetValue("OutputResult", 0)
        study.GetMaterial("Coil").SetValue("OutputResult", 0)
        # Rotor
        if True:
            cond = study.CreateCondition("Ironloss", "IronLossConRotor")
            cond.SetValue("BasicFrequencyType", 2)
            cond.SetValue("BasicFrequency", "freq")
                # cond.SetValue(u"BasicFrequency", u"slip*freq") # this require the signal length to be at least 1/4 of slip period, that's too long!
            cond.ClearParts()
            sel = cond.GetSelection()
            sel.SelectPartByPosition(-im_variant.Radius_Shaft-EPS, 0 ,0)
            cond.AddSelected(sel)
            # Use FFT for hysteresis to be consistent with FEMM's results
            cond.SetValue("HysteresisLossCalcType", 1)
            cond.SetValue("PresetType", 3)
            # Specify the reference steps yourself because you don't really know what JMAG is doing behind you
            cond.SetValue("StartReferenceStep", number_of_total_steps+1-number_of_steps_2ndTTS*0.5) # 1/4 period <=> number_of_steps_2ndTTS*0.5
            cond.SetValue("EndReferenceStep", number_of_total_steps)
            cond.SetValue("UseStartReferenceStep", 1)
            cond.SetValue("UseEndReferenceStep", 1)
            cond.SetValue("Cyclicity", 4) # specify reference steps for 1/4 period and extend it to whole period
            cond.SetValue("UseFrequencyOrder", 1)
            cond.SetValue("FrequencyOrder", "1-50") # Harmonics up to 50th orders 
        self.study_name = study_name
        return study

    # structural Study
    def add_structural_study(self, app, model, dir_csv_output_folder):
        study = model.CreateStudy(u"StructuralStatic2D", u"StructuralStatic")
        print('Csv output directory: %s'%(dir_csv_output_folder))
        study.GetStudyProperties().SetValue(u"CsvOutputPath", dir_csv_output_folder)
        study.GetStudyProperties().SetValue(u"CsvOutputResponseData", 1)
        study.GetStudyProperties().SetValue(u"CsvOutputProbe", 1)
        study.GetStudyProperties().SetValue(u"CsvOutputCalculation", 1)

        # Material
        if 'M19' in self.fea_config_dict['Steel']:
            study.SetMaterialByName("Rotor Core", "M-19 Steel Gauge-29")
            study.GetMaterial("Rotor Core").SetValue("Laminated", 1)
            study.GetMaterial("Rotor Core").SetValue("LaminationFactor", 95)

            study.SetMaterialByName("Shaft", "M-19 Steel Gauge-29")
            study.GetMaterial("Shaft").SetValue("Laminated", 1)
            study.GetMaterial("Shaft").SetValue("LaminationFactor", 95)
        else:
            raise Exception('Not implemented error.')

        study.SetMaterialByName("Cage", "Aluminium")
        study.GetMaterial("Cage").SetValue("EddyCurrentCalculation", 1)
        study.GetMaterial("Cage").SetValue("UserConductivityType", 1)
        study.GetMaterial("Cage").SetValue("UserConductivityValue", self.Bar_Conductivity)    

        # Conditions
        # Fixture
        study.CreateCondition(u"Constraint", u"FixtureConstraints")
        study.GetCondition(u"FixtureConstraints").SetValue(u"Type", 1) # Normal direction
        study.GetCondition(u"FixtureConstraints").ClearParts()
        study.GetCondition(u"FixtureConstraints").AddSet(model.GetSetList().GetSet(u"FixtureConstraints"), 0)

        # Contact
        study.CreateCondition(u"Contact", u"RotorCoreBar_Contact")
        study.GetCondition(u"RotorCoreBar_Contact").SetValue(u"ConsiderAdhesion", 1)
        study.GetCondition(u"RotorCoreBar_Contact").SetValue(u"LayerYoungsModulus", 0.025) # JAC223_IPM_Opti
        study.GetCondition(u"RotorCoreBar_Contact").SetValue(u"CalcTypeOfAdhesionLayer", 1)
        study.GetCondition(u"RotorCoreBar_Contact").SetValue(u"ThickOfAdhesionLayer", 0.001)
        study.GetCondition(u"RotorCoreBar_Contact").ClearParts()
        study.GetCondition(u"RotorCoreBar_Contact").AddSet(model.GetSetList().GetSet(u"RotorCoreShaft_Master"), 0)
        study.GetCondition(u"RotorCoreBar_Contact").AddSet(model.GetSetList().GetSet(u"RotorCoreShaft_Slave"), 1)

        study.CreateCondition(u"Contact", u"BarRotorCore_Contact")
        study.GetCondition(u"BarRotorCore_Contact").SetValue(u"ConsiderAdhesion", 1)
        study.GetCondition(u"BarRotorCore_Contact").SetValue(u"LayerYoungsModulus", 0.025)
        study.GetCondition(u"BarRotorCore_Contact").SetValue(u"CalcTypeOfAdhesionLayer", 1)
        study.GetCondition(u"BarRotorCore_Contact").SetValue(u"ThickOfAdhesionLayer", 0.001)
        study.GetCondition(u"BarRotorCore_Contact").ClearParts()
        study.GetCondition(u"BarRotorCore_Contact").AddSet(model.GetSetList().GetSet(u"BarRotorCore_Slave"), 0) # 这里注意一下，根据JAC223，我们认为 Contact 条件
        study.GetCondition(u"BarRotorCore_Contact").AddSet(model.GetSetList().GetSet(u"BarRotorCore_Master"), 1)

        # Centrifugal Force
        study.CreateCondition(u"CentrifugalForce", u"Centrifugal_Force")
        print(self.DriveW_Freq*60. / (0.5*self.DriveW_poles),'r/min')
        study.GetCondition(u"Centrifugal_Force").SetValue(u"AngularVelocity", self.DriveW_Freq*60. / (0.5*self.DriveW_poles))
        study.GetCondition(u"Centrifugal_Force").SetXYZPoint(u"Axis", 0, 0, 1)
        study.GetCondition(u"Centrifugal_Force").ClearParts()
        study.GetCondition(u"Centrifugal_Force").AddSet(model.GetSetList().GetSet(u"Motion_Region"), 0)

        # Mesh Size Control
        study.GetMeshControl().CreateCondition(u"Part", u"MeshSizeControl")
        study.GetMeshControl().GetCondition(u"MeshSizeControl").SetValue(u"Size", 0.100) # mm
        study.GetMeshControl().GetCondition(u"MeshSizeControl").ClearParts()
        study.GetMeshControl().GetCondition(u"MeshSizeControl").AddSet(model.GetSetList().GetSet(u"Motion_Region"), 0)
        study.CreateMesh()
        app.View().ShowMesh()

        return study

    # Static FEA
    def add_cad_parameters(self, study):

        ''' CAD Parameters '''
        study.SetCheckForTopologyChanges(False) # 

        # the order matters for param_no to begin at 0.
        study.GetDesignTable().AddEquation("RotorPosition")
        study.GetDesignTable().GetEquation("RotorPosition").SetType(0) # 0-value. 1-expression which cannot be modified for different cases.
        study.GetDesignTable().GetEquation("RotorPosition").SetExpression("0.0")
        study.GetDesignTable().GetEquation("RotorPosition").SetDescription("rotor position")
        # study.GetDesignTable().SetValue(0, 0, 0.0) # case_no, param_no, value

        # app.GetModel(u"EC_Rotate_32").RestoreCadLink()
        def add_vertex_xy_as_cad_param(list_vertex_names, sketch_name):
            for vertex_name in list_vertex_names:
                study.AddCadParameter("X@%s@%s"%(vertex_name,sketch_name))
                study.AddCadParameter("Y@%s@%s"%(vertex_name,sketch_name))
                study.GetDesignTable().AddCadParameterVariableName("X value@%s@%s"%(vertex_name,sketch_name))
                study.GetDesignTable().AddCadParameterVariableName("Y value@%s@%s"%(vertex_name,sketch_name))

        #print im.list_rotorCore_vertex_names
        add_vertex_xy_as_cad_param(self.list_rotorCore_vertex_names, 'Rotor Core')
        add_vertex_xy_as_cad_param(self.list_rotorCage_vertex_names, 'Cage')
            # add_vertex_xy_as_cad_param(self.list_rotorAirWithin_vertex_names, 'Air Within Rotor Slots')

        # total_number_of_cad_parameters = len(self.list_rotorCore_vertex_names) + len(self.list_rotorCage_vertex_names)

    def add_cases_rotate_rotor(self, study, total_number_of_cases):
        print('total_number_of_cases:', total_number_of_cases)
        if total_number_of_cases > 1:

            # add case label!
            study.GetDesignTable().AddCases(total_number_of_cases - 1)    

            def rotate_vertex_in_cad_param(theta, case_no_list, param_no):
                for case_no in case_no_list:
                    # print case_no, total_number_of_cases
                    X = study.GetDesignTable().GetValue(0, param_no)
                    Y = study.GetDesignTable().GetValue(0, param_no+1)
                    radian = case_no*theta
                    study.GetDesignTable().SetValue(case_no, begin_here-1, radian/pi*180)
                    # print case_no, radian/pi*180, case_no_list
                    # Park Transformation
                    NEW_X = cos(radian)*X - sin(radian)*Y # 注意转动的方向，要和磁场旋转的方向一致
                    NEW_Y = sin(radian)*X + cos(radian)*Y
                    study.GetDesignTable().SetValue(case_no, param_no, NEW_X)
                    study.GetDesignTable().SetValue(case_no, param_no+1, NEW_Y)

            begin_here = 4 # rotor_position, freq, slip, speed, 
            study.GetDesignTable().SetValue(0, begin_here-1, 0.0) # rotor_position for case#0 is 0.0 # BUG???
            end_here = study.GetDesignTable().NumParameters()
            # print 'This value should be 4+34:', end_here

            # self.theta = 1./180.0*pi

            for param_no in range(begin_here, end_here, 2): # begin at one for the RotorPosition Variable
                rotate_vertex_in_cad_param(self.theta, list(range(1,total_number_of_cases)), param_no) # case_no = 1
            study.ApplyCadParameters()
                    # study.GetDesignTable().SetActive(0, True)

    def add_rotor_current_condition(self, app, model, study, total_number_of_cases, eddy_current_circuit_current_csv_file): # r'D:\Users\horyc\OneDrive - UW-Madison\csv\Freq_#4_circuit_current.csv'

        # study.GetMaterial(u"Cage").SetValue(u"EddyCurrentCalculation", 0)
        rotor_phase_name_list = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

        # read from eddy current results
        dict_circuit_current_complex = {}
        with open(eddy_current_circuit_current_csv_file, 'r') as f:
            for row in self.csv_row_reader(f):
                try: 
                    float(row[0])
                except:
                    continue
                else:
                    if np.abs(self.slip_freq_breakdown_torque - float(row[0])) < 1e-3:                
                        beginning_column = 1 + 2*3*2 # title + drive/bearing * 3 phase * real/imag
                        for i in range(0, int(self.no_slot_per_pole)):
                            natural_i = i+1
                            current_phase_column = beginning_column + i * int(self.DriveW_poles) * 2
                            for j in range(int(self.DriveW_poles)):
                                natural_j = j+1
                                re = float(row[current_phase_column+2*j])
                                im = float(row[current_phase_column+2*j+1])
                                dict_circuit_current_complex["%s%d"%(rotor_phase_name_list[i], natural_j)] = (re, im)
        dict_circuit_current_amp_and_phase = {}
        for key, item in dict_circuit_current_complex.items():
            amp = np.sqrt(item[1]**2 + item[0]**2)
            phase = np.arctan2(item[0], -item[1]) # atan2(y, x), y=a, x=-b
            dict_circuit_current_amp_and_phase[key] = (amp, phase)


        # Link the FEMCoil with the conditions
        begin_parameters = study.GetDesignTable().NumParameters()
        # print 'num param:', begin_parameters
        count_parameters = 0
        for i in range(0, int(self.no_slot_per_pole)):

            # values of CurCon are determined by equations/variables: so buld variables for that
            study.GetDesignTable().AddEquation( 'var' + "RotorCurCon%s"%(rotor_phase_name_list[i]) )
            study.GetDesignTable().GetEquation( 'var' + "RotorCurCon%s"%(rotor_phase_name_list[i]) ).SetType(0)
            study.GetDesignTable().GetEquation( 'var' + "RotorCurCon%s"%(rotor_phase_name_list[i]) ).SetExpression("0")
            study.GetDesignTable().GetEquation( 'var' + "RotorCurCon%s"%(rotor_phase_name_list[i]) ).SetDescription("")

            # now, assign different values to these variables w.r.t. different cases 
            # self.theta = 1./180.0*pi
            t = 0.0
            time_one_step = self.theta / (2*pi * self.the_speed) * 60 # sec
            # ScriptComments
            amp, phase = dict_circuit_current_amp_and_phase[rotor_phase_name_list[i]+'1']
            for case_no in range(total_number_of_cases):
                current_value = amp * sin(2*pi*(self.the_slip*self.DriveW_Freq)*t + phase)
                study.GetDesignTable().SetValue(case_no, begin_parameters + count_parameters, current_value)
                t += time_one_step
            count_parameters += 1

            natural_i = i+1
            for index, j in enumerate([ natural_i,
                                        natural_i + self.no_slot_per_pole]):

                # for Static FEA: amp * sin(2*fpi*f*t + phase)
                current_condition = study.CreateCondition("Current", "RotorCurCon%s%d"%(rotor_phase_name_list[i],j))
                current_condition.SetValue("Current", 'var' + "RotorCurCon%s"%(rotor_phase_name_list[i]))
                    # current_condition.SetValue(u"XType", 1)
                if index == 0:
                    current_condition.SetValue("Direction2D", 0)
                else:
                    current_condition.SetValue("Direction2D", 1) # 结果不正常？可能是转子电流相位差了180度？ Check.

                current_condition.ClearParts()
                current_condition.AddSet(model.GetSetList().GetSet("Bar %d"%(j)), 0)
                current_condition.AddSet(model.GetSetList().GetSet("Bar %d"%(j+2*self.no_slot_per_pole)), 0)

        # for key, item in dict_circuit_current_complex.iteritems():
        #     print key, dict_circuit_current_complex[key],
        #     print dict_circuit_current_amp_and_phase[key]

    def add_stator_current_condition(self, app, model, study, total_number_of_cases, eddy_current_circuit_current_csv_file):

        # read from eddy current results
        dict_circuit_current_complex = {}
        with open(eddy_current_circuit_current_csv_file, 'r') as f:
            for row in self.csv_row_reader(f):
                try: 
                    float(row[0])
                except:
                    continue
                else:
                    if '%g'%(self.slip_freq_breakdown_torque) in row[0]:
                        beginning_column = 1 # title column is not needed
                        for i, phase in zip(list(range(0,12,2)), ['2A','2B','2C','4A','4B','4C']): # 3 phase
                            natural_i = i+1
                            current_phase_column = beginning_column + i 
                            re = float(row[current_phase_column])
                            im = float(row[current_phase_column+1])
                            dict_circuit_current_complex[phase] = (re, im)
        dict_circuit_current_amp_and_phase = {}
        for key, item in dict_circuit_current_complex.items():
            amp = np.sqrt(item[1]**2 + item[0]**2)
            phase = np.arctan2(item[0], -item[1]) # atan2(y, x), y=a, x=-b
            dict_circuit_current_amp_and_phase[key] = (amp, phase)

        # for key, item in dict_circuit_current_complex.iteritems():
        #     print key, dict_circuit_current_complex[key],
        #     print dict_circuit_current_amp_and_phase[key]


        ''' Create Variables for being used for different cases '''
        begin_parameters = study.GetDesignTable().NumParameters()
        count_parameters = 0
        print('num param begins at:', begin_parameters)
        for phase_name in ['2A','2B','2C','4A','4B','4C']:
            # ScriptComments
            amp, phase = dict_circuit_current_amp_and_phase[phase_name]

            # values of CurCon are determined by equations/variables: so buld variables for that
            study.GetDesignTable().AddEquation( 'var' + phase_name )
            study.GetDesignTable().GetEquation( 'var' + phase_name ).SetType(0)
            study.GetDesignTable().GetEquation( 'var' + phase_name ).SetExpression("0")

            # now, assign different values to these variables w.r.t. different cases 
            # self.theta = 1./180.0*pi
            t = 0.0
            time_one_step = self.theta / (2*pi * self.the_speed) * 60 # sec
            for case_no in range(total_number_of_cases):
                current_value = amp * sin(2*pi*self.DriveW_Freq*t + phase)
                study.GetDesignTable().SetValue(case_no, begin_parameters + count_parameters, current_value)
                t += time_one_step 
            count_parameters += 1
        print('num param ends at:', study.GetDesignTable().NumParameters())



        def create_stator_current_conditions(turns,condition_name_list):

            set_list = model.GetSetList()

            for condition_name in condition_name_list:
                condition = study.CreateCondition("Current", condition_name)
                condition.SetValue("Turn", turns)
                condition.SetValue("Current", 'var' + condition_name[4:6])
                if '+' in condition_name:
                    condition.SetValue("Direction2D", 1) # 1=x
                elif '-' in condition_name:
                    condition.SetValue("Direction2D", 0) # 0=.
                else:
                    raise Exception('cannot find + or i in current condition name')

                condition.ClearParts()
                set_name_list = [set_list.GetSet(i).GetName() for i in range(set_list.NumSet()) if condition_name[4:7] in set_list.GetSet(i).GetName()]
                print(set_name_list)
                for set_name in set_name_list:
                    condition.AddSet(model.GetSetList().GetSet(set_name), 0) # 0: group
        create_stator_current_conditions(self.DriveW_zQ, ["Coil4A-","Coil4B-","Coil4C-","Coil4A+","Coil4B+","Coil4C+"])
        create_stator_current_conditions(self.DriveW_zQ, ["Coil2A-","Coil2B-","Coil2C-","Coil2A+","Coil2B+","Coil2C+"])

    def add_rotor_current_condition_obsolete_slow_version(self, app, model, study, total_number_of_cases, eddy_current_circuit_current_csv_file): # r'D:\Users\horyc\OneDrive - UW-Madison\csv\Freq_#4_circuit_current.csv'

        # study.GetMaterial(u"Cage").SetValue(u"EddyCurrentCalculation", 0)
        rotor_phase_name_list = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

        # read from eddy current results
        dict_circuit_current_complex = {}
        with open(eddy_current_circuit_current_csv_file, 'r') as f:
            for row in self.csv_row_reader(f):
                try: 
                    float(row[0])
                except:
                    continue
                else:
                    if '%g'%(self.slip_freq_breakdown_torque) in row[0]:
                        beginning_column = 1 + 2*3*2 # title + drive/bearing * 3 phase * real/imag
                        for i in range(0, int(self.no_slot_per_pole)):
                            natural_i = i+1
                            current_phase_column = beginning_column + i * int(self.DriveW_poles) * 2
                            for j in range(int(self.DriveW_poles)):
                                natural_j = j+1
                                re = float(row[current_phase_column+2*j])
                                im = float(row[current_phase_column+2*j+1])
                                dict_circuit_current_complex["%s%d"%(rotor_phase_name_list[i], natural_j)] = (re, im)
        dict_circuit_current_amp_and_phase = {}
        for key, item in dict_circuit_current_complex.items():
            amp = np.sqrt(item[1]**2 + item[0]**2)
            phase = np.arctan2(item[0], -item[1]) # atan2(y, x), y=a, x=-b
            dict_circuit_current_amp_and_phase[key] = (amp, phase)


        # Link the FEMCoil with the conditions
        begin_parameters = study.GetDesignTable().NumParameters()
        # print 'num param:', begin_parameters
        count_parameters = 0
        for i in range(0, int(self.no_slot_per_pole)):
            natural_i = i+1
            for index0123, j in enumerate([  natural_i,
                                natural_i +   self.no_slot_per_pole,
                                natural_i + 2*self.no_slot_per_pole,
                                natural_i + 3*self.no_slot_per_pole]):
                # study.GetCondition(u"BarCoilCon %d"%(natural_j)).SetLink(u"ConductorCircuit %d"%(j))

                # ScriptComments
                amp, phase = dict_circuit_current_amp_and_phase["%s%d"%(rotor_phase_name_list[i], index0123+1)]

                # values of CurCon are determined by equations/variables: so buld variables for that
                study.GetDesignTable().AddEquation( 'var' + "RotorCurCon%s%d"%(rotor_phase_name_list[i],j) )
                study.GetDesignTable().GetEquation( 'var' + "RotorCurCon%s%d"%(rotor_phase_name_list[i],j) ).SetType(0)
                study.GetDesignTable().GetEquation( 'var' + "RotorCurCon%s%d"%(rotor_phase_name_list[i],j) ).SetExpression("0")
                study.GetDesignTable().GetEquation( 'var' + "RotorCurCon%s%d"%(rotor_phase_name_list[i],j) ).SetDescription("")

                # for Static FEA: amp * sin(2*fpi*f*t + phase)
                current_condition = study.CreateCondition("Current", "RotorCurCon%s%d"%(rotor_phase_name_list[i],j))
                current_condition.SetValue("XType", 1)
                current_condition.SetValue("Current", 'var' + "RotorCurCon%s%d"%(rotor_phase_name_list[i],j))
                current_condition.ClearParts()
                current_condition.AddSet(model.GetSetList().GetSet("Bar %d"%(j)), 0)

                # now, assign different values to these variables w.r.t. different cases 
                # self.theta = 1./180.0*pi
                t = 0.0
                time_one_step = self.theta / (2*pi * self.the_speed) * 60 # sec
                for case_no in range(total_number_of_cases):
                    t += time_one_step
                    current_value = amp * sin(2*pi*(self.the_slip*self.DriveW_Freq)*t + phase)
                    study.GetDesignTable().SetValue(case_no, begin_parameters + count_parameters, current_value)
                count_parameters += 1


        # for key, item in dict_circuit_current_complex.iteritems():
        #     print key, dict_circuit_current_complex[key],
        #     print dict_circuit_current_amp_and_phase[key]

class VanGogh_JMAG(VanGogh):
    def __init__(self, im, child_index=1, doNotRotateCopy=False):
        super(VanGogh_JMAG, self).__init__(im, child_index)
        self.doNotRotateCopy = doNotRotateCopy

        self.SketchName = None
        self.dict_count_arc = {}
        self.dict_count_region = {}

        self.count = 0 # counter of region

    def mirror_and_copyrotate(self, Q, Radius, fraction, 
                                edge4ref=None,
                                symmetry_type=None,
                                merge=True,
                                do_you_have_region_in_the_mirror=False
                                ):

        region = self.create_region([art.GetName() for art in self.artist_list]) 

        self.region_mirror_copy(region, edge4ref=edge4ref, symmetry_type=symmetry_type, merge=merge)
        self.count+=1
        # if self.count == 4: # debug
            # raise Exception
            # merge = True # When overlap occurs between regions because of copying, a boolean operation (sum) is executed and they are merged into 1 region.
        if not self.doNotRotateCopy:
            self.region_circular_pattern_360_origin(region, float(Q), merge=merge,
                                                    do_you_have_region_in_the_mirror=do_you_have_region_in_the_mirror)
        # print self.artist_list
        self.sketch.CloseSketch()

    def draw_arc(self, center, p1, p2, **kwarg):
        art = self.sketch.CreateArc(center[0], center[1], p1[0], p1[1], p2[0], p2[1])
        self.artist_list.append(art)

    def draw_line(self, p1, p2):
        # return self.line(p1[0],p1[1],p2[0],p2[1])
        art = self.sketch.CreateLine(p1[0],p1[1],p2[0],p2[1])
        self.artist_list.append(art)

    def add_line(self, p1, p2):
        draw_line(p1, p2)

    def plot_sketch_shaft(self):
        self.SketchName = "Shaft"
        SketchName = self.SketchName
        sketch = self.create_sketch(SketchName, "#D1B894")

        if self.doNotRotateCopy:
            Rotor_Sector_Angle = 2*pi/self.im.Qr*0.5
            PA = [ self.im.Radius_Shaft*-cos(Rotor_Sector_Angle), self.im.Radius_Shaft*sin(Rotor_Sector_Angle) ]
            PB = [ PA[0], -PA[1] ]
            self.draw_arc([0,0], PA, PB)
            self.draw_line(PA, [0,0])
            self.draw_line(PB, [0,0])
            self.doc.GetSelection().Clear()
            self.doc.GetSelection().Add(sketch.GetItem("Arc"))
            self.doc.GetSelection().Add(sketch.GetItem("Line"))
            self.doc.GetSelection().Add(sketch.GetItem("Line.2"))
        else:
            self.circle(0, 0, self.im.Radius_Shaft)

            self.doc.GetSelection().Clear()
            self.doc.GetSelection().Add(sketch.GetItem("Circle"))
        sketch.CreateRegions()

        sketch.CloseSketch()
        # sketch.SetProperty(u"Visible", 0)

    def init_sketch_statorCore(self):
        self.SketchName="Stator Core"
        sketch = self.create_sketch(self.SketchName, "#E8B5CE")
        return 

    def init_sketch_coil(self):
        self.SketchName="Coil"
        sketch = self.create_sketch(self.SketchName, "#EC9787")
        return

    def init_sketch_rotorCore(self):
        self.SketchName="Rotor Core"
        sketch = self.create_sketch(self.SketchName, "#FE840E")
        return

    def init_sketch_cage(self):
        self.SketchName="Cage"
        sketch = self.create_sketch(self.SketchName, "#8D9440")
        return


    # obsolete
    def draw_arc_using_shapely(self, p1, p2, angle, maxseg=1): # angle in rad
        center = self.find_center_of_a_circle_using_2_points_and_arc_angle(p1, p2, angle) # ordered p1 and p2 are
        art = self.sketch.CreateArc(center[0], center[1], p1[0], p1[1], p2[0], p2[1])
        self.artist_list.append(art)

    def add_arc_using_shapely(self, p1, p2, angle, maxseg=1): # angle in rad
        self.draw_arc_using_shapely(p1, p2, angle, maxseg)


    # Utility wrap function for JMAG
    def create_sketch(self, SketchName, color):
        self.artist_list = []

        try:self.dict_count_arc[SketchName]
        except: self.dict_count_arc[SketchName] = 0
        try:self.dict_count_region[SketchName]
        except: self.dict_count_region[SketchName] = 0
        ref1 = self.ass.GetItem("XY Plane")
        ref2 = self.doc.CreateReferenceFromItem(ref1)
        self.sketch = self.ass.CreateSketch(ref2)
        self.sketch.OpenSketch()
        self.sketch.SetProperty("Name", SketchName)
        self.sketch.SetProperty("Color", color)
        return self.sketch
    def circle(self, x,y,r):
        # SketchName = self.SketchName
        self.sketch.CreateVertex(x, y)
        # return self.circle(x, y, r)
        return self.sketch.CreateCircle(x, y, r)
    def line(self, x1,y1,x2,y2):
        # SketchName = self.SketchName
        self.sketch.CreateVertex(x1,y1)
        self.sketch.CreateVertex(x2,y2)
        # return self.line(x1,y1,x2,y2)
        return self.sketch.CreateLine(x1,y1,x2,y2)
    def create_region(self, l):
        SketchName = self.SketchName
        self.doc.GetSelection().Clear()
        for art_name in l:
            self.doc.GetSelection().Add(self.sketch.GetItem(art_name))
            # self.doc.GetSelection().Add(el)
        self.sketch.CreateRegions() # this returns None
        # self.sketch.CreateRegionsWithCleanup(0.05, True) # mm. difference at stator outter radius is up to 0.09mm! This turns out to be neccessary for shapely to work with JMAG. Shapely has poor 

        self.dict_count_region[SketchName] += 1
        if self.dict_count_region[SketchName]==1:
            return self.sketch.GetItem("Region")
        else:
            return self.sketch.GetItem("Region.%d"%(self.dict_count_region[SketchName]))
    def region_mirror_copy(self, region, edge4ref=None, symmetry_type=None, merge=True):
        mirror = self.sketch.CreateRegionMirrorCopy()
        mirror.SetProperty("Merge", merge)
        ref2 = self.doc.CreateReferenceFromItem(region)
        mirror.SetPropertyByReference("Region", ref2)

        # å¯¹ç§°è½´
        if edge4ref == None:
            if symmetry_type == None:
                print("At least give one of edge4ref and symmetry_type")
                raise Exception
            else:
                mirror.SetProperty("SymmetryType", symmetry_type)
        else:
            ref1 = self.sketch.GetItem(edge4ref.GetName()) # e.g., u"Line"
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            mirror.SetPropertyByReference("Symmetry", ref2)

        # print region
        # print ass.GetItem(u"Region.1")
        if merge == False and region.GetName()=="Region":
            return self.ass.GetItem("Region.1")
    def region_circular_pattern_360_origin(self, region, Q_float, merge=True, do_you_have_region_in_the_mirror=False):
        circular_pattern = self.sketch.CreateRegionCircularPattern()
        circular_pattern.SetProperty("Merge", merge)

        ref2 = self.doc.CreateReferenceFromItem(region)
        circular_pattern.SetPropertyByReference("Region", ref2)
        face_region_string = circular_pattern.GetProperty("Region")
        if isinstance(face_region_string, tuple): #type(face_region_string) == type(tuple()):
            print('[DEBUG] TUPLE type face_region_string', face_region_string)
            face_region_string = face_region_string[0]
        else:
            print('[DEBUG] INT type face_region_string', face_region_string)
            # face_region_string = face_region_string[0]

        if do_you_have_region_in_the_mirror == True:
            # è¿™é‡Œå‡è®¾face_region_stringæœ€åŽä¸¤ä½æ˜¯æ•°å­—
            if face_region_string[-7:-3] == 'Item':
                number_plus_1 = str(int(face_region_string[-3:-1]) + 1)
                refarray = [0 for i in range(2)]
                refarray[0] = "faceregion(TRegionMirrorPattern%s+%s_2)" % (number_plus_1, face_region_string)
                refarray[1] = face_region_string
                circular_pattern.SetProperty("Region", refarray)
                # print refarray[0]
                # print refarray[1]
            elif face_region_string[-6:-2] == 'Item':
                # è¿™é‡Œå‡è®¾face_region_stringæœ€åŽä¸€ä½æ˜¯æ•°å­—
                number_plus_1 = str(int(face_region_string[-2:-1]) + 1)
                refarray = [0 for i in range(2)]
                refarray[0] = "faceregion(TRegionMirrorPattern%s+%s_2)" % (number_plus_1, face_region_string)
                refarray[1] = face_region_string
                circular_pattern.SetProperty("Region", refarray)
            elif face_region_string[-8:-4] == 'Item':
                # è¿™é‡Œå‡è®¾face_region_stringæœ€åŽä¸‰ä½æ˜¯æ•°å­—
                number_plus_1 = str(int(face_region_string[-4:-1]) + 1)
                refarray = [0 for i in range(2)]
                refarray[0] = "faceregion(TRegionMirrorPattern%s+%s_2)" % (number_plus_1, face_region_string)
                refarray[1] = face_region_string
                circular_pattern.SetProperty("Region", refarray)



        circular_pattern.SetProperty("CenterType", 2)
        circular_pattern.SetProperty("Angle", "360/%d"%(int(Q_float)))
        circular_pattern.SetProperty("Instance", int(Q_float))

# if __name__ == '__main__':
#     vg_jmag = VanGogh_JMAG(None)
#     print vg_jmag.find_center_of_a_circle_using_2_points_and_arc_angle

class TrimDrawer(object):
    ''' JMAG's geometry editor works like SolidWorks. '''
    def __init__(self, im):
        self.SketchName = None
        self.trim_a = self.trim_l
        self.dict_count_arc = {}
        self.dict_count_region = {}

        self.im = im

    ''' Basic Functions '''
    def create_sketch(self, SketchName, color):
        try:self.dict_count_arc[SketchName]
        except: self.dict_count_arc[SketchName] = 0
        try:self.dict_count_region[SketchName]
        except: self.dict_count_region[SketchName] = 0
        ref1 = self.ass.GetItem("XY Plane")
        ref2 = self.doc.CreateReferenceFromItem(ref1)
        self.sketch = self.ass.CreateSketch(ref2)
        self.sketch.OpenSketch()
        self.sketch.SetProperty("Name", SketchName)
        self.sketch.SetProperty("Color", color)
        return self.sketch

    def circle(self, x,y,r):
        # SketchName = self.SketchName
        # self.sketch.CreateVertex(x, y)
        # return self.circle(x, y, r)
        return self.sketch.CreateCircle(x, y, r)

    def line(self, x1,y1,x2,y2):
        # SketchName = self.SketchName
        # self.sketch.CreateVertex(x1,y1)
        # self.sketch.CreateVertex(x2,y2)
        # return self.line(x1,y1,x2,y2)
        return self.sketch.CreateLine(x1,y1,x2,y2)

    def trim_l(self, who,x,y):
        # SketchName = self.SketchName
        self.doc.GetSelection().Clear()
        ref1 = self.sketch.GetItem(who.GetName())
        self.doc.GetSelection().Add(ref1)
        self.doc.GetSketchManager().SketchTrim(x,y)
        # l1 trim 完以后还是l1，除非你切中间，这样会多生成一个Line，你自己捕捉一下吧

    def trim_c(self, who,x,y):
        SketchName = self.SketchName
        self.doc.GetSelection().Clear()
        ref1 = self.sketch.GetItem(who.GetName())
        self.doc.GetSelection().Add(ref1)
        self.doc.GetSketchManager().SketchTrim(x,y)

        # print who 
        self.dict_count_arc[SketchName] += 1
        if self.dict_count_arc[SketchName]==1:
            return self.sketch.GetItem("Arc")
        else:
            return self.sketch.GetItem("Arc.%d"%(self.dict_count_arc[SketchName]))

    def create_region(self, l):
        SketchName = self.SketchName
        self.doc.GetSelection().Clear()
        for string in l:
            self.doc.GetSelection().Add(self.sketch.GetItem(string))
            # self.doc.GetSelection().Add(el)
        self.sketch.CreateRegions() # this returns None

        self.dict_count_region[SketchName] += 1
        if self.dict_count_region[SketchName]==1:
            return self.sketch.GetItem("Region")
        else:
            return self.sketch.GetItem("Region.%d"%(self.dict_count_region[SketchName]))

    def region_mirror_copy(self, region, edge4ref=None, symmetry_type=None, merge=True):
        mirror = self.sketch.CreateRegionMirrorCopy()
        mirror.SetProperty("Merge", merge)
        ref2 = self.doc.CreateReferenceFromItem(region)
        mirror.SetPropertyByReference("Region", ref2)

        # å¯¹ç§°è½´
        if edge4ref == None:
            if symmetry_type == None:
                print("At least give one of edge4ref and symmetry_type")
            else:
                mirror.SetProperty("SymmetryType", symmetry_type)
        else:
            ref1 = self.sketch.GetItem(edge4ref.GetName()) # e.g., u"Line"
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            mirror.SetPropertyByReference("Symmetry", ref2)

        # print region
        # print ass.GetItem(u"Region.1")
        if merge == False and region.GetName()=="Region":
            return self.ass.GetItem("Region.1")
            # å…¶ä»–æƒ…å†µè¿˜æ²¡å†™å“¦ï¼Œæœ€å\½ä¿è¯ä½ å°±è¿™ä¸€ä¸ªregion

            # Region_Mirror_Copy_2 = sketch.CreateRegionMirrorCopy()
            # sketch.GetItem(u"Region Mirror Copy").SetProperty(u"Merge", 1)
            # refarray = [0 for i in range(1)]
            # refarray[0] = u"faceregion(TRegionItem68)" # åœ¨å½•è¿™ä¸€æ­\å‰ï¼Œå¿…é¡»æŠŠGeometry Editorç»™å…³äº†ï¼Œé‡æ–°è·‘ä¸€éå‰é¢çš„ä»£ç ï¼Œè¿™ä¸ªæ•°å­—æ‰å¯¹ï¼
            # sketch.GetItem(u"Region Mirror Copy").SetProperty(u"Region", refarray)
            # ref1 = sketch.GetItem(u"Line")
            # ref2 = self.doc.CreateReferenceFromItem(ref1)
            # sketch.GetItem(u"Region Mirror Copy").SetPropertyByReference(u"Symmetry", ref2)

    def region_circular_pattern_360_origin(self, region, Q_float, merge=True, do_you_have_region_in_the_mirror=False):
        circular_pattern = self.sketch.CreateRegionCircularPattern()
        circular_pattern.SetProperty("Merge", merge)

        ref2 = self.doc.CreateReferenceFromItem(region)
        circular_pattern.SetPropertyByReference("Region", ref2)
        face_region_string = circular_pattern.GetProperty("Region")
        if isinstance(face_region_string, tuple): #type(face_region_string) == type(tuple()):
            print('[DEBUG] TUPLE type face_region_string', face_region_string)
            face_region_string = face_region_string[0]
        else:
            print('[DEBUG] INT type face_region_string', face_region_string)
            # face_region_string = face_region_string[0]

        if do_you_have_region_in_the_mirror == True:
            # 这里假设face_region_string最后两位是数字
            # 乱码的原因是因为我有一次手贱，用DOS-CMD把所有文件的大小扫了一遍还是怎么的，中文就乱码了。
            # 20181114 (Designer from scratch) 没有乱码哦（好像是从Onedrive上找回的）
            # 总结JMAG的代码生成规律……
            # è¿™é‡Œå‡è®¾face_region_stringæœ€åŽä¸¤ä½æ˜¯æ•°å­—
            if face_region_string[-7:-3] == 'Item':
                number_plus_1 = str(int(face_region_string[-3:-1]) + 1)
                refarray = [0 for i in range(2)]
                refarray[0] = "faceregion(TRegionMirrorPattern%s+%s_2)" % (number_plus_1, face_region_string)
                refarray[1] = face_region_string
                circular_pattern.SetProperty("Region", refarray)
                # print refarray[0]
                # print refarray[1]
            elif face_region_string[-6:-2] == 'Item':
                # 这里假设face_region_string最后一位是数字
                # è¿™é‡Œå‡è®¾face_region_stringæœ€åŽä¸€ä½æ˜¯æ•°å­—
                number_plus_1 = str(int(face_region_string[-2:-1]) + 1)
                refarray = [0 for i in range(2)]
                refarray[0] = "faceregion(TRegionMirrorPattern%s+%s_2)" % (number_plus_1, face_region_string)
                refarray[1] = face_region_string
                circular_pattern.SetProperty("Region", refarray)
            elif face_region_string[-8:-4] == 'Item':
                # 这里假设face_region_string最后三位是数字
                # è¿™é‡Œå‡è®¾face_region_stringæœ€åŽä¸‰ä½æ˜¯æ•°å­—
                number_plus_1 = str(int(face_region_string[-4:-1]) + 1)
                refarray = [0 for i in range(2)]
                refarray[0] = "faceregion(TRegionMirrorPattern%s+%s_2)" % (number_plus_1, face_region_string)
                refarray[1] = face_region_string
                circular_pattern.SetProperty("Region", refarray)




        circular_pattern.SetProperty("CenterType", 2)
        circular_pattern.SetProperty("Angle", "360/%d"%(int(Q_float)))
        circular_pattern.SetProperty("Instance", int(Q_float))

        # if merge == False:


    # I add constraint only for rotate the model for static FEA in JMAG
    def constraint_fixture_circle_center(self, c):
        sketch = self.ass.GetItem(self.SketchName) # ass is global
        ref1 = c.GetCenterVertex()
        ref2 = self.doc.CreateReferenceFromItem(ref1)
        constraint = sketch.CreateMonoConstraint("fixture", ref2)
        # constraint.SetProperty(u"Name", constraint_name)

    def constraint_radius_arc(self, arc, radius_value, constraint_name):
        sketch = self.ass.GetItem(self.SketchName)
        ref1 = arc
        ref2 = self.doc.CreateReferenceFromItem(ref1)
        constraint = sketch.CreateMonoConstraint("radius", ref2)
        constraint.SetProperty("Radius", radius_value)
        constraint.SetProperty("Name", constraint_name)

    def rigidset(self, vertex_name_list): # this won't work, reason is mystery
        sketch = self.ass.GetItem(self.SketchName)
        sketch.CreateConstraint("rigidset")
        # var_list = ['ref%d'%(i+1) for i in range(2*len(vertex_name_list))]
        var_list = [None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None]
        for ind, vtx_name in enumerate(vertex_name_list):
            var_list[ind]   = sketch.GetItem(vtx_name)
            var_list[ind+1] = self.doc.CreateReferenceFromItem(ref1)
            sketch.GetItem("Relative Fixation").SetPropertyByReference("TargetList", var_list[ind+1])


    # global SketchName 
    ''' Parts to Plot '''
    def plot_shaft(self, name=None):
        if name == None:
            self.SketchName="Shaft"
        else:
            self.SketchName=name
        SketchName = self.SketchName
        sketch = self.create_sketch(SketchName, "#D1B894")

        self.circle(0, 0, self.im.Radius_Shaft)

        self.doc.GetSelection().Clear()
        self.doc.GetSelection().Add(sketch.GetItem("Circle"))
        sketch.CreateRegions()

        sketch.CloseSketch()
        # sketch.SetProperty(u"Visible", 0)
    def plot_rotorCore(self, name=None):
        if name == None:
            self.SketchName="Rotor Core"
        else:
            self.SketchName=name
        SketchName = self.SketchName
        sketch = self.create_sketch(SketchName, "#FE840E")

        c1=self.circle(0, 0, self.im.Radius_Shaft) # Circle.1
        c2=self.circle(0, 0, self.im.Radius_OuterRotor) # Circle.2

        Vertex_RotorBarCenter = sketch.CreateVertex(-self.im.Location_RotorBarCenter, 0)
        c3=self.circle(-self.im.Location_RotorBarCenter, 0, self.im.Radius_of_RotorSlot) # Circle.3

        l1=self.line(-5.5-self.im.Radius_OuterRotor, 0.5*self.im.Width_RotorSlotOpen, -self.im.Location_RotorBarCenter, 0.5*self.im.Width_RotorSlotOpen) # Line.1 # -5.5 is arbitrary float <0

        ref1 = sketch.GetItem("Line")
        ref2 = self.doc.CreateReferenceFromItem(ref1)
        sketch.CreateMonoConstraint("horizontality", ref2) # how to set constraint

        if self.im.use_drop_shape_rotor_bar == True:
            l2=self.line(0, 0, -self.im.Location_RotorBarCenter2+self.im.Radius_of_RotorSlot2, 0) # Line.2
            ref1 = sketch.GetItem("Line.2")
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            sketch.CreateMonoConstraint("horizontality", ref2)    
        else:
            l2=self.line(0, 0, -self.im.Location_RotorBarCenter+self.im.Radius_of_RotorSlot, 0) # Line.2
            ref1 = sketch.GetItem("Line.2")
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            sketch.CreateMonoConstraint("horizontality", ref2)

        R = self.im.Radius_OuterRotor
        THETA = (180-0.5*360/self.im.Qr)/180.*pi
        X = R*cos(THETA)
        Y = R*sin(THETA)
        l3 = self.line(0, 0, X, Y) # Line.3

        # raise Exception('before trim')

        # trim the lines first, because there is a chance the inner circle of rotor slot can intesect l1
        self.trim_l(l1,-EPS-self.im.Radius_OuterRotor, 0.5*self.im.Width_RotorSlotOpen)
        self.trim_l(l1,-EPS-self.im.Location_RotorBarCenter, 0.5*self.im.Width_RotorSlotOpen)


        if self.im.use_drop_shape_rotor_bar == True:
            # the inner rotor slot for drop shape rotor suggested by Gerada2011
            c4=self.circle(-self.im.Location_RotorBarCenter2, 0, self.im.Radius_of_RotorSlot2) # Circle.4

            l4 = self.line(-self.im.Location_RotorBarCenter-0.5*self.im.Radius_of_RotorSlot, c3.GetRadius(), -self.im.Location_RotorBarCenter2+0.5*self.im.Radius_of_RotorSlot2, c4.GetRadius())

            # Constraint to fix c4's center
            ref1 = c4.GetCenterVertex()
            ref2 = self.doc.CreateReferenceFromItem(ref1)
                # sketch.CreateMonoConstraint(u"distancefromxaxis", ref2)
                # sketch.GetItem(u"Distance From X Axis").SetProperty(u"Distance", 0)
            # Constrants to avoid moving of circle center
            # self.constraint_fixture_circle_center(c1)
            self.constraint_fixture_circle_center(c4) # Fixture constraint


                # # Constraint to fix c3's center
                # ref1 = c3.GetCenterVertex()
                # ref2 = self.doc.CreateReferenceFromItem(ref1)
                # sketch.CreateMonoConstraint(u"distancefromxaxis", ref2)
                # sketch.GetItem(u"Distance From X Axis").SetProperty(u"Distance", 0)

            ref1 = c4
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            ref3 = l4
            ref4 = self.doc.CreateReferenceFromItem(ref3)
            sketch.CreateBiConstraint("tangency", ref2, ref4)
            # sketch.GetItem(u"Vertex.15").SetProperty(u"Y", 2.345385)
            # sketch.GetItem(u"Vertex.16").SetProperty(u"Y", 2.345385)


            ref1 = c3
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            ref3 = l4
            ref4 = self.doc.CreateReferenceFromItem(ref3)
            sketch.CreateBiConstraint("tangency", ref2, ref4)
            # sketch.GetItem(u"Vertex.7").SetProperty(u"Y", 2.993642)
            # sketch.GetItem(u"Vertex.8").SetProperty(u"Y", 2.993642)

            # we won't need this fixture constraint anymore
            self.doc.GetSelection().Clear()
            self.doc.GetSelection().Add(sketch.GetItem("Fixture")) # Fixture.1
            self.doc.GetSelection().Delete()

            if -self.im.Location_RotorBarCenter + self.im.Radius_of_RotorSlot > -self.im.Location_RotorBarCenter2 - self.im.Radius_of_RotorSlot2:
                'Two circles have intersections.'
                # delete c4
                self.doc.GetSelection().Clear()
                self.doc.GetSelection().Add(c4)
                self.doc.GetSelection().Delete()

                # raise Exception('Is c4 still there? If c4 is not deleted, you will get Unexpected Part Number error afterwards. So fix it here now.')


        # Trim the Sketch Object!
        arc2 = self.trim_c(c2,0, c2.GetRadius()) # or self.trim_c(c2,0, self.im.Radius_OuterRotor)

        X = float(c3.GetCenterVertex().GetX()) # the returned type is long, which will cause failing to trim. Convert to float. Or simply add EPS to it, so Python will do the conversion for you.
        arc3 = self.trim_c(c3, X, -c3.GetRadius()+1e-6) # or self.trim_c(c3,-self.im.Location_RotorBarCenter, -self.im.Radius_of_RotorSlot)


        # åˆ°è¿™ä¸€æ­\ï¼Œä¸èƒ½å…ˆåˆ‡Line.3ï¼Œå¦åˆ™Circle.1ä¼šå\½åƒä¸å­˜åœ¨ä¸€æ ·ï¼Œå¯¼è‡´Line.3æ•´æ ¹æ¶ˆå¤±ï¼æ‰€ä»\ï¼Œå…ˆåˆ‡Circle.1
        arc1 = self.trim_c(c1, 0, c1.GetRadius()) # or self.trim_c(c1, 0, self.im.Radius_Shaft)

        self.trim_l(l2,-0.5*self.im.Radius_Shaft, 0)

        # 上面说的“导致Line.3整根消失！”的原因是直接剪在了原点(0,0)上，所以把整根线都删掉了，稍微偏一点(-0.1,0.1)操作即可
        # ä¸Šé¢è¯´çš„â€œå¯¼è‡´Line.3æ•´æ ¹æ¶ˆå¤±ï¼â€çš„åŽŸå› æ˜¯ç›´æŽ\å‰ªåœ¨äº†åŽŸç‚¹(0,0)ä¸Šï¼Œæ‰€ä»\æŠŠæ•´æ ¹çº¿éƒ½åˆ æŽ‰äº†ï¼Œç¨å¾®åä¸€ç‚¹(-0.1,0.1)æ“ä½œå³å¯
            # self.doc.GetSketchManager().SketchTrim(0,0) # BUG - delete the whole Line.3
        self.trim_l(l3,-EPS, EPS)

        if self.im.use_drop_shape_rotor_bar == True:

            if -self.im.Location_RotorBarCenter + self.im.Radius_of_RotorSlot > -self.im.Location_RotorBarCenter2 - self.im.Radius_of_RotorSlot2:
                'Two circles have intersections.'
                # re-create c4
                c4=self.circle(-self.im.Location_RotorBarCenter2, 0, self.im.Radius_of_RotorSlot2) # Circle.4

            arc4 = self.trim_c(c4, -self.im.Location_RotorBarCenter2, -c4.GetRadius())

            self.trim_l(l4, l4.GetStartVertex().GetX()+EPS, l4.GetStartVertex().GetY()+EPS)
            self.trim_l(l4, l4.GetEndVertex().GetX()-EPS, l4.GetEndVertex().GetY()-EPS)

            # Mirror and Duplicate
            region = self.create_region(["Arc","Arc.2","Arc.3","Arc.4","Line","Line.2","Line.3","Line.4"])
        else:
            # Mirror and Duplicate
            region = self.create_region(["Arc","Arc.2","Arc.3","Line","Line.2","Line.3"])


        # This is necessary if you want to do MODEL_ROTATE
        if self.im.MODEL_ROTATE:
            self.doc.GetSelection().Clear()
            self.doc.GetSelection().Add(sketch.GetItem("Horizontality"))
            self.doc.GetSelection().Add(sketch.GetItem("Horizontality.2"))
            self.doc.GetSelection().Add(sketch.GetItem("Tangency"))
            self.doc.GetSelection().Add(sketch.GetItem("Tangency.2"))
            self.doc.GetSelection().Add(sketch.GetItem("Coincident"))
            self.doc.GetSelection().Add(sketch.GetItem("Coincident.2"))
            self.doc.GetSelection().Delete()

            self.im.list_rotorCore_vertex_names = [ arc1.GetCenterVertex().GetName(),
                                                    arc2.GetCenterVertex().GetName(),
                                                    arc3.GetCenterVertex().GetName(),
                                                    arc4.GetCenterVertex().GetName(),
                                                    l1.GetStartVertex().GetName(),
                                                    l1.GetEndVertex().GetName(),
                                                    l2.GetStartVertex().GetName(),
                                                    l2.GetEndVertex().GetName(),
                                                    l3.GetStartVertex().GetName(),
                                                    l3.GetEndVertex().GetName(),
                                                    l4.GetStartVertex().GetName(),
                                                    l4.GetEndVertex().GetName()]
            self.im.list_rotorCore_vertex_names = list(dict.fromkeys(self.im.list_rotorCore_vertex_names).keys()) #https://stackoverflow.com/questions/7961363/removing-duplicates-in-lists?page=1&tab=votes#tab-top

        self.region_mirror_copy(region, l3)
        self.region_circular_pattern_360_origin(region, self.im.Qr)


        sketch.CloseSketch()
        # sketch.SetProperty(u"Visible", 0)
    def plot_statorCore(self, name=None):
        if name == None:
            self.SketchName="Stator Core"
        else:
            self.SketchName=name
        SketchName = self.SketchName
        sketch = self.create_sketch(SketchName, "#E8B5CE")

        self.im.Radius_InnerStator = self.im.Length_AirGap + self.im.Radius_OuterRotor
        sketch.CreateVertex(0., 0.)
        c1=self.circle(0., 0., self.im.Radius_InnerStator)
        c2=self.circle(0., 0., self.im.Radius_InnerStatorYoke)
        c3=self.circle(0., 0., self.im.Radius_OuterStatorYoke)
        c4=self.circle(0., 0., self.im.Radius_InnerStator + self.im.Width_StatorTeethHeadThickness)

        l1=self.line(-self.im.Radius_OuterStatorYoke, 0., 0., 0.)
        l2=self.line(-0.5*(self.im.Radius_OuterStatorYoke+self.im.Radius_InnerStatorYoke), 0.5*self.im.Width_StatorTeethBody, \
                -(self.im.Radius_InnerStator + (self.im.Width_StatorTeethHeadThickness+self.im.Width_StatorTeethNeck)), 0.5*self.im.Width_StatorTeethBody) # legacy 1. NEW: add neck here. Approximation is adopted: the accurate results should be (self.im.Width_StatorTeethHeadThickness+self.im.Width_StatorTeethNeck) * cos(6 deg) or so. 6 deg = 360/24/2 - 3/2

        R = self.im.Radius_OuterStatorYoke
        THETA = (180-0.5*360/self.im.Qs)/180.*pi
        X = R*cos(THETA)
        Y = R*sin(THETA)
        l3 = StatorCore_Line_3 = self.line(0., 0., X, Y) # Line.3

        R = self.im.Radius_InnerStator + self.im.Width_StatorTeethHeadThickness
        THETA = (180.-(0.5*360./self.im.Qs-0.5*self.im.Angle_StatorSlotOpen))/180.*pi # BUG is found here, 3 is instead used for Angle_StatorSlotOpen
        X = R*cos(THETA)
        Y = R*sin(THETA)
        l4 = self.line(0., 0., X, Y) # Line.4


        # Trim for Arcs

        # 为了避免数Arc，我们应该在Circle.4时期就把Circle.4上的小片段优先剪了！
        arc4 = self.trim_c(c4, X+1e-6, Y+1e-6) # Arc.1

        self.trim_a(arc4, 0., self.im.Radius_InnerStator + self.im.Width_StatorTeethHeadThickness)
        self.trim_a(arc4, -(self.im.Radius_InnerStator + self.im.Width_StatorTeethHeadThickness), 1e-6)

        arc1 = self.trim_c(c1, 0., self.im.Radius_InnerStator) # self.trim_c(c1,0, self.im.Radius_InnerStator)

        R = self.im.Radius_InnerStator
        THETA = (180-0.5*360./self.im.Qs)/180.*pi
        X = R*cos(THETA)
        Y = R*sin(THETA)
        self.trim_a(arc1, X, Y-1e-6)

        # 为了避免数Arc，我们应该在Circle.2时期就把Circle.2上的小片段优先剪了！
        arc2 = self.trim_c(c2,-self.im.Radius_InnerStatorYoke, 0.25*self.im.Width_StatorTeethBody)

        self.trim_a(arc2, 0., self.im.Radius_InnerStatorYoke)

        self.trim_c(c3,0., self.im.Radius_OuterStatorYoke)

        self.trim_l(l1,-0.1, 0.)


        self.trim_l(l2,1e-6 -0.5*(self.im.Radius_OuterStatorYoke+self.im.Radius_InnerStatorYoke), 0.5*self.im.Width_StatorTeethBody)
            # self.trim_l(l2,1e-6 -1e-6 -(self.im.Radius_InnerStator + 0.5*(self.im.Width_StatorTeethHeadThickness+self.im.Width_StatorTeethNeck)), 0.5*self.im.Width_StatorTeethBody)

        self.trim_l(l3,-1e-6, 1e-6)

        self.trim_l(l4,-1e-6, 1e-6) # float number

        # This error is caused because Y=2 is an integer!
        # 2019-01-15 17:48:55,243 - population - ERROR - The drawing is terminated. Please check whether the specified bounds are proper.
        # Traceback (most recent call last):
        #   File "D:/OneDrive - UW-Madison/c/codes/population.py", line 816, in draw_jmag_model
        #     d.plot_statorCore(u"Stator Core")
        #   File "D:/OneDrive - UW-Madison/c/codes/population.py", line 4066, in plot_statorCore
        #     l5_start_vertex = sketch.CreateVertex(X, Y)
        #   File "<string>", line 43, in CreateVertex
        # ValueError: Error calling the remote method

        # Similaryly, this error is caused because self.l5_start_vertex_y=2 is an integer!
        # 2019-01-15 19:00:33,286 - population - ERROR - The drawing is terminated. Please check whether the specified bounds are proper.
        # Traceback (most recent call last):
        #   File "D:/OneDrive - UW-Madison/c/codes/population.py", line 818, in draw_jmag_model
        #     d.plot_statorCore(u"Stator Core")
        #   File "D:/OneDrive - UW-Madison/c/codes/population.py", line 4106, in plot_statorCore
        #     raise e
        # ValueError: Error calling the remote method

        # we forget to plot the neck of stator tooth
        self.l5_start_vertex_x = l2.GetEndVertex().GetX()        # used later for plot_coil()
        self.l5_start_vertex_y = float(l2.GetEndVertex().GetY()) # used later for plot_coil()        
        X = arc4.GetStartVertex().GetX()
        Y = arc4.GetStartVertex().GetY()
        try:
            l5 = self.line( self.l5_start_vertex_x, self.l5_start_vertex_y, X, Y) # Line.5
        except Exception as e:
            print(l2.GetEndVertex().GetX(), l2.GetEndVertex().GetY(), arc4.GetStartVertex().GetX(), arc4.GetStartVertex().GetY())
            logger = logging.getLogger(__name__)
            logger.error('Draw Line.5 for Stator Core Failed, because integer cannot be passed as paramters to JMAG remote mathod. It is wierd, but just follow the rule.', exc_info=True)
            raise e


        self.doc.GetSelection().Clear()
        self.doc.GetSelection().Add(sketch.GetItem(arc4.GetName())) # arc4 corresponds to u"Arc".
        self.doc.GetSelection().Delete()

            # self.trim_l(l2, l2_end_vertex.GetX()-1e-3, l2_end_vertex.GetY()+1e-3) # convert to float (it returns long int) # legacy 3


        region = self.create_region(["Arc.2","Arc.3","Arc.4","Line","Line.2","Line.3","Line.4","Line.5"])

        self.region_mirror_copy(region, l1)

        self.region_circular_pattern_360_origin(region, self.im.Qs)

        ''' Constraints Stator Core '''
        if 0:
            ref1 = sketch.GetItem("Line.2")
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            ref3 = sketch.GetItem("Line")
            ref4 = self.doc.CreateReferenceFromItem(ref3)
            sketch.CreateBiConstraint("distance", ref2, ref4)
            sketch.GetItem("Distance").SetProperty("Distance", 0.5*self.im.Width_StatorTeethBody)

        sketch.CloseSketch()
        # sketch.SetProperty(u"Visible", 0)
    def plot_cage(self, name=None):
        if name == None:
            self.SketchName="Cage"
        else:
            self.SketchName=name
        SketchName = self.SketchName
        sketch = self.create_sketch(SketchName, "#8D9440")

        c3=self.circle(-self.im.Location_RotorBarCenter, 0., self.im.Radius_of_RotorSlot)

        if self.im.use_drop_shape_rotor_bar == True:
            # the inner rotor slot for drop shape rotor suggested by Gerada2011
            c4=self.circle(-self.im.Location_RotorBarCenter2, 0., self.im.Radius_of_RotorSlot2) # Circle.4
            # Constraint to fix c4's center
            ref1 = c4.GetCenterVertex()
            ref2 = self.doc.CreateReferenceFromItem(ref1)
                # sketch.CreateMonoConstraint(u"distancefromxaxis", ref2)
                # sketch.GetItem(u"Distance From X Axis").SetProperty(u"Distance", 0)
            # Constrants to avoid moving of circle center
            self.constraint_fixture_circle_center(c4)


            l41 = self.line(-self.im.Location_RotorBarCenter-0.5*self.im.Radius_of_RotorSlot, c3.GetRadius(), -self.im.Location_RotorBarCenter2+0.5*self.im.Radius_of_RotorSlot2, c4.GetRadius())
            l42 = self.line(-self.im.Location_RotorBarCenter-0.5*self.im.Radius_of_RotorSlot, -c3.GetRadius(), -self.im.Location_RotorBarCenter2+0.5*self.im.Radius_of_RotorSlot2, -c4.GetRadius())

            ref1 = c4
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            ref3 = l41
            ref4 = self.doc.CreateReferenceFromItem(ref3)
            sketch.CreateBiConstraint("tangency", ref2, ref4)
            ref1 = c3
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            ref3 = l41
            ref4 = self.doc.CreateReferenceFromItem(ref3)
            sketch.CreateBiConstraint("tangency", ref2, ref4)

            ref1 = c4
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            ref3 = l42
            ref4 = self.doc.CreateReferenceFromItem(ref3)
            sketch.CreateBiConstraint("tangency", ref2, ref4)
            ref1 = c3
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            ref3 = l42
            ref4 = self.doc.CreateReferenceFromItem(ref3)
            sketch.CreateBiConstraint("tangency", ref2, ref4)

            if -self.im.Location_RotorBarCenter + self.im.Radius_of_RotorSlot > -self.im.Location_RotorBarCenter2 - self.im.Radius_of_RotorSlot2:
                'Two circles have intersections.'
                # delete c4
                self.doc.GetSelection().Clear()
                self.doc.GetSelection().Add(c4)
                self.doc.GetSelection().Delete()
            arc3 = self.trim_c(c3, c3.GetCenterVertex().GetX()+EPS+c3.GetRadius(), 0) # make sure it is float number

            if -self.im.Location_RotorBarCenter + self.im.Radius_of_RotorSlot > -self.im.Location_RotorBarCenter2 - self.im.Radius_of_RotorSlot2:
                'Two circles have intersections.'
                # re-create c4
                c4=self.circle(-self.im.Location_RotorBarCenter2, 0, self.im.Radius_of_RotorSlot2) # Circle.4
            arc4 = self.trim_c(c4, -self.im.Location_RotorBarCenter2-c4.GetRadius(), 0)

            self.trim_l(l41, l41.GetStartVertex().GetX()+EPS, l41.GetStartVertex().GetY()+EPS)
            self.trim_l(l41, l41.GetEndVertex().GetX()-EPS, l41.GetEndVertex().GetY()-EPS)

            self.trim_l(l42, l42.GetStartVertex().GetX()+EPS, l42.GetStartVertex().GetY()+EPS)
            self.trim_l(l42, l42.GetEndVertex().GetX()-EPS, l42.GetEndVertex().GetY()-EPS)

            region = self.create_region(["Arc","Arc.2","Line","Line.2"])
        else:
            print('Round bar rotor slot is being plotted\n'*3)
            region = self.create_region(["Circle"])

        if self.im.MODEL_ROTATE:
            self.doc.GetSelection().Clear()
            self.doc.GetSelection().Add(sketch.GetItem("Fixture"))
            self.doc.GetSelection().Add(sketch.GetItem("Tangency"))
            self.doc.GetSelection().Add(sketch.GetItem("Tangency.2"))
            self.doc.GetSelection().Add(sketch.GetItem("Tangency.3"))
            self.doc.GetSelection().Add(sketch.GetItem("Tangency.4"))
            self.doc.GetSelection().Delete()

            self.im.list_rotorCage_vertex_names = [  arc3.GetCenterVertex().GetName(),
                                                arc4.GetCenterVertex().GetName(),
                                                l41.GetStartVertex().GetName(),
                                                l41.GetEndVertex().GetName(),
                                                l42.GetStartVertex().GetName(),
                                                l42.GetEndVertex().GetName()]
            self.im.list_rotorCage_vertex_names = list(dict.fromkeys(self.im.list_rotorCage_vertex_names).keys())


        self.region_circular_pattern_360_origin(region, self.im.Qr)
        # raise

        sketch.CloseSketch()
    def plot_coil(self, name=None):
        if name == None:
            self.SketchName="Coil"
        else:
            self.SketchName=name
        SketchName = self.SketchName
        sketch = self.create_sketch(SketchName, "#EC9787")


        # 方案一 画中线，然后镜像到对面，得到另外两个点！
        # 方案二 为什么不直接画圆然后Trim哦！
        if 1:
            self.im.Radius_InnerStator = self.im.Length_AirGap + self.im.Radius_OuterRotor
            # sketch.CreateVertex(0, 0)
            c1=self.circle(0, 0, self.im.Radius_InnerStatorYoke)
            c2=self.circle(0, 0, self.im.Radius_InnerStator + self.im.Width_StatorTeethHeadThickness)


            l2l5 = self.line(-0.5*(self.im.Radius_OuterStatorYoke+self.im.Radius_InnerStatorYoke), 0.5*self.im.Width_StatorTeethBody, \
                             self.l5_start_vertex_x, self.l5_start_vertex_y)

            # l2 = self.line(-0.5*(self.im.Radius_OuterStatorYoke+self.im.Radius_InnerStatorYoke), 0.5*self.im.Width_StatorTeethBody, \
            #                -(self.im.Radius_InnerStator + 1.0*self.im.Width_StatorTeethHeadThickness) - self.im.Width_StatorTeethNeck, 0.5*self.im.Width_StatorTeethBody)

            R = self.im.Radius_InnerStatorYoke
            THETA = (180-0.5*360/self.im.Qs)/180.*pi
            X = R*cos(THETA)
            Y = R*sin(THETA)
            l3 = self.line(0, 0, X, Y) # Line.3


            # Trim Lines
            self.trim_l(l2l5,1e-6 -0.5*(self.im.Radius_OuterStatorYoke+self.im.Radius_InnerStatorYoke), 0.5*self.im.Width_StatorTeethBody)
                # self.trim_l(l2l5,-1e-6 -(self.im.Radius_InnerStator + 0.5*self.im.Width_StatorTeethHeadThickness), 0.5*self.im.Width_StatorTeethBody)
            self.trim_l(l3,-1e-6, 1e-6)


            l6 = self.line(self.l5_start_vertex_x, self.l5_start_vertex_y, \
                            l3.GetStartVertex().GetX(), l3.GetStartVertex().GetY())


            from CrossSectStator import get_area_polygon
            # print(  sketch.GetItem("Line").GetStartVertex().GetX(), sketch.GetItem("Line").GetStartVertex().GetY(),
            #         sketch.GetItem("Line").GetEndVertex().GetX(), sketch.GetItem("Line").GetEndVertex().GetY(),
            #         sketch.GetItem("Line.2").GetStartVertex().GetX(), sketch.GetItem("Line.2").GetStartVertex().GetY(),
            #         sketch.GetItem("Line.2").GetEndVertex().GetX(), sketch.GetItem("Line.2").GetEndVertex().GetY()
            #         )
            self.mm2_slot_area = 2 * get_area_polygon(
                                                        [sketch.GetItem("Line").GetStartVertex().GetX(), sketch.GetItem("Line").GetStartVertex().GetY()],
                                                        [sketch.GetItem("Line").GetEndVertex().GetX(), sketch.GetItem("Line").GetEndVertex().GetY()],
                                                        [sketch.GetItem("Line.2").GetStartVertex().GetX(), sketch.GetItem("Line.2").GetStartVertex().GetY()],
                                                        [sketch.GetItem("Line.2").GetEndVertex().GetX(), sketch.GetItem("Line.2").GetEndVertex().GetY()]
                                                        )
            print('Stator Slot Area: %g mm^2'%(self.mm2_slot_area))
            # raise KeyboardInterrupt

            # Trim Circles
            self.trim_c(c1,0, self.im.Radius_InnerStatorYoke)
            self.trim_c(c2,0, self.im.Radius_InnerStator + self.im.Width_StatorTeethHeadThickness)


        # Mirror and Duplicate
        region = self.create_region(["Line","Line.2","Line.3","Arc"])

        # region_mirror_pattern_which_is_ItemObject = self.region_mirror_copy(region, l3, merge=False)
        # region_in_the_mirror = RegionItem(region_mirror_pattern_which_is_ItemObject)
        region_in_the_mirror = self.region_mirror_copy(region, l3, merge=False)
        # print region_in_the_mirror # it's ItemObject rather than RegionObject! We are fucked, because we cannot generate reference for ItemObject. We can generate ref for RegionObject.

        self.region_circular_pattern_360_origin(region, self.im.Qs, merge=False, do_you_have_region_in_the_mirror=True)

        ''' Constraints Coil '''
        if 0:
            ref1 = sketch.GetItem("Line") # or l1
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            sketch.CreateMonoConstraint("horizontality", ref2)

            ref1 = sketch.GetItem("Line").GetEndVertex()
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            sketch.CreateMonoConstraint("distancefromxaxis", ref2)
            sketch.GetItem("Distance From X Axis").SetProperty("Distance", 0.5*self.im.Width_StatorTeethBody)

        sketch.CloseSketch()
    def plot_airWithinRotorSlots(self, name=None):
        if name == None:
            self.SketchName="Air Within Rotor Slots"
        else:
            self.SketchName=name
        SketchName = self.SketchName
        sketch = self.create_sketch(SketchName, "#BFD641")


        c1=self.circle(0,0,self.im.Radius_OuterRotor)
        c2=self.circle(-self.im.Location_RotorBarCenter,0,self.im.Radius_of_RotorSlot)

        l1=self.line(-5.5-self.im.Radius_OuterRotor, 0.5*self.im.Width_RotorSlotOpen, -self.im.Location_RotorBarCenter, 0.5*self.im.Width_RotorSlotOpen)
        l2=self.line(-5.5-self.im.Radius_OuterRotor,                       0, -self.im.Location_RotorBarCenter,                       0)

        self.trim_l(l1,-5-self.im.Radius_OuterRotor, 0.5*self.im.Width_RotorSlotOpen)
        self.trim_l(l2,-5-self.im.Radius_OuterRotor, 0)

        self.trim_l(l1,-1e-6-self.im.Location_RotorBarCenter, 0.5*self.im.Width_RotorSlotOpen)
        self.trim_l(l2,-1e-6-self.im.Location_RotorBarCenter, 0)

        arc1 = self.trim_c(c1,0,c1.GetRadius())
        arc2 = self.trim_c(c2,0,c2.GetRadius())


        if self.im.MODEL_ROTATE:
            self.doc.GetSelection().Clear()
            self.doc.GetSelection().Add(sketch.GetItem("Coincident"))
            self.doc.GetSelection().Add(sketch.GetItem("Fixture"))
            self.doc.GetSelection().Delete()

            self.im.list_rotorAirWithin_vertex_names = [ arc1.GetCenterVertex().GetName(),
                                                    arc2.GetCenterVertex().GetName(),
                                                    l1.GetStartVertex().GetName(),
                                                    l1.GetEndVertex().GetName(),
                                                    l2.GetStartVertex().GetName(),
                                                    l2.GetEndVertex().GetName()]
            self.im.list_rotorAirWithin_vertex_names = list(dict.fromkeys(self.im.list_rotorAirWithin_vertex_names).keys())

        region = self.create_region(["Arc.2","Arc","Line.2","Line"])

        self.region_mirror_copy(region, l2)
        self.region_circular_pattern_360_origin(region, self.im.Qr)

        ''' Constraint Air within Rotor Slots '''
        if 0:
            ref1 = sketch.GetItem("Vertex.2")
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            sketch.CreateMonoConstraint("distancefromyaxis", ref2)
            sketch.GetItem("Distance From Y Axis").SetProperty("Distance", self.im.Location_RotorBarCenter)
            sketch.GetItem("Distance From Y Axis").SetProperty("Name", "_Air_LocationBar")

            ref1 = sketch.GetItem("Arc.2")
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            sketch.CreateMonoConstraint("radius", ref2)
            sketch.GetItem("Radius/Diameter").SetProperty("Radius", self.im.Radius_of_RotorSlot)
            sketch.GetItem("Radius/Diameter").SetProperty("Name", "_Air_RadiusOfRotorSlot")

            ref1 = sketch.GetItem("Vertex")
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            sketch.CreateMonoConstraint("fixture", ref2)

            ref1 = sketch.GetItem("Arc")
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            sketch.CreateMonoConstraint("radius", ref2)
            sketch.GetItem("Radius/Diameter.2").SetProperty("Radius", self.im.Radius_OuterRotor)
            sketch.GetItem("Radius/Diameter.2").SetProperty("Name", "_Air_RadiusOuterRotor")

            ref1 = sketch.GetItem("Line")
            ref2 = self.doc.CreateReferenceFromItem(ref1)
            ref3 = sketch.GetItem("Line.2")
            ref4 = self.doc.CreateReferenceFromItem(ref3)
            sketch.CreateBiConstraint("distance", ref2, ref4)
            sketch.GetItem("Distance").SetProperty("Distance", 0.5*self.im.Width_RotorSlotOpen)
            sketch.GetItem("Distance").SetProperty("Name", "_Air_HalfSlotOpenRotor")

        sketch.CloseSketch()
    def plot_airHugRotor(self, name=None):
        if name == None:
            self.SketchName="Air Hug Rotor"
        else:
            self.SketchName=name
        SketchName = self.SketchName
        sketch = self.create_sketch(SketchName, "#4F84C4")

        R = self.im.Radius_OuterRotor+0.5*self.im.Length_AirGap
        c1=self.circle(0,0,self.im.Radius_OuterRotor)
        c2=self.circle(0,0,R)

        l1=self.line(0,-R,0,R)
        self.trim_l(l1,0,0)

        self.trim_c(c1,c1.GetRadius(),0)
        self.trim_c(c2,c2.GetRadius(),0)


        region = self.create_region(["Arc","Arc.2","Line","Line.2"])

        self.region_mirror_copy(region, symmetry_type=3) # y-axis æ˜¯å¯¹ç§°è½´
        # sketch.CreateRegionMirrorCopy()
        # sketch.GetItem(u"Region Mirror Copy").SetProperty(u"Merge", 1)
        # refarray = [0 for i in range(1)]
        # refarray[0] = u"faceregion(TRegionItem126)"
        # sketch.GetItem(u"Region Mirror Copy").SetProperty(u"Region", refarray)
        # sketch.GetItem(u"Region Mirror Copy").SetProperty(u"SymmetryType", 3)

        sketch.CloseSketch()

def add_M1xSteel(app, dir_parent, steel_name="M-19 Steel Gauge-29"):

    if '19' in steel_name:
        BH = np.loadtxt(dir_parent + 'codes3/M-19-Steel-BH-Curve-afterJMAGsmooth.BH', unpack=True, usecols=(0,1)) # after JMAG smooth, it beomces HB rather than BH
    elif '15' in steel_name:
        BH = np.loadtxt(dir_parent + './M-15-Steel-BH-Curve.txt', unpack=True, usecols=(0,1))


    app.GetMaterialLibrary().CreateCustomMaterial(steel_name, "Custom Materials")
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("Density", 7.85)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("MagneticSteelPermeabilityType", 2)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("CoerciveForce", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).GetTable("BhTable").SetName("Untitled")

    refarray = BH.T.tolist()

    app.GetMaterialLibrary().GetUserMaterial(steel_name).GetTable("BhTable").SetTable(refarray)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("DemagnetizationCoerciveForce", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("MagnetizationSaturated", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("MagnetizationSaturated2", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("YoungModulus", 210000)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("ShearModulus", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("YoungModulusX", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("YoungModulusY", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("YoungModulusZ", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("ShearModulusXY", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("ShearModulusYZ", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("ShearModulusZX", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G11", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G12", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G13", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G14", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G15", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G16", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G22", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G23", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G24", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G25", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G26", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G33", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G34", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G35", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G36", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G44", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G45", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G46", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G55", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G56", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("G66", 0)

    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("MagnetizationSaturated2", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("MagnetizationSaturatedMakerValue", 0)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("Loss_Type", 1)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("LossConstantKhX", 143.)
    app.GetMaterialLibrary().GetUserMaterial(steel_name).SetValue("LossConstantKeX", 0.530)

def add_Arnon5(app, dir_parent):
    app.GetMaterialLibrary().CreateCustomMaterial("Arnon5-final", "Custom Materials")
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("Density", 7.85)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("MagneticSteelPermeabilityType", 2)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("CoerciveForce", 0)
    # app.GetMaterialLibrary().GetUserMaterial(u"Arnon5-final").GetTable("BhTable").SetName(u"SmoothZeroPointOne")

    BH = np.loadtxt(dir_parent + 'Arnon5/Arnon5-final.txt', unpack=True, usecols=(0,1))
    refarray = BH.T.tolist()

    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").GetTable("BhTable").SetTable(refarray)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("DemagnetizationCoerciveForce", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("MagnetizationSaturated", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("MagnetizationSaturated2", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("ExtrapolationMethod", 1)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("YoungModulus", 210000)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("ShearModulus", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("YoungModulusX", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("YoungModulusY", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("YoungModulusZ", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("ShearModulusXY", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("ShearModulusYZ", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("ShearModulusZX", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G11", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G12", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G13", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G14", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G15", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G16", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G22", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G23", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G24", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G25", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G26", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G33", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G34", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G35", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G36", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G44", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G45", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G46", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G55", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G56", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("G66", 0)

    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("MagnetizationSaturated2", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("MagnetizationSaturatedMakerValue", 0)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("Loss_Type", 1)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("LossConstantKhX", 186.6)
    app.GetMaterialLibrary().GetUserMaterial("Arnon5-final").SetValue("LossConstantKeX", 0.07324)

def add_Arnon7(app, dir_parent):
    pass

if __name__ == '__main__':
    pass

