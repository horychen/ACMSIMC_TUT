from time import time as clock_time
from pylab import plt, np
import os
import win32com.client
import logging
import utility
import pyrhonen_procedure_as_function
import population
import FEMM_Solver

# BPMSM codes
# import bearingless_spmsm_design
# import JMAG

class swarm_data_container(object):
    def __init__(self, swarm_data_raw, fea_config_dict, bound_filter):
        self.swarm_data_raw = swarm_data_raw
        self.fea_config_dict = fea_config_dict

        # x, f(x)
        self.swarm_data_xf = []
        self.machine_data = []
        self.rated_data = []
        for el in self.swarm_data_raw:
            design_parameters_denorm = [float(x) for x in el[5].split(',')]
            loc1 = el[2].find('f1')
            loc2 = el[2].find('f2')
            loc3 = el[2].find('f3')
            f1 = float(el[2][loc1+3:loc2-1])
            f2 = float(el[2][loc2+3:loc3-1])
            f3 = float(el[2][loc3+3:])

            x_denorm = self.get_x_denorm_from_design_parameters(design_parameters_denorm, bound_filter)
            self.swarm_data_xf.append(x_denorm + [f1, f2, f3])
            # print(design_parameters_denorm, f1, f2, f3)

            self.machine_data.append([float(x) for x in el[3].split(',')])
            self.rated_data.append(  [float(x) for x in el[4].split(',')])

        self.number_of_free_variables = len(x_denorm)
        print('\tCount of individuals:', len(self.swarm_data_raw))

        self.l_OA = [el[-3] for el in self.swarm_data_xf]
        self.l_OB = [el[-2] for el in self.swarm_data_xf]
        self.l_OC = [el[-1] for el in self.swarm_data_xf]
        self.l_design_parameters = [el[:-3] for el in self.swarm_data_xf]

        # [power_factor, efficiency, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle]
        self.l_power_factor                     = [el[0] for el in self.machine_data]
        self.l_efficiency                       = [el[1] for el in self.machine_data]
        self.l_torque_average                   = [el[2] for el in self.machine_data]
        self.l_normalized_torque_ripple         = [el[3] for el in self.machine_data]
        self.l_ss_avg_force_magnitude           = [el[4] for el in self.machine_data]
        self.l_normalized_force_error_magnitude = [el[5] for el in self.machine_data]
        self.l_force_error_angle                = [el[6] for el in self.machine_data]

        self.l_rated_shaft_power                    = [el[0] for el in self.rated_data]
        self.l_rated_efficiency                     = [el[1] for el in self.rated_data]
        self.l_rated_total_loss                     = [el[2] for el in self.rated_data]
        self.l_rated_stator_copper_loss_along_stack = [el[3] for el in self.rated_data]
        self.l_rated_rotor_copper_loss_along_stack  = [el[4] for el in self.rated_data]
        self.l_stator_copper_loss_in_end_turn       = [el[5] for el in self.rated_data]
        self.l_rotor_copper_loss_in_end_turn        = [el[6] for el in self.rated_data]
        self.l_rated_iron_loss                      = [el[7] for el in self.rated_data]
        self.l_rated_windage_loss                   = [el[8] for el in self.rated_data]
        self.l_rated_rotor_volume                   = [el[9] for el in self.rated_data]
        self.l_rated_rotor_weight                   = [(V*8050*9.8) for V in self.l_rated_rotor_volume] # density of rotor is estimated to be that of steel of 8050 g/cm^3
        self.l_rated_stack_length                   = [el[10] for el in self.rated_data] # new!
        self.l_original_stack_length                = [el[11] for el in self.rated_data] # new!
        self.l_original_rotor_weight                = [el/rated*ori for el, ori, rated in zip(self.l_rated_rotor_weight, self.l_rated_stack_length, self.l_original_stack_length)]

        required_torque = 50e3 / (30000/60*2*np.pi)
        self.l_TRV = [required_torque/el for el in self.l_rated_rotor_volume]

    def get_list_y_data(self):

        list_y_data = [ self.l_TRV, ##self.l_rated_stack_length,
                        [100*el for el in self.l_OB], 
                        self.l_force_error_angle,
                        ]
        return list_y_data

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Utility
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def sensitivity_bar_charts(self):
        number_of_variant = self.fea_config_dict['local_sensitivity_analysis_number_of_variants'] + 1
        number_of_free_variables = self.number_of_free_variables

        from pylab import subplots, mpl, plt
        mpl.style.use('classic')
        mpl.rcParams['legend.fontsize'] = 12
        # mpl.rcParams['legend.family'] = 'Times New Roman'
        mpl.rcParams['font.family'] = ['Times New Roman']
        # mpl.rcParams['font.size'] = 15.0
        font = {'family' : 'Times New Roman', #'serif',
                'color' : 'darkblue',
                'weight' : 'normal',
                'size' : 14,}
        textfont = {'family' : 'Times New Roman', #'serif',
                    'color' : 'darkblue',
                    'weight' : 'normal',
                    'size' : 11.5,}

        fig, axeses = subplots(4, 2, sharex=True, dpi=150, figsize=(16*0.75, 8*0.75), facecolor='w', edgecolor='k', constrained_layout=True)
        ax_list = []
        for i in range(4):
            ax_list.extend(axeses[i].tolist())
        # O2_prototype_ax.plot(O2_prototype_data[1], 'o-', lw=0.75, alpha=0.5, label=r'$\delta$'         )
        # O2_prototype_ax.plot(O2_prototype_data[0], 'v-', lw=0.75, alpha=0.5, label=r'$b_{\rm tooth,s}$')
        # O2_prototype_ax.plot(O2_prototype_data[3], 's-', lw=0.75, alpha=0.5, label=r'$b_{\rm tooth,r}$')
        # O2_prototype_ax.plot(O2_prototype_data[5], '^-', lw=0.75, alpha=0.5, label=r'$w_{\rm open,s}$')
        # O2_prototype_ax.plot(O2_prototype_data[2], 'd-', lw=0.75, alpha=0.5, label=r'$w_{\rm open,r}$')
        # O2_prototype_ax.plot(O2_prototype_data[6], '*-', lw=0.75, alpha=0.5, label=r'$h_{\rm head,s}$')
        # O2_prototype_ax.plot(O2_prototype_data[4], 'X-', lw=0.75, alpha=0.5, label=r'$h_{\rm head,r}$')

        # Extract data
        free_param_list = [
        r'$L_g$',        
        r'$w_{st}$',     
        r'$w_{rt}$',     
        r'$\theta_{so}$',
        r'$w_{ro}$',
        r'$d_{so}$',
        r'$d_{ro}$']
        y_label_list = ['PF', r'$\eta$ [100%]', r'$T_{em} [N]$', r'$T_{rip}$ [100%]', r'$|F|$ [N]', r'$E_m$ [100%]', r'$E_a$ [deg]', 
                        r'$P_{Cu,s,JMAG}$', r'$P_{Cu,r,JMAG}$', r'$P_{Fe}$ [W]', r'$P_{eddy}$', r'$P_{hyst}$', r'$P_{Cu,s,FEMM}$', r'$P_{Cu,r,FEMM}$', 
                        r'Windage loss', r'Total loss']

        list_y_label = [r'$O_A$ [$\rm kNm/m^3$]', 
                         '$O_C$ [1]', 
                         'FRW [p.u.]',
                         '$E_a$ [deg]', 
                         '$O_B$ [%]', 
                         '$E_m$ [%]', 
                         r'$P_{\rm loss}$ [W]',
                         r'$T_{\rm rip}$ [%]',
                         # 'Rotor Weight [N]', #'Power Factor [1]',
                         ]
        list_y_data_max = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
        list_y_data_min = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

        list_y_data = [ [el/1e3 for el in self.l_OA], 
                        self.l_OC,
                        [F/W for W, F in zip(self.l_original_rotor_weight, self.l_ss_avg_force_magnitude)], # FRW
                        self.l_force_error_angle,
                        [100*el for el in self.l_OB], 
                        [100*el for el in self.l_normalized_force_error_magnitude],
                        self.l_rated_total_loss,
                        [100*el for el in self.l_normalized_torque_ripple],
                        # self.l_rated_rotor_weight,
                        ]
        for i in range(len(list_y_label)):
            ax = ax_list[i]
            y_data = list_y_data[i]
            y_value_reference = y_data[0]
            ax.plot(y_value_reference*np.ones(len(y_data)), '-k', alpha=1, zorder=10)
            y_data = y_data[1:]
            number_of_points_per_geometry_variable = len(y_data)/number_of_free_variables
            for idx, part_of_y_data in enumerate([y_data[int(number_of_points_per_geometry_variable*_)\
                                                        :int(number_of_points_per_geometry_variable*(_+1))]\
                                                        for _ in range(number_of_free_variables)]):
                if idx%2 == 0:
                    line_style = '--bo'
                else:
                    line_style = '--ro'
                ax.plot(list(range(len(y_data)))[int(number_of_points_per_geometry_variable*idx)\
                                                :int(number_of_points_per_geometry_variable*(idx+1))], 
                                                part_of_y_data, line_style, alpha=0.33)

            low, high = ax.get_ylim()
            # ax.legend()
            ax.grid()
            ax.set_ylabel(list_y_label[i], **font)
            ax.set_xlim([0,140])
            for j in range(number_of_free_variables):
                if j%2==0:
                    alpha = 0.05
                else:
                    alpha = 0.15
                ax.axvspan(j*number_of_variant-0.5, (j+1)*number_of_variant-0.5, facecolor='k', alpha=alpha)
                ax.text(0.33*number_of_free_variables+j*number_of_variant, high-(high-low)*0.125, free_param_list[j])

            if i == 0:
                ax.set_yticks([-24, -22, -20, -18, -16])
            list_y_data_max[i].append(max(y_data))
            list_y_data_min[i].append(min(y_data))

        ax_list[-2].set_xlabel('Count of design variant', **font)
        ax_list[-1].set_xlabel('Count of design variant', **font)
        fig.savefig(r'C:\Users\horyc\Desktop/'+ 'LSA_curves.png', dpi=300)
        # plt.show()
        return













        INDEX_TOTAL_LOSS = 15 + 4 # index of total loss in the machine_data list

        # ------------------------------------ Sensitivity Analysis Bar Chart Scripts

        # print next(self.get_list_objective_function())
        data_max = []
        data_min = []
        eta_at_50kW_max = []
        eta_at_50kW_min = []
        O1_max   = []
        O1_min   = []
        for ind, i in enumerate(list(range(7))+[INDEX_TOTAL_LOSS]):
            print('\n-----------', y_label_list[i])
            l = list(self.get_certain_objective_function(i))
            y = l
            print('ind=', ind, 'i=', i, 'len(y)=', len(y))

            data_max.append([])
            data_min.append([])

            for j in range(int(len(y)/number_of_variant)): # iterate design parameters
                y_vs_design_parameter = y[j*number_of_variant:(j+1)*number_of_variant]

                try:
                    # if j == 6:
                    ax_list[ind].plot(y_vs_design_parameter, 'o-', lw=0.75, label=str(j)+' '+param_list[j], alpha=0.5)
                except IndexError as e:
                    print('Check the length of y should be 7*(%d+1)=%d, or else you should remove the redundant results in swarm_data.txt (they are produced because of the interrupted/resumed script run.)'%(number_of_variant, 7*number_of_variant))
                    raise e
                print('\tj=', j, param_list[j], '\t\t Max-Min:', max(y_vs_design_parameter) - min(y_vs_design_parameter))

                data_max[ind].append(max(y_vs_design_parameter))
                data_min[ind].append(min(y_vs_design_parameter))            

            if i==1:
                ax_list[ind].legend(prop={'family':'Times New Roman'})
            ax_list[ind].grid()
            ax_list[ind].set_ylabel(y_label_list[i], **font)

        print('\nObjectives vs. geometry variables:')
        for ind, el in enumerate(data_max):
            print(ind, 'Max', el)
        print('\nObjectives vs. geometry variables:')
        for ind, el in enumerate(data_min):
            print(ind, 'Min', el)

        if self.reference_design is not None:
            print('\n-------------------- Here goes the reference design:')
            for el in self.reference_design[1:]:
                print(el, end=' ')
            self.reference_data = [float(el) for el in self.reference_design[3].split(',')]
            O2_ref = fobj_scalar(self.reference_data[2],
                                 self.reference_data[4],
                                 self.reference_data[3],
                                 self.reference_data[5],
                                 self.reference_data[6],
                                 self.reference_data[INDEX_TOTAL_LOSS],
                                 weights=use_weights('O2'), rotor_volume=self.rotor_volume, rotor_weight=self.rotor_weight)
            O1_ref = fobj_scalar(self.reference_data[2],
                                 self.reference_data[4],
                                 self.reference_data[3],
                                 self.reference_data[5],
                                 self.reference_data[6],
                                 self.reference_data[INDEX_TOTAL_LOSS],
                                 weights=use_weights('O1'), rotor_volume=self.rotor_volume, rotor_weight=self.rotor_weight)
        else:
            raise Exception('self.reference_design is None.')

        print('Objective function 1')
        O1 = fobj_list( list(self.get_certain_objective_function(2)), 
                        list(self.get_certain_objective_function(4)), 
                        list(self.get_certain_objective_function(3)), 
                        list(self.get_certain_objective_function(5)), 
                        list(self.get_certain_objective_function(6)), 
                        np.array(list(self.get_certain_objective_function(9))) + np.array(list(self.get_certain_objective_function(12))) + np.array(list(self.get_certain_objective_function(13))),
                        weights=use_weights('O1'), rotor_volume=self.rotor_volume, rotor_weight=self.rotor_weight)
        O1_max = []
        O1_min = []
        from pylab import figure
        O1_ax  = figure().gca()
        O2_prototype_data = []
        results_for_refining_bounds = {}
        results_for_refining_bounds['O1'] = []
        for j in range(int(len(O1)/number_of_variant)): # iterate design parameters
            O1_vs_design_parameter = O1[j*number_of_variant:(j+1)*number_of_variant]
            O2_prototype_data.append(O1_vs_design_parameter)

            O1_ax.plot(O1_vs_design_parameter, label=str(j)+' '+param_list[j], alpha=0.5)
            print('\t', j, param_list[j], '\t\t max O1 - min O1:', max(O1_vs_design_parameter) - min(O1_vs_design_parameter), '\t\t', end=' ')

            # narrow bounds (refine bounds)
            results_for_refining_bounds['O1'].append( [ind for ind, el in enumerate(O1_vs_design_parameter) if el < O1_ref*1.0] )
            print(results_for_refining_bounds['O1'][j]) #'<- to derive new original_bounds.'

            O1_max.append(max(O1_vs_design_parameter))
            O1_min.append(min(O1_vs_design_parameter))            
        O1_ax.legend()
        O1_ax.grid()
        O1_ax.set_ylabel('O1 [1]', **font)
        O1_ax.set_xlabel('Count of design variants', **font)

        # fig_prototype = figure(500, figsize=(10, 5), facecolor='w', edgecolor='k')
        # O2_prototype_ax = fig_prototype.gca()
        # O2_prototype_ax.plot(list(range(-1, 22)), O1_ref*np.ones(23), 'k--', label='Reference design')
        # O2_prototype_ax.plot(O2_prototype_data[1], 'o-', lw=0.75, alpha=0.5, label=r'$L_g$')
        # O2_prototype_ax.plot(O2_prototype_data[0], 'v-', lw=0.75, alpha=0.5, label=r'$w_{st}$')
        # O2_prototype_ax.plot(O2_prototype_data[3], 's-', lw=0.75, alpha=0.5, label=r'$w_{rt}$')
        # O2_prototype_ax.plot(O2_prototype_data[5], '^-', lw=0.75, alpha=0.5, label=r'$\theta_{so}$')
        # O2_prototype_ax.plot(O2_prototype_data[2], 'd-', lw=0.75, alpha=0.5, label=r'$w_{ro}$')
        # O2_prototype_ax.plot(O2_prototype_data[6], '*-', lw=0.75, alpha=0.5, label=r'$d_{so}$')
        # O2_prototype_ax.plot(O2_prototype_data[4], 'X-', lw=0.75, alpha=0.5, label=r'$d_{ro}$')
        # O2_prototype_ax.legend()
        # O2_prototype_ax.set_ylabel('$O_2(x)$ [1]', **font)

        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        # O2
        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        print('Objective function 2')
        O2 = fobj_list( list(self.get_certain_objective_function(2)), 
                        list(self.get_certain_objective_function(4)), 
                        list(self.get_certain_objective_function(3)), 
                        list(self.get_certain_objective_function(5)), 
                        list(self.get_certain_objective_function(6)), 
                        np.array(list(self.get_certain_objective_function(9))) + np.array(list(self.get_certain_objective_function(12))) + np.array(list(self.get_certain_objective_function(13))),
                        weights=use_weights('O2'), rotor_volume=self.rotor_volume, rotor_weight=self.rotor_weight )
        O2_max = []
        O2_min = []
        O2_ax  = figure().gca()
        O2_ecce_data = []
        results_for_refining_bounds['O2'] = []
        for j in range(int(len(O2)/number_of_variant)): # iterate design parameters: range(7)
            O2_vs_design_parameter = O2[j*number_of_variant:(j+1)*number_of_variant]
            O2_ecce_data.append(O2_vs_design_parameter)

            # narrow bounds (refine bounds)
            O2_ax.plot(O2_vs_design_parameter, 'o-', label=str(j)+' '+param_list[j], alpha=0.5)
            print('\t', j, param_list[j], '\t\t max O2 - min O2:', max(O2_vs_design_parameter) - min(O2_vs_design_parameter), '\t\t', end=' ')
            results_for_refining_bounds['O2'].append( [ind for ind, el in enumerate(O2_vs_design_parameter) if el < O2_ref*1.0] )
            print(results_for_refining_bounds['O2'][j]) #'<- to derive new original_bounds.'

            O2_max.append(max(O2_vs_design_parameter))
            O2_min.append(min(O2_vs_design_parameter))
        O2_ax.legend()
        O2_ax.grid()
        O2_ax.set_ylabel('O2 [1]', **font)
        O2_ax.set_xlabel('Count of design variants', **font)

        # for ecce digest
        fig_ecce = figure(figsize=(10, 5), facecolor='w', edgecolor='k')
        O2_ecce_ax = fig_ecce.gca()
        O2_ecce_ax.plot(list(range(-1, 22)), O2_ref*np.ones(23), 'k--', label='Reference design')
        O2_ecce_ax.plot(O2_ecce_data[1], 'o-', lw=0.75, alpha=0.5,      label=r'$L_g$')
        O2_ecce_ax.plot(O2_ecce_data[0], 'v-', lw=0.75, alpha=0.5,      label=r'$w_{st}$')
        O2_ecce_ax.plot(O2_ecce_data[3], 's-', lw=0.75, alpha=0.5,      label=r'$w_{rt}$')
        O2_ecce_ax.plot(O2_ecce_data[5], '^-', lw=0.75, alpha=0.5,      label=r'$\theta_{so}$')
        O2_ecce_ax.plot(O2_ecce_data[2], 'd-', lw=0.75, alpha=0.5,      label=r'$w_{ro}$')
        O2_ecce_ax.plot(O2_ecce_data[6], '*-', lw=0.75, alpha=0.5,      label=r'$d_{so}$')
        O2_ecce_ax.plot(O2_ecce_data[4], 'X-', lw=0.75, alpha=0.5,      label=r'$d_{ro}$')

        myfontsize = 12.5
        from pylab import plt
        plt.rcParams.update({'font.size': myfontsize})


        # Reference candidate design
        ref = np.zeros(8)
            # ref[0] = 0.635489                                   # PF
            # ref[1] = 0.963698                                   # eta
            # ref[1] = efficiency_at_50kW(1817.22+216.216+224.706)# eta@50kW

        if self.reference_design is not None:
            list_plotting_weights = [8, 3, self.required_torque, 0.1, self.rotor_weight, 0.2, 10, 2500]
            ref[0] = O2_ref                  / list_plotting_weights[0] 
            ref[1] = O1_ref                  / list_plotting_weights[1] 
            ref[2] = self.reference_data[2]  / list_plotting_weights[2]  # 100%
            ref[3] = self.reference_data[3]  / list_plotting_weights[3]  # 100%
            ref[4] = self.reference_data[4]  / list_plotting_weights[4]  # 100% = FRW
            ref[5] = self.reference_data[5]  / list_plotting_weights[5]  # 100%
            ref[6] = self.reference_data[6]  / list_plotting_weights[6]  # deg
            ref[7] = self.reference_data[INDEX_TOTAL_LOSS] / list_plotting_weights[7]  # W

        O1_ax.plot(list(range(-1, 22)), O1_ref*np.ones(23), 'k--')
        O2_ax.plot(list(range(-1, 22)), O2_ref*np.ones(23), 'k--')
        O2_ecce_ax.legend()
        O2_ecce_ax.grid()
        O2_ecce_ax.set_xticks(list(range(21)))
        O2_ecce_ax.annotate('Lower bound', xytext=(0.5, 5.5), xy=(0, 4), xycoords='data', arrowprops=dict(arrowstyle="->"))
        O2_ecce_ax.annotate('Upper bound', xytext=(18.0, 5.5),  xy=(20, 4), xycoords='data', arrowprops=dict(arrowstyle="->"))
        O2_ecce_ax.set_xlim((-0.5,20.5))
        O2_ecce_ax.set_ylim((0,14)) # 4,14
        O2_ecce_ax.set_xlabel(r'Number of design variant', **font)
        O2_ecce_ax.set_ylabel(r'$O_2(x)$ [1]', **font)
        fig_ecce.tight_layout()
        # fig_ecce.savefig(r'D:\OneDrive\[00]GetWorking\32 blimopti\p2019_ecce_bearingless_induction_full_paper\images\O2_vs_params.png', dpi=150)
        # plt.show()
        # quit() ###################################


        # Maximum
        data_max = np.array(data_max)
        O1_max   = np.array(O1_max)
        O2_max   = np.array(O2_max)
            # data_max[0] = (data_max[0])                   # PF
            # data_max[1] = (data_max[1])                   # eta
            # data_max[1] = efficiency_at_50kW(data_max[7]) # eta@50kW # should use data_min[7] because less loss, higher efficiency
        data_max[0] = O2_max       / list_plotting_weights[0]  
        data_max[1] = O1_max       / list_plotting_weights[1]  
        data_max[2] = (data_max[2])/ list_plotting_weights[2]  # 100%
        data_max[3] = (data_max[3])/ list_plotting_weights[3]  # 100%
        data_max[4] = (data_max[4])/ list_plotting_weights[4]  # 100% = FRW
        data_max[5] = (data_max[5])/ list_plotting_weights[5]  # 100%
        data_max[6] = (data_max[6])/ list_plotting_weights[6]  # deg
        data_max[7] = (data_max[7])/ list_plotting_weights[7]  # W
        y_max_vs_design_parameter_0 = [el[0] for el in data_max]
        y_max_vs_design_parameter_1 = [el[1] for el in data_max]
        y_max_vs_design_parameter_2 = [el[2] for el in data_max]
        y_max_vs_design_parameter_3 = [el[3] for el in data_max]
        y_max_vs_design_parameter_4 = [el[4] for el in data_max]
        y_max_vs_design_parameter_5 = [el[5] for el in data_max]
        y_max_vs_design_parameter_6 = [el[6] for el in data_max]

        # Minimum
        data_min = np.array(data_min)
        O1_min   = np.array(O1_min)
        O2_min   = np.array(O2_min)
            # data_min[0] = (data_min[0])                    # PF
            # data_min[1] = (data_min[1])                    # eta
            # data_min[1] = efficiency_at_50kW(data_min[7])  # eta@50kW
        data_min[0] = O2_min        / list_plotting_weights[0] 
        data_min[1] = O1_min        / list_plotting_weights[1] 
        data_min[2] = (data_min[2]) / list_plotting_weights[2] # 100%
        data_min[3] = (data_min[3]) / list_plotting_weights[3] # 100%
        data_min[4] = (data_min[4]) / list_plotting_weights[4] # 100% = FRW
        data_min[5] = (data_min[5]) / list_plotting_weights[5] # 100%
        data_min[6] = (data_min[6]) / list_plotting_weights[6] # deg
        data_min[7] = (data_min[7]) / list_plotting_weights[7] # W
        y_min_vs_design_parameter_0 = [el[0] for el in data_min]
        y_min_vs_design_parameter_1 = [el[1] for el in data_min]
        y_min_vs_design_parameter_2 = [el[2] for el in data_min]
        y_min_vs_design_parameter_3 = [el[3] for el in data_min]
        y_min_vs_design_parameter_4 = [el[4] for el in data_min]
        y_min_vs_design_parameter_5 = [el[5] for el in data_min]
        y_min_vs_design_parameter_6 = [el[6] for el in data_min]

        count = np.arange(len(y_max_vs_design_parameter_0))  # the x locations for the groups
        width = 1.0  # the width of the bar

        fig = figure(dpi=150, figsize=(16, 8), facecolor='w', edgecolor='k')
        ax = fig.gca()
        # fig, ax = plt.subplots(dpi=150, figsize=(16, 8), facecolor='w', edgecolor='k')                                      #  #1034A
        rects1 = ax.bar(count - 3*width/8, y_min_vs_design_parameter_0, width/8, alpha=0.5, label=r'$L_g$, Air gap length', color='#6593F5')
        rects2 = ax.bar(count - 2*width/8, y_min_vs_design_parameter_1, width/8, alpha=0.5, label=r'$b_{st}$, Stator tooth width', color='#1D2951') # https://digitalsynopsis.com/design/beautiful-color-palettes-combinations-schemes/
        rects3 = ax.bar(count - 1*width/8, y_min_vs_design_parameter_2, width/8, alpha=0.5, label=r'$b_{rt}$, Rotor tooth width', color='#03396c')
        rects4 = ax.bar(count - 0*width/8, y_min_vs_design_parameter_3, width/8, alpha=0.5, label=r'$\theta_{so}$, Stator open width', color='#6497b1')
        rects5 = ax.bar(count + 1*width/8, y_min_vs_design_parameter_4, width/8, alpha=0.5, label=r'$w_{ro}$, Rotor open width',  color='#0E4D92')
        rects6 = ax.bar(count + 2*width/8, y_min_vs_design_parameter_5, width/8, alpha=0.5, label=r'$d_{so}$, Stator open depth', color='#005b96')
        rects7 = ax.bar(count + 3*width/8, y_min_vs_design_parameter_6, width/8, alpha=0.5, label=r'$d_{ro}$, Rotor open depth', color='#b3cde0') 
        print('ylim=', ax.get_ylim())
        autolabel(ax, rects1, bias=-0.10, textfont=textfont)
        autolabel(ax, rects2, bias=-0.10, textfont=textfont)
        autolabel(ax, rects3, bias=-0.10, textfont=textfont)
        autolabel(ax, rects4, bias=-0.10, textfont=textfont)
        autolabel(ax, rects5, bias=-0.10, textfont=textfont)
        autolabel(ax, rects6, bias=-0.10, textfont=textfont)
        autolabel(ax, rects7, bias=-0.10, textfont=textfont)
        one_one = np.array([1, 1])
        minus_one_one = np.array([-1, 1])
        ax.plot(rects4[0].get_x() + 0.5*width*minus_one_one, ref[0]*one_one, 'k--', lw=1.0, alpha=0.6, label='Reference design' )
        ax.plot(rects4[1].get_x() + 0.5*width*minus_one_one, ref[1]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[2].get_x() + 0.5*width*minus_one_one, ref[2]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[3].get_x() + 0.5*width*minus_one_one, ref[3]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[4].get_x() + 0.5*width*minus_one_one, ref[4]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[5].get_x() + 0.5*width*minus_one_one, ref[5]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[6].get_x() + 0.5*width*minus_one_one, ref[6]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.plot(rects4[7].get_x() + 0.5*width*minus_one_one, ref[7]*one_one, 'k--', lw=1.0, alpha=0.6 )
        ax.legend(loc='upper right', prop={'family':'Times New Roman'})
        # text for indicating reference values
        ax.text(rects4[0].get_x() - 3.5/8*width, ref[0]*1.01, '%.2f'%(ref[0]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[1].get_x() - 3.5/8*width, ref[1]*1.01, '%.2f'%(ref[1]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[2].get_x() - 3.5/8*width, ref[2]*1.01, '%.2f'%(ref[2]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[3].get_x() - 3.5/8*width, ref[3]*1.01, '%.2f'%(ref[3]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[4].get_x() - 3.5/8*width, ref[4]*1.01, '%.2f'%(ref[4]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[5].get_x() - 3.5/8*width, ref[5]*1.01, '%.2f'%(ref[5]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[6].get_x() - 3.5/8*width, ref[6]*1.01, '%.2f'%(ref[6]), ha='center', va='bottom', rotation=90, **textfont)
        ax.text(rects4[7].get_x() - 3.5/8*width, ref[7]*1.01, '%.2f'%(ref[7]), ha='center', va='bottom', rotation=90, **textfont)

        rects1 = ax.bar(count - 3*width/8, y_max_vs_design_parameter_0, width/8, alpha=0.5, label=r'$L_g$,         Air gap length', color='#6593F5')    # bottom=y_min_vs_design_parameter_0, 
        rects2 = ax.bar(count - 2*width/8, y_max_vs_design_parameter_1, width/8, alpha=0.5, label=r'$b_{st}$, Stator tooth width', color='#1D2951')     # bottom=y_min_vs_design_parameter_1, 
        rects3 = ax.bar(count - 1*width/8, y_max_vs_design_parameter_2, width/8, alpha=0.5, label=r'$b_{rt}$, Rotor tooth width', color='#03396c')      # bottom=y_min_vs_design_parameter_2, 
        rects4 = ax.bar(count - 0*width/8, y_max_vs_design_parameter_3, width/8, alpha=0.5, label=r'$\theta_{so}$, Stator open width', color='#6497b1') # bottom=y_min_vs_design_parameter_3, 
        rects5 = ax.bar(count + 1*width/8, y_max_vs_design_parameter_4, width/8, alpha=0.5, label=r'$w_{ro}$, Rotor open width',  color='#0E4D92')      # bottom=y_min_vs_design_parameter_4, 
        rects6 = ax.bar(count + 2*width/8, y_max_vs_design_parameter_5, width/8, alpha=0.5, label=r'$d_{so}$, Stator open depth', color='#005b96')      # bottom=y_min_vs_design_parameter_5, 
        rects7 = ax.bar(count + 3*width/8, y_max_vs_design_parameter_6, width/8, alpha=0.5, label=r'$d_{ro}$, Rotor open depth', color='#b3cde0')       # bottom=y_min_vs_design_parameter_6, 
        autolabel(ax, rects1, textfont=textfont)
        autolabel(ax, rects2, textfont=textfont)
        autolabel(ax, rects3, textfont=textfont)
        autolabel(ax, rects4, textfont=textfont)
        autolabel(ax, rects5, textfont=textfont)
        autolabel(ax, rects6, textfont=textfont)
        autolabel(ax, rects7, textfont=textfont)

        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylabel('Normalized Objective Functions', **font)
        ax.set_xticks(count)
        # ax.set_xticklabels(('Power Factor [100%]', r'$\eta$@$T_{em}$ [100%]', r'$T_{em}$ [15.9 N]', r'$T_{rip}$ [10%]', r'$|F|$ [51.2 N]', r'    $E_m$ [20%]', r'      $E_a$ [10 deg]', r'$P_{\rm Cu,Fe}$ [2.5 kW]')))
        # ax.set_xticklabels(('Power Factor [100%]', r'$O_1$ [3]', r'$T_{em}$ [15.9 N]', r'$T_{rip}$ [10%]', r'$|F|$ [51.2 N]', r'    $E_m$ [20%]', r'      $E_a$ [10 deg]', r'$P_{\rm Cu,Fe}$ [2.5 kW]'))
        ax.set_xticklabels(('$O_2$ [%g]'               %(list_plotting_weights[0]), 
                            '$O_1$ [%g]'               %(list_plotting_weights[1]), 
                            '$T_{em}$ [%g Nm]'         %(list_plotting_weights[2]), 
                            '$T_{rip}$ [%g%%]'         %(list_plotting_weights[3]*100), 
                            '$|F|$ [%g N]'             %(list_plotting_weights[4]), 
                            '    $E_m$ [%g%%]'         %(list_plotting_weights[5]*100), 
                            '      $E_a$ [%g deg]'     %(list_plotting_weights[6]), 
                            '$P_{\\rm Cu,Fe}$ [%g kW]' %(list_plotting_weights[7]*1e-3) ), **font)
        ax.grid()
        ax.set_ylim([0,4])
        # fig.tight_layout()
        # fig.savefig(r'D:\OneDrive\[00]GetWorking\32 blimopti\p2019_ecce_bearingless_induction\images\sensitivity_results.png', dpi=150)

        # plt.show()
        return results_for_refining_bounds


    def get_x_denorm_from_design_parameters(self, design_parameters, bound_filter=None):
        if bound_filter is None:
            x_denorm = design_parameters
            return x_denorm

        if len(bound_filter) == 13:
            # step 1: get free_variables from design_parameters
            free_variables = [None]*13
            free_variables[0]  = design_parameters[0] # spmsm_template.deg_alpha_st 
            free_variables[1]  = design_parameters[3] # spmsm_template.mm_d_so         
            free_variables[2]  = design_parameters[5] # spmsm_template.mm_d_st
            free_variables[3]  = sum([design_parameters[i] for i in (2,4,5,6)]) # spmsm_template.mm_r_si + spmsm_template.mm_d_sp + spmsm_template.mm_d_st + spmsm_template.mm_d_sy # stator outer radius
            free_variables[4]  = design_parameters[7] # spmsm_template.mm_w_st         
            free_variables[5]  = design_parameters[12] # spmsm_template.sleeve_length   
            free_variables[6]  = design_parameters[14] # spmsm_template.mm_d_pm         
            free_variables[7]  = design_parameters[15] # spmsm_template.deg_alpha_rm    
            free_variables[8]  = design_parameters[16] # spmsm_template.deg_alpha_rs    
            free_variables[9]  = design_parameters[17] # spmsm_template.mm_d_ri         
            free_variables[10] = sum([design_parameters[i] for i in (18,17,19)]) # spmsm_template.mm_r_ri + spmsm_template.mm_d_ri + spmsm_template.mm_d_rp
            free_variables[11] = design_parameters[19] # spmsm_template.mm_d_rp         
            free_variables[12] = design_parameters[20] # spmsm_template.mm_d_rs         
        elif len(bound_filter) == 9:
            free_variables = design_parameters # For IM, free_variables are design_parameters (even always having the same length)

        # step 2: get x_denorm from free_variables
        x_denorm = []
        for idx, boo in enumerate(bound_filter):
            if boo == 1:
                x_denorm.append( free_variables[idx] )
        return x_denorm

class FEA_Solver:
    def __init__(self, fea_config_dict):
        self.fea_config_dict = fea_config_dict
        self.app = None

        self.output_dir = self.fea_config_dict['dir_parent'] + self.fea_config_dict['run_folder']
        self.dir_csv_output_folder = self.output_dir + 'csv/'
        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)
        if not os.path.isdir(self.dir_csv_output_folder):
            os.makedirs(self.dir_csv_output_folder)

        # post-process feature
        self.fig_main, self.axeses = plt.subplots(2, 2, sharex=True, dpi=150, figsize=(16, 8), facecolor='w', edgecolor='k')
        utility.pyplot_clear(self.axeses)

        self.folder_to_be_deleted = None

        # if os.path.exists(self.output_dir+'swarm_MOO_log.txt'):
        #     os.rename(self.output_dir+'swarm_MOO_log.txt', self.output_dir+'swarm_MOO_log_backup.txt')
        open(self.output_dir+'swarm_MOO_log.txt', 'a').close()

    def read_swarm_survivor(self, popsize):
        if not os.path.exists(self.output_dir + 'swarm_survivor.txt'):
            return None

        with open(self.output_dir + 'swarm_survivor.txt', 'r') as f:
            buf = f.readlines()
            survivor_data_raw = buf[-popsize:]
            survivor_data = [[float(s) for s in line.split(',')] for line in survivor_data_raw]

        self.survivor_title_number = float(buf[-popsize-1][len('---------'):])

        # for el in survivor_data:
        #     print('\t', el)
        # quit()
        return survivor_data

    def write_swarm_survivor(self, pop, counter_fitness_return):
        with open(self.output_dir + 'swarm_survivor.txt', 'a') as f:
            f.write('\n---------%d\n'%(counter_fitness_return) \
                    + '\n'.join(','.join('%.16f'%(x) for x in el[0].tolist() + el[1].tolist() ) for el in zip(pop.get_x(), pop.get_f()) )) # convert 2d array to string

    def read_swarm_data(self, bound_filter=None):
        if not os.path.exists(self.output_dir + 'swarm_data.txt'):
            print('\tNo file!')
            return None

        print('\tRead in', self.output_dir + 'swarm_data.txt')
        with open(self.output_dir + 'swarm_data.txt', 'r') as f:
            buf = f.readlines()
            buf = buf[1:]
            length_buf = len(buf) 

            if length_buf % 21 == 0:
                pass
            else:
                raise Exception('Invalid swarm_data.txt!')

            number_of_chromosome = length_buf / 21
            if number_of_chromosome == 0:
                return None

            self.swarm_data_raw = [buf[i:i+21] for i in range(0, len(buf), 21)]
            self.swarm_data_container = swarm_data_container(self.swarm_data_raw, self.fea_config_dict, bound_filter)
            self.swarm_data = self.swarm_data_container.swarm_data_xf
            # for el in self.swarm_data:
            #     print(el)
            # quit()
            return int(number_of_chromosome)

            # while True:
            #     try:
            #         if 'Extra Info:' in buf.pop():
            #             info_is_at_this_line = buf[-10]
            #             loc_first_comma = info_is_at_this_line.find(',') + 1
            #             loc_second_comma = info_is_at_this_line.find(',', loc_first_comma+1)
            #             counter = int(info_is_at_this_line[loc_first_comma+1, loc_second_comma])
            #             print(counter)
            #             quit()
            #             break
            #     except:
            #         print('swarm_data.txt is empty')
            #         return None

    def fea_bearingless_spmsm(self, spmsm_template, x_denorm, counter, counter_loop, bool_re_evaluate=False):
        logger = logging.getLogger(__name__)
        msg = 'SPMSM: Run FEA for individual #%d'%(counter)
        logger.info(msg)
        print(msg)

        # get local design variant
        print('#', counter, 'x_denorm', x_denorm)
        spmsm_variant = bearingless_spmsm_design.bearingless_spmsm_design(
                                            spmsm_template=spmsm_template,
                                            x_denorm=x_denorm,
                                            counter=counter,
                                            counter_loop=counter_loop
                                            )

        spmsm_variant.spec = spmsm_template.spec
        self.spmsm_variant = spmsm_variant

        # print('------------------Js = %g'%(spmsm_variant.Js))

        # project name
        if counter_loop == 1:
            self.project_name          = 'proj%d'%(counter)
        else:
            self.project_name          = 'proj%d-redo%d'%(counter, counter_loop)
        self.expected_project_file = self.output_dir + "%s.jproj"%(self.project_name)

        # study name
        study_name = spmsm_variant.name + "TranPMSM" # Change here and there 
        # ind1Tran_torque

        # Leave the solving task to JMAG
        if bool_re_evaluate==False:

            # def draw_jmag_bpmsm():
            import JMAG
            toolJd = JMAG.JMAG(self.fea_config_dict)

            toolJd.open(self.expected_project_file)
            DRAW_SUCCESS = spmsm_variant.draw_spmsm(toolJd)
            if DRAW_SUCCESS != 1:
                raise Exception('Drawer failed.')

            app = toolJd.app

            # JMAG
            if app.NumModels()>=1:
                model = app.GetModel(spmsm_variant.name)
            else:
                logger.error('there is no model yet for %s'%(spmsm_variant.name))
                raise Exception('why is there no model yet? %s'%(spmsm_variant.name))

            spmsm_variant.pre_process(app, model)
            study = spmsm_variant.add_magnetic_transient_study(app, model, self.dir_csv_output_folder, study_name) # Change here and there 
            self.mesh_study(spmsm_variant, app, model, study)
            # raise KeyboardInterrupt
            self.run_study(spmsm_variant, app, study, clock_time())


            # export Voltage if field data exists.
            if self.fea_config_dict['delete_results_after_calculation'] == False:
                # Export Circuit Voltage
                ref1 = app.GetDataManager().GetDataSet("Circuit Voltage")
                app.GetDataManager().CreateGraphModel(ref1)
                app.GetDataManager().GetGraphModel("Circuit Voltage").WriteTable(self.dir_csv_output_folder + spmsm_variant.name + "_EXPORT_CIRCUIT_VOLTAGE.csv")
        else:

            THE_mm2_magnet_area = self.spmsm_variant.rotorMagnet.draw(None, bool_re_evaluate=True)
            THE_mm2_slot_area = self.spmsm_variant.coils.draw(None, bool_re_evaluate=True)

            CurrentAmp_in_the_slot = THE_mm2_slot_area * spmsm_variant.fill_factor * spmsm_variant.Js*1e-6 * np.sqrt(2) #/2.2*2.8
            CurrentAmp_per_conductor = CurrentAmp_in_the_slot / spmsm_variant.DriveW_zQ
            CurrentAmp_per_phase = CurrentAmp_per_conductor * spmsm_variant.wily.number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
            variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor
            spmsm_variant.CurrentAmp_per_phase = CurrentAmp_per_phase
            spmsm_variant.DriveW_CurrentAmp = spmsm_variant.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
            spmsm_variant.BeariW_CurrentAmp = spmsm_variant.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp

            slot_area_utilizing_ratio = (spmsm_variant.DriveW_CurrentAmp + spmsm_variant.BeariW_CurrentAmp) / spmsm_variant.CurrentAmp_per_phase
            print('---Heads up! slot_area_utilizing_ratio is', slot_area_utilizing_ratio)

            print('---Variant CurrentAmp_in_the_slot =', CurrentAmp_in_the_slot)
            print('---variant_DriveW_CurrentAmp = CurrentAmp_per_phase =', variant_DriveW_CurrentAmp)
            print('---spmsm_variant.DriveW_CurrentAmp =', spmsm_variant.DriveW_CurrentAmp)
            print('---spmsm_variant.BeariW_CurrentAmp =', spmsm_variant.BeariW_CurrentAmp)
            print('---TORQUE_CURRENT_RATIO:', spmsm_variant.fea_config_dict['TORQUE_CURRENT_RATIO'])
            print('---SUSPENSION_CURRENT_RATIO:', spmsm_variant.fea_config_dict['SUSPENSION_CURRENT_RATIO'])


        ################################################################
        # Load data for cost function evaluation
        ################################################################
        results_to_be_unpacked = utility.build_str_results(self.axeses, spmsm_variant, self.project_name, study_name, self.dir_csv_output_folder, self.fea_config_dict, femm_solver=None, machine_type='PMSM')
        if results_to_be_unpacked is not None:
            self.fig_main.savefig(self.output_dir + spmsm_variant.name + 'results.png', dpi=150)
            utility.pyplot_clear(self.axeses)
            # show()
            return results_to_be_unpacked # spmsm_variant.stack_length, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle, jmag_loss_list, femm_loss_list, power_factor, total_loss
        else:
            raise Exception('results_to_be_unpacked is None.')

    def fea_bearingless_induction(self, im_template, x_denorm, counter, counter_loop, bool_re_evaluate=False):
        logger = logging.getLogger(__name__)
        print('Run FEA for individual #%d'%(counter))

        # get local design variant
        im_variant = population.bearingless_induction_motor_design.local_design_variant(im_template, 0, counter, x_denorm)
        if counter_loop == 1:
            im_variant.name = 'ind%d'%(counter)
        else:
            im_variant.name = 'ind%d-redo%d'%(counter, counter_loop)
        im_variant.spec = im_template.spec
        self.im_variant = im_variant
        self.femm_solver = FEMM_Solver.FEMM_Solver(self.im_variant, flag_read_from_jmag=False, freq=50) # eddy+static
        im = None

        if counter_loop == 1:
            self.project_name          = 'proj%d'%(counter)
        else:
            self.project_name          = 'proj%d-redo%d'%(counter, counter_loop)
        self.expected_project_file = self.output_dir + "%s.jproj"%(self.project_name)

        original_study_name = im_variant.name + "Freq"
        tran2tss_study_name = im_variant.name + 'Tran2TSS'

        self.dir_femm_temp         = self.output_dir + 'femm_temp/'
        self.femm_output_file_path = self.dir_femm_temp + original_study_name + '.csv'

        # self.jmag_control_state = False

        # local scripts
        def open_jmag(expected_project_file_path):
            if self.app is None:
                app = win32com.client.Dispatch('designer.Application.171')
                if self.fea_config_dict['designer.Show'] == True:
                    app.Show()
                else:
                    app.Hide()
                # app.Quit()
                self.app = app # means that the JMAG Designer is turned ON now.
                self.bool_run_in_JMAG_Script_Editor = False

                def add_steel(self):
                    print('[First run on this computer detected]', self.fea_config_dict['Steel'], 'is added to jmag material library.')
                    import population
                    if 'M15' in self.fea_config_dict['Steel']:
                        population.add_M1xSteel(self.app, self.fea_config_dict['dir_parent'], steel_name="M-15 Steel")
                    elif 'M19' in self.fea_config_dict['Steel']:
                        population.add_M1xSteel(self.app, self.fea_config_dict['dir_parent'])
                    elif 'Arnon5' == self.fea_config_dict['Steel']:
                        population.add_Arnon5(self.app, self.fea_config_dict['dir_parent'])        

                # too avoid tons of the same material in JAMG's material library
                fname = self.fea_config_dict['dir_parent'] + '.jmag_state.txt'
                if not os.path.exists(fname):
                    with open(fname, 'w') as f:
                        f.write(self.fea_config_dict['pc_name'] + '/' + self.fea_config_dict['Steel'] + '\n')
                    add_steel(self)
                else:
                    with open(fname, 'r') as f:
                        for line in f.readlines():
                            if self.fea_config_dict['pc_name'] + '/' + self.fea_config_dict['Steel'] not in line:
                                add_steel(self)

            else:
                app = self.app

            print(expected_project_file_path)
            if os.path.exists(expected_project_file_path):
                print('JMAG project exists already. I learned my lessions. I will NOT delete it but create a new one with a different name instead.')
                # os.remove(expected_project_file_path)
                attempts = 2
                temp_path = expected_project_file_path[:-len('.jproj')] + 'attempts%d.jproj'%(attempts)
                while os.path.exists(temp_path):
                    attempts += 1
                    temp_path = expected_project_file_path[:-len('.jproj')] + 'attempts%d.jproj'%(attempts)

                expected_project_file_path = temp_path

            app.NewProject("Untitled")
            app.SaveAs(expected_project_file_path)
            logger.debug('Create JMAG project file: %s'%(expected_project_file_path))

            return app

        def draw_jmag(app):
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # Draw the model in JMAG Designer
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            DRAW_SUCCESS = self.draw_jmag_induction(app,
                                                counter, 
                                                im_variant,
                                                im_variant.name)
            if DRAW_SUCCESS == 0:
                # TODO: skip this model and its evaluation
                cost_function = 99999 # penalty
                logging.getLogger(__name__).warn('Draw Failed for %s-%s\nCost function penalty = %g.%s', self.project_name, im_variant.name, cost_function, self.im_variant.show(toString=True))
                raise Exception('Draw Failed: Are you working on the PC? Sometime you by mistake operate in the JMAG Geometry Editor, then it fails to draw.')
                return None
            elif DRAW_SUCCESS == -1:
                raise Exception(' DRAW_SUCCESS == -1:')

            # JMAG
            if app.NumModels()>=1:
                model = app.GetModel(im_variant.name)
            else:
                logger.error('there is no model yet for %s'%(im_variant.name))
                raise Exception('why is there no model yet? %s'%(im_variant.name))
            return model

        if bool_re_evaluate==False:
            # this should be summoned even before initializing femm, and it will decide whether the femm results are reliable
            app = open_jmag(self.expected_project_file) # will set self.jmag_control_state to True

            ################################################################
            # Begin from where left: Frequency Study
            ################################################################
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # Eddy Current Solver for Breakdown Torque and Slip
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            if self.fea_config_dict['jmag_run_list'][0] == 0:
                # check for existing results
                if os.path.exists(self.femm_output_file_path):
                    # for file in os.listdir(self.dir_femm_temp):
                    #     if original_study_name in file:
                    #         print('----------', original_study_name, file)
                    print('Remove legacy femm output files @ %s'%(self.femm_output_file_path))
                    os.remove(self.femm_output_file_path)
                    os.remove(self.femm_output_file_path[:-4]+'.fem')
                    # quit()
                    # quit()


                # At this point, no results exist from femm.
                print('Run greedy_search_for_breakdown_slip...')
                femm_tic = clock_time()
                # self.femm_solver.__init__(im_variant, flag_read_from_jmag=False, freq=50.0)
                if im_variant.DriveW_poles == 2:
                    self.femm_solver.greedy_search_for_breakdown_slip( self.dir_femm_temp, original_study_name, 
                                                                        bool_run_in_JMAG_Script_Editor=self.bool_run_in_JMAG_Script_Editor, fraction=1) # 转子导条必须形成通路
                else:
                    self.femm_solver.greedy_search_for_breakdown_slip( self.dir_femm_temp, original_study_name, 
                                                                        bool_run_in_JMAG_Script_Editor=self.bool_run_in_JMAG_Script_Editor, fraction=2)
            else:
                raise Exception("if self.fea_config_dict['jmag_run_list'][0] == 0:")

            ################################################################
            # Begin from where left: Transient Study
            ################################################################
            model = draw_jmag(app)

            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # TranFEAwi2TSS for ripples and iron loss
            #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # add or duplicate study for transient FEA denpending on jmag_run_list
            if self.fea_config_dict['jmag_run_list'][0] == 0:
                # FEMM+JMAG
                study = im_variant.add_TranFEAwi2TSS_study( 50.0, app, model, self.dir_csv_output_folder, tran2tss_study_name, logger)
                self.mesh_study(im_variant, app, model, study)

                # wait for femm to finish, and get your slip of breakdown
                slip_freq_breakdown_torque, breakdown_torque, breakdown_force = self.femm_solver.wait_greedy_search(femm_tic)

                # Now we have the slip, set it up!
                im_variant.update_mechanical_parameters(slip_freq_breakdown_torque) # do this for records only
                if im_variant.the_slip != slip_freq_breakdown_torque / im_variant.DriveW_Freq:
                    raise Exception('Check update_mechanical_parameters().')
                study.GetDesignTable().GetEquation("slip").SetExpression("%g"%(im_variant.the_slip))

                self.run_study(im_variant, app, study, clock_time())
            else:
                raise Exception("self.fea_config_dict['jmag_run_list'][0] == 0: TranFEAwi2TSS")

            # export Voltage if field data exists.
            if self.fea_config_dict['delete_results_after_calculation'] == False:
                # Export Circuit Voltage
                ref1 = app.GetDataManager().GetDataSet("Circuit Voltage")
                app.GetDataManager().CreateGraphModel(ref1)
                app.GetDataManager().GetGraphModel("Circuit Voltage").WriteTable(self.dir_csv_output_folder + im_variant.name + "_EXPORT_CIRCUIT_VOLTAGE.csv")
        else:
            # FEMM的转子电流，哎，是个麻烦事儿。。。

            self.femm_solver.vals_results_rotor_current = []

            new_fname = self.dir_femm_temp + original_study_name + '.csv'
            try:
                with open(new_fname, 'r') as f:
                    buf = f.readlines()
                    for idx, el in enumerate(buf):
                        if idx == 0:
                            slip_freq_breakdown_torque = float(el)
                            continue
                        if idx == 1:
                            breakdown_torque = float(el)
                            continue
                        if idx == 2:
                            self.femm_solver.stator_slot_area = float(el)
                            continue
                        if idx == 3:
                            self.femm_solver.rotor_slot_area = float(el)
                            continue
                        # print(el)
                        # print(el.split(','))
                        temp = el.split(',')
                        self.femm_solver.vals_results_rotor_current.append( float(temp[0])+ 1j*float(temp[1]) )
                        # print(self.femm_solver.vals_results_rotor_current)
                self.dirty_backup_stator_slot_area = self.femm_solver.stator_slot_area
                self.dirty_backup_rotor_slot_area = self.femm_solver.rotor_slot_area
                self.dirty_backup_vals_results_rotor_current = self.femm_solver.vals_results_rotor_current
            except FileNotFoundError as error:
                print(error)
                print('Use dirty_backup to continue...') # 有些时候，不知道为什么femm的结果文件（.csv）没了，这时候曲线救国，凑活一下吧
                self.femm_solver.stator_slot_area = self.dirty_backup_stator_slot_area           
                self.femm_solver.rotor_slot_area = self.dirty_backup_rotor_slot_area            
                self.femm_solver.vals_results_rotor_current = self.dirty_backup_vals_results_rotor_current 


            # 电机的电流值取决于槽的面积。。。。
            THE_mm2_slot_area = self.femm_solver.stator_slot_area*1e6

            if 'VariableStatorSlotDepth' in self.fea_config_dict['which_filter']:
                # set DriveW_CurrentAmp using the calculated stator slot area.
                print('[A]: DriveW_CurrentAmp is updated.')

                # 槽深变化，电密不变，所以电流也会变化。
                CurrentAmp_in_the_slot = THE_mm2_slot_area * im_variant.fill_factor * im_variant.Js*1e-6 * np.sqrt(2)
                CurrentAmp_per_conductor = CurrentAmp_in_the_slot / im_variant.DriveW_zQ
                CurrentAmp_per_phase = CurrentAmp_per_conductor * im_variant.wily.number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
                variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor

                im_variant.CurrentAmp_per_phase = CurrentAmp_per_phase

                im_variant.DriveW_CurrentAmp = self.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
                im_variant.BeariW_CurrentAmp = self.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp

                slot_area_utilizing_ratio = (im_variant.DriveW_CurrentAmp + im_variant.BeariW_CurrentAmp) / im_variant.CurrentAmp_per_phase
                print('---Heads up! slot_area_utilizing_ratio is', slot_area_utilizing_ratio)
                
                print('---Variant CurrentAmp_in_the_slot =', CurrentAmp_in_the_slot)
                print('---variant_DriveW_CurrentAmp = CurrentAmp_per_phase =', variant_DriveW_CurrentAmp)
                print('---im_variant.DriveW_CurrentAmp =', im_variant.DriveW_CurrentAmp)
                print('---im_variant.BeariW_CurrentAmp =', im_variant.BeariW_CurrentAmp)
                print('---TORQUE_CURRENT_RATIO:', self.fea_config_dict['TORQUE_CURRENT_RATIO'])
                print('---SUSPENSION_CURRENT_RATIO:', self.fea_config_dict['SUSPENSION_CURRENT_RATIO'])



        ################################################################
        # Load data for cost function evaluation
        ################################################################
        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        # Load Results for Tran2TSS
        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        results_to_be_unpacked = utility.build_str_results(self.axeses, im_variant, self.project_name, tran2tss_study_name, self.dir_csv_output_folder, self.fea_config_dict, self.femm_solver, machine_type='IM')
        if results_to_be_unpacked is not None:
            self.fig_main.savefig(self.output_dir + im_variant.name + 'results.png', dpi=150)
            utility.pyplot_clear(self.axeses)
            # show()
            return results_to_be_unpacked # im_variant.stack_length, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle, jmag_loss_list, femm_loss_list, power_factor, total_loss
        else:
            raise Exception('results_to_be_unpacked is None.')
        # winding analysis? 之前的python代码利用起来啊
        # 希望的效果是：设定好一个设计，马上进行运行求解，把我要看的数据都以latex报告的形式呈现出来。
        # OP_PS_Qr36_M19Gauge29_DPNV_NoEndRing.jproj

    def draw_jmag_induction(self, app, individual_index, im_variant, model_name, bool_trimDrawer_or_vanGogh=True, doNotRotateCopy=False):

        if individual_index == -1: # 后处理是-1
            print('Draw model for post-processing')
            if individual_index+1 + 1 <= app.NumModels():
                logger = logging.getLogger(__name__)
                logger.debug('The model already exists for individual with index=%d. Skip it.', individual_index)
                return -1 # the model is already drawn

        elif individual_index+1 <= app.NumModels(): # 一般是从零起步
            logger = logging.getLogger(__name__)
            logger.debug('The model already exists for individual with index=%d. Skip it.', individual_index)
            return -1 # the model is already drawn

        # open JMAG Geometry Editor
        app.LaunchGeometryEditor()
        geomApp = app.CreateGeometryEditor()
        # geomApp.Show()
        geomApp.NewDocument()
        doc = geomApp.GetDocument()
        ass = doc.GetAssembly()

        # draw parts
        try:
            if bool_trimDrawer_or_vanGogh:
                d = population.TrimDrawer(im_variant) # 传递的是地址哦
                d.doc, d.ass = doc, ass
                d.plot_shaft("Shaft")

                d.plot_rotorCore("Rotor Core")
                d.plot_cage("Cage")

                d.plot_statorCore("Stator Core")

                d.plot_coil("Coil")
                # d.plot_airWithinRotorSlots(u"Air Within Rotor Slots")

                if 'VariableStatorSlotDepth' in self.fea_config_dict['which_filter']:
                    # set DriveW_CurrentAmp using the calculated stator slot area.
                    print('[A]: DriveW_CurrentAmp is updated.')

                    # 槽深变化，电密不变，所以电流也会变化。
                    CurrentAmp_in_the_slot = d.mm2_slot_area * im_variant.fill_factor * im_variant.Js*1e-6 * np.sqrt(2)
                    CurrentAmp_per_conductor = CurrentAmp_in_the_slot / im_variant.DriveW_zQ
                    CurrentAmp_per_phase = CurrentAmp_per_conductor * im_variant.wily.number_parallel_branch # 跟几层绕组根本没关系！除以zQ的时候，就已经变成每根导体的电流了。
                    variant_DriveW_CurrentAmp = CurrentAmp_per_phase # this current amp value is for non-bearingless motor

                    im_variant.CurrentAmp_per_phase = CurrentAmp_per_phase

                    im_variant.DriveW_CurrentAmp = self.fea_config_dict['TORQUE_CURRENT_RATIO'] * variant_DriveW_CurrentAmp 
                    im_variant.BeariW_CurrentAmp = self.fea_config_dict['SUSPENSION_CURRENT_RATIO'] * variant_DriveW_CurrentAmp

                    slot_area_utilizing_ratio = (im_variant.DriveW_CurrentAmp + im_variant.BeariW_CurrentAmp) / im_variant.CurrentAmp_per_phase
                    print('---Heads up! slot_area_utilizing_ratio is', slot_area_utilizing_ratio)
                    
                    print('---Variant CurrentAmp_in_the_slot =', CurrentAmp_in_the_slot)
                    print('---variant_DriveW_CurrentAmp = CurrentAmp_per_phase =', variant_DriveW_CurrentAmp)
                    print('---im_variant.DriveW_CurrentAmp =', im_variant.DriveW_CurrentAmp)
                    print('---im_variant.BeariW_CurrentAmp =', im_variant.BeariW_CurrentAmp)
                    print('---TORQUE_CURRENT_RATIO:', self.fea_config_dict['TORQUE_CURRENT_RATIO'])
                    print('---SUSPENSION_CURRENT_RATIO:', self.fea_config_dict['SUSPENSION_CURRENT_RATIO'])

            else:
                d = VanGogh_JMAG(im_variant, doNotRotateCopy=doNotRotateCopy) # 传递的是地址哦
                d.doc, d.ass = doc, ass
                d.draw_model()
            self.d = d
        except Exception as e:
            print('See log file to plotting error.')
            logger = logging.getLogger(__name__)
            logger.error('The drawing is terminated. Please check whether the specified bounds are proper.', exc_info=True)

            raise e

            # print 'Draw Failed'
            # if self.pc_name == 'Y730':
            #     # and send the email to hory chen
            #     raise e

            # or you can skip this model and continue the optimization!
            return False # indicating the model cannot be drawn with the script.

        # Import Model into Designer
        doc.SaveModel(True) # True=on : Project is also saved. 
        model = app.GetCurrentModel() # model = app.GetModel(u"IM_DEMO_1")
        model.SetName(model_name)
        model.SetDescription(im_variant.model_name_prefix + '\n' + im_variant.show(toString=True))

        if doNotRotateCopy:
            im_variant.pre_process_structural(app, d.listKeyPoints)
        else:
            im_variant.pre_process(app)

        model.CloseCadLink() # this is essential if you want to create a series of models
        return True

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
        app.ExportImageWithSize(self.output_dir + model.GetName() + '.png', 2000, 2000)
        app.View().ShowModel() # 1st btn. close mesh view, and note that mesh data will be deleted if only ouput table results are selected.

class acm_designer(object):
    def __init__(self, fea_config_dict, spec):

        self.spec = spec
        self.solver = FEA_Solver(fea_config_dict)
        self.fea_config_dict = fea_config_dict

        self.flag_do_not_evaluate_when_init_pop = False

    def init_logger(self, prefix='pygmo_'):
        self.logger = utility.myLogger(self.fea_config_dict['dir_codes'], prefix=prefix+self.fea_config_dict['run_folder'][:-1])

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Automatic Report Generation
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def build_oneReport(self):
        if '730' in self.fea_config_dict['pc_name']:
            os.system('cd /d "'+ self.fea_config_dict['dir_parent'] + 'release/OneReport/OneReport_TEX" && z_nul"')

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Talk to Database
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def talk_to_mysql_database(self):
        if self.spec.bool_bad_specifications:
            print('\nThe specifiaction can not be fulfilled. Read script log or OneReport.pdf for information and revise the specifiaction for $J_r$ or else your design name is wrong.')
        else:
            print('\nThe specifiaction is meet. Now check the database of blimuw if on Y730.')
            if '730' in self.fea_config_dict['pc_name']:
                utility.communicate_database(self.spec)

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Automatic Performance Evaluation (This is just a wraper)
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def evaluate_design(self, acm_template, x_denorm, counter, counter_loop=1):
        # print(dir(acm_template))
        if 'SM' in acm_template.name:
            function = self.solver.fea_bearingless_spmsm
        else:
            function = self.solver.fea_bearingless_induction

        return function(acm_template, x_denorm, counter, counter_loop, bool_re_evaluate=self.bool_re_evaluate)

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 1. Bounds for DE optimiazation
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # def get_original_bounds(self, Jr_max=8e6):

        # from math import tan, pi
        # 定子齿宽最小值 = 1
        # 定子齿宽最大值 = tan(2*pi/self.spec.Qs*0.5)*self.spec.Radius_OuterRotor * 2 # 圆（半径为Radius_OuterRotor）的外接正多边形（Regular polygon）的边长
        # # print(定子齿宽最大值, 2*pi*self.spec.Radius_OuterRotor/self.spec.Qs) # 跟弧长比应该比较接近说明对了

        # 转子齿宽最小值 = 1
        # 内接圆的半径 = self.spec.Radius_OuterRotor - (self.spec.d_ro + self.spec.Radius_of_RotorSlot)
        # 转子齿宽最大值 = tan(2*pi/self.spec.Qr*0.5)*内接圆的半径 * 2 
        # # print(转子齿宽最大值, 2*pi*内接圆的半径/self.spec.Qr) # 跟弧长比应该比较接近说明对了

        # self.original_bounds = [ [           1.0,                3],          # air_gap_length_delta
        #                          [定子齿宽最小值,   定子齿宽最大值],#--# stator_tooth_width_b_ds
        #                          [转子齿宽最小值,   转子齿宽最大值],#--# rotor_tooth_width_b_dr
        #                          [             1, 360/self.spec.Qs],           # Angle_StatorSlotOpen
        #                          [          5e-1,                3],           # Width_RotorSlotOpen 
        #                          [          5e-1,                3],           # Width_StatorTeethHeadThickness
        #                          [          5e-1,                3] ]          # Length_HeadNeckRotorSlot
        # # 定子齿范围检查
        # # 定子齿再宽，都可以无限加长轭部来满足导电面积。
        # # stator_inner_radius_r_is_eff = stator_inner_radius_r_is + (width_statorTeethHeadThickness + width_StatorTeethNeck)
        # # temp = (2*pi*stator_inner_radius_r_is_eff - self.Qs*stator_tooth_width_b_ds)
        # # stator_tooth_height_h_ds = ( sqrt(temp**2 + 4*pi*area_stator_slot_Sus*self.Qs) - temp ) / (2*pi)

        # def check_valid_rotor_slot_height(rotor_tooth_width_b_dr, Jr_max):

        #     area_conductor_rotor_Scr = self.spec.rotor_current_actual / Jr_max
        #     area_rotor_slot_Sur = area_conductor_rotor_Scr
            
        #     rotor_outer_radius_r_or_eff = 1e-3*(self.spec.Radius_OuterRotor - self.spec.d_ro)

        #     slot_height, _, _ = pyrhonen_procedure_as_function.get_parallel_tooth_height(area_rotor_slot_Sur, rotor_tooth_width_b_dr, self.spec.Qr, rotor_outer_radius_r_or_eff)
        #     return np.isnan(slot_height)


        # # 转子齿范围检查
        # 下界, 上界 = self.original_bounds[2][0], self.original_bounds[2][1]
        # 步长 = (上界-下界)*0.05
        # list_valid_tooth_width = []
        # for rotor_tooth_width_b_dr in np.arange(下界, 上界, 步长):
        #     # print('b_dr =', rotor_tooth_width_b_dr)
        #     list_valid_tooth_width.append( check_valid_rotor_slot_height(rotor_tooth_width_b_dr, Jr_max) ) # 8e6 from Pyrhonen's book for copper
        # # print(list_valid_tooth_width)
        # 有效上界 = 下界
        # for ind, el in enumerate(list_valid_tooth_width):
        #     if el == True:
        #         break
        #     else:
        #         有效上界 += 步长
        # self.original_bounds[2][1] = 有效上界

        # return self.original_bounds

    # def get_classic_bounds(self):
        # if 'SM' in self.spec.acm_template.name:
        #     Q = self.spec.pmsm_template.Q
        #     s = self.spec.pmsm_template.s
        #     p = self.spec.pmsm_template.p
        #     PMSM = self.spec.pmsm_template
        #     self.classic_bounds =  [ 
        #                         [ 0.35*360/Q, 0.9*360/Q],    # deg_alpha_st        = free_variables[0]
        #                         [  0.5,   5],                # mm_d_so             = free_variables[1]
        #                         [0.8*PMSM.mm_d_st,                1.2*PMSM.mm_d_st], # mm_d_st    = free_variables[2]
        #                         [0.8*PMSM.Radius_OuterStatorYoke, 1.2*PMSM.Radius_OuterStatorYoke], # stator_outer_radius = free_variables[3]
        #                         [0.8*PMSM.mm_w_st,                1.2*PMSM.mm_w_st], # mm_w_st    = free_variables[4]
        #                         [3,   4],                    # sleeve_length       = free_variables[5]
        #                         [2.5, 7],                    # mm_d_pm             = free_variables[6]
        #                         [0.6*360/(2*p), 1.0*360/(2*p)],      # deg_alpha_rm        = free_variables[7]
        #                         [0.8*360/(2*p)/s, 0.975*360/(2*p)/s], # deg_alpha_rs        = free_variables[8]
        #                         [0.8*PMSM.mm_d_ri,  1.2*PMSM.mm_d_ri], # mm_d_ri   = free_variables[9]
        #                         [0.8*PMSM.Radius_OuterRotor, 1.2*PMSM.Radius_OuterRotor], # rotor_outer_radius  = free_variables[10] 
        #                         [2.5,   6],                  # mm_d_rp             = free_variables[11]
        #                         [2.5,   6] ]                 # mm_d_rs             = free_variables[12]
        #     self.original_bounds = self.classic_bounds
        #     return self.classic_bounds
        # else:
        #     self.get_original_bounds()
        #     self.classic_bounds = [ [self.spec.delta*0.9, self.spec.delta*2  ],          # air_gap_length_delta
        #                             [self.spec.w_st *0.5, self.spec.w_st *1.5],          #--# stator_tooth_width_b_ds
        #                             [self.spec.w_rt *0.5, self.spec.w_rt *1.5],          #--# rotor_tooth_width_b_dr
        #                             [                1.5,                  12],           # Angle_StatorSlotOpen
        #                             [               5e-1,                   3],           # Width_RotorSlotOpen 
        #                             [               5e-1,                   3],           # Width_StatorTeethHeadThickness
        #                             [               5e-1,                   3] ]          # Length_HeadNeckRotorSlot
        #     # classic_bounds cannot be beyond original_bounds
        #     index = 0
        #     for A, B in zip(self.classic_bounds, self.original_bounds):
        #         if A[0] < B[0]:
        #             self.classic_bounds[index] = B[0]
        #         if A[1] > B[1]:
        #             self.classic_bounds[index] = B[1]
        #         index += 1
                
        #     return self.classic_bounds

    def get_de_config(self):

        self.de_config_dict = { 'original_bounds': self.get_original_bounds(),
                                'mut':        0.8,
                                'crossp':     0.7,
                                'popsize':    35, # 5~10 \times number of geometry parameters --JAC223
                                'iterations': 70,
                                'narrow_bounds_normalized':[[],
                                                            [],
                                                            [],
                                                            [],
                                                            [],
                                                            [],
                                                            [] ], # != []*7 （完全是两回事）
                                'bounds':None}
        return self.de_config_dict

    def run_local_sensitivity_analysis(self, the_bounds, design_denorm=None):
        # if design_denorm not in the_bounds: then raise 
        if design_denorm is not None:
            for ind, el in enumerate(design_denorm):
                if el < the_bounds[ind][0] or el > the_bounds[ind][1]:
                    raise Exception('给的设计不在边界的内部')

        self.logger.debug('---------\nBegin Local Sensitivity Analysis')

        # de_config_dict['bounds'] 还没有被赋值
        self.de_config_dict['bounds'] = the_bounds

        self.init_swarm() # define app.sw
        self.sw.generate_pop(specified_initial_design_denorm=design_denorm)
        de_generator = self.sw.de()
        for result in de_generator:
            print(result)

    def check_results_of_local_sensitivity_analysis(self):
        if self.fea_config_dict['local_sensitivity_analysis'] == True:
            run_folder = self.fea_config_dict['run_folder'][:-1] + 'lsa/'
        else:
            run_folder = self.fea_config_dict['run_folder']

        # Sensitivity Bar Charts
        return os.path.exists(fea_config_dict['dir_parent'] + 'pop/' + run_folder + 'swarm_data.txt')

    def collect_results_of_local_sensitivity_analysis(self):
        print('Start to collect results for local sensitivity analysis...')
        try:
            results_for_refining_bounds = utility.build_sensitivity_bar_charts(self.spec, self.sw)
            # quit()
        except Exception as e:
            raise e
            os.remove(fea_config_dict['dir_parent'] + 'pop/' + run_folder + 'swarm_data.txt')
            print('Remove ' + fea_config_dict['dir_parent'] + 'pop/' + run_folder + 'swarm_data.txt')
            print('Continue for sensitivity analysis...')

        self.results_for_refining_bounds = results_for_refining_bounds
        # return results_for_refining_bounds

    def build_refined_bounds(self, the_bounds):

        de_config_dict = self.de_config_dict
        number_of_variants= self.fea_config_dict['local_sensitivity_analysis_number_of_variants'] # 故意不加1的

        print('-'*20, 'results_for_refining_bounds (a.k.a. results_for_refining_bounds):')
        str_to_file = 'When you are done, replace this line with "Done", save file and close all figures.\n'
        for key, val in self.results_for_refining_bounds.items():
            print('---', key)
            for el in val:
                print('\t', el)
            str_to_file += key + '\n' + '\n'.join( ','.join(f'{x}' for x in y) for y in val ) + '\n'

        self.logger.debug('The default refining factors will be for: %s'%(self.fea_config_dict['use_weights']))

        with open('./refining_factors.txt', 'w') as f:
            f.write(str_to_file)
        file_backup_user_input = './refining_factors_%s.txt'%(self.fea_config_dict['run_folder'][:-1])
        if os.path.exists(file_backup_user_input):
            os.system('start %s'%(file_backup_user_input))
        os.system('start ./refining_factors.txt')
        from pylab import show
        show()
        from time import sleep
        while True:
            with open('./refining_factors.txt') as f:
                buf = f.read()
                if buf[:4] == 'Done' or buf[:4] == 'done':
                    print('Done')
                    break
                else:
                    print('.', end='')
                    sleep(1)
        if os.path.exists(file_backup_user_input):
            os.remove(file_backup_user_input)
        os.rename('./refining_factors.txt', file_backup_user_input)

        buf_list = buf.split('\n')
        # for el in buf_list:
        #     print(el)

        print('fea_config_dict - use_weights changed from %s to %s'%(self.fea_config_dict['use_weights'], buf_list[1]))
        self.fea_config_dict['use_weights'] = buf_list[1]

        user_input_for_refining_bounds = [[float(x) for x in y.split(',')] for y in buf_list[2:2+7]]
        print('user_input_for_refining_bounds')
        for el in user_input_for_refining_bounds:
            print('\t', el)

        for ind, bound in enumerate(user_input_for_refining_bounds):
            下界 = bound[0]
            if 下界 != 0:
                下界 -= 1
            上界 = bound[-1]
            if 上界 != number_of_variants:
                上界 += 1
            self.de_config_dict['narrow_bounds_normalized'][ind].append(下界/number_of_variants)
            self.de_config_dict['narrow_bounds_normalized'][ind].append(上界/number_of_variants)

        self.de_config_dict['bounds'] = []
        for bnd1, bnd2 in zip(self.de_config_dict['original_bounds'], self.de_config_dict['narrow_bounds_normalized']):
            diff = bnd1[1] - bnd1[0]
            self.de_config_dict['bounds'].append( [ bnd1[0]+diff*bnd2[0] , bnd1[0]+diff*bnd2[1] ]) # 注意，都是乘以original_bounds的上限哦！

        print('-'*40)
        print('narrow_bounds_normalized:')
        for el in self.de_config_dict['narrow_bounds_normalized']:
            print('\t', el)

        print('original_bounds:')
        for el in self.de_config_dict['original_bounds']:
            print('\t', el)

        print('refined bounds:')
        for el in self.de_config_dict['bounds']:
            print('\t', el)

        return de_config_dict['bounds']

    def build_local_bounds_from_best_design(self, best_design):
        raise Exception('build_local_bounds_from_best_design')
        return local_bounds

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 2. Initilize Swarm and Initial Pyrhonen's Design (Run this part in JMAG) and femm solver (if required by run_list)
    #    Bounds: 1e-1也还是太小了（第三次报错），至少0.5mm长吧 
    #    # 1e-1 is the least geometry value. 
    #    A 1e-2 will leads to：转子闭口槽极限，会导致edge过小，从而报错：small arc entity exists.png
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def init_swarm(self):
        self.sw = population.swarm(self.fea_config_dict, de_config_dict=self.de_config_dict)
        # sw.show(which='all')
        # print sw.im.show(toString=True)
        # quit()

        if self.fea_config_dict['jmag_run_list'][0] == 0:
            self.sw.init_femm_solver() # define app.sw.femm_solver        

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 3. Run DE Optimization
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def run_de(self):
        sw = self.sw
        logger = self.logger

        count_abort = 0
        logger.debug('-------------------------count_abort=%d' % (count_abort))

        # if optimization_flat == True:
        # generate the initial generation
        sw.generate_pop(specified_initial_design_denorm=None)
        # quit()

        # add initial_design of Pyrhonen09 to the initial generation
        if sw.fea_config_dict['local_sensitivity_analysis'] == False:
            if count_abort == 0:
                utility.add_Pyrhonen_design_to_first_generation(sw, self.de_config_dict, logger)

        # write FEA config to disk
        sw.write_to_file_fea_config_dict()

        if True:
            try:
                de_generator = sw.de()
                # run
                # result = list(de_generator)
                for result in de_generator:
                    print(result)
            except Exception as e:
                print('See log file for the error msg.')
                logger.error('Optimization aborted.', exc_info=True)

                raise e
                # 避免死循环
                count_abort+1
                if count_abort > 10:
                    quit()
            else:
                logger.info('Done.')
                utility.send_notification('Done.')

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 4. Post-processing
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def best_design_by_weights(self, use_weights):
        if use_weights != self.fea_config_dict['use_weights']:
            raise Exception('Not implemented')

        # import utility
        swda = utility.build_Pareto_plot(self.spec, self.sw)
        # sw.fobj(999, individual_denorm)

        # Final LaTeX Report
        print('According to Pyrhonen09, 300 MPa is the typical yield stress of iron core.')
        initial_design = utility.Pyrhonen_design(self.sw.im, self.de_config_dict['bounds'])
        print (initial_design.design_parameters_denorm)
        print (swda.best_design_denorm)
        print (self.de_config_dict['bounds'])

        best_report_dir_prefix = '../release/OneReport/BestReport_TEX/contents/'
        file_name = 'final_report'
        file_suffix = '.tex'
        fname = open(best_report_dir_prefix+file_name+'_s01'+file_suffix, 'w', encoding='utf-8')
        def combine_lists_alternating_2(list1, list2):
            if abs(len(list1) - len(list2))<=1:
                result = [None]*(len(list1)+len(list2))
                result[::2] = list1
                result[1::2] = list2
                return result
            else:
                raise Exception('Try this (not tested).') # https://stackoverflow.com/questions/3678869/pythonic-way-to-combine-two-lists-in-an-alternating-fashion
                import itertools 
                return [x for x in itertools.chain.from_iterable(itertools.izip_longest(list1,list2)) if x]
        def combine_lists_alternating_4(list1, list2, list3, list4):
            result = [None]*(len(list1)+len(list2)+len(list3)+len(list4))
            result[::4] = list1
            result[1::4] = list2
            result[2::4] = list3
            result[3::4] = list4
            return result
        lower_bounds = [el[0] for el in self.de_config_dict['bounds']]
        upper_bounds = [el[1] for el in self.de_config_dict['bounds']]
        design_data = combine_lists_alternating_4(  initial_design.design_parameters_denorm, 
                                                    swda.best_design_denorm,
                                                    lower_bounds,
                                                    upper_bounds,)
        latex_table = r'''
            \begin{table*}[!t]
              \caption{Comparison of the key geometry parameters between best design and initial design}
              \centering
                \begin{tabular}{ccccc}
                    \hline
                    \hline
                    \thead{Geometry parameters} &
                    \thead{Initial\\\relax design} &
                    \thead{The best\\\relax design} &
                    \thead{Lower\\\relax bounds} &
                    \thead{Upper\\\relax bounds} \\
                    \hline
                    Air gap length $L_g$                    & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\
                    Stator tooth width $w_{st}$             & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\
                    Rotor tooth width $w_{rt}$              & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\
                    Stator slot open angle $\theta_{so}$    & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\
                    Rotor slot open width $w_{ro}$          & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\
                    Stator slot open depth $d_{so}$         & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\
                    Rotor slot open depth $d_{ro}$          & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\
                    \hline
                    \vspace{-2.5ex}
                    \\
                    \multicolumn{3}{l}{*Note: Blah.}
                \end{tabular}
              \label{tab:001}
              \vspace{-3ex}
            \end{table*}
            ''' % tuple(design_data)
        print('''\\subsection{Best Design of Objective Function %s}\n\n
                    ''' % (self.fea_config_dict['use_weights']) + latex_table, file=fname)
        print(swda.str_best_design_details, file=fname)
        fname.close()

        os.system('cd /d '+ r'"D:\OneDrive - UW-Madison\c\release\OneReport\BestReport_TEX" && z_nul"') # 必须先关闭文件！否则编译不起来的
        # import subprocess
        # subprocess.call(r"D:\OneDrive - UW-Madison\c\release\OneReport\OneReport_TEX\z_nul", shell=True)

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 5. Check mechanical strength for the best design
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    def run_static_structural_fea(self, design_denorm):
        sw = self.sw
        im_best = population.bearingless_induction_motor_design.local_design_variant(sw.im, \
                     -1, -1, design_denorm)

        # initialize JMAG Designer
        sw.designer_init()
        sw.app.Show()
        project_name = fea_config_dict['run_folder'][:-1]+'_Best' # spec.build_name() is too long...
        expected_project_file = sw.dir_project_files + "%s.jproj"%(project_name)
        print(expected_project_file)
        if not os.path.exists(expected_project_file):
            sw.app.NewProject("Untitled")
            sw.app.SaveAs(expected_project_file)
            logger.debug('Create JMAG project file: %s'%(expected_project_file))
        else:
            sw.app.Load(expected_project_file)
            logger.debug('Load JMAG project file: %s'%(expected_project_file))
            logger.debug('Existing models of %d are found in %s', sw.app.NumModels(), sw.app.GetDefaultModelFolderPath())

        # draw the model in JMAG Designer
        DRAW_SUCCESS = sw.draw_jmag_induction( -1, 
                                             im_best,
                                             'Best %s %s'%(fea_config_dict['use_weights'], fea_config_dict['run_folder'][:-1]),
                                             bool_trimDrawer_or_vanGogh=False,
                                             doNotRotateCopy=True)
        if DRAW_SUCCESS == 0:
            raise Exception('Drawing failed')
        elif DRAW_SUCCESS == -1:
            print('Model Already Exists')
        # quit()

        model = sw.app.GetCurrentModel()
        if model.NumStudies() == 0:
            expected_csv_output_dir = sw.dir_csv_output_folder+'structural/'
            if not os.path.isdir(expected_csv_output_dir):
                os.makedirs(expected_csv_output_dir)
            study = im_best.add_structural_study(sw.app, model, expected_csv_output_dir) # 文件夹名应该与jproj同名
        else:
            study = model.GetStudy(0)


        if study.AnyCaseHasResult():
            pass
        else:
            # Add cases 
            study.GetDesignTable().AddParameterVariableName(u"Centrifugal_Force (CentrifugalForce2D): AngularVelocity")
            study.GetDesignTable().AddCase()
            study.GetDesignTable().SetValue(1, 0, 45000) # r/min
            study.GetDesignTable().AddCase()
            study.GetDesignTable().SetValue(2, 0, 15000) # r/min
            study.GetDesignTable().AddParameterVariableName(u"MeshSizeControl (ElementSizeOnPart): Size")
            study.GetDesignTable().AddCase()
            study.GetDesignTable().SetValue(3, 1, 0.5) # mm
            study.GetDesignTable().AddCase()
            study.GetDesignTable().SetValue(4, 1, 0.05) # mm

            # run (mesh is included in add_study_structural)
            study.RunAllCases()
            sw.app.Save()

            # Results-Graphs-Calculations-Add Part Calculation
            study.CreateCalculationDefinition(u"VonMisesStress")
            study.GetCalculationDefinition(u"VonMisesStress").SetResultType(u"MisesStress", u"")
            study.GetCalculationDefinition(u"VonMisesStress").SetResultCoordinate(u"Global Rectangular")
            study.GetCalculationDefinition(u"VonMisesStress").SetCalculationType(u"max")
            study.GetCalculationDefinition(u"VonMisesStress").ClearParts()
            study.GetCalculationDefinition(u"VonMisesStress").AddSet(model.GetSetList().GetSet(u"Motion_Region"), 0)

            # show graph, Tab File|Edit|Calculation, click on Calculation - Response Graph Data to register response data
            parameter = sw.app.CreateResponseDataParameter(u"VonMisesStress")
            parameter.SetCalculationType(u"SingleValue")
            parameter.SetStartValue(u"-1")
            parameter.SetEndValue(u"-1")
            parameter.SetUnit(u"s")
            parameter.SetVariable(u"VonMisesStress")
            parameter.SetAllLine(False)
            parameter.SetCaseRangeType(1)
            parameter.SetLine(u"Maximum Value")
            sw.app.GetDataManager().CreateParametricDataWithParameter(study.GetDataSet(u"VonMisesStress", 4), parameter)

            # Contour results
            study.CreateScaling(u"Scale100")
            study.GetScaling(u"Scale100").SetScalingFactor(100)
            # sw.app.View().SetOriginalModelView(True)
            sw.app.View().ShowMeshGeometry()
            sw.app.View().SetScaledDisplacementView(True)
            study.CreateContour(u"MisesElement")
            study.GetContour(u"MisesElement").SetResultType(u"MisesStress", u"")
            study.GetContour(u"MisesElement").SetResultCoordinate(u"Global Rectangular")
            study.GetContour(u"MisesElement").SetContourType(2)
            study.GetContour(u"MisesElement").SetDigitsNotationType(2)

            study.GetContour(u"MisesElement").SetLogScale(True)
            study.GetContour(u"MisesElement").SetNumLabels(u"11")
            study.GetContour(u"MisesElement").SetPrecision(u"1")
            study.GetContour(u"MisesElement").SetGradient(u"PurpleRed", u"11", False)

            sw.app.View().SetContourView(True)



            # study.CreateContour(u"PrincipleStressElement")
            # study.GetContour(u"PrincipleStressElement").SetResultType(u"PrincipalStress", u"")
            # study.GetContour(u"PrincipleStressElement").SetComponent(u"I")
            # study.GetContour(u"PrincipleStressElement").SetContourType(2)
            # study.GetContour(u"PrincipleStressElement").SetDigitsNotationType(2)

            sw.app.Save()


        from pylab import show
        show()

def get_bad_fintess_values(machine_type='IM', ref=False):
    if ref == False:
        if 'IM' in machine_type:
            return 0, 0, 99
        elif 'PMSM' in machine_type:
            return 9999, 0, 999
    else:
        if 'IM' in machine_type:
            return 1, 1, 100
        elif 'PMSM' in machine_type:
            return 10000, 1, 1000        


