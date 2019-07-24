import shutil
from utility import my_execfile
from utility_moo import *
from win32com.client import pywintypes
bool_re_evaluate = False # re-evaluate the designs using csv (without calling FEA softwares)

#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# 0. FEA Setting / General Information & Packages Loading
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
my_execfile('./default_setting.py', g=globals(), l=locals())
fea_config_dict
fea_config_dict['local_sensitivity_analysis'] = False
fea_config_dict['bool_refined_bounds'] = False
fea_config_dict['use_weights'] = 'O2' # this is not working
if 'Y730' in fea_config_dict['pc_name']:
    ################################################################
    # Y730
    ################################################################
    # Combined winding IM
    fea_config_dict['TORQUE_CURRENT_RATIO'] = 0.975
    fea_config_dict['SUSPENSION_CURRENT_RATIO'] = 0.025
    fea_config_dict['which_filter'] = 'FixedStatorSlotDepth' # 'VariableStatorSlotDepth'
    run_folder = r'run#1000/'

fea_config_dict['run_folder'] = run_folder

# spec's
my_execfile('./spec_Prototype2poleOD150mm500Hz_SpecifyTipSpeed.py', g=globals(), l=locals())
spec
fea_config_dict['Active_Qr'] = spec.Qr
fea_config_dict['use_drop_shape_rotor_bar'] = spec.use_drop_shape_rotor_bar
build_model_name_prefix(fea_config_dict) # rebuild the model name for fea_config_dict
spec.build_im_template(fea_config_dict)

# select motor type ehere
print('Build ACM template...')
spec.acm_template = spec.im_template


import acm_designer
global ad
ad = acm_designer.acm_designer(fea_config_dict, spec)
if 'Y730' in fea_config_dict['pc_name']:
    ad.build_oneReport() # require LaTeX
    # ad.talk_to_mysql_database() # require MySQL
    # quit()
ad.init_logger()
ad.bool_re_evaluate = bool_re_evaluate

ad.bounds_denorm = spec.get_im_classic_bounds(which_filter=fea_config_dict['which_filter'])
ad.bound_filter  = spec.bound_filter
print('---------------------\nBounds:')
idx_ad = 0
for idx, f in enumerate(ad.bound_filter):
    if f == True:
        print(idx, f, '[%g,%g]'%tuple(spec.original_template_neighbor_bounds[idx]), '[%g,%g]'%tuple(ad.bounds_denorm[idx_ad]))
        idx_ad += 1
    else:
        print(idx, f, '[%g,%g]'%tuple(spec.original_template_neighbor_bounds[idx]))
print('-'*20)
# quit()
counter_fitness_called, counter_loop = 0, 0

cost_function, f1, f2, f3, FRW, \
normalized_torque_ripple, \
normalized_force_error_magnitude, \
force_error_angle = ad.evaluate_design(ad.spec.im_template, ad.spec.im_template.spec.x_denorm, counter_fitness_called, counter_loop=counter_loop)

print(cost_function, f1, f2, f3, FRW, \
        normalized_torque_ripple, \
        normalized_force_error_magnitude, \
        force_error_angle)
