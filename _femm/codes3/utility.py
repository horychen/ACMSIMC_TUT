# coding:u8
import os
import datetime
import itertools


def my_execfile(filename, g=None, l=None):
    # g=globals(), l=locals()
    exec(compile(open(filename, "rb").read(), filename, 'exec'), g, l)

class ExceptionReTry(Exception):
    """Exception for requiring a second attemp."""
    def __init__(self, message, payload=None):
        self.message = message
        self.payload = payload # you could add more args
    def __str__(self):
        return str(self.message)    

class ExceptionBadNumberOfParts(Exception):
    """Exception for unexpected number of parts in JMAG Designer."""
    def __init__(self, message, payload=None):
        self.message = message
        self.payload = payload # you could add more args
    def __str__(self):
        return str(self.message)

        
def communicate_database(spec):
    try:
        import mysql.connector
    except:
        print('MySQL python connector is not installed. Skip database communication.')
    else:
        db = mysql.connector.connect(
            host ='localhost',
            user ='root',
            passwd ='password123',
            database ='blimuw',
            )
        cursor = db.cursor()
        cursor.execute('SELECT name FROM designs')
        result = cursor.fetchall()
        if spec.build_name() not in [row[0] for row in result]:
            def sql_add_one_record(spec):
                # Add one record
                sql = "INSERT INTO designs " \
                    + "(" \
                        + "name, " \
                            + "PS_or_SC, " \
                            + "DPNV_or_SEPA, " \
                            + "p, " \
                            + "ps, " \
                            + "MecPow, " \
                            + "Freq, " \
                            + "Voltage, " \
                            + "TanStress, " \
                            + "Qs, " \
                            + "Qr, " \
                            + "Js, " \
                            + "Jr, " \
                            + "Coil, " \
                            + "kCu, " \
                            + "Condct, " \
                            + "kAl, " \
                            + "Temp, " \
                            + "Steel, " \
                            + "kFe, " \
                            + "Bds, " \
                            + "Bdr, " \
                            + "Bys, " \
                            + "Byr, " \
                            + "G_b, " \
                            + "G_eta, " \
                            + "G_PF, " \
                            + "debug, " \
                            + "Sskew, " \
                            + "Rskew, " \
                            + "Pitch " \
                    + ") VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
                record = (  spec.build_name(), 
                            'PS' if spec.PS_or_SC else 'SC',
                            'DPNV' if spec.DPNV_or_SEPA else 'SEPA',
                            spec.p,
                            spec.ps,
                            spec.mec_power,
                            spec.ExcitationFreq,
                            spec.VoltageRating,
                            spec.TangentialStress,
                            spec.Qs,
                            spec.Qr,
                            spec.Js,
                            spec.Jr,
                            spec.Coil,
                            spec.space_factor_kCu,
                            spec.Conductor,
                            spec.space_factor_kAl,
                            spec.Temperature,
                            spec.Steel,
                            spec.lamination_stacking_factor_kFe,
                            spec.stator_tooth_flux_density_B_ds,
                            spec.rotor_tooth_flux_density_B_dr,
                            spec.stator_yoke_flux_density_Bys,
                            spec.rotor_yoke_flux_density_Byr,
                            spec.guess_air_gap_flux_density,
                            spec.guess_efficiency,
                            spec.guess_power_factor,
                            spec.debug_or_release,
                            spec.bool_skew_stator,
                            spec.bool_skew_rotor,
                            spec.winding_layout.coil_pitch
                        )
                cursor.execute(sql, record)
                db.commit()
            sql_add_one_record(spec)
            'A new record is added to table named designs.'
        else:
            'Record already exists, skip database communication.'

import operator
def get_index_and_max(the_list):
    return max(enumerate(the_list), key=operator.itemgetter(1))
    # index, max_value

def get_index_and_min(the_list):
    return min(enumerate(the_list), key=operator.itemgetter(1))
    # index, min_value

def gcd(a,b):
    while b:
        a,b = b, a%b
    return a
# Now find LCM using GCD
def lcm(a,b):
    return a*b // gcd(a,b)
# print lcm(8,7.5)

import logging
def myLogger(dir_codes, prefix='default_prefix_'): # This works even when the module is reloaded (which is not the case of the other answers) https://stackoverflow.com/questions/7173033/duplicate-log-output-when-using-python-logging-module

    # logging.getLogger("imported_module").setLevel(logging.WARNING) # disable logging from matplotlib and others
    logging.getLogger('matplotlib').setLevel(logging.WARNING)

    logger=logging.getLogger()
    if not len(logger.handlers):
        logger.setLevel(logging.DEBUG)
        now = datetime.datetime.now()

        if not os.path.isdir(dir_codes + 'log/'):
            os.makedirs(dir_codes + 'log/')

        # create a file handler
        handler=logging.FileHandler(dir_codes + 'log/' + prefix + '-' + now.strftime("%Y-%m-%d") +'.log')
        handler.setLevel(logging.DEBUG)

        # create a logging format
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)

        # add the handlers to the logger
        logger.addHandler(handler)
    return logger

def logger_init(): # This will lead to duplicated logging output
    # logger = logging.getLogger(__name__) # this is used in modules 
    logger = logging.getLogger() # use this (root) in the main executable file
    logger.setLevel(logging.DEBUG)

    # create a file handler
    now = datetime.datetime.now()
    if not os.path.isdir(dir_codes + 'log/'):
        os.makedir(dir_codes + 'log/')
    handler = logging.FileHandler(dir_codes + r'opti_script.log')
    handler.setLevel(logging.DEBUG)

    # create a logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)

    # add the handlers to the logger
    logger.addHandler(handler)
    return logger



import sys
# For a better solution for print, see https://stackoverflow.com/questions/4230855/why-am-i-getting-ioerror-9-bad-file-descriptor-error-while-making-print-st/4230866
# decorater used to block function printing to the console
def blockPrinting(func):
    def func_wrapper(*args, **kwargs):
        # block all printing to the console
        sys.stdout = open(os.devnull, 'w')
        # call the method in question
        func(*args, **kwargs)
        # enable all printing to the console
        sys.stdout = sys.__stdout__

    return func_wrapper
# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')
# Restore
def enablePrint():
    sys.stdout = sys.__stdout__



import math
def to_precision(x,p=4):
    """ http://randlet.com/blog/python-significant-figures-format/
    returns a string representation of x formatted with a precision of p

    Based on the webkit javascript implementation taken from here:
    https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
    """

    x = float(x)

    if x == 0.:
        return "0." + "0"*(p-1)

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(math.log10(x))
    tens = math.pow(10, e - p + 1)
    n = math.floor(x/tens)

    if n < math.pow(10, p - 1):
        e = e -1
        tens = math.pow(10, e - p+1)
        n = math.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens -x):
        n = n + 1

    if n >= math.pow(10,p):
        n = n / 10.
        e = e + 1

    m = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p -1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e+1])
        if e+1 < len(m):
            out.append(".")
            out.extend(m[e+1:])
    else:
        out.append("0.")
        out.extend(["0"]*-(e+1))
        out.append(m)

    return "".join(out)



import numpy as np
def singleSidedDFT(signal, samp_freq):
    NFFT = len(signal)
    dft_complex = np.fft.fft(signal,NFFT) # y is a COMPLEX defined in numpy
    dft_single_sided = [2 * abs(dft_bin) / NFFT for dft_bin in dft_complex][0:int(NFFT/2)+1] # /NFFT for spectrum aplitude consistent with actual signal. 2* for single-sided. abs for amplitude of complem number.
    dft_single_sided[0] *= 0.5 # DC bin in single sided spectrem does not need to be times 2
    return np.array(dft_single_sided)
    dft_freq = 0.5*samp_freq*np.linspace(0,1,NFFT/2+1) # unit is Hz # # f = np.fft.fftfreq(NFFT, Ts) # for double-sided

def basefreqDFT(signal, samp_freq, ax_time_domain=None, ax_freq_domain=None, base_freq=1):
    NFFT = len(signal)

    dft_complex = np.fft.fft(signal,NFFT) # y is a COMPLEX defined in numpy
    dft_double_sided = [abs(dft_bin) / NFFT for dft_bin in dft_complex]
    dft_single_sided = [2 * abs(dft_bin) / NFFT for dft_bin in dft_complex][0:int(NFFT/2)+1] # /NFFT for spectrum aplitude consistent with actual signal. 2* for single-sided. abs for amplitude of complem number.
    dft_single_sided[0] *= 0.5 # DC bin in single sided spectrem does not need to be times 2

    if base_freq==None:
        dft_freq = 0.5*samp_freq*np.linspace(0,1,NFFT/2+1) # unit is Hz # # f = np.fft.fftfreq(NFFT, Ts) # for double-sided
    else: # 0.5*samp_freq is the Nyquist frequency
        dft_freq = 0.5*samp_freq/base_freq*np.linspace(0,1,NFFT/2+1) # unit is per base_freq 

    if ax_time_domain != None:
        Ts = 1.0/samp_freq
        t = [el*Ts for el in range(0,NFFT)]
        ax_time_domain.plot(t, signal, '*', alpha=0.4)
        ax_time_domain.set_xlabel('time [s]')
        ax_time_domain.set_ylabel('B [T]')
    if ax_freq_domain != None:
        ax_freq_domain.plot(dft_freq, dft_single_sided, '*',alpha=0.4)
        ax_freq_domain.set_xlabel('per base frequency')
        ax_time_domain.set_ylabel('B [T]')

class Pyrhonen_design(object):
    def __init__(self, im, bounds=None):
        ''' Determine bounds for these parameters:
            stator_tooth_width_b_ds              = design_parameters[0]*1e-3 # m                       # stator tooth width [mm]
            air_gap_length_delta                 = design_parameters[1]*1e-3 # m                       # air gap length [mm]
            b1                                   = design_parameters[2]*1e-3 # m                       # rotor slot opening [mm]
            rotor_tooth_width_b_dr               = design_parameters[3]*1e-3 # m                       # rotor tooth width [mm]
            self.Length_HeadNeckRotorSlot        = design_parameters[4]      # mm       # rotor tooth head & neck length [mm]
            self.Angle_StatorSlotOpen            = design_parameters[5]      # mm       # stator slot opening [deg]
            self.Width_StatorTeethHeadThickness  = design_parameters[6]      # mm       # stator tooth head length [mm]
        '''
        # rotor_slot_radius = (2*pi*(Radius_OuterRotor - Length_HeadNeckRotorSlot)*1e-3 - rotor_tooth_width_b_dr*Qr) / (2*Qr+2*pi)
        # => rotor_tooth_width_b_dr = ( 2*pi*(Radius_OuterRotor - Length_HeadNeckRotorSlot)*1e-3  - rotor_slot_radius * (2*Qr+2*pi) ) / Qr
        # Verified in Tran2TSS_PS_Opti.xlsx.
            # ( 2*pi*(im.Radius_OuterRotor - im.Length_HeadNeckRotorSlot)  - im.Radius_of_RotorSlot * (2*Qr+2*pi) ) / Qr
            # = ( 2*PI()*(G4 - I4) - J4 * (2*C4+2*PI()) ) / C4
        from math import pi
        Qr = im.Qr

        # unit: mm 
        self.air_gap_length_delta           = im.Length_AirGap
        self.stator_tooth_width_b_ds        = im.Width_StatorTeethBody
        self.rotor_tooth_width_b_dr         = ( 2*pi*(im.Radius_OuterRotor - im.Length_HeadNeckRotorSlot)  - im.Radius_of_RotorSlot * (2*Qr+2*pi) ) / Qr
        self.Angle_StatorSlotOpen           = im.Angle_StatorSlotOpen # deg
        self.b1                             = im.Width_RotorSlotOpen
        self.Width_StatorTeethHeadThickness = im.Width_StatorTeethHeadThickness
        self.Length_HeadNeckRotorSlot       = im.Length_HeadNeckRotorSlot

        self.design_parameters_denorm = [   self.air_gap_length_delta,
                                            self.stator_tooth_width_b_ds,
                                            self.rotor_tooth_width_b_dr,
                                            self.Angle_StatorSlotOpen,
                                            self.b1,
                                            self.Width_StatorTeethHeadThickness,
                                            self.Length_HeadNeckRotorSlot ]

        if bounds is None:
            self.design_parameters_denorm
        else:
            self.show_norm(bounds, self.design_parameters_denorm)


    def show_denorm(self, bounds, design_parameters_norm):
        pop = design_parameters_norm
        min_b, max_b = np.asarray(bounds).T 
        diff = np.fabs(min_b - max_b)
        pop_denorm = min_b + pop * diff
        print('[De-normalized]:', end=' ')
        print(pop_denorm.tolist())
        
    def show_norm(self, bounds, design_parameters_denorm):
        min_b, max_b = np.asarray(bounds).T 
        diff = np.fabs(min_b - max_b)
        print(design_parameters_denorm)
        print(min_b)
        print(bounds)
        self.design_parameters_norm = (design_parameters_denorm - min_b)/diff #= pop
        # print type(self.design_parameters_norm)
        print('[Normalized]:', end=' ')
        print(self.design_parameters_norm.tolist())

        # pop = design_parameters_norm
        # min_b, max_b = np.asarray(bounds).T 
        # diff = np.fabs(min_b - max_b)
        # pop_denorm = min_b + pop * diff
        # print '[De-normalized:]---------------------------Are these two the same?'
        # print pop_denorm.tolist()
        # print design_parameters_denorm

def add_Pyrhonen_design_to_first_generation(sw, de_config_dict, logger):
    initial_design = Pyrhonen_design(sw.im, de_config_dict['bounds'])
    # print 'SWAP!'
    # print initial_design.design_parameters_norm.tolist()
    # print '\nInitial Population:'
    # for index in range(len(sw.init_pop)):
    #     # print sw.init_pop[index].tolist()
    #     print index,
    #     initial_design.show_denorm(de_config_dict['bounds'], sw.init_pop[index])
    # print sw.init_pop[0].tolist()
    sw.init_pop_denorm[0] = initial_design.design_parameters_denorm
    # print sw.init_pop[0].tolist()
    with open(sw.get_gen_file(0), 'r') as f:
        before = f.read()
    with open(sw.get_gen_file(0), 'w') as f:
        f.write('\n'.join(','.join('%.16f'%(x) for x in y) for y in sw.init_pop_denorm)) # convert 2d array to string
    logger.info('Initial design from Pyrhonen09 is added to the first generation of pop (i.e., gen#0000ind#0000).')
    # print('Before:', before)
    # print('After:', '\n'.join(','.join('%.16f'%(x) for x in y) for y in sw.init_pop_denorm))
    # quit()

from smtplib import SMTP
def send_notification(text='Hello', subject='AUTO BLIM OPTIMIZATION NOTIFICATION'):
    content = 'Subject: %s\n\n%s' % (subject, text)
    mail = SMTP('smtp.gmail.com', 587)
    mail.ehlo()
    mail.starttls()
    with open('temp.txt', 'r') as f:
        pswd = f.read()
    mail.login('horyontour@gmail.com', pswd)
    mail.sendmail('horyontour@gmail.com', 'jiahao.chen@wisc.edu', content) 
    mail.close()
    print("Notificaiont sent.")

def get_windage_loss(im_variant, mm_stack_length, TEMPERATURE_OF_AIR=75):

    # %Air friction loss calculation
    nu_0_Air  = 13.3e-6#;  %[m^2/s] kinematic viscosity of air at 0
    rho_0_Air = 1.29#;     %[kg/m^3] Air density at 0
    Shaft = [mm_stack_length,                               #1;         %End position of the sections mm (Absolut)
             im_variant.Radius_OuterRotor+im_variant.Length_AirGap, #1;         %Inner Radius in mm
             1,                                                     #0;         %Shrouded (1) or free surface (0)
             im_variant.Length_AirGap]                              #0];        %Airgap in mm
    Num_shaft_section = 1
    T_Air = TEMPERATURE_OF_AIR #20:(120-20)/((SpeedMax-SpeedMin)/SpeedStep):120         #; % Air temperature []
    
    nu_Air  = nu_0_Air*((T_Air+273)/(0+273))**1.76
    rho_Air = rho_0_Air*(0+273)/(T_Air+273)
    windage_loss_radial = 0 

    # Calculation of the section length ...
    L     = Shaft[0]*1e-3 # in meter
    R     = Shaft[1]*1e-3 # radius of air gap
    delta = Shaft[3]*1e-3 # length of air gap
    
    Omega = 2*np.pi*im_variant.the_speed/60.
    if abs(Omega - im_variant.Omega) < 0.1:
        pass
    else:
        print(Omega, im_variant.Omega, im_variant.the_speed)
        raise Exception('Check speed calc. resutls.')

    # Reynolds number
    Rey = R**2 * (Omega)/nu_Air

    if Shaft[2] == 0: # free running cylinder
        if Rey <= 170:
            c_W = 8. / Rey
        elif Rey>170 and Rey<4000:
            c_W = 0.616*Rey**(-0.5)
        else:
            c_W = 6.3e-2*Rey**(-0.225)
        windage_loss_radial = c_W*np.pi*rho_Air* Omega**3 * R**5 * (1.+L/R)

    else: # shrouded cylinder by air gap from <Loss measurement of a 30 kW High Speed Permanent Magnet Synchronous Machine with Active Magnetic Bearings>
        Tay = R*(Omega)*(delta/nu_Air)*np.sqrt(delta/R) # Taylor number 
        if Rey <= 170:
            c_W = 8. / Rey
        elif Rey>170 and Tay<41.3:
            c_W = 1.8 * Rey**(-1) * delta/R**(-0.25) * (R+delta)**2 / ((R+delta)**2 - R**2)
        else:
            c_W = 7e-3
        windage_loss_radial = c_W*np.pi*rho_Air* Omega**3 * R**4 * L
        
    # end friction loss added - 05192018.yegu
    # the friction coefficients from <Rotor Design of a High-Speed Permanent Magnet Synchronous Machine rating 100,000 rpm at 10 kW>
    Rer = rho_Air * (im_variant.Radius_OuterRotor * 1e-3)**2 * Omega/nu_Air
    if Rer <= 30:
        c_f = 64/3. / Rer
    elif Rer>30 and Rer<3*10**5:
        c_f = 3.87 * Rer**(-0.5)
    else:
        c_f = 0.146 * Rer**(-0.2)

    windage_loss_axial = 0.5 * c_f * rho_Air * Omega**3 * (im_variant.Radius_OuterRotor*1e-3)**5
    
    windage_loss_total = windage_loss_radial + windage_loss_axial
    return windage_loss_total


#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# Evaluation Utility
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
def add_plots(axeses, dm, title=None, label=None, zorder=None, time_list=None, sfv=None, torque=None, range_ss=None, alpha=0.7):

    info = '%s' % (title)
    torque_average = sum(torque[-range_ss:])/len(torque[-range_ss:])
    info += '\nAverage Torque: %g Nm' % (torque_average)
    # torque error = torque - avg. torque
    torque_error = np.array(torque) - torque_average
    ss_max_torque_error = max(torque_error[-range_ss:]), min(torque_error[-range_ss:])
    # we use  half of peak-to-peak value to compute error rather than use peak-to-peak value
    normalized_torque_ripple   = 1.0*(ss_max_torque_error[0] - ss_max_torque_error[1]) / torque_average
    info += '\nNormalized Torque Ripple: %g %%' % (normalized_torque_ripple*100)

    info += '\nAverage Force Mag: %g N'% (sfv.ss_avg_force_magnitude)
    # we use half of peak-to-peak value to compute error rather than use peak-to-peak value
    normalized_force_error_magnitude = sfv.normalized_force_error_magnitude
    info += '\nNormalized Force Error Mag: %g%%, (+)%g%% (-)%g%%' % (normalized_force_error_magnitude*100,
                                                                  sfv.ss_max_force_err_abs[0]/sfv.ss_avg_force_magnitude*100,
                                                                  sfv.ss_max_force_err_abs[1]/sfv.ss_avg_force_magnitude*100)
    # we use peak value to compute error rather than use peak-to-peak value
    # 跟Eric讨论过后，确定了悬浮力的角度误差不能用峰峰值的一半，而是要用最大值和最小值中绝对值更大的那一个。
    force_error_angle = sfv.force_error_angle
    info += '\nMaximum Force Error Angle: %g [deg], (+)%g deg (-)%g deg' % (force_error_angle,
                                                                 sfv.ss_max_force_err_ang[0],
                                                                 sfv.ss_max_force_err_ang[1])
    info += '\nExtra Info:'
    info += '\n\tAverage Force Vector: (%g, %g) N' % (sfv.ss_avg_force_vector[0], sfv.ss_avg_force_vector[1])
    info += '\n\tTorque Ripple (Peak-to-Peak): %g Nm'% ( max(torque[-range_ss:]) - min(torque[-range_ss:]))
    info += '\n\tForce Mag Ripple (Peak-to-Peak): %g N'% (sfv.ss_max_force_err_abs[0] - sfv.ss_max_force_err_abs[1])

    if axeses is not None:
        # plot for torque and force
        ax = axeses[0][0]; ax.plot(time_list, torque,                                           alpha=alpha, label=label, zorder=zorder)
        ax = axeses[0][1]; ax.plot(time_list, sfv.force_abs,                                    alpha=alpha, label=label, zorder=zorder)
        ax = axeses[1][0]; ax.plot(time_list, 100*sfv.force_err_abs/sfv.ss_avg_force_magnitude, label=label, alpha=alpha, zorder=zorder)
        ax = axeses[1][1]; ax.plot(time_list, sfv.force_ang - sfv.ss_avg_force_angle,           label=label, alpha=alpha, zorder=zorder)

    # plot for visialization of power factor 
    # dm.get_voltage_and_current(range_ss)
    # ax = axeses[2][0]; ax.plot(dm.mytime, dm.myvoltage, label=label, alpha=alpha, zorder=zorder)
    # ax = axeses[2][0]; ax.plot(dm.mytime, dm.mycurrent, label=label, alpha=alpha, zorder=zorder)

    return info, torque_average, normalized_torque_ripple, sfv.ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle

def build_str_results(axeses, acm_variant, project_name, tran_study_name, dir_csv_output_folder, fea_config_dict, femm_solver=None, machine_type=None):
    # originate from fobj

    try:
        dm = read_csv_results_4_general_purpose(tran_study_name, dir_csv_output_folder, fea_config_dict, femm_solver, machine_type=machine_type, acm_variant=acm_variant)
    except Exception as e:
        print(e)
        logging.getLogger(__name__).error('Error when loading csv results for Tran2TSS. Check the Report of JMAG Designer. (Maybe Material is not added.)', exc_info=True)
        陌生感 = 'CSV results are not found. Will re-build and re-run the JMAG project...' 
        raise ExceptionReTry(陌生感)
        # return None
        # raise e

    basic_info, time_list, TorCon_list, ForConX_list, ForConY_list, ForConAbs_list = dm.unpack()
    sfv = suspension_force_vector(ForConX_list, ForConY_list, range_ss=fea_config_dict['number_of_steps_2ndTTS']) # samples in the tail that are in steady state
    str_results, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle = \
        add_plots( axeses, dm,
                      title=tran_study_name,
                      label='Transient FEA w/ 2 Time Step Sections',
                      zorder=8,
                      time_list=time_list,
                      sfv=sfv,
                      torque=TorCon_list,
                      range_ss=sfv.range_ss)
    str_results += '\n\tbasic info:' +   ''.join(  [str(el) for el in basic_info])

    if dm.jmag_loss_list is None:
        raise Exception('Loss data is not loaded?')
    else:
        str_results += '\n\tjmag loss info: '  + ', '.join(['%g'%(el) for el in dm.jmag_loss_list]) # dm.jmag_loss_list = [stator_copper_loss, rotor_copper_loss, stator_iron_loss, stator_eddycurrent_loss, stator_hysteresis_loss]

    if fea_config_dict['jmag_run_list'][0] == 0:
        str_results += '\n\tfemm loss info: '  + ', '.join(['%g'%(el) for el in dm.femm_loss_list])

    if fea_config_dict['delete_results_after_calculation'] == False:
        power_factor = dm.power_factor(fea_config_dict['number_of_steps_2ndTTS'], targetFreq=acm_variant.DriveW_Freq)
        str_results += '\n\tPF: %g' % (power_factor)

    # compute the fitness 
    rotor_volume = acm_variant.get_rotor_volume() 
    rotor_weight = acm_variant.get_rotor_weight()
    shaft_power  = acm_variant.Omega * torque_average # make sure update_mechanical_parameters is called so that Omega corresponds to slip_freq_breakdown_torque

    if 'IM' in machine_type:
        if False: # fea_config_dict['jmag_run_list'][0] == 0
            # by JMAG only
            copper_loss  = dm.jmag_loss_list[0] + dm.jmag_loss_list[1] 
            iron_loss    = dm.jmag_loss_list[2] 
        else:
            # by JMAG for iron loss and FEMM for copper loss
            if dm.femm_loss_list[0] is None: # this will happen for running release_design.py
                copper_loss  = dm.jmag_loss_list[0] + dm.jmag_loss_list[1]
            else:
                copper_loss  = dm.femm_loss_list[0] + dm.femm_loss_list[1]
            iron_loss = dm.jmag_loss_list[2] 
    elif 'PMSM' in machine_type:
        # Rotor magnet loss by JMAG
        magnet_Joule_loss = dm.jmag_loss_list[1]
                      # Stator copper loss by Bolognani  
        copper_loss = dm.femm_loss_list[0] + magnet_Joule_loss
        iron_loss = dm.jmag_loss_list[2] 

    windage_loss = get_windage_loss(acm_variant, acm_variant.stack_length)

    # 这样计算效率，输出转矩大的，铁耗大一倍也没关系了，总之就是气隙变得最小。。。要不就不要优化气隙了。。。
    total_loss   = copper_loss + iron_loss + windage_loss
    efficiency   = shaft_power / (total_loss + shaft_power)  # 效率计算：机械功率/(损耗+机械功率)
    str_results  += '\n\teta, windage, total_loss: %g, %g, %g' % (efficiency, windage_loss, total_loss)

    # for easy access to codes
    machine_results = [power_factor, efficiency, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle]
    machine_results.extend(dm.jmag_loss_list)
    if dm.femm_loss_list is None:
        raise
    machine_results.extend(dm.femm_loss_list)
    machine_results.extend([windage_loss, total_loss])

    str_machine_results = ','.join('%g'%(el) for el in machine_results if el is not None) # note that femm_loss_list can be None called by release_design.py
    
    cost_function_O1, list_cost_O1 = compute_list_cost(use_weights('O1'), rotor_volume, rotor_weight, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle, dm.jmag_loss_list, dm.femm_loss_list, power_factor, total_loss)
    cost_function_O2, list_cost_O2 = compute_list_cost(use_weights('O2'), rotor_volume, rotor_weight, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle, dm.jmag_loss_list, dm.femm_loss_list, power_factor, total_loss)

    ################################################################
    # NEW CODES for rated performance
    ################################################################
    # caculate the fitness
    print('-'*40)
    print('Calculate the fitness for', acm_variant.name)

    # LOSS
    if 'IM' in machine_type:
        stator_copper_loss_along_stack = dm.femm_loss_list[2]
        rotor_copper_loss_along_stack  = dm.femm_loss_list[3]

        stator_copper_loss_in_end_turn = dm.femm_loss_list[0] - stator_copper_loss_along_stack 
        rotor_copper_loss_in_end_turn  = dm.femm_loss_list[1] - rotor_copper_loss_along_stack

    elif 'PMSM' in machine_type:
        stator_copper_loss_along_stack = dm.femm_loss_list[2]
        rotor_copper_loss_along_stack  = magnet_Joule_loss

        stator_copper_loss_in_end_turn = dm.femm_loss_list[0] - stator_copper_loss_along_stack 
        rotor_copper_loss_in_end_turn  = 0

    rated_ratio                          = acm_variant.spec.required_torque / torque_average 
    rated_stack_length_mm                = rated_ratio * acm_variant.stack_length
    rated_stator_copper_loss_along_stack = rated_ratio * stator_copper_loss_along_stack
    rated_rotor_copper_loss_along_stack  = rated_ratio * rotor_copper_loss_along_stack
    rated_iron_loss                      = rated_ratio * dm.jmag_loss_list[2]
    rated_windage_loss                   = get_windage_loss(acm_variant, rated_stack_length_mm)

    # total_loss   = copper_loss + iron_loss + windage_loss
    rated_total_loss =  rated_stator_copper_loss_along_stack \
                      + rated_rotor_copper_loss_along_stack \
                      + stator_copper_loss_in_end_turn \
                      + rotor_copper_loss_in_end_turn \
                      + rated_iron_loss \
                      + rated_windage_loss


    # THERMAL
    if 'IM' in machine_type:
        stator_current_density = dm.femm_loss_list[4]
        rotor_current_density  = dm.femm_loss_list[5]
        # print('Current density [Arms/m^2]:', stator_current_density, rotor_current_density, sep='\n')
        # if rotor_current_density > 8e6:
        #     print('rotor_current_density is over 8e6 Arms/m^2')
    else:
                                 # 基波电流幅值（在一根导体里的电流，六相逆变器中的4W相的电流，所以相当于已经考虑了并联支路数了）
        stator_current_density = dm.ui_info[2] / 1.4142135623730951 / (acm_variant.coils.mm2_slot_area*1e-6/acm_variant.DriveW_zQ)
        print('Data Magager: stator_current_density (4W) = %g Arms/m^2'%(stator_current_density))
        rotor_current_density = 0

    print('Required torque: %g Nm'%(acm_variant.spec.required_torque))
    print('acm_variant.Omega: %g rad/s'%(acm_variant.Omega))
    rated_shaft_power  = acm_variant.Omega * acm_variant.spec.required_torque
    rated_efficiency   = rated_shaft_power / (rated_total_loss + rated_shaft_power)  # 效率计算：机械功率/(损耗+机械功率)

    rated_rotor_volume = np.pi*(acm_variant.Radius_OuterRotor*1e-3)**2 * (rated_stack_length_mm*1e-3)
    print('rated_stack_length_mm =', rated_stack_length_mm)

    # This weighted list suggests that peak-to-peak torque ripple of 5% is comparable with Em of 5% or Ea of 1 deg. Ref: Ye gu ECCE 2018
    # Eric suggests Ea is 1 deg. But I think this may be too much emphasis on Ea so large Trip does not matter anymore (not verified yet).
    list_weighted_ripples = [normalized_torque_ripple/0.05, normalized_force_error_magnitude/0.05, force_error_angle]

    if 'IM' in machine_type:
        # - Torque per Rotor Volume
        f1_IM = - acm_variant.spec.required_torque / rated_rotor_volume
        f1 = f1_IM
    elif 'PMSM' in machine_type:
        # - Cost
        price_per_volume_steel    = 0.28  * 61023.744 # $/in^3 (M19 Gauge26) # 0.23 for low carbon, semi-processed 24 Gauge electrical steel
        price_per_volume_copper   = 1.2   * 61023.744 # $/in^3 wire or bar or end-ring
        price_per_volume_magnet   = 11.61 * 61023.744 # $/in^3 NdFeB PM
        # price_per_volume_aluminum = 0.88  / 16387.064 # $/in^3 wire or cast Al
        Vol_Fe = (acm_variant.Radius_OuterStatorYoke*1e-3) ** 2 * (rated_stack_length_mm*1e-3) # 注意，硅钢片切掉的方形部分全部消耗了。
        Vol_PM = (acm_variant.rotorMagnet.mm2_magnet_area*1e-6) * (rated_stack_length_mm*1e-3)
        print('Area_Fe', (acm_variant.Radius_OuterStatorYoke*1e-3) ** 2)
        print('Area_Cu (est.)', dm.Vol_Cu/(rated_stack_length_mm*1e-3))
        print('Area_PM', (acm_variant.rotorMagnet.mm2_magnet_area*1e-6))
        print('Volume_Fe',    Vol_Fe)
        print('Volume_Cu', dm.Vol_Cu)
        print('Volume_PM',    Vol_PM)
        f1_PMSM =    Vol_Fe * price_per_volume_steel \
                + dm.Vol_Cu * price_per_volume_copper\
                +    Vol_PM * price_per_volume_magnet
        print('Cost Fe:',   Vol_Fe * price_per_volume_steel )
        print('Cost Cu:',dm.Vol_Cu * price_per_volume_copper)
        print('Cost PM:',   Vol_PM * price_per_volume_magnet)
        f1 = f1_PMSM

    # - Efficiency @ Rated Power
    f2 = - rated_efficiency
    # Ripple Performance (Weighted Sum)
    f3 = sum(list_weighted_ripples)

    FRW = ss_avg_force_magnitude / rotor_weight
    print('FRW:', FRW, 'Rotor weight:', rotor_weight, 'Stack length:', acm_variant.stack_length, 'Rated stack length:', rated_stack_length_mm)
    rated_rotor_volume = acm_variant.get_rotor_volume(stack_length=rated_stack_length_mm) 
    rated_rotor_weight = acm_variant.get_rotor_weight(stack_length=rated_stack_length_mm)
    print('rated_rotor_volume:', rated_rotor_volume, 'rated_rotor_weight:', rated_rotor_weight)

    rated_results = [   rated_shaft_power, 
                        rated_efficiency,
                        rated_total_loss, 
                        rated_stator_copper_loss_along_stack, 
                        rated_rotor_copper_loss_along_stack, 
                        stator_copper_loss_in_end_turn, 
                        rotor_copper_loss_in_end_turn, 
                        rated_iron_loss, 
                        rated_windage_loss,
                        rated_rotor_volume,
                        rated_stack_length_mm,  # new!
                        acm_variant.stack_length]           # new! 在计算FRW的时候，我们只知道原来的叠长下的力，所以需要知道原来的叠长是多少。

    str_results = '\n-------\n%s-%s\n%d,%d,O1=%g,O2=%g,f1=%g,f2=%g,f3=%g\n%s\n%s\n%s\n' % (
                    project_name, acm_variant.get_individual_name(), 
                    acm_variant.number_current_generation, acm_variant.individual_index, cost_function_O1, cost_function_O2, f1, f2, f3,
                    str_machine_results,
                    ','.join(['%g'%(el) for el in rated_results]), # 改为输出 rated_results
                    ','.join(['%g'%(el) for el in acm_variant.design_parameters]) ) + str_results
 
    # str_results, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle, jmag_loss_list, femm_loss_list, power_factor, total_loss, cost_function = results_to_be_unpacked

    # write design evaluation data to file
    with open(dir_csv_output_folder[:-4] + 'swarm_data.txt', 'a') as f:
        f.write(str_results)

    if fea_config_dict['use_weights'] == 'O1':
        cost_function = cost_function_O1
    elif fea_config_dict['use_weights'] == 'O2':
        cost_function = cost_function_O2
    else:
        raise Exception('Not implemented error.')

    return cost_function, f1, f2, f3, FRW, normalized_torque_ripple, normalized_force_error_magnitude, force_error_angle
 # str_results, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle, dm.jmag_loss_list, dm.femm_loss_list, power_factor, total_loss, cost_function

class suspension_force_vector(object):
    """docstring for suspension_force_vector"""
    def __init__(self, force_x, force_y, range_ss=None): # range_ss means range_steadystate
        super(suspension_force_vector, self).__init__()
        self.force_x = force_x
        self.force_y = force_y
        self.force_ang = []
        temp_force_ang = np.arctan2(force_y, force_x) / np.pi * 180 # [deg]
        for angle in temp_force_ang:
            if angle < 0:
                angle += 360
            self.force_ang.append(angle)
        # print('-'*40+'\nsfv:', self.force_ang)
        self.force_abs = np.sqrt(np.array(force_x)**2 + np.array(force_y)**2 )

        if range_ss == None:
            range_ss = len(force_x)
        self.range_ss = range_ss

        self.ss_avg_force_vector    = np.array([sum(force_x[-range_ss:]), sum(force_y[-range_ss:])]) / range_ss #len(force_x[-range_ss:])
        self.ss_avg_force_angle     = np.arctan2(self.ss_avg_force_vector[1], self.ss_avg_force_vector[0]) / np.pi * 180
        if self.ss_avg_force_angle < 0:
            self.ss_avg_force_angle += 360
        # print('sfv:', self.ss_avg_force_angle)
        self.ss_avg_force_magnitude = np.sqrt(self.ss_avg_force_vector[0]**2 + self.ss_avg_force_vector[1]**2)

        self.force_err_ang = self.force_ang - self.ss_avg_force_angle
        # print('sfv:', self.force_err_ang)
        self.force_err_abs = self.force_abs - self.ss_avg_force_magnitude

        self.ss_max_force_err_ang = max(self.force_err_ang[-range_ss:]), min(self.force_err_ang[-range_ss:])
        self.ss_max_force_err_abs = max(self.force_err_abs[-range_ss:]), min(self.force_err_abs[-range_ss:])

        # method 1 
        # self.force_error_angle = 0.5*(self.ss_max_force_err_ang[0]-self.ss_max_force_err_ang[1])
        # normalized_force_error_magnitude = 0.5*(sfv.ss_max_force_err_abs[0]-sfv.ss_max_force_err_abs[1])/sfv.ss_avg_force_magnitude
        
        # method 2 suggested by Eric 
        self.force_error_angle = max( [ abs(self.ss_max_force_err_ang[0]), 
                                        abs(self.ss_max_force_err_ang[1]) ] )
        self.normalized_force_error_magnitude = max( [  abs(self.ss_max_force_err_abs[0]), 
                                                        abs(self.ss_max_force_err_abs[1]) ] ) / self.ss_avg_force_magnitude

def pyplot_clear(axeses):
    # self.fig_main.clf()
    # axeses = self.axeses
    # for ax in [axeses[0][0],axeses[0][1],axeses[1][0],axeses[1][1],axeses[2][0],axeses[2][1]]:
    for ax in [axeses[0][0],axeses[0][1],axeses[1][0],axeses[1][1]]:
        ax.cla()
        ax.grid()
    ax = axeses[0][0]; ax.set_xlabel('(a)',fontsize=14.5); ax.set_ylabel('Torque [Nm]',fontsize=14.5)
    ax = axeses[0][1]; ax.set_xlabel('(b)',fontsize=14.5); ax.set_ylabel('Force Amplitude [N]',fontsize=14.5)
    ax = axeses[1][0]; ax.set_xlabel('Time [s]\n(c)',fontsize=14.5); ax.set_ylabel('Norm. Force Error Mag. [%]',fontsize=14.5)
    ax = axeses[1][1]; ax.set_xlabel('Time [s]\n(d)',fontsize=14.5); ax.set_ylabel('Force Error Angle [deg]',fontsize=14.5)



# IEMDC 2019
def read_csv_results_4_comparison__transient(study_name, path_prefix):
    print('look into:', path_prefix)

    # Torque
    basic_info = []
    time_list = []
    TorCon_list = []
    with open(path_prefix + study_name + '_torque.csv', 'r') as f:
        count = 0
        for row in csv_row_reader(f):
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
        for row in csv_row_reader(f):
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
        for row in csv_row_reader(f):
            count +=1
            if count<=8:
                if 'Time' in row[0]:
                    for key in row:
                        key_list.append(key)
                        Current_dict[key] = []
                else:
                    continue
            else:
                for ind, val in enumerate(row):
                    Current_dict[key_list[ind]].append(float(val))

    dm = data_manager()
    dm.basic_info     = basic_info
    dm.time_list      = time_list
    dm.TorCon_list    = TorCon_list
    dm.ForConX_list   = ForConX_list
    dm.ForConY_list   = ForConY_list
    dm.ForConAbs_list = ForConAbs_list
    dm.Current_dict   = Current_dict
    dm.key_list       = key_list
    return dm

def read_csv_results_4_comparison_eddycurrent(study_name, path_prefix):
    # path_prefix = self.dir_csv_output_folder
    # Torque
    TorCon_list = []
    with open(path_prefix + study_name + '_torque.csv', 'r') as f:
        count = 0
        for row in csv_row_reader(f):
            count +=1
            if count<=5:
                continue
            else:
                TorCon_list.append(float(row[1]))

    # Force
    ForConX_list = []
    ForConY_list = []
    with open(path_prefix + study_name + '_force.csv', 'r') as f:
        count = 0
        for row in csv_row_reader(f):
            count +=1
            if count<=5:
                continue
            else:
                ForConX_list.append(float(row[1]))
                ForConY_list.append(float(row[2]))
    ForConAbs_list = np.sqrt(np.array(ForConX_list)**2 + np.array(ForConY_list)**2 )

    dm = data_manager()
    dm.basic_info     = None
    dm.time_list      = None
    dm.TorCon_list    = TorCon_list
    dm.ForConX_list   = ForConX_list
    dm.ForConY_list   = ForConY_list
    dm.ForConAbs_list = ForConAbs_list
    return dm

def collect_jmag_Tran2TSSProlong_results(im_variant, path_prefix, fea_config_dict, axeses, femm_solver_data=None):

    data_results = []

    ################################################################
    # TranRef
    ################################################################   
    study_name = im_variant.get_individual_name() + 'Tran2TSS'
    dm = read_csv_results_4_comparison__transient(study_name, path_prefix)
    basic_info, time_list, TorCon_list, ForConX_list, ForConY_list, ForConAbs_list = dm.unpack()
    end_time = time_list[-1]
    
    if fea_config_dict['number_cycles_prolonged'] != 0:
        sfv = suspension_force_vector(ForConX_list, ForConY_list, range_ss=400) # samples in the tail that are in steady state
        str_results, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle = \
            add_plots( axeses, dm,
                          title='TranRef',#tran_study_name,
                          label='Reference Transient FEA',
                          zorder=8,
                          time_list=time_list,
                          sfv=sfv,
                          torque=TorCon_list,
                          range_ss=sfv.range_ss)
        str_results += '\n\tbasic info:' +   ''.join(  [str(el) for el in basic_info])
        # print str_results
        data_results.extend([torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle])
    else:
        data_results.extend([0, 0, 0, 0, 0])
        

    ################################################################
    # Tran2TSS
    ################################################################
    number_of_total_steps = 1 + fea_config_dict['number_of_steps_1stTTS'] + fea_config_dict['number_of_steps_2ndTTS']
    time_list, TorCon_list, ForConX_list, ForConY_list, ForConAbs_list = time_list[:number_of_total_steps], \
                                                                        TorCon_list[:number_of_total_steps], \
                                                                        ForConX_list[:number_of_total_steps], \
                                                                        ForConY_list[:number_of_total_steps], \
                                                                        ForConAbs_list[:number_of_total_steps]

    sfv = suspension_force_vector(ForConX_list, ForConY_list, range_ss=fea_config_dict['number_of_steps_2ndTTS']) # samples in the tail that are in steady state
    str_results, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle = \
        add_plots( axeses, dm,
                      title='Tran2TSS',#tran_study_name,
                      label='Transient FEA w/ 2 Time Step Sections',
                      zorder=8,
                      time_list=time_list,
                      sfv=sfv,
                      torque=TorCon_list,
                      range_ss=sfv.range_ss)
    str_results += '\n\tbasic info:' +   ''.join(  [str(el) for el in basic_info])
    # print str_results
    data_results.extend([torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle])


    ################################################################
    # FEMM
    ################################################################
    if femm_solver_data is not None:
        study_name = 'FEMM'
        rotor_position_in_deg = femm_solver_data[0]*0.1 
        time_list = rotor_position_in_deg/180.*math.pi / im_variant.Omega
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

        sfv = suspension_force_vector(femm_force_x, femm_force_y, range_ss=len(rotor_position_in_deg)) # samples in the tail that are in steady state
        str_results, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle = \
            add_plots( axeses, dm,
                  title=study_name,
                  label='Static FEA', #'StaticFEAwiRR',
                  zorder=3,
                  time_list=time_list,
                  sfv=sfv,
                  torque=femm_torque,
                  range_ss=sfv.range_ss,
                  alpha=0.5) 
        # print str_results
        data_results.extend([torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle])


    # # for easy access to codes
    # machine_results = [torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle]
    # str_machine_results = ','.join('%g'%(el) for el in machine_results if el is not None) # note that femm_loss_list can be None called by release_design.py
    # return str_results, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle

    return data_results


class data_manager(object):

    def __init__(self):
        self.basic_info = []
        self.time_list = []
        self.TorCon_list = []
        self.ForConX_list = []
        self.ForConY_list = []
        self.ForConAbs_list = []

        self.jmag_loss_list = None
        self.femm_loss_list = None

    def unpack(self):
        return self.basic_info, self.time_list, self.TorCon_list, self.ForConX_list, self.ForConY_list, self.ForConAbs_list

    def terminal_voltage(self, which='4W'): # 2A 2B 2C 4A 4B 4C
        return self.Current_dict['Terminal%s [Case 1]'%(which)]
        # 端点电压是相电压吗？应该是，我们在中性点设置了地电位

    def circuit_current(self, which='4W'): # 2A 2B 2C 4A 4B 4C
        return self.Current_dict['CircuitCoil%s'%(which)]

    def get_voltage_and_current(self, number_of_steps_2ndTTS):

        # 4C <- the C-phase of the 4 pole winding
        mytime  = self.Current_dict['Time(s)'][-number_of_steps_2ndTTS:]
        voltage =      self.terminal_voltage()[-number_of_steps_2ndTTS:]
        current =       self.circuit_current()[-number_of_steps_2ndTTS:]

        # if len(mytime) > len(voltage):
        #     mytime = mytime[:len(voltage)]

        # print len(mytime), len(voltage), number_of_steps_2ndTTS
        # print len(mytime), len(voltage)
        # print len(mytime), len(voltage)

        # for access to plot
        self.myvoltage = voltage
        self.mycurrent = current
        self.mytime    = mytime

    def power_factor(self, number_of_steps_2ndTTS, targetFreq=1e3, numPeriodicalExtension=1000):
        # number_of_steps_2ndTTS: steps corresponding to half the period 

        # for key, val in self.Current_dict.iteritems():
        #     if 'Terminal' in key:
        #         print key, val
        # quit()

        self.get_voltage_and_current(number_of_steps_2ndTTS)
        mytime  = self.mytime
        voltage = self.myvoltage
        current = self.mycurrent

        # from pylab import *
        # print len(mytime), len(voltage), len(current)
        # figure()
        # plot(mytime, voltage)
        # plot(mytime, current)
        # show()
        power_factor, u, i, phase_diff_ui = compute_power_factor_from_half_period(voltage, current, mytime, targetFreq=targetFreq, numPeriodicalExtension=numPeriodicalExtension)
        self.ui_info = [power_factor, u, i, phase_diff_ui]
        return power_factor

def csv_row_reader(handle):
    from csv import reader
    read_iterator = reader(handle, skipinitialspace=True)
    return whole_row_reader(read_iterator)        

def whole_row_reader(reader):
    for row in reader:
        yield row[:]



def get_copper_loss_Bolognani(stator_slot_area, rotor_slot_area=None, STATOR_SLOT_FILL_FACTOR=0.5, ROTOR_SLOT_FILL_FACTOR=1.0, TEMPERATURE_OF_COIL=75, copper_loss_parameters=None): 
    # make sure these two values 
    # space_factor_kCu = SLOT_FILL_FACTOR in Pyrhonen09 design
    # space_factor_kAl = 1 in Pyrhonen09 design
    # as well as temperature
    #   TEMPERATURE_OF_COIL: temperature increase, degrees C
    #   Here, by conductor it means the parallel branch (=1) is considered, and number_of_coil_per_slot = 8 (and will always be 8, see example below).
    #   Just for example, if parallel branch is 2, number_of_coil_per_slot will still be 8, even though we cut the motor in half and count there will be 16 wire cross-sections,
    #   because the reduction in resistance due to parallel branch will be considered into the variable resistance_per_coil.
    #   Remember that the number_of_coil_per_slot (in series) is limited by the high back EMF rather than slot area.
    rho_Copper = (3.76*TEMPERATURE_OF_COIL+873)*1e-9/55. # resistivity

    air_gap_length_delta     = copper_loss_parameters[0]*1e-3 # m (including sleeve, for PMSM)

    # http://127.0.0.1:4000/tech/ECCE-2019-Documentation/

    ################################################################
    # Stator Copper Loss 
    ################################################################
    tooth_width_w_t          = copper_loss_parameters[1]*1e-3 # m
    Area_S_slot              = stator_slot_area
    area_copper_S_Cu         = STATOR_SLOT_FILL_FACTOR * Area_S_slot
    a                        = copper_loss_parameters[2]
    zQ                       = copper_loss_parameters[3]
    coil_pitch_yq            = copper_loss_parameters[4]
    Q                        = copper_loss_parameters[5]
    stack_length_m           = 1e-3*copper_loss_parameters[6]
    current_rms_value        = copper_loss_parameters[7] / 1.4142135623730951 # for one phase
    # Area_conductor_Sc        = Area_S_slot * STATOR_SLOT_FILL_FACTOR / zQ

    Js = (current_rms_value/a) * zQ / area_copper_S_Cu # 逆变器电流current_rms_value在流入电机时，
                                                       # 受到并联支路分流，(current_rms_value/a)才是实际导体中流动的电流值，
                                                       # 这样的电流在一个槽内有zQ个，所以Islot=(current_rms_value/a) * zQ
                                                       # 槽电流除以槽内铜的面积，就是电流密度

    stator_inner_diameter_D = 2*(air_gap_length_delta + copper_loss_parameters[8]*1e-3)
    slot_height_h_t = 0.5*(copper_loss_parameters[9] - stator_inner_diameter_D)
    slot_pitch_pps = np.pi * (stator_inner_diameter_D + slot_height_h_t) / Q
    kov = 1.8 # \in [1.6, 2.0]
    end_winding_length_Lew = np.pi*0.5 * (slot_pitch_pps + tooth_width_w_t) + slot_pitch_pps*kov * (coil_pitch_yq - 1)

    Vol_Cu = area_copper_S_Cu * (stack_length_m + end_winding_length_Lew) * Q
    stator_copper_loss = rho_Copper * Vol_Cu * Js**2

    Vol_Cu_along_stack = area_copper_S_Cu * (end_winding_length_Lew) * Q
    stator_copper_loss_along_stack = rho_Copper * Vol_Cu_along_stack * Js**2

    print('Stator current [Arms]:', current_rms_value, 'Js:', Js)

    if rotor_slot_area is not None:
        raise Exception('Not supported')
        ################################################################
        # Rotor Copper Loss (Valid for IM only)
        ################################################################
        tooth_width_w_t          = self.im.design_parameters[2]*1e-3 # m
        Area_S_slot              = rotor_slot_area
        area_copper_S_Cu         = ROTOR_SLOT_FILL_FACTOR * Area_S_slot
        a                        = 1
        zQ                       = 1
        coil_pitch_yq            = self.im.Qr/self.im.DriveW_poles
        Q                        = self.im.Qr
        # the_radius_m             = 1e-3*(self.im.Radius_OuterRotor - self.im.Radius_of_RotorSlot)
        stack_length_m           = 1e-3*self.im.stack_length
        # number_of_phase          = self.im.Qr/self.im.DriveW_poles
        # Ns                       = zQ * self.im.Qr / (2 * number_of_phase * a) # turns in series = 2 for 4 pole winding; = 1 for 2 pole winding
        # density_of_copper        = 8960 
        # k_R                      = 1 # AC resistance factor
        current_rms_value = sum(self.list_rotor_current_amp) / ( 1.4142135623730951 * len(self.list_rotor_current_amp) )
        # Area_conductor_Sc        = Area_S_slot * ROTOR_SLOT_FILL_FACTOR / zQ        
        Jr = (current_rms_value/a) * zQ / area_copper_S_Cu # 逆变器电流current_rms_value在流入电机时，


        rotor_outer_diameter_Dor = 2*(self.im.Radius_OuterRotor*1e-3)
        slot_height_h_t = self.im.rotor_slot_height_h_sr
        slot_pitch_pps = np.pi * (rotor_outer_diameter_Dor - slot_height_h_t) / Q
        kov = 1.6 
        end_winding_length_Lew = np.pi*0.5 * (slot_pitch_pps + tooth_width_w_t) + slot_pitch_pps*kov * (coil_pitch_yq - 1)

        Vol_Cu = area_copper_S_Cu * (stack_length_m + end_winding_length_Lew) * Q
        rotor_copper_loss = rho_Copper * Vol_Cu * Jr**2

        Vol_Cu_along_stack = area_copper_S_Cu * (end_winding_length_Lew) * Q
        rotor_copper_loss_along_stack = rho_Copper * Vol_Cu_along_stack * Jr**2
        print('Rotor current [Arms]:', current_rms_value, 'Jr:', Jr)
    else:
        rotor_copper_loss, rotor_copper_loss_along_stack, Jr = 0, 0, 0

    return stator_copper_loss, rotor_copper_loss, stator_copper_loss_along_stack, rotor_copper_loss_along_stack, Js, Jr, Vol_Cu



def read_csv_results_4_general_purpose(study_name, path_prefix, fea_config_dict, femm_solver, machine_type=None, acm_variant=None):
    # Read TranFEAwi2TSS results
    
    # logging.getLogger(__name__).debug('Look into: ' + path_prefix)

    # Torque
    basic_info = []
    time_list = []
    TorCon_list = []
    with open(path_prefix + study_name + '_torque.csv', 'r') as f:
        count = 0
        for row in csv_row_reader(f):
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
        for row in csv_row_reader(f):
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
        for row in csv_row_reader(f):
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
    if fea_config_dict['delete_results_after_calculation'] == False:
        # file name is by individual_name like ID32-2-4_EXPORT_CIRCUIT_VOLTAGE.csv rather than ID32-2-4Tran2TSS_circuit_current.csv
        fname = path_prefix + study_name[:-8] + "_EXPORT_CIRCUIT_VOLTAGE.csv"
        # print 'Terminal Voltage - look into:', fname
        with open(fname, 'r') as f:
            count = 0
            for row in csv_row_reader(f):
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
        for row in csv_row_reader(f):
            count +=1
            if 'IM' in machine_type:
                if count>8:
                    rotor_iron_loss = float(row[2]) # Rotor Core
                    stator_iron_loss = float(row[3]) # Stator Core
                    print('Iron loss:', stator_iron_loss, rotor_iron_loss)
                    break
            elif 'PMSM' in machine_type:
                if count>7:
                    print('This should be 0:', float(row[0]))
                    rotor_iron_loss = float(row[1]) # Rotor Core
                    stator_iron_loss = float(row[4]) # Stator Core
                    print('Iron loss:', stator_iron_loss, rotor_iron_loss)
                    break                    
    with open(path_prefix + study_name + '_joule_loss_loss.csv', 'r') as f:
        count = 0
        for row in csv_row_reader(f):
            count +=1
            if 'IM' in machine_type:
                if count>8:
                    rotor_eddycurrent_loss  = float(row[2]) # Rotor Core
                    stator_eddycurrent_loss = float(row[3]) # Stator Core
                    print('Eddy current loss:', stator_eddycurrent_loss, rotor_eddycurrent_loss)
                    break
            elif 'PMSM' in machine_type:
                if count>7:
                    rotor_eddycurrent_loss  = float(row[1]) # Rotor Core
                    stator_eddycurrent_loss = float(row[4]) # Stator Core
                    print('Eddy current loss:', stator_eddycurrent_loss, rotor_eddycurrent_loss)
                    break                    
    with open(path_prefix + study_name + '_hysteresis_loss_loss.csv', 'r') as f:
        count = 0
        for row in csv_row_reader(f):
            count +=1
            if 'IM' in machine_type:
                if count>8:
                    rotor_hysteresis_loss = float(row[2]) # Rotor Core
                    stator_hysteresis_loss = float(row[3]) # Stator Core
                    print('Hysteresis loss:', stator_hysteresis_loss, rotor_hysteresis_loss)
                    break
            elif 'PMSM' in machine_type:
                if count>7:
                    rotor_hysteresis_loss  = float(row[1]) # Rotor Core
                    stator_hysteresis_loss = float(row[4]) # Stator Core
                    print('Hysteresis loss:', stator_hysteresis_loss, rotor_hysteresis_loss)
                    break

    # Joule Loss (Copper and Magnet)
    rotor_Joule_loss_list = []
    with open(path_prefix + study_name + '_joule_loss.csv', 'r') as f:
        count = 0
        for row in csv_row_reader(f):
            count +=1
            if 'IM' in machine_type:
                if count == 8:
                    headers = row
                    for idx_coil, h in enumerate(headers): # on server there are 3 air regions... while on PC there are 2...
                        if 'Coil' in h:
                            break

                if count>8:
                    if count==8+1:
                        if 'Coil' not in headers[idx_coil]:
                            print(headers)
                            raise Exception('Error when load csv data for Coil.')
                        stator_copper_loss = float(row[idx_coil]) # Coil # it is the same over time, this value does not account for end coil

                    if 'Cage' not in headers[idx_coil-1]:
                        print(headers)
                        raise Exception('Error when load csv data for Cage.')
                    rotor_Joule_loss_list.append(float(row[idx_coil-1])) # Cage

            elif 'PMSM' in machine_type:
                if count == 7: # 少一个slip变量，所以不是8，是7。
                    headers = row
                    for idx_coil, h in enumerate(headers):
                        if 'Coil' in h:
                            break

                if count>7:
                    if count==7+1:
                        if 'Coil' not in headers[idx_coil]:
                            print(headers)
                            raise Exception('Error when load csv data for Coil.')
                        stator_copper_loss = float(row[idx_coil]) # Coil # it is the same over time, this value does not account for end coil

                    if 'Magnet' not in headers[idx_coil-1]:
                        print(headers)
                        raise Exception('Error when load csv data for Magnet.')
                    rotor_Joule_loss_list.append(float(row[idx_coil-1])) # Magnet

    # use the last 1/4 period data to compute average copper loss of Tran2TSS rather than use that of Freq study
    effective_part = rotor_Joule_loss_list[-int(0.5*fea_config_dict['number_of_steps_2ndTTS']):] # number_of_steps_2ndTTS = steps for half peirod
    # print(rotor_Joule_loss_list)
    # print(effective_part)
    # quit()
    rotor_Joule_loss = sum(effective_part) / len(effective_part)
    if 'PMSM' in machine_type:
        print('Magnet Joule loss:', rotor_Joule_loss)

    if fea_config_dict['jmag_run_list'][0] == 0 and femm_solver is not None:
        # blockPrint()
        try:
            # convert rotor current results (complex number) into its amplitude
            femm_solver.list_rotor_current_amp = [abs(el) for el in femm_solver.vals_results_rotor_current] # el is complex number
            # settings not necessarily be consistent with Pyrhonen09's design: , STATOR_SLOT_FILL_FACTOR=0.5, ROTOR_SLOT_FILL_FACTOR=1., TEMPERATURE_OF_COIL=75

            # slot_area_utilizing_ratio = (acm_variant.DriveW_CurrentAmp + acm_variant.BeariW_CurrentAmp) / acm_variant.CurrentAmp_per_phase
            # if slot_area_utilizing_ratio < 1:
            #     print('Heads up! slot_area_utilizing_ratio is', slot_area_utilizing_ratio, 'which means you are simulating a separate winding? If not, contrats--you found a bug...')
            #     print('DW, BW, Total:', acm_variant.DriveW_CurrentAmp, acm_variant.BeariW_CurrentAmp, acm_variant.CurrentAmp_per_phase)

            _s, _r, _sAlongStack, _rAlongStack, _Js, _Jr = femm_solver.get_copper_loss_pyrhonen(acm_variant.slot_area_utilizing_ratio*femm_solver.stator_slot_area, 
                                                                                                                          femm_solver.rotor_slot_area, 
                                                                                                                          total_CurrentAmp=acm_variant.DriveW_CurrentAmp+acm_variant.BeariW_CurrentAmp)
            s, r, sAlongStack, rAlongStack, Js, Jr, Vol_Cu = femm_solver.get_copper_loss_Bolognani(acm_variant.slot_area_utilizing_ratio*femm_solver.stator_slot_area, 
                                                                                                                          femm_solver.rotor_slot_area, 
                                                                                                                          total_CurrentAmp=acm_variant.DriveW_CurrentAmp+acm_variant.BeariW_CurrentAmp)

            msg1 = 'Pyrhonen : %g, %g | %g, %g | %g, %g ' % (_s, _r, _sAlongStack, _rAlongStack, _Js, _Jr) 
            msg2 = 'Bolognani: %g, %g | %g, %g | %g, %g ' % (s, r, sAlongStack, rAlongStack, Js, Jr) 
            logger = logging.getLogger(__name__)
            logger.debug(msg1)
            logger.debug(msg2)
        except Exception as e:
            raise e
        # enablePrint()
    else:
        copper_loss_parameters = [acm_variant.sleeve_length + acm_variant.fixed_air_gap_length,
                             acm_variant.mm_w_st,
                             acm_variant.wily.number_parallel_branch,
                             acm_variant.DriveW_zQ,
                             acm_variant.wily.coil_pitch,
                             acm_variant.Q,
                             acm_variant.stack_length,
                             acm_variant.DriveW_CurrentAmp + acm_variant.BeariW_CurrentAmp, # total current amplitude
                             acm_variant.Radius_OuterRotor,
                             acm_variant.stator_yoke_diameter_Dsyi
                             ]
        # slot_area_utilizing_ratio = (acm_variant.DriveW_CurrentAmp + acm_variant.BeariW_CurrentAmp) / acm_variant.CurrentAmp_per_phase
        # if slot_area_utilizing_ratio < 1:
        #     print('Heads up! slot_area_utilizing_ratio is', slot_area_utilizing_ratio, 'which means you are simulating a separate winding? If not, contrats--you found a bug...')
        #     print('DW, BW, Total:', acm_variant.DriveW_CurrentAmp, acm_variant.BeariW_CurrentAmp, acm_variant.CurrentAmp_per_phase)
        s, r, sAlongStack, rAlongStack, Js, Jr, Vol_Cu = get_copper_loss_Bolognani(acm_variant.slot_area_utilizing_ratio*acm_variant.coils.mm2_slot_area*1e-6, copper_loss_parameters=copper_loss_parameters)
        # s, r, sAlongStack, rAlongStack, Js, Jr = 0, 0, 0, 0, 0, 0

    dm = data_manager()
    dm.basic_info     = basic_info
    dm.time_list      = time_list
    dm.TorCon_list    = TorCon_list
    dm.ForConX_list   = ForConX_list
    dm.ForConY_list   = ForConY_list
    dm.ForConAbs_list = ForConAbs_list
    dm.Current_dict   = Current_dict
    dm.key_list       = key_list
    dm.jmag_loss_list = [   stator_copper_loss, 
                            rotor_Joule_loss, 
                            stator_iron_loss+rotor_iron_loss, 
                            stator_eddycurrent_loss+rotor_eddycurrent_loss, 
                            stator_hysteresis_loss+rotor_hysteresis_loss ]
    dm.femm_loss_list = [s, r, sAlongStack, rAlongStack, Js, Jr ]
    dm.Vol_Cu         = Vol_Cu
    return dm

def check_csv_results_4_general_purpose(study_name, path_prefix, returnBoolean=False):
    # Check frequency analysis results

    print('Check:', path_prefix + study_name + '_torque.csv')

    if not os.path.exists(path_prefix + study_name + '_torque.csv'):
        if returnBoolean == False:
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

        # fitness_in_physics_data = []
        with open(path_prefix + study_name + '_torque.csv', 'r') as f: 
            for ind, row in enumerate(csv_row_reader(f)):
                if ind >= 5:
                    try:
                        float(row[0])
                    except:
                        continue
                    l_slip_freq.append(float(row[0]))
                    l_TorCon.append(float(row[1]))

        with open(path_prefix + study_name + '_force.csv', 'r') as f: 
            for ind, row in enumerate(csv_row_reader(f)):
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

        index, breakdown_torque = get_index_and_max(l_TorCon)
        slip_freq_breakdown_torque = l_slip_freq[index]
        return slip_freq_breakdown_torque, breakdown_torque, breakdown_force
    except NameError as e:
        logger = logging.getLogger(__name__)
        logger.error('No CSV File Found.', exc_info=True)
        raise e


# https://dsp.stackexchange.com/questions/11513/estimate-frequency-and-peak-value-of-a-signals-fundamental
#define N_SAMPLE ((long int)(1.0/(0.1*TS))) // Resolution 0.1 Hz = 1 / (N_SAMPLE * TS)
class Goertzel_Data_Struct(object):
    """docstring for Goertzel_Data_Struct"""
    def __init__(self, id=None, ):
        self.id = id    
        self.bool_initialized = False
        self.sine = None
        self.cosine = None
        self.coeff = None
        self.scalingFactor = None
        self.q = None
        self.q2 = None
        self.count = None
        self.k = None #// k is the normalized target frequency 
        self.real = None
        self.imag = None
        self.accumSquaredData = None
        self.ampl = None
        self.phase = None


    # /************************************************
    #  * Real time implementation to avoid the array of input double *data[]
    #  * with Goertzel Struct to store the variables and the output values
    #  *************************************************/
    def goertzel_realtime(gs, targetFreq, numSamples, samplingRate, data):
        # gs is equivalent to self
        try:
            len(data)
        except:
            pass
        else:
            raise Exception('This is for real time implementation of Goertzel, hence data must be a scalar rather than array-like object.')

        if gs.bool_initialized == False:
            gs.bool_initialized = True

            gs.count = 0
            gs.k = (0.5 + ((numSamples * targetFreq) / samplingRate))
            omega = (2.0 * math.pi * gs.k) / numSamples
            gs.sine = sin(omega)
            gs.cosine = cos(omega)
            gs.coeff = 2.0 * gs.cosine
            gs.q1=0
            gs.q2=0
            gs.scalingFactor = 0.5 * numSamples
            gs.accumSquaredData = 0.0

        q0 = gs.coeff * gs.q1 - gs.q2 + data
        gs.q2 = gs.q1
        gs.q1 = q0 # // q1 is the newest output vk[N], while q2 is the last output vk[N-1].

        gs.accumSquaredData += data*data

        gs.count += 1
        if gs.count>=numSamples:
            # // calculate the real and imaginary results with scaling appropriately
            gs.real = (gs.q1 * gs.cosine - gs.q2) / gs.scalingFactor #// inspired by the python script of sebpiq
            gs.imag = (gs.q1 * gs.sine) / gs.scalingFactor

            # // reset
            gs.bool_initialized = False
            return True
        else:
            return False

    def goertzel_offline(gs, targetFreq, samplingRate, data_list):
        # gs is equivalent to self

        numSamples = len(data_list)

        if gs.bool_initialized == False:
            gs.bool_initialized = True

            gs.count = 0
            gs.k = (0.5 + ((numSamples * targetFreq) / samplingRate))
            omega = (2.0 * math.pi * gs.k) / numSamples
            gs.sine = math.sin(omega)
            gs.cosine = math.cos(omega)
            gs.coeff = 2.0 * gs.cosine
            gs.q1=0
            gs.q2=0
            gs.scalingFactor = 0.5 * numSamples
            gs.accumSquaredData = 0.0

        for data in data_list:
            q0 = gs.coeff * gs.q1 - gs.q2 + data
            gs.q2 = gs.q1
            gs.q1 = q0 # // q1 is the newest output vk[N], while q2 is the last output vk[N-1].

            gs.accumSquaredData += data*data

            gs.count += 1
            if gs.count>=numSamples:

                # // calculate the real and imaginary results with scaling appropriately
                gs.real = (gs.q1 * gs.cosine - gs.q2) / gs.scalingFactor #// inspired by the python script of sebpiq
                gs.imag = (gs.q1 * gs.sine) / gs.scalingFactor

                # // reset
                gs.bool_initialized = False
                return True
        print(data_list)
        print(gs.count)
        print(numSamples)
        return None

def compute_power_factor_from_half_period(voltage, current, mytime, targetFreq=1e3, numPeriodicalExtension=1000): # 目标频率默认是1000Hz

    gs_u = Goertzel_Data_Struct("Goertzel Struct for Voltage\n")
    gs_i = Goertzel_Data_Struct("Goertzel Struct for Current\n")

    TS = mytime[-1] - mytime[-2]

    if type(voltage)!=type([]):
        voltage = voltage.tolist() + (-voltage).tolist()
        current = current.tolist() + (-current).tolist()
    else:
        voltage = voltage + [-el for el in voltage]
        current = current + [-el for el in current]

    voltage *= numPeriodicalExtension
    current *= numPeriodicalExtension

    N_SAMPLE = len(voltage)
    gs_u.goertzel_offline(targetFreq, 1./TS, voltage)
    gs_i.goertzel_offline(targetFreq, 1./TS, current)

    gs_u.ampl = math.sqrt(gs_u.real*gs_u.real + gs_u.imag*gs_u.imag) 
    gs_u.phase = math.atan2(gs_u.imag, gs_u.real)

    gs_i.ampl = math.sqrt(gs_i.real*gs_i.real + gs_i.imag*gs_i.imag) 
    gs_i.phase = math.atan2(gs_i.imag, gs_i.real)

    print('gs_u.ampl', gs_u.ampl)
    print('gs_i.ampl', gs_i.ampl)

    phase_difference_in_deg = ((gs_i.phase-gs_u.phase)/math.pi*180)
    power_factor = math.cos(gs_i.phase-gs_u.phase)
    return power_factor, gs_u.ampl, gs_i.ampl, phase_difference_in_deg
    




def max_indices_2(arr, k):
    '''
    Returns the indices of the k first largest elements of arr
    (in descending order in values)
    '''
    assert k <= arr.size, 'k should be smaller or equal to the array size'
    arr_ = arr.astype(float)  # make a copy of arr
    max_idxs = []
    for _ in range(k):
        max_element = np.max(arr_)
        if np.isinf(max_element):
            break
        else:
            idx = np.where(arr_ == max_element)
        max_idxs.append(idx)
        arr_[idx] = -np.inf
    return max_idxs

def min_indices(arr, k):
    if type(arr) == type([]):
        arr = np.array(arr)
    indices = np.argsort(arr)[:k]
    items   = arr[indices] # arr and indices must be np.array
    return indices, items

def max_indices(arr, k):
    if type(arr) == type([]):
        arr = np.array(arr)
    arr_copy = arr[::]
    indices = np.argsort(-arr)[:k]
    items   = arr_copy[indices]
    return indices, items



def autolabel(ax, rects, xpos='center', bias=0.0, textfont=None):
    """
    Attach a text label above each bar in *rects*, displaying its height.

    *xpos* indicates which side to place the text w.r.t. the center of
    the bar. It can be one of the following {'center', 'right', 'left'}.
    """

    xpos = xpos.lower()  # normalize the case of the parameter
    ha = {'center': 'center', 'right': 'left', 'left': 'right'}
    offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height+bias,
                '%.2f'%(height), ha=ha[xpos], va='bottom', rotation=90, **textfont)

def efficiency_at_50kW(total_loss):
    return 50e3 / (np.array(total_loss) + 50e3)  # 用50 kW去算才对，只是这样转子铜耗数值会偏大哦。 效率计算：机械功率/(损耗+机械功率)        

def use_weights(which='O1'):
    if which == 'O1':
        return [ 1, 0.1,   1, 0.1, 0.1,   0 ]
    if which == 'O2':
        return [ 1,1,1,1,1,  0 ]
    if which == 'O3':
        return [ 1,0,1,0,0,  1 ]
    if which == 'O1R':
        return [ 1, 0.2,   1, 0.1, 0.1,   0 ]
    if which == 'O4':
        return [ 2,2,1,1,1,  0 ]
    return None

def compute_list_cost(weights, rotor_volume, rotor_weight, torque_average, normalized_torque_ripple, ss_avg_force_magnitude, normalized_force_error_magnitude, force_error_angle, jmag_loss_list, femm_loss_list, power_factor, total_loss):
    # O2
    # weights = [ 1, 1.0,   1, 1.0, 1.0,   0 ]
    # O1
    # weights = [ 1, 0.1,   1, 0.1, 0.1,   0 ]
    list_cost = [   30e3 / ( torque_average/rotor_volume ),
                    normalized_torque_ripple         *  20, 
                    1.0 / ( ss_avg_force_magnitude/rotor_weight ),
                    normalized_force_error_magnitude *  20, 
                    force_error_angle                * 0.2, # [deg] 
                    total_loss                       / 2500. ] 
    cost_function = np.dot(np.array(list_cost), np.array(weights))
    return cost_function, list_cost

    # The weight is [TpRV=30e3, FpRW=1, Trip=50%, FEmag=50%, FEang=50deg, eta=sqrt(10)=3.16]
    # which means the FEang must be up to 50deg so so be the same level as TpRV=30e3 or FpRW=1 or eta=316%
    # list_weighted_cost = [  30e3 / ( torque_average/rotor_volume ),
    #                         1.0 / ( ss_avg_force_magnitude/rotor_weight ),
    #                         normalized_torque_ripple         *   2, #       / 0.05 * 0.1
    #                         normalized_force_error_magnitude *   2, #       / 0.05 * 0.1
    #                         force_error_angle * 0.2          * 0.1, # [deg] /5 deg * 0.1 is reported to be the base line (Yegu Kang) # force_error_angle is not consistent with Yegu Kang 2018-060-case of TFE
    #                         2*total_loss/2500., #10 / efficiency**2,
    #                         im_variant.thermal_penalty ] # thermal penalty is evaluated when drawing the model according to the parameters' constraints (if the rotor current and rotor slot size requirement does not suffice)
    # cost_function = sum(list_weighted_cost)


def fobj_scalar(torque_average, ss_avg_force_magnitude, normalized_torque_ripple, normalized_force_error_magnitude, force_error_angle, total_loss, 
                weights=None, rotor_volume=None, rotor_weight=None):

    list_cost = [   30e3 / ( torque_average/rotor_volume ),
                    normalized_torque_ripple         *  20, #       / 0.05 * 0.1
                    1.0 / ( ss_avg_force_magnitude/rotor_weight ),
                    normalized_force_error_magnitude *  20, #       / 0.05 * 0.1
                    force_error_angle                * 0.2, # [deg] /5 deg * 0.1 is reported to be the base line (Yegu Kang) # force_error_angle is not consistent with Yegu Kang 2018-060-case of TFE
                    total_loss                       / 2500. ] #10 / efficiency**2,
    cost_function = np.dot(np.array(list_cost), np.array(weights))
    return cost_function

def fobj_list(l_torque_average, l_ss_avg_force_magnitude, l_normalized_torque_ripple, l_normalized_force_error_magnitude, l_force_error_angle, l_total_loss,
                weights=None, rotor_volume=None, rotor_weight=None):

    l_cost_function = []
    for torque_average, ss_avg_force_magnitude, normalized_torque_ripple, normalized_force_error_magnitude, force_error_angle, total_loss in zip(l_torque_average, l_ss_avg_force_magnitude, l_normalized_torque_ripple, l_normalized_force_error_magnitude, l_force_error_angle, l_total_loss):
        list_cost = [   30e3 / ( torque_average/rotor_volume ),
                        normalized_torque_ripple         *  20,
                        1.0 / ( ss_avg_force_magnitude/rotor_weight ),
                        normalized_force_error_magnitude *  20,
                        force_error_angle                * 0.2,
                        total_loss                       / 2500. ]
        cost_function = np.dot(np.array(list_cost), np.array(weights))
        l_cost_function.append(cost_function)
    return np.array(l_cost_function)



class SwarmDataAnalyzer(object):
    """docstring for SwarmDataAnalyzer"""
    def __init__(self, sw, spec, dir_run=None, run_integer=None, bool_sensitivity_analysis=True):

        number_of_variant = sw.fea_config_dict['local_sensitivity_analysis_number_of_variants'] + 1
        self.sw = sw
        self.spec = spec

        self.dir_run = dir_run
        self.run_integer = run_integer
        with open(dir_run+'swarm_data.txt', 'r') as f:
            self.buf = f.readlines()[1:]

        self.number_of_designs = len(self.buf) / 21 # 此21是指每个个体的结果占21行，非彼20+1哦。

        if bool_sensitivity_analysis:
            if self.number_of_designs == number_of_variant*7:
                print('These are legacy results without the initial design within the swarm_data.txt. ')

            elif self.number_of_designs == number_of_variant*7 + 1:
                print('Initial design is among the pop and used as reference design.')
                self.reference_design = self.buf[:21] # 此21是指每个个体的结果占21行，非彼20+1哦。
                self.buf = self.buf[21:] # 此21是指每个个体的结果占21行，非彼20+1哦。
            else:
                raise Exception('Please remove duplicate results in swarm_data.txt for run#', run_integer)

        # now we have the actual pop size
        self.buf_length = len(self.buf)
        self.number_of_designs = len(self.buf) / 21

        msg = 'self.buf_length %% 21 = %d, / 21 = %d' % (self.buf_length % 21, self.number_of_designs)
        logger = logging.getLogger(__name__)
        logger.debug(msg)

        self.build_basic_info()

    def design_display_generator(self):
        for i in range(int(self.number_of_designs)):
            yield ''.join(self.buf[i*21:(1+i)*21])

    def design_parameters_generator(self):
        for i in range(int(self.number_of_designs)):
            yield [float(el) for el in self.buf[i*21:(1+i)*21][5].split(',')]

    def list_generations(self):

        the_dict = {}
        for i in range(int(self.number_of_designs)):
            the_row = [float(el) for el in self.buf[i*21:(1+i)*21][2].split(',')]        
            the_dict[int(the_row[0])] = (the_row[1], the_row[2])
            # the_row[0] # generation
            # the_row[1] # individual
            # the_row[2] # cost_function
        return the_dict #每代只留下最后一个个体的结果，因为字典不能重复key

    def list_cost_function(self):
        l = []
        for i in range(self.number_of_designs):
            l.append( float( self.buf[i*21:(1+i)*21][2].split(',')[2] ) )
        return l

    def find_individual(self, generation_index, individual_index):
        # for i in range(self.number_of_designs):
        for i in range(self.number_of_designs)[::-1]:
            # print i
            the_row = [int(el) for el in self.buf[i*21:(1+i)*21][2].split(',')[:2]]
            if the_row[0] == generation_index:
                if the_row[1] == individual_index:
                    # found it, return the design parameters and cost_function
                    return [float(el) for el in self.buf[i*21:(1+i)*21][5].split(',')], float(self.buf[i*21:(1+i)*21][2].split(',')[2])
        return None, None

    def get_best_generation(self, popsize=50, generator=None, returnMore=False):

        if generator is None:
            generator = self.design_parameters_generator()

        cost = self.list_cost_function()
        indices, items = min_indices(cost, popsize)
        print(indices)
        print(items)

        # for ind, el in enumerate(cost):
        #     print ind, el,
        # quit()

        # this is in wrong order as inidices, so use dict instead
        # gen_best = [design for index, design in enumerate(generator) if index in indices]
        gen_best = {}
        for index, design in enumerate(generator):
            if index in indices:
                gen_best[index] = design
        gen_best = [gen_best[index] for index in indices] # now it is in order

        if returnMore == False:
            return gen_best
        else:
            return gen_best, indices, items

    def get_list_objective_function(self):
        for i in range(self.number_of_designs):
            individual = self.buf[i*21:(1+i)*21]
            yield     [float(el) for el in individual[3].split(',')] #\
                    # + [float(el) for el in individual[17][16:].split(',')] \
                    # + [float(el) for el in individual[18][16:].split(',')]

    def get_certain_objective_function(self, which):
        for i in range(int(self.number_of_designs)):
            individual = self.buf[i*21:(1+i)*21]
            try:
                yield [float(el) for el in individual[3].split(',')][which]
            except Exception as e:
                print([float(el) for el in individual[3].split(',')], which, i, individual)
                raise e

    def get_windage_loss(self, which=1):
        for i in range(int(self.number_of_designs)):
            individual = self.buf[i*21:(1+i)*21]
            try:
                loc_data_begin = individual[-1].find(':') + 1
                yield [float(el) for el in individual[-1][loc_data_begin:].split(',')][which]
            except Exception as e:
                print([float(el) for el in individual[-1][loc_data_begin:].split(',')], which, i, individual)
                raise e

    def my_population_distribution_plots(self, de_config_dict):
        from pylab import subplots, subplots_adjust, show, legend
        fig, axes = subplots(7, 1, sharex=True, dpi=150, figsize=(12, 6), facecolor='w', edgecolor='k')
        subplots_adjust(left=None, bottom=0.25, right=None, top=None, wspace=0.8, hspace=None)
        list_label = [  'Stator tooth width $w_{st}$         ',
                        'Air gap length $L_g$                ',
                        'Rotor slot open width $w_{ro}$      ',
                        'Rotor tooth width $w_{rt}$          ',
                        'Rotor slot open depth $d_{ro}$      ',
                        'Stator slot open angle $\theta_{so}$',
                        'Stator slot open depth $d_{so}$     ',]
        geo_param = {}
        geo_param[list_label[0]] = []
        geo_param[list_label[1]] = []
        geo_param[list_label[2]] = []
        geo_param[list_label[3]] = []
        geo_param[list_label[4]] = []
        geo_param[list_label[5]] = []
        geo_param[list_label[6]] = []
        for el in self.design_parameters_generator():
            geo_param[list_label[0]].append(el[0])
            geo_param[list_label[1]].append(el[1])
            geo_param[list_label[2]].append(el[2])
            geo_param[list_label[3]].append(el[3])
            geo_param[list_label[4]].append(el[4])
            geo_param[list_label[5]].append(el[5])
            geo_param[list_label[6]].append(el[6])
        from collections import OrderedDict
        ordered_geo_param = [(list_label[0], geo_param[list_label[0]]),
                             (list_label[1], geo_param[list_label[1]]),
                             (list_label[2], geo_param[list_label[2]]),
                             (list_label[3], geo_param[list_label[3]]),
                             (list_label[4], geo_param[list_label[4]]),
                             (list_label[5], geo_param[list_label[5]]),
                             (list_label[6], geo_param[list_label[6]])]
        ordered_geo_param = OrderedDict(ordered_geo_param)

        
        popsize = de_config_dict['popsize']
        no_generations = len(ordered_geo_param[list_label[0]]) // popsize
        list_avg = [[], [], [], [], [], [], []] # != [[]]*7
        for number_current_generation in range(no_generations):
            for ind, key in enumerate(list_label):
                item = ordered_geo_param[key]
                item_in_question = item[number_current_generation*popsize:(number_current_generation+1)*popsize]

                avg = sum(item_in_question) / len(item_in_question)
                list_avg[ind].append(avg)
                # print(list_avg)

                axes[ind].plot( number_current_generation*np.ones(len(item_in_question)), 
                                item_in_question, 
                                label=list_label[ind].strip(), marker='_', color=None, alpha=1.00, lw=0.2 )
        for ind, key in enumerate(list_label):
            axes[ind].plot(list(range(no_generations)), list_avg[ind], color='k', lw=1.5)

        print ('Distribution of total individuals:', popsize*no_generations, 'No. generations', no_generations)
        # legend()
        fig.tight_layout()
        # show()

    def my_scatter_plot(self, x, y, O, fig=None, ax=None, s=15, marker='.', index_list=None): # index_list is for filtered case
        O1, O2 = None, None
        # O is a copy of your list rather than array or the adress of the list
        # O is a copy of your list rather than array or the adress of the list
        # O is a copy of your list rather than array or the adress of the list

        if ax is None or fig is None:
            fig = figure()
            ax = fig.gca()
        scatter_handle = ax.scatter(x, y, c=O, s=s, alpha=0.5, cmap='viridis', marker=marker)
        # print('---------------------For plotting in Matlab:')
        # for a, b in zip(x, y):
        #     print(a, b)
        #     # O2_mix = np.concatenate([[O_ref], O2], axis=0) # # https://stackoverflow.com/questions/46106912/one-colorbar-for-multiple-scatter-plots
        #     # min_, max_ = O2_mix.min(), O2_mix.max()
        #     # scatter(*xy_ref, marker='s', c=O_ref, s=20, alpha=0.75, cmap='viridis')
        #     # clim(min_, max_)
        #     # clim(min_, max_)
        #     # ax.grid()
        # print('---------------------End for plotting in Matlab:')

        if True:
            # Decide which one is the best in this Pareto plot
            best_index, best_O = get_index_and_min(O)
            # best_index, best_O = get_index_and_min(O[:best_index]+O[best_index+1:])
 
            # index_list being not None means filtered data are used. 注意，部分个体被滤除掉了，所以index_list是断断续续的，所以是必需的。
            if index_list is not None:
                print('run_integer:', self.run_integer)
                print('BEST metrics:', best_index, index_list[best_index], best_O)
                print('Best design:', x[best_index], y[best_index])
                print(u'目标函数值对不上，是因为跑完优化以后，我修改了钢的密度，所以后处理中重新计算的目标函数会有所不同。')
                print(u'此外，铜耗按叠长放大，其实是保守估计了，因为端部的铜耗不会因为叠长变化而变化。')

                self.best_design_display = next(itertools.islice(self.design_display_generator(),    index_list[best_index], None))
                self.best_design_denorm = next(itertools.islice(self.design_parameters_generator(), index_list[best_index], None))
                print( self.best_design_display )
                    # for ind, el in enumerate(self.design_display_generator()):
                    #     if ind == index_list[best_index]:
                    #         print(ind, el)
                    #         break
                    # print (next(itertools.islice(self.design_parameters_generator(), index_list[best_index], index_list[best_index]+1)))
                    # print (list(self.design_parameters_generator())[index_list[best_index]])
                data = [float(el) for el in self.best_design_display.split('\n')[3].split(',')]
                print(data)
                # unpacking
                PF, _, torque, Trip, Fmag, Em, Ea, \
                    jmag_stator_copper_loss, jmag_rotor_copper_loss, iron_loss, eddy_loss, hyst_loss, \
                    femm_stator_copper_loss, femm_rotor_copper_loss, windage_loss, total_loss = data
                self.speed_rpm        
                self.Omega            
                self.mec_power        
                self.required_torque  
                self.Radius_OuterRotor
                self.stack_length     
                self.Qs               
                self.Qr               
                self.rotor_volume # 这个体积是对应原来的叠长的！
                self.rotor_weight     
                self.weights_name     
                self.weights_used     
                self.stack_length     
                self.stack_length_max 
                rated_stack_length = self.stack_length / torque * self.required_torque
                rated_total_loss   = total_loss / self.stack_length * rated_stack_length # 这样简单计算，端部的损耗就多计入了。
                eta = 1 - total_loss/self.mec_power
                self.str_best_design_details = f'''
                \\begin{{table*}}[!t]
                  \\caption{{Design Release Data}}
                  \\centering
                    \\begin{{tabular}}{{rl}}
                      \\hline
                      \\hline
                      \\thead{{Item}} &
                      \\thead{{Value}} \\\\
                      \\hline
                      TRV [kNm/$\\rm m^3$]                                  & {torque/self.rotor_volume*1e-3          :.1f} \\\\
                      FRW [N/kg]                                            & {Fmag/self.rotor_weight                 :.1f} \\\\
                      Torque ripple [\\%]                                   & {Trip*100                               :.1f} \\\\
                      Force Error angle [deg]                               & {Ea                                     :.1f} \\\\
                      Force Error magnitude [\\%]                           & {Em*100                                 :.1f} \\\\
                      Efficiency at rated load (include windage loss) [\\%] & {eta*100                                :.1f} \\\\
                      Power Factor at rated load                            & {PF                                     :.2f} \\\\
                      Stator OD [mm]                                        & {2*self.sw.im.Radius_OuterStatorYoke    :.1f} \\\\
                      Rotor OD [mm]                                         & {2*self.sw.im.Radius_OuterRotor         :.1f} \\\\
                      Stack length [mm]                                     & {rated_stack_length                     :.1f} \\\\
                      Rotor volume [$\\rm m^3$]                             & {self.rotor_volume                        :g} \\\\
                      Weight of the rotor [N]                               & {self.rotor_weight                        :g} \\\\
                      Airgap length [mm]                                    & {self.best_design_denorm[0]             :.2f} \\\\
                      Sleeve thickness [mm]                                 & {0                                      :.1f} \\\\
                      Suspension poles                                      & {int(self.sw.im.BeariW_poles)             :g} \\\\
                      Torque poles                                          & {int(self.sw.im.DriveW_poles)             :g} \\\\
                      Motor Electric Frequency [Hz]                         & {self.ExcitationFreqSimulated             :g} \\\\
                      Motor mechanical speed (r/min)                        & {self.speed_rpm                           :g} \\\\
                      Rated Power [W]                                       & {self.mec_power                           :g} \\\\
                      Stator slots                                          & {int(self.Qs)                             :g} \\\\
                      Stator slot fill factor                               & {0.5                                    :.1f} \\\\
                      Rotor slots                                           & {int(self.Qr)                             :g} \\\\
                      Rotor slot fill factor                                & {1.0                                    :.1f} \\\\
                      Stator Conductor current density (Arms/$\\rm mm^2$)   & {self.spec.Js                             :g} \\\\
                      Rotor Conductor current density (Arms/$\\rm mm^2$)    & {self.spec.Jr                             :g} \\\\
                      Lamination Material                                   & {self.spec.Steel                            } \\\\
                      Motor Phase Current (Arms)                            & {self.spec.stator_phase_current_rms     :.1f} \\\\
                      Motor Phase Voltage (Vrms)                            & {self.spec.VoltageRating                :.1f} \\\\
                      Criteria for candidate design:                        & {self.sw.fea_config_dict['use_weights']     } \\\\
                    \\hline
                    \\vspace{{-2.5ex}}
                    \\\\
                \\end{{tabular}}
              \\label{{tab:001}}
              \\vspace{{-3ex}}
            \\end{{table*}}
            '''
            # quit()
            xy_best = (x[best_index], y[best_index])
            handle_best = ax.scatter(*xy_best, s=s*3, marker='s', facecolors='none', edgecolors='r')
            ax.legend((handle_best,), ('Best',))

        return scatter_handle

    def pareto_plot_eta_vs_stack_legnth(self, fig, ax, marker='.', bool_filtered=True):
        ''' Pareto Plot for efficiency v.s. statck length @ rated power '''

        O = fobj_list(  list(self.get_certain_objective_function(2)), #torque_average, 
                        list(self.get_certain_objective_function(4)), #ss_avg_force_magnitude, 
                        list(self.get_certain_objective_function(3)), #normalized_torque_ripple, 
                        list(self.get_certain_objective_function(5)), #normalized_force_error_magnitude, 
                        list(self.get_certain_objective_function(6)), #force_error_angle, 
                        list(self.get_certain_objective_function(15)), #total_loss, 
                        # np.array(  list(self.get_certain_objective_function(9)) ) \
                        # + np.array(list(self.get_certain_objective_function(12))) \
                        # + np.array(list(self.get_certain_objective_function(13))) \
                        # + self.get_windage_loss(), # or windage is the 14
                        weights=self.weights_used, rotor_volume=self.rotor_volume, rotor_weight=self.rotor_weight)
        # print O
        # print np.array(swda.list_cost_function()) - np.array(O) # they are the same
        O = O.tolist()

        # Efficiency vs Rated power stack length
        def get_rated_values(l_torque_average, l_total_loss):
            l_rated_stack_length = []
            l_rated_total_loss = []
            for torque_average, total_loss in zip(l_torque_average, l_total_loss):
                rated_stack_length = self.stack_length / torque_average * self.required_torque
                l_rated_stack_length.append(rated_stack_length)

                rated_total_loss   = total_loss / self.stack_length * rated_stack_length
                l_rated_total_loss.append(rated_total_loss)

            return l_rated_stack_length, l_rated_total_loss

            # a, b = get_rated_values([torque_average], [total_loss])
            # rated_stack_length = a[0]
            # rated_total_loss = b[0]
            # print('stack_length=', stack_length, 'mm, rated_stack_length=', rated_stack_length, 'mm')
            # print('total_loss=', total_loss, 'W, rated_total_loss=', rated_total_loss, 'W')
            # xy_ref = (rated_stack_length, 1 - rated_total_loss/self.mec_power)

        if not bool_filtered: # no filter
            x, y = get_rated_values(list(self.get_certain_objective_function(2)), 
                                    list(self.get_certain_objective_function(15)) )
            y = 1 - np.array(y)/self.mec_power
            y = y.tolist()
            self.my_scatter_plot(x, y, O[::], fig=fig, ax=ax, marker=marker)
        else:
            filtered_torque_list = []
            filtered_loss_list = []
            filtered_O_list = []
            index_list = []
            for torque, loss, err_angle, err_mag, o, index in zip(
                                                        list(self.get_certain_objective_function(2)),
                                                        list(self.get_certain_objective_function(15)),
                                                        list(self.get_certain_objective_function(6)),
                                                        list(self.get_certain_objective_function(5)),
                                                        O,
                                                        list(range(len(O)))):
                if err_angle<5 and err_mag<0.25:
                    # print(err_angle, err_mag)
                    filtered_torque_list.append(torque)
                    filtered_loss_list.append(loss)
                    filtered_O_list.append(o)
                    index_list.append(index)
            print('Max index is', index)

            x, y = get_rated_values(filtered_torque_list, filtered_loss_list)
            filtered_x = []
            filtered_y = []
            filtered_filtered_O_list = []
            filtered_index_list = []
            for x_el, y_el, O_el, index in zip(x, y, filtered_O_list, list(range(len(filtered_O_list)))):
                if x_el < self.stack_length_max*2: ######################################################################## *2
                    filtered_x.append(x_el)
                    filtered_y.append(y_el)
                    filtered_filtered_O_list.append(O_el)
                    filtered_index_list.append(index_list[index])

            filtered_y = 1 - np.array(filtered_y)/self.mec_power
            filtered_y = filtered_y.tolist()
            if len(filtered_y) == 0:
                print(filtered_y)
                raise Exception('After filtering, nothing is left')
            self.my_scatter_plot(filtered_x, filtered_y, filtered_filtered_O_list[::], 
                                    fig=fig, ax=ax, marker=marker,index_list=filtered_index_list)

        ax.set_xlabel('Stack length [mm]')
        ax.set_ylabel(r'Efficiency at %d kW [1]'%(self.mec_power*1e-3))
        ax.grid(True)

    def pareto_plot_torque_force(self, fig2, axeses, marker='.'):

        O = fobj_list(  list(self.get_certain_objective_function(2)),
                        list(self.get_certain_objective_function(4)),
                        list(self.get_certain_objective_function(3)),
                        list(self.get_certain_objective_function(5)),
                        list(self.get_certain_objective_function(6)),
                        list(self.get_certain_objective_function(15)),
                        weights=self.weights_used, rotor_volume=self.rotor_volume, rotor_weight=self.rotor_weight).tolist()

        # fig, axeses = subplots(2, 2, sharex=False, dpi=150, figsize=(12, 6), facecolor='w', edgecolor='k')
        # # fig, axeses = subplots(2, 2, sharex=False, dpi=150, figsize=(10, 8), facecolor='w', edgecolor='k')
        # fig.subplots_adjust(right=0.9, hspace=0.21, wspace=0.11) # won't work after I did something. just manual adjust!

        # Use FRW and TRV
        if True:
            # TRV vs Torque Ripple
            ax = axeses[0][0]
            x, y = np.array(list(self.get_certain_objective_function(2)))/self.rotor_volume/1e3, list(self.get_certain_objective_function(3))
            x = x.tolist()
            self.my_scatter_plot(x,y,O[::], fig=fig2, ax=ax, marker=marker)
            ax.set_xlabel('TRV [kNm/m^3]\n(a)')
            ax.set_ylabel(r'$T_{\rm rip}$ [100%]')

            # FRW vs Ea
            ax = axeses[0][1]
            x, y = np.array(list(self.get_certain_objective_function(4)))/self.rotor_weight, list(self.get_certain_objective_function(6))
            x = x.tolist()
            self.my_scatter_plot(x,y,O[::], fig=fig2, ax=ax, marker=marker)
            ax.set_xlabel('FRW [1]\n(b)')
            ax.set_ylabel(r'$E_a$ [deg]')

            # FRW vs Em
            ax = axeses[1][0]
            x, y = np.array(list(self.get_certain_objective_function(4)))/self.rotor_weight, list(self.get_certain_objective_function(5))
            x = x.tolist()
            self.my_scatter_plot(x,y,O[::], fig=fig2, ax=ax, marker=marker)
            ax.set_xlabel('FRW [1]\n(c)')
            ax.set_ylabel(r'$E_m$ [100%]')

            # Em vs Ea
            ax = axeses[1][1]
            x, y = list(self.get_certain_objective_function(5)), list(self.get_certain_objective_function(6))
            scatter_handle = self.my_scatter_plot(x,y,O[::], fig=fig2, ax=ax, marker=marker)
            ax.set_xlabel('$E_m$ [100%]\n(d)')
            ax.set_ylabel(r'$E_a$ [deg]')

            # fig.subplots_adjust(right=0.8)
            # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7]) # left, bottom, width, height
            # fig.subplots_adjust(right=0.9)
            # cbar_ax = fig.add_axes([0.925, 0.15, 0.025, 0.7])
            cbar_ax = fig2.add_axes([0.925, 0.15, 0.02, 0.7])
            cbar_ax.get_yaxis().labelpad = 10
            clb = fig2.colorbar(scatter_handle, cax=cbar_ax)
            clb.ax.set_ylabel('Cost function $%s$'%(self.weights_name), rotation=270)
            # clb.ax.set_title(r'Cost function $O_2$', rotation=0)

        # Use Torque and Force
        if False:

            # Torque vs Torque Ripple
            ax = axeses[0][0]
            xy_ref = (19.1197, 0.0864712) # from run#117
            x, y = list(self.get_certain_objective_function(2)), list(self.get_certain_objective_function(3))
            my_scatter_plot(x,y,O[::],xy_ref,O_ref, fig=fig2, ax=ax)
            ax.set_xlabel('$T_{em}$ [Nm]\n(a)')
            ax.set_ylabel(r'$T_{\rm rip}$ [100%]')

            # Force vs Ea
            ax = axeses[0][1]
            xy_ref = (96.9263, 6.53137)
            x, y = list(self.get_certain_objective_function(4)), list(self.get_certain_objective_function(6))
            my_scatter_plot(x,y,O[::],xy_ref,O_ref, fig=fig2, ax=ax)
            ax.set_xlabel('$|F|$ [N]\n(b)')
            ax.set_ylabel(r'$E_a$ [deg]')

            # Force vs Em
            ax = axeses[1][0]
            xy_ref = (96.9263, 0.104915)
            x, y = list(self.get_certain_objective_function(4)), list(self.get_certain_objective_function(5))
            my_scatter_plot(x,y,O[::],xy_ref,O_ref, fig=fig2, ax=ax)
            ax.set_xlabel('$|F|$ [N]\n(c)')
            ax.set_ylabel(r'$E_m$ [100%]')

            # Em vs Ea
            ax = axeses[1][1]
            xy_ref = (0.104915, 6.53137)
            x, y = list(self.get_certain_objective_function(5)), list(self.get_certain_objective_function(6))
            scatter_handle = my_scatter_plot(x,y,O[::],xy_ref,O_ref, fig=fig2, ax=ax)
            ax.set_xlabel('$E_m$ [100%]\n(d)')
            ax.set_ylabel(r'$E_a$ [deg]')

            # fig.subplots_adjust(right=0.8)
            # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7]) # left, bottom, width, height
            # fig.subplots_adjust(right=0.9)
            # cbar_ax = fig.add_axes([0.925, 0.15, 0.025, 0.7])
            cbar_ax = fig2.add_axes([0.925, 0.15, 0.02, 0.7])
            cbar_ax.get_yaxis().labelpad = 10
            clb = fig2.colorbar(scatter_handle, cax=cbar_ax)
            clb.ax.set_ylabel('Cost function $%s$'%(self.weights_name), rotation=270)
            # clb.ax.set_title(r'Cost function $O_2$', rotation=0)


            fig2.tight_layout()
            # fig.savefig(r'D:\OneDrive\[00]GetWorking\32 blimopti\p2019_ecce_bearingless_induction\images\pareto_plot.png', dpi=150, bbox_inches='tight')
            show()

            # Loss vs Ea
            xy_ref = ((1817.22+216.216+224.706), 6.53137)
            x, y = np.array(list(self.get_certain_objective_function(9))) + np.array(list(self.get_certain_objective_function(12))) + np.array(list(self.get_certain_objective_function(13))), list(self.get_certain_objective_function(6))
            x = x.tolist()
            my_scatter_plot(x,y,O[::],xy_ref,O_ref)
            xlabel(r'$P_{\rm Cu,Fe}$ [W]')
            ylabel(r'$E_a$ [deg]')

            quit()   

    def sensitivity_bar_charts(self):
        INDEX_TOTAL_LOSS = 15 + 4 # index of total loss in the machine_data list

        # ------------------------------------ Sensitivity Analysis Bar Chart Scripts
        number_of_variant = self.sw.fea_config_dict['local_sensitivity_analysis_number_of_variants'] + 1

        from pylab import subplots, mpl
        mpl.rcParams['legend.fontsize'] = 12
        # mpl.rcParams['legend.family'] = 'Times New Roman'
        mpl.rcParams['font.family'] = ['Times New Roman']
        font = {'family' : 'Times New Roman', #'serif',
                'color' : 'darkblue',
                'weight' : 'normal',
                'size' : 14,}
        textfont = {'family' : 'Times New Roman', #'serif',
                    'color' : 'darkblue',
                    'weight' : 'normal',
                    'size' : 11.5,}
        mpl.style.use('classic')

        fig, axeses = subplots(4, 2, sharex=True, dpi=150, figsize=(16*0.75, 8*0.75), facecolor='w', edgecolor='k')
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
        param_list = [
        r'$L_g$',# 'air_gap_length_delta',
        r'$w_{st}$',# 'stator_tooth_width_b_ds',
        r'$w_{rt}$',# 'rotor_tooth_width_b_dr',
        r'$\theta_{so}$',# 'Angle_StatorSlotOpen',
        r'$w_{ro}$',# 'Width_RotorSlotOpen ',
        r'$d_{so}$',# 'Width_StatorTeethHeadThickness',
        r'$d_{ro}$']# 'Length_HeadNeckRotorSlot']
        y_label_list = ['PF', r'$\eta$ [100%]', r'$T_{em} [N]$', r'$T_{rip}$ [100%]', r'$|F|$ [N]', r'$E_m$ [100%]', r'$E_a$ [deg]', 
                        r'$P_{Cu,s,JMAG}$', r'$P_{Cu,r,JMAG}$', r'$P_{Fe}$ [W]', r'$P_{eddy}$', r'$P_{hyst}$', r'$P_{Cu,s,FEMM}$', r'$P_{Cu,r,FEMM}$', 
                        r'Windage loss', r'Total loss']
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

    def build_basic_info(self):
        spec = self.spec
        sw = self.sw
        # Basic information
        self.ExcitationFreqSimulated = spec.ExcitationFreqSimulated
        self.speed_rpm         = self.ExcitationFreqSimulated * 60 / spec.p # rpm
        self.Omega             = self.speed_rpm / 60. * 2*np.pi
        self.mec_power         = spec.mec_power / spec.ExcitationFreq * self.ExcitationFreqSimulated
        self.required_torque   = self.mec_power / self.Omega # Nm
        self.Radius_OuterRotor = sw.im.Radius_OuterRotor
        self.stack_length      = sw.im.stack_length
        self.Qs                = sw.im.Qs
        self.Qr                = sw.im.Qr
        self.rotor_volume      = sw.im.get_rotor_volume()
        self.rotor_weight      = sw.im.get_rotor_weight()
        self.weights_name      = sw.fea_config_dict['use_weights']
        self.weights_used      = use_weights(which=sw.fea_config_dict['use_weights'])
        self.stack_length      = sw.im.stack_length
        self.stack_length_max  = spec.Stack_Length_Max
        # self.spec              = spec
        # self.sw                = sw
        print('-'*50 + '\nutility.py')
        print('Qs=%d, rotor_volume=%g'%(self.Qs, self.rotor_volume), 'm^3')
        print('Qr=%d, rotor_weight=%g'%(self.Qr, self.rotor_weight), 'N')
        # O1_weights = use_weights(which='O1') # [ 1, 0.1,   1, 0.1, 0.1,   0 ]
        # O2_weights = use_weights(which='O2') # [ 1, 1.0,   1, 1.0, 1.0,   0 ]


#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# Post-processing
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~

def build_sensitivity_bar_charts(spec, sw):
    run_integer = int(sw.fea_config_dict['run_folder'][4:-1])
    swda = SwarmDataAnalyzer(sw, spec, dir_run=sw.dir_run, run_integer=run_integer, bool_sensitivity_analysis=True)
    return swda.sensitivity_bar_charts()

def build_Pareto_plot(spec, sw):
    run_integer = int(sw.fea_config_dict['run_folder'][4:-1])
    swda = SwarmDataAnalyzer(sw, spec, dir_run=sw.dir_run, run_integer=run_integer, bool_sensitivity_analysis=False) 

    # Plotting the Pareto plot
    from pylab import plt, subplots#, show
    plt.rcParams["font.family"] = "Times New Roman"

    fig, ax = subplots(1, 1, sharex=False, dpi=150, figsize=(12, 6), facecolor='w', edgecolor='k')
    swda.pareto_plot_eta_vs_stack_legnth(fig, ax, marker='o', bool_filtered=True) # Prototype DPNV O1
    fig.tight_layout()

    fig2, axeses = subplots(2, 2, sharex=False, dpi=150, figsize=(12, 6), facecolor='w', edgecolor='k')
    fig2.subplots_adjust(right=0.9, hspace=0.21, wspace=0.11) # won't work after I did something. just manual adjust!
    swda.pareto_plot_torque_force(fig2, axeses, marker='o') # Prototype DPNV O1

    swda.my_population_distribution_plots(sw.de_config_dict)

        # add initial design to the plot
        # ax.annotate('Initial design', xytext=(xy_ref[0]*0.95, xy_ref[1]*1.0), xy=xy_ref, xycoords='data', arrowprops=dict(arrowstyle="->"))

    # if run_integer == 501:
    #     swda.__init__(run_integer=500, bool_sensitivity_analysis=True) # the sensitivity results
    #     swda.pareto_plot_eta_vs_stack_legnth(fig, ax, marker='^', bool_filtered=True) # Prototype DPNV O1

    # show()
    return swda # swda.best_design_denorm




# test for power factor (Goertzel Algorithm with periodic extension 1000), etc.
if __name__ == '__main__':
    swda = SwarmDataAnalyzer(run_integer=121)

    # Pseudo Pareto Optimal Front
    gen_best = swda.get_best_generation()
    with open('d:/gen#0000.txt', 'w') as f:
        f.write('\n'.join(','.join('%.16f'%(x) for x in y) for y in gen_best)) # convert 2d array to string            

    design_parameters_norm = (gen_best - min_b) / diff


    for el in design_parameters_norm:
        print(','.join('%.4f'%(_) for _ in el.tolist()))
    print('airgap length\n', [el[1] for el in gen_best])

    print('Average value of design parameters')
    avg_design_parameters = []
    for el in design_parameters_norm.T:
        avg_design_parameters.append(sum(el)/len(el))
    print(avg_design_parameters)
    avg_design_parameters_denorm = min_b + avg_design_parameters * diff
    print(avg_design_parameters_denorm)

    # for design in swda.get_best_generation(generator=swda.design_display_generator()):
    #     print ''.join(design),
    quit()

    cost = swda.list_cost_function()
    indices, items = min_indices(cost, 50)
    print(indices)
    print(items)

    stop = max(indices)
    start = min(indices)
    print(start, stop)

    gene = swda.design_parameters_generator()
    gen_best = []
    for index, design in enumerate(gene):
        if index in indices:
            gen_best.append(design)
    print('there are', index, 'designs')

    # print  len(list(gene))
    # for index in indices:
    #     start = index
    # print next(itertools.islice(gene, start, None)) 


    quit()
    # print min_indices([3,5,7,4,2], 5)
    # print max_indices([3,5,7,4,2], 5)
    # quit()



    print(''.join(swda.buf[:21]), end=' ')

    print(swda.find_individual(14, 0))


    for design in swda.design_display_generator():
        print(design, end=' ')
        break

    for design in swda.design_parameters_generator():
        print(design)
        break

    print() 

    # for generation in range(5):
    #     print '----------gen#%d'%(generation)
    #     generation_file_path = r'D:\OneDrive - UW-Madison\c\pop\run#107/' + 'gen#%04d.txt'%(generation)
    #     print generation_file_path
    #     if os.path.exists( generation_file_path ):
    #         with open(generation_file_path, 'r') as f:
    #             for el in f.readlines():
    #                 print el[:-1]

    # read voltage and current to see the power factor!
    # read voltage and current to see the power factor!
    # read voltage and current to see the power factor!
    
    # 绘制损耗图形。
    # 绘制损耗图形。
    # 绘制损耗图形。

    if False:
        from pylab import *
        gs_u = Goertzel_Data_Struct("Goertzel Struct for Voltage\n")
        gs_i = Goertzel_Data_Struct("Goertzel Struct for Current\n")

        phase = arccos(0.6) # PF=0.6
        targetFreq = 1000.
        TS = 1.5625e-5

        if False: # long signal
            time = arange(0./targetFreq, 100.0/targetFreq, TS)
            voltage = 3*sin(targetFreq*2*pi*time + 30/180.*pi)
            current = 2*sin(targetFreq*2*pi*time + 30/180.*pi - phase)
        else: # half period
            time = arange(0.5/targetFreq, 1.0/targetFreq, TS)
            voltage = 3*sin(targetFreq*2*pi*time + 30/180.*pi)
            current = 2*sin(targetFreq*2*pi*time + 30/180.*pi - phase)

        N_SAMPLE = len(voltage)
        noise = ( 2*rand(N_SAMPLE) - 1 ) * 0.233
        voltage += noise

        # test
        print('PF=', compute_power_factor_from_half_period(voltage, current, time, targetFreq=1e3)[0])
        quit()

        # plot(voltage+0.5)
        # plot(current+0.5)

        # Check for resolution
        end_time = N_SAMPLE * TS
        resolution = 1./end_time
        print(resolution, targetFreq)
        if resolution > targetFreq:

            if True: # for half period signal
                print('Data length (N_SAMPLE) too short or sampling frequency too high (1/TS too high).')
                print('Periodical extension is applied. This will not really increase your resolution. It is a visual trick for Fourier Analysis.')
                voltage = voltage.tolist() + (-voltage).tolist() #[::-1]
                current = current.tolist() + (-current).tolist() #[::-1]

            voltage *= 1000
            current *= 1000

            N_SAMPLE = len(voltage)
            end_time = N_SAMPLE * TS
            resolution = 1./end_time
            print(resolution, targetFreq, 'Hz')
            if resolution <= targetFreq:
                print('Now we are good.')


            # 目前就写了二分之一周期的处理，只有四分之一周期的数据，可以用反对称的方法周期延拓。

        print(gs_u.goertzel_offline(targetFreq, 1./TS, voltage))
        print(gs_i.goertzel_offline(targetFreq, 1./TS, current))

        gs_u.ampl = sqrt(gs_u.real*gs_u.real + gs_u.imag*gs_u.imag) 
        gs_u.phase = arctan2(gs_u.imag, gs_u.real)

        print('N_SAMPLE=', N_SAMPLE)
        print("\n")
        print(gs_u.id)
        print("RT:\tAmplitude: %g, %g\n" % (gs_u.ampl, sqrt(2.0*gs_u.accumSquaredData/N_SAMPLE)))
        print("RT:\tre=%g, im=%g\n\tampl=%g, phase=%g\n" % (gs_u.real, gs_u.imag, gs_u.ampl, gs_u.phase*180/math.pi))

        # // do the analysis
        gs_i.ampl = sqrt(gs_i.real*gs_i.real + gs_i.imag*gs_i.imag) 
        gs_i.phase = arctan2(gs_i.imag, gs_i.real)
        print(gs_i.id)
        print("RT:\tAmplitude: %g, %g\n" % (gs_i.ampl, sqrt(2.0*gs_i.accumSquaredData/N_SAMPLE)))
        print("RT:\tre=%g, im=%g\n\tampl=%g, phase=%g\n" % (gs_i.real, gs_i.imag, gs_i.ampl, gs_i.phase*180/math.pi))

        print("------------------------")
        print("\tPhase Difference of I (version 1): %g\n" % ((gs_i.phase-gs_u.phase)*180/math.pi), 'PF=', cos(gs_i.phase-gs_u.phase))

        # // Take reference to the voltage phaser
        Ireal = sqrt(2.0*gs_i.accumSquaredData/N_SAMPLE) * cos(gs_i.phase-gs_u.phase)
        Iimag = sqrt(2.0*gs_i.accumSquaredData/N_SAMPLE) * sin(gs_i.phase-gs_u.phase)
        print("\tAmplitude from accumSquaredData->%g\n" % (sqrt(2.0*gs_u.accumSquaredData/N_SAMPLE)))
        print("\tPhaser I:\tre=%g, im=%g\n" % (Ireal, Iimag))
        print("\tPhase Difference of I (version 2): %g\n" % (arctan2(Iimag,Ireal)*180/math.pi), 'PF=', cos(arctan2(Iimag,Ireal)))

        plot(voltage)
        plot(current)

        show()

        # it works!
        # send_notification('Test email')

# find best individual
if __name__ == '__main__':
    swda = SwarmDataAnalyzer(run_integer=142)
    gen_best = swda.get_best_generation(popsize=30)

    with open('d:/Qr16_gen_best.txt', 'w') as f:
        f.write('\n'.join(','.join('%.16f'%(x) for x in y) for y in gen_best)) # convert 2d array to string            
    quit()

