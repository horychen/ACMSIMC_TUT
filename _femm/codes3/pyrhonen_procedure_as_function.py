from pylab import *
import sys
import os
import logging

# General Information
# Steel
    # Yield Stress of M19 Steel:
        # 350 MPa \cite{lovelace2004mechanical}
        # % The typical yield stress of steel is 300 MPa \cite{pyrhonen2009design}.
    # Density:
        # Silicon steels around 2 to 3% silicon composition have a density about 7.65 gr/cc.  The cold rolling process is volume conserving.  Just currious why you care about the density. from https://www.quora.com/What-is-the-material-density-for-M19-silicon-steel

# print('-'*20+'\nPyrhonen 2009 Chapter 7.')
one_report_dir_prefix = '../release/OneReport/OneReport_TEX/contents/'
file_name = 'pyrhonen_procedure'
file_suffix = '.tex'

def get_parallel_tooth_height(area_rotor_slot_Sur, rotor_tooth_width_b_dr, Qr, rotor_outer_radius_r_or_eff):
    # rotor slot height depends on rotor_tooth_width_b_dr and rotor current (power factor)
    temp = (2*pi*rotor_outer_radius_r_or_eff - Qr*rotor_tooth_width_b_dr)
    rotor_Delta = temp**2 - 4*pi*area_rotor_slot_Sur*Qr
    rotor_tooth_height_h_dr_plus  = ( +sqrt(rotor_Delta) + temp ) / (2*pi)
    rotor_tooth_height_h_dr_minus = ( -sqrt(rotor_Delta) + temp ) / (2*pi)
    rotor_tooth_height_h_dr = rotor_tooth_height_h_dr_minus
    return rotor_tooth_height_h_dr, rotor_tooth_height_h_dr_plus, rotor_Delta

class geometry_data(object):
    def __init__(self):
        pass

from winding_layout import winding_layout
class desgin_specification(object):
    def __init__(self,
                    PS_or_SC = None,
                    DPNV_or_SEPA = None,
                    p  = None,
                    ps = None,
                    mec_power = None,
                    ExcitationFreq = None,
                    ExcitationFreqSimulated = None,
                    VoltageRating = None,
                    TangentialStress = None,
                    Qs = None,
                    Qr = None,
                    Js = None,
                    Jr = None,
                    Coil = None,
                    space_factor_kCu = None,
                    Conductor = None,
                    space_factor_kAl = None,
                    Temperature = None,
                    Steel = None,
                    lamination_stacking_factor_kFe = None,
                    stator_tooth_flux_density_B_ds = None,
                    rotor_tooth_flux_density_B_dr = None,
                    stator_yoke_flux_density_Bys = None,
                    rotor_yoke_flux_density_Byr = None,
                    guess_air_gap_flux_density = None,
                    guess_efficiency = None,
                    guess_power_factor = None,
                    safety_factor_to_yield = None,
                    safety_factor_to_critical_speed = None,
                    use_drop_shape_rotor_bar = None,
                    no_segmented_magnets = None,
                    tip_speed = None,
                    debug_or_release= None,
                    bool_skew_stator = None,
                    bool_skew_rotor = None,
                ):
        self.PS_or_SC = PS_or_SC
        self.DPNV_or_SEPA = DPNV_or_SEPA
        self.p = p
        self.ps = ps
        self.mec_power = mec_power
        self.ExcitationFreq = ExcitationFreq
        self.ExcitationFreqSimulated = ExcitationFreqSimulated
        self.VoltageRating = VoltageRating
        self.TangentialStress = TangentialStress
        self.Qs = Qs
        self.Qr = Qr
        self.Js = Js
        self.Jr = Jr
        self.Jr_backup = self.Jr
        self.Coil = Coil
        self.space_factor_kCu = space_factor_kCu
        self.Conductor = Conductor
        self.space_factor_kAl = space_factor_kAl
        self.Temperature = Temperature
        self.Steel = Steel
        self.lamination_stacking_factor_kFe = lamination_stacking_factor_kFe
        self.stator_tooth_flux_density_B_ds = stator_tooth_flux_density_B_ds
        self.rotor_tooth_flux_density_B_dr = rotor_tooth_flux_density_B_dr
        self.stator_yoke_flux_density_Bys = stator_yoke_flux_density_Bys
        self.rotor_yoke_flux_density_Byr = rotor_yoke_flux_density_Byr
        self.guess_air_gap_flux_density = guess_air_gap_flux_density
        self.guess_efficiency = guess_efficiency
        self.guess_power_factor = guess_power_factor
        self.safety_factor_to_yield = safety_factor_to_yield
        self.safety_factor_to_critical_speed = safety_factor_to_critical_speed
        self.use_drop_shape_rotor_bar = use_drop_shape_rotor_bar
        self.no_segmented_magnets = no_segmented_magnets
        self.tip_speed = tip_speed
        self.debug_or_release = debug_or_release
        self.bool_skew_stator = bool_skew_stator
        self.bool_skew_rotor  = bool_skew_rotor 

        self.winding_layout = winding_layout(self.DPNV_or_SEPA, self.Qs, self.p)

        self.bool_high_speed_design = self.tip_speed is not None

        self.geometry = geometry_data() # not used

        if not os.path.isdir('../' + 'pop/'):
            os.mkdir('../' + 'pop/')
        self.loc_txt_file = '../' + 'pop/' + r'initial_design.txt'
        open(self.loc_txt_file, 'w').close() # clean slate to begin with

    def build_name(self, bLatex=False):
        u'''转子类型+基本信息+槽配合+电负荷+导电材料+温度+导磁材料+磁负荷+叠压系数+假设+目的'''
        name = '%s%sp%dps%d_%dkW%dHz%dVtan%d_Qs%dQr%dJs%gJr%g%s%d%s%d@%dK_%s@%d_Bt%gs%grBy%gs%gr_b%.2fEff%gPf%g' % (
                'PS' if self.PS_or_SC else 'SC',
                'DPNV' if self.DPNV_or_SEPA else 'SEPA',
                self.p,
                self.ps,
                self.mec_power*1e-3,
                self.ExcitationFreq,
                self.VoltageRating,
                self.TangentialStress*1e-3,
                self.Qs,
                self.Qr,
                self.Js*1e-6,
                self.Jr*1e-6,
                self.Coil,
                self.space_factor_kCu*100,
                self.Conductor,
                self.space_factor_kAl*100,
                self.Temperature,
                self.Steel,
                self.lamination_stacking_factor_kFe*100,
                self.stator_tooth_flux_density_B_ds,
                self.rotor_tooth_flux_density_B_dr,
                self.stator_yoke_flux_density_Bys,
                self.rotor_yoke_flux_density_Byr,
                self.guess_air_gap_flux_density,
                self.guess_efficiency,
                self.guess_power_factor
                )
        if self.tip_speed is None:
            name_part2 = "SFy%gSFcs%g"%(
                    self.safety_factor_to_yield,
                    self.safety_factor_to_critical_speed,
                    )
        else:
            name_part2 = "tip%.0f"%(
                    self.tip_speed
                    )
        name_part3 = "_%s%s%s%s%s"%(
                    'drop'  if self.use_drop_shape_rotor_bar else 'rod',
                    'Pitch%d'%(self.winding_layout.coil_pitch),
                    'Sskew' if self.bool_skew_stator is not None else '',
                    'Rskew' if self.bool_skew_rotor  is not None else '',
                    'debug' if self.debug_or_release else 'release'
                    )

        name += name_part2 + name_part3

        if bLatex:
            return name.replace('_', '\\_')
        else:
            return name

    def show(self, toString=False):
        attrs = list(vars(self).items())
        key_list = [el[0] for el in attrs]
        val_list = [el[1] for el in attrs]
        the_dict = dict(list(zip(key_list, val_list)))
        sorted_key = sorted(key_list, key=lambda item: (int(item.partition(' ')[0]) if item[0].isdigit() else float('inf'), item)) # this is also useful for string beginning with digiterations '15 Steel'.
        tuple_list = [(key, the_dict[key]) for key in sorted_key]
        if toString==False:
            print('- Bearingless Induction Motor Design Specs\n\t', end=' ')
            print(', \n\t'.join("%s = %s" % item for item in tuple_list))
            return ''
        else:
            return '\n- Bearingless Induction Motor Design SPecs\n\t' + ', \n\t'.join("%s = %s" % item for item in tuple_list)

    # @blockPrinting
    def pyrhonen_procedure(self, pc_name, bool_loop_Jr=True):
        # for self.TangentialStress in arange(12000, 33001, 3000): # self.TangentialStress - this value will change the motor size, so it is not suited for loop for optimization bounds

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 1. Initial Design Parameters
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        no_phase_m = 3
        speed_rpm = self.ExcitationFreq * 60 / self.p # rpm
        U1_rms = self.VoltageRating / sqrt(3) # V - Wye-connect #480 V is standarnd # 电压越高，意味着越厚的绝缘占去槽空间（Lipo2017书）
            # U1_rms = 500. / sqrt(3) # The design used in ECCE
        stator_phase_voltage_rms = U1_rms
        bool_we_have_plenty_voltage = True
        print('bool_we_have_plenty_voltage is %s' % 'True' if bool_we_have_plenty_voltage else 'False')
        if math.isclose(480, self.VoltageRating):
            bool_standard_voltage_rating = True
        else:
            bool_standard_voltage_rating = False

        fname = open(one_report_dir_prefix+file_name+'_s01'+file_suffix, 'w', encoding='utf-8') if 'Y730' in pc_name else None
        print(r'''\subsubsection{Initial Design Parameters}''', file=fname)
        print('\nDesign Name:\\\\{\\tiny %s}' % self.build_name(bLatex=True), file=fname)
        if self.PS_or_SC:
            print('\nRotor type: pole-specific rotor', file=fname)
        else:
            print('\nRotor type: squirrel cage rotor', file=fname)
        print('\nTorque winding poles: $2p=%d$' %(2*self.p), file=fname)
        print('\nSuspension winding poles: $2p_s=%d$' %(2*self.ps), file=fname)
        print('\nPower: $P_{mec}=%g$ kW' %(self.mec_power*1e-3), file=fname)
        print('\nMaximum frequency: %g Hz' % (self.ExcitationFreq), file=fname)
        print('\nRated frequency (for FEA study): %g Hz' % (self.ExcitationFreqSimulated), file=fname)
        print('\nVoltage rating: %g Vrms (line-to-line, Wye-Connect)' % self.VoltageRating, file=fname)
        print('\nTangential stress: $\\sigma_{\\rm tan}=%g$ Pa' % self.TangentialStress, file=fname)
        print('\nRotor speed: $n_{mec}=%g$ r/min' % speed_rpm, file=fname)
        print('\nNumber of stator slot: $Q_s=%d$'%(self.Qs), file=fname)
        print('\nNumber of rotor slot: $Q_r=%d$'%(self.Qr), file=fname)
        print('\nStator current density: $J_s=%g$ Arms/$\\rm m^2$'%(self.Js), file=fname)
        print('\nRotor current density: $J_r=%g$ Arms/$\\rm m^2$'%(self.Jr), file=fname)
        print('\nMagnetic material: %s with a stack factor of %g\\%%'  %(self.Steel, 100*self.lamination_stacking_factor_kFe), file=fname)
        print('\nCoil material: %s with a fill/packing factor %g\\%%'      %(self.Coil, 100*self.space_factor_kCu), file=fname)
        print('\nConductor material: %s with a fill factor %g\\%%' %(self.Conductor, 100*self.space_factor_kAl), file=fname)
        print('\nCoil/conductor temperature: %g deg Celcius' %(self.Temperature), file=fname)
        print('\nStator tooth flux density: $B_{ds}=%g$ T'    % (self.stator_tooth_flux_density_B_ds), file=fname)
        print('\nRotor tooth flux density: $B_{dr}=%g$ T'     % (self.rotor_tooth_flux_density_B_dr), file=fname)
        print('\nStator yoke flux density: $B_{ys}=%g$ T'     % (self.stator_yoke_flux_density_Bys), file=fname)
        print('\nRotor yoke flux density: $B_{yr}=%g$ T'      % (self.rotor_yoke_flux_density_Byr), file=fname)
        print('\nGuess air gap flux density: $B_\\delta=%g$ T'% (self.guess_air_gap_flux_density), file=fname)
        print('\nGuess efficiency: $\\eta=%g$\\%%'% (100*self.guess_efficiency), file=fname)
        print('\nGuess efficiency: ${\\rm PF}=%g$'% (self.guess_efficiency), file=fname)
        print('\nDebug or release: %s'% (self.debug_or_release), file=fname)
        print('\n*Note: The default option for voltage availability is plenty. That is, the dc bus can be high enough (higher than the specified votlage rating) for allowing more conductors in stator slots.', file=fname)



    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 2. Machine Constants
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # if self.PS_or_SC == False:
            #     self.TangentialStress = 21500 # 12000 ~ 33000
            # else: # for wound rotor, the required rotor slot is really too large to find a solution, we must reduce self.TangentialStress to allow a larger rotor.
            #     self.TangentialStress = 12000 # larger rotor to increase rotor slot area and to reduce rotor heating.
        if self.p == 1:
            machine_constant_Cmec = 150 # kW s / m^3. Figure 6.3 <- This is what Eric calls the penalty on the 2 pole IM---you loss power density.
        else:
            machine_constant_Cmec = 250 # kW s / m^3. 这里的电机常数都是按100kW去查的表哦，偏大了。不过电机常数在这里只用于检查AJ值的合理性。
        fname = open(one_report_dir_prefix+file_name+'_s02'+file_suffix, 'w', encoding='utf-8') if 'Y730' in pc_name else None
        print(r'''\subsubsection{Machines Constant and Tangential Stress}
                Tangential stress.
                \[{\sigma _{F\tan }} = \frac{{\hat A{{\hat B}_\delta }\cos \varphi }}{2} = \frac{{A{{\hat B}_\delta }\cos \varphi }}{{\sqrt 2 }} \in [12000-21500] {\rm Pa}\]

                Machine constant (mechanical version).
                \[{C_{{\rm{mec}}}} = \frac{{{\pi ^2}}}{{\sqrt 2 }}{k_{ws1}}A{{\hat B}_\delta } \Rightarrow {V_r} = \frac{\pi }{4}\frac{{{S_i}}}{{C{n_{syn}}}} ~(\text{Not used.})\]

                The former is used to determine motor size, and the latter is for determining linear current density $A$ according to the air gap $B_\delta$. We choose $\sigma_{F\tan}$ to be 12000 Pa to reduce current density on the rotor bars.
                ''', file=fname)
        print('\nTangential stress: $\\sigma_{F{\\rm tan}}=%g$ Pa' % self.TangentialStress, file=fname)
        print('\nMachine constant: $C_{\\rm mec}=%g$ kWs/$\\rm m^3$'% machine_constant_Cmec, file=fname)



    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 3. Machine Sizing
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        required_torque = self.mec_power/(2*pi*speed_rpm)*60
        rotor_volume_Vr = required_torque/(2*self.TangentialStress)
        self.required_torque = required_torque    

        if self.bool_high_speed_design == True:
            length_ratio_chi = None
            rotor_outer_radius_r_or = eric_specify_tip_speed_get_radius(self.tip_speed, speed_rpm)
            rotor_outer_diameter_Dr = rotor_outer_radius_r_or*2
            stack_length = rotor_volume_Vr / (pi * rotor_outer_radius_r_or**2)
        else:
            #     Usually, the ratio of the length of the machine to the air-gap diameter χ = l/D is
            # selected for operation of the rotor below the first critical rotation speed. --P296
            length_ratio_chi = pi/(2*self.p) * self.p**(1/3.) # Table 6.5

            # Considering only electromagnetic performances
            rotor_outer_diameter_Dr = (4/pi*rotor_volume_Vr*length_ratio_chi)**(1/3.)
            rotor_outer_radius_r_or = 0.5 * rotor_outer_diameter_Dr
            stack_length = rotor_outer_diameter_Dr * length_ratio_chi
            print('rotor_outer_radius_r_or is', rotor_outer_radius_r_or)

            # Considering machanical loading
            rotor_radius_max = get_outer_rotor_radius_yield(speed_rpm, safety_factor_to_yield=self.safety_factor_to_yield)
            if rotor_outer_radius_r_or>rotor_radius_max:
                rotor_outer_radius_r_or = rotor_radius_max
                rotor_outer_diameter_Dr = 2 * rotor_outer_radius_r_or
                stack_length = rotor_volume_Vr / (pi * rotor_outer_radius_r_or**2)
                bool_mechanical_dominate = True
            print('rotor_outer_radius_r_or is', rotor_outer_radius_r_or, 'under mechanical limit')

        stack_length_max = get_stack_length_critical_speed(speed_rpm, rotor_outer_radius_r_or, safety_factor_to_critical_speed=self.safety_factor_to_critical_speed)
        self.Stack_Length_Max = stack_length_max*1e3
        if stack_length>stack_length_max:
            raise Exception('For current safety factor and power rating, no design can be sought under first critical speed.')
        
        # rotor_volume_Vr = pi * rotor_outer_radius_r_or**2 * stack_length, so
        # rotor_outer_radius_r_or = 28.8e-3 # based on the mechanical check with a safety factor of 1.5 according to the ECCE paper of Yegu Kang
        # rotor_outer_diameter_Dr = 2*rotor_outer_radius_r_or
        # # We have: rotor_volume_Vr = pi * rotor_outer_radius_r_or**2 * stack_length, so

        fname = open(one_report_dir_prefix+file_name+'_s03'+file_suffix, 'w', encoding='utf-8') if 'Y730' in pc_name else None
        print(r'''\subsubsection{Machine Sizing}\label{subsubsec:machine_sizing}
                Required torque.
                \[{T_{em}} = \frac{{{P_{mec}}}}{{2\pi {n_{mec}} \times 60}}\]
                \[{T_{em}} = {\rm Radius} \times {\rm Tangential~Force}  = {r_r}{\sigma _{F\tan }}{S_r} = {r_r}{\sigma _{F\tan }}\left( {2\pi {r_r}l'} \right)\]

                Required volume of rotor.
                \[{V_r} = \frac{{{T_{em}}}}{{2{\sigma _{\tan }}}},{V_r} = \frac{\pi }{4}D_r^2l'\]

                Ratio between stack length and diameter.
                \[\chi  = l'/D \Rightarrow l' = \chi D \Rightarrow \frac{\pi }{4}D_r^2\chi D = {V_r},\]
                \[let\,{D_r} \approx D \Rightarrow {D_r} = {\left( {\frac{4}{\pi }\frac{{{V_r}}}{\chi }} \right)^{\frac{1}{3}}}\]

                Limit by centrifugal force.
                \[{r_{r,\max }} = \sqrt {\frac{1}{{{k_\sigma }}}\frac{{{\sigma _{yield}}}}{{C'\rho {\Omega ^2}}}} \]

                Limit by natural frequency.
                \[l_{\max }^2 = {n^2}\frac{{{\pi ^2}}}{{k\Omega }}\sqrt {\frac{{EI}}{{\rho S}}} \]

                 ''', file=fname)
        print('\nRequired Torque: %g Nm'% required_torque, file=fname)

        if self.bool_high_speed_design == True:
            print('\nTip speed: %g m/s (specified)' % (self.tip_speed), file=fname)
        else:
            tip_speed = get_tip_speed(speed_rpm, rotor_outer_radius_r_or)
            print('\nTip speed: %g m/s with a safety factor to yield of %g' %(tip_speed, self.safety_factor_to_yield), file=fname)

            self.tip_speed = tip_speed

            # print 'Yegu Kang: rotor_outer_diameter_Dr, 95 mm'

            if bool_mechanical_dominate:
                print('\nRadially restricted by the machanical loading, and the maximal rotor radius is %g mm ($k_\\sigma=%g$)' % (rotor_radius_max*1e3, self.safety_factor_to_yield), file=fname)
            else:
                print('\nUse recommended length ratio: $\\chi=%g$' % length_ratio_chi, file=fname)

        print('\nCentrifugal stress: %g MPa' %(1e-6*check_stress_due_to_centrifugal_force_Pyrhonen14(speed_rpm, rotor_outer_radius_r_or)), file=fname)
        print('\nRotor outer diameter $D_{or}=%g$ mm'% (rotor_outer_diameter_Dr*1e3), file=fname)
        print('\nRotor outer radius $r_{or}=%g$ mm'% (rotor_outer_radius_r_or*1e3), file=fname)

        print('\nStack length: $L_{st}=%g$ mm'% (stack_length*1e3), file=fname)
        print('\nStack length (max): $L_{st,\\max}=%g$ mm ($k=%g$)'% (stack_length_max*1e3, self.safety_factor_to_critical_speed), file=fname)



    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 4. Air Gap Length
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        if self.p == 1:
            air_gap_length_delta = (0.2 + 0.01*self.mec_power**0.4) / 1000
        else:
            air_gap_length_delta = (0.18 + 0.006*self.mec_power**0.4) / 1000
        if self.p == 2:
            if self.mec_power <= 75e3:
                air_gap_length_delta *= 2 # *=3 will not converge
            else:
                air_gap_length_delta *= 1.5
        elif self.p == 1:
            if self.mec_power <= 75e3:
                air_gap_length_delta *= 1.5
            else:
                air_gap_length_delta *= 1.25

        if self.tip_speed>100:
            air_gap_length_delta_high_speed = 0.001 + (rotor_outer_diameter_Dr / 0.07 + self.tip_speed/400) * 1e-3 # 第二版(6.25)
            print('High-Speed motor detected. The air_gap_length_delta_high_speed =', air_gap_length_delta_high_speed)

        stack_length_eff = stack_length + 2 * air_gap_length_delta

        air_gap_diameter_D = 1*air_gap_length_delta + rotor_outer_diameter_Dr
        stator_inner_diameter_Dis = 2*air_gap_length_delta + rotor_outer_diameter_Dr
        stator_inner_radius_r_is = 0.5*stator_inner_diameter_Dis

        fname = open(one_report_dir_prefix+file_name+'_s04'+file_suffix, 'w', encoding='utf-8') if 'Y730' in pc_name else None
        print(r'''\subsubsection{Air Gap}
                For 50 Hz machines, we have
                \[\delta  = \frac{{0.2 + 0.01P_{mec}^{0.4}}}{{1000}},\,p = 1\]
                \[\delta  = \frac{{0.18 + 0.006P_{mec}^{0.4}}}{{1000}},\,p > 1\]
                We double $\delta$ for high speed machines, as suggested by Pyrhonen---increase air gap length by 50\%--100\%.
                % (看integrated box 硕士论文 'Kevin S. Campbell: this is too small. 3.5 mm is good! 3.1.3，Pyrhonen09给的只适用50Hz电机，其实pyrhonen自己也有提到气隙要加大100%哦) ''', file=fname)
        print('\nAir gap length: $\\delta=%g$ mm' % (air_gap_length_delta*1e3), file=fname)
        if self.tip_speed>100:
            print('\nAir gap length (high speed): $\\delta_{hs} = %g$ mm' % (air_gap_length_delta_high_speed*1e3), file=fname)
        print('\nEffective stack length: $L_{st,{\\rm eff}}=%g$ mm'% (stack_length_eff*1e3), file=fname)



    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 5. Stator Winding and Slots
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        '''为了方便散热，应该加大槽数！
            为了获得正弦磁势，应该加大槽数！
            缺点，就是基波绕组系数降低。A high number of slots per pole and phase lead to a situation in which the influence of the fifth and seventh harmonics is almost insignificant. Py09

            The magnitude of the harmonic torques depends greatly on the ratio of the slot numbers of
            the stator and the rotor. The torques can be reduced by skewing the rotor slots with reselft to
            the stator. In that case, the slot torques resulting from the slotting are damped. In the design
            of an induction motor, selfial attention has to be paid to the elimination of harmonic torques. Py09

            注意，其实50kW、30krpm的电机和5kW、3krpm的电机产生的转矩是一样大的，区别就在产生的电压，如果限制电压，那就只能降低匝数（线圈变粗，电阻变小），方便提高额定电流！
            把代码改改，弄成双层绕组的（短距系数）。

            按照公式计算一下短距、分布、斜槽系数啊！

            Lahne 2015【这一篇文章顶十篇】Comparison of State-of-the-Art High-Speed High-Power Machines.pdf cite Lipo2017
            Why Qs=24 is better than Qs=12? 
            Answer: Because it gives higher breakdown torque as the leakage inductance is smaller.
            '''
        print('Regular available Qs choice is: k * 2*p * no_phase_m. However, we have to make sure it is integral slot for the bearing winding that has two pole pairs')
        for i in range(1, 5):
            if self.p == 1:
                print(i * 2*(self.p+1)*no_phase_m, end=' ')
            else:
                print(i * 2*(self.p)*no_phase_m, end=' ')

        distribution_q = self.Qs / (2*self.p*no_phase_m)
        pole_pitch_tau_p = pi*air_gap_diameter_D/(2*self.p) # (3/17)

        if self.winding_layout.no_winding_layer == 1:
            # full pitch - easy
            coil_span_W = pole_pitch_tau_p
        else: 
            # short pitch (not tested)
            stator_slot_pitch_tau_us = pi * air_gap_diameter_D / self.Qs
            coil_span_W = self.winding_layout.coil_pitch * stator_slot_pitch_tau_us
                # for 2 pole motor, the recommended short pitch is 0.7. --p76

        kd1 = 2*sin(1/no_phase_m*pi/2)/(self.Qs/(no_phase_m*self.p)*sin(1*pi*self.p/self.Qs))
        kq1 = sin(1*coil_span_W/pole_pitch_tau_p*pi/2)        
        ksq1 = 1
        if self.bool_skew_rotor:
            # (not tested)
            ksq1 = sin(pi/(2*no_phase_m*distribution_q)) / (pi/(2*no_phase_m*distribution_q)) # (4.77), skew_s = slot_pitch_tau_u
        if self.bool_skew_stator:
            raise Exception('Not implemented error.')
        kw1 = kd1 * kq1 * ksq1

        fname = open(one_report_dir_prefix+file_name+'_s05'+file_suffix, 'w', encoding='utf-8') if 'Y730' in pc_name else None
        print(r'''\subsubsection{Stator Winding and Slots}
                \[{k_{wv}} = {k_{pv}}{k_{dv}}{k_{sqv}} = \sin \left( {v\frac{W}{{{\tau _p}}}\frac{\pi }{2}} \right)\frac{{2\sin \left( {\frac{v}{m}\frac{\pi }{2}} \right)}}{{\frac{Q}{{mp}}\sin \left( {v\pi \frac{p}{Q}} \right)}}\frac{{\sin \left( {v\frac{s}{{{\tau _p}}}\frac{\pi }{2}} \right)}}{{v\frac{s}{{{\tau _p}}}\frac{\pi }{2}}}\]
                 ''', file=fname)

        if self.winding_layout.no_winding_layer == 1:
            # full pitch - easy
            print('\nSingle layer torque winding')
        else: 
            # short pitch (not tested)
            print('\nDouble layer torque winding')

        print('\nSlot per phase per pole: $q=\\frac{%d}{2\\times%d\\times%d}=%g$' % (self.Qs, self.p, no_phase_m, distribution_q), file=fname)
        print('\nCoil pitch in count is: %d'% self.winding_layout.coil_pitch, file=fname)
        print('\nCoil span/pitch in fraction: $\\frac{W}{\\tau_p}=%g$' % (coil_span_W/pole_pitch_tau_p), file=fname)
        print('\nPole pitch $\\tau_p =%g$ m'% pole_pitch_tau_p, file=fname)
        print('\nWinding factor: $k_{w1}=k_{d1}k_{q1}k_{sq1} = %g\\times%g\\times%g = %g$' % (kd1, kq1, ksq1, kw1), file=fname)
            # print('$pole_pitch_tau_p / (pi*air_gap_diameter_D)=%g$'% (pole_pitch_tau_p / (pi*air_gap_diameter_D)), file=fname)
            # self.Qr = 30 # Table 7.5. If self.Qs=24, then self.Qr!=28, see (7.113).



    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 6. Air gap flux density
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        self.guess_air_gap_flux_density
        linear_current_density_A = machine_constant_Cmec / (pi**2/sqrt(2)*kw1*self.guess_air_gap_flux_density)
        
        fname = open(one_report_dir_prefix+file_name+'_s06'+file_suffix, 'w', encoding='utf-8') if 'Y730' in pc_name else None
        print(r'''\subsubsection{Air Gap Flux Density}
                    Recall machine constant.
                    \[{C_{{\rm{mec}}}} = \frac{{{\pi ^2}}}{{\sqrt 2 }}{k_{ws1}}A{{\hat B}_\delta } = 150{\rm{kWs/}}{{\rm{m}}^{\rm{3}}}\]

                    Flux density in the air gap can be inferred from machine constant.
                    \[ \Rightarrow {{\hat B}_\delta } = \frac{{{C_{{\rm{mec}}}}}}{{\frac{{{\pi ^2}}}{{\sqrt 2 }}{k_{ws1}}A}}\]
                    where linear current density @(6.31) is 
                    \[A = \frac{{{J_{rms}}{S_{Cu,u}}}}{{{\tau _u}}}\]
                    Reversely, knowing air gap flux density gives you the linear current density.
                 ''', file=fname)

        print('\n[Guess] Air gap flux density (amplitude): $\\hat B=%g$ T' % self.guess_air_gap_flux_density, file=fname)
        if linear_current_density_A<65 and linear_current_density_A>30: # Example 6.4
            print('\nLinear current density: $A=%g$ kA/m'% linear_current_density_A, file=fname)
        else:
            print('\n[Warning] Bad linear current density: $A=%g$ kA/m' % linear_current_density_A, file=fname)



    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 7. Number of Coil Turns
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        desired_emf_Em = 0.95 * stator_phase_voltage_rms # 0.96~0.98, high speed motor has higher leakage reactance hence 0.95
        flux_linkage_Psi_m = sqrt(2)*desired_emf_Em / (2*pi*self.ExcitationFreq)

        alpha_i = 2/pi # ideal sinusoidal flux density distribusion, when the saturation happens in teeth, alpha_i becomes higher.
        air_gap_flux_Phi_m = alpha_i * self.guess_air_gap_flux_density * pole_pitch_tau_p * stack_length_eff
        no_series_coil_turns_N = sqrt(2)*desired_emf_Em / (2*pi*self.ExcitationFreq * kw1 * air_gap_flux_Phi_m) # p306

        fname = open(one_report_dir_prefix+file_name+'_s07'+file_suffix, 'w', encoding='utf-8') if 'Y730' in pc_name else None
        print(r'''\subsubsection{Number of Coil Turns}
                \[{\alpha _i} > \pi /2\]

                Turns per slot.
                \[{z_Q} = \frac{{{aN_s}}}{{pq}} = \frac{{{aN_s}}}{{p\frac{{{Q_s}}}{{2pm}}}} = \frac{{2{aN_s}m}}{{{Q_s}}}\]

                Total turns in series $N_s$ is derived from desired EMF.
                \[N = \frac{{\sqrt 2 {E_m}}}{{\omega {k_{w1}}{{\hat \Phi }_m}}} = \frac{{\sqrt 2 {E_m}}}{{\omega {k_{w1}}{\alpha _i}{{\hat B}_\delta }{\tau _p}l'}} \Rightarrow {{\hat B}_\delta } = \frac{{\sqrt 2 {E_m}}}{{\omega {k_{w1}}{\alpha _i}N{\tau _p}l'}}\]
                ''', file=fname)
        print('\n[EMF] Air gap flux linkage: $\\Psi_m=%g$ Wb'% (flux_linkage_Psi_m), file=fname)
        print('\n[Guess] Air gap flux: $\\Phi_m=%g$ Wb'% (air_gap_flux_Phi_m), file=fname)
        print('\n$\\Psi_m/\\Phi_m=%g$'% (flux_linkage_Psi_m/air_gap_flux_Phi_m), file=fname)
        print('\n$\\Psi_m/\\Phi_m/k_{w1}=%g$'% (flux_linkage_Psi_m/air_gap_flux_Phi_m/kw1), file=fname)
        print('\nPhase stator voltage: $U_1=%g$ Vrms'% (stator_phase_voltage_rms), file=fname)
        print('\nEMF to terminal voltage coefficient: %g' % (0.95), file=fname)
        print('\nDesired EMF $E_m=%g$ Vrms'% desired_emf_Em, file=fname)
        print('\nNumber of series coil turns (EMF): $N=%g$'% no_series_coil_turns_N, file=fname)

        no_series_coil_turns_N = round(no_series_coil_turns_N)
        print('\nNumber of series coil turns (round-up): $N=%g$'% no_series_coil_turns_N, file=fname)
        print('\nDivide it by: $pq=%d$' % (self.p*distribution_q), file=fname)
        print('\nThe remainder is: %g' % ( (no_series_coil_turns_N) % (self.p*distribution_q) ), file=fname)



    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 8. Round Up Number of Coil Turns
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        backup = no_series_coil_turns_N
        if bool_we_have_plenty_voltage:
            no_series_coil_turns_N = min([self.p*distribution_q*i for i in range(100,0,-1)], key=lambda x:abs(x - no_series_coil_turns_N)) # using larger turns value has priority
        else:
            no_series_coil_turns_N = min([self.p*distribution_q*i for i in range(100)], key=lambda x:abs(x - no_series_coil_turns_N))  # using lower turns value has priority # https://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value

        fname = open(one_report_dir_prefix+file_name+'_s08'+file_suffix, 'w', encoding='utf-8') if 'Y730' in pc_name else None
        print(r'''\subsubsection{Round up Number of Coil Turns}
                    This $N$ must be multiples of $pq$ and it affects air gap B.
                    For example of $pq=4$, you get an $N$ value as 18,
                    if you choose $N=20$, then the air gap B will reduce.
                    and if you choose $N$ value as 16, then the air gap B will increase.
                    %air_gap_flux_density_B=', air_gap_flux_density_B, 'T', '变得比0.8T大了，是因为你减少了匝数取整，反之亦然。'
                     ''', file=fname) # 2nd ed.@(6.12)
        print('\nSince voltage is%s plenty, the suggested N is %g---modify this manually, if necessary.' \
            % (' NOT' if not bool_we_have_plenty_voltage else '', no_series_coil_turns_N), file=fname)
        # if bool_standard_voltage_rating == True:
        #     # reduce the magnetic load on rotor tooth by allowing the back EMF to become larger.
        #     if no_series_coil_turns_N < backup:
        #         no_series_coil_turns_N += self.p*distribution_q
        #         print('\bAuto use a larger no_series_coil_turns_N=%d,\n\tto make sure the air gap B is lower than 0.8T so as to limit the rotor tooth width and the solution of a rotor slot height exist.' % (no_series_coil_turns_N), file=fname)
        #         print('\b\tThis means that you sacrifice your power factor (directly related to air gap B only) to make sure you have enough rotor slot area for your rotor current.', file=fname)
        #     # quit()
        # else:
        #     if no_series_coil_turns_N < backup:
        #         print('\nSince your voltage is not fixed, you can vary it to make the selection of no_series_coil_turns_N works without increasing your air gap B.', file=fname)
        #         quit()
        print('\nN=%d means desired EMF $E_m=%g$ Vrms becomes %g Vrms.'                              % (no_series_coil_turns_N, desired_emf_Em,         no_series_coil_turns_N * (2*pi*self.ExcitationFreq * kw1 * air_gap_flux_Phi_m) / sqrt(2)), file=fname)
        print('\nN=%d means desired line-to-line EMF $E_m\\times\\sqrt{3}=%g$ Vrms becomes %g Vrms.' % (no_series_coil_turns_N, desired_emf_Em*sqrt(3), no_series_coil_turns_N * (2*pi*self.ExcitationFreq * kw1 * air_gap_flux_Phi_m) / sqrt(2)*sqrt(3)), file=fname)
        print('\nEMF induced by one turn is: %g Vrms'% ((2*pi*self.ExcitationFreq * kw1 * air_gap_flux_Phi_m) / sqrt(2)), file=fname)

        if self.DPNV_or_SEPA:
            number_parallel_branch = 2
        else:
            number_parallel_branch = 1
        # print('\n''In some cases, eselfially in low-voltage, high-power machines, there may be a need to change the stator slot number, the number of parallel paths or even the main dimensions of the machine in order to find the appropriate number of conductors in a slot.''')
        no_conductors_per_slot_zQ = 2* no_phase_m * no_series_coil_turns_N /self.Qs * number_parallel_branch # (7.8)

        print('\nNumber of parallel branch: $a=%d$' % number_parallel_branch, file=fname)
        print('\nNumber of conductors per slot: $z_Q=%g$'%no_conductors_per_slot_zQ, file=fname)
        # print('\nThe actual number of series coil turns $N_a$ is %d' % (no_series_coil_turns_N/number_parallel_branch), file=fname) # this is wrong concept.
        print('\nAir gap length: $\\delta=%g$ mm' % (air_gap_length_delta*1e3), file=fname)
        print('\n[Guess] Air gap flux density: $B_\\delta=%g$ T' % (self.guess_air_gap_flux_density), file=fname)

        if no_conductors_per_slot_zQ % 2 != 0:
            raise Exception('This zQ does not suit for two layer winding.')

        loop_count = 0
        global BH_lookup, bdata, hdata
        while True:
            if loop_count>25:
                raise Exception("Abort. This script may not converge anyway.")
            loop_count += 1



    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 9. (Loop ON) Recalculate Air Gap Flux Density
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        #
            self.guess_air_gap_flux_density = (sqrt(2)*desired_emf_Em) / (2*pi*self.ExcitationFreq * kw1 *  alpha_i * no_series_coil_turns_N * pole_pitch_tau_p * stack_length_eff) # p306
            fname = open(one_report_dir_prefix+file_name+'_s09'+file_suffix, 'w', encoding='utf-8') if 'Y730' in pc_name else None
            print(r'''\subsubsection{Re-calculate Air Gap Flux Density and Beginning of Loop for $\alpha_i$}
                    \[{{\hat \Phi }_m} = {\alpha _i}{{\hat B}_\delta }{\tau _p}l'\]
                    \[{{\hat B}_\delta } = \frac{{\sqrt 2 {E_m}}}{{\omega {k_{w1}}{\alpha _i}N{\tau _p}l'}}\]
                        ''', file=fname) # 反电势 sqrt(2)*desired_emf_Em 用幅值哦！
            print('\nField sinusoidal coefficient (converged): $\\alpha_i=%g>\\frac{2}{\\pi}=0.6366$'%(alpha_i), file=fname)
            print('\nNew air gap flux density: $\\hat B_\\delta=%g$ T' % self.guess_air_gap_flux_density, file=fname) # '变得比0.8T大了，是因为你减少了匝数取整，反之亦然。'

            # re-compute other B_delta related values
            if True:
                air_gap_flux_Phi_m = alpha_i * self.guess_air_gap_flux_density * pole_pitch_tau_p * stack_length_eff
                no_series_coil_turns_N = sqrt(2)*desired_emf_Em / (2*pi*self.ExcitationFreq * kw1 * air_gap_flux_Phi_m) # p306 # 2ndEdition p316 (6.32)
                print('\nNew air gap flux: $\\Phi_m=%g$ Wb'% (air_gap_flux_Phi_m), file=fname)

                print('\nPhase stator voltage: $U_1=%g$ Vrms'% (stator_phase_voltage_rms), file=fname)
                print('\nEMF to terminal voltage coefficient: %g' % (0.95), file=fname)
                print('\nDesired EMF $E_m=%g$ Vrms'% desired_emf_Em, file=fname)
                print('\nNumber of series coil turns (EMF): $N=%g$'% no_series_coil_turns_N, file=fname)

                print('\nNumber of series coil turns (round-up): $N=%g$'% no_series_coil_turns_N, file=fname)
                print('\nDivide it by: $pq=%d$' % (self.p*distribution_q), file=fname)
                print('\nThe remainder is: %g' % ( (no_series_coil_turns_N) % (self.p*distribution_q) ), file=fname)
        
                linear_current_density_A = machine_constant_Cmec / (pi**2/sqrt(2)*kw1*self.guess_air_gap_flux_density)
                if linear_current_density_A<65 and linear_current_density_A>30: # Example 6.4
                    print('\nLinear current density: $A=%g$ kA/m'% linear_current_density_A, file=fname)
                else:
                    print('\n[Warning] Bad linear current density: $A=%g$ kA/m' % linear_current_density_A, file=fname)

                print('\nN=%g means desired EMF $E_m=%g$ Vrms becomes %g Vrms.'                              % ((no_series_coil_turns_N), desired_emf_Em,         no_series_coil_turns_N * (2*pi*self.ExcitationFreq * kw1 * air_gap_flux_Phi_m) / sqrt(2)), file=fname)
                print('\nN=%g means desired line-to-line EMF $E_m\\times\\sqrt{3}=%g$ Vrms becomes %g Vrms.' % ((no_series_coil_turns_N), desired_emf_Em*sqrt(3), no_series_coil_turns_N * (2*pi*self.ExcitationFreq * kw1 * air_gap_flux_Phi_m) / sqrt(2)*sqrt(3)), file=fname)
                print('\nEMF induced by one turn is: %g Vrms'% ((2*pi*self.ExcitationFreq * kw1 * air_gap_flux_Phi_m) / sqrt(2)), file=fname)

                if self.DPNV_or_SEPA:
                    number_parallel_branch = 2
                else:
                    number_parallel_branch = 1
                no_conductors_per_slot_zQ = 2* no_phase_m * no_series_coil_turns_N /self.Qs * number_parallel_branch # (7.8) # 2ndEdition (4.63)!!! (2.9) (6.52)

                print('\nNumber of parallel branch: $a=%d$' % number_parallel_branch, file=fname)
                print('\nNumber of conductors per slot: $z_Q=%g$'%no_conductors_per_slot_zQ, file=fname)


    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 10. Tooth Flux Density
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        #       
            # self.stator_tooth_flux_density_B_ds = 1.4 #1.4–2.1 (stator) 
            # rotor_tooth_flux_density_B_dr = 1.5 #1.8 #1.5–2.2 (rotor)
            stator_slot_pitch_tau_us = pi * air_gap_diameter_D / self.Qs
            stator_tooth_apparent_flux_over_slot_pitch_Phi_ds = stack_length_eff * stator_slot_pitch_tau_us * self.guess_air_gap_flux_density
            stator_tooth_width_b_ds = stack_length_eff*stator_slot_pitch_tau_us*self.guess_air_gap_flux_density / (self.lamination_stacking_factor_kFe*stack_length*self.stator_tooth_flux_density_B_ds) + 0.1e-4 # (3.42)

            rotor_slot_pitch_tau_ur = pi * air_gap_diameter_D / self.Qr
            rotor_tooth_apparent_flux_over_slot_pitch_Phi_dr = stack_length_eff*rotor_slot_pitch_tau_ur*self.guess_air_gap_flux_density
            rotor_tooth_width_b_dr = stack_length_eff*rotor_slot_pitch_tau_ur*self.guess_air_gap_flux_density / (self.lamination_stacking_factor_kFe*stack_length*self.rotor_tooth_flux_density_B_dr) + 0.1e-4

            fname = open(one_report_dir_prefix+file_name+'_s10'+file_suffix, 'w', encoding='utf-8') if 'Y730' in pc_name else None
            print(r'''\subsubsection{Tooth Flux Density}
                    Stator tooth.
                    \[\begin{array}{l}
                    {I_u} = {I_s}{z_Q}\\
                    A = {I_u}/{\tau _s} = \frac{{2{I_s}{N_s}m}}{{\pi D}}\\
                    {\tau _s} = \pi D/{Q_s}\\
                    {{\hat B'}_d} = \frac{{{{\hat \Phi '}_d}}}{{{S_d}}} = \frac{{l'{\tau _u}}}{{{k_{Fe}}l\left( {{b_d} - {\rm{0}}{\rm{.1e - 4}}} \right)}}{{\hat B}_\delta }\\
                     \Rightarrow {b_d} = \frac{{l'{\tau _u}}}{{{k_{Fe}}l{{\hat B'}_d}}}{{\hat B}_\delta } + {\rm{0}}{\rm{.1e - 4}}
                    \end{array}\]

                    Similar for rotor tooth.
                    \[{b_{dr}} = \frac{{l'{\tau _{ur}}}}{{{k_{Fe}}l{{\hat B'}_{dr}}}}{{\hat B}_\delta } + {\rm{0}}{\rm{.1e - 4}}\]
                    ''', file=fname)
            print('\nStator tooth flux density $B_{ds}=%g$ T' % (self.stator_tooth_flux_density_B_ds), file=fname)
            print('\nStator tooth width $b_{ds}=%g$ mm---Here, we neglect the flux in stator slot' % (stator_tooth_width_b_ds*1e3), file=fname)
            print('\nStator slot pitch $\\tau_{us}=%g$ mm'    % (stator_slot_pitch_tau_us*1e3), file=fname)
            print('\nRotor tooth flux density $B_{dr}=%g$ T'  % (self.rotor_tooth_flux_density_B_dr), file=fname)
            print('\nRotor tooth width $b_{dr}=%g$ mm---Here, we neglect the flux in rotor slot' % (rotor_tooth_width_b_dr*1e3), file=fname)
            print('\nRotor slot pitch $\\tau_{ur}=%g$ mm'      % (rotor_slot_pitch_tau_ur*1e3), file=fname)



    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 11. Dimension of Slots
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # if bool_standard_voltage_rating:
            #     efficiency = 0.9375 # 500/1.732 * 0.9 / (480/1.732)
            # else:
            #     efficiency = 0.9 # design used in ECCE
            # power_factor = 0.6  #0.85

            stator_phase_current_rms = self.mec_power / (no_phase_m*self.guess_efficiency*stator_phase_voltage_rms*self.guess_power_factor)
            self.stator_phase_current_rms = stator_phase_current_rms

            rotor_current_referred = stator_phase_current_rms * self.guess_power_factor
            rotor_current_actual = no_conductors_per_slot_zQ / number_parallel_branch * self.Qs / self.Qr * rotor_current_referred
            self.rotor_current_actual = rotor_current_actual

            # Stator current density (Pyrhonen09@Example6.4)
            # Js = 3.7e6 # typical value is 3.7e6 Arms/m^2 = Js
            my_AJ_value = self.Js*(linear_current_density_A*1e3)   # \in [9e10, 52e10]
            area_one_conductor_stator_Scs = stator_phase_current_rms / (number_parallel_branch * self.Js) # (7.13)

            # space factor or slot packing factor
            # space_factor_kCu = 0.50 # 不计绝缘的导体填充率：也就是说，一般来说槽满率能达到60%-66%，但是，这里用于下式计算的，要求考虑导体，而槽满率已经考虑了细导线的绝缘了，所以space factor会比槽满率更小，一般在0.5-0.6，低压电机则取下限0.5。
            area_stator_slot_Sus = no_conductors_per_slot_zQ * area_one_conductor_stator_Scs / self.space_factor_kCu 
            print('area_stator_slot_Sus:', area_stator_slot_Sus)

            # guess these local design values or adapt from other designs
            width_statorTeethHeadThickness = 1e-3 # m
            width_StatorTeethNeck = 0.5 * width_statorTeethHeadThickness

            # stator slot height
            stator_inner_radius_r_is_eff = stator_inner_radius_r_is + (width_statorTeethHeadThickness + width_StatorTeethNeck)
            temp = (2*pi*stator_inner_radius_r_is_eff - self.Qs*stator_tooth_width_b_ds)
            stator_tooth_height_h_ds = ( sqrt(temp**2 + 4*pi*area_stator_slot_Sus*self.Qs) - temp ) / (2*pi)
            # The slot depth equals the tooth height hs = hd Pyrhonen09@p309
            stator_slot_height_h_ss = stator_tooth_height_h_ds

            fname = open(one_report_dir_prefix+file_name+'_s11'+file_suffix, 'w', encoding='utf-8') if 'Y730' in pc_name else None
            print(r'''\subsubsection{Dimension of Slots}
                    Guess efficiency ($\eta=0.9$) and power factor (${\rm PF}=0.6$). According to my simulation, it is important that you take a good guess at power factor, or else your design can be low in average torque. For example, if ${\rm PF}$ is 0.8, then the average torque is only 11 Nm, while according to my FEA simulation, a ${\rm PF}=0.6$ is reported and if this value is used, then the average torque reaches 21 Nm. The key physical quantity related to power factor is the rated stator current. That is, PF decreases from 0.8 down to 0.6 and in the meantime stator current 107 A increases up to 151 A.
                    \emph{P.S.: you can take an estimate of PF from the ratio of your designed magnetizing current to the stator current.}

                    The required slot area is derived from the required conductor area and slot packing factor.
                    the required conductor area depends on current density $J$.
                    Particularly, rotor slot height depends on rotor tooth width $b_{dr}$ and rotor current (related to power factor).

                    Define effective radius of inner stator
                    \[{r_{is,{\rm eff}}} = {r_{is}} + {h_{neck}} + {h_{head}}\]
                    \[\text{Slot~Area} = \frac{{\pi {{\left( {{r_{is,{\rm eff}}} + {h_{slot}}} \right)}^2} - \pi {r_{is,{\rm eff}}}^2 - {Q_s}{b_{ds}}{h_{slot}}}}{{{Q_s}}} = {S_{us}}\]

                    The end of a slot means the terminal of a tooth.
                    \[{h_{slot}} = {h_{tooth}}\]

                    Stator slot height.
                    \[\begin{array}{l}
                     \Rightarrow {S_{us}}{Q_s} = 2\pi {r_{is}}{h_{slot}} + \pi h_{slot}^2 - {Q_s}{b_{ds}}{h_{slot}}\\
                     \Rightarrow 0 = \pi h_{slot}^2 + \left( {2\pi {r_{is}} - {Q_s}{b_{ds}}} \right){h_{slot}} - {S_{us}}{Q_s}\\
                    {h_{slot}} = \frac{{ - \left( {2\pi {r_{is}} - {Q_s}{b_{ds}}} \right) \pm \sqrt {{{\left( {2\pi {r_{is}} - {Q_s}{b_{ds}}} \right)}^2} + 4\pi {S_{us}}{Q_s}} }}{{2\pi }}
                    \end{array}\]

                    Similar for rotor.
                    \[{r_{os,{\rm eff}}} = {r_{os}} - {h_{neckhead}}\]
                    \[\begin{array}{l}
                    Rotor\,Slot\,Area = \frac{{\pi r_{os,{\rm eff}}^2 - \pi {{\left( {{r_{or,{\rm eff}}} - {h_{slot,r}}} \right)}^2} - {Q_r}{b_{dr}}{h_{slot,r}}}}{{{Q_r}}} = {S_{ur}}\\
                     \Rightarrow \pi r_{os}^2 - \pi {\left( {{r_{or}} - {h_{slot,r}}} \right)^2} - {Q_r}{b_{dr}}{h_{slot,r}} = {S_{ur}}{Q_r}\\
                     \Rightarrow  - \pi h_{slot,r}^2 + \left( {2\pi {r_{or}} - {Q_r}{b_{dr}}} \right){h_{slot,r}} = {S_{ur}}{Q_r}\\
                     \Rightarrow \pi h_{slot,r}^2 - \left( {2\pi {r_{or}} - {Q_r}{b_{dr}}} \right){h_{slot,r}} + {S_{ur}}{Q_r} = 0\\
                    {h_{slot,r}} = \frac{{\left( {2\pi {r_{or}} - {Q_r}{b_{dr}}} \right) \pm \sqrt {{{\left( {2\pi {r_{or}} - {Q_r}{b_{dr}}} \right)}^2} - 4\pi {S_{ur}}{Q_r}} }}{{2\pi }}
                    \end{array}\]

                    Radius of outer rotor slot.
                    \[\begin{array}{l}
                    2\pi \left( {{r_{or}} - {l_{rhn}} - {r_{rslot}}} \right) - {Q_r}{b_{dr}} = 2{Q_r}{r_{rslot}}\\
                     \Rightarrow {r_{rslot}} = \frac{{2\pi \left( {{r_{or}} - {h_{headneck}}} \right) - {Q_r}{b_{dr}}}}{{\left( {2{Q_r} + 2\pi } \right)}}
                    \end{array}\]

                    Radius of inner rotor slot.
                    \[\begin{array}{l}
                    2\pi \left( {{r_{or}} - {l_{rhn}} + {r_{rslot2}} - {h_{slot,r}}} \right) - {Q_r}{b_{dr}} = 2{Q_r}{r_{rslot2}}\\
                     \Rightarrow {r_{rslot2}} = \frac{{2\pi \left( {{r_{or}} - {l_{rhn}} - {h_{slot,r}}} \right) - {Q_r}{b_{dr}}}}{{2{Q_r} - 2\pi }}
                    \end{array}\]

                    Check in FEMM by integrating block for the rotor slot and rotor current.
                    ''', file=fname)
            print('\n[Guess] Efficiency: $\\eta=%g$'% self.guess_efficiency, file=fname)
            print('\n[Guess] Power factor: $PF=%g$'% self.guess_power_factor, file=fname)
            print('\nStator phase current: $I_1=%g$ Arms'% stator_phase_current_rms, file=fname)
            print('\nRotor current: $I_2^\\prime=%g$ Arms (referred)'% rotor_current_referred, file=fname)
            print('\nRotor current: $I_2^\\prime=%g$ Arms (actual)'% rotor_current_actual, file=fname)
            # quit()
        
            print('\nStator current density: $J_s=%g$ Arms/$\\rm m^2$' % self.Js, file=fname)
            if my_AJ_value*1e-10 < 52 and my_AJ_value*1e-10 > 9:
                print('\nMy $AJ$ value is %ge10 $\\rm A^2/m^3$. It is $\\in$ [9e10, 52e10].'% (my_AJ_value*1e-10), file=fname)
            else:
                print('\nMy $AJ$ value is %ge10 $\\rm A^2/m^3$.'% (my_AJ_value*1e-10), file=fname)
                print('\n[Warning] The liear current density or slot current density is bad.', file=fname)
                # raise Exception('The liear current density or slot current density is bad.')
            print('\nStator one conductor area: $S_{cs}=%g$ ${\\rm mm^2}$'% (area_one_conductor_stator_Scs * 1e6), file=fname)
            print('\nStator space factor: $k_{Cu}=%g$'% self.space_factor_kCu, file=fname)
            print('\nStator slot area: $S_{us}=%g$ ${\\rm mm^2}$ (required)'%( area_stator_slot_Sus*1e6), file=fname)
            print('\nStator tooth height: $h_{ds}=%g$ mm'%(stator_tooth_height_h_ds*1e3), file=fname)
            print('\nStator slot height: $h_{ss}=%g$ mm'% (stator_slot_height_h_ss*1e3), file=fname)

            # Rotor current density
            # we move this loop outside
            # for self.Jr in arange(3e6, 8e6+1, 1e6):
            # if True:
            count_Jr_search = 0
            while True:
                print('-'*20)
                count_Jr_search += 1
                if count_Jr_search>100:
                    if bool_enough_rotor_slot_space == False:
                        print('Lower your self.TangentialStress. If that is not possible, increase your B in the rotor tooth because rotor iron loss is of slip freqnecy---it is not that 怖い.')
                    raise Exception('After 100 search for Jr, there is still fit.')
                    # self.Jr = 8e6 # 8e6 for Cu # 6.5e6 for Al # 减少电流密度有效减少转子欧姆热 # 修正4*pi*area_rotor_slot_Sur*self.Qr处的BUG前，这里设置的转子电流密度都是没用的，和FEA的结果对不上，现在能对上了（用FEMM积分电流和面积验证）！
                area_conductor_rotor_Scr = rotor_current_actual / (1 * self.Jr) # number_parallel_branch_of_rotor_winding=1

                if self.PS_or_SC == True:
                    # space_factor_kAl = 0.6 # We use wound-rotor to implement Chiba's pole selfific rotor. That is, it is coils.
                    # space_factor_kAl = 0.95 # (debug)
                    self.space_factor_kAl = 1.0
                else:
                    self.space_factor_kAl = 1.0 # no clearance for die cast aluminium bar
                     # However, if a cage winding is produced from copper bars by soldering, a clearance of about 0.4mm in width and 1mm in height has to be left in the rotor slot. This clearance also decreases the space factor.
                # no_conductors_per_slot_zQ (z_self.Qr) for cage rotor is clearly 1. 
                # Even though you use coils in rotor (wound rotor), since all the coils within a slot are parallel connected.
                # This will give give a number_parallel_branch equal to z_self.Qr, so they will cancel out, anyway. 
                # Just let both number_parallel_branch and z_self.Qr to be 1 and don't bother.
                area_rotor_slot_Sur = 1 * area_conductor_rotor_Scr / self.space_factor_kAl 

                # guess this local design values or adapt from other designs
                length_headNeckRotorSlot = 1e-3 # for 2 pole motor # 1e-3 m for 4 pole motor

                rotor_outer_radius_r_or_eff = rotor_outer_radius_r_or - length_headNeckRotorSlot
                rotor_tooth_height_h_dr, rotor_tooth_height_h_dr_plus, rotor_Delta = get_parallel_tooth_height(area_rotor_slot_Sur, rotor_tooth_width_b_dr, self.Qr, rotor_outer_radius_r_or_eff)

                if isnan(rotor_tooth_height_h_dr) == True:
                    bool_enough_rotor_slot_space = False
                    print('[%d]'%count_Jr_search, 'Loop for working self.Jr, not', self.Jr, 'Arms/(m^2)')
                    if bool_loop_Jr == True:
                        self.Jr += 0.25e6
                        continue
                    else:
                        return
                    # You can also increase your magnetic load at rotor tooth to create larger rotor slot area.
                    # raise Exception('There are no space on the rotor to fulfill the required rotor slot to reach your J_r and self.space_factor_kAl. Reduce your self.TangentialStress and run the script again.')
                else:
                    bool_enough_rotor_slot_space = True
                    print('[%d]'%count_Jr_search, 'The working self.Jr is', self.Jr)
                    break

            print('\nLoop for searching valid $J_r$...', file=fname)
            print('\nRotor current density (after loop): $J_r=%g$ Arms/${\\rm m^2}$'%self.Jr, file=fname)
            print('\nRotor conductor area: $S_{cr}=%g$ $\\rm mm^2$'% (area_conductor_rotor_Scr*1e6), file=fname)
            print('\nRorot space factor: $k_{Al}=%g$' % self.space_factor_kAl, file=fname)
            print('\nRotor slot area: $S_{ur}=%g$ $\\rm mm^2$'% (area_rotor_slot_Sur*1e6), file=fname)
            print('\nRotor tooth height $h_{dr}(+)=%g$ mm' % (rotor_tooth_height_h_dr_plus*1e3), file=fname) # too large to be right answer
            print('\nRotor tooth height $h_{dr}(-)=%g$ mm' % (rotor_tooth_height_h_dr*1e3), file=fname)
            print('\nQuadratic equation Delta: $\\Delta=%g$'%rotor_Delta, file=fname)
            # The slot depth equals the tooth height hs = hd Pyrhonen09@p309
            rotor_slot_height_h_sr = rotor_tooth_height_h_dr
            print('\nRotor slot height: $h_{sr}=%g$ mm'% (rotor_slot_height_h_sr*1e3), file=fname)

            maximum__Jr = 8e6
            minimum__area_rotor_slot_Sur = rotor_current_actual / (maximum__Jr) / self.space_factor_kAl
            minimum__rotor_tooth_height_h_dr = ( -sqrt(temp**2 - 4*pi*minimum__area_rotor_slot_Sur*self.Qr) + temp ) / (2*pi)
            print('\nMinimum rotor slot area ($J_r=%g$): $S_{ur}=%g$ $\\rm m^2$ $<$ %g $\\rm m^2$'   % (maximum__Jr, minimum__area_rotor_slot_Sur, area_rotor_slot_Sur), file=fname)
            print('\nMinimum rotor tooth height ($J_r=%g$): $h_{dr}=%g$ $\\rm m^2$ $<$ %g $\\rm m^2$'% (maximum__Jr, minimum__rotor_tooth_height_h_dr, rotor_tooth_height_h_dr), file=fname)
            # print('rotor_tooth_width_b_dr', rotor_tooth_width_b_dr*1e3, 'mm', file=fname)
            # quit()


    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 12. Magnetic Voltage
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        #
            if 'M19' in self.Steel:
                # M19-Gauge29 (from FEMM@spmloss example)
                hdata, bdata = np.loadtxt('./M-19-Steel-BH-Curve-afterJMAGsmooth.BH', unpack=True, usecols=(0,1))
                print('The magnetic material is M19-Gauge29.', file=fname)
            elif 'Arnon' in self.Steel and '7' in self.Steel:
                # Arnon-7 (from ELS)
                hdata, bdata = np.loadtxt('../Arnon5/Arnon-7-from-ELS.txt', unpack=True, usecols=(0,1))
                print('The magnetic material is Arnon7.', file=fname)
            else:
                # Arnon-5 (from Ye gu Kang)
                hdata = [0, 9.51030700000000, 11.2124700000000, 13.2194140000000, 15.5852530000000, 18.3712620000000, 21.6562210000000, 25.5213000000000, 30.0619920000000, 35.3642410000000, 41.4304340000000, 48.3863030000000, 56.5103700000000, 66.0660360000000, 77.3405760000000, 90.5910260000000, 106.212089000000, 124.594492000000, 146.311191000000, 172.062470000000, 202.524737000000, 238.525598000000, 281.012026000000, 331.058315000000, 390.144609000000, 459.695344000000, 541.731789000000, 638.410494000000, 752.333643000000, 886.572927000000, 1044.77299700000, 1231.22308000000, 1450.53867000000, 1709.16554500000, 2013.86779200000, 2372.52358500000, 2795.15968800000, 3292.99652700000, 3878.92566000000, 4569.10131700000, 5382.06505800000, 6339.70069300000, 7465.56316200000, 8791.72220000000, 10352.2369750000, 12188.8856750000, 14347.8232500000, 16887.9370500000, 19872.0933000000, 23380.6652750000, 27504.3713250000, 32364.9650250000, 38095.3408000000, 44847.4916750000, 52819.5656250000, 62227.2176750000, 73321.1169500000]
                bdata = [0, 0.0654248125493027, 0.0748613131259592, 0.0852200097732390, 0.0964406675582732, 0.108404414030963, 0.120978202862830, 0.133981410558774, 0.147324354453074, 0.161128351463696, 0.175902377184132, 0.193526821151857, 0.285794748353625, 0.411139883513949, 0.532912618951425, 0.658948953940289, 0.787463844307836, 0.911019620277348, 1.01134216103736, 1.09097860155578, 1.15946725009315, 1.21577636425715, 1.26636706123955, 1.29966244236095, 1.32941739086224, 1.35630922421149, 1.37375630182574, 1.39003487040401, 1.41548927346395, 1.43257623013269, 1.44423937756642, 1.45969672805890, 1.47405771023894, 1.48651531058339, 1.49890498452922, 1.51343941451204, 1.52867783835158, 1.54216506561365, 1.55323686869400, 1.56223503150867, 1.56963683394210, 1.57600636116484, 1.58332795425880, 1.59306861236599, 1.60529276088440, 1.61939615147952, 1.63357053682375, 1.64622605475232, 1.65658227422276, 1.66426678010510, 1.66992280459884, 1.67585542605930, 1.68316554465867, 1.69199548893857, 1.70235212334602, 1.71387033561736, 1.72578827760282]
                print('The magnetic material is Arnon5.', file=fname)
            # for B, H in zip(bdata, hdata):
            #     print B, H, file=fname

            def BH_lookup(B_list, H_list, your_B):
                if your_B<=0:
                    print('positive B only', file=fname)
                    return None
                for ind, B in enumerate(B_list):
                    if your_B > B:
                        continue
                    elif your_B <= B:
                        return (your_B - B_list[ind-1]) / (B-B_list[ind-1]) * (H_list[ind] - H_list[ind-1]) + H_list[ind-1]

                    if ind == len(B_list)-1:
                        slope = (H_list[ind]-H_list[ind-1]) / (B-B_list[ind-1]) 
                        return (your_B - B) * slope + H_list[ind]

            stator_tooth_field_strength_H = BH_lookup(bdata, hdata, self.stator_tooth_flux_density_B_ds)
            rotor_tooth_field_strength_H = BH_lookup(bdata, hdata, self.rotor_tooth_flux_density_B_dr)
            mu0 = 4*pi*1e-7
            air_gap_field_strength_H = self.guess_air_gap_flux_density / mu0

            stator_tooth_magnetic_voltage_Um_ds = stator_tooth_field_strength_H*stator_tooth_height_h_ds
            rotor_tooth_magnetic_voltage_Um_dr = rotor_tooth_field_strength_H*rotor_tooth_height_h_dr

            angle_stator_slot_open = 0.2*(360/self.Qs) / 180 * pi # 参考Chiba的电机槽的开口比例
            b1 = angle_stator_slot_open * stator_inner_radius_r_is
            kappa = b1/air_gap_length_delta / (5 + b1/air_gap_length_delta)
            Carter_factor_kCs = stator_slot_pitch_tau_us / (stator_slot_pitch_tau_us - kappa * b1)

            angle_rotor_slop_open = 0.1*(360/self.Qr) / 180 * pi # 参考Chiba的电机槽的开口比例。Gerada11建议不要开口。
            b1 = angle_rotor_slop_open * rotor_outer_radius_r_or
            kappa = b1/air_gap_length_delta / (5+b1/air_gap_length_delta) 
            Carter_factor_kCr = rotor_slot_pitch_tau_ur / (rotor_slot_pitch_tau_ur - kappa * b1)

            Carter_factor_kC = Carter_factor_kCs * Carter_factor_kCr
            air_gap_length_delta_eff = Carter_factor_kC * air_gap_length_delta
            air_gap_magnetic_voltage_Um_delta = air_gap_field_strength_H * air_gap_length_delta_eff

            fname = open(one_report_dir_prefix+file_name+'_s12'+file_suffix, 'w', encoding='utf-8') if 'Y730' in pc_name else None
            print(r'''\subsubsection{Magnetic Voltage}
                        Look up in BH table for
                        \[\begin{array}{l}
                        {U_{m,ds}} = \\
                        {U_{m,dr}} = \\
                        {U_{m,\delta }} = \frac{{{{\hat H}_\delta }}}{{{\delta _{eff}}}}\\
                        {{\hat H}_\delta } = \frac{{{{\hat B}_\delta }}}{{{\mu _0}}}\\
                        {\delta _{\rm eff}} \buildrel \Delta \over = {\delta _{es}} + {\delta _{er}} = {k_{Cs}}\delta  + {k_{Cr}}\delta \\
                        {k_{Cs}} = \frac{{{\tau _{us}}}}{{{\tau _{us}} - \kappa {b_1}}} = \frac{{{\tau _{us}}}}{{{\tau _{us}} - \frac{{{b_1}/\delta }}{{5 + {b_1}/\delta }}{b_1}}}\\
                        Stator slot open width {b_1} = {\rm{Angle\_StatorSlotOpen}} \times {r_{is}}
                        \end{array}\]''', file=fname)
            print('\nStator tooth field strength: $H=%g$ A/m'% stator_tooth_field_strength_H, file=fname)
            print('\nRotor tooth field strength: $H=%g$ A/m'% rotor_tooth_field_strength_H, file=fname)
            print('\nAir gap field strength: $H_\\delta=%g$ A/m'% air_gap_field_strength_H, file=fname)
            print('\nstator tooth magnetic voltage: $U_{m,ds}=%g$ A'%stator_tooth_magnetic_voltage_Um_ds, file=fname)
            print('\nrotor tooth magnetic voltage: $U_{m,dr}=%g$ A'%rotor_tooth_magnetic_voltage_Um_dr, file=fname)
            print('\nCarter factor: $k_C=k_{Cs} \\times k_{Cr}=%g\\times%g=%g$'% (Carter_factor_kCs, Carter_factor_kCr, Carter_factor_kC), file=fname)
            print('\nEffective air gap length: $\\delta_{\\rm eff}=%g$ mm $>%g$ mm'% (air_gap_length_delta_eff*1e3, air_gap_length_delta*1e3), file=fname)
            print('\nMagnetic voltage drop at air gap: $U_{m\\delta}=%g$ A'% air_gap_magnetic_voltage_Um_delta, file=fname)

            if 'M19' in self.Steel:
                #M19
                coef_eddy = 0.530 # % Eddy current coefficient in (Watt/(meter^3 * T^2 * Hz^2)
                coef_hysteresis = 143.  # % Hysteresis coefficient in (Watts/(meter^3 * T^2 * Hz)
            elif 'Arnon' in self.Steel and '7' in self.Steel:
                #Arnon7
                coef_eddy = 0.07324; # Eddy current coefficient in (Watt/(meter^3 * T^2 * Hz^2)
                coef_hysteresis = 187.6; # Hysteresis coefficient in (Watts/(meter^3 * T^2 * Hz)
            else:
                raise Exception('No loss data for Arnon5.')

            volume_stator_tooth = stator_tooth_width_b_ds * stator_tooth_height_h_ds * stack_length
            volume_rotor_tooth = rotor_tooth_width_b_dr * rotor_tooth_height_h_dr * stack_length

                # 这里在瞎算什么，人家Example 7.4用的是loss table的结果P15，不要混淆了！还有，这个2pi是什么鬼？虽然写的是omega，但是人家是指频率Hz啊！
                # stator_tooth_core_loss_power_per_meter_cube = 1.8 * coef_hysteresis* 2*pi*self.ExcitationFreq * self.stator_tooth_flux_density_B_ds**2 \
                #                                                   + coef_eddy* (2*pi*self.ExcitationFreq * self.stator_tooth_flux_density_B_ds) **2  # table 3.2
                # stator_tooth_core_loss_power = self.Qs* volume_stator_tooth * stator_tooth_core_loss_power_per_meter_cube
                # print 'stator_tooth_core_loss_power=', stator_tooth_core_loss_power, 'W', file=fname

                # rotor_tooth_core_loss_power_per_meter_cube = 1.8 * coef_hysteresis* 2*pi*self.ExcitationFreq * rotor_tooth_flux_density_B_dr**2 \
                #                                                  + coef_eddy* (2*pi*self.ExcitationFreq * rotor_tooth_flux_density_B_dr) **2 
                # rotor_tooth_core_loss_power = self.Qr* volume_rotor_tooth * rotor_tooth_core_loss_power_per_meter_cube
                # print 'rotor_tooth_core_loss_power=', rotor_tooth_core_loss_power, 'W', file=fname



    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 13. (End Loop) Saturation Factor
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        #
            saturation_factor_k_sat = (stator_tooth_magnetic_voltage_Um_ds + rotor_tooth_magnetic_voltage_Um_dr) / air_gap_magnetic_voltage_Um_delta

            tick = 1/60. # Figure 7.2
            k_sat_list =   [   0, 2.5*tick, 5.5*tick, 9.5*tick, 14.5*tick, 20*tick, 28*tick, 37.5*tick, 52.5*tick, 70*tick, 100*tick]
            alpha_i_list = [2/pi,     0.66,     0.68,     0.70,      0.72,    0.74,    0.76,      0.78,      0.80,    0.82,     0.84]

            def alpha_i_lookup(k_sat_list, alpha_i_list, your_k_sat):
                if your_k_sat<=0:
                    print('positive k_sat only', file=fname)
                    return None 
                for ind, k_sat in enumerate(k_sat_list):
                    if your_k_sat > k_sat:
                        if ind == len(k_sat_list)-1: # it is the last one
                            slope = (alpha_i_list[ind]-alpha_i_list[ind-1]) / (k_sat-k_sat_list[ind-1]) 
                            return (your_k_sat - k_sat) * slope + alpha_i_list[ind]
                        else:
                            continue
                    elif your_k_sat <= k_sat:
                        return (your_k_sat - k_sat_list[ind-1]) / (k_sat-k_sat_list[ind-1]) * (alpha_i_list[ind] - alpha_i_list[ind-1]) + alpha_i_list[ind-1]

                # these will not be reached
                print('Reach the end of the curve of k_sat.', file=fname)
                print('End of Loop Error.\n'*3, file=fname)
                return None
            alpha_i_next = alpha_i_lookup(k_sat_list, alpha_i_list, saturation_factor_k_sat)



            fname = open(one_report_dir_prefix+file_name+'_s13'+file_suffix, 'w', encoding='utf-8') if 'Y730' in pc_name else None
            print(r'''\subsubsection{(End Loop) Saturation Factor}
                        \[{k_{sat}} = \frac{{{{\hat U}_{m,ds}} + {{\hat U}_{m,dr}}}}{{{{\hat U}_{m,\delta }}}}\]
                        ''', file=fname)
            print('\nSaturation factor (converged): $k_{sat}=%g$'% saturation_factor_k_sat, file=fname)



            print('alpha_i_next=', alpha_i_next, 'alpha_i=', alpha_i)
            if abs(alpha_i_next - alpha_i)< 1e-3:
                print('alpha_i converges for %d loop' % (loop_count))
                break
            else:
                print('loop for alpha_i_next')
                alpha_i = alpha_i_next
                continue



    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 14. Yoke Geometry (its Magnetic Voltage cannot be caculated because Dse is still unknown)
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
            # if self.p == 1:
            #     stator_yoke_flux_density_Bys = 1.2
            #     rotor_yoke_flux_density_Byr = 1.1 + 0.3
            # else:
            #     stator_yoke_flux_density_Bys = 1.2
            #     rotor_yoke_flux_density_Byr = 1.1

        # compute this again for the new alpha_i and new self.guess_air_gap_flux_density
        air_gap_flux_Phi_m = alpha_i * self.guess_air_gap_flux_density * pole_pitch_tau_p * stack_length_eff
        stator_yoke_height_h_ys = 0.5*air_gap_flux_Phi_m / (self.lamination_stacking_factor_kFe*stack_length*self.stator_yoke_flux_density_Bys)
        rotor_yoke_height_h_yr = 0.5*air_gap_flux_Phi_m / (self.lamination_stacking_factor_kFe*stack_length*self.rotor_yoke_flux_density_Byr)

        fname = open(one_report_dir_prefix+file_name+'_s14'+file_suffix, 'w', encoding='utf-8') if 'Y730' in pc_name else None
        print(r'''\subsubsection{Yoke Geometry}
                \[\begin{array}{l}
                Select\,{{\hat B}_{ys}},{{\hat B}_{yr}}\,in\,Table6.1\\
                {{\hat B}_{ys}} = \frac{{0.5{{\hat \Phi }_m}}}{{{S_{ys}}}} = \frac{{0.5{{\hat \Phi }_m}}}{{{k_{Fe}}\left( {l - {n_v}{b_v}} \right){h_{ys}}}},\,\left( {3.48} \right)\\
                 \Rightarrow {h_{ys}} = \frac{{0.5{{\hat \Phi }_m}}}{{{k_{Fe}}\left( {l - {n_v}{b_v}} \right){{\hat B}_{ys}}}} = \frac{{0.5{{\hat \Phi }_m}}}{{{k_{Fe}}l{{\hat B}_{ys}}}}\\
                {{\hat \Phi }_m} = {\alpha _i}{{\hat B}_\delta }{\tau _p}l',\,\left( {7.5a} \right)\\
                Similarly,\,{h_{yr}} = \frac{{0.5{{\hat \Phi }_m}}}{{{k_{Fe}}\left( {l - {n_v}{b_v}} \right){{\hat B}_{yr}}}} = \frac{{0.5{{\hat \Phi }_m}}}{{{k_{Fe}}l{{\hat B}_{yr}}}}\\
                {\tau _{ys}} = \frac{{\pi \left( {0.5{D_{se}} - 0.5{h_{ys}}} \right)}}{p} = \frac{{\pi \left( {{D_{se}} - {h_{ys}}} \right)}}{{2p}}
                \end{array}\]
                The flux density maxima $\hat B_{ys}$ and $\hat B_{yr}$ of the stator and rotor yokes are selected according to Table 6.1. With the peak value of the flux of the machine, together with the flux density peaks $\hat B_{ys}$ and $\hat B_{yr}$, we are able to determine the heights $h_{ys}$ and $h_{yr}$ of the rotor and stator yokes that realize the selected flux density maxima.
                ''', file=fname)
        print('Stator yoke height: $h_{ys}=%g$ mm'% (stator_yoke_height_h_ys*1e3), file=fname)
        print('Rotor yoke height: $h_{yr}=%g$ mm'% (rotor_yoke_height_h_yr*1e3), file=fname)




    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 15. Machine Main Geometry and Total Magnetic Voltage
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        stator_yoke_diameter_Dsyi = stator_inner_diameter_Dis + 2*stator_tooth_height_h_ds
        stator_outer_diameter_Dse = stator_yoke_diameter_Dsyi + 2*stator_yoke_height_h_ys

        rotor_yoke_diameter_Dryi = rotor_outer_diameter_Dr - 2*rotor_tooth_height_h_dr
        rotor_inner_diameter_Dri = rotor_yoke_diameter_Dryi - 2*rotor_yoke_height_h_yr

        print('Gerada 2011 TIE:')
        rin = rotor_inner_diameter_Dri * 0.5
        rout = rotor_outer_radius_r_or
        r_list = np.arange(rin, rout, 2*1e-3)
        stress_list = check_stress_due_to_centrifugal_force_Gerada11(r_list*1e3, speed_rpm, rout*1e3, rin*1e3)
        for r, stress in zip(r_list, stress_list):
            print('\t%g mm, %g MPa'% (r*1e3, stress))


        fname = open(one_report_dir_prefix+file_name+'_s15'+file_suffix, 'w', encoding='utf-8') if 'Y730' in pc_name else None
        print(r'''\subsubsection{Total Magnetic Voltage \& Machine Geometry}
                    When the air-gap diameter $D_s$, the heights $h_{ds}$ and $h_{dr}$ of the teeth, and the heights $h_{ys}$ and
                    $h_{yr}$ of the stator and rotor yokes are known, we obtain the outer diameter $D_{se}$ of the stator
                    and the inner diameter $D_{ri}$ of the machine; cf. Figures 3.1 and 3.2.
                    \[\begin{array}{l}
                    {D_{se}} - 2{h_{ys}} = {D_{syi}}\\
                    {D_{syi}} - 2{h_{ds}} = {D_s}\\
                    {D_s} - 2\delta  = {D_r}\\
                    {D_r} - 2{h_{dr}} = {D_{ryi}}\\
                    {D_{ryi}} - 2{h_{yr}} = {D_{ri}}
                    \end{array}\]
                    ''', file=fname)
        print('\nStator outer diameter: $Dse=%g$ mm '% (stator_outer_diameter_Dse*1e3), file=fname)
        print('\nRotor inner diameter: $Dri=%g$ mm '% (rotor_inner_diameter_Dri*1e3), file=fname)
        print('\nRotor outer diameter: $D_{or}=%g$ mm'% (1e3*2*rotor_outer_radius_r_or), file=fname)
        # print('Yegu Kang: stator_outer_diameter fixed at 250 mm.')

        # B@yoke
        By_list     = [    0,  0.5,   0.6, 0.7,   0.8,  0.9, 0.95,  1.0, 1.08,   1.1, 1.2,   1.3,   1.4, 1.5,  1.6,  1.7,   1.8,  1.9,   2.0] # Figure 3.17
        coef_c_list = [ 0.72, 0.72, 0.715, 0.7, 0.676, 0.62,  0.6, 0.56,  0.5, 0.485, 0.4, 0.325, 0.255, 0.2, 0.17, 0.15, 0.142, 0.13, 0.125]


        def coef_c_lookup(By_list, coef_c_list, your_By):
            if your_By<=0:
                print('positive By only', file=fname)
                return None 
            for ind, By in enumerate(By_list):
                if your_By > By:
                    continue
                elif your_By <= By:
                    return (your_By - By_list[ind-1]) / (By-By_list[ind-1]) * (coef_c_list[ind] - coef_c_list[ind-1]) + coef_c_list[ind-1]

                if ind == len(By_list)-1:
                    slope = (coef_c_list[ind]-coef_c_list[ind-1]) / (By-By_list[ind-1]) 
                    return (your_By - By) * slope + coef_c_list[ind]

        stator_yoke_field_strength_Hys = BH_lookup(bdata, hdata, self.stator_yoke_flux_density_Bys)
        stator_yoke_middle_pole_pitch_tau_ys = pi*(stator_outer_diameter_Dse - stator_yoke_height_h_ys) / (2*self.p)
        coef_c = coef_c_lookup(By_list, coef_c_list, self.stator_yoke_flux_density_Bys)
        stator_yoke_magnetic_voltage_Um_ys = coef_c * stator_yoke_field_strength_Hys * stator_yoke_middle_pole_pitch_tau_ys

        print('\nStator coeffience c: c=%g'% coef_c, file=fname)
        print('\nStator yoke magnetic voltage: $U_{m,ys}=%g$ A'% stator_yoke_magnetic_voltage_Um_ys, file=fname)

        rotor_yoke_field_strength_Hys = BH_lookup(bdata, hdata, self.rotor_yoke_flux_density_Byr)
        rotor_yoke_middle_pole_pitch_tau_yr = pi*(rotor_yoke_diameter_Dryi - rotor_yoke_height_h_yr) / (2*self.p) # (3.53b)
        coef_c = coef_c_lookup(By_list, coef_c_list, self.rotor_yoke_flux_density_Byr)
        rotor_yoke_magnetic_voltage_Um_yr = coef_c * rotor_yoke_field_strength_Hys * rotor_yoke_middle_pole_pitch_tau_yr

        print('\nRotor coeffience c: c=%g'% coef_c, file=fname)
        print('\nRotor yoke magnetic voltage: $U_{m,yr}=%g$ A'% rotor_yoke_magnetic_voltage_Um_yr, file=fname)



    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 16. Magnetizing Current Providing Total Magnetic Voltage
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        total_magnetic_voltage_Um_tot = air_gap_magnetic_voltage_Um_delta + stator_tooth_magnetic_voltage_Um_ds + rotor_tooth_magnetic_voltage_Um_dr \
                                        + 0.5*stator_yoke_magnetic_voltage_Um_ys + 0.5*rotor_yoke_magnetic_voltage_Um_yr
        stator_magnetizing_current_Is_mag = total_magnetic_voltage_Um_tot * pi* self.p / (no_phase_m * kw1 * no_series_coil_turns_N * sqrt(2))

        fname = open(one_report_dir_prefix+file_name+'_s16'+file_suffix, 'w', encoding='utf-8') if 'Y730' in pc_name else None
        print(r'''\subsubsection{Magnetizing Current Providing Total Magnetic Voltage}
                    \[q = \frac{{{Q_s}}}{{2pm}},q > 1 \Rightarrow {Q_s} = 12,18,24,30,36,42,48\]
                    \[\begin{array}{l}
                    {N_s} = N = pq{z_Q}/a = \frac{{{Q_s}}}{{2m}}{z_Q}/a,\,\left( {2.7} \right),{z_Q} = no.\,conductors\,in\,a\,slot\\
                    {{\hat \Theta }_{s1}} = m{k_{w1}}{N_s}\sqrt 2 {I_s}/\left( {\pi p} \right){\rm{,}}\,{\rm{Equation (2}}{\rm{.15),}}v = 1{\rm{ }}\\
                    {{\hat U}_{m,tot}} = {U_{m,\delta e}} + {{\hat U}_{m,ds}} + {{\hat U}_{m,dr}} + \frac{1}{2}{{\hat U}_{m,ys}} + \frac{1}{2}{{\hat U}_{m,yr}}
                    \end{array}\]
                    \[Let\,{{\hat U}_{m,tot}} = {{\hat \Theta }_{s1}} \Rightarrow {I_{s,mag}} = \frac{{{{\hat U}_{m,tot}}\pi p}}{{m{k_{w1}}{N_s}\sqrt 2 }}\]
                    ''', file=fname)
        print('\nTotal magnetic voltage (half magnetic circuit): $U_{m,tot} =%g$ A'% total_magnetic_voltage_Um_tot, file=fname)
        print('\nStator magnetizing current: $I_{1,mag}=%g$ Arms'% stator_magnetizing_current_Is_mag, file=fname)
        print('\nStator phase current: $I_1=%g$ Arms'% stator_phase_current_rms, file=fname)
        print('\nRotor current: $I_2^\\prime=%g$ Arms (referred)'% rotor_current_referred, file=fname)
        print('\nRotor current: $I_2^\\prime=%g$ Arms (actual)'% rotor_current_actual, file=fname)



    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 17. Losses Computation and Efficiency
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        fname = open(one_report_dir_prefix+file_name+'_s17'+file_suffix, 'w', encoding='utf-8') if 'Y730' in pc_name else None
        print(r'''\subsubsection{Efficiency and Equivalent Circuit Parameters}
                    Since the dimensions have been defined and the winding has been selected, the resistances and inductances of the machine are now calculated. 
                    The magnetizing inductance was discussed in Chapter 3, the leakage inductances in Chapter 4 and the resistances in Chapter 5. 
                    With them, the equivalent circuit parameters of the machine per phase is obtained. 
                    Or, we can resort to eddy current FEA to obtain the equivalent circuit parameters.

                    Now the losses, efficiency, temperature rise and torques of the machine can be determined.
                    ''', file=fname)


    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # 18. Collect Key Geometry Parameters for Performance Evaluation
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        if 'Y730' in pc_name:
            fname.close() # 别忘了关闭文件！
        
        # fig, axes = subplots(1,3, dpi=80)
        # ax = axes[0]
        # ax.plot(hdata, bdata)
        # ax.set_xlabel(r'$H [A/m]$')
        # ax.set_ylabel(r'$B [T]$')
        # ax.grid()
        # ax = axes[1]
        # ax.plot(k_sat_list, alpha_i_list)
        # ax.set_xlabel(r'$k_{sat}$')
        # ax.set_ylabel(r'$\alpha_i$')
        # ax.grid()
        # ax = axes[2]
        # ax.plot(By_list, coef_c_list)
        # ax.set_xlabel(r'$B_y$ [T]')
        # ax.set_ylabel(r'$c$')
        # ax.grid()
        # show()
        print('-'*20)
        fname = None

        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
        # Those variables are used for geometry generation so they are in mm
        #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~

        self.Qs # number of stator slots
        self.Qr # number of rotor slots
        Angle_StatorSlotSpan = 360 / self.Qs # in deg.
        Angle_RotorSlotSpan = 360 / self.Qr # in deg.

        Radius_OuterStatorYoke  = 0.5*stator_outer_diameter_Dse * 1e3
        Radius_InnerStatorYoke  = 0.5*stator_yoke_diameter_Dsyi * 1e3
        Length_AirGap           = air_gap_length_delta * 1e3
        Radius_OuterRotor       = 0.5*rotor_outer_diameter_Dr * 1e3
        Radius_Shaft            = 0.5*rotor_inner_diameter_Dri * 1e3

        Length_HeadNeckRotorSlot = length_headNeckRotorSlot *1e3 # mm # 这里假设与HeadNeck相对的槽的部分也能放导体了。准确来说应该有：rotor_inner_diameter_Dri = rotor_yoke_diameter_Dryi - 2*rotor_yoke_height_h_yr - 2*1e-3*Length_HeadNeckRotorSlot

        rotor_slot_radius = (2*pi*(Radius_OuterRotor - Length_HeadNeckRotorSlot)*1e-3 - rotor_tooth_width_b_dr*self.Qr) / (2*self.Qr+2*pi)
        print('Rotor Slot Radius=', rotor_slot_radius * 1e3, 'mm', file=fname)
        print('rotor_tooth_width_b_dr=', rotor_tooth_width_b_dr * 1e3, 'mm', file=fname)

        Radius_of_RotorSlot = rotor_slot_radius*1e3 
        Location_RotorBarCenter = Radius_OuterRotor - Length_HeadNeckRotorSlot - Radius_of_RotorSlot
        Width_RotorSlotOpen = b1*1e3 # 10% of 360/self.Qr

        # compute Radius_of_RotorSlot2 via the new exact method
        Radius_of_RotorSlot2 = 1e3 * (2*pi*(Radius_OuterRotor - Length_HeadNeckRotorSlot - rotor_slot_height_h_sr*1e3)*1e-3 - rotor_tooth_width_b_dr*self.Qr) / (2*self.Qr-2*pi)
        print('Radius_of_RotorSlot2=', Radius_of_RotorSlot2, file=fname)
        Location_RotorBarCenter2 = Radius_OuterRotor - Length_HeadNeckRotorSlot - rotor_slot_height_h_sr*1e3 + Radius_of_RotorSlot2 

        Angle_StatorSlotOpen = angle_stator_slot_open / pi *180 # in deg.
        Width_StatorTeethBody = stator_tooth_width_b_ds*1e3
        
        Width_StatorTeethHeadThickness = width_statorTeethHeadThickness*1e3 # mm # 这里假设与齿头、齿脖相对的槽的部分也能放导体了。准确来说应该stator_yoke_diameter_Dsyi = stator_inner_diameter_Dis + 2*stator_tooth_height_h_ds + 2*1e-3*(Width_StatorTeethHeadThickness+Width_StatorTeethNeck)
        Width_StatorTeethNeck = 0.5*Width_StatorTeethHeadThickness # mm 



        DriveW_poles=self.p*2
        DriveW_zQ=no_conductors_per_slot_zQ
        DriveW_CurrentAmp = stator_phase_current_rms * sqrt(2); 
        DriveW_Freq = self.ExcitationFreqSimulated
        if 'Cu' not in self.Coil:
            raise Exception('Not implemented error.')
        else:
            print('self.Temperature=', self.Temperature, 'deg Celcius', file=fname)
            rho_Copper = (3.76*self.Temperature+873)*1e-9/55.
            # 参考FEMM_Solver.py 或 按照书上的公式算一下
            stator_slot_area = area_stator_slot_Sus # FEMM是对槽积分得到的，更准哦
            SLOT_FILL_FACTOR = self.space_factor_kCu
            coil_pitch_slot_count = self.winding_layout.coil_pitch # 短距 # self.Qs / DriveW_poles # 整距！
            length_endArcConductor = coil_pitch_slot_count/self.Qs * (0.5*(Radius_OuterRotor + Length_AirGap + Radius_InnerStatorYoke)) * 2*pi # [mm] arc length = pi * diameter  
            length_conductor = (stack_length*1e3 + length_endArcConductor) * 1e-3 # mm to m  ## imagine: two conductors + two end conducotors = one loop (in and out)
            area_conductor   = (stator_slot_area) * SLOT_FILL_FACTOR / DriveW_zQ # TODO: 这里绝缘用槽满率算进去了，但是没有考虑圆形导体之间的空隙？槽满率就是空隙，这里没有考虑绝缘的面积占用。
            # number_parallel_branch = 1
            resistance_per_conductor = rho_Copper * length_conductor / (area_conductor * number_parallel_branch)
        DriveW_Rs = resistance_per_conductor * DriveW_zQ * self.Qs / no_phase_m # resistance per phase
        print('DriveW_Rs=', DriveW_Rs, 'Ohm', file=fname)

        if True:
            with open(self.loc_txt_file, 'a') as f:
                f.write('%d, ' %(self.Qr))
                f.write('%d, %d, ' %(self.Qs, self.Qr))
                f.write('%f, %f, %f, %f, %f, '      % (Radius_OuterStatorYoke, Radius_InnerStatorYoke, Length_AirGap, Radius_OuterRotor, Radius_Shaft))
                f.write('%f, %f, %f, %f, %f, %f, '  % (Length_HeadNeckRotorSlot, Radius_of_RotorSlot, Location_RotorBarCenter, Width_RotorSlotOpen, Radius_of_RotorSlot2, Location_RotorBarCenter2))
                f.write('%f, %f, %f, %f, '          % (Angle_StatorSlotOpen, Width_StatorTeethBody, Width_StatorTeethHeadThickness, Width_StatorTeethNeck))
                f.write('%f, %f, %f, %f, %f, %f, '  % (DriveW_poles, DriveW_zQ, DriveW_Rs, DriveW_CurrentAmp, DriveW_Freq, stack_length*1000))
                f.write('%.14f, %.14f, %.14f, '     % (area_stator_slot_Sus, area_rotor_slot_Sur, minimum__area_rotor_slot_Sur)) # this line exports values need to impose constraints among design parameters for the de optimization
                f.write('%g, %g, %g, %g\n'          % (self.rotor_tooth_flux_density_B_dr, self.stator_tooth_flux_density_B_ds, self.Jr, rotor_tooth_width_b_dr))

        ''' Determine bounds for these parameters:
            stator_tooth_width_b_ds        = design_parameters[0]*1e-3 # m                # stator tooth width [mm]
            air_gap_length_delta           = design_parameters[1]*1e-3 # m                # air gap length [mm]
            b1                             = design_parameters[2]*1e-3 # m                # rotor slot opening [mm]
            rotor_tooth_width_b_dr         = design_parameters[3]*1e-3 # m                # rotor tooth width [mm]
            Length_HeadNeckRotorSlot       = design_parameters[4]             # [4]       # rotor tooth head & neck length [mm]
            Angle_StatorSlotOpen           = design_parameters[5]             # [5]       # stator slot opening [deg]
            Width_StatorTeethHeadThickness = design_parameters[6]             # [6]       # stator tooth head length [mm]
        '''

        self.delta      = 1e3*air_gap_length_delta      # delta or L_g
        self.w_st       = 1e3*stator_tooth_width_b_ds   # w_st
        self.w_rt       = 1e3*rotor_tooth_width_b_dr    # w_rt
        self.theta_so   = Angle_StatorSlotOpen          # theta_so
        self.w_ro       = 1e3*b1                        # w_ro
        self.d_so       = Width_StatorTeethHeadThickness# d_so
        self.d_ro       = Length_HeadNeckRotorSlot      # d_ro
        self.d_st       = 1e3*stator_slot_height_h_ss   # d_st
        self.d_sy       = 1e3*stator_yoke_height_h_ys   # d_sy

        print('-'*20+'\ndelta', '%g mm'%self.delta, file=fname)
        print('w_st',           '%g mm'%self.w_st, file=fname)
        print('w_rt',           '%g mm'%self.w_rt, file=fname)
        print('theta_so',       '%g deg'%self.theta_so, file=fname)
        print('w_ro',           '%g mm'%self.w_ro, file=fname)
        print('d_so',           '%g mm'%self.d_so, file=fname)
        print('d_ro',           '%g mm'%self.d_ro, '\n'+'-'*20, file=fname)
        
        self.x_denorm = [self.delta, self.w_st, self.w_rt, self.theta_so, self.w_ro, self.d_so, self.d_ro]

        self.Radius_of_RotorSlot = Radius_of_RotorSlot # 外圆 of 转子槽

        self.Radius_OuterRotor = Radius_OuterRotor
        self.Radius_OuterStatorYoke = Radius_OuterStatorYoke

        return True if self.Jr_backup < self.Jr else False # bool_bad_specifications

    def build_im_template(self, fea_config_dict):
        import utility
        import population
        # get rid of Swarm class
        self.initial_design_file = '../' + 'pop/' + 'initial_design.txt'

        # load initial design using the obsolete class bearingless_induction_motor_design
        self.im_list = []
        with open(self.initial_design_file, 'r') as f: 
            print(self.initial_design_file)
            for row in utility.csv_row_reader(f):
                im = population.bearingless_induction_motor_design([row[0]]+[float(el) for el in row[1:]], fea_config_dict, model_name_prefix='MOO')
                self.im_list.append(im)
        for im in self.im_list:
            if im.Qr == fea_config_dict['Active_Qr']:
                self.im_template = im
                print('Take the first Active_Qr match as initial design')
                break
        try:
            self.im_template
            self.im_template.Js = self.Js
            self.im_template.fill_factor = self.space_factor_kCu
        except Exception as e:
            print('There is no design matching Active_Qr.')
            msg = 'Please activate one initial design. Refer %s.' % (self.initial_design_file)
            logger = logging.getLogger(__name__)
            logger.warn(msg+str(e))
            print(e)
            raise Exception('no match for Active_Qr: %d'%(fea_config_dict['Active_Qr']))

        # 让儿子能访问爸爸
        self.im_template.spec = self

    def build_pmsm_template(self, fea_config_dict, im_template=None):
        import bearingless_spmsm_design
        self.pmsm_template = bearingless_spmsm_design.bearingless_spmsm_template(fea_config_dict=fea_config_dict)

        if im_template is not None:
            Q = im_template.Qs
            p = im_template.DriveW_poles/2
            self.pmsm_template.deg_alpha_st         = 360/Q - im_template.Angle_StatorSlotOpen
            self.pmsm_template.deg_alpha_so         =                                   self.pmsm_template.deg_alpha_st/2 # im_template uses alpha_so as 0.
            self.pmsm_template.mm_r_si              = im_template.Radius_OuterRotor + im_template.Length_AirGap
            self.pmsm_template.mm_d_so              = im_template.Width_StatorTeethHeadThickness
            self.pmsm_template.mm_d_sp              =                                   1.5*self.pmsm_template.mm_d_so
            self.pmsm_template.mm_d_st              = im_template.Radius_InnerStatorYoke - self.pmsm_template.mm_r_si - self.pmsm_template.mm_d_sp
            self.pmsm_template.mm_d_sy              = im_template.Radius_OuterStatorYoke - im_template.Radius_InnerStatorYoke
            self.pmsm_template.mm_w_st              = self.im_template.Width_StatorTeethBody
            self.pmsm_template.mm_r_st              = 0
            self.pmsm_template.mm_r_sf              = 0
            self.pmsm_template.mm_r_sb              = 0
            self.pmsm_template.Q                    = Q
            self.pmsm_template.sleeve_length        = 3  # mm
            self.pmsm_template.fixed_air_gap_length = 0.75 # mm
            self.pmsm_template.mm_d_pm              = 5  # mm
            self.pmsm_template.deg_alpha_rm         = 90/90*360/(2*p) # deg
            self.pmsm_template.deg_alpha_rs         = self.pmsm_template.deg_alpha_rm # if s=1
            self.pmsm_template.mm_d_ri              = im_template.Location_RotorBarCenter2 - im_template.Radius_of_RotorSlot2 - im_template.Radius_Shaft
            self.pmsm_template.mm_r_ri              = im_template.Radius_Shaft
            self.pmsm_template.mm_d_rp              = 4  # mm
            self.pmsm_template.mm_d_rs              = 0*3
            self.pmsm_template.p                    = p
            self.pmsm_template.s                    = 1

            # Those are some obsolete variables that are convenient to have.
            self.pmsm_template.Radius_OuterStatorYoke = im_template.Radius_OuterStatorYoke
            self.pmsm_template.Radius_OuterRotor      = im_template.Radius_OuterRotor

            # Those variables are for PMSM convenience
            self.rotor_steel_outer_radius             = im_template.Radius_OuterRotor

            # Excitation Properties
            self.pmsm_template.DriveW_Freq       = im_template.DriveW_Freq      
            self.pmsm_template.DriveW_Rs         = im_template.DriveW_Rs        
            self.pmsm_template.DriveW_zQ         = im_template.DriveW_zQ        
            self.pmsm_template.DriveW_CurrentAmp = im_template.DriveW_CurrentAmp
            self.pmsm_template.DriveW_poles      = im_template.DriveW_poles

            print('Pyrhonen TODO')
            self.pmsm_template.Js                = im_template.Js # 4e6 # Arms/mm^2 im_template.Js 
            self.pmsm_template.fill_factor       = im_template.fill_factor # 0.5 # im_template.fill_factor 

            self.pmsm_template.stack_length      = im_template.stack_length 
            self.pmsm_template.wily              = im_template.wily # winding_layout.winding_layout(fea_config_dict['DPNV'], self.Qs, self.p)

            # Specification details:
            #         # 让儿子能访问爸爸
            self.pmsm_template.spec = self


            # free_variables = [None]*13
            # # inherit the stator of IM for the PMSM
            # deg_alpha_st             = free_variables[0]  = 
            # mm_d_so                  = free_variables[1]  = 
            # mm_d_st                  = free_variables[2]  = 
            # stator_outer_radius      = free_variables[3]  = 
            # mm_w_st                  = free_variables[4]  = 
            # sleeve_length            = free_variables[5]  = 
            # mm_d_pm                  = free_variables[6]  = 
            # deg_alpha_rm             = free_variables[7]  = 
            # deg_alpha_rs             = free_variables[8]  = 
            # mm_d_ri                  = free_variables[9]  = 
            # rotor_steel_outer_radius = free_variables[10] =  # the outer radius of the rotor without magnet / the steel
            # mm_d_rp                  = free_variables[11] = 
            # mm_d_rs                  = free_variables[12] = 
            # self.free_variables = free_variables
        else:
            raise Exception('Not implemented error.')

    def get_im_classic_bounds(self, which_filter='FixedStatorSlotDepth', user_bound_filter=None):
        spec = self
        # bound_filter is used to filter out some free_variables that are not going to be optimized.
        self.bound_filter = [ 1,                     # air_gap_length_delta
                              1,                  #--# stator_tooth_width_b_ds
                              1,                  #--# rotor_tooth_width_b_dr
                              1,                     # Angle_StatorSlotOpen
                              1,                     # Width_RotorSlotOpen 
                              1,                     # Width_StatorTeethHeadThickness
                              1,                     # Length_HeadNeckRotorSlot
                              0,                     # Depth_StatorSlot
                              0 ]                    # Depth_StatorYoke (not supported yet)
        if 'FixedStatorSlotDepth' in which_filter:
            pass
        elif 'VariableStatorSlotDepth' in which_filter:
            self.bound_filter[7] = 1
        else:
            if len(user_bound_filter) != 9:
                raise Exception('Invalid bound_filter for bounds. Should be 9 length but %d.'%(len(user_bound_filter)))
            self.bound_filter = user_bound_filter

        self.original_template_neighbor_bounds = [ 
                                [     spec.delta*0.9,       spec.delta*2  ],           # air_gap_length_delta
                                [     spec.w_st *0.5,       spec.w_st *1.5],        #--# stator_tooth_width_b_ds
                                [     spec.w_rt *0.5,       spec.w_rt *1.5],        #--# rotor_tooth_width_b_dr
                                [    0.2*360/spec.Qs,      0.8*360/spec.Qs],           # Angle_StatorSlotOpen
                                [               5e-1,                    3],           # Width_RotorSlotOpen 
                                [               5e-1,                    3],           # Width_StatorTeethHeadThickness
                                [               5e-1,                    3],           # Length_HeadNeckRotorSlot
                                [     spec.d_st *0.8,       spec.d_st *1.2],           # Depth_StatorSlot
                                [     spec.d_sy *0.8,       spec.d_sy *1.2]]           # Depth_StatorYoke

        index_not_included = [idx for idx, el in enumerate(self.bound_filter) if el==0]
        self.filtered_template_neighbor_bounds = [bound for idx, bound in enumerate(self.original_template_neighbor_bounds) if idx not in index_not_included]
        return self.filtered_template_neighbor_bounds


#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# Play with this Pyrhonen procedure
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~

def loop_for_bounds(spec, run_folder):
    spec.loc_txt_file = '../' + 'pop/' + 'loop_for_bounds_%s.txt'%(run_folder[:-1])
    open(spec.loc_txt_file, 'w').close()
    import utility

    # for THE_IM_DESIGN_ID, spec.Qr in enumerate([16,20,28,32,36]): # any spec.Qr>36 will not converge (for alpha_i and k_sat)
    # for THE_IM_DESIGN_ID, spec.Qr in enumerate([32,36]): # any spec.Qr>36 will not converge (for alpha_i and k_sat) with Arnon5 at least
    # for THE_IM_DESIGN_ID, spec.Qr in enumerate([32]):
    bool_run_for_bounds = True
    for rotor_tooth_flux_density_B_dr      in np.arange(1.2, 1.8+0.01, 0.1): #1.5–2.2 (rotor) 
        for stator_tooth_flux_density_B_ds in np.arange(1.2, 1.8+0.01, 0.1): #1.4–2.1 (stator) # too large you will get End of Loop Error (Fixed by extropolating the k_sat vs alpha_i curve.)

            # for spec.Jr in arange(3e6, 8e6+1, 1e6):
            for Jr in np.arange(5.5e6, 8e6+0.01, 0.25e6):
                print(rotor_tooth_flux_density_B_dr, stator_tooth_flux_density_B_ds, Jr)
                utility.blockPrint()

                if not bool_run_for_bounds:
                    rotor_tooth_flux_density_B_dr = 1.5
                    stator_tooth_flux_density_B_ds = 1.4
                    if spec.Qr == 32:
                        Jr = 6.4e6
                    if spec.Qr == 16:
                        Jr = 10.25e6 # for p=1 blim design
                    if spec.Qr == 10:
                        Jr = 9e6 # for p=1 blim design
                try:
                    spec.rotor_tooth_flux_density_B_dr  = rotor_tooth_flux_density_B_dr
                    spec.stator_tooth_flux_density_B_ds = stator_tooth_flux_density_B_ds
                    spec.Jr = Jr
                    # print 
                    bool_bad_specifications = spec.pyrhonen_procedure(bool_loop_Jr=False)
                except Exception as e:
                    raise e
                utility.enablePrint()

                if not bool_run_for_bounds:
                    break
            if not bool_run_for_bounds:
                break
        if not bool_run_for_bounds:
            break

    return spec.loc_txt_file

#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
# Utility for analytical mechanical check
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~

def get_material_data():
    material_density_rho = 7860 # kg/m^3
    Poisson_ratio_nu = 0.29 # (i.e. the ratio of lateral contraction to longitudinal extension in the direction of the stretching force)
    Youngs_modulus_of_elasticity = 190 * 1e9 # Young's modulus for steel is in [190, 210] GPa
    return material_density_rho, Poisson_ratio_nu, Youngs_modulus_of_elasticity

def get_C_prime(Poisson_ratio_nu, bBore):
    if bBore:
        return (3+Poisson_ratio_nu) / 4
    else:
        return (3+Poisson_ratio_nu) / 8
    # Note: C_prime=1 for thin cynlinder
    # -----------------------------------
    # Table 6.4 Poisson’s ratios for certain pure metals
    # Metal ν Metal ν
    # Aluminium Al 0.34 Nickel Ni 0.30
    # Copper Cu 0.34 Titanium Ti 0.34
    # Iron Fe 0.29 Cobalt Co 0.31

def check_stress_due_to_centrifugal_force_Pyrhonen14(speed_rpm, rotor_radius, material_density_rho=None, Poisson_ratio_nu=None, bBore=True):
    Omega = speed_rpm/(60)*2*pi

    if material_density_rho is None: # Example 6.3
        material_density_rho, Poisson_ratio_nu, _ = get_material_data()

    return get_C_prime(Poisson_ratio_nu, bBore) * material_density_rho * rotor_radius**2 * Omega**2

def check_stress_due_to_centrifugal_force_Gerada11(R_list, speed_rpm, Rout, Rin, material_density_rho=None, Poisson_ratio_nu=None):
    Omega = speed_rpm/(60)*2*pi

    if material_density_rho is None: # Example 6.3
        material_density_rho, Poisson_ratio_nu, _ = get_material_data()

    m = inverse_Poisson_ratio_m = 1.0 / Poisson_ratio_nu

    sigma_at_R_list = [ material_density_rho * Omega**2 / (8e12)* ( 
                                                                    (3*m-2)/(m-1) * (Rin**2 + Rout**2 + Rin**2 * Rout**2/R**2) - (m+2)/(m-1)*R**2
                                                                  ) for R in R_list
                      ]

    return sigma_at_R_list

def get_outer_rotor_radius_yield(speed_rpm, yield_stress=None, material_density_rho=None, Poisson_ratio_nu=None, bBore=True, safety_factor_to_yield=1.5): # centrifugal force
    if yield_stress is None: # Example 6.3
        yield_stress = 300e6 # Pa

    yield_stress = yield_stress / safety_factor_to_yield

    Omega = speed_rpm/(60)*2*pi

    if material_density_rho is None:
        material_density_rho, Poisson_ratio_nu, _ = get_material_data()

    rotor_radius_max = sqrt( yield_stress / (get_C_prime(Poisson_ratio_nu, bBore) * material_density_rho * Omega**2) )
    return rotor_radius_max

def eric_specify_tip_speed_get_radius(tip_speed, speed_rpm):
    Omega = speed_rpm/(60)*2*pi
    rotor_radius = tip_speed/Omega
    return rotor_radius

def get_tip_speed(speed_rpm, rotor_radius):
    Omega = speed_rpm/(60)*2*pi
    print('Induction motors with a laminated squirrel cage rotor, rotor surface speed <=200 m/s. --Figure 6.4')
    return rotor_radius * Omega

def get_stack_length_critical_speed(speed_rpm, rotor_radius, material_density_rho=None, Youngs_modulus_of_elasticity=None, safety_factor_to_critical_speed=1.5): # natural frequency
    Omega = speed_rpm/(60)*2*pi

    if material_density_rho is None:
        material_density_rho, _, Youngs_modulus_of_elasticity = get_material_data()

    D_out = rotor_radius * 2
    D_in = 0 # Modifying this is not allowed or else cross section area formula needs to be updated.
    cross_section_area_S = pi*(D_out/2)**2
    second_moment_of_inertia_of_area_I = pi*(D_out**4 - D_in**4) / 64 
    order_critical = 1
    stack_length_max =  sqrt(
                              (order_critical**2 * pi**2) / (safety_factor_to_critical_speed*Omega) \
                              * sqrt( Youngs_modulus_of_elasticity * second_moment_of_inertia_of_area_I \
                                      / (material_density_rho* cross_section_area_S) ) )
    return stack_length_max
    # print('stack_length_max=', stack_length_max*1e3, 'mm', file=fname)

