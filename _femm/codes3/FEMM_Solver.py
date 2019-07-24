#coding:utf-8

import femm
from math import tan, pi, atan, cos, sin, sqrt, copysign, exp
import numpy as np
from csv import reader as csv_reader

import logging
import os
from collections import OrderedDict

import sys
import subprocess

import utility
 # will not create new list as zip does

from time import sleep
from time import time as clock_time

from VanGogh import VanGogh

SELECT_ALL = 4
EPS = 1e-2 # unit mm

class VanGogh_FEMM(VanGogh):
    def __init__(self, im, child_index=0):
        super(VanGogh_FEMM, self).__init__(im, child_index)

    @staticmethod
    def mirror_and_copyrotate(Q, Radius, fraction):
        # Mirror
        femm.mi_selectcircle(0,0,Radius+EPS,SELECT_ALL) # this EPS is sometime necessary to selece the arc at Radius.
        femm.mi_mirror2(0,0,-Radius,0, SELECT_ALL)

        # Rotate
        femm.mi_selectcircle(0,0,Radius+EPS,SELECT_ALL)
        femm.mi_copyrotate2(0, 0, 360./Q, int(Q)/fraction, SELECT_ALL)

    @staticmethod
    def draw_arc(p1, p2, angle, maxseg=1, center=None, **kwarg):
        femm.mi_drawarc(p1[0],p1[1],p2[0],p2[1],angle/pi*180,maxseg) # [deg]
    
    @staticmethod
    def add_arc(p1, p2, angle, maxseg=1, center=None, **kwarg):
        femm.mi_addarc(p1[0],p1[1],p2[0],p2[1],angle/pi*180,maxseg) # [deg]

    @staticmethod
    def draw_line(p1, p2):
        femm.mi_drawline(p1[0],p1[1],p2[0],p2[1])

    @staticmethod
    def add_line(p1, p2):
        femm.mi_addsegment(p1[0],p1[1],p2[0],p2[1])

    def some_solver_related_operations_rotor_before_mirror_rotation(self, im, P6, P8):

        if im.use_drop_shape_rotor_bar == True:
            # constraint to reduce element number @rotor-P6
            femm.mi_selectarcsegment(P6[0], P6[1])
            femm.mi_setarcsegmentprop(8, "<None>", False, 100)
            femm.mi_clearselected()

            # constraint to reduce element number @rotor-P8
            femm.mi_selectarcsegment(P8[0], P8[1])
            femm.mi_setarcsegmentprop(8, "<None>", False, 100)
            femm.mi_clearselected()
        else:
            # constraint to reduce element number @rotor-P8
            femm.mi_selectarcsegment(P8[0], P8[1])
            femm.mi_setarcsegmentprop(8, "<None>", False, 100)
            femm.mi_clearselected()


    def some_solver_related_operations_fraction(self, im, fraction):
        # Boundary
        if fraction == 1:
            femm.mi_drawarc(im.Radius_Shaft,0, -im.Radius_Shaft,0, 180, 20) # 边界不要用太小的segment咯！避免剖分过细（这里设置无效）
            femm.mi_drawarc(-im.Radius_Shaft,0, im.Radius_Shaft,0, 180, 20)
            femm.mi_drawarc(im.Radius_OuterStatorYoke,0, -im.Radius_OuterStatorYoke,0, 180, 20)
            femm.mi_drawarc(-im.Radius_OuterStatorYoke,0, im.Radius_OuterStatorYoke,0, 180, 20)
        elif fraction == 4:
            femm.mi_drawarc(-im.Radius_Shaft,0, 0, -im.Radius_Shaft, 90, 10)
            femm.mi_drawarc(-im.Radius_OuterStatorYoke,0, 0, -im.Radius_OuterStatorYoke, 90, 10)
            femm.mi_selectrectangle(-EPS-im.Radius_Shaft,EPS,EPS-im.Radius_OuterStatorYoke,im.Radius_OuterStatorYoke,SELECT_ALL)
            femm.mi_selectrectangle(EPS,-EPS-im.Radius_Shaft,im.Radius_OuterStatorYoke,EPS-im.Radius_OuterStatorYoke,SELECT_ALL)
            femm.mi_deleteselected()

            # between 2rd and 3th quarters
            p1 = (-im.Location_RotorBarCenter2+im.Radius_of_RotorSlot2, 0)
            p2 = (-im.Radius_Shaft, 0)
            self.add_line(p1, p2)
            p2 = (-im.Location_RotorBarCenter-im.Radius_of_RotorSlot, 0)
            self.add_line(p1, p2)
            p1 = (-im.Radius_OuterRotor-0.5*im.Length_AirGap, 0) # for later extending for moverotate with anti-periodic boundary condition
            self.draw_line(p1, p2)
            p2 = (-im.Radius_OuterRotor-im.Length_AirGap, 0)
            self.draw_line(p1, p2)
            p1 = (-im.Radius_OuterStatorYoke, 0)
            self.add_line(p1, p2)

            # between 3rd and 4th quarters
            p1 = (0, -im.Location_RotorBarCenter2+im.Radius_of_RotorSlot2)
            p2 = (0, -im.Radius_Shaft)
            self.add_line(p1, p2)
            p2 = (0, -im.Location_RotorBarCenter-im.Radius_of_RotorSlot)
            self.add_line(p1, p2)
            p1 = (0, -im.Radius_OuterRotor-0.5*im.Length_AirGap)
            self.draw_line(p1, p2)
            p2 = (0, -im.Radius_OuterRotor-im.Length_AirGap)
            self.draw_line(p1, p2)
            p1 = (0, -im.Radius_OuterStatorYoke)
            self.add_line(p1, p2)
        elif fraction == 2:
            femm.mi_drawarc(-im.Radius_Shaft,0, im.Radius_Shaft,0, 180, 15)
            femm.mi_drawarc(-im.Radius_OuterStatorYoke,0, im.Radius_OuterStatorYoke,0, 180, 15)
            femm.mi_selectrectangle(EPS-im.Radius_OuterStatorYoke,EPS, -EPS+im.Radius_OuterStatorYoke,EPS+im.Radius_OuterStatorYoke, SELECT_ALL)
            femm.mi_deleteselected()

            # between 2rd and 3th quarters
            p1 = (-im.Location_RotorBarCenter2+im.Radius_of_RotorSlot2, 0)
            p2 = (-im.Radius_Shaft, 0)
            self.add_line(p1, p2)
            p2 = (-im.Location_RotorBarCenter-im.Radius_of_RotorSlot, 0)
            self.add_line(p1, p2)
            p1 = (-im.Radius_OuterRotor-0.5*im.Length_AirGap, 0) # for later extending for moverotate with anti-periodic boundary condition
            self.draw_line(p1, p2)
            p2 = (-im.Radius_OuterRotor-im.Length_AirGap, 0)
            self.draw_line(p1, p2)
            p1 = (-im.Radius_OuterStatorYoke, 0)
            self.add_line(p1, p2)

            # between 1rd and 4th quarters
            p1 = (+im.Location_RotorBarCenter2-im.Radius_of_RotorSlot2, 0)
            p2 = (+im.Radius_Shaft, 0)
            self.add_line(p1, p2)
            p2 = (+im.Location_RotorBarCenter+im.Radius_of_RotorSlot, 0)
            self.add_line(p1, p2)
            p1 = (+im.Radius_OuterRotor+0.5*im.Length_AirGap, 0) # for later extending for moverotate with anti-periodic boundary condition
            self.draw_line(p1, p2)
            p2 = (+im.Radius_OuterRotor+im.Length_AirGap, 0)
            self.draw_line(p1, p2)
            p1 = (+im.Radius_OuterStatorYoke, 0)
            self.add_line(p1, p2)
        else:
            raise Exception('not supported fraction = %d' % (fraction))
        # Air Gap Boundary for Rotor Motion #1
        # R = im.Radius_OuterRotor+0.6*im.Length_AirGap
        # femm.mi_drawarc(R,0, -R,0, 180, 5)
        # femm.mi_drawarc(-R,0, R,0, 180, 5)
        # R = im.Radius_OuterRotor+0.4*im.Length_AirGap
        # femm.mi_drawarc(R,0, -R,0, 180, 5)
        # femm.mi_drawarc(-R,0, R,0, 180, 5)

class FEMM_Solver(object):

    def __init__(self, im, flag_read_from_jmag=True, freq=0, individual_index=None, bool_static_fea_loss=False):
        self.bool_static_fea_loss = bool_static_fea_loss
        self.individual_index = individual_index
        self.vangogh = VanGogh_FEMM(im)

        self.deg_per_step = im.fea_config_dict['femm_deg_per_step'] # deg, we need this for show_results
        self.flag_read_from_jmag = flag_read_from_jmag # read the pre-determined rotor currents from the eddy current FEA of jmag or femm

        self.freq = freq
        if freq == 0:
            self.flag_eddycurrent_solver = False
            self.flag_static_solver = not self.flag_eddycurrent_solver
            self.fraction = 1
        else:
            self.flag_eddycurrent_solver = True
            self.flag_static_solver = not self.flag_eddycurrent_solver
            self.fraction = 2

        self.stack_length = im.stack_length

        self.im = im
        self.dir_codes = im.fea_config_dict['dir_codes']

        # evaluate initial design only or optimize
        if im.bool_initial_design == True:
            self.dir_run = im.fea_config_dict['dir_femm_files'] + im.fea_config_dict['model_name_prefix'] + '/'

            if not os.path.exists(self.dir_run):
                logging.getLogger(__name__).debug('FEMM: There is no run yet. Generate the run folder under %s.', self.dir_run)
                os.makedirs(self.dir_run)

            if flag_read_from_jmag == True:
                if self.individual_index is not None:
                    self.dir_run += 'static-jmag/' + 'ind#%04d/'%(self.individual_index)
                else:
                    self.dir_run += 'static-jmag/'
                if not os.path.exists(self.dir_run):
                    os.makedirs(self.dir_run)
            else:
                if self.individual_index is not None:
                    self.dir_run += 'static-femm/' + 'ind#%04d/'%(self.individual_index)
                else:
                    self.dir_run += 'static-femm/'
                if not os.path.exists(self.dir_run):
                    os.makedirs(self.dir_run)
        else:
            if self.individual_index is not None:
                self.dir_run = im.fea_config_dict['dir_femm_files'] + im.fea_config_dict['run_folder'] +  'ind#%04d/'%(self.individual_index)
            else:
                self.dir_run = im.fea_config_dict['dir_femm_files'] + im.fea_config_dict['run_folder']

            if not os.path.exists(self.dir_run):
                logger = logging.getLogger(__name__)
                logger.debug('FEMM: There is no run yet. Generate the run folder as %s.', self.dir_run)
                os.makedirs(self.dir_run)

        self.dir_run_sweeping = self.dir_run + 'sweeping/'
        if not os.path.isdir(self.dir_run_sweeping):
            os.makedirs(self.dir_run_sweeping)

        self.output_file_name = self.get_output_file_name()
        self.rotor_slot_per_pole = int(im.Qr/im.DriveW_poles)
        self.rotor_phase_name_list = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    def add_block_labels(self, fraction=1):
        im = self.im

        SERIES_CONNECTED = 1
        PARALLEL_CONNECTED = 0

        if self.im.fea_config_dict['FEMM_Coarse_Mesh']==True: # Coarse mesh
            MESH_SIZE_ALUMINUM = 2 * 6    # 3
            MESH_SIZE_STEEL    = 2 * 6    # 4
            MESH_SIZE_AIR      = 2 * 0.75 # 0.5 
            MESH_SIZE_COPPER   = 2 * 10   # 8
        else:
            MESH_SIZE_ALUMINUM = 3
            MESH_SIZE_STEEL    = 4
            MESH_SIZE_AIR      = 0.5 
            MESH_SIZE_COPPER   = 8

        def block_label(group_no, material_name, p, meshsize_if_no_automesh, incircuit='<None>', turns=0, automesh=True, magdir=0):
            femm.mi_addblocklabel(p[0],p[1])
            femm.mi_selectlabel(p[0],p[1])
            femm.mi_setblockprop(material_name, automesh, meshsize_if_no_automesh, incircuit, magdir, group_no, turns)
            femm.mi_clearselected()

        # air region @225deg
        X = Y = -(im.Radius_OuterRotor+0.5*im.Length_AirGap) / 1.4142135623730951
        block_label(9, 'Air', (X, Y), MESH_SIZE_AIR, automesh=self.bool_automesh)

        # # Air Gap Boundary for Rotor Motion #2
        # block_label(9, '<No Mesh>',   (0, im.Radius_OuterRotor+0.5*im.Length_AirGap), 5, automesh=self.bool_automesh)
        # block_label(9, 'Air',         (0, im.Radius_OuterRotor+0.7*im.Length_AirGap), 0.5, automesh=self.bool_automesh)
        # block_label(9, 'Air',         (0, im.Radius_OuterRotor+0.3*im.Length_AirGap), 0.5, automesh=self.bool_automesh)


        # shaft
        if fraction == 1:
            block_label(102, '<No Mesh>',         (0, 0),  20)
            # block_label(101, 'Air',         (0, 0),  10, automesh=True) # for deeply-saturated rotor yoke

        # Iron Core @225deg
        if 'M19' in self.im.fea_config_dict['Steel']:
            steel_name = 'M19Gauge29'
        elif self.im.fea_config_dict['Steel'] == 'M15':
            steel_name = 'My M-15 Steel' 
        elif self.im.fea_config_dict['Steel'] == 'Arnon5':
            steel_name = 'Arnon5-final'
        X = Y = -(im.Radius_Shaft+EPS*10) / 1.4142135623730951
        block_label(100, steel_name, (X, Y), MESH_SIZE_STEEL, automesh=self.bool_automesh)
        X = Y = -(0.5*(im.Radius_InnerStatorYoke+im.Radius_OuterStatorYoke)) / 1.4142135623730951
        block_label(10, steel_name, (X, Y), MESH_SIZE_STEEL, automesh=self.bool_automesh)

        # Circuit Configuration
        # Rotor Winding
        if fraction == 1:
            # Pole-Specific Rotor Winding
            # R = 0.5*(im.Location_RotorBarCenter + im.Location_RotorBarCenter2)
            R = im.Location_RotorBarCenter # Since 5/23/2019
            angle_per_slot = 2*pi/im.Qr
            THETA_BAR = pi - angle_per_slot

            for i in range(self.rotor_slot_per_pole):
                circuit_name = 'r%s'%(self.rotor_phase_name_list[i])

                if self.flag_static_solver == True: #self.freq == 0: # Static FEA
                    femm.mi_addcircprop(circuit_name, self.dict_rotor_current_function[i](0.0), SERIES_CONNECTED)
                    # print self.dict_rotor_current_function[i](0.0)
                else:  # Eddy Current FEA (with multi-phase 4-bar cage... haha this is practically nothing)
                    femm.mi_addcircprop(circuit_name, 0, PARALLEL_CONNECTED) # PARALLEL for PS circuit

                THETA_BAR += angle_per_slot

                THETA = THETA_BAR
                X = R*cos(THETA); Y = R*sin(THETA)
                block_label(101, 'Aluminum', (X, Y), MESH_SIZE_ALUMINUM, automesh=self.bool_automesh, incircuit=circuit_name, turns=1)

                THETA = THETA_BAR + angle_per_slot*self.rotor_slot_per_pole
                X = R*cos(THETA); Y = R*sin(THETA)
                block_label(101, 'Aluminum', (X, Y), MESH_SIZE_ALUMINUM, automesh=self.bool_automesh, incircuit=circuit_name, turns=-1)

                THETA = THETA_BAR + angle_per_slot*2*self.rotor_slot_per_pole
                X = R*cos(THETA); Y = R*sin(THETA)
                block_label(101, 'Aluminum', (X, Y), MESH_SIZE_ALUMINUM, automesh=self.bool_automesh, incircuit=circuit_name, turns=1)

                THETA = THETA_BAR + angle_per_slot*3*self.rotor_slot_per_pole
                X = R*cos(THETA); Y = R*sin(THETA)
                block_label(101, 'Aluminum', (X, Y), MESH_SIZE_ALUMINUM, automesh=self.bool_automesh, incircuit=circuit_name, turns=-1)
        elif fraction == 4 or fraction == 2:
            # poly-four-bar-Cage + no bearing current excitated <=> pole specific winding 
            # R = 0.5*(im.Location_RotorBarCenter + im.Location_RotorBarCenter2)
            R = im.Location_RotorBarCenter # Since 5/23/2019
            angle_per_slot = 2*pi/im.Qr
            THETA_BAR = pi - angle_per_slot + EPS # add EPS for the half bar

            for i in range(self.rotor_slot_per_pole):
                circuit_name = 'r%s'%(self.rotor_phase_name_list[i])
                # Eddy Current FEA (with multi-phase 4-bar cage behave the same with PS rotor winding when no bearing current is excited!
                femm.mi_addcircprop(circuit_name, 0, PARALLEL_CONNECTED) # PARALLEL for PS circuit (valid only if there is no 2-pole field)

                THETA_BAR += angle_per_slot

                THETA = THETA_BAR
                X = R*cos(THETA); Y = R*sin(THETA)
                block_label(101, 'Aluminum', (X, Y), MESH_SIZE_ALUMINUM, automesh=self.bool_automesh, incircuit=circuit_name, turns=1) # rA+ ~ rH+

                if fraction == 2:
                    THETA = THETA_BAR + angle_per_slot*self.rotor_slot_per_pole
                    X = R*cos(THETA); Y = R*sin(THETA)
                    block_label(101, 'Aluminum', (X, Y), MESH_SIZE_ALUMINUM, automesh=self.bool_automesh, incircuit=circuit_name, turns=-1) # rA- However, this turns=-1 is not effective for PARALLEL_CONNECTED circuit

            # the other half bar 
            # THETA_BAR += angle_per_slot
            THETA = THETA + angle_per_slot - 2*EPS
            X = R*cos(THETA); Y = R*sin(THETA)
            block_label(101, 'Aluminum', (X, Y), MESH_SIZE_ALUMINUM, automesh=self.bool_automesh, incircuit='r%s'%(self.rotor_phase_name_list[0]), turns=-1) # However, this turns=-1 is not effective for PARALLEL_CONNECTED circuit

        # Stator Winding
        npb = im.wily.number_parallel_branch
        nwl = im.wily.no_winding_layer # number of windign layers 
        if self.flag_static_solver == True: #self.freq == 0: 
            # static solver
            femm.mi_addcircprop('dU', self.dict_stator_current_function[3](0.0), SERIES_CONNECTED)
            femm.mi_addcircprop('dV', self.dict_stator_current_function[4](0.0), SERIES_CONNECTED)
            femm.mi_addcircprop('dW', self.dict_stator_current_function[5](0.0), SERIES_CONNECTED)
            femm.mi_addcircprop('bU', self.dict_stator_current_function[0](0.0), SERIES_CONNECTED)
            femm.mi_addcircprop('bV', self.dict_stator_current_function[1](0.0), SERIES_CONNECTED)
            femm.mi_addcircprop('bW', self.dict_stator_current_function[2](0.0), SERIES_CONNECTED)
        else: # eddy current solver
            # if im.fea_config_dict['DPNV_separate_winding_implementation'] == True or im.fea_config_dict['DPNV'] == False:
            if im.fea_config_dict['DPNV'] == False:
                # either a separate winding or a DPNV winding implemented as a separate winding
                ampD = im.DriveW_CurrentAmp/npb
                ampB = im.BeariW_CurrentAmp
            else:
                # case: DPNV as an actual two layer winding
                ampD = im.DriveW_CurrentAmp/npb
                ampB = ampD

            if im.wily.CommutatingSequenceD == 1:
                MyCommutatingSequence = ['-', '+'] # 2 pole
            else:
                raise
                MyCommutatingSequence = ['+', '-'] # 4 pole legacy

            femm.mi_addcircprop('dU', '%g'                             %(ampD), SERIES_CONNECTED)
            femm.mi_addcircprop('dV', '%g*(-0.5%sI*0.8660254037844386)'%(ampD, MyCommutatingSequence[0]), SERIES_CONNECTED)
            femm.mi_addcircprop('dW', '%g*(-0.5%sI*0.8660254037844386)'%(ampD, MyCommutatingSequence[1]), SERIES_CONNECTED)
            femm.mi_addcircprop('bU', '%g'                             %(ampB), SERIES_CONNECTED)
            femm.mi_addcircprop('bV', '%g*(-0.5%sI*0.8660254037844386)'%(ampB, MyCommutatingSequence[0]), SERIES_CONNECTED)
            femm.mi_addcircprop('bW', '%g*(-0.5%sI*0.8660254037844386)'%(ampB, MyCommutatingSequence[1]), SERIES_CONNECTED)

            # if fraction == 1: # I thought PS can be realized in FEMM but I was wrong, this fraction==1 case should be deleted!
            #     # femm.mi_addcircprop('bA', '%g'                            %(im.BeariW_CurrentAmp), SERIES_CONNECTED)
            #     # femm.mi_addcircprop('bB', '%g*(-0.5+I*0.8660254037844386)'%(im.BeariW_CurrentAmp), SERIES_CONNECTED)
            #     # femm.mi_addcircprop('bC', '%g*(-0.5-I*0.8660254037844386)'%(im.BeariW_CurrentAmp), SERIES_CONNECTED)
            #     femm.mi_addcircprop('bU', '%g'                             %(ampB), SERIES_CONNECTED)
            #     femm.mi_addcircprop('bV', '%g*(-0.5%sI*0.8660254037844386)'%(ampB, MyCommutatingSequence[0]), SERIES_CONNECTED)
            #     femm.mi_addcircprop('bW', '%g*(-0.5%sI*0.8660254037844386)'%(ampB, MyCommutatingSequence[1]), SERIES_CONNECTED)
            # elif fraction == 4 or fraction == 2: # no bearing current
            #     femm.mi_addcircprop('bU', 0, SERIES_CONNECTED)
            #     femm.mi_addcircprop('bV', 0, SERIES_CONNECTED)
            #     femm.mi_addcircprop('bW', 0, SERIES_CONNECTED)

        # dict_dir = {'+':1, '-':-1} # wrong (not consistent with JMAG)
        dict_dir = {'+':-1, '-':1, 'o':0}
        R = 0.5*(im.Radius_OuterRotor + im.Radius_InnerStatorYoke)
        angle_per_slot = 2*pi/im.Qs

        # torque winding's blocks
        THETA = - angle_per_slot + 0.5*angle_per_slot - 3.0/360 # This 3 deg must be less than 360/Qs/2，取这么大是为了在GUI上看得清楚点。
        count = 0
        # for phase, up_or_down in zip(im.l_rightlayer1,im.l_rightlayer2):
        for phase, up_or_down in zip(im.wily.l41,im.wily.l42):
            circuit_name = 'd' + phase
            THETA += angle_per_slot
            X = R*cos(THETA); Y = R*sin(THETA)
            count += 1
            if fraction == 4:
                if not (count > im.Qs*0.5+EPS and count <= im.Qs*0.75+EPS): 
                    continue
            if fraction == 2:
                if not (count > im.Qs*0.5+EPS): 
                    continue
            block_label(11, 'Copper', (X, Y), MESH_SIZE_COPPER, automesh=self.bool_automesh, incircuit=circuit_name, turns=im.DriveW_zQ/nwl*dict_dir[up_or_down])

        # bearing winding's blocks
        if fraction == 1:
            THETA = - angle_per_slot + 0.5*angle_per_slot + 3.0/360

            # for phase, up_or_down in zip(im.l_leftlayer1,im.l_leftlayer2):
            for phase, up_or_down in zip(im.wily.l21,im.wily.l22):
                circuit_name = 'b' + phase
                THETA += angle_per_slot
                X = R*cos(THETA); Y = R*sin(THETA)

                # if self.im.fea_config_dict['DPNV'] == True: 
                # else： # separate winding (e.g., Chiba's)
                block_label(11, 'Copper', (X, Y), MESH_SIZE_COPPER, automesh=self.bool_automesh, incircuit=circuit_name, turns=im.BeariW_turns/nwl*dict_dir[up_or_down])

        elif fraction == 4 or fraction == 2:
            # 危险！FEMM默认把没有设置incircuit的导体都在无限远短接在一起——也就是说，你可能把定子悬浮绕组也短接到鼠笼上去了！
            # 所以，一定要设置好悬浮绕组，而且要用serial-connected，电流给定为 0 A。
            THETA = - angle_per_slot + 0.5*angle_per_slot + 3.0/360
            count = 0
            # for phase, up_or_down in zip(im.l_leftlayer1,im.l_leftlayer2):
            for phase, up_or_down in zip(im.wily.l21,im.wily.l22):
                circuit_name = 'b' + phase
                THETA += angle_per_slot
                X = R*cos(THETA); Y = R*sin(THETA)
                count += 1
                if fraction == 4:
                    if not (count > im.Qs*0.5+EPS and count <= im.Qs*0.75+EPS): 
                        continue
                elif fraction == 2:
                    if not (count > im.Qs*0.5+EPS): 
                        continue
                block_label(11, 'Copper', (X, Y), MESH_SIZE_COPPER, automesh=self.bool_automesh, incircuit=circuit_name, turns=im.BeariW_turns/nwl*dict_dir[up_or_down])

        # Boundary Conditions 
        # femm.mi_makeABC() # open boundary
        if fraction == 1:
            femm.mi_addboundprop('BC:A=0', 0,0,0, 0,0,0,0,0,0,0,0)

            femm.mi_selectarcsegment(0,-im.Radius_OuterStatorYoke)
            femm.mi_setarcsegmentprop(20, "BC:A=0", False, 10) # maxseg = 20 deg (only this is found effective)
            femm.mi_clearselected()
            femm.mi_selectarcsegment(0,im.Radius_OuterStatorYoke)
            femm.mi_setarcsegmentprop(20, "BC:A=0", False, 10)
            femm.mi_clearselected()

            femm.mi_selectarcsegment(0,-im.Radius_Shaft)
            femm.mi_setarcsegmentprop(20, "BC:A=0", False, 100)
            femm.mi_clearselected()
            femm.mi_selectarcsegment(0,im.Radius_Shaft)
            femm.mi_setarcsegmentprop(20, "BC:A=0", False, 100)
            femm.mi_clearselected()
        elif fraction == 4:
            femm.mi_addboundprop('BC:A=0', 0,0,0, 0,0,0,0,0,0,0,0)

            X = Y = -(im.Radius_OuterStatorYoke) / 1.4142135623730951
            femm.mi_selectarcsegment(X, Y)
            femm.mi_setarcsegmentprop(10, "BC:A=0", False, 10) # maxseg = 20 deg (only this is found effective)
            femm.mi_clearselected()

            X = Y = -(im.Radius_Shaft) / 1.4142135623730951
            femm.mi_selectarcsegment(0,-im.Radius_Shaft)
            femm.mi_setarcsegmentprop(10, "BC:A=0", False, 100)
            femm.mi_clearselected()

            femm.mi_addboundprop('apbc1', 0,0,0, 0,0,0,0,0, 5, 0,0)
            femm.mi_addboundprop('apbc2', 0,0,0, 0,0,0,0,0, 5, 0,0)
            femm.mi_addboundprop('apbc3', 0,0,0, 0,0,0,0,0, 5, 0,0)
            femm.mi_addboundprop('apbc5', 0,0,0, 0,0,0,0,0, 5, 0,0)
            femm.mi_addboundprop('apbc6', 0,0,0, 0,0,0,0,0, 5, 0,0) # http://www.femm.info/wiki/periodicboundaries

            R = im.Radius_Shaft+EPS
            femm.mi_selectsegment(0,-R)
            femm.mi_selectsegment(-R,0)
            femm.mi_setsegmentprop("apbc1", 4, False, False, 100)
            femm.mi_clearselected()

            R = im.Location_RotorBarCenter
            femm.mi_selectsegment(0,-R)
            femm.mi_selectsegment(-R,0)
            femm.mi_setsegmentprop("apbc2", 3, False, False, 100)
            femm.mi_clearselected()

            R = im.Radius_OuterRotor + 0.25*im.Length_AirGap
            femm.mi_selectsegment(0,-R)
            femm.mi_selectsegment(-R,0)
            femm.mi_setsegmentprop("apbc3", 0.5, False, False, 9)
            femm.mi_clearselected()

            R = im.Radius_OuterRotor + 0.75*im.Length_AirGap
            femm.mi_selectsegment(0,-R)
            femm.mi_selectsegment(-R,0)
            femm.mi_setsegmentprop("apbc5", 0.5, False, False, 9)
            femm.mi_clearselected()

            R = im.Radius_OuterStatorYoke-EPS
            femm.mi_selectsegment(0,-R)
            femm.mi_selectsegment(-R,0)
            femm.mi_setsegmentprop("apbc6", 4, False, False, 10)
            femm.mi_clearselected()
        elif fraction == 2:
            femm.mi_addboundprop('BC:A=0', 0,0,0, 0,0,0,0,0,0,0,0)

            X = Y = -(im.Radius_OuterStatorYoke) / 1.4142135623730951
            femm.mi_selectarcsegment(X, Y)
            femm.mi_setarcsegmentprop(10, "BC:A=0", False, 10) # maxseg = 20 deg (only this is found effective)
            femm.mi_clearselected()

            X = Y = -(im.Radius_Shaft) / 1.4142135623730951
            femm.mi_selectarcsegment(0,-im.Radius_Shaft)
            femm.mi_setarcsegmentprop(10, "BC:A=0", False, 100)
            femm.mi_clearselected()

            femm.mi_addboundprop('pbc1', 0,0,0, 0,0,0,0,0, 4, 0,0)
            femm.mi_addboundprop('pbc2', 0,0,0, 0,0,0,0,0, 4, 0,0)
            femm.mi_addboundprop('pbc3', 0,0,0, 0,0,0,0,0, 4, 0,0)
            femm.mi_addboundprop('pbc5', 0,0,0, 0,0,0,0,0, 4, 0,0)
            femm.mi_addboundprop('pbc6', 0,0,0, 0,0,0,0,0, 4, 0,0)

            R = im.Radius_Shaft+EPS
            femm.mi_selectsegment(-R,0)
            femm.mi_selectsegment(+R,0)
            femm.mi_setsegmentprop("pbc1", 4, False, False, 100)
            femm.mi_clearselected()

            R = im.Location_RotorBarCenter
            femm.mi_selectsegment(+R,0)
            femm.mi_selectsegment(-R,0)
            femm.mi_setsegmentprop("pbc2", 3, False, False, 100)
            femm.mi_clearselected()

            R = im.Radius_OuterRotor + 0.25*im.Length_AirGap
            femm.mi_selectsegment(+R,0)
            femm.mi_selectsegment(-R,0)
            femm.mi_setsegmentprop("pbc3", 0.5, False, False, 9)
            femm.mi_clearselected()

            R = im.Radius_OuterRotor + 0.75*im.Length_AirGap
            femm.mi_selectsegment(+R,0)
            femm.mi_selectsegment(-R,0)
            femm.mi_setsegmentprop("pbc5", 0.5, False, False, 9)
            femm.mi_clearselected()

            R = im.Radius_OuterStatorYoke-EPS
            femm.mi_selectsegment(+R,0)
            femm.mi_selectsegment(-R,0)
            femm.mi_setsegmentprop("pbc6", 4, False, False, 10)
            femm.mi_clearselected()

        # Air Gap Boundary for Rotor Motion #3
        # inner_angle = 0; outer_angle = 0
        # femm.mi_addboundprop('AGB4RM', 0,0,0, 0,0,0,0,0, 6, inner_angle, outer_angle)
        # R = im.Radius_OuterRotor+0.6*im.Length_AirGap
        # femm.mi_selectarcsegment(0,-R)
        # femm.mi_setarcsegmentprop(5, "AGB4RM", False, 9)
        # femm.mi_clearselected()
        # femm.mi_selectarcsegment(0,R)
        # femm.mi_setarcsegmentprop(5, "AGB4RM", False, 9)
        # femm.mi_clearselected()
        # R = im.Radius_OuterRotor+0.4*im.Length_AirGap
        # femm.mi_selectarcsegment(0,-R)
        # femm.mi_setarcsegmentprop(5, "AGB4RM", False, 9)
        # femm.mi_clearselected()
        # femm.mi_selectarcsegment(0,R)
        # femm.mi_setarcsegmentprop(5, "AGB4RM", False, 9)
        # femm.mi_clearselected()

        # Other arc-segment-specific mesh constraints are already done in draw_model()

    def draw_model(self, fraction=1):

        from shapely.geometry import LineString
        from shapely.geometry import Point

        im = self.im

        origin = Point(0,0)
        Stator_Sector_Angle = 2*pi/im.Qs*0.5
        Rotor_Sector_Angle = 2*pi/im.Qr*0.5

        def mirror_and_copyrotate(Q, Radius, fraction):
            # Mirror
            femm.mi_selectcircle(0,0,Radius+EPS,SELECT_ALL) # this EPS is sometime necessary to selece the arc at Radius.
            femm.mi_mirror2(0,0,-Radius,0, SELECT_ALL)

            # Rotate
            femm.mi_selectcircle(0,0,Radius+EPS,SELECT_ALL)
            femm.mi_copyrotate2(0, 0, 360./Q, int(Q)/fraction, SELECT_ALL)

        def create_circle(p, radius):
            return p.buffer(radius).boundary

        def get_node_at_intersection(c,l): # this works for c and l having one intersection only
            i = c.intersection(l)
            # femm.mi_addnode(i.coords[0][0], i.coords[0][1])
            return i.coords[0][0], i.coords[0][1]

        def draw_arc(p1, p2, angle, maxseg=1):
            femm.mi_drawarc(p1[0],p1[1],p2[0],p2[1],angle/pi*180,maxseg) # [deg]
        def add_arc(p1, p2, angle, maxseg=1):
            femm.mi_addarc(p1[0],p1[1],p2[0],p2[1],angle/pi*180,maxseg) # [deg]
        def draw_line(p1, p2):
            femm.mi_drawline(p1[0],p1[1],p2[0],p2[1])
        def add_line(p1, p2):
            femm.mi_addsegment(p1[0],p1[1],p2[0],p2[1])
        def get_postive_angle(p, origin=(0,0)):
            # using atan loses info about the quadrant
            return atan(abs((p[1]-origin[1]) / (p[0]-origin[0])))


        ''' Part: Stator '''
        # Draw Points as direction of CCW
        # P1
        P1 = (-im.Radius_OuterRotor-im.Length_AirGap, 0)

        # P2
        # Parallel to Line? No they are actually not parallel
        P2_angle = Stator_Sector_Angle -im.Angle_StatorSlotOpen*0.5/180*pi
        k = -tan(P2_angle) # slope
        l_sector_parallel = LineString([(0,0), (-im.Radius_OuterStatorYoke, -im.Radius_OuterStatorYoke*k)])
        c = create_circle(origin, im.Radius_OuterRotor+im.Length_AirGap)
        P2 = get_node_at_intersection(c,l_sector_parallel)
        draw_arc(P2, P1, get_postive_angle(P2))

        # P3
        c = create_circle(origin, im.Radius_OuterRotor+im.Length_AirGap+im.Width_StatorTeethHeadThickness)
        P3 = get_node_at_intersection(c,l_sector_parallel)
        draw_line(P2, P3)

        # P4
        c = create_circle(origin, im.Radius_OuterRotor+im.Length_AirGap+im.Width_StatorTeethHeadThickness+im.Width_StatorTeethNeck)
        l = LineString([(0, 0.5*im.Width_StatorTeethBody), (-im.Radius_OuterStatorYoke, 0.5*im.Width_StatorTeethBody)])
        P4 = get_node_at_intersection(c,l)
        draw_line(P3, P4)

        # P5
        c = create_circle(origin, im.Radius_InnerStatorYoke)
        P5 = get_node_at_intersection(c,l)
        draw_line(P4, P5)

        # P6
        k = -tan(Stator_Sector_Angle)
        l_sector = LineString([(0,0), (-im.Radius_OuterStatorYoke, -im.Radius_OuterStatorYoke*k)])
        P6 = get_node_at_intersection(c,l_sector)
        draw_arc(P6, P5, Stator_Sector_Angle - get_postive_angle(P5))

        # P7
        c = create_circle(origin, im.Radius_OuterStatorYoke)
        P7 = get_node_at_intersection(c,l_sector)
        # draw_line(P6, P7)

        # P8
        P8 = (-im.Radius_OuterStatorYoke, 0)
        # draw_arc(P7, P8, Stator_Sector_Angle)
        # draw_line(P8, P1)

        # P_Coil
        l = LineString([(P3[0], P3[1]), (P3[0], im.Radius_OuterStatorYoke)])
        P_Coil = get_node_at_intersection(l_sector, l)
        draw_line(P4, P_Coil)
        draw_line(P6, P_Coil)

        mirror_and_copyrotate(im.Qs, im.Radius_OuterStatorYoke, fraction)



        ''' Part: Rotor '''
        # Draw Points as direction of CCW
        # P1
        # femm.mi_addnode(-im.Radius_Shaft, 0)
        P1 = (-im.Radius_Shaft, 0)

        # P2
        c = create_circle(origin, im.Radius_Shaft)
        # Line: y = k*x, with k = -tan(2*pi/im.Qr*0.5)
        P2_angle = P3_anlge = Rotor_Sector_Angle
        k = -tan(P2_angle)
        l_sector = LineString([(0,0), (-im.Radius_OuterStatorYoke, -im.Radius_OuterStatorYoke*k)])
        P2 = get_node_at_intersection(c,l_sector)
        # draw_arc(P2, P1, P2_angle)

        # P3
        c = create_circle(origin, im.Radius_OuterRotor)
        P3 = get_node_at_intersection(c,l_sector)
        # draw_line(P2, P3)


        # P4
        l = LineString([(-im.Location_RotorBarCenter, 0.5*im.Width_RotorSlotOpen), (-im.Radius_OuterRotor, 0.5*im.Width_RotorSlotOpen)])
        P4 = get_node_at_intersection(c,l)
        draw_arc(P3, P4, P3_anlge - get_postive_angle(P4))

        # P5
        p = Point(-im.Location_RotorBarCenter, 0)
        c = create_circle(p, im.Radius_of_RotorSlot)
        P5 = get_node_at_intersection(c,l)
        draw_line(P4, P5)

        # P6
        # femm.mi_addnode(-im.Location_RotorBarCenter, im.Radius_of_RotorSlot)
        P6 = (-im.Location_RotorBarCenter, im.Radius_of_RotorSlot)
        draw_arc(P6, P5, 0.5*pi - get_postive_angle(P5, c.centroid.coords[0]))
        # constraint to reduce element number
        femm.mi_selectarcsegment(P6[0], P6[1])
        femm.mi_setarcsegmentprop(8, "<None>", False, 100)
        femm.mi_clearselected()


        # P7
        # femm.mi_addnode(-im.Location_RotorBarCenter2, im.Radius_of_RotorSlot2)
        P7 = (-im.Location_RotorBarCenter2, im.Radius_of_RotorSlot2)
        draw_line(P6, P7)

        # P8
        # femm.mi_addnode(-im.Location_RotorBarCenter2+im.Radius_of_RotorSlot2, 0)
        P8 = (-im.Location_RotorBarCenter2+im.Radius_of_RotorSlot2, 0)
        draw_arc(P8, P7, 0.5*pi)
        # draw_line(P8, P1)
        # constraint to reduce element number
        femm.mi_selectarcsegment(P8[0], P8[1])
        femm.mi_setarcsegmentprop(8, "<None>", False, 100)
        femm.mi_clearselected()

        # P_Bar
        P_Bar = (-im.Location_RotorBarCenter-im.Radius_of_RotorSlot, 0)
        draw_arc(P5, P_Bar, get_postive_angle(P5))
        # add_line(P_Bar, P8)

        mirror_and_copyrotate(im.Qr, im.Radius_OuterRotor, fraction)

        # Boundary
        if fraction == 1:
            femm.mi_drawarc(im.Radius_Shaft,0, -im.Radius_Shaft,0, 180, 20) # 边界不要用太小的segment咯！避免剖分过细（这里设置无效）
            femm.mi_drawarc(-im.Radius_Shaft,0, im.Radius_Shaft,0, 180, 20)
            femm.mi_drawarc(im.Radius_OuterStatorYoke,0, -im.Radius_OuterStatorYoke,0, 180, 20)
            femm.mi_drawarc(-im.Radius_OuterStatorYoke,0, im.Radius_OuterStatorYoke,0, 180, 20)
        elif fraction == 4:
            femm.mi_drawarc(-im.Radius_Shaft,0, 0, -im.Radius_Shaft, 90, 10)
            femm.mi_drawarc(-im.Radius_OuterStatorYoke,0, 0, -im.Radius_OuterStatorYoke, 90, 10)
            femm.mi_selectrectangle(-EPS-im.Radius_Shaft,EPS,EPS-im.Radius_OuterStatorYoke,im.Radius_OuterStatorYoke,SELECT_ALL)
            femm.mi_selectrectangle(EPS,-EPS-im.Radius_Shaft,im.Radius_OuterStatorYoke,EPS-im.Radius_OuterStatorYoke,SELECT_ALL)
            femm.mi_deleteselected()

            # between 2rd and 3th quarters
            p1 = (-im.Location_RotorBarCenter2+im.Radius_of_RotorSlot2, 0)
            p2 = (-im.Radius_Shaft, 0)
            add_line(p1, p2)
            p2 = (-im.Location_RotorBarCenter-im.Radius_of_RotorSlot, 0)
            add_line(p1, p2)
            p1 = (-im.Radius_OuterRotor-0.5*im.Length_AirGap, 0) # for later extending for moverotate with anti-periodic boundary condition
            draw_line(p1, p2)
            p2 = (-im.Radius_OuterRotor-im.Length_AirGap, 0)
            draw_line(p1, p2)
            p1 = (-im.Radius_OuterStatorYoke, 0)
            add_line(p1, p2)

            # between 3rd and 4th quarters
            p1 = (0, -im.Location_RotorBarCenter2+im.Radius_of_RotorSlot2)
            p2 = (0, -im.Radius_Shaft)
            add_line(p1, p2)
            p2 = (0, -im.Location_RotorBarCenter-im.Radius_of_RotorSlot)
            add_line(p1, p2)
            p1 = (0, -im.Radius_OuterRotor-0.5*im.Length_AirGap)
            draw_line(p1, p2)
            p2 = (0, -im.Radius_OuterRotor-im.Length_AirGap)
            draw_line(p1, p2)
            p1 = (0, -im.Radius_OuterStatorYoke)
            add_line(p1, p2)
        elif fraction == 2:
            femm.mi_drawarc(-im.Radius_Shaft,0, im.Radius_Shaft,0, 180, 15)
            femm.mi_drawarc(-im.Radius_OuterStatorYoke,0, im.Radius_OuterStatorYoke,0, 180, 15)
            femm.mi_selectrectangle(EPS-im.Radius_OuterStatorYoke,EPS, -EPS+im.Radius_OuterStatorYoke,EPS+im.Radius_OuterStatorYoke, SELECT_ALL)
            femm.mi_deleteselected()

            # between 2rd and 3th quarters
            p1 = (-im.Location_RotorBarCenter2+im.Radius_of_RotorSlot2, 0)
            p2 = (-im.Radius_Shaft, 0)
            add_line(p1, p2)
            p2 = (-im.Location_RotorBarCenter-im.Radius_of_RotorSlot, 0)
            add_line(p1, p2)
            p1 = (-im.Radius_OuterRotor-0.5*im.Length_AirGap, 0) # for later extending for moverotate with anti-periodic boundary condition
            draw_line(p1, p2)
            p2 = (-im.Radius_OuterRotor-im.Length_AirGap, 0)
            draw_line(p1, p2)
            p1 = (-im.Radius_OuterStatorYoke, 0)
            add_line(p1, p2)

            # between 1rd and 4th quarters
            p1 = (+im.Location_RotorBarCenter2-im.Radius_of_RotorSlot2, 0)
            p2 = (+im.Radius_Shaft, 0)
            add_line(p1, p2)
            p2 = (+im.Location_RotorBarCenter+im.Radius_of_RotorSlot, 0)
            add_line(p1, p2)
            p1 = (+im.Radius_OuterRotor+0.5*im.Length_AirGap, 0) # for later extending for moverotate with anti-periodic boundary condition
            draw_line(p1, p2)
            p2 = (+im.Radius_OuterRotor+im.Length_AirGap, 0)
            draw_line(p1, p2)
            p1 = (+im.Radius_OuterStatorYoke, 0)
            add_line(p1, p2)
        else:
            raise Exception('not supported fraction = %d' % (fraction))
        # Air Gap Boundary for Rotor Motion #1
        # R = im.Radius_OuterRotor+0.6*im.Length_AirGap
        # femm.mi_drawarc(R,0, -R,0, 180, 5)
        # femm.mi_drawarc(-R,0, R,0, 180, 5)
        # R = im.Radius_OuterRotor+0.4*im.Length_AirGap
        # femm.mi_drawarc(R,0, -R,0, 180, 5)
        # femm.mi_drawarc(-R,0, R,0, 180, 5)

    def model_rotor_rotate(self, time):
        if self.deg_per_step != 0.0:
            # 之前用的方法是打开上一个FEM文件，然后将其模型旋转deg_per_step，用不到rotor_position_in_deg的！
            # 当然，我们也试过AirGapBoundary（David Meeker推荐的），转动转子不需要重复剖分，但是计算出来的力不准（转矩是准的）
            # 现在，我们打开在0位置的fem文件，然后转动，savea。这样，就不用不断地打开文件了
            femm.mi_selectgroup(100) # this only select the block labels
            femm.mi_selectgroup(101)
            femm.mi_selectcircle(0,0,self.im.Radius_OuterRotor+EPS,SELECT_ALL) # this selects the nodes, segments, arcs
            femm.mi_moverotate(0,0, self.deg_per_step)
            # femm.mi_zoomnatural()

        # rotor current
        for i in range(self.rotor_slot_per_pole):
            circuit_name = 'r%s'%(self.rotor_phase_name_list[i])
            femm.mi_modifycircprop(circuit_name, 1, self.dict_rotor_current_function[i](time))

        # stator current
        femm.mi_modifycircprop('dU', 1, self.dict_stator_current_function[3](time))
        femm.mi_modifycircprop('dV', 1, self.dict_stator_current_function[4](time))
        femm.mi_modifycircprop('dW', 1, self.dict_stator_current_function[5](time))
        femm.mi_modifycircprop('bU', 1, self.dict_stator_current_function[0](time))
        femm.mi_modifycircprop('bV', 1, self.dict_stator_current_function[1](time))
        femm.mi_modifycircprop('bW', 1, self.dict_stator_current_function[2](time))

    def run_rotating_static_FEA(self): # deg_per_step is key parameter for this function
        # STATIC_RUN_PERIOD = 180 # deg
        STATIC_RUN_PERIOD = 90 # deg

        self.flag_static_solver = True
        self.flag_eddycurrent_solver = False

        femm.openfemm(True) # bHide # False for debug
        femm.newdocument(0) # magnetic
        self.freq = 0 # static 
        self.stack_length = self.im.stack_length * 1
        self.probdef()

        if self.deg_per_step == 0.0:
            print('Locked Rotor! Run 40 stEPS for one slip period.')
            self.im.update_mechanical_parameters(syn_freq=0.0)
        # read currents from previous ec solve
        self.dict_rotor_current_function, self.dict_stator_current_function = self.read_current_from_EC_FEA() # DriveW_Freq and slip_freq_breakdown_torque are used here

        # self.time = 0.0
        # self.rotor_position_in_deg = 0.0
        self.add_material()
        # self.draw_model()
        self.vangogh.draw_model()
        self.add_block_labels()

        # # debug here
        # femm.mi_maximize()
        # femm.mi_zoomnatural()
        # return

        self.output_file_name = self.get_output_file_name()

        if self.deg_per_step == 0.0:
            for i in range(40): # don't forget there
                time += 1.0/self.im.DriveW_Freq / 40. # don't forget here
                # self.rotor_position_in_deg = i # used in file naming
                print(i, time, 's')

                last_out_file_name = output_file_name
                output_file_name = self.output_file_name + '%04d'%(i)
                if os.path.exists(output_file_name + '.ans'):
                    print('.ans file exists. skip this fem file: %s' % (output_file_name))
                    continue
                if last_out_file_name != None:
                    femm.opendocument(last_out_file_name + '.fem')
                    self.model_rotor_rotate(0.0)

                femm.mi_saveas(output_file_name + '.fem') # locked-rotor test
        else: # rotating static FEA

            self.list_rotor_position_in_deg = np.arange(0, STATIC_RUN_PERIOD, self.deg_per_step)
            self.list_name = ['%04d'%(10*el) for el in self.list_rotor_position_in_deg] # with no suffix

            femm.mi_saveas(self.output_file_name + self.list_name[0] + '.fem')
            for rotor_position_in_deg, name in zip(self.list_rotor_position_in_deg[1:], # skip the intial position
                                                    self.list_name[1:]):
                fem_file = self.output_file_name + name + '.fem'
                time = np.abs(rotor_position_in_deg/180*pi / self.im.Omega) # DEBUG: 查了这么久的BUG，原来就是转速用错了！应该用机械转速啊！
                self.model_rotor_rotate(time)
                femm.mi_saveas(fem_file)
                print(time*1e3, 'ms', rotor_position_in_deg, 'deg')
        femm.closefemm()

    def parallel_solve(self, dir_run=None, number_of_instantces=5, bool_watchdog_postproc=True, bool_run_in_JMAG_Script_Editor=False):
        '''[并行求解] 当初没想好，旋转转子竟然不是并行的。。。
        Keyword Arguments:
            dir_run {[type]} -- [静态场用self.dir_run，涡流场用self.dir_run_sweeping] (default: {None})
                                Not not use space in dir_run!!!
                                Not not use space in dir_run!!!
                                Not not use space in dir_run!!!
            number_of_instantces {number} -- [几个？] (default: {5})
            bool_watchdog_postproc {bool} -- [有些时候我们不喜欢用看门狗，后面看了就知道] (default: {True})
        '''
        if dir_run == None:
           dir_run = self.dir_run

        if bool_run_in_JMAG_Script_Editor: # for running script in JMAG, we need a .bat file wrapper for the subprocess calls.
            # os.system('python "%smethod_parallel_solve_4jmag.py" %s' % (self.dir_codes, dir_run))
            with open('temp.bat', 'w') as f:
                if '01' in self.im.fea_config_dict['pc_name']: # python is not added to path in Severson01
                    f.write('"C:/Program Files/JMAG-Designer17.1/python2.7/python" "%smethod_parallel_solve_4jmag.py" %s %d' % (self.dir_codes, dir_run, number_of_instantces))
                else:
                    f.write('python "%smethod_parallel_solve_4jmag.py" %s %d' % (self.dir_codes, dir_run, number_of_instantces))
            os.startfile('temp.bat')
            # os.remove('temp.bat')

        else: # run script in other platforms such as command prompt
            raise Exception('Please explicitly specify bool_run_in_JMAG_Script_Editor.')
            procs = []
            for i in range(number_of_instantces):
                # proc = subprocess.Popen([sys.executable, 'parasolve.py', '{}in.csv'.format(i), '{}out.csv'.format(i)], bufsize=-1)
                proc = subprocess.Popen([sys.executable, 'parasolve.py', str(i), str(number_of_instantces), dir_run], bufsize=-1)
                procs.append(proc)

            for proc in procs:
                proc.wait() # return exit code

        # To collct static results, use while instead, it is more straightforward
        if self.flag_static_solver == True:
            if self.im.fea_config_dict['flag_optimization'] == False: # 优化的话，还是不要用这种看门狗的后处理了，直接求解完就并行后处理。
                self.keep_collecting_static_results_for_optimization(self.list_name, self.list_rotor_position_in_deg)
        return 

        # TODO: loop for post_process results
        # search for python: constently check for new ans file
        if self.im.fea_config_dict['pc_name'] == 'Y730':
            print('Initialize watchdog...')
        else:
            print('watchdog is not installed on servers, quit.')
            return

        return 'Testing JMAG with no watchdog'

        if bool_watchdog_postproc:
            import time
            from watchdog.observers import Observer
            from watchdog.events import FileSystemEventHandler
            class MyHandler(FileSystemEventHandler):
                def __init__(self, solver):
                    self.count_ans = 0
                    self.bool_stop = False
                    self.solver = solver
                    super(FileSystemEventHandler, self).__init__()
                def on_created(self, event):
                    if '.ans' in event.src_path:
                        self.count_ans += 1
                        if self.solver.flag_eddycurrent_solver == True:
                            if self.count_ans == len(self.solver.freq_range):
                                print('\t[Eddy Current Solver] count_ans matches.')
                                self.solver.show_results(bool_plot=True)
                                self.bool_stop = True

                        if self.solver.flag_static_solver == True:
                            # write data to file per 10 .ans files

                            if self.count_ans == self.solver.number_ans or self.count_ans == int(0.5*self.solver.number_ans):
                                print('[Static Solver] count_ans matches for %d' % (self.count_ans))
                                if self.solver.has_results():
                                    self.solver.show_results(bool_plot=True)
                                    self.bool_stop = True
                                else:
                                    self.solver.show_results(bool_plot=False)      
            event_handler = MyHandler(self)
            observer = Observer()
            observer.schedule(event_handler, path=dir_run, recursive=False)
            observer.start()
            try:
                while not event_handler.bool_stop:
                    time.sleep(1)
            except KeyboardInterrupt:
                observer.stop()
            else:
                print('after viewing the plot, watchdog is killed.')
                observer.stop()
            observer.join()

    def read_current_from_EC_FEA(self):
        print('Read rotor current conditions from %s...' % ('JMAG' if self.flag_read_from_jmag else 'FEMM'))
        if self.flag_read_from_jmag == True:
            return self.read_current_conditions_from_JMAG()
        else:
            return self.read_current_conditions_from_FEMM()

    def read_current_conditions_from_FEMM(self):
        self.list_rotor_current_amp = []

        # print 'femm_rotor_current_conditions000' # for noBar test (iron loss with only stator current excitation)
        data = np.loadtxt(self.dir_run_sweeping + 'femm_rotor_current_conditions.txt', unpack=True, usecols=(0,1))
        dict_rotor_current_function = []
        print('[FEMM] Rotor Current') 
        for item in zip(data[0], data[1]):
            item = item[0] + 1j * item[1]
            item *= -1j # -1j is added to be consistent with JMAG (whose Current Source uses sine function)
            amp = np.sqrt(item.imag**2 + item.real**2)
            phase = np.arctan2(item.real, -item.imag) # atan2(y, x), y=a, x=-b
            dict_rotor_current_function.append(lambda t, amp=amp, phase=phase: amp * sin(2*pi*self.im.slip_freq_breakdown_torque*t + phase))
            print('\t', item, amp, phase/pi*180)
            self.list_rotor_current_amp.append(amp)

        dict_stator_current_function = []
        print('[FEMM] Stator Current')                                # -1j is added to be consistent with JMAG
        self.dict_stator_current_from_EC_FEA = [ ('bU', complex(eval( '-1j*%g'                             %(self.im.BeariW_CurrentAmp) ))  ),
                                                 ('bV', complex(eval( '-1j*%g*(-0.5+1j*0.8660254037844386)'%(self.im.BeariW_CurrentAmp) ))  ),
                                                 ('bW', complex(eval( '-1j*%g*(-0.5-1j*0.8660254037844386)'%(self.im.BeariW_CurrentAmp) ))  ),
                                                 ('dU', complex(eval( '-1j*%g'                             %(self.im.DriveW_CurrentAmp) ))  ),
                                                 ('dV', complex(eval( '-1j*%g*(-0.5+1j*0.8660254037844386)'%(self.im.DriveW_CurrentAmp) ))  ),
                                                 ('dW', complex(eval( '-1j*%g*(-0.5-1j*0.8660254037844386)'%(self.im.DriveW_CurrentAmp) ))  )]

        self.dict_stator_current_from_EC_FEA = OrderedDict(self.dict_stator_current_from_EC_FEA)
        dict_stator_current_function = []
        for key, item in self.dict_stator_current_from_EC_FEA.items():
            amp = np.sqrt(item.imag**2 + item.real**2)
            phase = np.arctan2(item.real, -item.imag) # atan2(y, x), y=a, x=-b
            dict_stator_current_function.append(lambda t, amp=amp, phase=phase: amp * sin(2*pi*self.im.DriveW_Freq*t + phase))
            print('\t', key, item, amp, phase/pi*180)

        return dict_rotor_current_function, dict_stator_current_function

    def read_current_conditions_from_JMAG(self):
        try:
            print('The breakdown torque slip frequency is', self.im.slip_freq_breakdown_torque)
        except:
            raise Exception('no breakdown torque slip freqeuncy available.')

        self.dict_rotor_current_from_EC_FEA = []
        self.dict_stator_current_from_EC_FEA = []

        rotor_phase_name_list = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        with open(self.im.csv_previous_solve, 'r') as f:
            read_iterator = csv_reader(f, skipinitialspace=True)
            for row in self.whole_row_reader(read_iterator):
                try:
                    float(row[0])
                except:
                    continue
                else:
                    if np.abs(self.im.slip_freq_breakdown_torque - float(row[0])) < 1e-3:
                        # print row
                        ''' Rotor Current '''
                        beginning_column = 1 + 2*3*2 # title + drive/bearing * 3 phase * real/imag
                        for i in range(0, int(self.im.Qr/self.im.DriveW_poles)):
                            natural_i = i+1
                            current_phase_column = beginning_column + i * int(self.im.DriveW_poles) * 2
                            for j in range(int(self.im.DriveW_poles)):
                                natural_j = j+1
                                re = float(row[current_phase_column+2*j])
                                im = float(row[current_phase_column+2*j+1])
                                self.dict_rotor_current_from_EC_FEA.append( ("r%s%d"%(rotor_phase_name_list[i], natural_j), (re, im)) )

                        ''' Stator Current '''
                        beginning_column = 1 # title column is not needed
                        for i, str_phase in zip(list(range(0,12,2)), ['2A','2B','2C','4A','4B','4C']): # 3 phase
                            natural_i = i+1
                            current_phase_column = beginning_column + i
                            re = float(row[current_phase_column])
                            im = float(row[current_phase_column+1])
                            self.dict_stator_current_from_EC_FEA.append( (str_phase, (re, im)) )

        print('[JMAG] Rotor Current')
        self.list_rotor_current_amp = []
        self.dict_rotor_current_from_EC_FEA = OrderedDict(self.dict_rotor_current_from_EC_FEA)
        dict_rotor_current_function = []
        for key, item in self.dict_rotor_current_from_EC_FEA.items():
            amp = np.sqrt(item[1]**2 + item[0]**2)
            phase = np.arctan2(item[0], -item[1]) # atan2(y, x), y=a, x=-b

            if '1' in key:
                dict_rotor_current_function.append(lambda t, amp=amp, phase=phase: amp * sin(2*pi*self.im.slip_freq_breakdown_torque*t + phase))
                print('\t', key, item, amp, phase/pi*180)
            self.list_rotor_current_amp.append(amp)

        print('[JMAG] Stator Current')
        self.dict_stator_current_from_EC_FEA = OrderedDict(self.dict_stator_current_from_EC_FEA)
        self.dict_stator_current_function_wrong = []
        dict_stator_current_function = []
        for key, item in self.dict_stator_current_from_EC_FEA.items():
            amp = np.sqrt(item[1]**2 + item[0]**2)
            phase = np.arctan2(item[0], -item[1]) # atan2(y, x), y=a, x=-b
            self.dict_stator_current_function_wrong.append(lambda t, amp=amp, phase=phase: amp * sin(2*pi*self.im.DriveW_Freq*t + phase))
            # dict_stator_current_function.append(lambda t, amp=amp, phase=phase: amp * sin(2*pi*self.im.slip_freq_breakdown_torque*t + phase))
            dict_stator_current_function.append(lambda t, amp=amp, phase=phase: amp * sin(2*pi*self.im.DriveW_Freq*t + phase))
            print('\t', key, item, amp, phase/pi*180)
        # dict_stator_current_function =[self.dict_stator_current_function_wrong[0],
        #                                self.dict_stator_current_function_wrong[2],
        #                                self.dict_stator_current_function_wrong[1],
        #                                self.dict_stator_current_function_wrong[3],
        #                                self.dict_stator_current_function_wrong[5],
        #                                self.dict_stator_current_function_wrong[4],]
        # print dict_stator_current_function
        if False:

            from pylab import show, figure
            t = np.arange(0, 0.5, 1e-4)
            # t = np.arange(0, 0.5, 1e-3) # down sampling effect of TranFEAwi2TSS. 5e-2 is too less
            ax1 = figure(1).gca()
            ax2 = figure(2).gca()

            ax = ax1
            rotor_current_one_pole = np.zeros(t.__len__())
            for ind, func in enumerate(dict_rotor_current_function):
                ax.plot(t, [func(el) for el in t], label=ind)
                rotor_current_one_pole += np.array([func(el) for el in t])
            ax.plot(t, rotor_current_one_pole, label='one pole')

            ax = ax2
            iabc_wrong = []
            for ind, func in enumerate(self.dict_stator_current_function_wrong):
                if ind == 3 or ind == 0:
                    ax.plot(t, [-func(el) for el in t], label=str(ind)+'wrong-reversed')
                iabc_wrong.append(np.array([func(el) for el in t]))

            ax = ax1
            iabc = []
            for ind, func in enumerate(self.dict_stator_current_function):
                ax.plot(t, [func(el) for el in t], label=ind)
                iabc.append(np.array([func(el) for el in t]))
            # amplitude-invariant transform - the alpha-beta frame is actually rotating, because iabs is in rotor ref frame (slip frequency)
            ialbe = 2/3.* np.dot(np.array([ [1,      -0.5,       -0.5], 
                                            [0, sqrt(3)/2, -sqrt(3)/2] ]), np.array([iabc[3],iabc[4],iabc[5]]))
            print(np.shape(np.array([iabc[3],iabc[4],iabc[5]])))
            print(np.shape(ialbe))

            print(self.im.omega/2/pi)
            ids = []
            iqs = []
            ''' Speed negation is done by ? '''
            iabc_stationary_2_tor = []
            iabc_stationary_2_sus = []
            for i in range(len(t)):
                # theta = -self.im.omega * t[i]
                # theta = -self.im.DriveW_Freq*2*pi * t[i]
                # theta = self.im.omega * t[i]
                theta = self.im.DriveW_Freq*2*pi * t[i]

                # turn into stationary dq frame
                temp = np.dot( np.array([ [np.cos(theta), -np.sin(theta)], 
                                          [np.sin(theta),  np.cos(theta)] ]), np.array([ [ialbe[0][i]], 
                                                                                         [ialbe[1][i]] ]) )
                ids.append(temp[0])
                iqs.append(temp[1])

                iabc_stationary_2_tor.append(self.transformed_dict_stator_current_function(3, t[i], theta))
                iabc_stationary_2_sus.append(self.transformed_dict_stator_current_function(0, t[i], theta))

            ids = np.array(ids).T[0]
            iqs = np.array(iqs).T[0]
            # print ids
            # print iqs
            print('idq', np.shape(ids), np.shape(iqs))
            # ax_r.plot(t, idq[0], label='i_ds')
            # ax_r.plot(t, idq[1], label='i_qs')
            # ax_r.plot(t, ialbe[0], label='alpha')
            # ax_r.plot(t, ialbe[1], label='beta')




            # tansform to phase coordinates
            ax = ax2
            iabc_stationary = 1.5 * np.dot(np.array([ [ 2/3.,          0], 
                                                      [-1/3.,  sqrt(3)/3],
                                                      [-1/3., -sqrt(3)/3] ]), np.array([ids, iqs]))
            ax.plot(t, iabc_stationary[0], label='i_a')
            # ax.plot(t, iabc_stationary[1], label='i_b')
            # ax.plot(t, iabc_stationary[2], label='i_c')

            ax.plot(t, [el[0] for el in iabc_stationary_2_sus], label='i_a_2_sus')

            ax.plot(t, [el[0] for el in iabc_stationary_2_tor], label='i_a_2_tor')
            # ax.plot(t, [el[1]+5 for el in iabc_stationary_2_tor], label='i_b_2')
            # ax.plot(t, [el[2]+5 for el in iabc_stationary_2_tor], label='i_c_2')

            ax1.legend()
            ax2.legend()
            show()

            quit()

        return dict_rotor_current_function, dict_stator_current_function

    def transformed_dict_stator_current_function(self, index_phase_A, time, theta):
        ia_slip_freq = self.dict_stator_current_function[index_phase_A](time)
        ib_slip_freq = self.dict_stator_current_function[index_phase_A+1](time)
        ic_slip_freq = self.dict_stator_current_function[index_phase_A+2](time)

        iabc_vector_slip_freq = np.array([[ia_slip_freq, ib_slip_freq, ic_slip_freq]]).T
        ialbe = 2/3.* np.dot(np.array([ [1,      -0.5,       -0.5], 
                                        [0, sqrt(3)/2, -sqrt(3)/2] ]), iabc_vector_slip_freq)
        # print 'ialbe', ialbe

        # turn into stationary dq frame
        idq = np.dot( np.array([ [np.cos(theta), -np.sin(theta)], 
                                 [np.sin(theta),  np.cos(theta)] ]), ialbe )
        # print 'idq', idq

        iabc_stationary = 1.5 * np.dot(np.array([ [ 2/3.,          0], 
                                                  [-1/3.,  sqrt(3)/3],
                                                  [-1/3., -sqrt(3)/3] ]), idq)

        return iabc_stationary[0][0], iabc_stationary[1][0], iabc_stationary[2][0]

    def get_air_gap_B(self, number_of_points=360):
        im = self.im
        femm.opendocument(self.output_file_name + '.fem')
        femm.mi_loadsolution()

        list_B_magitude = []
        R = im.Radius_OuterRotor + 0.25*im.Length_AirGap
        for i in range(number_of_points):
            THETA = i / 180.0 * pi
            X = R*cos(THETA)
            Y = R*sin(THETA)
            B_vector_complex = femm.mo_getb(X, Y)
            B_X_complex = B_vector_complex[0]
            B_Y_complex = B_vector_complex[1]
            B_X_real = np.real(B_X_complex)
            B_Y_real = np.real(B_Y_complex)
            # Assume the magnitude is all due to radial component
            B_magitude = sqrt(B_X_real**2 + B_Y_real**2)
            inner_product = B_X_real * X + B_Y_real *Y
            list_B_magitude.append( B_magitude * copysign(1, inner_product) )
        return list_B_magitude

    def probdef(self):
        # femm.smartmesh(False) <- This will not work due to bug of femm.__init__  # let mi_smartmesh deside. You must turn it off in parasolver.py
        femm.callfemm_noeval('smartmesh(0)') # call this after probdef
        femm.mi_probdef(self.freq, 'millimeters', 'planar', 1e-8, # must < 1e-8
                        self.stack_length, 18, 1) # The acsolver parameter (default: 0) specifies which solver is to be used for AC problems: 0 for successive approximation, 1 for Newton.
                                             # 1 for 'I intend to try the acsolver of Newton, as this is the default for JMAG@[Nonlinear Calculation] Setting Panel in the [Study Properties] Dialog Box'
        # femm.callfemm_noeval('mi_smartmesh(0)') # call this after probdef
        self.bool_automesh = False # setting to false gives no effect?

        # femm.smartmesh(True) # let mi_smartmesh deside. You must turn it off in parasolver.py
        # self.bool_automesh = True # setting to false gives no effect?

    def add_material(self):
        # mi_addmaterial('matname', mu x, mu y, H c, J, Cduct, Lam d, Phi hmax, lam fill, LamType, Phi hx, Phi hy, nstr, dwire)
        femm.mi_getmaterial('Air')
        femm.mi_getmaterial('Copper') # for coil
        # femm.mi_getmaterial('18 AWG') # for coil
        # femm.mi_getmaterial('Aluminum, 1100') # for bar?
        # femm.mi_getmaterial('304 Stainless Steel') # for shaft?

        # femm.mi_addmaterial('Air', 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0);
        # femm.mi_addmaterial('Aluminum', 1, 1, 0, 0, 35, 0, 0, 1, 0, 0, 0)
        femm.mi_addmaterial('Aluminum', 1, 1, 0, 0, self.im.fea_config_dict['Bar_Conductivity']*1e-6, 0, 0, 1, 0, 0, 0) # [MS/m]
        # femm.mi_addmaterial('Aluminum', 1, 1, 0, 0, 1/1.673e-2, 0, 0, 1, 0, 0, 0)

        # femm.mi_addmaterial('LinearIron', 2000, 2000, 0, 0, 0, 0, 0, 1, 0, 0, 0);

        if self.im.fea_config_dict['Steel'] == 'M19Gauge29':
            # femm.mi_getmaterial('M-19 Steel') # for Stator & Rotor Iron Cores (Nonlinear with B-H curve)
            femm.mi_addmaterial('M19Gauge29',0,0, 0,0, 0,0.3556,0, 0.95) # no lamination for testing consistency with JMAG
            hdata, bdata = np.loadtxt(self.dir_codes + './M-19-Steel-BH-Curve-afterJMAGsmooth.BH', unpack=True, usecols=(0,1))
            for n in range(0,len(bdata)):
                femm.mi_addbhpoint('M19Gauge29', bdata[n], hdata[n])

        elif self.im.fea_config_dict['Steel'] == 'Arnon5':
            # Arnon5 is 1/5 thick as M15, which is too thin to use and it is expensive as well
            femm.mi_addmaterial('Arnon5-final',0,0, 0,0, 0.0,0.127,0, 0.96)
            BH = np.loadtxt(self.dir_codes + '../Arnon5/Arnon5-final.txt', unpack=True, usecols=(0,1))
            bdata = BH[1][1:] # if not skip the first point, there will be two (0,0) in FEMM software, reason unknown.
            hdata = BH[0][1:] # if not skip the first point, there will be two (0,0) in FEMM software, reason unknown.
            for n in range(0,len(bdata)):
                femm.mi_addbhpoint('Arnon5-final', bdata[n], hdata[n])

        elif self.im.fea_config_dict['Steel'] == 'M15':
            femm.mi_addmaterial('My M-15 Steel',0,0, 0,0, 0,0.635,0, 0.98)
            BH = np.loadtxt(self.dir_codes + '../Arnon5/M-15-Steel-BH-Curve.txt', unpack=True, usecols=(0,1))
            bdata = BH[1]
            hdata = BH[0]
            for n in range(0,len(bdata)):
                femm.mi_addbhpoint('My M-15 Steel', bdata[n], hdata[n])


        if False:
            # A more interesting material to add is the iron with a nonlinear
            # BH curve.  First, we create a material in the same way as if we
            # were creating a linear material, except the values used for
            # permeability are merely placeholders.
            femm.mi_addmaterial('Arnon5',0,0,0,0, 0.0,0.127,0,0.96)
            # A set of points defining the BH curve is then specified.
            BH = np.loadtxt(self.dir_codes + 'Arnon5_Kang_after_JMAG_Smoothed.txt', unpack=True, usecols=(0,1))
            bdata = BH[1]
            hdata = BH[0]
            for n in range(0,len(bdata)):
                femm.mi_addbhpoint('Arnon5', bdata[n], hdata[n])

    def get_output_file_name(self, booL_dir=True):
        fname ='%s-%gHz'%(self.im.ID, self.freq)
        if booL_dir == True:
            self.output_file_name = self.dir_run + fname
            return self.output_file_name
        else:
            return fname

    def whole_row_reader(self, reader):
        for row in reader:
            yield row[:]

    def has_results(self, dir_run=None):
        # print 'self.freq', self.freq
        if dir_run == None:
            dir_run = self.dir_run

        a = [f for f in os.listdir(dir_run) if '.ans' in f].__len__()
        b = [f for f in os.listdir(dir_run) if '.fem' in f].__len__()
        if a == 0:
            if 'no_duplicates.txt' in os.listdir(dir_run):
                return True # 直接把femm的结果从服务器上拷过来也行
            else:
                return False
        print('[FEMM.has_results] ans count: %d. fem count: %d.' % (a, b))
        return a == b

    def show_results(self, bool_plot=True):

        if self.flag_eddycurrent_solver == True:
            print('show results for eddy current solver')
            return self.show_results_eddycurrent(bool_plot=bool_plot)

        if self.flag_static_solver == True:
            print('show results for static solver')
            return self.show_results_static(bool_plot=bool_plot)

        return None

    def show_results_eddycurrent(self, bool_plot):
        if self.fraction == 1:
            self.fraction = 2 # this fixes a bug calling femm_integrate_4_current without calling run_frequency_sweeping before
            raise Exception('You should initialize FEMM Solver with freq!=0')
    
        ans_file_list = os.listdir(self.dir_run_sweeping)
        ans_file_list = [f for f in ans_file_list if '.ans' in f]

        femm.openfemm(True)
        list_ans_file = []
        list_torque = []
        # write results to a data file
        str_results = ''
        for ind, f in enumerate( ans_file_list ):
            # femm.opendocument(self.dir_run_sweeping + f[:-4] + '.fem')
            # femm.mi_loadsolution()
            femm.opendocument(self.dir_run_sweeping + f)

            # physical amount on rotor
            femm.mo_groupselectblock(100)
            femm.mo_groupselectblock(101)
            Fx = femm.mo_blockintegral(18) #-- 18 x (or r) part of steady-state weighted stress tensor force
            Fy = femm.mo_blockintegral(19) #--19 y (or z) part of steady-state weighted stress tensor force
            torque = femm.mo_blockintegral(22) #-- 22 = Steady-state weighted stress tensor torque
            femm.mo_clearblock()

            # rotor current
            # _ = self.femm_integrate_4_current()
            # # TODO: this is for testing phase
            # if float(f[3:-6]) == 3:
            #     print '\n', f[3:-6], 'Hz'
            #     for el in vals_results_rotor_current:
            #         print abs(el)

            str_results += "%s %g %g %g\n" % ( f[3:-6], torque, Fx, Fy ) 
            list_ans_file.append(f)
            list_torque.append(torque)
            femm.mo_close() 
        with open(self.dir_run_sweeping + "eddycurrent_results.txt", "w") as stream:
            stream.write(str_results)

        # find breakdown torque and slip frequency that we are interested
        index, breakdown_torque = utility.get_max_and_index(list_torque)
        slip_freq_breakdown_torque = list_ans_file[index][3:-6]
        print("FEMM's breakdown data: %s Hz, %g Nm" % (slip_freq_breakdown_torque, breakdown_torque))
        self.im.update_mechanical_parameters(float(slip_freq_breakdown_torque))
        # if self.im.slip_freq_breakdown_torque != float(slip_freq_breakdown_torque):
        #     raise Exception('[DEBUG] JMAG disagrees with FEMM (!= 3 Hz).')

        # write rotor currents to file
        self.femm_integrate_4_current(self.dir_run_sweeping + list_ans_file[index], self.fraction)

        femm.closefemm()

    def femm_integrate_4_current(self, fname, fraction, dir_output=None, returnData=False):
        '''Make sure femm is opened
        Returns:
            [type] -- [list of complex number of rotor currents from FEMM]
        '''

        # get corresponding rotor current conditions for later static FEA
        femm.opendocument(fname)

        if True:
            # physical amount of Cage
            im = self.im
            vals_results_rotor_current = []
            # R = 0.5*(im.Location_RotorBarCenter + im.Location_RotorBarCenter2)
            R = im.Location_RotorBarCenter # Since 5/23/2019
            angle_per_slot = 2*pi/im.Qr
            THETA_BAR = pi - angle_per_slot + EPS # add EPS for the half bar
            # print 'number of rotor_slot per partial model', self.rotor_slot_per_pole * int(4/fraction)
            for i in range(self.rotor_slot_per_pole * int(4/fraction)):
                THETA_BAR += angle_per_slot
                THETA = THETA_BAR
                X = R*cos(THETA); Y = R*sin(THETA)
                femm.mo_selectblock(X, Y) # or you can select circuit rA rB ...
                vals_results_rotor_current.append(femm.mo_blockintegral(7)) # integrate for current
                femm.mo_clearblock()
            # the other half bar of rA
            THETA_BAR += angle_per_slot
            THETA = THETA_BAR - 2*EPS
            X = R*cos(THETA); Y = R*sin(THETA)
            femm.mo_selectblock(X, Y)
            vals_results_rotor_current.append(femm.mo_blockintegral(7)) # integrate for current
            femm.mo_clearblock()

            ################################################################
            # Also collect slot area information for loss evaluation in JMAG optimization 
            ################################################################
            if True:
                # get stator slot area for copper loss calculation
                femm.mo_groupselectblock(11)
                stator_slot_area = femm.mo_blockintegral(5) / (im.Qs/fraction) # unit: m^2 (verified by GUI operation)
                femm.mo_clearblock()

                # get rotor slot area for copper loss calculation
                femm.mo_groupselectblock(101)
                rotor_slot_area = femm.mo_blockintegral(5) / (im.Qs/fraction)
                femm.mo_clearblock()

            femm.mo_close()
            # return [-el for el in vals_results_rotor_current[self.rotor_slot_per_pole:2*self.rotor_slot_per_pole]] # 用第四象限的转子电流，因为第三象限的被切了一半，麻烦！
            # vals_results_rotor_current[self.rotor_slot_per_pole:2*self.rotor_slot_per_pole]这里用的都是第四象限的转子电流了，我们后面默认用的是第三象限的转子电流，即rA1 rB1 ...，所以要反相一下(-el)
            vals_results_rotor_current = [-el for el in vals_results_rotor_current[self.rotor_slot_per_pole:2*self.rotor_slot_per_pole]]
        # vals_results_rotor_current = self.femm_integrate_4_current(self.fraction)

        if dir_output is None:
            dir_output = self.dir_run_sweeping

        if returnData == False: # no return then write to file
            with open(dir_output + "femm_rotor_current_conditions.txt", "w") as stream:
                str_results = ''
                for el in vals_results_rotor_current:
                    stream.write("%g %g \n" % (el.real, el.imag))
            print('done. append to eddycurrent_results.txt.')
            return None
        else:
            return vals_results_rotor_current, stator_slot_area, rotor_slot_area


    def show_results_static(self, bool_plot=True):
        # Collect results from all the .ans file in the dir_run folder of FEMM.

        # recall that for static FEA, you call show_results once when half .ans files are generated from watchdog
        self.freq = 0 # needed for 

        # TODO 判断，如果文件存在，且不是空的！
        # if exists .txt file, then load it
        missed_ans_file_list = []
        if os.path.exists(self.dir_run + "static_results.txt"):
            data = np.loadtxt(self.dir_run + "static_results.txt", unpack=True, usecols=(0,1,2,3))

            # use dict to eliminate duplicates
            results_dict = {}
            for i in range(len(data[0])):
                results_dict[data[0][i]] = (data[1][i], data[2][i], data[3][i]) 
            keys_without_duplicates = list(OrderedDict.fromkeys(data[0])) # remove duplicated item because it is a dict now!
            keys_without_duplicates.sort()

            # check for missed .ans files
            if len( np.arange(0, 180, self.deg_per_step) ) == len(keys_without_duplicates):
                pass
            else:
                for self.rotor_position_in_deg in np.arange(0, 180, self.deg_per_step):
                    flag_missed = True
                    for key in keys_without_duplicates:
                        if int('%04d'%(10*self.rotor_position_in_deg)) == key:
                            flag_missed = False
                            break
                    if flag_missed == True:
                        missed_ans_file_list.append( self.get_output_file_name(booL_dir=False) + '%04d'%(10*self.rotor_position_in_deg) + '.ans' )
            print('missed:', missed_ans_file_list)

            # typical print gives: 5 1813 1795 1799.0
            # print len(missed_ans_file_list), len(data[0]), len(keys_without_duplicates), keys_without_duplicates[-1]
            # quit()

            # write data without duplicates to file
            with open(self.dir_run + "static_results_no_duplicates.txt", 'w') as f:
                for key in keys_without_duplicates:
                    f.writelines('%g %g %g %g\n' % (key, results_dict[key][0], results_dict[key][1], results_dict[key][2]))
                print('[FEMM.show_results_static] the last key is', max(keys_without_duplicates), '[begin from 0]. the length of keys is', len(keys_without_duplicates))

            data = np.loadtxt(self.dir_run + "static_results_no_duplicates.txt", unpack=True, usecols=(0,1,2,3))

            last_index = len(data[0])
        else:
            last_index = 0

        ans_file_list = os.listdir(self.dir_run)
        ans_file_list = [f for f in ans_file_list if '.ans' in f]

        # last_index == 0 则说明是第一次运行后处理
        if last_index > 0:
            if len(ans_file_list) <= last_index:
                if bool_plot == True:
                    self.plot_results(data)
                return data
            else:
                print('There are new .ans files. Now append them')

        # iter ans_file_list and missed_ans_file_list, and write to .txt 
        femm.openfemm(True) # bHide
        print('there are total %d .ans files'%(len(ans_file_list)))
        print('I am going to append the rest %d ones.'%(len(ans_file_list) - last_index))
        for ind, f in enumerate( ans_file_list[last_index:] + missed_ans_file_list ):
            if ind >= len(ans_file_list[last_index:]):
                print('...open missed .ans files')
                if os.path.exists(self.dir_run + f) == False:
                    print('run mi_analyze for %s' % (f))
                    femm.opendocument(self.dir_run + f[:-4] + '.fem')
                    femm.mi_analyze(1)
                else:
                    femm.opendocument(self.dir_run + f[:-4] + '.fem')                    
            else:
                print(last_index + ind, end=' ')
                femm.opendocument(self.dir_run + f[:-4] + '.fem')

            # load solution (if corrupted, re-run)
            try:
                femm.mi_loadsolution()
            except Exception as e:
                logger = logging.getLogger(__name__)
                logger.error('The .ans file to this .fem file is corrupted. re-run the .fem file %s'%(f), exc_info=True)
                femm.opendocument(self.dir_run + f[:-4] + '.fem')
                femm.mi_analyze(1)
                femm.mi_loadsolution()

            # get the physical amounts on the rotor
            try:
                femm.mo_groupselectblock(100)
                femm.mo_groupselectblock(101)
                Fx = femm.mo_blockintegral(18) #-- 18 x (or r) part of steady-state weighted stress tensor force
                Fy = femm.mo_blockintegral(19) #--19 y (or z) part of steady-state weighted stress tensor force
                torque = femm.mo_blockintegral(22) #-- 22 = Steady-state weighted stress tensor torque
                femm.mo_clearblock()

                # Air Gap Boundary for Rotor Motion #5
                # gap_torque = femm.mo_gapintegral("AGB4RM", 0)
                # gap_force = femm.mo_gapintegral("AGB4RM", 1)
                # print gap_force, gap_torque, torque, Fx, Fy

                # write results to a data file
                with open(self.dir_run + "static_results.txt", "a") as stream:
                    stream.write("%s %g %g %g\n" % ( f[-8:-4], torque, Fx, Fy ))
            except Exception as e:
                logger = logging.getLogger(__name__)
                logger.error('Encounter error while post-processing (integrating, etc.).', exc_info=True)
                raise e

            # avoid run out of RAM when there are a thousand of ans files loaded into femm...
            # if ind % 10 == 0:
            #     femm.closefemm() 
            #     femm.openfemm(True)
            femm.mo_close()  # use mo_ to close .ans file
            femm.mi_close()  # use mi_ to close .fem file

        print('done. append to static_results.txt.')
        femm.closefemm()

        try:
            data
        except:
            print('call this method again to plot...')
            return None
        if bool_plot == True:
            self.plot_results(data)
        return data

    def write_physical_data(self, results_list):
        with open(self.dir_run + "static_results.txt", "a") as f:
            results_rotor = ''
            for ind, row in enumerate(results_list):
                results_rotor += "%s %g %g %g\n" \
                    % ( row[0][-8:-4], row[1], row[2], row[3] )
            f.write(results_rotor)        

    def plot_results(self, data):
        from pylab import subplots, legend, show
        try: 
            self.fig
        except:
            fig, axes = subplots(3, 1, sharex=True)
            self.fig = fig
            self.axes = axes
        else:
            fig = self.fig
            axes = self.axes

        if self.flag_eddycurrent_solver:
            ax = axes[0]; ax.plot(data[0], data[1],label='torque'); ax.legend(); ax.grid()
            ax = axes[1]; ax.plot(data[0], data[2],label='Fx'); ax.legend(); ax.grid()
            ax = axes[2]; ax.plot(data[0], data[3],label='Fy'); ax.legend(); ax.grid()

        if self.flag_static_solver:
            ax = axes[0]; ax.plot(data[0]*0.1, data[1],label='torque'); ax.legend(); ax.grid()
            ax = axes[1]; ax.plot(data[0]*0.1, data[2],label='Fx'); ax.legend(); ax.grid()
            ax = axes[2]; ax.plot(data[0]*0.1, data[3],label='Fy'); ax.legend(); ax.grid()

    def run_frequency_sweeping(self, freq_range, fraction=2):

        if self.has_results(dir_run=self.dir_run_sweeping):
            return

        self.flag_static_solver = False
        self.flag_eddycurrent_solver = True
        self.fraction = fraction

        for f in os.listdir(self.dir_run_sweeping):
            os.remove(self.dir_run_sweeping + f)

        femm.openfemm(True) # bHide # False for debug
        femm.newdocument(0) # magnetic
        self.freq_range = freq_range
        self.freq = freq_range[0]
        # Alternatively, the length of the machine could be scaled by the number of segments to make this correction automatically.
        self.stack_length = self.im.stack_length * fraction
        self.probdef()

        # is coarse mesh causing biased rotor current?
        # femm.smartmesh(True)
        # self.bool_automesh = True

        self.add_material()
        # self.draw_model(fraction=fraction)
        self.vangogh.draw_model(fraction=fraction)
        self.add_block_labels(fraction=fraction)

        # # debug here
        # femm.mi_maximize()
        # femm.mi_zoomnatural()
        # return

        list_ans_file = []
        for freq in freq_range:
            self.freq = freq
            temp = self.get_output_file_name(booL_dir=False)
            list_ans_file.append(temp+'.ans')
            self.output_file_name = self.dir_run_sweeping + temp            
            print(temp)
            if os.path.exists(self.output_file_name + '.ans'):
                continue
            self.probdef()        
            femm.mi_saveas(self.output_file_name + '.fem')

        self.parallel_solve(dir_run=self.dir_run_sweeping, number_of_instantces=5) # subprocess will wait for cmd but not the pytho script
        self.wait(list_ans_file)

        # flux and current of circuit can be used for parameter identification
        if False:
            dict_circuits = {}
            # i = femm.mo_getprobleminfo()
            logging.getLogger().info('Sweeping: %g Hz.'%(self.freq))
            femm.mi_analyze(1) # None for inherited. 1 for a minimized window,
            femm.mi_loadsolution()

            # circuit
                # i1_re,i1_im, v1_re,v1_im, flux1_re,flux1_im = femm.mo_getcircuitproperties("dA")
                # i2_re,i2_im, v2_re,v2_im, flux2_re,flux2_im = femm.mo_getcircuitproperties("bA")
                # i3_re,i3_im, v3_re,v3_im, flux3_re,flux3_im = femm.mo_getcircuitproperties("rA")
            dict_circuits['dU'] = femm.mo_getcircuitproperties("dU")
            dict_circuits['dV'] = femm.mo_getcircuitproperties("dV")
            dict_circuits['dW'] = femm.mo_getcircuitproperties("dW")
            dict_circuits['bU'] = femm.mo_getcircuitproperties("bU")
            dict_circuits['bV'] = femm.mo_getcircuitproperties("bV")
            dict_circuits['bW'] = femm.mo_getcircuitproperties("bW")
            for i in range(self.rotor_slot_per_pole):
                circuit_name = 'r%s'%(self.rotor_phase_name_list[i])
                dict_circuits[circuit_name] = femm.mo_getcircuitproperties(circuit_name)

            # write results to a data file, multiplying by ? to get
            # the results for all ? poles of the machine. # 如果只分析了对称场，记得乘回少掉的那部分。
            with open(self.dir_run_sweeping + "results.txt", "a") as f:

                results_circuits = "[DW] %g + j%g A. %g + j%g V. %g + j%g Wb. [BW] %g + j%g A. %g + j%g V. %g + j%g Wb. [BAR] %g + j%g A. %g + j%g V. %g + j%g Wb. " \
                    % ( np.real(dict_circuits['dU'][0]), np.imag(dict_circuits['dU'][0]),
                        np.real(dict_circuits['dU'][1]), np.imag(dict_circuits['dU'][1]),
                        np.real(dict_circuits['dU'][2]), np.imag(dict_circuits['dU'][2]),
                        np.real(dict_circuits['bU'][0]), np.imag(dict_circuits['bU'][0]),
                        np.real(dict_circuits['bU'][1]), np.imag(dict_circuits['bU'][1]),
                        np.real(dict_circuits['bU'][2]), np.imag(dict_circuits['bU'][2]),
                        np.real(dict_circuits['rA'][0]), np.imag(dict_circuits['rA'][0]),
                        np.real(dict_circuits['rA'][1]), np.imag(dict_circuits['rA'][1]),
                        np.real(dict_circuits['rA'][2]), np.imag(dict_circuits['rA'][2]))

                results_rotor = "%g Hz. [Rotor] %g Nm. %g N. %g N. " \
                    % ( freq, torque, Fx, Fy)

                f.write(results_rotor + '\n')

    def wait(self, list_ans_file):
        # 删掉FEMM涡流场求解完的结果！
        flag_pso_running = True
        current_index = 0
        print('\nbegin waiting for eddy current solver...')
        while flag_pso_running:
            for fname in os.listdir(self.dir_run_sweeping):
                if fname == list_ans_file[current_index]:
                    current_index += 1
                    print(current_index)
                    if current_index == len(list_ans_file):
                        flag_pso_running = False
                        break
            sleep(1)

    def get_rotor_current_function(self, i=0):
        '''[For post-process of plotting rotor current waveform]
        
        Keyword Arguments:
            i {number} -- [which rotor bar? i=0 means rA1] (default: {0})
        
        Returns:
            [type] -- [function of rotor current]
        '''
        try:
            self.dict_rotor_current_function
        except:
            print('\n\n\nrotor current function does not exist. build it now...')
            self.dict_rotor_current_function, _ = self.read_current_conditions_from_JMAG()
        return self.dict_rotor_current_function[i]

    def keep_collecting_static_results_for_optimization(self, list_name=None, list_rotor_position_in_deg=None):
        if list_name is None:
            list_rotor_position_in_deg = np.arange(0, 180, self.deg_per_step)
            list_name = ['%04d'%(10*el) for el in list_rotor_position_in_deg] # with no suffix 

        self.freq = 0
        prefix = self.get_output_file_name(booL_dir=False)
        # print prefix + list_name[0] + '.ans'
        self.number_ans = len(list_name)
        handle_torque = open(self.dir_run + "static_results.txt", 'w')
        femm.openfemm(True)

        time_init = clock_time()
        current_index = 0
        flag_run_while = True
        # x_division_stator = 
        # y_division_stator = 
        while flag_run_while:
            for fname in os.listdir(self.dir_run):
                current_ans_file = prefix + list_name[current_index] + '.ans'
                if fname == current_ans_file:
                    print('[debug]', current_index, list_name[current_index], list_rotor_position_in_deg[current_index])
                    # initialization
                    if current_index == 0:
                        femm.opendocument(self.dir_run + current_ans_file)

                        # get slot area for copper loss calculation
                        femm.mo_groupselectblock(11) # fraction is 1 
                        self.stator_slot_area = femm.mo_blockintegral(5) / self.im.Qs # unit: m^2 (verified by GUI operation)
                        femm.mo_clearblock()

                        femm.mo_groupselectblock(101)
                        self.rotor_slot_area = femm.mo_blockintegral(5) / self.im.Qr
                        femm.mo_clearblock()

                        self.number_of_elements = femm.mo_numelements()
                        self.stator_Area_data = []
                        self.stator_xy_complex_data = []
                        self.rotor_Area_data = []
                        self.rotor_xy_complex_data = []
                        for id_element in range(1, self.number_of_elements+1):
                            _, _, _, x, y, area, group = femm.mo_getelement(id_element)
                            # consider 1/4 model for loss (this is valid if we presume the suspension two pole field is weak)
                            # Use the mesh info of the initial rotor position 
                            if group == 10: # stator iron
                                new_xy_complex = (x+1j*y) * np.exp(1j* pi/self.im.Qs) # 分割线经过槽
                                if new_xy_complex.real>0 and new_xy_complex.imag>0:
                                    self.stator_Area_data.append(area)
                                    self.stator_xy_complex_data.append(x+1j*y)
                            elif group == 100: # rotor iron
                                if y>0 and x>0:
                                    self.rotor_Area_data.append(area)
                                    self.rotor_xy_complex_data.append(x+1j*y)
                        self.stator_Bx_data = []
                        self.stator_By_data = []
                        for i in range(len(self.stator_xy_complex_data)):
                            self.stator_Bx_data.append([])
                            self.stator_By_data.append([])
                        self.rotor_Bx_data = []
                        self.rotor_By_data = []
                        for i in range(len(self.rotor_xy_complex_data)):
                            self.rotor_Bx_data.append([])
                            self.rotor_By_data.append([])
                        femm.mo_close()
                        # debug
                        from pylab import scatter, figure
                        figure()
                        for ind, complex_number in enumerate(self.stator_xy_complex_data):
                            scatter( complex_number.real, complex_number.imag)
                        figure()
                        for ind, complex_number in enumerate(self.rotor_xy_complex_data):
                            scatter( complex_number.real, complex_number.imag)
                        # show()
                        # print 'stator_xy_complex_data', len(self.stator_xy_complex_data)
                        # print 'number_of_elements', self.number_of_elements
                        # group10_area = 0.0
                        # for el in self.stator_Area_data:
                        #     group10_area += el
                        # group100_area = 0.0
                        # for el in self.rotor_Area_data:
                        #     group100_area += el
                        # print group10_area*1e-6, group100_area*1e-6
                        # quit()
                    # read field data and write torque data
                    current_rotor_position = list_rotor_position_in_deg[current_index] / 180. * pi
                    self.read_Torque_and_B_data(self.dir_run + current_ans_file,
                                                np.exp(1j*current_rotor_position),
                                                handle_torque)
                    current_index += 1
                    if current_index == self.number_ans:
                        flag_run_while = False
                        break
                    continue
            # no file matches, wait a while.
            print('Sleep', clock_time() - time_init, 's')
            if current_index != self.number_ans:
                sleep(3)

        femm.closefemm()
        handle_torque.close()
        logger.debug('FEMM with %g deg per .fem file takes %g sec.'%(self.deg_per_step, clock_time() - time_init))

        if self.bool_static_fea_loss == True:
            if False: 
                from pylab import subplots, show
                def test(Bx_data, base_freq):
                    print('There are in total', len(Bx_data), 'elements per step.')
                    print('There are in total', len(Bx_data[0]), 'steps.')

                    fig_dft, axes_dft = subplots(2,1)
                    fig, ax = subplots()
                    for id_element in range(0, len(Bx_data), 200): # typical id
                        Bx_waveform = Bx_data[id_element]
                        ax.plot(np.arange(0, 180, self.deg_per_step), Bx_waveform, label=id_element, alpha=0.3)

                        # DFT
                        samp_freq = 1. / (self.deg_per_step/180.*pi / self.im.Omega)
                        utility.basefreqDFT(Bx_waveform, samp_freq, ax_time_domain=axes_dft[0], ax_freq_domain=axes_dft[1], base_freq=base_freq)
                    ax.legend()
                test(self.stator_Bx_data, 500)
                test(self.stator_By_data, 500)
                test(self.rotor_Bx_data, 1)
                test(self.rotor_By_data, 1)
                show()

            temp = self.get_iron_loss()
            print('iron loss', temp)
            temp = self.get_copper_loss(self.stator_slot_area, self.rotor_slot_area)
            print('copper loss', temp)

            # np.savetxt(f, np.c_[self.stator_Bx_data])
            np.savetxt(self.dir_run + 'stator_Bx_data.txt', self.stator_Bx_data)
            np.savetxt(self.dir_run + 'stator_By_data.txt', self.stator_By_data)
            np.savetxt(self.dir_run + 'stator_Area_data.txt', self.stator_Area_data)
            np.savetxt(self.dir_run + 'rotor_Bx_data.txt', self.rotor_Bx_data)
            np.savetxt(self.dir_run + 'rotor_By_data.txt', self.rotor_By_data)
            np.savetxt(self.dir_run + 'rotor_Area_data.txt', self.rotor_Area_data)


    def read_Torque_and_B_data(self, ans_file, rotation_operator, handle_torque):
        # (str_rotor_position, rotation_operator):
        femm.opendocument(ans_file)

        # Physical Amount on the Rotor
        femm.mo_groupselectblock(100) # rotor iron
        femm.mo_groupselectblock(101) # rotor bars
        Fx = femm.mo_blockintegral(18) #-- 18 x (or r) part of steady-state weighted stress tensor force
        Fy = femm.mo_blockintegral(19) #--19 y (or z) part of steady-state weighted stress tensor force
        torque = femm.mo_blockintegral(22) #-- 22 = Steady-state weighted stress tensor torque
        femm.mo_clearblock()
        # write results to a data file (write to partial files to avoid compete between parallel instances)
        handle_torque.write("%s %g %g %g\n" % ( ans_file[-8:-4], torque, Fx, Fy ))    

        # stator iron (group==10)
        for ind, stator_xy_complex_number in enumerate(self.stator_xy_complex_data):
            # 1. What we need for iron loss evaluation is the B waveform at a fixed point (x,y). 
            #    For example, (x,y) is the centeroid of element in stator tooth.
            Bx, By = femm.mo_getb( stator_xy_complex_number.real,
                                   stator_xy_complex_number.imag)
            self.stator_Bx_data[ind].append(Bx)
            self.stator_By_data[ind].append(By)

        # rotor iron (group==100)
        for ind, rotor_xy_complex_number in enumerate(self.rotor_xy_complex_data):
            # 2. The element at (x,y) is no longer the same element from last rotor position.
            #    To find the exact element from last rotor position,
            #    we rotate the (x,y) forward as we rotate the model (rotor), get the B value there: (x,y)*rotation_operator, and correct the (Bx,By)/rotation_operator
            new_xy_complex = rotor_xy_complex_number * rotation_operator
            Bx, By = femm.mo_getb( new_xy_complex.real, 
                                   new_xy_complex.imag )
            new_BxBy_complex = (Bx + 1j*By) / rotation_operator
            self.rotor_Bx_data[ind].append(new_BxBy_complex.real)
            self.rotor_By_data[ind].append(new_BxBy_complex.imag)

        femm.mo_close()

    def get_copper_loss_Bolognani(self, stator_slot_area, rotor_slot_area, STATOR_SLOT_FILL_FACTOR=0.5, ROTOR_SLOT_FILL_FACTOR=1.0, TEMPERATURE_OF_COIL=75, total_CurrentAmp=None):
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

        air_gap_length_delta     = self.im.design_parameters[0]*1e-3 # m

        # http://127.0.0.1:4000/tech/ECCE-2019-Documentation/

        ################################################################
        # Stator Copper Loss 
        ################################################################
        tooth_width_w_t          = self.im.design_parameters[1]*1e-3 # m
        Area_S_slot              = stator_slot_area
        area_copper_S_Cu         = STATOR_SLOT_FILL_FACTOR * Area_S_slot
        a                        = self.im.wily.number_parallel_branch
        zQ                       = self.im.DriveW_zQ
        coil_pitch_yq            = self.im.wily.coil_pitch
        Q                        = self.im.Qs
        # the_radius_m             = 1e-3*(0.5*(self.im.Radius_OuterRotor + self.im.Length_AirGap + self.im.Radius_InnerStatorYoke))
        stack_length_m           = 1e-3*self.im.stack_length
        # number_of_phase          = 3
        # Ns                       = zQ * self.im.Qs / (2 * number_of_phase * a) # 3 phase winding
        # density_of_copper        = 8960 
        # k_R                      = 1 # AC resistance factor
        current_rms_value        = total_CurrentAmp / 1.4142135623730951 # for one phase
        # Area_conductor_Sc        = Area_S_slot * STATOR_SLOT_FILL_FACTOR / zQ

        Js = (current_rms_value/a) * zQ / area_copper_S_Cu # 逆变器电流current_rms_value在流入电机时，
                                                           # 受到并联支路分流，(current_rms_value/a)才是实际导体中流动的电流值，
                                                           # 这样的电流在一个槽内有zQ个，所以Islot=(current_rms_value/a) * zQ
                                                           # 槽电流除以槽内铜的面积，就是电流密度

        stator_inner_diameter_D = 2*(air_gap_length_delta + self.im.Radius_OuterRotor*1e-3)
        slot_height_h_t = 0.5*(self.im.stator_yoke_diameter_Dsyi - stator_inner_diameter_D)
        slot_pitch_pps = np.pi * (stator_inner_diameter_D + slot_height_h_t) / Q
        kov = 1.8 # \in [1.6, 2.0]
        end_winding_length_Lew = np.pi*0.5 * (slot_pitch_pps + tooth_width_w_t) + slot_pitch_pps*kov * (coil_pitch_yq - 1)

        Vol_Cu = area_copper_S_Cu * (stack_length_m + end_winding_length_Lew) * Q
        stator_copper_loss = rho_Copper * Vol_Cu * Js**2

        Vol_Cu_along_stack = area_copper_S_Cu * (end_winding_length_Lew) * Q
        stator_copper_loss_along_stack = rho_Copper * Vol_Cu_along_stack * Js**2

        print('Stator current [Arms]:', current_rms_value, 'Js:', Js)


        ################################################################
        # Rotor Copper Loss
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

        return stator_copper_loss, rotor_copper_loss, stator_copper_loss_along_stack, rotor_copper_loss_along_stack, Js, Jr, Vol_Cu

    def get_copper_loss_pyrhonen(self, stator_slot_area, rotor_slot_area, STATOR_SLOT_FILL_FACTOR=0.5, ROTOR_SLOT_FILL_FACTOR=1.0, TEMPERATURE_OF_COIL=75, total_CurrentAmp=None):

        # 这里有一个问题：就是用的几何参数是im_template的值，而不是来自 self.im.design_parameter
        # 这里有一个问题：就是用的几何参数是im_template的值，而不是来自 self.im.design_parameter
        # 这里有一个问题：就是用的几何参数是im_template的值，而不是来自 self.im.design_parameter

        # make sure these two values 
        # space_factor_kCu = SLOT_FILL_FACTOR in Pyrhonen09 design
        # space_factor_kAl = 1 in Pyrhonen09 design
        # as well as temperature
        #   TEMPERATURE_OF_COIL: temperature increase, degrees C
        #   Here, by conductor it means the parallel branch (=1) is considered, and number_of_coil_per_slot = 8 (and will always be 8, see example below).
        #   Just for example, if parallel branch is 2, number_of_coil_per_slot will still be 8, even though we cut the motor in half and count there will be 16 wire cross-sections,
        #   because the reduction in resistance due to parallel branch will be considered into the variable resistance_per_coil.
        #   Remember that the number_of_coil_per_slot (in series) is limited by the high back EMF rather than slot area.
        im = self.im
        rho_Copper = (3.76*TEMPERATURE_OF_COIL+873)*1e-9/55. # resistivity


        # Stator Copper Loss 
        Area_slot                = stator_slot_area
        a                        = im.wily.number_parallel_branch
        zQ                       = im.DriveW_zQ
        coil_pitch_by_slot_count = im.wily.coil_pitch
        Q                        = im.Qs
        the_radius_m             = 1e-3*(0.5*(im.Radius_OuterRotor + im.Length_AirGap + im.Radius_InnerStatorYoke))
        stack_length_m           = 1e-3*im.stack_length
        number_of_phase          = 3
        Ns                       = zQ * im.Qs / (2 * number_of_phase * a) # 3 phase winding
        density_of_copper        = 8960 
        k_R                      = 1 # AC resistance factor
        current_rms_value        = total_CurrentAmp / 1.4142135623730951 # for one phase

        Area_conductor_Sc        = Area_slot * STATOR_SLOT_FILL_FACTOR / zQ
        Js = current_rms_value / (a * Area_conductor_Sc)

        coil_span_W = coil_pitch_by_slot_count / Q * (the_radius_m * 2*pi)
        lav = 2*stack_length_m + 2.4*coil_span_W + 0.1
        mass_Cu             = density_of_copper * Ns * lav * Area_conductor_Sc
        mass_Cu_along_stack = density_of_copper * Ns * (2*stack_length_m) * Area_conductor_Sc

        stator_copper_loss             = number_of_phase * rho_Copper / density_of_copper * k_R * Js**2 * mass_Cu
        stator_copper_loss_along_stack = number_of_phase * rho_Copper / density_of_copper * k_R * Js**2 * mass_Cu_along_stack # this part of loss will scale as the depth of the 2D model is scaled


        # Rotor Copper Loss
        Area_slot                = rotor_slot_area
        a                        = 1
        zQ                       = 1
        coil_pitch_by_slot_count = im.Qr/im.DriveW_poles
        Q                        = im.Qr
        the_radius_m             = 1e-3*(im.Radius_OuterRotor - im.Radius_of_RotorSlot)
        stack_length_m           = 1e-3*im.stack_length
        number_of_phase          = im.Qr/im.DriveW_poles
        Ns                       = zQ * im.Qr / (2 * number_of_phase * a) # turns in series = 2 for 4 pole winding; = 1 for 2 pole winding
        density_of_copper        = 8960 
        k_R                      = 1 # AC resistance factor
        current_rms_value        = None

        Area_conductor_Sc        = Area_slot * ROTOR_SLOT_FILL_FACTOR / zQ        
        # print('list_rotor_current_amp', self.list_rotor_current_amp) # self.list_rotor_current_amp is defined in population.py
        rotor_copper_loss             = 0.0
        rotor_copper_loss_along_stack = 0.0
        # sum_rotor_current_density     = 0.0
        list_Jr = []
        for amp in self.list_rotor_current_amp:
            current_rms_value              = amp / 1.4142135623730951
            list_Jr.append(current_rms_value / (a * Area_conductor_Sc))
        Jr = sum(list_Jr) / len(list_Jr) # take average for Jr
        # print('Jr=%g Arms/m^2'%(Jr))

        coil_span_W = coil_pitch_by_slot_count / Q * (the_radius_m * 2*pi)
        lav = 2*stack_length_m + 2.4*coil_span_W + 0.1
        mass_Cu             = density_of_copper * Ns * lav * Area_conductor_Sc
        mass_Cu_along_stack = density_of_copper * Ns * (2*stack_length_m) * Area_conductor_Sc

        rotor_copper_loss             = number_of_phase * rho_Copper / density_of_copper * k_R * Jr**2 * mass_Cu
        rotor_copper_loss_along_stack = number_of_phase * rho_Copper / density_of_copper * k_R * Jr**2 * mass_Cu_along_stack # this part of loss will scale as the depth of the 2D model is scaled

        # print('stator slot area', stator_slot_area, 'm^2')
        # print('rotor slot area', rotor_slot_area, 'm^2')
        return stator_copper_loss, rotor_copper_loss, stator_copper_loss_along_stack, rotor_copper_loss_along_stack, Js, Jr

    def get_iron_loss(self, MAX_FREQUENCY=50e3, SLOT_FILL_FACTOR=0.5, TEMPERATURE_OF_COIL=75):
        # http://www.femm.info/wiki/SPMLoss
        # % Now, total core loss can be computed in one fell swoop...
        im = self.im
        samp_freq = 1. / (self.deg_per_step/180.*pi / self.im.Omega)
        print('Sampling frequency is', samp_freq, 'Hz', '<=> deg_per_step is', self.deg_per_step)

        # Iron Loss
        # % Dividing the result by cs corrects for the lamination stacking factor
        if 'M19' in self.im.fea_config_dict['Steel'] or 'M15' in self.im.fea_config_dict['Steel']:
            # M-19 Steel
            ce = 0.530 # % Eddy current coefficient in (Watt/(meter^3 * T^2 * Hz^2)
            ch = 143.  # % Hysteresis coefficient in (Watts/(meter^3 * T^2 * Hz)
            cs = 0.95  # % Lamination stacking factor (nondimensional)
        elif self.im.fea_config_dict['Steel'] == 'Arnon5':
            # %Arnon7
            ce = 0.07324 # % Eddy current coefficient in (Watt/(meter^3 * T^2 * Hz^2)
            ch = 187.6   # % Hysteresis coefficient in (Watts/(meter^3 * T^2 * Hz)
            cs = 0.96    # % Lamination stacking factor (nondimensional)

        # % Get parameters for proximity effect loss computation for phase windings
        # AWG     = 25       # % Magnet wire gauge used in winding
        # dwire   = 0.324861*0.0254*exp(-0.115942*AWG)   # % wire diameter in meters as a function of AWG
        # owire   = (58.*1e6) / (1+TEMPERATURE_OF_COIL*0.004) # % conductivity of the wire in S/m at prescribed deltaT
        # cePhase = SLOT_FILL_FACTOR * (pi**2/8.) * dwire**2 *owire

        # dff = MyLowestHarmonic*thisFrequency*w.*(w<(ns/2));  
        try:
            NFFT = self.number_ans
        except:
            NFFT = len(np.arange(0, 180, self.deg_per_step))

        dft_freq = 0.5*samp_freq*np.linspace(0,1,NFFT/2+1)
        print(dft_freq[1])
        dft_freq = dft_freq + im.slip_freq_breakdown_torque
        print(dft_freq[1])

        stator_eddycurrent_loss = 0.0
        stator_hysteresis_loss = 0.0
        rotor_eddycurrent_loss = 0.0
        rotor_hysteresis_loss = 0.0
        stator_volume = 0.0
        rotor_volume = 0.0
        # prox_loss = 0.0

        # DEBUG:
        # for id_element in range(self.number_ans):
        # for id_element in range(self.number_of_elements):

        # Element-wise calculation
        try:
            self.stator_Area_data
        except:
            # def what_loadtxt_should_be_doing(fname):
            #     data = []
            #     with open(fname, 'r') as f:
            #         read_iterator = csv_reader(f, delimiter =' ')
            #         for ind, row in enumerate(self.whole_row_reader(read_iterator)):
            #             print row
            #             data.append([float(el) for el in row])
            #     return data
            # self.stator_Bx_data   = what_loadtxt_should_be_doing(self.dir_run + 'stator_Bx_data.txt')
            # self.stator_By_data   = what_loadtxt_should_be_doing(self.dir_run + 'stator_By_data.txt')
            # self.stator_Area_data = what_loadtxt_should_be_doing(self.dir_run + 'stator_Area_data.txt')
            # self.rotor_Bx_data    = what_loadtxt_should_be_doing(self.dir_run + 'rotor_Bx_data.txt')
            # self.rotor_By_data    = what_loadtxt_should_be_doing(self.dir_run + 'rotor_By_data.txt')
            # self.rotor_Area_data  = what_loadtxt_should_be_doing(self.dir_run + 'rotor_Area_data.txt')

            self.stator_Bx_data   = np.loadtxt(self.dir_run + 'stator_Bx_data.txt'  )
            self.stator_By_data   = np.loadtxt(self.dir_run + 'stator_By_data.txt'  )
            self.stator_Area_data = np.loadtxt(self.dir_run + 'stator_Area_data.txt')
            self.rotor_Bx_data    = np.loadtxt(self.dir_run + 'rotor_Bx_data.txt'   )
            self.rotor_By_data    = np.loadtxt(self.dir_run + 'rotor_By_data.txt'   )
            self.rotor_Area_data  = np.loadtxt(self.dir_run + 'rotor_Area_data.txt' )
    
        # print np.shape(self.stator_Bx_data)
        # print np.shape(self.stator_Bx_data[0])
        # quit()

        from pylab import subplots, show
        if False:
            def test(Bx_data, base_freq):
                print('There are in total', len(Bx_data), 'elements per step.')
                print('There are in total', len(Bx_data[0]), 'steps.')

                fig_dft, axes_dft = subplots(2,1)
                fig, ax = subplots()
                for id_element in range(0, len(Bx_data), 200): # typical id
                    Bx_waveform = Bx_data[id_element]
                    ax.plot(np.arange(0, 180, self.deg_per_step), Bx_waveform, label=id_element, alpha=0.3)

                    # DFT
                    samp_freq = 1. / (self.deg_per_step/180.*pi / self.im.Omega)
                    utility.basefreqDFT(Bx_waveform, samp_freq, ax_time_domain=axes_dft[0], ax_freq_domain=axes_dft[1], base_freq=base_freq)
                ax.legend()
            test(self.stator_Bx_data, 500)
            test(self.stator_By_data, 500)
            test(self.rotor_Bx_data, 1)
            test(self.rotor_By_data, 1)
            show()


        # remove noises in frequency domain to estimate correct power spectrum
        # https://dsp.stackexchange.com/questions/9054/removing-noise-from-audio-using-fourier-transform-in-matlab
        # https://dsp.stackexchange.com/questions/6220/why-is-it-a-bad-idea-to-filter-by-zeroing-out-fft-bins
        # https://www.mathworks.com/help/matlab/math/fourier-transforms.html

        global threshold 
        def remove_noises(bxfft, threshold_tunner=0.3): # square window in freq domain
            global threshold 
            threshold = threshold_tunner * np.mean(bxfft)
            noises_places = np.where(bxfft<threshold, 0, 1)
            bxfft *= noises_places
            # print 'threshold', threshold
            return bxfft

        # Stator iron loss
        stator_eddycurrent_loss_dft = np.zeros(len(dft_freq))
        stator_hysteresis_loss_dft = np.zeros(len(dft_freq))
        for id_element, area_element in enumerate(self.stator_Area_data):             
            stator_Bx_waveform = self.stator_Bx_data[id_element]
            stator_By_waveform = self.stator_By_data[id_element]
            bxfft = utility.singleSidedDFT(stator_Bx_waveform, samp_freq)
            byfft = utility.singleSidedDFT(stator_By_waveform, samp_freq)

            # # test remove noises in spectrum
            # index_enough = -1
            # fig_dft, axes_dft = subplots(4,1, sharex=True)
            # axes_dft[3].plot(dft_freq[:index_enough], bxfft[:index_enough], '>',alpha=0.4)

            bxfft = remove_noises(bxfft)
            byfft = remove_noises(byfft)

            bsq = (bxfft*bxfft) + (byfft*byfft) # the sqaure of the amplitude of each harmonic component of flux density
            volume_element = (area_element*1e-6) * (im.stack_length*1e-3) # # Compute the volume of each element in units of meter^3            

            stator_eddycurrent_loss_dft += (ce*dft_freq**2 * volume_element/cs ) * bsq
            stator_hysteresis_loss_dft  += (ch*dft_freq    * volume_element/cs ) * bsq
            stator_volume += volume_element

            # # test remove noises in spectrum
            # axes_dft[0].plot(dft_freq[:index_enough], stator_eddycurrent_loss_dft[:index_enough], '+',alpha=0.4)
            # axes_dft[0].set_xlabel('Frequency [Hz]')
            # axes_dft[0].set_ylabel('Stator\nEddy\nCurrent\nLoss [W]')
            # axes_dft[1].plot(dft_freq[:index_enough], stator_hysteresis_loss_dft[:index_enough], '^',alpha=0.4)
            # axes_dft[1].set_xlabel('Frequency [Hz]')
            # axes_dft[1].set_ylabel('Stator\nHysteresis\nLoss [W]')
            # axes_dft[2].plot(dft_freq[:index_enough], bxfft[:index_enough], '>',alpha=0.4)
            # axes_dft[2].plot(dft_freq[:index_enough], np.ones(len(dft_freq[:index_enough]))*threshold)
            # axes_dft[3].plot(dft_freq[:index_enough], np.ones(len(dft_freq[:index_enough]))*threshold)
            # show()

        # find the index of dft_freq that corresponds to 50e3 Hz, because higher results are contaminated by noises in fourier analysis
        index_enough = next(ind for ind, el in enumerate(dft_freq) if el > MAX_FREQUENCY) # 50e3 Hz is enough 10 times base frequency 500 Hz
        stator_eddycurrent_loss = sum( stator_eddycurrent_loss_dft[:index_enough] ) 
        stator_hysteresis_loss  = sum( stator_hysteresis_loss_dft[:index_enough] )

        fig_dft, axes_dft = subplots(4,1, sharex=True)
        axes_dft[0].plot(dft_freq[:index_enough], stator_eddycurrent_loss_dft[:index_enough], '+',alpha=0.4)
        axes_dft[0].set_xlabel('Frequency [Hz]')
        axes_dft[0].set_ylabel('Stator\nEddy\nCurrent\nLoss [W]')
        axes_dft[1].plot(dft_freq[:index_enough], stator_hysteresis_loss_dft[:index_enough], '^',alpha=0.4)
        axes_dft[1].set_xlabel('Frequency [Hz]')
        axes_dft[1].set_ylabel('Stator\nHysteresis\nLoss [W]')

        # Rotor iron loss
        rotor_eddycurrent_loss_dft = np.zeros(len(dft_freq))
        rotor_hysteresis_loss_dft = np.zeros(len(dft_freq))
        for id_element, area_element in enumerate(self.rotor_Area_data):
            rotor_Bx_waveform = self.rotor_Bx_data[id_element]
            rotor_By_waveform = self.rotor_By_data[id_element]
            bxfft = utility.singleSidedDFT(rotor_Bx_waveform, samp_freq)
            byfft = utility.singleSidedDFT(rotor_By_waveform, samp_freq)

            bxfft = remove_noises(bxfft)
            byfft = remove_noises(byfft)

            bsq = (bxfft*bxfft) + (byfft*byfft) # the sqaure of the amplitude of each harmonic component of flux density
            volume_element = (area_element*1e-6) * (im.stack_length*1e-3) # # Compute the volume of each element in units of meter^3

            rotor_eddycurrent_loss_dft += (ce*dft_freq**2 * volume_element/cs ) * bsq
            rotor_hysteresis_loss_dft  += (ch*dft_freq    * volume_element/cs ) * bsq
            rotor_volume += volume_element
        rotor_eddycurrent_loss = sum( rotor_eddycurrent_loss_dft[:index_enough] ) 
        rotor_hysteresis_loss  = sum( rotor_hysteresis_loss_dft[:index_enough]  )

        axes_dft[2].plot(dft_freq[:index_enough], rotor_eddycurrent_loss_dft[:index_enough], '+',alpha=0.4)
        axes_dft[2].set_xlabel('Frequency [Hz]')
        axes_dft[2].set_ylabel('\nRotor\nEddy Current\nLoss [W]')
        axes_dft[3].plot(dft_freq[:index_enough], rotor_hysteresis_loss_dft[:index_enough], '^',alpha=0.4)
        axes_dft[3].set_xlabel('Frequency [Hz]')
        axes_dft[3].set_ylabel('\nRotor\nHysteresis\nLoss [W]')


        # # Stator iron loss
        # for id_element, area_element in enumerate(self.stator_Area_data):             
        #     stator_Bx_waveform = self.stator_Bx_data[id_element]
        #     stator_By_waveform = self.stator_By_data[id_element]
        #     bxfft = utility.singleSidedDFT(stator_Bx_waveform, samp_freq)
        #     byfft = utility.singleSidedDFT(stator_By_waveform, samp_freq)
        #     bsq = (bxfft*bxfft) + (byfft*byfft) # the sqaure of the amplitude of each harmonic component of flux density
        #     volume_element = (area_element*1e-6) * (im.stack_length*1e-3) # # Compute the volume of each element in units of meter^3
        #     stator_eddycurrent_loss += ( ce*np.dot(dft_freq**2, bsq) * volume_element ) / cs
        #     stator_hysteresis_loss  += ( ch*np.dot(dft_freq, bsq) * volume_element ) / cs
        #     stator_volume += volume_element

        # # Rotor iron loss
        # for id_element, area_element in enumerate(self.rotor_Area_data):
        #     rotor_Bx_waveform = self.rotor_Bx_data[id_element]
        #     rotor_By_waveform = self.rotor_By_data[id_element]
        #     bxfft = utility.singleSidedDFT(rotor_Bx_waveform, samp_freq)
        #     byfft = utility.singleSidedDFT(rotor_By_waveform, samp_freq)
        #     bsq = (bxfft*bxfft) + (byfft*byfft) # the sqaure of the amplitude of each harmonic component of flux density
        #     volume_element = (area_element*1e-6) * (im.stack_length*1e-3) # # Compute the volume of each element in units of meter^3
        #     rotor_eddycurrent_loss += ( ce*np.dot(dft_freq**2, bsq) * volume_element ) / cs
        #     rotor_hysteresis_loss  += ( ch*np.dot(dft_freq, bsq) * volume_element ) / cs
        #     rotor_volume += volume_element

            # Copper Loss - Stator Winding Proximity Effect (Rotor side is neglected because the slip frequency is low)
            # this should be done with g==2, i.e., field data on coil area
            # % and prox losses can be totalled up in a similar way as iron loss
            # prox_loss += np.dot(cePhase * dft_freq**2, bsq) * volume_element

        # we use only 1/4 model for the loss caculation
        number_of_fraction = 4
        return ( number_of_fraction*stator_eddycurrent_loss, 
                 number_of_fraction*stator_hysteresis_loss, 
                 number_of_fraction*rotor_eddycurrent_loss, 
                 number_of_fraction*rotor_hysteresis_loss, 
                 number_of_fraction*stator_volume, 
                 number_of_fraction*rotor_volume )


    # this is for JMAG calling to seek for breakdown torque and slip
    def greedy_search_for_breakdown_slip(self, dir_femm_temp, study_name, fraction=2, number_of_instantces=5, bool_run_in_JMAG_Script_Editor=False):
        if not os.path.isdir(dir_femm_temp):
            os.makedirs(dir_femm_temp)

        # wait_greedy_search detect femm_found.csv files to confirm the complete state of this greedy_search
        for suffix in ['.ans', '.csv', '.fem']:
            fname = dir_femm_temp+'femm_found'+suffix
            if os.path.exists(fname):
                os.remove(fname)

        self.flag_static_solver = False
        self.flag_eddycurrent_solver = True
        self.fraction = fraction

        femm.openfemm(True) # bHide # False for debug
        femm.newdocument(0) # magnetic
        
        self.freq = 1
        # Alternatively, the length of the machine could be scaled by the number of segments to make this correction automatically. -David Meeker
        self.stack_length = self.im.stack_length * fraction
        self.probdef()

        self.add_material()
        # print 'Call VanGogh to draw a FEMM model of fraction %d' % (fraction)
        self.vangogh.draw_model(fraction=fraction)
        self.add_block_labels(fraction=fraction)
        # print dir_femm_temp + 'femm_temp.fem'
        # print dir_femm_temp + 'femm_temp.fem'
        # print dir_femm_temp + 'femm_temp.fem'
        femm.mi_saveas(dir_femm_temp + 'femm_temp.fem')
        # quit()
        femm.mi_close()
        femm.closefemm()

        # save for wait_greedy_search()
        self.dir_femm_temp = dir_femm_temp
        self.study_name = study_name


        # 在搜索的过程中，显然是需要知道结果的，
        # 这说明，不仅仅要跑五个instance，
        # 而且还需要有一个总管，负责总调。
        if bool_run_in_JMAG_Script_Editor: # run inside JMAG then a wrapper batch file is required to work properly.
            with open('temp2.bat', 'w') as f:
                f.write('python "%smethod_parasolve_greedy_search.py" "%s" %d %.16f' % (self.dir_codes, 
                                                                                        dir_femm_temp, 
                                                                                        number_of_instantces,
                                                                                        self.stack_length))
            os.startfile('temp2.bat')
            # os.remove('temp2.bat')
        else:
            print('\n' + '-'*20, self.study_name)
            proc = subprocess.Popen([sys.executable, 'parasolve_greedy_search_manager.py', 
                                     str(number_of_instantces), self.dir_femm_temp, str(self.stack_length)], bufsize=-1)
            # proc.wait() # don't wait on femm solver, and let jmag plot the model and get ready for the breakdownd slip info.


    def wait_greedy_search(self, tic):
        while True:
            fname = self.dir_femm_temp + 'femm_found.csv'
            if os.path.exists(fname):
                with open(fname, 'r') as f:
                    data = f.readlines()
                    freq = float(data[0][:-1])
                    torque = float(data[1][:-1])

                femm.openfemm(True)
                vals_results_rotor_current, stator_slot_area, rotor_slot_area = \
                    self.femm_integrate_4_current(fname[:-4]+'.ans', self.fraction, dir_output=self.dir_femm_temp, returnData=True)
                femm.closefemm()
                
                new_fname = self.dir_femm_temp + self.study_name + '.csv'
                os.rename(fname, new_fname)
                with open(new_fname, 'w') as f:
                    str_results = "%g\n%g\n" % (freq, torque)
                    str_results += "%g\n%g\n" % (stator_slot_area, rotor_slot_area)
                    for el in vals_results_rotor_current:
                        str_results += "%g,%g\n" % (el.real, el.imag)
                    f.write(str_results)

                # also save for loss evaluation query # this is not good for restart, just recover your data from csv files
                self.vals_results_rotor_current = vals_results_rotor_current
                self.stator_slot_area = stator_slot_area
                self.rotor_slot_area = rotor_slot_area

                # leave the fem file for ease of reproduction
                # os.remove(fname[:-4]+'.fem')
                os.remove(fname[:-4]+'.ans')
                os.rename(fname[:-4]+'.fem', new_fname[:-4]+'.fem')
                break
            else:
                print('Wait for greedy search: sleep 1 sec...')
                sleep(1)
                # print clock_time() - tic, 's'
        toc = clock_time()
        logger = logging.getLogger(__name__)
        logger.debug('Time spent on femm frequency search is %g s.', toc-tic)
        return freq, torque, None


# def get_magnet_loss(self):
    # pass
    # % Add up eddy current losses in the magnets
    # % Magnet properties
    # RotorMagnets = Num_pole;
    # omag = 0.556*10^6;                  % conductivity of sintered NdFeB in S/m
    # % compute fft of A at the center of each element
    # Jm=fft(A)*(2/ns);
    # for k=1:RotorMagnets
    #     g3=(g==(10+k));
    #     % total volume of the magnet under consideration;
    #     vmag=v'*g3;
    #     % average current in the magnet for each harmonic
    #     Jo=(Jm*(v.*g3))/vmag;
    #     % subtract averages off of each each element in the magnet
    #     Jm = Jm - Jo*g3';
    # end
    # %     magnet_loss = (1/2)*((omag*(2*pi*w).^2)*(abs(Jm).^2)*v);
    # magnet_loss = (1/2)*((omag*(w).^2)*(abs(Jm).^2)*v);

    # total_loss = rotor_loss + stator_loss + prox_loss + SlotOhmic + magnet_loss + Air_friction_loss;
    # results = [results; thisSpeed, rotor_loss, stator_loss, magnet_loss, SlotOhmic + prox_loss, Air_friction_loss, total_loss];

