# -*- coding: utf-8 -*-

# from shapely.geometry import LineString, Point
from math import tan, pi, atan, sqrt, sin, cos, copysign, atan2, asin, acos
import utility
from pylab import plt, np

CUSTOM = 2
JMAG = 1
FEMM = 0

class VanGogh(object):
    """One VanGogh for both FEMM and JMAG"""
    def __init__(self, im, child_index):
        self.im = im
        self.child_index = child_index
        self.plot_object_list = []

    def draw_model(self, fraction=1):
        # utility functions than wrap shapely functions

        if self.child_index == JMAG: # for being consistent with obselete codes
            self.plot_sketch_shaft() # Shaft if any
            self.draw_rotor_eMach(fraction)
            # self.draw_stator_without_non_accurate_shapely(fraction)

        elif self.child_index == FEMM: # for easy selecting of objects
            utility.blockPrint()
            self.draw_stator_without_non_accurate_shapely(fraction)
            self.draw_rotor_without_non_accurate_shapely(fraction)
            # self.draw_rotor_eMach(fraction)
            utility.enablePrint()

            # self.draw_stator(fraction)
            # try:
            #     self.draw_rotor(fraction)
            # except Exception as e:
            #     raise e

        elif self.child_index == CUSTOM: 
            self.draw_rotor_without_non_accurate_shapely()
            self.draw_stator_without_non_accurate_shapely()

    def draw_rotor_without_non_accurate_shapely(self, fraction=1):
        # Shapely is very poor in accuracy, use your high school geometry knowledge to derive the coordinates!

        im = self.im

        # origin = Point(0,0)
        Stator_Sector_Angle = 2*pi/im.Qs*0.5
        Rotor_Sector_Angle = 2*pi/im.Qr*0.5

        ''' Part: Rotor '''
        if self.child_index == JMAG:
            self.init_sketch_rotorCore()
        # Draw Points as direction of CCW
        # P1
        P1 = (-im.Radius_Shaft, 0)

        # P2
        P2 = ( -im.Radius_Shaft*cos(Rotor_Sector_Angle), im.Radius_Shaft*sin(Rotor_Sector_Angle) )

        # P3
        P2_angle = P3_angle = Rotor_Sector_Angle
        P3 = ( -im.Radius_OuterRotor*cos(Rotor_Sector_Angle), im.Radius_OuterRotor*sin(Rotor_Sector_Angle) )

        # P4
        P4 = (-sqrt(im.Radius_OuterRotor**2 - (0.5*im.Width_RotorSlotOpen)**2), 0.5*im.Width_RotorSlotOpen)
        self.draw_arc(P3, P4, P3_angle - self.get_postive_angle(P4), center=(0,0))

        # P5
        P5 = (-im.Location_RotorBarCenter-sqrt(im.Radius_of_RotorSlot**2 - (0.5*im.Width_RotorSlotOpen)**2), 0.5*im.Width_RotorSlotOpen)
        self.draw_line(P4, P5)

        if im.use_drop_shape_rotor_bar == True:
            # P6
            P6 = (-im.Location_RotorBarCenter, im.Radius_of_RotorSlot)
            self.draw_arc(P6, P5, 0.5*pi - self.get_postive_angle(P5, (-im.Location_RotorBarCenter, 0)), center=(-im.Location_RotorBarCenter, 0))

            # P7
            P7 = (-im.Location_RotorBarCenter2, im.Radius_of_RotorSlot2)
            self.draw_line(P6, P7)

            # P8
            P8 = (-im.Location_RotorBarCenter2+im.Radius_of_RotorSlot2, 0)
            self.draw_arc(P8, P7, 0.5*pi, center=(-im.Location_RotorBarCenter2, 0))

        else:
            P6 = P7 = None
            P8 = (-im.Location_RotorBarCenter+im.Radius_of_RotorSlot, 0)
            self.draw_arc(P8, P5, pi - self.get_postive_angle(P5, (-im.Location_RotorBarCenter, 0)), center=(-im.Location_RotorBarCenter, 0))

        if self.child_index == FEMM:
            self.some_solver_related_operations_rotor_before_mirror_rotation(im, P6, P8) # call this before mirror_and_copyrotate

        if self.child_index == JMAG:
            self.draw_line(P8, P1)
            self.draw_arc(P2, P1, P2_angle)
            self.draw_line(P2, P3)

            self.mirror_and_copyrotate(im.Qr, im.Radius_OuterRotor, fraction,
                                        symmetry_type=2) 
            self.init_sketch_cage()

        if self.child_index == CUSTOM:
            # self.draw_line(P8, P1, ls='-.')
            self.draw_arc(P2, P1, P2_angle)
            # self.draw_line(P2, P3, ls='-.')

        # 导条
        # P_Bar
        P_Bar = (-im.Location_RotorBarCenter-im.Radius_of_RotorSlot, 0)
        self.draw_arc(P5, P_Bar, self.get_postive_angle(P5, (-im.Location_RotorBarCenter, 0)), center=(-im.Location_RotorBarCenter, 0), ls=':')

        if self.child_index == JMAG:
            self.add_line(P_Bar, P8)

            if im.use_drop_shape_rotor_bar == True:

                # draw the outline of rotor core for coil to form a region in JMAG
                self.draw_arc(P6, P5, 0.5*pi - self.get_postive_angle(P5, (-im.Location_RotorBarCenter, 0)), center=(-im.Location_RotorBarCenter, 0))
                self.draw_arc(P8, P7, 0.5*pi, center=(-im.Location_RotorBarCenter2, 0))

            else:
                raise Exception('Not implemented error')

            self.mirror_and_copyrotate(im.Qr, None, fraction,
                                        symmetry_type=2
                                        # merge=False, # bars are not connected to each other, so you don't have to specify merge=False, they will not merge anyway...
                                        # do_you_have_region_in_the_mirror=True # In short, this should be true if merge is false...
                                        )


        if self.child_index == FEMM:
            self.mirror_and_copyrotate(im.Qr, im.Radius_OuterRotor, fraction)
            self.some_solver_related_operations_fraction(im, fraction)

        if self.child_index == CUSTOM:
            # self.draw_line(P8, P_Bar, ls='-.')
            self.mirror_and_copyrotate(im.Qr, im.Radius_OuterRotor, fraction)

            # self.Pr_list = [P1,P2,P3,P4,P5,P6,P7,P8,P_Bar]
            self.Pr_list = [np.array(P) for P in [P1,P2,P3,P4,P5,P6,P7,P8]]
            for P in self.Pr_list:
                self.ax.scatter(*P, c='k', marker='None')

            # rotor objects
            self.rotor_object_list = self.plot_object_list
            self.plot_object_list = []

            for ind, P in enumerate([P1, P2, P3, P4, P5, P6, P7, P8]):
                print('rP'+str(ind+1), P) 

    def draw_stator_without_non_accurate_shapely(self, fraction=1):
        im = self.im

        # origin = Point(0,0)
        origin = [0, 0]
        Stator_Sector_Angle = 2*pi/im.Qs*0.5
        Rotor_Sector_Angle = 2*pi/im.Qr*0.5

        ''' Part: Stator '''
        if self.child_index == JMAG:
            self.init_sketch_statorCore()

        # Draw Points as direction of CCW
        # P1
        P1 = (-im.Radius_OuterRotor-im.Length_AirGap, 0)

        # P2
        Radius_InnerStator = im.Radius_OuterRotor+im.Length_AirGap
        P2_angle = im.Angle_StatorSlotOpen*0.5/180*pi
        P2_rot = (-Radius_InnerStator*cos(P2_angle), -Radius_InnerStator*sin(P2_angle))
        P2 = self.park_transform(P2_rot, Stator_Sector_Angle)
        self.draw_arc(P2, P1, self.get_postive_angle(P2))

        # P3
        P3_rot = (P2_rot[0]-im.Width_StatorTeethHeadThickness, P2_rot[1])
        P3 = self.park_transform(P3_rot, Stator_Sector_Angle)
        self.draw_line(P2, P3)

        # P4 (P4 is better to compute using intersection, error from shapely will not really cause a problem while modeling)
        c = (origin, im.Radius_OuterRotor+im.Length_AirGap+im.Width_StatorTeethHeadThickness+im.Width_StatorTeethNeck)
        l = [(0, 0.5*im.Width_StatorTeethBody), (-im.Radius_OuterStatorYoke, 0.5*im.Width_StatorTeethBody)]
        intersections = self.get_node_at_intersection(c,l)
        if intersections[0][0] < 0:
            P4 = intersections[0]
        else:
            P4 = intersections[1]
        self.draw_line(P3, P4)

        # P5
        P5 = (-sqrt(im.Radius_InnerStatorYoke**2 - (0.5*im.Width_StatorTeethBody)**2), 0.5*im.Width_StatorTeethBody)
        self.draw_line(P4, P5)

        # P6
        # k = -tan(Stator_Sector_Angle)
        # l_sector = LineString([(0,0), (-im.Radius_OuterStatorYoke, -im.Radius_OuterStatorYoke*k)])
        # P6 = self.get_node_at_intersection(c,l_sector)
        P6 = ( -im.Radius_InnerStatorYoke*cos(Stator_Sector_Angle), im.Radius_InnerStatorYoke*sin(Stator_Sector_Angle) )
        self.draw_arc(P6, P5, Stator_Sector_Angle - self.get_postive_angle(P5))


        # P7
        # c = self.create_circle(origin, im.Radius_OuterStatorYoke)
        # P7 = self.get_node_at_intersection(c,l_sector)
        P7 = [ -im.Radius_OuterStatorYoke*cos(Stator_Sector_Angle), im.Radius_OuterStatorYoke*sin(Stator_Sector_Angle) ]

        # P8
        P8 = (-im.Radius_OuterStatorYoke, 0)

        if self.child_index == JMAG:
            self.draw_line(P6, P7)
            self.draw_arc(P7, P8, Stator_Sector_Angle)
            self.draw_line(P8, P1)
            self.mirror_and_copyrotate(im.Qs, None, fraction,
                                        symmetry_type=2,  # 2: x-axis
                                        )
            self.init_sketch_coil()

        if self.child_index == CUSTOM:
            # self.draw_line(P6, P7, ls='-.')
            self.draw_arc(P7, P8, Stator_Sector_Angle)
            # self.draw_line(P8, P1, ls='-.')

        # P_Coil
        # l = LineString([(P3[0], P3[1]), (P3[0], im.Radius_OuterStatorYoke)])
        # P_Coil = self.get_node_at_intersection(l_sector, l)
        P_Coil = ( P3[0], abs(P3[0])*tan(Stator_Sector_Angle) )




        if self.child_index == JMAG: # F*ck you, Shapely for putting me through this!
            #     temp = 0.7071067811865476*(P_Coil[1] - P4[1])
            #     P4 = [                      P4[0], P4[1] + temp]
            #     P5 = [P5[0] + 0.7071067811865476*(P4[0] - P5[0]), P5[1] + temp]
            #     P6 = []
            # we use sin cos to find P6 and P7, now this suffices.
            P6[1] -= 0.01
            # Conclusion: do not use the intersection between l_sector and a circle!
            # Conclusion: do not use the intersection between l_sector and a circle!
            # Conclusion: do not use the intersection between l_sector and a circle!
            # 总之，如果overlap了，merge一下是没事的，麻烦就在Coil它不能merge，所以上层绕组和下层绕组之间就产生了很多Edge Parts。

        if self.child_index == CUSTOM:
            # self.Ps_list = [P1,P2,P3,P4,P5,P6,P7,P8,P_Coil]
            self.Ps_list = [np.array(P) for P in [P1,P2,P3,P4,P5,P6,P7,P8]]
            for P in self.Ps_list:
                self.ax.scatter(*P, c='k', marker='None')

            # stator objects
            self.stator_object_list = self.plot_object_list
            self.plot_object_list = []

            for ind, P in enumerate([P1, P2, P3, P4, P5, P6, P7, P8]):
                print('sP'+str(ind+1), P) 

        else:
            self.draw_line(P4, P_Coil)
            self.draw_line(P6, P_Coil)

        if self.child_index == JMAG:
            # draw the outline of stator core for coil to form a region in JMAG
            self.draw_line(P4, P5)
            self.draw_arc(P6, P5, Stator_Sector_Angle - self.get_postive_angle(P5))
            self.mirror_and_copyrotate(im.Qs, None, fraction,
                                        edge4ref=self.artist_list[1], #'Line.2' 
                                        # symmetry_type=2,
                                        merge=False, # two layers of windings
                                        do_you_have_region_in_the_mirror=True # In short, this should be true if merge is false...
                                        )
            # it is super wierd that use edge Line.2 as symmetry axis will lead to redundant parts imported into JMAG Designers (Extra Coil and Stator Core Parts)
            # symmetry_type=2 will not solve this problem either
            # This is actually caused by the overlap of different regions, because the precision of shapely is shit!

        # FEMM does not model coil and stator core separately        
        if self.child_index == FEMM:
            self.mirror_and_copyrotate(im.Qs, im.Radius_OuterStatorYoke, fraction)


    def draw_rotor_eMach(self, fraction):

        def utilityTangentPointsOfTwoCircles(C1, C2, r, R):
            x1 = C1[0]; y1 = C1[1]
            x2 = C2[0]; y2 = C2[1]
            gamma = -atan((y2-y1)/(x2-x1))
            distance = sqrt((x2-x1)**2+(y2-y1)**2)
            beta = asin((R-r)/distance)
            alpha = gamma - beta
            x3 = x1 + r*cos(0.5*pi - alpha)
            y3 = y1 + r*sin(0.5*pi - alpha)
            x4 = x2 + R*cos(0.5*pi - alpha)
            y4 = y2 + R*sin(0.5*pi - alpha) 
            # (x3,y3) and (x4,y4) are outer tangent points on one side.
            coord3 = [x3, y3]
            coord4 = [x4, y4]
            return coord3, coord4

        d_rs  = self.im.Location_RotorBarCenter - self.im.Location_RotorBarCenter2
        d_ro  = self.im.Length_HeadNeckRotorSlot
        w_ro  = self.im.Width_RotorSlotOpen
        w_rs1 = self.im.Radius_of_RotorSlot
        w_rs2 = self.im.Radius_of_RotorSlot2
        R_or  = self.im.Radius_OuterRotor
        R_ir  = self.im.Radius_Shaft
        Qr    = self.im.Qr

        # origin = [0,0]
        angleRotorSector = 2*pi/Qr*0.5

        ''' Part: Rotor '''
        if self.child_index == JMAG:
            self.init_sketch_rotorCore()

        # Draw Points as direction of CCW
        P1 = (-R_ir, 0)
        P2 = [-R_ir*cos(angleRotorSector), R_ir*sin(angleRotorSector)]
        P3 = [-R_or*cos(angleRotorSector), R_or*sin(angleRotorSector)]
        x, y = self.linecirc(0, 0.5*w_ro, 0, 0, R_or)
        if np.isnan(y).all(): # all element are True?
            raise Exception('Error: provided line and circle have no intersection.')
        elif x[0]<0:
            P4 = [x[0],y[0]]
        else:
            P4 = [x[1],y[1]]
        CenterP5P6 = [-(R_or - d_ro - w_rs1), 0]
        [x,y] = self.linecirc(0, 0.5*w_ro, CenterP5P6[0], CenterP5P6[1], w_rs1)
        if np.isnan(y).all():
            raise Exception('Error: provided line and circle have no intersection.')
        elif x[0] < CenterP5P6[0]:
            P5 = [x[0],y[0]]
        else:
            P5 = [x[1],y[1]]
        if self.im.use_drop_shape_rotor_bar == True:
            CenterP7P8 = [-(R_or - d_ro - w_rs1 - d_rs), 0]
            [P6, P7] = utilityTangentPointsOfTwoCircles(CenterP5P6, CenterP7P8, w_rs1, w_rs2)
            if P6[1] < 0:
                P6[1] = -1*P6[1]
                P7[1] = -1*P7[1]
            if P6[0] > P7[0]:
                P6, P7 = P7, P6
            P8 = [-(R_or - d_ro - w_rs1 - d_rs - w_rs2), 0]
        else:
            P6 = P7 = None
            P8 = [-(R_or - d_ro - 2*w_rs1), 0]

        self.listKeyPoints = [P1, P2, P3, P4, P5, P6, P7, P8];

        drawer = self
        arc21  = drawer.drawArc([0,0],P2,P1)
        line23 = drawer.drawLine(P2,P3)
        arc34  = drawer.drawArc([0,0],P3,P4)
        line45 = drawer.drawLine(P4,P5)
        line81 = drawer.drawLine(P8,P1)
        if self.im.use_drop_shape_rotor_bar == True:
            arc65  = drawer.drawArc(CenterP5P6,P6,P5)
            line67 = drawer.drawLine(P6,P7)
            arc87  = drawer.drawArc(CenterP7P8,P8,P7)
            segments = [arc21 ,line23,arc34 ,line45,arc65 ,line67,arc87 ,line81]
        else:
            arc85  = drawer.drawArc(CenterP5P6,P8,P5) # CenterP5P6 is the same as CenterP5P8
            segments = [arc21 ,line23,arc34 ,line45,arc85 ,line81]
        csToken = segments
        
        im = self.im
        if self.child_index == FEMM:
            self.some_solver_related_operations_rotor_before_mirror_rotation(im, P6, P8) # call this before mirror_and_copyrotate to asign segment size to the arc to reduce element number

        if self.child_index == JMAG:
            # JMAG needs to explictly draw the rotor slot.
            drawer.drawLine(P8, P1)
            drawer.drawArc([0,0],P2,P1)
            drawer.drawLine(P2, P3)

            self.mirror_and_copyrotate(im.Qr, im.Radius_OuterRotor, fraction,
                                        symmetry_type=2) 
            self.init_sketch_cage()

        if self.child_index == CUSTOM:
            # self.draw_line(P8, P1, ls='-.')
            self.draw_arc(P2, P1, P2_angle)
            # self.draw_line(P2, P3, ls='-.')

        # P_Bar
        P_Bar = (-im.Location_RotorBarCenter-im.Radius_of_RotorSlot, 0)
        drawer.drawArc((-im.Location_RotorBarCenter, 0), P5, P_Bar, ls=':')

        if self.child_index == JMAG:
            drawer.drawLine(P_Bar, P8)
            drawer.drawLine(P6, P7)
            # draw the outline of stator core for coil to form a region in JMAG
            drawer.drawArc((-im.Location_RotorBarCenter, 0), P6, P5)
            drawer.drawArc((-im.Location_RotorBarCenter2, 0), P8, P7)

            self.mirror_and_copyrotate(im.Qr, None, fraction,
                                        symmetry_type=2
                                        # merge=False, # bars are not connected to each other, so you don't have to specify merge=False, they will not merge anyway...
                                        # do_you_have_region_in_the_mirror=True # In short, this should be true if merge is false...
                                        )

        if self.child_index == FEMM:
            self.mirror_and_copyrotate(im.Qr, im.Radius_OuterRotor, fraction)
            self.some_solver_related_operations_fraction(im, fraction)

        if self.child_index == CUSTOM:
            # self.draw_line(P8, P_Bar, ls='-.')
            self.mirror_and_copyrotate(im.Qr, im.Radius_OuterRotor, fraction)

            # self.Pr_list = [P1,P2,P3,P4,P5,P6,P7,P8,P_Bar]
            self.Pr_list = [np.array(P) for P in [P1,P2,P3,P4,P5,P6,P7,P8]]
            for P in self.Pr_list:
                self.ax.scatter(*P, c='k', marker='None')

            # rotor objects
            self.rotor_object_list = self.plot_object_list
            self.plot_object_list = []

            for ind, P in enumerate([P1, P2, P3, P4, P5, P6, P7, P8]):
                print('rP'+str(ind+1), P)


    def drawLine(self, PA, PB):
        self.draw_line(PA, PB)

    def drawArc(self, CenterPAPB, PA, PB, **kwarg):

        if self.child_index == JMAG: # plot arc using center and two points
            self.draw_arc(CenterPAPB, PA, PB)

        elif self.child_index == FEMM or self.child_index == CUSTOM: # plot arc using arc angle and two points
            angle = abs(self.get_postive_angle(PA) - self.get_postive_angle(PB))
            self.draw_arc(PA, PB, angle, **kwarg) # angle in rad       


    @staticmethod
    def park_transform(P_rot, angle):
        return (  cos(angle)*P_rot[0] + sin(angle)*P_rot[1],
                 -sin(angle)*P_rot[0] + cos(angle)*P_rot[1] )


    # this method is found to be very low in accuracy, don't use it, just specify the center if know!
    @staticmethod
    def find_center_of_a_circle_using_2_points_and_arc_angle(p1, p2, angle):
        return find_center_and_radius_of_a_circle_using_2_points_and_arc_angle(p1, p2, angle)[0]

    @staticmethod
    def find_center_and_radius_of_a_circle_using_2_points_and_arc_angle(p1, p2, angle):
        distance = sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)
        radius   = distance*0.5 / sin(angle*0.5)
        c1 = Point(p1[0],p1[1]).buffer(radius).boundary
        c2 = Point(p2[0],p2[1]).buffer(radius).boundary
        i = c1.intersection(c2)
        # print len(i.geoms)
        # print i.geoms[0].coords[0]
        # print i.geoms[1].coords[0]

        # if len(i.geoms) > 2:
        #     print 'There are more than 2 intersections!'
        #     print p1, p2, angle
        #     print c1.coords[0], c1.coords[1], c1.coords[2]
        #     print c2.coords[0], c2.coords[1], c2.coords[2]
        #     for el in i.geoms:
        #         print '\t', el.coords[0]
        #     print 'Crazy Shapely.'
        #     # (0.013059471222689467, 0.46233200000038793)
        #     # (0.013059471222689467, 0.462332)
        #     # (-93.94908807906135, 4.615896979306289)
        #     # raise Exception('Too many intersections.')

        if i.geoms is None:
            raise Exception('There is no intersection.')

        for center in [ el.coords[0] for el in i.geoms ]: # shapely does not return intersections in a known order
            # print center
            determinant = (center[0]-p1[0]) * (p2[1]-p1[1]) - (center[1]-p1[1]) * (p2[0]-p1[0])
            if copysign(1, determinant) < 0: # CCW from p1 to p2
                return center, radius

    @staticmethod
    def create_circle(p, radius):
        return p.buffer(radius).boundary

    @staticmethod
    def get_postive_angle(p, center=(0,0)):
        # using atan loses info about the quadrant, so it is "positive"
        return atan(abs((p[1]-center[1]) / (p[0]-center[0])))

    @staticmethod
    def get_node_at_intersection(c,l):
        # this works for c and l having one or two intersections
        if c[0][0] != 0 or c[0][1] != 0:
            raise Exception('Not implemented for non-origin centered circle.')
        r = c[1]
        c = None
        x1, y1 = l[0][0], l[0][1]
        x2, y2 = l[1][0], l[1][1]
        if x1 == x2:
            raise Exception('Not implemented.')
        a = -(y2-y1)/(x2-x1)
        b = 1
        c = y1-(y2-y1)/(x2-x1)*x1
        Delta = sqrt(r**2*(a**2+b**2)-c**2)
        if Delta < 0:
            raise Exception('No intersection for given line and circle')
        x_solutions = (a*c + b*Delta)/(a**2+b**2), (a*c - b*Delta)/(a**2+b**2)
        y_solutions = (b*c - a*Delta)/(a**2+b**2), (b*c + a*Delta)/(a**2+b**2)
        return (x_solutions[0], y_solutions[0]), (x_solutions[1], y_solutions[1])

    @staticmethod
    def linecirc(slope,intercpt,centerx,centery,r):
        # Input data are in alpha-beta frame:
        #   -slope*alpha + beta = intercpt
        #   (alpha - c_x)**2 + (beta - c_y)**2 = r**2
        # Translation of coordinates:
        #   let x = alpha - c_x
        #   let y = beta - c_y
        # We have:
        #   a*x + b*y = c
        #   x**2 + y**2 = r**2
        # The solution in x-y coordinates should be transormed back:
        #   alpha = x + c_x
        #   beta = y + c_y
        a = -slope
        b = 1
        c = intercpt
        Delta = sqrt(r**2*(a**2+b**2)-c**2)
        if Delta < 0:
            raise Exception('No intersection for given line and circle')
        x_solutions = (a*c + b*Delta)/(a**2+b**2), (a*c - b*Delta)/(a**2+b**2)
        y_solutions = (b*c - a*Delta)/(a**2+b**2), (b*c + a*Delta)/(a**2+b**2)
        return (x_solutions[0]+centerx, x_solutions[1]+centerx), (y_solutions[0]+centery, y_solutions[1]+centery)

from utility import csv_row_reader

class VanGogh_pyPlotter(VanGogh):
    """VanGogh_pyPlotter is child class of VanGogh for plotting IM geometry in matplotlib."""
    def __init__(self, im, child_index=2):
        super(VanGogh_pyPlotter, self).__init__(im, child_index)

        self.add_line = self.draw_line
        self.add_arc  = self.draw_arc

        self.fig = plt.figure(figsize=(8, 8), facecolor='w', edgecolor='k')
        self.ax = self.fig.add_subplot(111, aspect='equal')
        # self.ax = self.fig.gca()
        # self.ax.set_xlim([-74,-40])
        self.ax.set_ylim([-10,25])

        # plt.gcf().gca().invert_yaxis()

        # Embeded TikZ function
        self.tikz = VanGogh_TikZPlotter()

    @staticmethod
    def mirror_and_copyrotate(Q, Radius, fraction):
        # Mirror
        # femm.mi_selectcircle(0,0,Radius+EPS,SELECT_ALL) # this EPS is sometime necessary to selece the arc at Radius.
        # femm.mi_mirror2(0,0,-Radius,0, SELECT_ALL)

        # Rotate
        # femm.mi_selectcircle(0,0,Radius+EPS,SELECT_ALL)
        # femm.mi_copyrotate2(0, 0, 360./Q, int(Q)/fraction, SELECT_ALL)

        return

    def draw_arc(self, p1, p2, angle, center=(0,0), **kwarg):
        # Embeded TikZ function
        self.tikz.draw_arc(p1, p2, center, relangle=-0.5*angle)

        # center, radius = self.find_center_and_radius_of_a_circle_using_2_points_and_arc_angle(p1, p2, angle) # ordered p1 and p2 are
        # print sqrt((p1[0]-center[0])**2 + (p1[1]-center[0])**2)
        # print sqrt((p2[0]-center[0])**2 + (p2[1]-center[0])**2)
        distance = sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)
        radius   = distance*0.5 / sin(angle*0.5)

        obj = self.pyplot_arc(radius, angle_span=angle, rotation=atan2(p1[1]-center[1],p1[0]-center[0]), center=center, **kwarg)
        self.plot_object_list.append(obj[0])
        return obj

    def pyplot_arc(self, radius, angle_span=3.14, center=(0,0), rotation=0, maxseg=0.1, ax=None, **kwarg):
        # Plot an arc starting from x-axis to axis aligned to angle_span. Use rotation to change the beginning axis.

        # make sure ax is not empty
        if ax is None:
            ax = self.ax
            # ax = plt.gcf().gca()

        # turn Polar into Cartesian
        def xy(radius, phi, center):
            return radius*np.cos(phi) + center[0], radius*np.sin(phi) + center[1]

        # get list of points
        # phis = np.arange(rotation, rotation+angle_span, 2*radius*np.pi / (360./maxseg))
        phis = np.linspace(rotation, rotation+angle_span, 360./maxseg)
        return ax.plot( *xy(radius, phis, center), c='k', **kwarg)
    
    def draw_line(self, p1, p2, ax=None, **kwarg):
        # Embeded TikZ function
        self.tikz.draw_line(p1, p2)

        # make sure ax is not empty
        if ax is None:
            ax = self.ax

        obj = ax.plot( [p1[0], p2[0]], [p1[1], p2[1]], c='k', **kwarg)
        self.plot_object_list.append(obj[0])
        return obj

import pyx
# from pyx import *
# text.set(mode="latex") 
# text.preamble(r"\usepackage{ulem}")  # underlining text...
# text.preamble(r"\usepackage{anyfontsize}")  # one way of scaling fonts (except math fonts)....

class VanGogh_TikZPlotter():
    """VanGogh_TikzPlotter is child class of VanGogh for plotting IM geometry in TikZ.
        The idea here is to incorporate this class with VanGogh_pyPlotter, 
        such that what you see in Matplotlib is what you get in TikZ."""
    def __init__(self):
        self.add_line = self.draw_line
        self.add_arc  = self.draw_arc
        self.file = open('VanGogh_TikZ.tex', 'w')
        self.c = pyx.canvas.canvas()

        self.track_path = []

    @staticmethod
    def mirror_and_copyrotate(Q, Radius, fraction):
        return

    def draw_line(self, p1, p2, ax=None, **kwarg):
        msg = '\\draw[black] (%12.8f,%12.8f) -- (%12.8f,%12.8f);\n' % (p1[0], p2[0], p1[1], p2[1])
        # print msg
        self.file.write(msg)

        # PyX
        self.c.stroke(pyx.path.line(p1[0], p1[1], p2[0], p2[1]), [pyx.color.rgb.black, pyx.style.linewidth.THICK])
        self.c.fill(pyx.path.circle(p1[0], p1[1], 0.075)) # use this if THICK is used.
        self.c.fill(pyx.path.circle(p2[0], p2[1], 0.075)) # use this if THICK is used.
            # p = pyx.path.path(pyx.path.moveto(p1[0], p1[1]), pyx.path.lineto(p2[0], p2[1]), pyx.path.closepath())
            # self.c.stroke(p, [pyx.color.rgb.black, pyx.style.linewidth.Thick])
        if 'untrack' not in kwarg:
            self.track_path.append([p1[0], p1[1], p2[0], p2[1]])

    def draw_arc(self, startxy, endxy, centerxy=(0,0), **kwarg):
        # %DRAWARC Draw an arc in the current TikZ graphic.
        # %   drawarc(mn, [center_x,_y], [start_x, _y], [end_x, _y])
        # %       draws an arc
        
        # % Extract and rename variables for processing        
        h = centerxy[0]
        k = centerxy[1]
        x1 = startxy[0]
        y1 = startxy[1]
        x2 = endxy[0]
        y2 = endxy[1]

        # % Need to calculate: x, y, start, stop, radius...
        x = x1
        y = y1
        radius = sqrt((x1 - h)**2 + (y1 - k)**2)
        
        # % These are both "correct" ways of determining start variable,
        # % but only one works depending on orientation of arc. Hacks
        # % below help solve this issue.
        start = acos((x1 - h) / radius)
        start2 = asin((y1 - k) / radius)
        temp = (x2 - h) / radius
        if temp>=1:
            # print 'temp=', temp
            temp = 1
        elif temp<=-1:
            # print 'temp=', temp
            temp = -1
        stop = acos(temp)
        
        # % ----------
        # % HACKS START
        # % ----------
        
        # % The code in this "hack" section is special cases for drawing
        # % an arc... This involves mostly cases where x1==x2 or
        # % y1==y2. I have not proved these are the correct solutions,
        # % but I have tried many arcs using various stator geometries
        # % and it seems to work for all cases...  �_(?)_/�            
        if (start > stop):
           start = -start
           stop = -stop
        
        # % Trying to check for start == stop.
        # % BUT, issues with rounding and floating numbers....
        # % Using this delta comparison fixes one silly case where
        # % rounding causes issues for direct comparison for equality...
        # % 
        # % Silly!
        # %
        if (abs(start - stop) < 0.0000001):
            if (y1-k < 0):
                start = -start
            elif (x1-h < 0):
                stop = start + (2 * start2)

        # % --------
        # % HACKS END
        # % --------
        def rad2deg(rad):
            return rad/pi*180.
        msg = '\\draw[black] (%12.8f,%12.8f) arc (%12.8f:%12.8f:%12.8f);\n' % (x, y, rad2deg(start), rad2deg(stop), radius)
        # print msg
        self.file.write(msg)

        # Below is PyX | Above is TikZ
        # PyX
        textattrs = [pyx.text.halign.center, pyx.text.vshift.middlezero]
        X = pyx.text.text(startxy[0], startxy[1], r"", textattrs) # must have textattrs or normsubpath cannot close erro
        Y = pyx.text.text(endxy[0], endxy[1], r"", textattrs)
        self.c.stroke(pyx.connector.arc(X, Y, boxdists=[0.0, 0.0], relangle=kwarg['relangle']/pi*180), [pyx.color.rgb.black, pyx.style.linewidth.THICK])
            # print '[%g,%g,%g,%g,%g],'%(startxy[0], startxy[1], endxy[0], endxy[1], 0.5*kwarg['relangle'])
        if 'untrack' not in kwarg:
            self.track_path.append([startxy[0], startxy[1], endxy[0], endxy[1], centerxy[0], centerxy[1], kwarg['relangle']])

# wiki: https://en.wikipedia.org/wiki/Tangent_lines_to_circles
def get_tangent_points_of_two_circles(center1, radius1, center2, radius2):
    x1, y1 = center1[0], center1[1]
    r      = radius1
    x2, y2 = center2[0], center2[1]
    R      = radius2

    gamma  = - atan( (y2-y1) / (x2-x1) )
    beta   = asin( (R-r) / sqrt((x2-x1)**2+(y2-y1)**2) )
    alpha  = gamma - beta

    x3, y3 = x1 + r*cos(0.5*pi - alpha), y1 + r*sin(0.5*pi - alpha), 
    x4, y4 = x2 + R*cos(0.5*pi - alpha), y2 + R*sin(0.5*pi - alpha), 

    return (x3, y3), (x4, y4)

# Eric required plot geometry with LaTeX labels
if __name__ == '!__main__':
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "Times New Roman"

    myfontsize = 13.5
    plt.rcParams.update({'font.size': myfontsize})

    from population import *
    im_list = []
    with open(r'D:\OneDrive - UW-Madison\c\pop\initial_design.txt', 'r') as f: 
        for row in csv_row_reader(f):
            # fea_config_dict = {}
            # fea_config_dict['DPNV'] = True
            # fea_config_dict['flag_optimization'] = False
            # fea_config_dict['End_Ring_Resistance'] = 0.0
            im = bearingless_induction_motor_design([row[0]]+[float(el) for el in row[1:]], None)
            im_list.append(im)
    # print im.show(toString=True)

    # 示意图而已，改改尺寸吧
    im.Radius_OuterStatorYoke -= 37
    im.Radius_InnerStatorYoke -= 20
    im.Radius_Shaft += 20
    # im.Location_RotorBarCenter2 += 5 # this will change the shape of rotor slot

    vg = VanGogh_pyPlotter(im, CUSTOM)
    vg.draw_model()

    # PyX
    vg.tikz.c = pyx.canvas.canvas() # clear the canvas because we want to redraw 90 deg with the data vg.tikz.track_path
    from copy import deepcopy
    def pyx_draw_path(vg, path, sign=1):
        if len(path) == 4:
            vg.tikz.draw_line(path[:2], path[2:4], untrack=True)
        else:
            vg.tikz.draw_arc(path[:2], path[2:4], path[4:6], relangle=sign*path[6], untrack=True)
    def rotate(_, x, y):
        return cos(_)*x + sin(_)*y, -sin(_)*x + cos(_)*y
    def is_at_stator(im, path):
        return sqrt(path[0]**2 + path[1]**2) > im.Radius_OuterRotor + 0.5*im.Length_AirGap
    # 整体转动90度。
    for path in vg.tikz.track_path:
        path[0], path[1] = rotate(0.5*pi, path[0], path[1])
        path[2], path[3] = rotate(0.5*pi, path[2], path[3])
        pyx_draw_path(vg, path)
    track_path_backup = deepcopy(vg.tikz.track_path)

    # Rotate
    for path in deepcopy(vg.tikz.track_path):
        if is_at_stator(im, path):
            Q = im.Qs
        else:
            Q = im.Qr
        _ = 2*pi/Q
        path[0], path[1] = rotate(_, path[0], path[1])
        path[2], path[3] = rotate(_, path[2], path[3])
        pyx_draw_path(vg, path)

    # Mirror
    for path in (vg.tikz.track_path): # track_path is passed by reference and is changed by mirror
        path[0] *= -1
        path[2] *= -1
        pyx_draw_path(vg, path, sign=-1)
    for path in (vg.tikz.track_path):
        if sqrt(path[0]**2 + path[1]**2) > im.Radius_OuterRotor + 0.5*im.Length_AirGap:
            Q = im.Qs
        else:
            Q = im.Qr
        _ = 2*pi/Q
        path[0], path[1] = rotate(_, path[0], path[1])
        path[2], path[3] = rotate(_, path[2], path[3])
        pyx_draw_path(vg, path, sign=-1)

    def pyx_label_dimention(track_path):
        def myarrow(c, PA, PB):
            c.stroke(path.line(PA[0], PA[1], PB[0], PB[1]),
                     [style.linewidth.Thick, style.linestyle.solid, color.rgb.black,
                      deco.earrow([deco.stroked([color.rgb.black, style.linejoin.round]),
                                   deco.filled([color.rgb.black])], size=0.5)])
        def mytext(c, PText, label):
            print(PText)
            textattrs = [text.halign.center, text.vshift.middlezero, text.size(5), trafo.scale(2)]
            handle = text.text(PText[0], PText[1], label, textattrs)
            c.insert(handle)

        def add_label_inside(c, label, PA, PB, PText=None):
            if PText is None:
                PText = 0.5*(PA + PB)
            mytext(c, PText, label)
            myarrow(c, PA, PB)
            myarrow(c, PB, PA)

        def add_label_outside(c, label, PA, PB, PText=None, ext_factor=1, left_or_right=True):
            vector_AB = np.array((PB[0]-PA[0], PB[1]-PA[1]))
            PB, PBExt = PB, PB + vector_AB * ext_factor
            myarrow(c, PBExt, PB)
            PA, PAExt = PA, PA - vector_AB * ext_factor
            myarrow(c, PAExt, PA)
            if PText is None:
                PText = 0.5*(PA + PB)
            mytext(c, PText, label)

        def mymarker(l_PA):
            for P in l_PA:
                vg.tikz.c.fill(pyx.path.circle(P[0], P[1], 0.25))

        count_stator = 0
        count_rotor = 0
        dictKeyPoints = {}
        for ind, el in enumerate(track_path):
            if is_at_stator(im, el):
                count_stator += 1
                print('S-Path'+str(count_stator), el)
                dictKeyPoints['S-Path'+str(count_stator)] = el
            else:
                count_rotor += 1
                print('R-Path'+str(count_rotor), el)
                dictKeyPoints['R-Path'+str(count_rotor)] = el
        # above SP means stator path, below SP means stator point

        d = dictKeyPoints # too long to type
        RP2 = array(d['R-Path6'][0:2])
        RP3 = array(d['R-Path1'][0:2])
        RP4 = array(d['R-Path2'][0:2])
        RP5 = array(d['R-Path7'][0:2])
        RP6 = array(d['R-Path3'][0:2])
        RP7 = array(d['R-Path4'][2:4])

        SP2 = array(d['S-Path1'][0:2])
        SP3 = array(d['S-Path3'][0:2])
        SP4 = array(d['S-Path3'][2:4])
        SP5 = array(d['S-Path4'][2:4])
        _s = 2*pi/im.Qs
        _r = 2*pi/im.Qr


        PA = 0.5 * (SP4 + SP5) 
        PB = -PA[0], PA[1]
        PText = 0.5*(PA[0]+PB[0]), PA[1] + 1.5
        add_label_inside(vg.tikz.c, r'$w_{st}$', PA, PB, PText=PText)

        PA = 0.5 * (SP2 + SP3) 
        PB = rotate(_s, -PA[0], PA[1])
        PText = 0.5*(PA[0]+PB[0]), PA[1] + 1.5
        add_label_inside(vg.tikz.c, r'$\theta_{so}$', PA, PB, PText=PText)

        PA = SP2
        PB = SP3
        PA = rotate(-0.9*_s, PA[0], PA[1])
        PB = rotate(-0.9*_s, PB[0], PB[1])
            # PText = 0.5*(PA[0]+PB[0]), PA[1] + 2.25
            # add_label_inside(vg.tikz.c, r'$d_{so}$', PA, PB, PText=PText)
        PText = 0.5*(PA[0]+PB[0]), PA[1] + 3
        add_label_outside(vg.tikz.c, r'$d_{so}$', PA, PB, PText=PText)
        # 辅助线
        vg.tikz.c.stroke(path.line(-SP2[0], SP2[1], PA[0], PA[1]),
                 [style.linewidth(0.08), style.linestyle.dashed, color.gray(0.5)])
        vg.tikz.c.stroke(path.line(-SP3[0], SP3[1], PB[0], PB[1]),
                 [style.linewidth(0.08), style.linestyle.dashed, color.gray(0.5)])


        PA = 0.5 * (RP6 + RP7) 
        PB = rotate(_r, -PA[0], PA[1])
        PText = 0.5*(PA[0]+PB[0]), PA[1] + 1.5
        add_label_inside(vg.tikz.c, r'$w_{rt}$', PA, PB, PText=PText)

        PA = 0.5 * (RP4 + RP5)
        PB = -PA[0], PA[1]
            # PText = 0.5*(PA[0]+PB[0]), PA[1] + 1.45
            # add_label_inside(vg.tikz.c, r'$w_{ro}$', PA, PB, PText=PText)
        PText = 0.5*(PA[0]+PB[0])+3, PA[1] - 0.1
        add_label_outside(vg.tikz.c, r'$w_{ro}$', PA, PB, PText=PText)

        PA = RP4
        PB = RP5
        PA = rotate(1.08*_s, PA[0], PA[1])
        PB = rotate(1.08*_s, PB[0], PB[1])
            # PText = 0.5*(PA[0]+PB[0]) - 1.5, PA[1] - 1.
            # add_label_inside(vg.tikz.c, r'$d_{ro}$', PA, PB, PText=PText)
        PText = 0.5*(PA[0]+PB[0]) - 1.5, PA[1] - 2
        add_label_outside(vg.tikz.c, r'$d_{ro}$', PA, PB, PText=PText)
        # 辅助线
        RP4 = rotate(_r, RP4[0], RP4[1])
        RP5 = rotate(_r, RP5[0], RP5[1])
        vg.tikz.c.stroke(path.line(RP4[0], RP4[1], PA[0], PA[1]),
                 [style.linewidth(0.08), style.linestyle.dashed, color.gray(0.5)])
        vg.tikz.c.stroke(path.line(RP5[0], RP5[1], PB[0], PB[1]),
                 [style.linewidth(0.08), style.linestyle.dashed, color.gray(0.5)])

        # Air gap
        PA = array([-cos(0.5*pi-0.5*_r), sin(0.5*pi-0.5*_r)])*im.Radius_OuterRotor
        PB = array([-cos(0.5*pi-0.5*_r), sin(0.5*pi-0.5*_r)])*(im.Radius_OuterRotor+im.Length_AirGap)
        PText = 0.5*(PA[0]+PB[0]), PA[1] - 1
        add_label_inside(vg.tikz.c, r'$L_{g}$', PA, PB, PText=PText)

        # Outer Rotor Diameter
        PA = array([cos(0.5*pi-0.5*_r), sin(0.5*pi-0.5*_r)])*im.Radius_OuterRotor
        PB = array([cos(0.5*pi-0.5*_r), sin(0.5*pi-0.5*_r)])*(im.Radius_OuterRotor*0.975)
        PText = PA[0] + 1.2, PA[1] - 1.5
        myarrow(vg.tikz.c, PB, PA)
        mytext(vg.tikz.c, PText, r'$r_{or}$')

        # Outer Sator Diameter
        PA = array([cos(0.5*pi-0.5*_r), sin(0.5*pi-0.5*_r)])*im.Radius_OuterStatorYoke
        PB = array([cos(0.5*pi-0.5*_r), sin(0.5*pi-0.5*_r)])*(im.Radius_OuterStatorYoke*0.975)
        PText = PA[0] - 1.5, PA[1] - 1.5
        myarrow(vg.tikz.c, PB, PA)
        mytext(vg.tikz.c, PText, r'$r_{os}$')

        # mymarker([PA, PB, PText])

        # print dir(color.rgb)
        # print dir(style.linestyle)
        # vg.tikz.c.stroke(path.line(0, 0, 0, 40),
        #          [style.linewidth.THICK, style.linestyle.solid, color.gray(0.5),
        #           deco.earrow([deco.stroked([color.rgb.black, style.linejoin.round]),
        #                        deco.filled([color.rgb.black])], size=0.5)])

    from pyx import *
    pyx_label_dimention(track_path_backup)
    vg.tikz.c.writePDFfile("pyx_output")
    # vg.tikz.c.writeEPSfile("pyx_output")
    print('Write to pdf file: pyx_output.pdf.')
    quit()

    # 以下应该是旧代码，用matplotlib画电机二维截面
    BIAS = 1
    Stator_Sector_Angle = 2*pi/im.Qs*0.5
    Rotor_Sector_Angle = 2*pi/im.Qr*0.5

    # draw the mirrored slot
    new_objs = []
    for obj in vg.rotor_object_list:
        xy_data = obj.get_xydata().T
        new_obj = vg.ax.plot(xy_data[0], -array(xy_data[1]), 'k-', zorder=0) # 'k--'
        new_objs.append(new_obj[0])
    vg.rotor_object_list += new_objs

    new_objs = []
    for obj in vg.stator_object_list:
        xy_data = obj.get_xydata().T
        new_obj = vg.ax.plot(xy_data[0], -array(xy_data[1]), 'k-', zorder=0) # 'k--'
        new_objs.append(new_obj[0])
    vg.stator_object_list += new_objs

    # draw the rotated slot
    # rotor 
    angle = Rotor_Sector_Angle*2
    TR = array( [ [cos(angle), sin(angle)],
                 [-sin(angle), cos(angle)] ] )
    for obj in vg.rotor_object_list:
        xy_data = obj.get_xydata().T
        xy_data_rot = np.dot( TR, xy_data )
        # xy_data_rot = xy_data_rot.T
        vg.ax.plot(xy_data_rot[0], xy_data_rot[1], 'k-', zorder=0) # 'k--'
    # stator
    angle = Stator_Sector_Angle*2
    TS = array( [ [cos(angle), sin(angle)],
                 [-sin(angle), cos(angle)] ] )
    for obj in vg.stator_object_list:
        xy_data = obj.get_xydata().T
        xy_data_rot = np.dot( TS, xy_data )
        # xy_data_rot = xy_data_rot.T
        vg.ax.plot(xy_data_rot[0], xy_data_rot[1], 'k-', zorder=0) # 'k--'

    ################################################################
    # add labels to geometry
    ################################################################
    # xy = (-0.5*(im.Location_RotorBarCenter+im.Location_RotorBarCenter2), 2.5)
    # vg.ax.annotate('Parallel tooth', xytext=(xy[0],-(xy[1]+5)), xy=(xy[0],-xy[1]), xycoords='data', arrowprops=dict(arrowstyle="->"),fontsize=24,rotation=90)

    xy = (-im.Radius_OuterStatorYoke, 0)
    vg.ax.annotate(r'$D_{os}$', xytext=(xy[0]+2.5,xy[1]-0.25), xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->"),fontsize=24,rotation=90)
    xy = (-im.Radius_OuterRotor, -4)
    vg.ax.annotate(r'$D_{or}$', xytext=(xy[0]+2.5,xy[1]-0.25), xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->"),fontsize=24,rotation=90)
    xy = (-im.Radius_Shaft, 0)
    vg.ax.annotate(r'$D_{ir}$', xytext=(xy[0]+2.5,xy[1]-0.25), xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->"),fontsize=24,rotation=90)
    vg.ax.set_xlim([-74,-45])

    def add_label_inside(vg, label, PA, PB):
        # vector_AB = (PB[0]-PA[0], PB[1]-PA[1])
        xy, xytext = PA, PB 
        vg.ax.annotate('', xytext=xytext, xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->"),fontsize=24,rotation=90)
        xy, xytext = xytext, xy
        vg.ax.annotate('', xytext=xytext, xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->"),fontsize=24,rotation=90)
        vg.ax.text(0.5+0.5*(xy[0]+xytext[0]), 0.+0.5*(xy[1]+xytext[1]), label,fontsize=24,rotation=90,
                    bbox=dict(facecolor='w', edgecolor='w',boxstyle='round,pad=0',alpha=0.9))


    def add_label_outside(vg, label, PA, PB, ext_factor=1, left_or_right=True):
        vector_AB = np.array((PB[0]-PA[0], PB[1]-PA[1]))
        xy, xytext = PB, PB + vector_AB * ext_factor
        vg.ax.annotate('', xytext=xytext, xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->"),fontsize=24,rotation=90)
        xy, xytext = PA, PA - vector_AB * ext_factor
        vg.ax.annotate('', xytext=xytext, xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->"),fontsize=24,rotation=90)
        if left_or_right:
            vg.ax.text(0.5+PB[0], -1.75+PB[1], label,fontsize=24,rotation=90,
                        bbox=dict(facecolor='w', edgecolor='w',boxstyle='round,pad=0',alpha=0.9))
        else:
            vg.ax.text(0.5+PB[0], +1.75+PB[1], label,fontsize=24,rotation=90,
                        bbox=dict(facecolor='w', edgecolor='w',boxstyle='round,pad=0',alpha=0.9))


    def add_extension_line(vg, P, vector_ext, distance_factor=1):

        distance = sqrt(vector_ext[0]**2 + vector_ext[1]**2)
        vector_ext /= distance/distance_factor # a unit vector * distance_factor

        x = [ P[0], P[0]+vector_ext[0] ]
        y = [ P[1], P[1]+vector_ext[1] ]
        vg.ax.plot(x, y, 'k-', lw=0.6)

        return vector_ext

    #1
    Ps2 = vg.Ps_list[2-BIAS] * array((1,-1)) # mirror
    Ptemp1 = Ps2[0], Ps2[1]-5
    Ptemp1 = vg.park_transform(Ptemp1, -Stator_Sector_Angle) # this is not exact but close enough, you can find the exact point by solve a triangle with Radius_OuterRotor + Length_AirGap from Pr3
    vector_ext = np.array((Ptemp1[0] - Ps2[0], Ptemp1[1] - Ps2[1]))
    vector_ext = add_extension_line(vg, Ps2, vector_ext, distance_factor=2.5)
    Ptemp1 = (Ps2 + 0.5*vector_ext)

    Pr3 = vg.Pr_list[3-BIAS] * array((1,-1)) # mirror
    Ptemp2 = array((Pr3[0], Pr3[1]-5))
    Ptemp2 = vg.park_transform(Ptemp2, -Rotor_Sector_Angle)
    vector_ext = np.array((Ptemp2[0] - Pr3[0], Ptemp2[1] - Pr3[1]))
    vector_ext = add_extension_line(vg, Pr3, vector_ext, distance_factor=2.5)
    Ptemp2 = (Pr3 + 0.8*vector_ext) # change this factor before vector_ext to rotate the label

    add_label_outside(vg, r'$L_g$', Ptemp1, Ptemp2, ext_factor=2)

    #2
    Ps4 = vg.Ps_list[4-BIAS]
    Ps5 = vg.Ps_list[5-BIAS]
    Ptemp = 0.5*(Ps4 + Ps5)
    add_label_inside(vg, r'$W_{\rm st}$', (Ptemp[0],-Ptemp[1]), Ptemp)

    #3
    Pr6 = array(vg.Pr_list[6-BIAS])
    Pr7 = array(vg.Pr_list[7-BIAS])
    Ptemp2 = 0.5*(Pr6 + Pr7) # middle point
    Ptemp1 = Ptemp2[0], -Ptemp2[1] # mirror
    Ptemp1 = np.dot( TR, Ptemp1 ) # rotate
    # vector_r67 = Pr7 - Pr6
    # vector_r67_perpendi = (-vector_r67[1], vector_r67[0])
    # some_angle = atan2(vector_r67[1], vector_r67[0])
    add_label_inside(vg, r'$W_{\rm rt}$', Ptemp1, Ptemp2)

    #4 (This one is angle, so we add text independently)
    Ps3 = vg.Ps_list[3-BIAS]
    Ptemp1 = Ps3 * np.array((1,-1)) # mirror
    Ptemp2 = np.dot( TS, Ptemp1 ) # rotate
    if False: # treat it as width
        add_label_inside(vg, r'$w_{\rm open,s}$', Ps3, Ptemp2)
    else: # it is an angle, ok.
        Ps2 = vg.Ps_list[2-BIAS]
        vector_ext = Ps3 - Ps2
        vector_ext = add_extension_line(vg, Ps3, vector_ext, distance_factor=3)
        Ptemp1 = (Ps3 + 0.5*vector_ext)

        vector_ext *= array([1,-1]) # mirror
        vector_ext = np.dot( TS, vector_ext ) # rotate
        vector_ext = add_extension_line(vg, Ptemp2, vector_ext, distance_factor=3)
        Ptemp2 = (Ptemp2 + 0.5*vector_ext)

        some_angle = atan2(vector_ext[1], vector_ext[0])
        print(some_angle/pi*180, 180 - Stator_Sector_Angle/pi*180)
        vg.pyplot_arc(im.Radius_OuterRotor+im.Length_AirGap+im.Width_StatorTeethHeadThickness+2, 
            angle_span=0.95*im.Angle_StatorSlotOpen/180*pi, rotation=some_angle-0.475*im.Angle_StatorSlotOpen/180*pi, center=(0,0), lw=0.6)
        # vg.ax.text(Ptemp1[0], Ptemp1[1], r'$w_{\rm open,s}$',fontsize=24,rotation=90,
        #             bbox=dict(facecolor='w', edgecolor='r',boxstyle='round,pad=0',alpha=0.9))
        vg.ax.text(0.5*(Ptemp1[0]+Ptemp2[0])-5, 0.5*(Ptemp1[1]+Ptemp2[1]), r'$\theta_{\rm so}$',fontsize=24,rotation=90,
                    bbox=dict(facecolor='w', edgecolor='w',boxstyle='round,pad=0',alpha=0.9))



    #5
    Pr4 = vg.Pr_list[4-BIAS]
    Pr5 = vg.Pr_list[5-BIAS]
    Ptemp1 = 0.5*(Pr4+Pr5)
    Ptemp2 = Ptemp1 * np.array((1,-1)) # mirror
    add_label_outside(vg, r'$W_{\rm ro}$', Ptemp1, Ptemp2, ext_factor=2, left_or_right=False)

    #6
    Ps2 = vg.Ps_list[2-BIAS]
    Ps3 = vg.Ps_list[3-BIAS] 
    Ptemp1 = np.dot( TS, Ps2 )
    Ptemp11 = np.dot( TS, Ps3 )
    vector_ext = -np.array((Ptemp11[0] - Ptemp1[0], Ptemp11[1] - Ptemp1[1]))
    vector_ext = np.array((-vector_ext[1], vector_ext[0]))
    vector_ext = add_extension_line(vg, Ptemp1, vector_ext, distance_factor=1)
    Ptemp1 = (Ptemp1 + 0.5*vector_ext)

    Ptemp2 = np.dot( TS, Ps3 )
    Ptemp22 = np.dot( TS, Ps2 )
    vector_ext = np.array((Ptemp22[0] - Ptemp2[0], Ptemp22[1] - Ptemp2[1]))
    vector_ext = np.array((-vector_ext[1], vector_ext[0]))
    vector_ext = add_extension_line(vg, Ptemp2, vector_ext, distance_factor=1)
    Ptemp2 = (Ptemp2 + 0.5*vector_ext)

    add_label_outside(vg, r'$d_{\rm so}$', Ptemp1, Ptemp2, ext_factor=3)


    #7
    Pr4 = vg.Pr_list[4-BIAS]
    Ptemp1 = np.dot( TR, Pr4 )
    Ptemp11 = np.array((Pr4[0], Pr4[1]+5))
    Ptemp11 = np.dot( TR, Ptemp11 )
    vector_ext = np.array((Ptemp11[0] - Ptemp1[0], Ptemp11[1] - Ptemp1[1]))
    vector_ext = add_extension_line(vg, Ptemp1, vector_ext, distance_factor=6)
    Ptemp1 = (Ptemp1 + 0.5*vector_ext)

    Pr5 = vg.Pr_list[5-BIAS] 
    Ptemp2 = np.dot( TR, Pr5 )
    Ptemp22 = np.array((Pr5[0], Pr5[1]+5))
    Ptemp22 = np.dot( TR, Ptemp22 )
    vector_ext = np.array((Ptemp22[0] - Ptemp2[0], Ptemp22[1] - Ptemp2[1]))
    vector_ext = add_extension_line(vg, Ptemp2, vector_ext, distance_factor=6)
    Ptemp2 = (Ptemp2 + 0.5*vector_ext)

    add_label_outside(vg, r'$d_{\rm ro}$', Ptemp1, Ptemp2, ext_factor=3)



    vg.ax.get_xaxis().set_visible(False)
    vg.ax.get_yaxis().set_visible(False)
    vg.fig.tight_layout()
    # vg.fig.savefig(r'D:\OneDrive\[00]GetWorking\32 blimopti\p2019_ecce_bearingless_induction_full_paper\images\CAD_CrossSection.png')
    vg.fig.savefig(r'D:\OneDrive\[00]GetWorking\31 Bearingless_Induction_FEA_Model\p2019_iemdc_bearingless_induction full paper\images\CAD_CrossSection_Rot.png')
    show()
    quit()


# My plot geometry with LaTeX labels
if __name__ == '!__main__':
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams["font.family"] = "Times New Roman"

    myfontsize = 13.5
    plt.rcParams.update({'font.size': myfontsize})

    from population import *
    im_list = []
    with open(r'D:\OneDrive - UW-Madison\c\pop\initial_design.txt', 'r') as f: 
        for row in csv_row_reader(f):
            im = bearingless_induction_motor_design([row[0]]+[float(el) for el in row[1:]], None)
            im_list.append(im)
    # print im.show(toString=True)

    # 示意图而已，改改尺寸吧
    im.Radius_OuterStatorYoke -= 25
    im.Radius_InnerStatorYoke -= 20
    im.Radius_Shaft += 10

    vg = VanGogh_pyPlotter(im, CUSTOM)
    vg.draw_model()

    # 
    BIAS = 1
    Stator_Sector_Angle = 2*pi/im.Qs*0.5
    Rotor_Sector_Angle = 2*pi/im.Qr*0.5

    # draw the mirrored slot
    new_objs = []
    for obj in vg.rotor_object_list:
        xy_data = obj.get_xydata().T
        new_obj = vg.ax.plot(xy_data[0], -array(xy_data[1]), 'k-', zorder=0) # 'k--'
        new_objs.append(new_obj[0])
    vg.rotor_object_list += new_objs

    new_objs = []
    for obj in vg.stator_object_list:
        xy_data = obj.get_xydata().T
        new_obj = vg.ax.plot(xy_data[0], -array(xy_data[1]), 'k-', zorder=0) # 'k--'
        new_objs.append(new_obj[0])
    vg.stator_object_list += new_objs

    # draw the rotated slot
    # rotor 
    angle = Rotor_Sector_Angle*2
    TR = array( [ [cos(angle), sin(angle)],
                 [-sin(angle), cos(angle)] ] )
    for obj in vg.rotor_object_list:
        xy_data = obj.get_xydata().T
        xy_data_rot = np.dot( TR, xy_data )
        # xy_data_rot = xy_data_rot.T
        vg.ax.plot(xy_data_rot[0], xy_data_rot[1], 'k-', zorder=0) # 'k--'
    # stator
    angle = Stator_Sector_Angle*2
    TS = array( [ [cos(angle), sin(angle)],
                 [-sin(angle), cos(angle)] ] )
    for obj in vg.stator_object_list:
        xy_data = obj.get_xydata().T
        xy_data_rot = np.dot( TS, xy_data )
        # xy_data_rot = xy_data_rot.T
        vg.ax.plot(xy_data_rot[0], xy_data_rot[1], 'k-', zorder=0) # 'k--'

    ################################################################
    # add labels to geometry
    ################################################################
    # xy = (-0.5*(im.Location_RotorBarCenter+im.Location_RotorBarCenter2), 2.5)
    # vg.ax.annotate('Parallel tooth', xytext=(xy[0],-(xy[1]+5)), xy=(xy[0],-xy[1]), xycoords='data', arrowprops=dict(arrowstyle="->"))

    xy = (-im.Radius_OuterStatorYoke, 0)
    vg.ax.annotate(r'$r_{so}=86.9$ mm', xytext=(xy[0]+2.5,xy[1]-0.25), xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->"))
    xy = (-im.Radius_OuterRotor, -4)
    vg.ax.annotate(r'$r_{ro}=47.1$ mm', xytext=(xy[0]+2.5,xy[1]-0.25), xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->"))

    def add_label_inside(vg, label, PA, PB):
        # vector_AB = (PB[0]-PA[0], PB[1]-PA[1])
        xy, xytext = PA, PB 
        vg.ax.annotate('', xytext=xytext, xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->"))
        xy, xytext = xytext, xy
        vg.ax.annotate('', xytext=xytext, xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->"))
        vg.ax.text(0.5+0.5*(xy[0]+xytext[0]), 0.+0.5*(xy[1]+xytext[1]), label,
                    bbox=dict(facecolor='w', edgecolor='r',boxstyle='round,pad=0',alpha=0.9))


    def add_label_outside(vg, label, PA, PB, ext_factor=1):
        vector_AB = np.array((PB[0]-PA[0], PB[1]-PA[1]))
        xy, xytext = PB, PB + vector_AB * ext_factor
        vg.ax.annotate('', xytext=xytext, xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->"))
        xy, xytext = PA, PA - vector_AB * ext_factor
        vg.ax.annotate('', xytext=xytext, xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->"))
        vg.ax.text(0.5+PB[0], -1.75+PB[1], label,
                    bbox=dict(facecolor='w', edgecolor='r',boxstyle='round,pad=0',alpha=0.9))


    def add_extension_line(vg, P, vector_ext, distance_factor=1):

        distance = sqrt(vector_ext[0]**2 + vector_ext[1]**2)
        vector_ext /= distance/distance_factor # a unit vector * distance_factor

        x = [ P[0], P[0]+vector_ext[0] ]
        y = [ P[1], P[1]+vector_ext[1] ]
        vg.ax.plot(x, y, 'k-', lw=0.6)

        return vector_ext

    #1
    Ps2 = vg.Ps_list[2-BIAS] * array((1,-1)) # mirror
    Ptemp1 = Ps2[0], Ps2[1]-5
    Ptemp1 = vg.park_transform(Ptemp1, -Stator_Sector_Angle) # this is not exact but close enough, you can find the exact point by solve a triangle with Radius_OuterRotor + Length_AirGap from Pr3
    vector_ext = np.array((Ptemp1[0] - Ps2[0], Ptemp1[1] - Ps2[1]))
    vector_ext = add_extension_line(vg, Ps2, vector_ext, distance_factor=2.5)
    Ptemp1 = (Ps2 + 0.5*vector_ext)

    Pr3 = vg.Pr_list[3-BIAS] * array((1,-1)) # mirror
    Ptemp2 = array((Pr3[0], Pr3[1]-5))
    Ptemp2 = vg.park_transform(Ptemp2, -Rotor_Sector_Angle)
    vector_ext = np.array((Ptemp2[0] - Pr3[0], Ptemp2[1] - Pr3[1]))
    vector_ext = add_extension_line(vg, Pr3, vector_ext, distance_factor=2.5)
    Ptemp2 = (Pr3 + 0.5*vector_ext)

    add_label_outside(vg, r'$\delta$', Ptemp1, Ptemp2, ext_factor=2)

    #2
    Ps4 = vg.Ps_list[4-BIAS]
    Ps5 = vg.Ps_list[5-BIAS]
    Ptemp = 0.5*(Ps4 + Ps5)
    add_label_inside(vg, r'$b_{\rm tooth,s}$', (Ptemp[0],-Ptemp[1]), Ptemp)

    #3
    Pr6 = array(vg.Pr_list[6-BIAS])
    Pr7 = array(vg.Pr_list[7-BIAS])
    Ptemp2 = 0.5*(Pr6 + Pr7) # middle point
    Ptemp1 = Ptemp2[0], -Ptemp2[1] # mirror
    Ptemp1 = np.dot( TR, Ptemp1 ) # rotate
    # vector_r67 = Pr7 - Pr6
    # vector_r67_perpendi = (-vector_r67[1], vector_r67[0])
    # some_angle = atan2(vector_r67[1], vector_r67[0])
    add_label_inside(vg, r'$b_{\rm tooth,r}$', Ptemp1, Ptemp2)

    #4 (This one is angle)
    Ps3 = vg.Ps_list[3-BIAS]
    Ptemp1 = Ps3 * np.array((1,-1)) # mirror
    Ptemp2 = np.dot( TS, Ptemp1 ) # rotate
    if False: # treat it as width
        add_label_inside(vg, r'$w_{\rm open,s}$', Ps3, Ptemp2)
    else: # it is an angle, ok.
        Ps2 = vg.Ps_list[2-BIAS]
        vector_ext = Ps3 - Ps2
        vector_ext = add_extension_line(vg, Ps3, vector_ext, distance_factor=3)
        Ptemp1 = (Ps3 + 0.5*vector_ext)

        vector_ext *= array([1,-1]) # mirror
        vector_ext = np.dot( TS, vector_ext ) # rotate
        vector_ext = add_extension_line(vg, Ptemp2, vector_ext, distance_factor=3)
        Ptemp2 = (Ptemp2 + 0.5*vector_ext)

        some_angle = atan2(vector_ext[1], vector_ext[0])
        print(some_angle/pi*180, 180 - Stator_Sector_Angle/pi*180)
        vg.pyplot_arc(im.Radius_OuterRotor+im.Length_AirGap+im.Width_StatorTeethHeadThickness+2, 
            angle_span=0.95*im.Angle_StatorSlotOpen/180*pi, rotation=some_angle-0.475*im.Angle_StatorSlotOpen/180*pi, center=(0,0), lw=0.6)
        # vg.ax.text(Ptemp1[0], Ptemp1[1], r'$w_{\rm open,s}$',
        #             bbox=dict(facecolor='w', edgecolor='r',boxstyle='round,pad=0',alpha=0.9))
        vg.ax.text(0.5*(Ptemp1[0]+Ptemp2[0])-5, 0.5*(Ptemp1[1]+Ptemp2[1]), r'$w_{\rm open,s}$',
                    bbox=dict(facecolor='w', edgecolor='r',boxstyle='round,pad=0',alpha=0.9))



    #5
    Pr4 = vg.Pr_list[4-BIAS]
    Pr5 = vg.Pr_list[5-BIAS]
    Ptemp1 = 0.5*(Pr4+Pr5)
    Ptemp2 = Ptemp1 * np.array((1,-1)) # mirror
    add_label_outside(vg, r'$w_{\rm open,r}$', Ptemp1, Ptemp2, ext_factor=2)

    #6
    Ps2 = vg.Ps_list[2-BIAS]
    Ps3 = vg.Ps_list[3-BIAS] 
    Ptemp1 = np.dot( TS, Ps2 )
    Ptemp11 = np.dot( TS, Ps3 )
    vector_ext = -np.array((Ptemp11[0] - Ptemp1[0], Ptemp11[1] - Ptemp1[1]))
    vector_ext = np.array((-vector_ext[1], vector_ext[0]))
    vector_ext = add_extension_line(vg, Ptemp1, vector_ext, distance_factor=1)
    Ptemp1 = (Ptemp1 + 0.5*vector_ext)

    Ptemp2 = np.dot( TS, Ps3 )
    Ptemp22 = np.dot( TS, Ps2 )
    vector_ext = np.array((Ptemp22[0] - Ptemp2[0], Ptemp22[1] - Ptemp2[1]))
    vector_ext = np.array((-vector_ext[1], vector_ext[0]))
    vector_ext = add_extension_line(vg, Ptemp2, vector_ext, distance_factor=1)
    Ptemp2 = (Ptemp2 + 0.5*vector_ext)

    add_label_outside(vg, r'$h_{\rm head,s}$', Ptemp1, Ptemp2, ext_factor=3)


    #7
    Pr4 = vg.Pr_list[4-BIAS]
    Ptemp1 = np.dot( TR, Pr4 )
    Ptemp11 = np.array((Pr4[0], Pr4[1]+5))
    Ptemp11 = np.dot( TR, Ptemp11 )
    vector_ext = np.array((Ptemp11[0] - Ptemp1[0], Ptemp11[1] - Ptemp1[1]))
    vector_ext = add_extension_line(vg, Ptemp1, vector_ext, distance_factor=6)
    Ptemp1 = (Ptemp1 + 0.5*vector_ext)

    Pr5 = vg.Pr_list[5-BIAS] 
    Ptemp2 = np.dot( TR, Pr5 )
    Ptemp22 = np.array((Pr5[0], Pr5[1]+5))
    Ptemp22 = np.dot( TR, Ptemp22 )
    vector_ext = np.array((Ptemp22[0] - Ptemp2[0], Ptemp22[1] - Ptemp2[1]))
    vector_ext = add_extension_line(vg, Ptemp2, vector_ext, distance_factor=6)
    Ptemp2 = (Ptemp2 + 0.5*vector_ext)

    add_label_outside(vg, r'$h_{\rm head,r}$', Ptemp1, Ptemp2, ext_factor=3)



    vg.ax.get_xaxis().set_visible(False)
    vg.ax.get_yaxis().set_visible(False)
    vg.fig.tight_layout()
    # vg.fig.savefig(r'D:\OneDrive\[00]GetWorking\32 blimopti\p2019_ecce_bearingless_induction_full_paper\images\CAD_CrossSection.png')
    show()
    quit()

    # obsolete
    def perpendi_label(vg, label, PA, PB, distance=3, distance_factor=1, mirror=False):
        if mirror == True:
            PA = PA[0], -PA[1]
            PB = PB[0], -PB[1]
            distance *= -1

        new_PA = PA[0], PA[1]+distance + PB[1]-PA[1]
        new_PB = PB[0], PB[1]+distance 

        x = [ PA[0], new_PA[0] ]
        y = [ PA[1], new_PA[1] ]
        # x = [ PA[0] ]*2
        # y = [ PA[1], PA[1]+distance + PB[1]-PA[1] ]
        vg.ax.plot(x, y, 'k', lw=0.5)

        x = [ PB[0], new_PB[0] ]
        y = [ PB[1], new_PB[1] ]
        # x = [ PB[0] ]*2
        # y = [ PB[1], PB[1]+distance ]
        vg.ax.plot(x, y, 'k', lw=0.5)

        distance = sqrt((PA[0]-new_PA[0])**2 + (PA[1]-new_PA[1])**2)*distance_factor
        xy = array([0.5*(PA[0]+new_PA[0]), 0.5*(PA[1]+new_PA[1])])
        vector = array([new_PA[0]-PA[0], new_PA[1]-PA[1]])
        xytext = array([-vector[1], vector[0]]) / distance # perpendicular vector
        xytext = xy + xytext
        vg.ax.annotate(label, xytext=xytext, xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->"))
    def parallel_label(vg, label, PA, PB, vector_extendion=None, translation=0, rotation=0, distance_factor=2, mirror=False):

        if vector_extendion is None:
            pass
        else:
            if mirror == True:
                PA = PA[0], -PA[1]
                PB = PB[0], -PB[1]
            if translation == 0:
                new_PA = vg.park_transform((PA[0],PA[1]), rotation)
                new_PB = vg.park_transform((PB[0],PB[1]), rotation)
            elif translation == 1:
                new_PA = PA[0], PA[1]+translation
                new_PB = PB[0], PB[1]+translation
            elif translation == -1:
                new_PA = PA[0]+translation, PA[1]
                new_PB = PB[0]+translation, PB[1]

            distance = sqrt((PA[0]-new_PA[0])**2 + (PA[1]-new_PA[1])**2)
            if distance == 0: # rotation == 0
                distance = 1
            unit_vector = array([new_PA[0]-PA[0], new_PA[1]-PA[1]])/distance
            new_PA = array([PA[0] + unit_vector[0] * distance_factor, PA[1] + unit_vector[1] * distance_factor])
            new_PB = array([PB[0] + unit_vector[0] * distance_factor, PB[1] + unit_vector[1] * distance_factor]) 

            # A
            x = [ PA[0], new_PA[0] ]
            y = [ PA[1], new_PA[1] ]
            vg.ax.plot(x, y, 'k-', lw=0.6)

            # B
            x = [ PB[0], new_PB[0] ]
            y = [ PB[1], new_PB[1] ]
            vg.ax.plot(x, y, 'k-', lw=0.6)

            xy = array([0.5*(PA[0]+new_PA[0]), 0.5*(PA[1]+new_PA[1])])
            xytext = distance_factor * array([-unit_vector[1], unit_vector[0]]) # perpendicular vector
            xytext = xy + xytext
            vg.ax.annotate('', xytext=xytext, xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->"))
            vg.ax.text(xytext[0]+0.5, xytext[1]-0.5, label,
                        bbox=dict(facecolor='w', edgecolor='w',boxstyle='round,pad=0',alpha=1))



    # print '---Unit Test of VanGogh.'
    # # >>> c2 = Point(0.5,0).buffer(1).boundary
    # # >>> i = c.intersection(c2)
    # # >>> i.geoms[0].coords[0]
    # # (0.24999999999999994, -0.967031122079967)
    # # >>> i.geoms[1].coords[0]
    # # (0.24999999999999997, 0.967031122079967)

    # vg = VanGogh(None, 1)
    # from numpy import arccos
    # theta = arccos(0.25) * 2
    # print theta / pi * 180, 'deg'
    # print vg.find_center_of_a_circle_using_2_points_and_arc_angle((0.24999999999999994, -0.967031122079967),
    #                                                               (0.24999999999999994, 0.967031122079967),
    #                                                               theta)
    # print vg.find_center_of_a_circle_using_2_points_and_arc_angle((0.24999999999999994, 0.967031122079967),
    #                                                               (0.24999999999999994, -0.967031122079967),
    #                                                               theta)


    # quit()


    # import itertools
    # import matplotlib.patches as patches
    # import matplotlib.pyplot as plt
    # import numpy as np
    # import sys

    # fig, ax = plt.subplots()

    # npoints = 5

    # # Calculate the xy coords for each point on the circle
    # s = 2 * np.pi / npoints
    # verts = np.zeros((npoints, 2))
    # for i in np.arange(npoints):
    #     angle = s * i
    #     x = npoints * np.cos(angle)
    #     y = npoints * np.sin(angle)
    #     verts[i] = [x, y]

    # # Plot the Bezier curves
    # numbers = [i for i in xrange(npoints)]
    # bezier_path = np.arange(0, 1.01, 0.01)
    # for a, b in itertools.product(numbers, repeat=2):
    #     if a == b:
    #         continue

    #     x1y1 = x1, y1 = verts[a]
    #     x2y2 = x2, y2 = verts[b]

    #     xbyb = xb, yb = [0, 0]

    #     # Compute and store the Bezier curve points
    #     x = (1 - bezier_path)** 2 * x1 + 2 * (1 - bezier_path) * bezier_path * xb + bezier_path** 2 * x2
    #     y = (1 - bezier_path)** 2 * y1 + 2 * (1 - bezier_path) * bezier_path * yb + bezier_path** 2 * y2

    #     ax.plot(x, y, 'k-')

    # x, y = verts.T
    # ax.scatter(x, y, marker='o', s=50, c='r')

    # ax.set_xlim(-npoints - 5, npoints + 6)
    # ax.set_ylim(-npoints - 5, npoints + 6)
    # ax.set(aspect=1)
    # plt.show()
    # quit()



# winding diagram for four pole dpnv motor
if __name__ == '__main__':

    bool_4pole_2pole = False

    # generate the dict of list: dl
    if bool_4pole_2pole: # 4 pole motor
        # ExampleQ24p2m3ps1y6
        l41 = ['C', 'C', 'A', 'A', 'B', 'B', 'C', 'C', 'A', 'A', 'B', 'B', 'C', 'C', 'A', 'A', 'B', 'B', 'C', 'C', 'A', 'A', 'B', 'B']
        l42 = ['+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-', '+', '+', '-', '-']
        # Eric Severson [2:12 PM]
        # How the coils are assigned to groups A-D for Q=24 example winding
        #     ```if (~bParallel)
        #         MachineParams.Group_UA = [-19 -20];
        #         MachineParams.Group_UB = [1 2];
        #         MachineParams.Group_UC = [13 14];
        #         MachineParams.Group_UD = [-7 -8];
        #     else
        #         MachineParams.Group_UA = [-19 -20 13 14];
        #         MachineParams.Group_UB = [1 2 -7 -8];
        #     end```
        # Eric Severson [2:14 PM]
        # These numbers are the coil numbers. Negative indicates a coil connected in reverse.
        l_rightlayer1 = l41[::]
        l_rightlayer2 = l42[::]
        l_leftlayer1  = l41[::]
        l_leftlayer2  = l42[::]

    else: # 2 pole motor
        # ExampleQ24p1m3ps2y9
        # -w -w  u  u  u u -v -v -v -v  w  w  w w -u -u -u -u  v  v  v v w w
        #  v -w -w -w -w u  u  u  u -v -v -v -v w  w  w  w -u -u -u -u v v v
        l_rightlayer1 = ['C', 'C', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'C', 'C'][::-1]
        l_rightlayer2 = ['-', '-', '+', '+', '+', '+', '-', '-', '-', '-', '+', '+', '+', '+', '-', '-', '-', '-', '+', '+', '+', '+', '-', '-'][::-1]
        l_leftlayer1  = ['B', 'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'C', 'C', 'C', 'C', 'A', 'A', 'A', 'A', 'B', 'B', 'B'][::-1]
        l_leftlayer2  = ['+', '-', '-', '-', '-', '+', '+', '+', '+', '-', '-', '-', '-', '+', '+', '+', '+', '-', '-', '-', '-', '+', '+', '+'][::-1]

    dl_rightlayer = {'A': [], 'B': [], 'C':[]}
    dl_leftlayer  = {'A': [], 'B': [], 'C':[]}
    count_slot = 0
    while True:
        try:
            count_slot += 1

            ABC     = l_rightlayer1.pop(0)
            up_down = l_rightlayer2.pop(0)
            for target in ['A', 'B', 'C']:
                if ABC == target:
                    dl_rightlayer[target].append(up_down+str(count_slot))
                    break

            ABC     = l_leftlayer1.pop(0)
            up_down = l_leftlayer2.pop(0)
            for target in ['A', 'B', 'C']:
                if ABC == target:
                    dl_leftlayer[target].append(up_down+str(count_slot))
                    break

        except IndexError as e: # nothing to pop
            # raise e
            break

    # for k, v in dl_rightlayer.iteritems():
    #     print k,v
    # for k, v in dl_leftlayer.iteritems():
    #     print k,v
    # quit()

    FOUR = 1
    SIX = 3
    NINE = 4

    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    # Motor Spec
    #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
    Qs = 24
    if bool_4pole_2pole: # 4 pole motor: full pitch is 6
        coil_pitch = 6 
    else: # 2 pole motor: full pitch is 12
        coil_pitch = 9
    fig_xlim = (0.25,Qs+1)
    fig_ylim = (-6.25,10)

    from pylab import *
    plt.rcParams["font.family"] = "Times New Roman"

    coil_bias = 0.3 # 绘制双层绕组的第二层的偏移量
    END   = Qs + 1 - 1/3.
    SHIFT = END - 1 + 1/3.
    # END = Qs+coil_pitch/2+coil_bias/2
    print("END=", END)
    print("SHIFT=", SHIFT)

    # x means coil location (meaning is equivalent to loc)
    def draw_vertical_line_and_arrow_head(x, x2):

        coil_spacing = 1
        x *= coil_spacing
        ax.plot([abs(x), abs(x)], 
                [-FOUR, FOUR], '-', color=color, lw=1,alpha=1.0)
        ax.plot([abs(x2)+coil_bias, abs(x2)+coil_bias], 
                [-FOUR, FOUR], '--', color=color, lw=1,alpha=1.0)

        annotate_bias = 0.1
        if x>0:
            xytext = [abs(x), (-FOUR+FOUR)/2 + annotate_bias + 0.0]
            xy     = [abs(x), (-FOUR+FOUR)/2 - annotate_bias + 0.0]
        else:
            xy     = [abs(x), (-FOUR+FOUR)/2 + annotate_bias + 0.25]
            xytext = [abs(x), (-FOUR+FOUR)/2 - annotate_bias + 0.25]
        ax.annotate('', xytext=xytext, xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->", color=color, alpha=1.0))
        ax.text(xy[0]+1.25*coil_bias, annotate_bias + 0.25, '%d'%(abs(x)), fontsize=6)

        if x2>0:
            xytext = [abs(x2), (-FOUR+FOUR)/2 + annotate_bias + 0.0]
            xy     = [abs(x2), (-FOUR+FOUR)/2 - annotate_bias + 0.0]
        else:
            xy     = [abs(x2), (-FOUR+FOUR)/2 + annotate_bias + 0.9]
            xytext = [abs(x2), (-FOUR+FOUR)/2 - annotate_bias + 0.9]
        ax.annotate('', xytext=(xytext[0]+coil_bias, xytext[1]), xy=(xy[0]+coil_bias, xy[1]), xycoords='data', arrowprops=dict(arrowstyle="->", color=color, alpha=1.0))

    def draw_end_turns(x, coil_pitch):

        def mirror_plot(ax, X, Y, *arg, **kwarg):
            ax.plot(X, Y, *arg, **kwarg)
            ax.plot(X, [-el for el in Y], *arg, **kwarg)

        slope = (SIX-FOUR) / (coil_pitch/2. + coil_bias/2.)

        # going up
        loc1 = abs(x) + coil_pitch/2. + coil_bias/2. # half pitch
        # going down
        loc2 = abs(x) + coil_pitch + coil_bias # full pitch

        if loc1 <= END:
            # going up
            mirror_plot(ax,  [ abs(x),
                               loc1   ],
                             [ FOUR,SIX ], '-', color=color, lw=0.5,alpha=1.0)
            if loc2 < END:
                # going down 
                mirror_plot(ax,  [ loc1,
                                   loc2 ],
                                 [ SIX,FOUR ], '--', color=color, lw=0.5,alpha=1.0)
            else:
                # going down (to be continued)
                mirror_plot(ax,  [ loc1,
                                   END],
                                 [ SIX,SIX-slope*(END-loc1) ], '--', color=color, lw=0.5,alpha=1.0)
                # going down (continued)
                mirror_plot(ax,  [ END-SHIFT,
                                   loc2-SHIFT],
                                 [ SIX-slope*(END-loc1),FOUR ], '--', color=color, lw=0.5,alpha=1.0)
        else:
            # going up (to be continued)
            mirror_plot(ax,  [ abs(x),
                               END],
                             [ FOUR,FOUR+slope*(END-abs(x)) ], '-', color=color, lw=0.5,alpha=1.0)
            # going up (continued)
            mirror_plot(ax,  [ END-SHIFT,
                               loc1-SHIFT],
                             [ FOUR+slope*(END-abs(x)),SIX ], '-', color=color, lw=0.5,alpha=1.0)
            # going down 
            mirror_plot(ax,  [ loc1-SHIFT,
                               loc2-SHIFT ],
                             [ SIX,FOUR ], '--', color=color, lw=0.5,alpha=1.0)

        # if loc2 <= END:
        # else: # loc2 超了
        #     if loc1 <= END: # loc1 没超
        #         mirror_plot(ax,  [ loc1,
        #                            END],
        #                          [ 6,6-slope*(END-loc1) ], '--'+color, lw=0.5,alpha=1.0)
        #     else: # loc1 也超了
        #         print '-------------------------AAAAAA'
        #         # partial plot on the left 
        #         mirror_plot(ax,  [ loc1-SHIFT,
        #                            loc2-SHIFT],
        #                          [ 6,4 ], '--'+color, lw=0.5,alpha=1.0)

        return slope

    def draw_terminals(list_terminals, slope, phase='u'):
        # xx means (x,x) <=> (loc, loc)
        for xx, letter in zip(list_terminals, ['b', 'a', 'd', 'c']): # ['a', 'b', 'c', 'd'] changed. 6/26/2019
            print(xx, letter)
            # ua+
            loc_plus = abs(xx[0])+coil_pitch/2+coil_bias/2 # top of the triangle
            if xx[0]>0: # solid line
                if loc_plus + 0.3 > END:
                    loc_plus -= SHIFT
                c1 = (loc_plus - 0.3, NINE)
                y_span = [c1[1]]+[SIX-0.3*slope]
                xy     = [c1[0], sum(y_span)/2]
                xytext = [c1[0], sum(y_span)/2 + 0.25]
                line_style = '-'
            else: # dashed line
                if loc_plus + 0.3 > END:
                    loc_plus -= SHIFT
                c1 = (loc_plus + 0.3, NINE)
                y_span = [c1[1]]+[SIX-0.3*slope]
                xy     = [c1[0], sum(y_span)/2]
                xytext = [c1[0], sum(y_span)/2 + 0.25]
                line_style = '--'
            ax.plot([c1[0]]*2, y_span, line_style, lw=0.5, alpha=1.0, color=color)
            ax.annotate('', xytext=xytext, xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->", color=color, lw=0.5, alpha=1.0))

            # ua-
            loc_minus = abs(xx[1])+coil_pitch/2+coil_bias/2 # top of the triangle
            if xx[1]>0: # solid line
                if loc_minus + 0.3 > END:
                    loc_minus -= SHIFT
                c2 = (loc_minus - 0.3, NINE)
                y_span = [c2[1]]+[SIX-0.3*slope]
                xytext = [c2[0], sum(y_span)/2]
                xy     = [c2[0], sum(y_span)/2 + 0.25]
                line_style = '-'
            else: # dashed line
                if loc_minus + 0.3 > END:
                    loc_minus -= SHIFT
                c2 = (loc_minus + 0.3, NINE)
                y_span = [c2[1]]+[SIX-0.3*slope]
                xytext = [c2[0], sum(y_span)/2]
                xy     = [c2[0], sum(y_span)/2 + 0.25]
                line_style = '--'
            ax.plot([c2[0]]*2, y_span, line_style, lw=0.5, alpha=1.0, color=color)
            ax.annotate('', xytext=xytext, xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->", color=color, lw=0.5, alpha=1.0))

            # circle and ua+/- labels
            terminal_symbols = (r'$%s_{%s}^+$'%(phase, letter), r'$%s_{%s}^-$'%(phase, letter))
            ax.text(c1[0]-0.5, c1[1]+0.5, terminal_symbols[0], fontsize=6, color=color)
            ax.text(c2[0]-0.5, c2[1]+0.5, terminal_symbols[1], fontsize=6, color=color)
            circle1 = plt.Circle(c1, 0.1, edgecolor=color, facecolor='w', lw=0.5, zorder=10)
            circle2 = plt.Circle(c2, 0.1, edgecolor=color, facecolor='w', lw=0.5, zorder=10)
            ax.add_artist(circle1)
            ax.add_artist(circle2)

    def draw_coil_connections(list_coil_connection, slope):

        # xx means (x,x) <=> (loc, loc)
        for xx in list_coil_connection:
            print(xx)
            loc1 = abs(xx[0])+coil_pitch/2+coil_bias/2 # top of the triangle
            loc2 = abs(xx[1])+coil_pitch/2+coil_bias/2 # top of the triangle
            flag_shift = False
            if loc1 > END:
                loc1 -= SHIFT
                flag_shift = not flag_shift
            if loc2 > END:
                loc2 -= SHIFT
                flag_shift = not flag_shift # two shifts are equivalent to no shift ~

            flag_swapped = False
            if loc1 > loc2:
                loc1, loc2 = loc2, loc1 # make sure that loc1/loc2 are the left/right location, respectively
                flag_swapped = True

            bias = 0.25
            # horizontal
            if flag_shift:
                bias *= -1
                # find the shortest connection (if shifted)
                ax.plot([loc2-bias, END], [SIX+.25, SIX+.25], '-', color=color, lw=0.5, alpha=1.0)
                ax.plot([END-SHIFT, loc1+bias], [SIX+.25, SIX+.25], '-', color=color, lw=0.5, alpha=1.0)
            else:
                # regular 
                ax.plot([loc1+bias, loc2-bias], [SIX+.25, SIX+.25], '-', color=color, lw=0.5, alpha=1.0)

            # vertical 1
            line_style = '-'
            if flag_swapped != (xx[0]<0): # XOR
                line_style = '-' + line_style
            ax.plot([loc1+bias]*2, [SIX+.25, SIX-abs(bias)*slope], line_style, color=color, lw=0.5, alpha=1.0)

            # vertical 2
            line_style = '-'
            if flag_swapped != (xx[1]<0): # XOR
                line_style = '-' + line_style
            ax.plot([loc2-bias]*2, [SIX+.25, SIX-abs(bias)*slope], line_style, color=color, lw=0.5, alpha=1.0)

    for phase in ['u', 'v', 'w']:
        fig = figure(dpi=300)
        ax = fig.add_subplot(111, aspect='equal') # fixed aspect ratio such that a circle is still a circle when you drag the tk plot window
        ax.axis('off')

        # draw coils in the slots as vertical lines
        for color, winding_layout_right, winding_layout_left, _ in zip( ['#5C9C31', '#E97675', '#6C8FAB'] , 
                                                                        [dl_rightlayer['A'], dl_rightlayer['B'], dl_rightlayer['C']],
                                                                        [dl_leftlayer['A'], dl_leftlayer['B'], dl_leftlayer['C']],
                                                                        list(range(3))):
            # if color != 'b':
            #     continue
            for el1, el2 in zip(winding_layout_right, winding_layout_left):
                x = float(el1)
                x2 = float(el2)
                print(_, color, x)
                draw_vertical_line_and_arrow_head(x, x2) # right layer, left layer
                slope = draw_end_turns(x, coil_pitch)
                # show()

        if bool_4pole_2pole: # 4 pole motor
            if phase == 'w': # Swap u and w. --6/26/2019
                # The key to determine this terminal location is to stick to the dual purpose feature 
                # meanwhile considering the fact that ub and ud will not change current direction, while ua and uc will change when under suspension excitation.
                # Minus sign means dashed line (the second/lower layer of the double layer winding),
                # while the first number of the tuple means the current is entering, and the second number of the tuple means the current is leaving.
                # E.g., for ub, current enters 4 and leaves 3. Particularly, the current first enters the lower layer of slot 4, hence -4.
                # Coil Group:     ua                       ub                           uc                         ud                         (see Severson15Dual@Fig.4b)
                color = '#5C9C31'
                list_terminals = [(9, -10),                (-4,  3),                    (-16, 15),                 (21, -22)] 
            if phase == 'v': 
                color = '#E97675'
                list_terminals = [(9+Qs/3, -(10+Qs/3)),    (-(4+Qs/3),     3+Qs/3),     (-(16+Qs/3), 15+Qs/3),     (21+Qs/3-Qs, -(22+Qs/3-Qs))]
            if phase == 'u': 
                color = '#6C8FAB'
                list_terminals = [(9-Qs/3, -(10-Qs/3)),    (-(4-Qs/3+Qs),  3-Qs/3+Qs),  (-(16-Qs/3), 15-Qs/3),     (21-Qs/3, -(22-Qs/3))]

        else: # 2 pole motor
            if phase == 'w': # Swap u and w. --7/15/2019
                color = '#5C9C31'
                list_terminals = [(19, -20),                (-10, 9),                    (21,  -22),                 (-8, 7)] # we have swapped the order of ub and ud
            if phase == 'v': 
                color = '#E97675'
                list_terminals = [(19+Qs/3, -(20+Qs/3)),    (-(10+Qs/3), 9+Qs/3),        (21+Qs/3,  -(22+Qs/3)),     (-(8+Qs/3), 7+Qs/3)] # we have swapped the order of ub and ud
            if phase == 'u': 
                color = '#6C8FAB'
                list_terminals = [(19-Qs/3, -(20-Qs/3)),    (-(10-Qs/3), 9-Qs/3),        (21-Qs/3,  -(22-Qs/3)),     (-(8-Qs/3+Qs), 7-Qs/3+Qs)] # we have swapped the order of ub and ud

        # draw terminals
        draw_terminals(list_terminals, slope, phase=phase)

        # draw connection lines between coil and coil
        list_coil_connection = [(-el[0], -el[1]) for el in list_terminals]
        if Qs!=24:
            raise Exception('this pattern only works for winding whose coil group has only two coils.')
        draw_coil_connections(list_coil_connection, slope)

        scatter(fig_xlim[0], fig_ylim[1], color='w') # 不然text要超出画框了
        ax.set_xlim(fig_xlim)
        ax.set_ylim(fig_ylim)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        fig.tight_layout()
        fig.savefig(r'D:\OneDrive - UW-Madison\c\release\winding_diagram_phase_%s.png'%(phase), dpi=300)


    # ROTOR WINDING (four pole 8 phase single layer wave winding)
    # ROTOR WINDING (four pole 8 phase single layer wave winding)
    # ROTOR WINDING (four pole 8 phase single layer wave winding)
    class one_phase_winding(object):
        def __init__(self, ax, list_slot_location, coil_pitch, list_labels=None, color='k'):
            self.ax = ax
            self.list_slot_location = list_slot_location
            self.coil_pitch = coil_pitch
            self.color=color
            if list_labels is None:
                self.list_labels = [abs(x) for x in list_slot_location]
            else:
                self.list_labels = list_labels

            self.draw_vertical_line()
            self.draw_arrow()
            self.draw_end_turns()

        def draw_vertical_line(self):
            for x in self.list_slot_location:
                self.ax.plot(   [abs(x), abs(x)], 
                                [-FOUR, FOUR], '-', color=self.color, lw=1,alpha=1.0)
                # ax.plot([abs(x2)+coil_bias, abs(x2)+coil_bias], 
                #         [-FOUR, FOUR], '--'+color, lw=1,alpha=1.0)

        def draw_arrow(self):
            annotate_bias = 0.1
            text_bias = 0.05
            for x, label in zip(self.list_slot_location, self.list_labels):
                if x>0:
                    xytext = [abs(x), (-FOUR+FOUR)/2 + annotate_bias + 0.0]
                    xy     = [abs(x), (-FOUR+FOUR)/2 - annotate_bias + 0.0]
                else:
                    xy     = [abs(x), (-FOUR+FOUR)/2 + annotate_bias + text_bias]
                    xytext = [abs(x), (-FOUR+FOUR)/2 - annotate_bias + text_bias]
                ax.annotate('', xytext=xytext, xy=xy, xycoords='data', arrowprops=dict(arrowstyle="->", color=self.color, alpha=1.0))
                ax.text(xy[0]+text_bias, annotate_bias, '%d'%(label), fontsize=24)

                # if x2>0:
                #     xytext = [abs(x2), (-FOUR+FOUR)/2 + annotate_bias + 0.0]
                #     xy     = [abs(x2), (-FOUR+FOUR)/2 - annotate_bias + 0.0]
                # else:
                #     xy     = [abs(x2), (-FOUR+FOUR)/2 + annotate_bias + 0.9]
                #     xytext = [abs(x2), (-FOUR+FOUR)/2 - annotate_bias + 0.9]
                # ax.annotate('', xytext=(xytext[0]+coil_bias, xytext[1]), xy=(xy[0]+coil_bias, xy[1]), xycoords='data', arrowprops=dict(arrowstyle="->", color=color, alpha=1.0))

        def draw_end_turns(self):

            count = 0
            for x, next_x in zip(self.list_slot_location, self.list_slot_location[1:]):
                count += 1
                loc1 = abs(x) + self.coil_pitch/2

                if count%2==1:
                    self.ax.plot(   [abs(x), loc1], 
                                    [FOUR, SIX], '-', color=self.color, lw=0.5, alpha=1.0)
                    self.ax.plot(   [loc1, abs(next_x)], 
                                    [SIX, FOUR], '-', color=self.color, lw=0.5, alpha=1.0)
                else:
                    self.ax.plot(   [abs(x), loc1], 
                                    [-FOUR, -SIX], '-', color=self.color, lw=0.5, alpha=1.0)
                    self.ax.plot(   [loc1, abs(next_x)], 
                                    [-SIX, -FOUR], '-', color=self.color, lw=0.5, alpha=1.0)

    FOUR = 0.5
    SIX = 0.75

    fig = figure(dpi=300)
    ax = fig.add_subplot(111, aspect='equal') # fixed aspect ratio such that a circle is still a circle when you drag the tk plot window
    ax.axis('off')

    if False:
        list_all_phase = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']*4
        list_all_direction = ['+', '+', '+', '+', '+', '+', '+', '+', '-', '-', '-', '-', '-', '-', '-', '-']*2

        list_phase_name = []
        count = 0
        for phase, direction in zip(list_all_phase, list_all_direction):
            count += 1
            if phase not in list_phase_name:
                list_phase_name.append(phase)
                exec(phase + ' = []')
                exec(phase + '_direction = []')
            exec(phase + '.append(count)')
            temp = 1 if list_all_direction[count-1] == '+' else -1
            exec(phase + '_direction.append(temp)')

        for phase in list_phase_name:
            print(eval(f'{phase}'))
            print(eval(f'{phase}_direction'))
            
            one_phase_winding(ax, [number*direction for number, direction in zip(eval(f'{phase}'), eval(f'{phase}_direction'))])
    else:
        list_labels = [1, 9, 17, 26]
        one_phase_winding(ax, [1, -2, 3, -4], 1, list_labels)

    ax.plot( [ 1, 1], 
             [-FOUR, -SIX*1.25], '-', color='k', lw=0.5, alpha=1.0)
    ax.plot( [ 1, 4], 
             [-SIX*1.25, -SIX*1.25], '-', color='k', lw=0.5, alpha=1.0)
    ax.plot( [ 4, 4], 
             [-SIX*1.25, -FOUR], '-', color='k', lw=0.5, alpha=1.0)

    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.tight_layout()
    fig.savefig(r'D:\OneDrive - UW-Madison\c\release\rotor_winding.png', dpi=300)

    show()
    quit()

# test tangent points of circles
if __name__ == '__main__':

    c1 = (-48,0)
    c2 = (-40,0)
    c2 = c1
    loc1, loc2 = get_tangent_points_of_two_circles(c1, 5, c2, 5)
    from pylab import *
    fig = figure()        
    ax = gcf().gca() # get current figure/axis
    ax = fig.add_subplot(111, aspect='equal')
    scatter(*loc1)
    scatter(*loc2)
    scatter(*c1)
    scatter(*c2)
    circle1 = plt.Circle(c1, 5, color='r')
    circle2 = plt.Circle(c2, 5, color='r')
    ax.add_artist(circle1)
    ax.add_artist(circle2)
    ax.plot([loc1[0], loc2[0]], [loc1[1], loc2[1]], 'k')
    show()
    quit()

