from pylab import np, cos, sin, arctan
class CrossSectInnerRotorStator:
    # CrossSectInnerRotorPMStator Describes the inner rotor PM stator.
    #    Properties are set upon class creation and cannot be modified.
    #    The anchor point for this is the center of the stator,
    #    with the x-axis directed down the center of one of the stator teeth.
    def __init__(self,
                    name = 'StatorCore',
                    color = '#BAFD01',
                    deg_alpha_st = 40, # span angle of tooth: class type DimAngular
                    deg_alpha_so = 20, # angle of tooth edge: class type DimAngular
                    mm_r_si      = 40,# inner radius of stator teeth: class type DimLinear
                    mm_d_so      = 5,# tooth edge length: class type DimLinear
                    mm_d_sp      = 10,# tooth tip length: class type DimLinear
                    mm_d_st      = 15,# tooth base length: class type DimLinear
                    mm_d_sy      = 15,# back iron thickness: class type DimLinear
                    mm_w_st      = 13,# tooth base width: class type DimLinear
                    mm_r_st      = 0,# fillet on outter tooth: class type DimLinear
                    mm_r_sf      = 0,# fillet between tooth tip and base: class type DimLinear
                    mm_r_sb      = 0,# fillet at tooth base: class type DimLinear
                    Q = 6,            # number of stator slots (integer)
                    location = None
                ):

        self.name = name
        self.color = color
        self.deg_alpha_st = deg_alpha_st
        self.deg_alpha_so = deg_alpha_so
        self.mm_r_si      = mm_r_si     
        self.mm_d_so      = mm_d_so     
        self.mm_d_sp      = mm_d_sp     
        self.mm_d_st      = mm_d_st     
        self.mm_d_sy      = mm_d_sy     
        self.mm_w_st      = mm_w_st     
        self.mm_r_st      = mm_r_st     
        self.mm_r_sf      = mm_r_sf     
        self.mm_r_sb      = mm_r_sb     
        self.Q = Q               
        self.location = location 

    def draw(self, drawer):

        drawer.getSketch(self.name, self.color)

        alpha_st = self.deg_alpha_st * np.pi/180
        alpha_so = -self.deg_alpha_so * np.pi/180
        r_si = self.mm_r_si
        d_so = self.mm_d_so
        d_sp = self.mm_d_sp
        d_st = self.mm_d_st
        d_sy = self.mm_d_sy
        w_st = self.mm_w_st
        r_st = self.mm_r_st
        r_sf = self.mm_r_sf
        r_sb = self.mm_r_sb
        Q    = self.Q

        alpha_slot_span = 360/Q * np.pi/180

        P1 = [r_si, 0]
        P2 = [r_si*cos(alpha_st*0.5), r_si*-sin(alpha_st*0.5)]
        P3_temp = [ d_so*cos(alpha_st*0.5), 
                    d_so*-sin(alpha_st*0.5)]
        P3_local_rotate = [  cos(alpha_so)*P3_temp[0] + sin(alpha_so)*P3_temp[1],
                             -sin(alpha_so)*P3_temp[0] + cos(alpha_so)*P3_temp[1] ]
        P3 = [  P3_local_rotate[0] + P2[0],
                P3_local_rotate[1] + P2[1] ]

        三角形的底 = r_si + d_sp
        三角形的高 = w_st*0.5
        三角形的角度 = arctan(三角形的高 / 三角形的底)
        P4 = [  三角形的底*cos(三角形的角度), 
                三角形的底*-sin(三角形的角度)]

        P5 = [ P4[0] + d_st, 
               P4[1]]

        P6 = [ (r_si+d_sp+d_st)*cos(alpha_slot_span*0.5),
               (r_si+d_sp+d_st)*-sin(alpha_slot_span*0.5) ]

        P7 = [ (r_si+d_sp+d_st+d_sy)*cos(alpha_slot_span*0.5),
               (r_si+d_sp+d_st+d_sy)*-sin(alpha_slot_span*0.5) ]
        P8 = [  r_si+d_sp+d_st+d_sy, 0]

        list_segments = []
        list_segments += drawer.drawArc([0,0], P2, P1)
        list_segments += drawer.drawLine(P2, P3)
        list_segments += drawer.drawLine(P3, P4)
        list_segments += drawer.drawLine(P4, P5)
        list_segments += drawer.drawArc([0,0], P6, P5)
        list_segments += drawer.drawLine(P6, P7)
        # l, vA = drawer.drawArc([0,0], P6, P5, returnVertexName=True)
        # list_segments += l
        # l, vB = drawer.drawLine(P6, P7, returnVertexName=True)
        # list_segments += l
        # drawer.addConstraintCocentricity(vA[0], vB[0])
        # raise

        list_segments += drawer.drawArc([0,0], P7, P8)
        list_segments += drawer.drawLine(P8, P1)

        return [list_segments]

def get_area_polygon(a,b,c,d):
    x1, x2, x3, x4 = a[0], b[0], c[0], d[0]
    y1, y2, y3, y4 = a[1], b[1], c[1], d[1]

    return 0.5*abs( x1*y2-y1*x2 + x2*y3-y2*x3 + x3*y4-y3*x4 + x4*y1-y4*x1 )

class CrossSectInnerRotorStatorWinding(object):
    def __init__(self, 
                    name = 'Coils',
                    color = '#3D9970',
                    stator_core = None
                    ):
        self.name = name
        self.color = color
        self.stator_core = stator_core

    def draw(self, drawer, bool_re_evaluate=False):

        if False == bool_re_evaluate:
            drawer.getSketch(self.name, self.color)

        alpha_st = self.stator_core.deg_alpha_st * np.pi/180
        alpha_so = self.stator_core.deg_alpha_so * np.pi/180
        r_si     = self.stator_core.mm_r_si
        d_so     = self.stator_core.mm_d_so
        d_sp     = self.stator_core.mm_d_sp
        d_st     = self.stator_core.mm_d_st
        d_sy     = self.stator_core.mm_d_sy
        w_st     = self.stator_core.mm_w_st
        r_st     = self.stator_core.mm_r_st
        r_sf     = self.stator_core.mm_r_sf
        r_sb     = self.stator_core.mm_r_sb
        Q        = self.stator_core.Q

        alpha_slot_span = 360/Q * np.pi/180

        P1 = [r_si, 0]

            # 乘以0.99或0.95避免上层导体和下层导体重合导致导入Designer时产生多余的Parts。
        POpen = [(r_si+d_sp)*cos(alpha_slot_span*0.5*1.00), (r_si+d_sp)*-sin(alpha_slot_span*0.5*1.00)]

        # P2 = [r_si*cos(alpha_st*0.5), r_si*-sin(alpha_st*0.5)]

        # P3_temp = [ d_so*cos(alpha_st*0.5), 
        #             d_so*-sin(alpha_st*0.5)]
        # P3_local_rotate = [  cos(alpha_so)*P3_temp[0] + sin(alpha_so)*P3_temp[1],
        #                      -sin(alpha_so)*P3_temp[0] + cos(alpha_so)*P3_temp[1] ]
        # P3 = [  P3_local_rotate[0] + P2[0],
        #         P3_local_rotate[1] + P2[1] ]

        三角形的底 = r_si + d_sp
        三角形的高 = w_st*0.5
        三角形的角度 = arctan(三角形的高 / 三角形的底)
        P4 = [  三角形的底*cos(三角形的角度), 
                三角形的底*-sin(三角形的角度) ]

        P5 = [ P4[0] + d_st, 
               P4[1]]

        PMiddle45 = [0.5*(P4[0] + P5[0]), P4[1]]
        TheRadius = (P5[0] - P4[0])*0.45

            # 为了使得槽和导体之间不要接触，试着添加5%的clearance？
        P6 = [ (r_si+d_sp+d_st)*cos(alpha_slot_span*0.5) *1.00,
               (r_si+d_sp+d_st)*-sin(alpha_slot_span*0.5) *1.00 ]

        self.mm2_slot_area = 2 * get_area_polygon(P4, P5, P6, POpen)
        print('Slot area is %g mm^2'%(self.mm2_slot_area))
        if bool_re_evaluate:
            return self.mm2_slot_area

        PMiddle6Open = [ 0.5*(P6[0]+POpen[0]), 0.5*(P6[1]+POpen[1])]
        self.PCoil = PCoil = [ 0.5*(PMiddle45[0]+PMiddle6Open[0]), 0.5*(PMiddle45[1]+PMiddle6Open[1])]

        # P7 = [ (r_si+d_sp+d_st+d_sy)*cos(alpha_slot_span*0.5),
        #        (r_si+d_sp+d_st+d_sy)*-sin(alpha_slot_span*0.5) ]
        # P8 = [  r_si+d_sp+d_st+d_sy, 0]

        # Compute the vector starting from PCoil to one of the corner of the polygon.
        def shrink(PC, P):
            vector = [ P[0] - PC[0], P[1] - PC[1]]
            return [ PC[0]+0.95*vector[0], PC[1]+0.95*vector[1] ]
        P6_Shrink = shrink(PCoil, P6)
        P5_Shrink = shrink(PCoil, P5)
        P4_Shrink = shrink(PCoil, P4)
        POpen_Shrink = shrink(PCoil, POpen)

        list_regions = []

        list_segments = []
        # list_segments += drawer.drawCircle(PCoil, TheRadius)
        list_segments += drawer.drawArc([0,0], P6_Shrink, P5_Shrink)
        list_segments += drawer.drawLine(P5_Shrink, P4_Shrink)
        list_segments += drawer.drawLine(P4_Shrink, POpen_Shrink)
        list_segments += drawer.drawLine(POpen_Shrink, P6_Shrink)
        list_regions.append(list_segments)

        PCoil[1] *= -1
        P6_Shrink[1] *= -1
        P5_Shrink[1] *= -1
        P4_Shrink[1] *= -1
        POpen_Shrink[1] *= -1
        list_segments = []
        # list_segments += drawer.drawCircle(PCoil, TheRadius)
        list_segments += drawer.drawArc([0,0], P5_Shrink, P6_Shrink)
        list_segments += drawer.drawLine(P5_Shrink, P4_Shrink)
        list_segments += drawer.drawLine(P4_Shrink, POpen_Shrink)
        list_segments += drawer.drawLine(POpen_Shrink, P6_Shrink)
        list_regions.append(list_segments)
        list_segments = []

        return list_regions

if __name__ == '__main__':
    import JMAG
    import Location2D
    if True:
        from utility import my_execfile
        my_execfile('./default_setting.py', g=globals(), l=locals())
        fea_config_dict

        toolJd = JMAG.JMAG(fea_config_dict)

        project_name          = 'proj%d'%(0)
        expected_project_file_path = './' + "%s.jproj"%(project_name)
        toolJd.open(expected_project_file_path)

    stator_core = CrossSectInnerRotorStator( name = 'StatorCore',
                                        deg_alpha_st = 40,
                                        deg_alpha_so = 20,
                                        mm_r_si = 40,
                                        mm_d_so = 5,
                                        mm_d_sp = 10,
                                        mm_d_st = 15,
                                        mm_d_sy = 15,
                                        mm_w_st = 13,
                                        mm_r_st = 0,
                                        mm_r_sf = 0,
                                        mm_r_sb = 0,
                                        Q = 6,
                                        location = Location2D.Location2D(anchor_xy=[0,0], deg_theta=0)
                                        )

    list_regions = stator_core.draw(toolJd)
    toolJd.bMirror = True
    toolJd.iRotateCopy = stator_core.Q
    region1 = toolJd.prepareSection(list_regions)

    coils = CrossSectInnerRotorStatorWinding(name = 'Coils',
                                                stator_core = stator_core)

    list_regions = coils.draw(toolJd)
    toolJd.bMirror = False
    toolJd.iRotateCopy = coils.stator_core.Q
    region2 = toolJd.prepareSection(list_regions)

    # Import Model into Designer
    toolJd.doc.SaveModel(False) # True: Project is also saved. 
    model = toolJd.app.GetCurrentModel()
    model.SetName('BPMSM Modeling')
    model.SetDescription('BPMSM Test')



