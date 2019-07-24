 #coding:u8
from pylab import *
mu0 = (4*pi*1e-7)

# fname = "../Arnon5/Arnon-5-NGOES-Magnetization-Curve_Redo.csv"
# data = loadtxt(fname, delimiter=',', skiprows=1) # https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.loadtxt.html
# x = data[:,0] * 1000/(4*np.pi) # Oesteds to A/m
# y = data[:,1] * 1e-4 # Gauss to Tesla

# fname = "../Arnon5/M-19-Steel-BH-Curve.txt"
# y, x = loadtxt(fname, unpack=True, usecols=(0,1))

fname = "../Arnon5/M-19-Steel-BH-Curve-afterJMAGsmooth.txt"
y, x = loadtxt(fname, unpack=True, usecols=(1,0))
y = y[1:] # skip the (0,0) point will be essential for pade approximate to work
x = x[1:]
print(fname)

def rational(x, p, q): # pade appximate
    """ https://stackoverflow.com/questions/29815094/rational-function-curve-fitting-in-python
    The general rational function description.
    p is a list with the polynomial coefficients in the numerator
    q is a list with the polynomial coefficients (except the first one)
    in the denominator
    The zeroth order coefficient of the denominator polynomial is fixed at 1.
    Numpy stores coefficients in [x**2 + x + 1] order, so the fixed
    zeroth order denominator coefficent must comes last. (Edited.)
    """
    return np.polyval(p, x) / np.polyval(q + [1.0], x)

def rational8_8(x, p0, p1, p2, p3, p4, p5, p6, p7, q1, q2, q3, q4, q5, q6, q7):
    '''88'''
    return rational(x, [p0, p1, p2, p3, p4, p5, p6, p7], [q1, q2, q3, q4, q5, q6, q7])
def rational8_7(x, p0, p1, p2, p3, p4, p5, p6, p7, q1, q2, q3, q4, q5, q6):
    '''87'''
    return rational(x, [p0, p1, p2, p3, p4, p5, p6, p7], [q1, q2, q3, q4, q5, q6])
def rational7_7(x, p0, p1, p2, p3, p4, p5, p6, q1, q2, q3, q4, q5, q6):
    '''77'''
    return rational(x, [p0, p1, p2, p3, p4, p5, p6], [q1, q2, q3, q4, q5, q6])
def rational7_6(x, p0, p1, p2, p3, p4, p5, p6, q1, q2, q3, q4, q5):
    '''76'''
    return rational(x, [p0, p1, p2, p3, p4, p5, p6], [q1, q2, q3, q4, q5])
def rational7_5(x, p0, p1, p2, p3, p4, p5, p6, q1, q2, q3, q4):
    '''75'''
    return rational(x, [p0, p1, p2, p3, p4, p5, p6], [q1, q2, q3, q4])
def rational6_6(x, p0, p1, p2, p3, p4, p5, q1, q2, q3, q4, q5):
    '''66'''
    return rational(x, [p0, p1, p2, p3, p4, p5], [q1, q2, q3, q4, q5])
def rational6_5(x, p0, p1, p2, p3, p4, p5, q1, q2, q3, q4):
    '''65'''
    return rational(x, [p0, p1, p2, p3, p4, p5], [q1, q2, q3, q4])
def rational6_4(x, p0, p1, p2, p3, p4, p5, q1, q2, q3):
    '''64'''
    return rational(x, [p0, p1, p2, p3, p4, p5], [q1, q2, q3])
def rational5_5(x, p0, p1, p2, p3, p4, q1, q2, q3, q4):
    '''55'''
    return rational(x, [p0, p1, p2, p3, p4], [q1, q2, q3, q4])
def rational5_3(x, p0, p1, p2, p3, p4, q1, q2):
    '''53'''
    return rational(x, [p0, p1, p2, p3, p4], [q1, q2])
def rational4_4(x, p0, p1, p2, p3, q1, q2, q3):
    '''44'''
    return rational(x, [p0, p1, p2, p3], [q1, q2, q3])
def rational4_3(x, p0, p1, p2, p3, q1, q2):
    '''43'''
    return rational(x, [p0, p1, p2, p3], [q1, q2])
def rational3_3(x, p0, p1, p2, q1, q2):
    '''33'''
    return rational(x, [p0, p1, p2], [q1, q2])
# x = np.linspace(0, 10, 100)  
# y = rational(x, [-0.2, 0.3, 0.5], [-1.0, 2.0])
# ynoise = y * (1.0 + np.random.normal(scale=0.1, size=x.shape))

ax1 = figure().gca()
ax2 = figure().gca()

from scipy.optimize import curve_fit
global y_fit
def fitting_test(ax1, ax2, x, y, rational_func, label):
    global y_fit

    popt, pcov = curve_fit(rational_func, x, y) #p0=(0.2, 0.3, 0.5, 1, -1.0, 2.0, 1), diag=(1./x.mean(),1./y.mean()))
    print('\n\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print(label)
    print(popt)
    y_fit = rational_func(x, *popt)

    ax1.plot(x, y_fit, label='fit'+label, alpha=0.33)
    ax2.plot(x, (y_fit - y) / y, label=label+'normalized fit error', alpha=0.33)

    print('accumulated error', sum((y_fit - y)/y)) 

    return y_fit

''' BH curve and its pade apprximate [no (0,0) allowed]
'''
# fitting_test(ax1, ax2, x, y, rational8_7, rational8_7.__doc__)
# fitting_test(ax1, ax2, x, y, rational7_7, rational7_7.__doc__)
fitting_test(ax1, ax2, x, y, rational7_6, rational7_6.__doc__) # after visual inspection, this is the best
# fitting_test(ax1, ax2, x, y, rational6_6, rational6_6.__doc__)
# fitting_test(ax1, ax2, x, y, rational6_5, rational6_5.__doc__)
# fitting_test(ax1, ax2, x, y, rational6_4, rational6_4.__doc__)
# fitting_test(ax1, ax2, x, y, rational5_5, rational5_5.__doc__)
# fitting_test(ax1, ax2, x, y, rational5_3, rational5_3.__doc__)
# fitting_test(ax1, ax2, x, y, rational4_4, rational4_4.__doc__)
# fitting_test(ax1, ax2, x, y, rational4_3, rational4_3.__doc__)
# fitting_test(ax1, ax2, x, y, rational3_3, rational3_3.__doc__)
# fitting_test(ax1, ax2, x, y, rational8_8, rational8_8.__doc__)

ax1.plot(x, y, '*', label='ori', alpha=0.33)
ax1.grid()
ax1.legend()
ax2.grid()
ax2.legend()

''' Dynamic Gradient/mu0 to H
'''
ax3 = figure().gca()
plot(x,     gradient(y) / gradient(x) / mu0, '*', label='ori', alpha=0.33)
plot(x, gradient(y_fit) / gradient(x) / mu0, 'o-', label='fit', alpha=0.33)
ax3.grid()
ax3.legend()

''' mu-H curve
'''
ax4 = figure().gca()
plot(x, y/x, '*', label='ori')
plot(x, y_fit/x, '^', label='fit')

ax4.grid()
ax4.legend()

''' B-mu_r curve
'''
ax5 = figure().gca()
plot(y,         y/x / mu0, '*', label='ori')
plot(y_fit, y_fit/x / mu0, '^', label='fit')

ax5.grid()
ax5.legend()

''' Staiti Gradient (What you see in JMAG)
'''
ax6 = figure().gca()
plot(x, gradient(y) / gradient(x), '*', label='ori', alpha=0.33)
plot(x, gradient(y_fit) / gradient(x), 'o-', label='fit', alpha=0.33)
ax6.grid()
ax6.legend()

''' Dynamic Gradient/mu0 to B
'''
ax7 = figure().gca()
plot(    y,     gradient(y) / gradient(x) / mu0, '*', label='ori', alpha=0.33)
plot(y_fit, gradient(y_fit) / gradient(x) / mu0, 'o-', label='fit', alpha=0.33)
ax7.grid()
ax7.legend()



show()

# from scipy.integrate import cumtrapz, simps
# def get_end_H(xM19, x): # xM19的数据足够长，而x的不够，
#     for ind, el in enumerate(xM19):
#         if el > x[-1]:
#             print el, x[-1], 'A/m'
#             print 'index is', ind
#             return ind
# M19_gradient = gradient(yM19) / gradient(xM19)
# # M19 appended to original
# end_index = get_end_H(xM19, x)
# y_fit_gradient_integrate = cumtrapz(M19_gradient[end_index:], xM19[end_index:], initial=0) + y[-1]
# ax1.plot(xM19[end_index:], y_fit_gradient_integrate, 'v-', label='M19-ori', alpha=0.33)



# from scipy import misc
# e_exp = [1.0, 1.0, 1.0/2.0, 1.0/6.0, 1.0/24.0, 1.0/120.0]
# p, q = misc.pade(e_exp, 2)
# '这里的misc.pade需要已知近似对象的泰勒级数，呵呵'

