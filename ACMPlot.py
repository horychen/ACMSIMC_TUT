#coding:u8
from pylab import plt, mpl, arange, show
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
# from pprint import pprint
from collections import OrderedDict as O
# plot style
plt.style.use('ggplot') 
# plt.style.use('grayscale') # print plt.style.available # get [u'dark_background', u'bmh', u'grayscale', u'ggplot', u'fivethirtyeight']
# plot setting
mpl.rcParams['legend.fontsize'] = 14
# fontdict
font = {'family' : 'Times New Roman', #'serif',
        'color' : 'darkblue',
        'weight' : 'normal',
        'size' : 14,}


######################
# Read in Data
import csv

try:
    f_name = './algorithm.dat'
    with open(f_name, mode='r') as f:
        print('found '+f_name)
except:
    f_name = '../algorithm.dat'    
print('[Python] Read in data...')
ll = [  [],[],[],[],[],[],[],[],[],[],
        [],[],[],[],[],[],[],[],[],[],
        [],[],[],[],[],[],[],[],[],[],
        [],[],[],[],[],[],[],[],[],[],
        [],[],[],[],[],[],[],[],[],[],
        [],[],[],[],[],[],[],[],[],[],
        [],[],[],[],[],[],[],[],[],[],
        [],[],[],[],[],[],[],[],[],[],
        [],[],[],[],[],[],[],[],[],[]]
with open(f_name, mode='r') as f:
    buf = f.readlines()
    reader = csv.reader(buf)
    for row in reader:
        # print (row)
        try:
            for ind, el in enumerate(row):
                ll[ind].append(float(el))
        except:
            break


for ind, l in enumerate(ll):
    # print len(l)
    if l==[]:
        ll_len = ind
        break
print('\tQuantities amount:', ll_len)
print('\tData Length:', ll[0].__len__())


######################
# Plotting

def get_axis(cNr):
    # fig, axes = plt.subplots(ncols=cNr[0], nrows=cNr[1], dpi=150, sharex=True);
    fig, axes = plt.subplots(ncols=cNr[0], nrows=cNr[1], sharex=True, figsize=(16*0.8, 9*0.8), dpi=80, facecolor='w', edgecolor='k');
    fig.subplots_adjust(right=0.95, bottom=0.1, top=0.95, hspace=0.2, wspace=0.02)    
    # fig.subplots_adjust(right=0.85, bottom=0.1, top=0.95, hspace=0.25)
    if sum(cNr)<=2:
        return axes
    else:
        return axes.ravel()

def plot_it(ax, ylabel, d):
    count = 0
    for k, v in d.items():
        if count == 0:
            count += 1
            # ax.plot(time, v, '--', lw=2, label=k)
            ax.plot(time, v, '-', lw=1)
        else:
            # ax.plot(time, v, '-', lw=2, label=k)
            ax.plot(time, v, '-', lw=1)

    # ax.legend(loc='lower right', shadow=True)
    # ax.legend(bbox_to_anchor=(1.08,0.5), borderaxespad=0., loc='center', shadow=True)
    ax.set_ylabel(ylabel, fontdict=font)
    # ax.set_xlim(0,35) # shared x
    # ax.set_ylim(0.85,1.45)



#################################
# Automatic Code Generation
time = arange(1,ll[0].__len__()+1,1) * 10/4000.000000
# title: Observer
ax_list = get_axis((1,3))
plot_it(ax_list[0], r'$i_s$ [V]', O([
                                             (r'0',   ll[0]),  
                                             (r'1',   ll[1]),  
                                             ]))
plot_it(ax_list[1], r'$\psi_\mu$ [A]', O([
                                             (r'0',   ll[2]),  
                                             (r'1',   ll[3]),  
                                             ]))
plot_it(ax_list[2], r'speed [rpm]', O([
                                             (r'0',   ll[4]),  
                                             ]))
# Automatic END
show()
# savefig(r'C:/Dr.H/(0) GET WORKING/05 3ph 3paramsId/3ph_TDDA_TEX/pic/'+where[:-1], dpi=300)





