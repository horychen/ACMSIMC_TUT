import ACMGUI
def restart():
    print('Restart?')

    # import ACMGUI
    # gui = ACMGUI.ACMGUI()
    # gui.app

    global ani, _index
    # root = Tkinter.Tk()
    # root.withdraw()
    # result = tkMessageBox.askyesno("Restart", "Do you want to restart animation?")
    gui = ACMGUI.ACMGUI()

    # print(int(gui.mb.buttonReply))
    if gui.mb.buttonReply == ACMGUI.QMessageBox.Yes:
        print('Yes clicked.')
        # ani.frame_seq = ani.new_frame_seq() 
        # ani.event_source.start()
        _index = 0
    if gui.mb.buttonReply == ACMGUI.QMessageBox.No:
        print('No clicked.')
        plt.close()
    if gui.mb.buttonReply == ACMGUI.QMessageBox.Cancel:
        print('Cancelled')
        raise KeyboardInterrupt
    # gui.app.exec_()

from pylab import subplots, mpl, plt, np
import pandas as pd
import time as time_package

def get_data(fname):
    try:
        data = pd.read_csv(fname) #, skiprows=2, nrows=7) # na_values = ['no info', '.', '']
    except pd.errors.EmptyDataError as error:
        # print(str(error))
        # time_package.sleep(1)
        raise error
    keys = data.keys()
    return data, keys
fname = './algorithm.dat'
data, keys = get_data(fname)
# print(keys)
# print(data)
# quit()

fig, ax_list = plt.subplots(len(keys), sharex=True, figsize=(16*0.8, 9*0.8), dpi=80, facecolor='w', edgecolor='k', constrained_layout=False)

SAMP_TS = 10/4000.000000 # TS=1/8000, down_exe=2, down-sampling 10: SAMP_TS=1/400
n = number_points_draw_at_once = 100 # plot n points at a time
POST_TS = SAMP_TS * n

_index = 0
def animate(_):
    global data, keys, _index
    if not bool_animate_pause:
        if _index > 0:
            time = np.arange(1, n*_index+1) * SAMP_TS
            for idx, ax in enumerate(ax_list):
                ax.cla()
                try:
                    ax.plot(time, data[keys[idx]][:n*_index])

                except ValueError as error:
                    # print(str(error))
                    data_backup = data
                    data, keys = get_data(fname)
                    print('Read data.')
                    if len(data[keys[0]]) == len(data_backup[keys[0]]):
                        print('Plotting data are exausted. Quit now.')
                        _index = -1
                        break
                else:
                    ax.set_ylabel(keys[idx])
            ax_list[-1].set_xlabel('Time [s]')
            # ax.set_title("Frame {}".format(i))
            # ax.relim()
            # ax.autoscale_view()

        if _index == -1:
            restart()
            # raise KeyboardInterrupt
        _index += 1

def onClick(event):
    global bool_animate_pause
    bool_animate_pause ^= True

import matplotlib.animation as animation
bool_animate_pause = False
ani = animation.FuncAnimation(fig, animate, blit=False, interval=10, repeat=False)
fig.canvas.mpl_connect('button_press_event', onClick)
time_package.sleep(0.1)

print('Animate starts...')
plt.show()

