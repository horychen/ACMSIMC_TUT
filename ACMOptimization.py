from ACMConfig import load_json, write_to_ACMConfig_header_file
import json, subprocess


def subprocess_cmd(command_string):
    process = subprocess.Popen(command_string,stdout=subprocess.PIPE, shell=True)
    streamdata = proc_stdout = process.communicate()[0].strip()
    # print(proc_stdout)
    return process

    # example:
    # subprocess_cmd('echo a; echo b')

if __name__ == '__main__':

    # select tempalte json config
    select_config = '#02 ACMSIMC IM FAA86'

    # load template json config
    raw_config_dicts = load_json('ACMConfig_backup.json', bool_print_out=False)

    # change json config for simulation variants
    NUMBER_OF_RUNS = 5
    for i in range(0, NUMBER_OF_RUNS):
        raw_config_dicts[select_config+f' - Variant Speed P {i:03d}'] = dict(raw_config_dicts[select_config]) # copy dict
        raw_config_dicts[select_config+f' - Variant Speed P {i:03d}']['SPEED_LOOP_PID_PROPORTIONAL_GAIN'] = 0.05 + 0.1*i
        raw_config_dicts[select_config+f' - Variant Speed P {i:03d}']['DATA_FILE_NAME'] = f"\"im_faa86_VSP{i:03d}.dat\""

    # dump json config's to file
    with open('ACMConfig.json', 'w') as f:
        json.dump(raw_config_dicts, f, indent=4)

    # with all those json config's, for each json config, re-write ACMConfig.h and compile main.c and run it.
    if True:
        for i in range(0, NUMBER_OF_RUNS):
            write_to_ACMConfig_header_file(raw_config_dicts[select_config+f' - Variant Speed P {i:03d}'])
            child = subprocess_cmd("gcc main.c controller.c observer.c -L. -o main && start cmd /c main")
            rc = child.returncode
            if rc == 0:
                continue
            else:
                raise Exception('The C codes end with returncode %d'%(rc))

    # now with a pile of .dat files, we are ready to plot
    import ACMPlot
    import pandas as pd
    import numpy as np
    from collections import OrderedDict as O
    df_info = pd.read_csv(r"./info.dat", na_values = ['1.#QNAN', '-1#INF00', '-1#IND00'])
    bool_initialized = False
    for i in range(0, NUMBER_OF_RUNS):

        df_profiles = pd.read_csv(f"im_faa86_VSP{i:03d}.dat", na_values = ['1.#QNAN', '-1#INF00', '-1#IND00'])

        if bool_initialized == False:
            bool_initialized = True
            no_samples = df_profiles.shape[0]
            no_traces  = df_profiles.shape[1]
            time = np.arange(1, no_samples+1) * df_info['DOWN_SAMPLE'].values[0] * df_info['TS'].values[0]
            ax_list = []
            for i in range(0, no_traces, 6):
                ax_list += list(ACMPlot.get_axis((1,6)))

        for idx, key in enumerate(df_profiles.keys()):
            ACMPlot.plot_it(ax_list[idx], key, O([
                                            (str(idx), df_profiles[key]),
                                            ]), time)
    from pylab import plt
    plt.show()
    quit()
