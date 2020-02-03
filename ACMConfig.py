
def load_json(fname='ACMConfig.json', bool_print_out=True):
    import json
    with open(fname, 'r') as f:
        raw_config_dicts = json.load(f)

    def decode_raw_fea_configs(raw_config_dicts):
        for key, val in raw_config_dicts.items():
            print('\n', key)
            for ke, va in val.items():
                print('\t', ke+':', va)
    if bool_print_out:
        decode_raw_fea_configs(raw_config_dicts)
    return raw_config_dicts

raw_config_dicts = load_json(bool_print_out=False); config_pairs = raw_config_dicts['#01 ACMSIMC PMSM EEMF']

str_c_defines = ''
for key, val in config_pairs.items():
    str_c_defines += f'#define {key} {val}\n'

bare_template = f'''
#ifndef ACMCONFIG_H
#define ACMCONFIG_H

#define NULL_D_AXIS_CURRENT_CONTROL -1
#define MTPA -2 // not supported

{str_c_defines}
#define TS (MACHINE_TS*DOWN_FREQ_EXE) //2.5e-4 
#define TS_INVERSE (MACHINE_TS_INVERSE*DOWN_FREQ_EXE_INVERSE) // 4000

#endif
'''

with open('ACMConfig.h', 'w') as f:
    f.write(bare_template)

