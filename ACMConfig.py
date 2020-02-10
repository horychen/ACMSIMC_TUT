import json
def load_json(fname='ACMConfig.json', bool_print_out=True):
    # In case there are comments in your json file, option 1 is to try https://pypi.org/project/jsmin/, option 2 is to remove it with GetJsonFromFileWithComments
    # def GetJsonFromFileWithComments(filePath):
    #     # https://stackoverflow.com/questions/29959191/how-to-parse-json-file-with-c-style-comments
    #     contents = ""
    #     fh = open(filePath)
    #     for line in fh:
    #         cleanedLine = line.split("//", 1)[0]
    #         if len(cleanedLine) > 0 and line.endswith("\n") and "\n" not in cleanedLine:
    #             cleanedLine += "\n"
    #         contents += cleanedLine
    #     fh.close
    #     while "/*" in contents:
    #         preComment, postComment = contents.split("/*", 1)
    #         contents = preComment + postComment.split("*/", 1)[1]
    #     return contents

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

def write_to_ACMConfig_header_file(config_pairs):
    str_c_defines = ''
    for key, val in config_pairs.items():
        str_c_defines += f'#define {key} {val}\n'

    bare_template = f'''
#ifndef ACMCONFIG_H
#define ACMCONFIG_H


#define VVVF_CONTROL 0
#define IFOC 1
#define DFOC 2

#define TAJIMA96 0 
#define TOTALLY_ADAPTIVE_OBSERVER 1 // Flux-MRAS
#define FIRST_APPROXIMATION_ANDERSON86 2

{str_c_defines}

#define MACHINE_TS         (TS*TS_UPSAMPLING_FREQ_EXE) //1.25e-4 
#define MACHINE_TS_INVERSE (TS_INVERSE*TS_UPSAMPLING_FREQ_EXE_INVERSE) // 8000

#endif
'''

    with open('ACMConfig.h', 'w') as f:
        f.write(bare_template)


if __name__ == '__main__':

    raw_config_dicts = load_json(fname='ACMConfig_backup.json', bool_print_out=False)
    config_pairs = raw_config_dicts['#02 ACMSIMC IM FAA86']
    write_to_ACMConfig_header_file(config_pairs)























readme = '''
"#02 ACMSIMC IM FAA86"

"NUMBER_OF_STEPS": 150000,                          <--- How many steps to simulate. NUMBER_OF_STEPS * MACHINE_TS = simulation time
"DOWN_SAMPLE": 1,                                   <--- For plot. Sometimes you may want to reduce the waveform sampling step to make the codes run faster.

"CONTROL_STRATEGY": "NULL_D_AXIS_CURRENT_CONTROL",  <--- Control method.
"OBSERVER_APPLIED": "TAJIMA96",                     <--- Observer applied.

"SENSORLESS_CONTROL": false,                        <--- Use reconstructed speed and position for feedback control.
"VOLTAGE_CURRENT_DECOUPLING_CIRCUIT": false,        <--- Decouple between d-axis and q-axis current loops.
"SATURATED_MAGNETIC_CIRCUIT": false,                <--- Saturation.
"INVERTER_NONLINEARITY": false,                     <--- Inverter nonlinearity.

"TS": 2.5e-4,                                       <--- Sampling time in DSP or controller. Make sure this TS matches your actual control codes execution frequency, e.g., PWM having a carrier frequency of 4 kHz (Additionally, if it interupts at both 波峰 and 波谷, you can have a interupt freuquency of 8 kHz).
"TS_INVERSE": 4000,                                 <--- TS * TS_INVERSE = 1
"TS_UPSAMPLING_FREQ_EXE": 0.5,                      <--- Up-sampling TS by a factor of this number to obtain MACHINE_TS. That is, the machine dynamics are executed in a higher frequency than your controller, to imitate a continuous time system
"TS_UPSAMPLING_FREQ_EXE_INVERSE": 2,                <--- 1/TS_UPSAMPLING_FREQ_EXE

"IM_NUMBER_OF_POLE_PAIRS": 2,
"IM_STATOR_RESISTANCE": 3.04,
"IM_INVERSE_GAMMA_CIRCUIT_MAGNETIZING_INDUCTANCE": 0.4482,
"IM_INVERSE_GAMMA_CIRCUIT_ROTOR_RESISTANCE": 1.69,
"IM_T_CIRCUIT_STATOR_LEAKAGE_INDUCTANCE": 0.0126,
"IM_T_CIRCUIT_ROTOR_LEAKAGE_INDUCTANCE": 0.0126,
"IM_SHAFT_INERTIA": 0.0636,

"SPEED_LOOP_PID_PROPORTIONAL_GAIN": 0.5,
"SPEED_LOOP_PID_INTEGRAL_TIME_CONSTANT": 1.05,
"SPEED_LOOP_PID_DIREVATIVE_TIME_CONSTANT": 0,
"SPEED_LOOP_LIMIT_NEWTON_METER": 8,

"SPEED_LOOP_CEILING": 40,                           <--- velocity control loop execution frequency is 40 times slower than current control loop execution frequency

"CURRENT_LOOP_PID_PROPORTIONAL_GAIN": 15,
"CURRENT_LOOP_PID_INTEGRAL_TIME_CONSTANT": 0.08,
"CURRENT_LOOP_PID_DIREVATIVE_TIME_CONSTANT": 0,
"CURRENT_LOOP_LIMIT_VOLTS": 650,

"DATA_FILE_NAME": "\"im_taa86.dat\"",

"PC_SIMULATION": true
'''
# print(readme)

