import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from AnalyseAdsorption import AnalyseAdsorption

def firstorder(t, q0, qe, k1):
    return qe - qe * np.exp(-k1 * t) + q0 * np.exp(-k1 * t)

def secondorder(t, q0, qe, k2):
    return (k2 * t * qe ^ 2 - k2 * q0 * t * qe + q0) / (k2 * qe * t - k2 * q0 * t + 1)

def elovich(t, q0, alpha, a):
    return -np.log(np.exp(-alpha * q0) - a * alpha * t) / alpha

def weber_morris(t, q0, alpha, a):
    pass

def langmuir(p,vl,pl):
    return (vl*p)/(pl+p)

input_values_L100_40 = {
    'adsorption_file': 'Data/SAMPLE_L100/CH4/40/experimental_raw.csv',
    'leak_calibration_file': 'Data/SAMPLE_L100/CH4/40/experimental_leak.csv',
    'sample_weight': 0.0232661,
    'eos': ['vanderWaals'],
    'specific_volume_calibration': {
        'vv40': 'Data/SAMPLE_L100/vv/300.csv'
    },
    'specific_volume': {
        #'net10': 0.0010,
        #'net7': 0.0007
    },
    'kinetic_models': {
        'firstorder': [firstorder, ['q01', 'qe1', 'k1']]
    },
    'isotherm_models': {
        'langmuir': [langmuir, 'qe1',['vl','pl']]
    },
    'results_location': 'Data/SAMPLE_L100/CH4/40/',
    'unit_system': {
        'pressure': 'bar_g',
        'temperature': 'degC'
    }
}
Lothian40 = AnalyseAdsorption(input_values_L100_40)