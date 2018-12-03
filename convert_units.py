def _get_pressure_index_():
    p_index = {
        'Pa':0,
        'bar':1,
        'atm':2,
        'psi':3,
        'torr':4,
        'g':5,
        'a':6
    }
    return p_index

def _get_pressure_matrix_():
    p_matrix = [
                [1, 1.00000000000000e-05, 9.86923266716000e-06, \
                 0.000145037891149100, 0.00750061682704200, 101325, 0],
                [100000, 1, 0.986923266716000, 14.5037891149100, \
                 750.061682704200, 1.01325000000000, 0],
                [101325, 1.01325000000000, 1, 14.6959643206800, 760, 1, 0],
                [6894.75000000000, 0.0689475000000000, 0.0680458919319000, \
                 1, 51.7148778682500, 14.6959500000000, 0],
                [133.322368421100, 0.00133322368421100, 0.00131578947368400, \
                 0.0193367951587900, 1, 760.002100000000, 0]
                ]
    return p_matrix

def convert_pressure(input_value, input_units, output_units):

    matrix = _get_pressure_matrix_()
    index = _get_pressure_index_()


    si = input_units.split("_");
    so = output_units.split("_");

    if not si[0] in index.keys() and so[0] in index.keys():
        print('invalid pressure units')

    absolute_input_value = input_value + matrix[index[si[0]]][index[si[1]]];
    converted_input_value = absolute_input_value * matrix[index[si[0]]][index[so[0]]];
    output_value = converted_input_value - matrix[index[so[0]]][index[so[1]]];
    return output_value

def _get_temperature_index_():
    t_index = {
        'degC':0,
        'degF':1,
        'K':2
    }
    return t_index

def convert_temperature(input_value, input_units, output_units):
    if input_units == 'degC':
        if output_units == 'degC':
            output_value =  input_value * 1
        elif output_units == 'degF':
            output_value = (input_value * 9 / 5) + 32
        elif output_units == 'K':
            output_value = input_value + 273.15
    elif input_units == 'degF':
        if output_units == 'degC':
            output_value =  (input_value - 32) * 5 / 9
        elif output_units == 'degF':
            output_value = input_value * 1
        elif output_units == 'K':
            output_value = (input_value + 459.67) * 5 / 9,
    elif input_units == 'K':
        if output_units == 'degC':
            output_value =  input_value - 273.15
        elif output_units == 'degF':
            output_value = (input_value * 9 / 5) - 459.67
        elif output_units == 'K':
            output_value = input_value * 1
    return output_value

def _get_volume_index_():
    v_index = {
        'm3':0,
        'l':1,
        'ml':2,
        'ft3':3
    }
    return v_index

def _get_volume_matrix_():
    v_matrix = [
                [1, 1000, 1000000, 35.3146667000000],
                [0.00100000000000000, 1, 1000, 0.0353146667000000],
                [1.00000000000000e-06, 0.00100000000000000, 1, 3.53146667000000e-05],
                [0.0283168466000000, 28.3168466000000, 28316.8466000000, 1]
               ]
    return v_matrix

def convert_volume(input_value, input_units, output_units):

    index = _get_volume_index_()
    matrix = _get_volume_matrix_()

    if not input_units in index.keys() and output_units in index.keys():
        print('invalid volume units')

    output_value = input_value * matrix[index[input_units]][index[output_units]]
    return output_value

def _get_weight_index_():
    w_index = {
        'kg': 0,
        'g': 1,
        'mg': 2,
        'lb': 3,
        'ton': 4,
        'USton': 5
    }
    return w_index

def _get_weight_matrix_():
    w_matrix = [
        [1,1000,1000000,2.20462,0.001,0.00110231],
        [0.001,1,1000,0.00220462,1.00e-06,1.10230e-06],
        [1.00e-06,0.001,1,2.20460e-6,1.00e-9,1.10230e-09],
        [0.453592,453.592,453592,1,0.000453592,0.0005],
        [1000,1.0e6,1.0e9,2204.62,1,1.10231],
        [907.185,907185,9.07185e8,2000,0.907185,1]
    ]
    return w_matrix

def convert_weight(input_value, input_units, output_units):

    index = _get_weight_index_()
    matrix = _get_weight_matrix_()

    if not input_units in index.keys() and output_units in index.keys():
        print('invalid units')

    output_value = input_value * matrix[index[input_units]][index[output_units]]
    return output_value

def _get_mole_index_():
    m_index = {
        'mol':0,
        'mmol':1,
        'SCCM':2,
        'SCF':3
    }
    return m_index

def convert_moles(input_value, input_units, output_units):
    if input_units == 'mol':
        if output_units == 'mol':
            output_value = input_value * 1
        elif output_units == 'mmol':
            output_value = input_value * 1000
        elif output_units == 'SCCM':
            output_value = input_value * 24710.480
        elif output_units == 'SCF':
            output_value = input_value * 1.883943106091513; # 0.986923266716*70*24.710480*0.03531467/32
    elif input_units == 'mmol':
        if output_units == 'mol':
            output_value = input_value * 1000;
        elif output_units == 'mmol':
            output_value = input_value * 1;
        elif output_units == 'SCCM':
            output_value = input_value * 24.710480;
        elif output_units == 'SCF':
            output_value = input_value * 0.001883943106091513;
    elif input_units == 'SCCM':
        if output_units == 'mol':
            output_value = input_value / 24710.480;
        elif output_units == 'mmol':
            output_value = input_value / 24.710480;
        elif output_units == 'SCCM':
            output_value = input_value * 1;
        elif output_units == 'SCF':
            output_value = input_value * 0.076240651986182; #0.986923266716*70*0.03531467/32
    elif input_units == 'SCF':
        if output_units == 'mol':
            output_value = input_value / 1.883943106091513; # 0.986923266716*70*24.710480*0.03531467/32
        elif output_units == 'mmol':
            output_value = input_value * 530.8015920261154;
        elif output_units == 'SCCM':
            output_value = input_value * 13.116362123729489;
        elif output_units == 'SCF':
            output_value = input_value * 1;
    return output_value

def _get_time_index_():
    t_index = {
        'sec': 0,
        'min' : 1,
        'hr' : 2,
        'day' :3
    }
    return t_index

def _get_time_matrix_():
    t_matrix = [
        [1,0.0166666666666667,2.777777777777778e-4,1.157407407407407e-5],
        [60,1,0.0166666666666667,6.944444444444444e-4],
        [3600,60,1,0.0416666666666667],
        [86400,1440,24,1]
    ]
    return t_matrix

def convert_time(input_value, input_units, output_units):

    index = _get_time_index_()
    matrix = _get_time_matrix_()

    if not input_units in index.keys() and output_units in index.keys():
        print('invalid time units')

    output_value = input_value * matrix[index[input_units]][index[output_units]]
    return output_value

