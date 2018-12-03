def calibrate_specific_volume(calibration_file, reference_cell_volume, sample_cell_volume):
    import pandas as pd
    import CoolProp.CoolProp as CP
    from convert_units import convert_pressure
    from convert_units import convert_temperature
    from convert_units import convert_weight

    specific_volume_dataframe = pd.read_excel(calibration_file)

    # convert pressure to SI units
    pressure_units = str(specific_volume_dataframe['pressure_units'][0])
    if pressure_units == 'Pa_a':
        specific_volume_dataframe['pressure_SI'] = specific_volume_dataframe['pressure']
    else:
        specific_volume_dataframe['pressure_SI'] = specific_volume_dataframe['pressure'] \
            .apply(convert_pressure, args=(pressure_units, 'Pa_a'))

    # convert temperature to SI units
    temperature_units = str(specific_volume_dataframe['temperature_units'][0])
    if temperature_units == 'K':
        specific_volume_dataframe['temperature_SI'] = specific_volume_dataframe['temperature']
    else:
        specific_volume_dataframe['temperature_SI'] = specific_volume_dataframe['temperature'] \
            .apply(convert_temperature, args=(temperature_units, 'K'))

    # convert weight to SI units
    weight_units = str(specific_volume_dataframe['weight_units'][0])
    if temperature_units == 'kg':
        specific_volume_dataframe['weight_SI'] = specific_volume_dataframe['weight']
    else:
        specific_volume_dataframe['weight_SI'] = specific_volume_dataframe['weight'] \
            .apply(convert_weight, args=(weight_units, 'kg'))

    # get mean of each pressure tag
    pressure_series = specific_volume_dataframe[['pressure_tag', 'pressure_SI']] \
        .groupby('pressure_tag').mean()['pressure_SI']

    # convert to array
    pressure_array = [pressure_series.loc['initial'], pressure_series.loc['final'],
                      pressure_series.loc['equilibrium']]

    # get mean for each temperature tag
    temperature_series = specific_volume_dataframe[['pressure_tag', 'temperature_SI']] \
        .groupby('pressure_tag').mean()['temperature_SI']

    # convert to array
    temperature_array = [temperature_series.loc['initial'], temperature_series.loc['final'],
                         temperature_series.loc['equilibrium']]

    # real gas law based on compressibility factors from coolprop is used as calibration is done with
    # helium at low pressures

    # get probe molecule from dataframe
    probe_molecule = str(specific_volume_dataframe['adsorptive_id'][0])

    # get sample weight from dataframe
    sample_weight = specific_volume_dataframe['weight'][0]

    # initialize phi = pressure/compressibility factor as phi as an empty array
    phi = []

    # loop through each pressure tag to get phi
    for i in range(3):
        z = CP.PropsSI('Z', 'P', pressure_array[i], 'T', temperature_array[i], probe_molecule)
        phi.append(pressure_array[i] / z)

    # calculate void volume = reference volume * (phi(equilibrium) - phi(final)) / (phi(final) - phi(initial))
    sample_cell_void_volume = reference_cell_volume * (phi[1] - phi[2]) / (phi[2] - phi[0])

    # calculate sample volume based on calculated void volume
    sample_volume = sample_cell_volume - sample_cell_void_volume

    # calculate specific volume based on weight
    sp_volume = sample_volume / sample_weight

    return sp_volume

def calibrate_leak_rate(calibration_file, reference_cell_volume, sample_cell_volume, sample_specific_volume, eos):
    from convert_units import convert_pressure
    from convert_units import convert_temperature
    from convert_units import convert_weight
    from calculate_moles import calculate_moles_row
    import pandas as pd
    from scipy.stats import linregress

    leak_data = pd.read_excel(calibration_file)

    # convert pressure to SI units
    pressure_units = str(leak_data['pressure_units'][0])
    if pressure_units == 'Pa_a':
        leak_data['pressure_SI'] = leak_data['pressure']
    else:
        leak_data['pressure_SI'] = leak_data['pressure'] \
            .apply(convert_pressure, args=(pressure_units, 'Pa_a'))

    # convert temperature to SI units
    temperature_units = str(leak_data['temperature_units'][0])
    if temperature_units == 'K':
        leak_data['temperature_SI'] = leak_data['temperature']
    else:
        leak_data['temperature_SI'] = leak_data['temperature'] \
            .apply(convert_temperature, args=(temperature_units, 'K'))

    # convert weight to SI units
    weight_units = str(leak_data['weight_units'][0])
    if temperature_units == 'kg':
        leak_data['weight_SI'] = leak_data['weight']
    else:
        leak_data['weight_SI'] = leak_data['weight'] \
            .apply(convert_weight, args=(weight_units, 'kg'))

    # get seconds since start
    leak_data['seconds_since_start'] = [(x - leak_data['time'][0]).total_seconds() for x in leak_data['time']]

    sample_weight = leak_data['weight'].mean()

    # insert columns for values required for leak calculation
    leak_data.insert(5, 'sample_specific_volume', sample_specific_volume)
    leak_data.insert(3, 'eos', eos)
    leak_data.insert(6, 'reference_cell_volume', reference_cell_volume)
    leak_data.insert(7, 'sample_cell_volume', sample_cell_volume)
    leak_data.insert(8, 'sample_weight', sample_weight)
    leak_data['void_volume'] = leak_data['reference_cell_volume'] \
                               + leak_data['sample_cell_volume'] - (leak_data['sample_specific_volume']
                                                                    * leak_data['sample_weight'])

    # calculate current moles
    leak_data['moles_present'] = leak_data.apply(lambda row: calculate_moles_row(row), axis=1)

    # fit regression between seconds and moles present
    x = leak_data['seconds_since_start']
    y = leak_data['moles_present']
    regression_result = linregress(x, y)

    # moles leaked is given by the slope
    linear_leak_rate = regression_result.slope

    # leak pressure is taken as the mean pressure for the whole experiment
    leak_pressure = leak_data['pressure_SI'].mean()

    return linear_leak_rate, leak_pressure
