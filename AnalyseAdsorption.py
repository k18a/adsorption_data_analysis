class PlotIsotherm:
    def __init__(self, isotherm_parameters, verbose=False):
        self.isotherm_parameters = isotherm_parameters
        self.verbose(verbose)

    def verbose(self, verbose=True):
        if verbose:
            def verboseprint(*args):
                # print each argument separately so caller doesn't need to
                # stuff everything to be printed into a single string
                for arg in args:
                    print(arg)
        else:
            def verboseprint(*args):
                # do nothing function
                pass
        self._verbose_print_ = verboseprint

    def plot_isotherm_fits(self, axis, label, model_name, isotherm_tags = ['adsorption','desorption'],
                           marker_formats = 'NaN', line_formats = 'NaN', colors = 'NaN',pressure_range=(0,6000000),
                           pressure_units = 'Pa_a', mole_units = 'mol', weight_units = 'kg'):
        from itertools import cycle
        import numpy as np
        from isotherm_models import get_isotherm_mapping_for_plots
        from process_data import define_plot_colors, define_plot_lines, define_plot_markers
        from convert_units import convert_pressure, convert_moles, convert_weight

        # define isotherm functions and relevant column names
        isotherm_function_mapping = get_isotherm_mapping_for_plots()

        # set plot formats
        if marker_formats == 'NaN':
            marker_formats = define_plot_markers()
        elif type(marker_formats) == list:
            marker_formats = cycle(marker_formats)
        else:
            print('marker formats not recognized; please input list')
        if line_formats =='NaN':
            line_formats = define_plot_lines()
        elif type(line_formats) == list:
            line_formats = cycle(line_formats)
        else:
            print('line formats not recognized; please input list')
        if colors == 'NaN':
            colors = define_plot_colors()
        elif type(colors) == list:
            colors = cycle(colors)
        else:
            print('colors not recognized; please input list')

        # iterate through adsorption and desorption isotherm tags
        for isotherm_tag in isotherm_tags:
            # check if model name is valid
            if not model_name in  isotherm_function_mapping.keys():
                print('model name not recognized')
            # get fitting function
            fitting_function = isotherm_function_mapping[model_name][0]
            self._verbose_print_('fitting_function', fitting_function)
            parameters = isotherm_function_mapping[model_name][1]
            self._verbose_print_('parameter names', parameters)
            pressure = np.linspace(pressure_range[0], pressure_range[1], 500)
            fitting_parameters = self.isotherm_parameters[self.isotherm_parameters['isotherm_tag']
                                                          == isotherm_tag][parameters].values[0]
            self._verbose_print_('parameters', fitting_parameters)
            adsorbed = fitting_function(pressure, *fitting_parameters)
            x = [convert_pressure(p,'Pa_a',pressure_units) for p in pressure]
            y = [q*convert_moles(1,'mol',mole_units)/convert_weight(1,'kg',weight_units)
                 for q in adsorbed]
            axis.plot(x, y, linestyle=line_formats.next(), marker = marker_formats.next(),
                      color=colors.next(), label=label)
        axis.set_xlabel('equilibrium pressure ({})'.format(pressure_units))
        axis.set_ylabel('equilibrium adsorbed ({}/{})'.format(mole_units,weight_units))

class AnalyseIsotherm(PlotIsotherm):
    def __init__(self, isotherm_data, verbose=False):
        from numpy import NaN
        import pandas as pd
        self.isotherm_data = isotherm_data
        self.isotherm_parameters = pd.DataFrame()
        self.pseudo_saturation_pressure = NaN
        self.micropore_volume = NaN
        self.adsorbed_phase_density = NaN
        self.verbose(verbose)
        self._initialize_data_()

    def verbose(self, verbose=True):
        if verbose:
            def verboseprint(*args):
                # print each argument separately so caller doesn't need to
                # stuff everything to be printed into a single string
                for arg in args:
                    print(arg)
        else:
            def verboseprint(*args):
                # do nothing function
                pass
        self._verbose_print_ = verboseprint

    def _initialize_data_(self):
        """
        checks if input data files are valid
        accepted data types include int, float, data frame, and file paths
        :return:
        """
        from pandas.core.frame import DataFrame
        assert(type(self.isotherm_data)==DataFrame)

        # read adsorbate id, adsorptive id, isotherm id, specific volume id, and eos from data frame
        self.adsorbent_id = self.isotherm_data['adsorbent_id'][0]
        self.adsorptive_id = self.isotherm_data['adsorptive_id'][0]
        self.isotherm_id = self.isotherm_data['isotherm_id'][0]
        self.specific_volume_id = self.isotherm_data['specific_volume_id'][0]
        self.eos = self.isotherm_data['eos'][0]
        # print assignments if verbose is true
        self._verbose_print_('adsorbent_id', self.adsorbent_id, 'adsorptive_id', self.adsorptive_id,
                             'specific_volume_id', self.specific_volume_id, 'eos', self.eos)

        self.pressure_steps = list(self.isotherm_data['pressure_step'])

    def fit_isotherm_data(self):
        from scipy.optimize import curve_fit
        import numpy as np
        from isotherm_models import linear, langmuir, freundlich, dubinin_radushkevich
        from thermodynamic_calculations import calculate_pseudo_saturation_pressure
        import warnings
        warnings.filterwarnings("ignore")

        isotherm_function_mapping = {
            'linear': [linear, ['henry_constant']],
            'langmuir': [langmuir, ['langmuir_volume','langmuir_pressure']],
            'freundlich': [freundlich, ['freundlich_k','freundlich_n']],
            'dubinin_radushkevich': [lambda P,E:
                                     dubinin_radushkevich(P,T,P0,V0,E),
                                     ['dr_E']]
        }

        # set pseudo saturation pressure
        self.pseudo_saturation_pressure = calculate_pseudo_saturation_pressure(
            self.pseudo_saturation_pressure,
            self.adsorptive_id,
            self.isotherm_data['equilibrium_temperature'].mean()
        )
        # set micro pore volume
        if self.micropore_volume == np.NaN:
            self.micropore_volume = float(input('enter micro-pore volume as int or float:'))
        # set adsorbed phase density
        if self.adsorbed_phase_density == np.NaN:
            self.adsorbed_phase_density = float(input('enter adsorbed phase density as int or float:'))

        print('Isotherm parameters:')
        templateh = '{:^15} {:^25} {:^50}'
        print(templateh.format('isotherm_tag', 'fitting function', 'fitting parameters'))
        print(templateh.replace(':', ':-').format('', '', '', ''))
        templateb = '{:^15} {:^25} {:^50}'

        # set temporary langmuir volume to a high value so you'll know if it changes
        temp_langmuir_volume = 9999999999

        # define isotherm constants
        T = self.isotherm_data['equilibrium_temperature'].mean()
        P0 = self.pseudo_saturation_pressure
        V0 = self.micropore_volume*self.adsorbed_phase_density

        isotherm_tags = set(self.isotherm_data['isotherm_tag'])

        for isotherm_tag in isotherm_tags:
            isotherm = self.isotherm_data[self.isotherm_data['isotherm_tag'] == isotherm_tag]

            # fix modelled and constant variables
            x_var = isotherm['equilibrium_pressure']
            y_var = isotherm['equilibrium_adsorbed']
            # iterate through user input models
            for model_name, model_values in isotherm_function_mapping.iteritems():
                try:
                    # define fitting function and fitting parameters
                    fitting_function = model_values[0]
                    fitting_parameters = model_values[1]

                    # insert columns into data frame for all model parameters
                    for parameter in fitting_parameters:
                        if not parameter in list(self.isotherm_parameters):
                            self.isotherm_parameters.insert(len(self.isotherm_parameters.columns), parameter, np.NaN)
                    self._verbose_print_('isotherm parameters before fitting', self.isotherm_parameters)

                    self._verbose_print_('x_var', x_var, 'y_var', y_var)

                    if model_name == 'langmuir':
                        parameter_bounds = ((0,0),(np.inf,np.inf))
                        if temp_langmuir_volume != 9999999999:
                            parameter_bounds = ((temp_langmuir_volume*0.99999999,0),
                                                 (temp_langmuir_volume*1.00000001,np.inf))
                        popt = curve_fit(fitting_function, x_var, y_var, bounds=parameter_bounds)
                        temp_langmuir_volume = popt[0][0]
                    elif model_name == 'dubinin_radushkevich':
                        # perform curve fitting for DR
                        parameter_bounds = ((0),(np.inf))
                        popt = curve_fit(fitting_function, x_var, y_var, p0=(20000), bounds=parameter_bounds)
                    else:
                        # perform curve fitting to get fitting variables for all other
                        popt = curve_fit(fitting_function, x_var, y_var)
                except:
                    popt = [[np.NaN]*len(fitting_parameters)]
                finally:
                    self.isotherm_parameters.loc[
                        self.isotherm_parameters['isotherm_tag']==isotherm_tag, fitting_parameters
                    ] = popt[0]
                    print(templateb.format(isotherm_tag,model_name,str(popt[0])))

        self.isotherm_parameters.insert(6,'adsorption_temperature',T)
        self.isotherm_parameters.insert(7,'micro_pore_volume',self.micropore_volume)
        self.isotherm_parameters.insert(8,'adsorbed_phase_density',self.adsorbed_phase_density)
        self.isotherm_parameters.insert(9,'pseudo_saturation_pressure',self.pseudo_saturation_pressure)
        self.isotherm_parameters.insert(10,'V0',V0)

        self._verbose_print_('isotherm parameters after fitting', self.isotherm_parameters)

    def plot_isotherm(self, axis, label, isotherm_tags = ['adsorption','desorption'],
                      marker_formats = 'NaN', line_formats = 'NaN', colors = 'NaN',
                      pressure_units = 'Pa_a', mole_units = 'mol', weight_units = 'kg'):
        from itertools import cycle
        import numpy as np
        from process_data import define_plot_colors, define_plot_lines, define_plot_markers
        from convert_units import convert_pressure, convert_moles, convert_weight

        # set plot formats
        if marker_formats == 'NaN':
            marker_formats = define_plot_markers()
        elif type(marker_formats) == list:
            marker_formats = cycle(marker_formats)
        else:
            print('marker formats not recognized; please input list')
        if line_formats =='NaN':
            line_formats = define_plot_lines()
        elif type(line_formats) == list:
            line_formats = cycle(line_formats)
        else:
            print('line formats not recognized; please input list')
        if colors == 'NaN':
            colors = define_plot_colors()
        elif type(colors) == list:
            colors = cycle(colors)
        else:
            print('colors not recognized; please input list')

        # iterate through adsorption and desorption isotherm tags
        for isotherm_tag in isotherm_tags:
            isotherm = self.isotherm_data[self.isotherm_data['isotherm_tag'] == isotherm_tag]
            pressure = isotherm['equilibrium_pressure']
            adsorbed = isotherm['equilibrium_adsorbed']
            x = [convert_pressure(p,'Pa_a',pressure_units) for p in pressure]
            y = [q*convert_moles(1,'mol',mole_units)/convert_weight(1,'kg',weight_units)
                 for q in adsorbed]
            axis.plot(x, y, linestyle=line_formats.next(), marker = marker_formats.next(),
                      color=colors.next(), label=label)
        axis.set_xlabel('equilibrium pressure ({})'.format(pressure_units))
        axis.set_ylabel('equilibrium adsorbed ({}/{})'.format(mole_units,weight_units))

class PlotAdsorption(AnalyseIsotherm):
    def __init__(self, adsorption_data, isotherm_data, verbose=False):
        self.adsorption_data = adsorption_data
        self.isotherm_data = isotherm_data
        self.verbose(verbose)

    def verbose(self, verbose=True):
        if verbose:
            def verboseprint(*args):
                # print each argument separately so caller doesn't need to
                # stuff everything to be printed into a single string
                for arg in args:
                    print(arg)
        else:
            def verboseprint(*args):
                # do nothing function
                pass
        self._verbose_print_ = verboseprint

    def _initialize_data_(self):
        """
        checks if input data files are valid
        accepted data types include int, float, data frame, and file paths
        :return:
        """
        from pandas.core.frame import DataFrame
        assert(type(self.isotherm_data)==DataFrame)

        # read adsorbate id, adsorptive id, isotherm id, specific volume id, and eos from data frame
        self.adsorbent_id = self.adsorption_data['adsorbent_id'][0]
        self.adsorptive_id = self.adsorption_data['adsorptive_id'][0]
        self.isotherm_id = self.adsorption_data['isotherm_id'][0]
        self.specific_volume_id = self.adsorption_data['specific_volume_id'][0]
        self.eos = self.adsorption_data['eos'][0]
        # print assignments if verbose is true
        self._verbose_print_('adsorbent_id', self.adsorbent_id, 'adsorptive_id', self.adsorptive_id,
                             'specific_volume_id', self.specific_volume_id, 'eos', self.eos)

        self.pressure_steps = list(set(self.adsorption_data['pressure_step']))

    def plot_kinetics_fits(self, axis, label, model_name, pressure_step,
                           marker_formats = 'NaN', line_formats = 'NaN', colors = 'NaN',time_range=(0,500),
                           time_units = 'sec', mole_units = 'mol', weight_units = 'kg'):
        from itertools import cycle
        import numpy as np
        from process_data import define_plot_colors, define_plot_lines, define_plot_markers
        from convert_units import convert_time, convert_moles, convert_weight
        from kinetic_models import get_kinetic_mapping_for_plots

        # create function mapping from name to function and parameters as dictionary key value pairs

        kinetic_model_function = get_kinetic_mapping_for_plots()

        # set plot formats
        if marker_formats == 'NaN':
            marker_formats = define_plot_markers()
        elif type(marker_formats) == list:
            marker_formats = cycle(marker_formats)
        else:
            print('marker formats not recognized; please input list')
        if line_formats =='NaN':
            line_formats = define_plot_lines()
        elif type(line_formats) == list:
            line_formats = cycle(line_formats)
        else:
            print('line formats not recognized; please input list')
        if colors == 'NaN':
            colors = define_plot_colors()
        elif type(colors) == list:
            colors = cycle(colors)
        else:
            print('colors not recognized; please input list')

        fitting_function = kinetic_model_function[model_name][0]
        self._verbose_print_('fitting_function', fitting_function)
        parameters = kinetic_model_function[model_name][1]
        self._verbose_print_('parameter names', parameters)
        time = np.linspace(time_range[0], time_range[1], 500)
        fitting_parameters = self.isotherm_data[self.isotherm_data['pressure_step']
                                                == pressure_step][parameters].values[0]
        self._verbose_print_('parameters', fitting_parameters)
        adsorbed = fitting_function(time, *fitting_parameters)

        x = [convert_time(t,'sec',time_units) for t in time]
        y = [q*convert_moles(1,'mol',mole_units)/convert_weight(1,'kg',weight_units)
             for q in adsorbed]
        axis.plot(x, y, linestyle=line_formats.next(), marker = marker_formats.next(),
                  color=colors.next(), label=label)
        axis.set_xlabel('time since start ({})'.format(time_units))
        axis.set_ylabel('equilibrium adsorbed ({}/{})'.format(mole_units,weight_units))

    def plot_kinetics(self, axis, label, pressure_step,
                          marker_formats = 'NaN', line_formats = 'NaN', colors = 'NaN',
                          time_units = 'sec', mole_units = 'mol', weight_units = 'kg'):
        from itertools import cycle
        from process_data import define_plot_colors, define_plot_lines, define_plot_markers
        from convert_units import convert_time, convert_moles, convert_weight

        # set plot formats
        if marker_formats == 'NaN':
            marker_formats = define_plot_markers()
        elif type(marker_formats) == list:
            marker_formats = cycle(marker_formats)
        else:
            print('marker formats not recognized; please input list')
        if line_formats =='NaN':
            line_formats = define_plot_lines()
        elif type(line_formats) == list:
            line_formats = cycle(line_formats)
        else:
            print('line formats not recognized; please input list')
        if colors == 'NaN':
            colors = define_plot_colors()
        elif type(colors) == list:
            colors = cycle(colors)
        else:
            print('colors not recognized; please input list')

        data = self.adsorption_data[self.adsorption_data['pressure_step'] == pressure_step]
        time = data['seconds_since_start']
        adsorbed = data['adsorbed_moles']
        # plot data
        x = [convert_time(t,'sec',time_units) for t in time]
        y = [q*convert_moles(1,'mol',mole_units)/convert_weight(1,'kg',weight_units)
             for q in adsorbed]
        axis.plot(x, y, linestyle=line_formats.next(), marker = marker_formats.next(),
                  color=colors.next(), label=label)
        axis.set_xlabel('time since start ({})'.format(time_units))
        axis.set_ylabel('equilibrium adsorbed ({}/{})'.format(mole_units,weight_units))

    def plot_kinetics_rate(self, axis, label, pressure_step,
                          marker_formats = 'NaN', line_formats = 'NaN', colors = 'NaN',
                          time_units = 'sec', mole_units = 'mol', weight_units = 'kg'):
        """
        plots rate of sorption for user specified pressure steps
        :return:
        """
        from itertools import cycle
        from process_data import define_plot_colors, define_plot_lines, define_plot_markers
        from convert_units import convert_time, convert_moles, convert_weight

        # set plot formats
        if marker_formats == 'NaN':
            marker_formats = define_plot_markers()
        elif type(marker_formats) == list:
            marker_formats = cycle(marker_formats)
        else:
            print('marker formats not recognized; please input list')
        if line_formats =='NaN':
            line_formats = define_plot_lines()
        elif type(line_formats) == list:
            line_formats = cycle(line_formats)
        else:
            print('line formats not recognized; please input list')
        if colors == 'NaN':
            colors = define_plot_colors()
        elif type(colors) == list:
            colors = cycle(colors)
        else:
            print('colors not recognized; please input list')

        # iterate through user specified pressure steps
        self._verbose_print_('pressure_step', pressure_step)
        # get sub-section of adsorption data
        data = self.adsorption_data[self.adsorption_data['pressure_step'] == pressure_step].dropna()
        time = data['seconds_since_start']
        adsorbed = data['rate_of_sorption']
        x = [convert_time(t,'sec',time_units) for t in time]
        mole_conv = convert_moles(1,'mol',mole_units)
        weight_conv = convert_weight(1,'kg',weight_units)
        time_conv = convert_time(1,'sec',time_units)
        y = [((q*mole_conv)/(weight_conv*time_conv)) for q in adsorbed]
        # plot data
        axis.semilogy(x, y, linestyle=line_formats.next(), marker = marker_formats.next(),
                  color=colors.next(), label=label)
        axis.set_xlabel('time since start ({})'.format(time_units))
        axis.set_ylabel('adsorption rate ({}/{}*{})'.\
                        format(mole_units,weight_units,time_units))

class AnalyseAdsorption(PlotAdsorption):

    def __init__(self, input_dictionary, verbose=False):
        import pandas as pd

        # read values from input dictionary
        self.raw_data = input_dictionary['adsorption_file']
        self.adsorbent_id = ''
        self.adsorptive_id = ''
        self.isotherm_id = ''
        self.probe_molecule = 'helium'
        self.eos = input_dictionary['eos']
        self.reference_cell_volume = input_dictionary['reference_cell_volume']
        self.sample_cell_volume = input_dictionary['sample_cell_volume']
        self.specific_volume = input_dictionary['sample_specific_volume']
        self.specific_volume_id = input_dictionary['specific_volume_id']
        self.leak_data = input_dictionary['leak_calibration_file']
        self.linear_leak_rate = 0
        self.poiseuille_leak_resistance = 0
        self.leak_calibration_function = input_dictionary['leak_calibration_function']
        self.equilibrium_rate_of_sopriton = input_dictionary['equilibrium_rate_of_sorption']
        self.equilibrium_rolling_window = input_dictionary['equilibrium_rolling_window']
        self.equilibrium_slope_points_required = input_dictionary['equilibrium_slope_points_required']
        self.ignore_non_equilibrium_points = input_dictionary['ignore_non_equilibrium_points']
        self.micropore_volume = input_dictionary['micropore_volume']
        self.pseudo_saturation_pressure = input_dictionary['pseudo_saturation_pressure']
        self.adsorbed_phase_density = input_dictionary['adsorbed_phase_density']
        self.absolute_adsorbed_method = input_dictionary['absolute_adsorbed_method']
        self.results_location = input_dictionary['results_location']
        self.collated_results_location = input_dictionary['collated_results_location']
        self.injected_data = pd.DataFrame()
        self.verbose(verbose)

        # perform data checks on input values read
        self._perform_data_checks_()

        # perform calibrations
        self._perform_calibrations_()

        # prepare data for calculations
        self._initialize_raw_data_()
        self._initialize_isotherm_data_()

        # calculate amount adsorbed
        self._calculate_adsorption_data_()

        # calculate equilibrium adsorbed
        self._calculate_equilibrium_data_()

        # prepare isotherm data for curve fitting
        self._initialize_isotherm_parameters_()

        # perform curve fits to find model parameters
        self.fit_kinetics_data_()
        self.fit_isotherm_data()

    def _perform_data_checks_(self):
        """
        checks if input datafiles are valid
        accepted data types include int, float, data frame, and file paths
        :return:
        """
        from os.path import exists
        from pandas.core.frame import DataFrame
        print('performing data checks\r'),
        variables_to_check = [
            self.reference_cell_volume,
            self.sample_cell_volume,
            self.specific_volume,
            self.leak_data,
            self.raw_data
        ]

        for variable in variables_to_check:
            # check data type for given variable and raise error if it's not numeric or string
            if not(type(variable) == float
                   or type(variable) == int
                   or type(variable) == DataFrame
                   or type(variable) == str):
                print('{} input type not recognized'.format(variable))
                # if reference cell volume input is string, check if file exists
            if type(self.reference_cell_volume) == str and not(exists(self.reference_cell_volume)):
                print('{} file not found'.format(variable))

    def _perform_calibrations_(self):
        """
        perform calibrations
        :return:
        """
        self._calibrate_reference_cell_volume_()
        self._calibrate_sample_cell_volume_()
        self._calibrate_sample_specific_volume_()
        self._calibrate_leak_rate_()

    def _calibrate_reference_cell_volume_(self):
        """
        calibrates reference cell volume if calibration file is given
        :return:
        """
        print('calibrating sample cell volume\r'),
        if type(self.reference_cell_volume) == str:
            pass

    def _calibrate_sample_cell_volume_(self):
        """
        calibrates sample cell volume if calibration file is given in self.sample_cell_volume
        :return:
        """
        print('calibrating sample cell volume\r'),
        if type(self.sample_cell_volume) == str:
            pass

    def _calibrate_sample_specific_volume_(self):
        """
        calibrates sample specific volume if calibration file is given in self.specific_volume
        :return:
        """
        if type(self.specific_volume) == str:
            from rig_calibration import calibrate_specific_volume
            print('calibrating sample specific volume\r'),

            # call function to calibrate specific volume
            sp_volume = calibrate_specific_volume(self.specific_volume,
                                                  self.reference_cell_volume,
                                                  self.sample_cell_volume)

            # update specific volume of class
            self.specific_volume = sp_volume

            print('Sample Specific Volume: {}'.format(self.specific_volume))

    def _calibrate_leak_rate_(self):
        """
        calibrates leak rates if calibration file is present in self.leak_data
        :return:
        """

        if type(self.leak_data) == int or type(self.leak_data) == float:
            # if leak data was given or if it was 0, read it into leak resistance
            self.linear_leak_rate = self.leak_data
            leak_pressure = 10101325 # assumed to be equal to 100 bars

        else:
            print('Leak Rate Calibration:')
            if type(self.leak_data) == str:
                from rig_calibration import calibrate_leak_rate
                self.linear_leak_rate, leak_pressure = calibrate_leak_rate(self.leak_data,self.reference_cell_volume,
                                    self.sample_cell_volume,self.specific_volume,self.eos)

        # calculate leak resistance using Poiseuille's Law Q = (pi R^4 (P1^2-P2^2))/(16 eta l P2)
        self.poiseuille_leak_resistance = -self.linear_leak_rate / abs(leak_pressure ** 2 - 101325 ** 2)

        # print fitted parameters
        templateh = '{:^30} {:^14}'
        print(templateh.format('leak parameter', 'value'))
        print(templateh.replace(':', ':-').format('', '', ''))
        templateb = '{:^30} {:^14.2e}'
        print(templateb.format('linear leak rate', self.linear_leak_rate))
        print(templateb.format('poiseuille leak resistance', self.poiseuille_leak_resistance))

    def _initialize_raw_data_(self):
        """
        process raw data for further calculations
        :return:
        """
        from convert_units import convert_pressure
        from convert_units import convert_temperature
        from convert_units import convert_weight
        import pandas as pd
        import numpy as np

        # if input is given as a file read file into dataframe
        if type(self.raw_data) == str:
            print('parsing raw data\r'),
            self.raw_data = pd.read_excel(self.raw_data)

        # reset adsorbate id, adsorptive id, and isotherm id of class from raw data
        self.adsorbent_id = self.raw_data['adsorbent_id'][0]
        self.adsorptive_id = self.raw_data['adsorptive_id'][0]
        self.isotherm_id = self.raw_data['isotherm_id'][0]

        # get seconds since start

        # get start times as the minimum reading for each pressure step
        start_times = self.raw_data[self.raw_data['pressure_tag'] == 'equilibrium'][['pressure_step', 'time']]\
            .groupby('pressure_step').min()

        self._verbose_print_('start times',start_times)

        # get first reading as the second reading recorded for each pressure step
        first_reading = self.raw_data[self.raw_data['pressure_tag'] == 'equilibrium'][['pressure_step', 'time']]\
            .groupby('pressure_step')['time'].nth(1)

        self._verbose_print_('first reading time', first_reading)

        # create a dataframe with start time and first reading
        times = pd.concat([start_times, first_reading], axis=1, sort=False)
        times.columns = ['start_times', 'first_reading']
        # get delta t for initial readings as the mean time difference between start time and first reading
        times['mean_time'] = times['first_reading'] - times['start_times']

        delta_t = times.loc[:, 'mean_time'].mean().total_seconds()

        self._verbose_print_('mean frequency for first readings', delta_t)

        # assuming adsorption begins between start time and first reading
        # calculate seconds since start as current start - start time + delta_t
        self.raw_data['seconds_since_start'] = [(t - times.loc[p, 'start_times']).total_seconds() + delta_t
                                                for (t, p) in zip(self.raw_data['time'],
                                                                  self.raw_data['pressure_step'])]

        # get delta t for all readings between 2 readings to calculate leak
        self.raw_data['delta_t'] = [x.total_seconds() for x in self.raw_data['time'].diff()]

        # convert_pressure to SI units
        print('parsing pressure\r'),
        pressure_units = self.raw_data['pressure_units'][0]
        if pressure_units == 'Pa_a':
            self.raw_data['pressure_SI'] = self.raw_data['pressure']
        else:
            self.raw_data['pressure_SI'] = self.raw_data['pressure'].apply(convert_pressure,
                                                                           args=(pressure_units, 'Pa_a'))

        # get a moving average with a window of 2 to calculate leak
        self.raw_data['mav2P'] = self.raw_data['pressure_SI'].rolling(2).mean()

        # convert_temperature to SI units
        print('parsing temperature\r'),
        temperature_units = self.raw_data['temperature_units'][0]
        if temperature_units == 'K':
            self.raw_data['temperature_SI'] = self.raw_data['temperature']
        else:
            self.raw_data['temperature_SI'] = self.raw_data['temperature']\
                .apply(convert_temperature, args=(temperature_units, 'K'))

        # convert weight to SI units
        weight_units = self.raw_data['weight_units'][0]
        if weight_units == 'kg':
            self.sample_weight = self.raw_data['weight'][0]
        else:
            self.sample_weight = convert_weight(self.raw_data['weight'][0],weight_units,'kg')

        # insert values into data frame required for moles calculation
        self.raw_data.insert(4, 'specific_volume_id', self.specific_volume_id)
        self.raw_data.insert(5, 'sample_specific_volume', self.specific_volume)
        self.raw_data.insert(3, 'eos', self.eos)
        self.raw_data.insert(6, 'reference_cell_volume', self.reference_cell_volume)
        self.raw_data.insert(7, 'sample_cell_volume', self.sample_cell_volume)
        self.raw_data.insert(8, 'sample_weight', self.sample_weight)
        self.raw_data['void_volume'] = self.raw_data['reference_cell_volume'] + self.raw_data['sample_cell_volume'] \
                                       - (self.raw_data['sample_specific_volume'] * self.raw_data['sample_weight'])

        # add leak resistance to dataframe
        self.raw_data.insert(9, 'leak_rate', self.linear_leak_rate)
        self.raw_data.insert(10, 'leak_resistance', self.poiseuille_leak_resistance)

        # calculate leak between current reading and previous reading

        if self.leak_calibration_function == 'Poiseuille':
            self.raw_data['leaked_moles'] = [delta_t * R * abs(P ** 2 - 101325 ** 2)
                                         for (delta_t, R, P) in
                                         zip(self.raw_data['delta_t'],
                                             self.raw_data['leak_resistance'],
                                             self.raw_data['pressure_SI'])]
        elif self.leak_calibration_function == 'linear':
            self.raw_data['leaked_moles'] = [-delta_t * R for (delta_t, R) in
                                             zip(self.raw_data['delta_t'],
                                                 self.raw_data['leak_rate'])]
        else:
            print('invalid leak calibration function')

        # take cumulative sum of leaked moles to get leaked moles from the start of experiment
        self.raw_data['leaked_moles'] = self.raw_data['leaked_moles'].cumsum()

        # insert rate of sorption to dataframe as NaN values
        self.raw_data.insert(10, 'rate_of_sorption', np.NaN)

    def _initialize_isotherm_data_(self):
        """
        process isotherm dataframe for further calculations
        :return:
        """
        import numpy as np
        from calculate_moles import calculate_injected_moles
        from calculate_moles import calculate_moles

        print('initializing isotherm calculation\r'),

        # only select data with initial and final pressure tags
        self.isotherm_data = self.raw_data[self.raw_data['pressure_tag'] != 'equilibrium']
        # get rid of data not used in calculation
        self.isotherm_data = self.isotherm_data[['adsorbent_id', 'adsorptive_id', 'isotherm_id', 'eos',
                                                 'adsorption_tag', 'desorption_tag',
                                                 'pressure_step', 'specific_volume_id', 'sample_specific_volume',
                                                 'void_volume', 'reference_cell_volume',
                                                 'pressure_tag', 'pressure_SI', 'temperature_SI']]
        # get mean values of pressure, temperature
        self.isotherm_data = self.isotherm_data.groupby(['adsorbent_id', 'adsorptive_id', 'isotherm_id',
                                                         'pressure_step',
                                                         'reference_cell_volume', 'specific_volume_id',
                                                         'eos', 'sample_specific_volume', 'void_volume',
                                                         'pressure_tag']).mean()
        self.isotherm_data = self.isotherm_data.unstack()
        self.isotherm_data = self.isotherm_data.reset_index()

        self.isotherm_data.columns = [''.join(col).strip() for col in self.isotherm_data.columns.values]

        self.isotherm_data['adsorption_tag'] = self.isotherm_data['adsorption_taginitial']
        self.isotherm_data['desorption_tag'] = self.isotherm_data['desorption_taginitial']

        # calculate injected moles
        self._verbose_print_('prepared_isotherm_data',self.isotherm_data)
        print('calculating injected moles\r'),
        self.isotherm_data['temp_injected_moles'] = self.isotherm_data.apply(lambda row:
                                                                             calculate_injected_moles(row),
                                                                             axis=1)
        self._verbose_print_('temp_injected_moles',self.isotherm_data['temp_injected_moles'])
        # calculate other moles present in rig before start of experiment

        other_moles = calculate_moles(self.isotherm_data['pressure_SIinitial'][0],
                                      self.isotherm_data['void_volume'].mean(),
                                      self.isotherm_data['temperature_SIinitial'][0],
                                      self.probe_molecule, 'ideal')
        self._verbose_print_('other moles', other_moles)
        # include other moles in injected moles to get gross injected moles
        self.isotherm_data['injected_moles'] = self.isotherm_data[['temp_injected_moles', 'specific_volume_id', 'eos']]\
                                                   .groupby(['specific_volume_id', 'eos']).cumsum() + other_moles

        # save injected moles dataframe
        self.injected_data = self.isotherm_data[['pressure_step', 'injected_moles']]
        self._verbose_print_('injected_data', self.injected_data)
        self._verbose_print_('injected data', self.injected_data)

        # transfer injected moles to raw data
        for index, row in self.injected_data.iterrows():
            # get index of rows in raw data where pressure step equals current pressure step
            ps = self.raw_data['pressure_step'] == row['pressure_step']
            # add injected moles to current pressure step of raw data
            self.raw_data.loc[ps, 'injected_moles'] = row['injected_moles']

        # find net injected moles as difference between gross injected moles and leaked moles
        self.raw_data['net_injected_moles'] = self.raw_data['injected_moles'] - self.raw_data['leaked_moles']

        # remove unwanted columns to initialize isotherm dataframe
        self.isotherm_data = self.isotherm_data[['adsorbent_id', 'adsorptive_id', 'isotherm_id', 'pressure_step',
                                                 'adsorption_tag', 'desorption_tag',
                                                 'specific_volume_id', 'eos']]
        self.isotherm_data.insert(8, 'equilibrium_conditions', self.equilibrium_rate_of_sopriton)
        self.isotherm_data.insert(9, 'equilibrium_reached', True)
        self.isotherm_data.insert(10, 'absolute_adsorbed_method', self.absolute_adsorbed_method)
        self.isotherm_data.insert(11, 'equilibrium_temperature', np.NaN)
        self.isotherm_data.insert(12, 'equilibrium_pressure', np.NaN)
        self.isotherm_data.insert(13, 'equilibrium_adsorbed', np.NaN)
        self.isotherm_data.insert(14, 'initial_adsorbed', np.NaN)
        self.pressure_steps = list(self.isotherm_data['pressure_step'])
        self.isotherm_data = self.isotherm_data.reset_index().set_index('pressure_step')

    def _calculate_adsorption_data_(self):
        """
        calculate amount adsorbed
        :return:
        """
        from calculate_moles import calculate_moles_row

        # only keep rows with equilibrium pressure tag
        print('calculating adsorption\r'),
        self.raw_data = self.raw_data[self.raw_data['pressure_tag'] == 'equilibrium']

        # calculate current moles
        self.raw_data['current_moles'] = self.raw_data.apply(lambda row: calculate_moles_row(row), axis=1)

        # calculate adsorbed moles as difference between injected moles and current moles
        self.raw_data['adsorbed_moles'] = (self.raw_data['net_injected_moles'] - self.raw_data[
            'current_moles']) / self.sample_weight

        self.raw_data = self.raw_data.reset_index()

    def _calculate_equilibrium_data_(self):
        """
        calculate equilibrium data
        :return:
        """
        import pandas as pd
        import statsmodels.api as sm
        from CoolProp import CoolProp as cp

        # calculate equilibrium based on experimental data

        # prepare table to print in
        print('Equilibrium Adsorbed:')
        templateh = '{:^5} {:^20} {:^20} {:^20}'
        print(templateh.format('step', 'pressure (Pa)', 'adsorbed (mol/kg)', 'equilibrium'))
        print(templateh.replace(':', ':-').format('', '', '', ''))
        templateb = '{:^5} {:^20.2e} {:^20.4f} {:^20}'

        # prepare to iterate through pressure steps
        pressure_steps = self.isotherm_data.index
        window = self.equilibrium_rolling_window
        points_needed = self.equilibrium_slope_points_required
        index_var = 'pressure_step'
        x_var = 'seconds_since_start'
        y_var = 'adsorbed_moles'

        for pressure_step in pressure_steps:
            min_index = self.raw_data[self.raw_data[index_var] == pressure_step].index.min()
            max_index = self.raw_data[self.raw_data[index_var] == pressure_step].index.max()
            num_points = 0
            equilibrium_reached = False
            equilibrium_pressure = 0
            equilibrium_adsorbed = 0
            self._verbose_print_('pressure_step',pressure_step,'min_index',min_index,'max_index',max_index)

            # fit rolling regression
            for i in range(window + min_index, max_index):
                # get y variable
                y = self.raw_data[self.raw_data[index_var] == pressure_step].loc[i - window:i][y_var]
                # get x variable
                x = self.raw_data[self.raw_data[index_var] == pressure_step].loc[i - window:i][x_var]
                # add constant to x variable
                x = sm.add_constant(x)
                # fit linear regression model
                model = sm.OLS(y, x)
                # get model results
                results = model.fit()
                # add predicted equilibrium adsorbed to dataframe
                slope = results.params[1]
                equilibrium_adsorbed = y.mean()
                equilibrium_pressure = self.raw_data['pressure_SI'].loc[i-window:i].mean()
                equilibrium_temperature = self.raw_data['temperature_SI'].loc[i-window:i].mean()
                # add slope to data frame
                self.raw_data.at[i,'rate_of_sorption']=slope
                if results.params[1] < self.equilibrium_rate_of_sopriton:
                    num_points += 1
                    equilibrium_reached = num_points > points_needed
                if equilibrium_reached:
                    break

            # add equilibrium parameters to isotherm dataframe
            self.isotherm_data.loc[pressure_step, 'equilibrium_pressure'] = equilibrium_pressure
            self.isotherm_data.loc[pressure_step, 'equilibrium_temperature'] = equilibrium_temperature
            if self.absolute_adsorbed_method == 'excess':
                self.isotherm_data.loc[pressure_step, 'equilibrium_adsorbed'] = equilibrium_adsorbed
            elif self.absolute_adsorbed_method == 'constant_volume':
                self.isotherm_data.loc[pressure_step, 'equilibrium_adsorbed'] = \
                    equilibrium_adsorbed + self.micropore_volume*self.sample_weight\
                                                *cp.PropsSI('DMOLAR','T',equilibrium_temperature,
                                                            'P',equilibrium_pressure,self.adsorptive_id)
            elif self.absolute_adsorbed_method == 'constant_density':
                self.isotherm_data.loc[pressure_step, 'equilibrium_adsorbed'] = \
                    equilibrium_adsorbed / (
                        1 - (cp.PropsSI('DMASS','T',equilibrium_temperature,'P',equilibrium_pressure,self.adsorptive_id)/
                                       self.adsorbed_phase_density)
                    )
            if pressure_step == 1:
                self.isotherm_data.loc[pressure_step, 'initial_adsorbed'] = 0
            else:
                self.isotherm_data.loc[pressure_step, 'initial_adsorbed'] = \
                    self.isotherm_data.loc[pressure_step-1, 'equilibrium_adsorbed']

            self.isotherm_data.loc[pressure_step, 'equilibrium_reached'] = equilibrium_reached

            # print equilibrium values
            print(templateb.format(pressure_step,equilibrium_pressure,equilibrium_adsorbed,str(equilibrium_reached)))

        # define adsorption data with only certain columns from raw data
        self.adsorption_data = self.raw_data[['adsorbent_id', 'adsorptive_id', 'isotherm_id', 'pressure_step',
                                              'specific_volume_id', 'eos', 'seconds_since_start', 'pressure_SI',
                                              'temperature_SI','adsorbed_moles', 'rate_of_sorption']]
        self.adsorption_data.insert(7, 'pressure', self.adsorption_data['pressure_SI'])
        self.adsorption_data.insert(7, 'temperature', self.adsorption_data['temperature_SI'])
        self.adsorption_data.drop('pressure_SI', axis='columns')
        self.adsorption_data.drop('temperature_SI', axis='columns')

        # reset index for adsorption dataframe
        self.adsorption_data.reset_index().set_index('pressure_step')

        # reset index
        # re-shape isotherm to give adsorption and desorption isotherms
        adsorption_isotherm = self.isotherm_data[self.isotherm_data['adsorption_tag'] == 1]
        adsorption_isotherm.loc[:,'adsorption_tag'] = 'adsorption'
        desorption_isotherm = self.isotherm_data[self.isotherm_data['desorption_tag'] == 1]
        desorption_isotherm.loc[:,'adsorption_tag'] = 'desorption'
        isotherm = pd.concat([adsorption_isotherm, desorption_isotherm])
        self.isotherm_data = isotherm.drop(columns=['desorption_tag'])
        self.isotherm_data = self.isotherm_data.rename(columns={'adsorption_tag': 'isotherm_tag'})
        self.isotherm_data = self.isotherm_data.drop(columns=['index']).reset_index()

    def _initialize_isotherm_parameters_(self):
        # create isotherm parameters data frame
        self.isotherm_parameters = self.isotherm_data.groupby(['adsorbent_id','adsorptive_id','isotherm_id','isotherm_tag',
                                                               'specific_volume_id','eos']).mean()
        self.isotherm_parameters = self.isotherm_parameters.reset_index()
        self.isotherm_parameters = self.isotherm_parameters[['adsorbent_id','adsorptive_id','isotherm_id','isotherm_tag',
                                                             'specific_volume_id','eos']]

    def fit_kinetics_data_(self):
        from scipy.optimize import curve_fit
        import numpy as np
        from kinetic_models import firstorder
        from kinetic_models import secondorder
        import warnings
        warnings.filterwarnings("ignore")

        print('Kinetic parameters:')
        templateh = '{:^5} {:^20} {:^25}'
        print(templateh.format('step', 'fitting function', 'fitting parameters'))
        print(templateh.replace(':', ':-').format('', '', '', ''))
        templateb = '{:^5} {:^20} {:^25}'

        kinetic_function_mapping = {
            'firstorder': [lambda t, k1: firstorder(t, q0, qe, k1), ['first_order_rate_constant']],
            'secondorder': [lambda t, k2: secondorder(t, q0, qe, k2), ['second_order_rate_constant']]
        }

        # copy row indices from isotherm data frame
        row_indices = self.pressure_steps

        self._verbose_print_('isotherm data before fitting', self.isotherm_data)

        # iterate through user input models
        for model_name in kinetic_function_mapping:
            # define fitting function and fitting parameters
            fitting_function = kinetic_function_mapping[model_name][0]
            fitting_parameters = kinetic_function_mapping[model_name][1]

            # insert columns into data frame for all model parameters
            for parameter in fitting_parameters:
                if not parameter in list(self.isotherm_data):
                    self.isotherm_data.insert(len(self.isotherm_data.columns), parameter, np.NaN)

            # iterate through pressure steps
            for row_index in row_indices:
                pressure_step = self.isotherm_data.loc[row_index]['pressure_step']
                self._verbose_print_('pressure_step', pressure_step)

                # fix modelled and constant variables
                index_boolean = self.adsorption_data['pressure_step'] == pressure_step
                x_var = self.adsorption_data[index_boolean]['seconds_since_start']
                y_var = self.adsorption_data[index_boolean]['adsorbed_moles']
                qe = self.isotherm_data.loc[row_index]['equilibrium_adsorbed']
                q0 = self.isotherm_data.loc[row_index]['initial_adsorbed']


                self._verbose_print_('index_boolean', index_boolean, 'x_var', x_var, 'y_var', y_var, 'qe', qe, 'q0', q0)

                # perform curve fitting to get fitting variables
                popt = curve_fit(fitting_function, x_var, y_var)

                self._verbose_print_('popt', popt)

                # populate data frame with fitting variables
                self.isotherm_data.loc[row_index, fitting_parameters] = popt[0]

                print(templateb.format(pressure_step,model_name,str(popt[0])))

        self._verbose_print_('isotherm data after fitting', self.isotherm_data)

    def write_to_file(self):
        import os
        import pandas as pd
        # check if director already exists, and make new directory
        if not os.path.exists(self.results_location):
            os.makedirs(self.results_location)
        # get results file name
        results_file_name = str(
            '{}_{}_{}_{}_{}_{}_{}.xlsx'.format(
                self.adsorbent_id,
                self.adsorptive_id,
                self.isotherm_id,
                self.specific_volume_id,
                self.eos,
                self.equilibrium_rate_of_sopriton,
                self.absolute_adsorbed_method
            )
        )
        results_file_name = self.results_location+results_file_name
        print('writing results to {}'.format(str(results_file_name)))
        # create excel writer
        writer = pd.ExcelWriter(results_file_name)
        self.isotherm_data.to_excel(writer,'Isotherm_data',index=False)
        self.isotherm_parameters.to_excel(writer,'Isotherm_parameters',index=False)
        self.adsorption_data.to_excel(writer,'Adsorption_data',index=False)
        self.adsorption_data.to_excel(writer,'Raw_data',index=False)
        writer.save()
        # append collated results to csv files
        # set file paths
        #raw_data_location = self.results_location+'calculated_raw_data.csv'
        #adsorption_data_location = self.results_location+'calculated_adsorption_data.csv'
        #isotherm_data_location = self.results_location+'calculated_isotherm_data.csv'
        #isotherm_parameters_location = self.results_location+'calculated_isotherm_parameters.csv'

