def get_isotherm_mapping_for_plots():
    isotherm_function_mapping = {
        'linear': [linear, ['henry_constant']],
        'langmuir': [langmuir, ['langmuir_volume','langmuir_pressure']],
        'freundlich': [freundlich, ['freundlich_k','freundlich_n']],
        'dubinin_radushkevich': [dubinin_radushkevich,
                                 ['adsorption_temperature','pseudo_saturation_pressure',
                                  'V0','dr_E']]
    }
    return isotherm_function_mapping

def get_isotherm_mapping_for_fits():
    isotherm_function_mapping = {
        'linear': [linear, ['henry_constant']],
        'langmuir': [langmuir, ['langmuir_volume','langmuir_pressure']],
        'freundlich': [freundlich, ['freundlich_k','freundlich_n']],
        'dubinin_radushkevich': [lambda P,E:
                                 dubinin_radushkevich(P,T,P0,V0,E),
                                 ['dr_E']]
    }
    return isotherm_function_mapping

def langmuir(P, VL, PL):
    """
    langmuir isotherm
    :param t: time
    :param q0: initial adsorbed at time t = 0
    :param qe: equilibrium adsorbed at time t = infinity
    :param k2: second order fitting parameter
    :return: amount adsorbed at time t
    """
    return (VL*P)/(PL+P)

def linear(P, K):
    """
    linear isotherm
    :param t: time
    :param q0: initial adsorbed at time t = 0
    :param qe: equilibrium adsorbed at time t = infinity
    :param k2: second order fitting parameter
    :return: amount adsorbed at time t
    """
    return K*P

def freundlich(P,K,n):
    """
    freundlich isotherm
    :param P:
    :param K:
    :param n:
    :return:
    """
    return K*P**(1/float(n))

def dubinin_radushkevich(P,T,P0,V0,E):
    """
    DR isotherm
    :param P:
    :param K:
    :param n:
    :return:
    """
    import numpy as np
    return V0*np.exp(-((8.314*T/E)*np.log(P0/P))**2)
