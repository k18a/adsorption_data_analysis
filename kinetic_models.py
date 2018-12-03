def get_kinetic_mapping_for_plots():
    kinetic_model_function = {
        'firstorder': [self._firstorder_,
                       ['initial_adsorbed', 'equilibrium_adsorbed',
                        'first_order_rate_constant']],
        'secondorder': [self._secondorder_,
                        ['initial_adsorbed', 'equilibrium_adsorbed',
                         'second_order_rate_constant']]
    }
    return kinetic_model_function

def get_kinetic_mapping_for_fits():
    kinetic_function_mapping = {
        'firstorder': [lambda t, k1: firstorder(t, q0, qe, k1), ['first_order_rate_constant']],
        'secondorder': [lambda t, k2: secondorder(t, q0, qe, k2), ['second_order_rate_constant']]
    }
    return kinetic_function_mapping

def firstorder(t, q0, qe, k1):
    """
    first order fitting function
    :param t: time
    :param q0: initial adsorbed at t = 0
    :param qe: equilibrium adsorbed at t = infinity
    :param k1: first order fitting parameter
    :return: amount adsorbed at time t
    """
    import numpy as np
    return qe - qe * np.exp(-k1 * t) + q0 * np.exp(-k1 * t)

def secondorder(t, q0, qe, k2):
    """
    second order fitting function
    :param t: time
    :param q0: initial adsorbed at time t = 0
    :param qe: equilibrium adsorbed at time t = infinity
    :param k2: second order fitting parameter
    :return: amount adsorbed at time t
    """
    return (qe * k2 * t * abs(qe - q0) + q0)/(1 + k2 * t * abs(qe - q0))

def elovich(t, q0, qe, k2):
    """
    second order fitting function
    :param t: time
    :param q0: initial adsorbed at time t = 0
    :param qe: equilibrium adsorbed at time t = infinity
    :param k2: second order fitting parameter
    :return: amount adsorbed at time t
    """
    return 1#(-log(exp))
