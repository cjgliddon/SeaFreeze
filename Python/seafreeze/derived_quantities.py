# Methods for computing additional thermodynamic quantities not directly represented by the Gibbs function

from collections import namedtuple
import numpy as np
import os.path as op
import seafreeze as sf

defpath = op.join(op.dirname(op.abspath(__file__)), 'SeaFreeze_Gibbs_VII_NaCl.mat')

def getTheta(PTm, phase, path=defpath, *tdvSpec):
    ''' Calculates the potential temperature under given conditions of pressure, absolute
        temperature, and salt molality at the specified phase.

    :param PTm:         The pressure (MPa), temperature (K), and molality (mol solute/kg solvent) at which potential
                        temperature should be calculated. May have either scatter-type or grid-type input, as in the
                        main seafreeze methods (see seafreeze.py for details). If phase is not a salt, molality is
                        not a required parameter.
    :param phase:       The phase of the substance in question. Should be one of the keys for the "phases" dictionary
                        in seafreeze.py (as a string). 

    '''

def _find_root_MNR(func, func_deriv, x0, tolerance='default'):
    ''' Helper function for root-finding according to the modified Newton-Raphson method described in McDougall &
        Wotherspoon (2014), which is used to solve for the potential temperature.

    :param func:        A function that accepts a single real number (i.e., float) as input and returns a single real
                        number as output. The function should (in principle) be continuously defined on the interval
                        in question.
    :param func_deriv:  Also a function that accepts and returns a single float, representing the derivative of
                        `func`.  
    :param x0:          The initial guess for the value of the function root.
    :param tolerance:   Optional parameter that specifies the tolerance for the algorithm; i.e., the smallest value
                        of the difference between root estimates at successive iterations at which the algorithm
                        will terminate. If set to 'default', the tolerance will be the machine "resolution" of the
                        datatype (which is 1E-15 for `float64`).
    '''

    if tolerance == 'default':
        tolerance = np.finfo(np.float64).resolution

    x0_int = x0 - func(x0)/func_deriv(x0)


