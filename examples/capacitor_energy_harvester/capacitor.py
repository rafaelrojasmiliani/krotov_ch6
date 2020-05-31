'''
    In this file we create the c++ class which represents a capacitor
'''
import os
import unittest
import functools
import traceback
import sys
import pdb
from sympy import symbols
from sympy import Function
from sympy import sympify
from sympy import sqrt
from sympy import solve
from sympy import Matrix

import pathlib

MOD_DIR = pathlib.Path(__file__).parent.absolute().parents[1]

sys.path.append(str(MOD_DIR))
from modelconstructor.affine import modelConstructorAffine


def debug_on(*exceptions):
    ''' Decorator for entering in debug mode after exceptions '''
    if not exceptions:
        exceptions = (Exception, )

    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            try:
                return f(*args, **kwargs)
            except exceptions:
                info = sys.exc_info()
                traceback.print_exception(*info)
                pdb.post_mortem(info[2])
                sys.exit(1)

        return wrapper

    return decorator


class cMyTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(cMyTest, self).__init__(*args, **kwargs)
        self.dir_path_ = os.path.dirname(os.path.realpath(__file__))+'/'

    @debug_on()
    def test(self):
        mc = modelConstructorAffine('EhCapacitor', 3, 2, _path=self.dir_path_)

        af = {}
        state = mc.getStateVar()
        u = mc.getConstrolVar()
        t = mc.getTimeVar()

        x = state[0]
        v = state[1]
        q = state[2]

        zeta0,  ebar, Abar, dbar =\
            symbols('zeta0 ebar Abar dbar', positive=True)
        y = symbols('y', cls=Function)
        mc.append_aux_func_t('y', y(t))
        f = [0, 0, 0]

        f[0] = v
        f[1] = -x - 2.0 * zeta0 * v + q**2 / (ebar * Abar) - y(t)
        f[2] = -(u[0] + u[1]) * q * (1.0 - x) / (Abar * ebar) + u[0]

        Pd = (1.0 - q * (1.0 - x) /
              (Abar * ebar))**2 * u[0] - ((1.0 - x) /
                                          (Abar * ebar) * q)**2 * u[1]

        rc = Pd

        fc = sympify(0)

        l = symbols('l')
        u0, u1 = symbols(r'u[0] u[1]', positive=True)
        eq = (u0 + u1) / sqrt(ebar * Abar) * l**3 - (
            u0 + u1) / sqrt(ebar * Abar) * l + u0
        xsqrt = solve(eq, l)[2]
        mc.init(f, rc, fc)

        mc.appendAuxFunction('equilibrium_points',
                             Matrix([xsqrt**2, 0.0,
                                     sqrt(Abar * ebar) * xsqrt]))
        mc.appendAuxVar(zeta0)
        mc.appendAuxVar(ebar)
        mc.appendAuxVar(Abar)
        mc.create()


if __name__ == '__main__':
    unittest.main()
