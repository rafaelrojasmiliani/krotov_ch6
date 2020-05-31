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


def main():
    mc = modelConstructorAffine('Eh1DofMb', 2, 1)

    state = mc.getStateVar()
    u = mc.getConstrolVar()
    t = mc.getTimeVar()

    xi = state[0]
    xid = state[1]

    zeta0 = symbols('zeta0')
    omega = symbols('Omega')
    y = symbols('etadd', cls=Function)
    f = [0, 0]

    f[0] = xid
    f[1] = -xi - 2 * (u[0] + zeta0) * xid + y(t)

    power = xid**2 * u[0]

    rc = -power

    fc = 0

    mc.init(f, rc, fc)

    mc.append_aux_func_t('etadd', y(t))
    mc.appendAuxVar(zeta0)
    mc.appendAuxVar(omega)
    mc.create()


if __name__ == '__main__':
    main()
