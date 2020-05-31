'''
    Model constructor of an affine model
'''

from sympy import sympify  # Concert numbers in symbpy objects
from sympy import printing  # print sympy expressions in c++
from sympy import diff  # differentiation
from sympy import Poly  # Polynomials
from sympy import polys  # Polynomials
from .modelconstructor import modelConstructorSolver


class modelConstructorAffine(modelConstructorSolver):
    def __init__(self, mn, n, m, _path='./'):
        modelConstructorSolver.__init__(self, mn, n, m, _path)
        self.aus = []
        self.gs = []
        self.hs = []
        self.g0 = sympify(0)
        self.h0s = []
        self.adim = 0

    def init(self, f, rc, fc):
        modelConstructorSolver.init(self, f, rc, fc)
        res = []

        # Here we compute in which components of the
        # control u the vector field and the running cost
        # are linear, the vector res will have
        # the value of one for each affine control
        for u in self.control_symbols_:
            res.append(1)
            c0 = Poly(self.f0, u).all_coeffs()
            for fi in self.funcion_symbols_:
                try:
                    cs = Poly(fi, u).all_coeffs()
                except polys.polyerrors.PolynomialError:
                    cs = [0, 0, 0, 0]
                n = max(len(cs), len(c0))
                if (n > 2):
                    res[-1] = 0

        for i in res:
            if (i == 1):
                self.adim = self.adim + 1

        print('affine vector: ' + str(res))

        self.gs = self.funcion_symbols_[:]
        for i, gi in zip(range(0, self.state_dim_), self.gs):
            for j, uj in zip(
                    range(0, self.control_dim_), self.control_symbols_):
                if (res[j]):
                    cs = Poly(gi, uj).all_coeffs()
                    gi = cs[-1]
                    if (len(cs) == 2):
                        self.hs.append(cs[0])
                    else:
                        self.hs.append(sympify(0))
            self.gs[i] = gi

        self.g0 = self.f0
        for j, uj in zip(range(0, self.control_dim_), self.control_symbols_):
            if (res[j]):
                cs = Poly(self.g0, uj).all_coeffs()
                self.g0 = cs[-1]
                if (len(cs) == 2):
                    self.h0s.append(cs[0])
                else:
                    self.h0s.append(sympify(0))

    def writeAffineFuncs(self):
        self.output_file_.write("""
    void dynSysAffinePart
        (double t,const double x[],const double u2[],double g[]){\n""")
        i = 0
        for gi in self.gs:
            self.output_file_.write("\t\tg[%(n)i]= %(c)s;\n" % {
                'n': i,
                'c': printing.ccode(gi)
            })
            i = i + 1

        self.output_file_.write("\t}\n")

        self.output_file_.write("""
        void dynSysLinearPart(double t,const double x[],double h[]){\n""")
        counter = 0
        for ex in self.hs:
            self.output_file_.write("\t\th[%(n)i]= %(c)s;" % {
                'n': counter,
                'c': printing.ccode(ex.evalf())
            })
            if (counter % self.control_dim_):
                self.output_file_.write("// ---- end of row\n")
            else:
                self.output_file_.write("\n")

            counter = counter + 1
        self.output_file_.write("\t}\n")

        self.output_file_.write("""
        void dynSysLinearPart_1G(double t,const double *x,double dhdx[]){\n""")
        counter = 0
        for x in self.state_symbols_:
            for ex in self.hs:
                self.output_file_.write(
                    "\t\tdhdx[%(n)i]= %(c)s;" % {
                        'n': counter,
                        'c': printing.ccode(diff(ex, x).evalf())
                    })
                if (counter % self.adim):
                    self.output_file_.write("// ----- end of row\n")
                else:
                    self.output_file_.write("\n")
                counter = counter + 1
            self.output_file_.write("//--end of matrix\n\n")
        self.output_file_.write("\t}\n")

        self.output_file_.write("""
        void dynSysLinearPart_dt(double t,const double *x,double dhdt[]){\n""")
        counter = 0
        for ex in self.hs:
            self.output_file_.write(
                "\t\tdhdt[%(n)i]= %(c)s;\n" % {
                    'n': counter,
                    'c': printing.ccode(diff(ex, self.time_symbol_).evalf())
                })
            counter = counter + 1
        self.output_file_.write("\t}\n")

        self.output_file_.write("""
        double runningCostAffinePart(double t,const double *x,const double u2[]){
               return  %s;
        }\n""" % printing.ccode(self.g0.evalf()))

        self.output_file_.write("""
    void runningCostLinearPart(double t,const double x[],double h0[]){\n""")
        counter = 0
        for ex in self.h0s:
            self.output_file_.write("\t\th0[%(n)i]= %(c)s;\n" % {
                'n': counter,
                'c': printing.ccode(ex.evalf())
            })
            counter = counter + 1
        self.output_file_.write("\t}\n")

        self.output_file_.write("""
    void runningCostLinearPart_1G(double t,const double x[],double dh0dx[]){\n
""")
        counter = 0
        for x in self.state_symbols_:
            for ex in self.h0s:
                self.output_file_.write(
                    "\t\tdh0dx[%(n)i]= %(c)s;\n" % {
                        'n': counter,
                        'c': printing.ccode(diff(ex, x).evalf())
                    })
                counter = counter + 1
        self.output_file_.write("\t}\n")

        self.output_file_.write("""
    void runningCostLinearPart_dt(double t,const double *x,double h0t[]){\n""")
        counter = 0
        for ex in self.h0s:
            self.output_file_.write(
                "\t\th0t[%(n)i]= %(c)s;\n" % {
                    'n': counter,
                    'c': printing.ccode(diff(ex, self.time_symbol_).evalf())
                })
            counter = counter + 1
        self.output_file_.write("\t}\n")

    def create(self):
        self.writeClassHeader()
        self.writeSolverFuncs()
        self.writeAuxFuncs()
        self.writeAffineFuncs()
        self.writeClassClosure()

    def writeClassHeader(self):
        self.output_file_.write("""
#ifndef """ + self.class_name_.upper() + """_H
#define """ + self.class_name_.upper() + """_H

#include<affinesys.h>
class c""" + self.class_name_ + """:public affinesys{
    public:
    double x0[%(dim)i];
    double umin[%(adim)i];
    double umax[%(adim)i];\n""" % {
            'dim': self.state_dim_,
            'adim': self.adim
        })

        for i in self.aux_vars_symbols_:
            self.output_file_.write("\tdouble %s;\n" % i.name)
        self.output_file_.write(
            """
    c%(cn)s():\n\taffinesys(%(dim)i,%(cdim)i,%(adim)i,1.0e-3,20.0,1.0e-7,x0)"""
            % {
                'cn': self.class_name_,
                'dim': self.state_dim_,
                'cdim': self.control_dim_,
                'adim': self.adim
            })
        for i in self.aux_vars_symbols_:
            self.output_file_.write(",\n\t\t%s ( 1.0 ) " % i.name)

        self.output_file_.write("\n\t{\n")
        for i in range(self.state_dim_):
            self.output_file_.write("\t\tx0[%i]= 0.0;\n" % i)
        self.output_file_.write("\n\t}\n")
