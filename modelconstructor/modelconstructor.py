'''
    A class for automatic differentiation and generation of
    c++ class with the derivatives of a system
'''
from sympy import sympify  # Concert numbers in symbpy objects
from sympy import symbols  # create sympy symbols
from sympy import printing  # print sympy expressions in c++
from sympy import diff  # differentiation
from sympy import MatrixSymbol  # Matrices


class modelConstructorSolver(object):
    ''' Class intended to construct the partial derivative of the
    vector field required for the Krotov's method
    '''

    def __init__(self, _class_name, _n, _m, _path='./'):
        self.state_dim_ = _n
        self.control_dim_ = _m
        self.class_name_ = _class_name
        self.state_symbols_ = []
        self.control_symbols_ = []
        self.control0_values_ = []
        self.control_max_ = []
        self.control_min_ = []
        self.funcion_symbols_ = []
        self.lagrange_symbols_ = []
        self.aux_func_symbols_ = {}
        self.aux_func_t_symbols_ = {}
        self.aux_vars_symbols_ = []
        self.aux_vars_values_ = []

        f0 = sympify(0)
        F = sympify(0)
        self.hamiltonian_ = sympify(0)
        fileName = _class_name.lower() + '.h'
        self.output_file_ = open(_path + fileName.lower(), "w")
        self.time_symbol_ = symbols('t')
        for i in range(0, self.state_dim_):
            self.state_symbols_.append(symbols('x[' + str(i) + ']', real=True))
            self.funcion_symbols_.append(sympify(0))
            self.lagrange_symbols_.append(
                symbols('lambda[' + str(i) + ']', real=True))
        for i in range(0, self.control_dim_):
            self.control_symbols_.append(
                symbols('u[' + str(i) + ']', real=True))
            self.control0_values_.append(sympify(0.0))
            self.control_min_.append(sympify(0.0))
            self.control_max_.append(sympify(0.0))

    def set_control0_val(self, _idx, _ex):
        self.control0_values_[_idx] = _ex

    def set_control_limits(self, _idx, _min, _max):
        self.control_max_[_idx] = _max
        self.control_min_[_idx] = _min

    def getStateVar(self):
        return self.state_symbols_

    def getConstrolVar(self):
        return self.control_symbols_

    def getTimeVar(self):
        return self.time_symbol_

    def init(self, f, rc, fc):
        self.funcion_symbols_ = sympify(f)
        self.f0 = sympify(rc)
        self.F = sympify(fc)
        for i in range(0, self.state_dim_):
            self.hamiltonian_ = self.hamiltonian_ + \
                self.lagrange_symbols_[i] * self.funcion_symbols_[i]
        self.hamiltonian_ = self.hamiltonian_ - self.f0

    def getHamiltonian(self):
        return self.hamiltonian_

    def appendAuxVar(self, var, _val=sympify(0.0)):
        self.aux_vars_symbols_.append(var)
        self.aux_vars_values_.append(var)

    def appendAuxFunction(self, name, ex):
        self.aux_func_symbols_[name] = ex

    def append_aux_func_t(self, name, ex):
        self.aux_func_t_symbols_[name] = ex

    def create(self):
        self.writeClassHeader()
        self.writeSolverFuncs()
        self.write_void_max_hamiltonian()
        self.writeAuxFuncs()
        self.writeClassClosure()

    def writeClassHeader(self):
        self.output_file_.write("""
#ifndef """ + self.class_name_.upper() + """_H
#define """ + self.class_name_.upper() + """_H

#include"../solver.h"
class c""" + self.class_name_ + """:public solver{
    public:\n""")

        for i in self.aux_vars_symbols_:
            self.output_file_.write("\tdouble %s;\n" % i.name)
        self.output_file_.write(
            """
    c%(cn)s():\n\t\tsolver(%(dim)i,%(cdim)i,1.0e-3,20.0,x0)""" % {
                'cn': self.class_name_,
                'dim': self.state_dim_,
                'cdim': self.control_dim_
            })
        self.output_file_.write("\nu0(0.01)")
        for i in self.aux_vars_symbols_:
            self.output_file_.write(",\n\t\t%s ( 1.0 ) " % i.name)

        self.output_file_.write("\n\t{\n\t}\n")

    def writeClassClosure(self):
        self.output_file_.write("""
};
#endif
    """)
        self.output_file_.close()

    def write_void_max_hamiltonian(self):
        self.output_file_.write("""
    void maxHamiltonian(double t,const double x[],
        const double lambda[],double u[]){
    }\n""")

    def writeSolverFuncs(self):
        self.output_file_.write("""
    virtual void controller0(double t,const double x[],double u[]){
""")
        for u, uval in zip(self.control_symbols_, self.control0_values_):
            self.output_file_.write('\t\t{:s} = {:s};'.format(
                printing.ccode(u), printing.ccode(uval)))

        self.output_file_.write("""
    }\n\n""")
        self.output_file_.write("\tint dynSystem(double t,const double x[],"
                                " const double u[], double xd[]){\n")

        for i, ex in zip(range(0, self.state_dim_), self.funcion_symbols_):
            self.output_file_.write("\t\txd[%(n)i]= %(c)s;\n" % {
                'n': i,
                'c': printing.ccode((ex.evalf()))
            })
        self.output_file_.write("\t\treturn GSL_SUCCESS;\n\t}\n")

        self.output_file_.write("\tint dynSys_1Gx(double t,const double x[],"
                                "const double u[],double dfdx[]){\n")

        for i, fi in zip(range(0, self.state_dim_), self.funcion_symbols_):
            for j, sj in zip(range(0, self.state_dim_), self.state_symbols_):
                self.output_file_.write(
                    "\t\tdfdx[%(n)i]= %(c)s;\n" % {
                        'n': i * self.state_dim_ + j,
                        'c': printing.ccode((diff(fi, sj)).evalf())
                    })
        self.output_file_.write("\t\treturn GSL_SUCCESS;\n\t}")

        self.output_file_.write("""
    void dynSys_2Gxx(double t,const double x[],const double u[],
                                                             double dfdxx[]){\n"""
                                )
        counter = 0
        for i in range(0, self.state_dim_):
            for j in range(0, self.state_dim_):
                for k in range(j, self.state_dim_):
                    self.output_file_.write(
                        '\t\t dfdxx[' + str(counter) + ']= ')
                    self.output_file_.write(
                        printing.ccode((diff(
                            self.funcion_symbols_[i], self.state_symbols_[j],
                            self.state_symbols_[k])).evalf()) + ';\n')
                    counter = counter + 1
            self.output_file_.write("\n")

        self.output_file_.write("\t}\n")
        self.output_file_.write("""
    double runningCost(double t,const double x[],const double u[]){
        return %s;
    }
    """ % printing.ccode(self.f0.evalf()))

        self.output_file_.write("""
    void runningCost_1Gx(double t,const double x[],const double u[],
                                double df0dx[]){\n""")
        for i in range(0, self.state_dim_):
            self.output_file_.write(
                "\t\tdf0dx[%(n)d]= %(c)s;\n" % {
                    'n':
                    i,
                    'c':
                    printing.ccode(
                        (diff(self.f0, self.state_symbols_[i])).evalf())
                })
        self.output_file_.write("\t}\n")
        self.output_file_.write("""
    void runningCost_2Gxx(double t,const double x[],const double u[],
                                double df0dxx[]){\n""")

        counter = 0
        for i in range(0, self.state_dim_):
            for j in range(i, self.state_dim_):
                self.output_file_.write(
                    "\t\tdf0dxx[%(n)d]= %(c)s;\n" % {
                        'n':
                        counter,
                        'c':
                        printing.ccode(
                            diff(self.f0, self.state_symbols_[i],
                                 self.state_symbols_[j]).evalf())
                    })
                counter = counter + 1
        self.output_file_.write("\t}\n")

        self.output_file_.write("""
    double finalCost(double t,const double x[],const double u[]){
        return %s;
    }
    """ % printing.ccode(self.F.evalf()))

        self.output_file_.write("""
    void finalCost_1Gx(double t,const double x[],const double u[],
                                                                    double *dFdx){\n"""
                                )
        for i in range(0, self.state_dim_):
            self.output_file_.write(
                "\t\tdFdx[%(n)d]= %(c)s;\n" % {
                    'n':
                    i,
                    'c':
                    printing.ccode(
                        diff(self.F, self.state_symbols_[i]).evalf())
                })
        self.output_file_.write("\t}\n")

        self.output_file_.write("""
    void finalCost_2Gxx(double t,const double x[],const double u[],
                               double *dFdxx){\n""")
        counter = 0
        for i in range(0, self.state_dim_):
            for j in range(i, self.state_dim_):
                self.output_file_.write(
                    "\t\tdFdxx[%(n)d]= %(c)s;\n" % {
                        'n':
                        counter,
                        'c':
                        printing.ccode((diff(self.F, self.state_symbols_[i],
                                             self.state_symbols_[j])).evalf())
                    })
                counter = counter + 1
        self.output_file_.write("\t}\n")

        self.output_file_.write("""
    double hamiltonian(double t,const double x[], const double lambda[],const double u[]){
        return %s;
    }\n""" % printing.ccode(self.hamiltonian_))

    def writeAuxFuncs(self):
        for name, ex in self.aux_func_symbols_.items():
            self.output_file_.write("""
    void %s(double t,const double x[],const double lambda[],
                                        const double u[],double res[]){
            """ % name)
            self.output_file_.write(
                printing.ccode(ex, MatrixSymbol('res', ex.shape[0],
                                                ex.shape[1])))
            self.output_file_.write('\n\t}\n')

        for name, ex in self.aux_func_t_symbols_.items():
            self.output_file_.write("""
    double {:s}(double t){{
        return {:s};
            """.format(name, printing.ccode(ex)))
            self.output_file_.write('\n\t}\n')
