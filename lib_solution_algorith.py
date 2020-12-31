from sympy import *
from lib_subclass import *
from lib_mechanics import *

def main_solve(op1,op2,kval,kmetodo='default'):
    kd=str(kval)
    sR=kd+' ='
    if kmetodo=='simple_ac':
        Eq1=op2.simple_ac()-op1.simple_ac()
        kres=csolve(Eq1,kval,kdisp=False)
        display(Math(sR+latex(kres)))
        return(kres)
    
    display(Math(sR+latex(kval)))
    return(kval)
    
    
    