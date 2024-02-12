

from  libaldo_math2 import *
from  libaldo_algorith import *
from  mathexponencial  import *
from  mathbasic import * 
from sympy import *
## PARTIAL FRACTION...

def primesymboldiff(var):
    svar=str(var)
    dsvar1=svar+"'"
    dsvar2=svar+"''"
    kres=  dsvar1+' '+dsvar2
    dsvar1=symbols(dsvar1)
    dsvar2=symbols(dsvar2)
    return dsvar1,dsvar2
    
def functiondiff(QQ,var,var2):
    ksym=QQ.ksym
    f=Function(str(var2))(var)
    ksym2=subsubs(ksym,var2,f)
    kres=diff(ksym2,var)
    return kres

def sympydiff2prime(expr,var,var2):
    Y=symbols(str(var2)) 
    p1=str(expr)   
    svar=str(var2) 
    dy,dy2=primesymboldiff(var2)

    p2='d'+str(var2)
    p3=p2+'2'

    F=Function(str(var2))(var)
    df=F.diff(var)
    df2=F.diff(var,var)
    p4=str(df)
    p5=str(df2)
    p1=p1.replace(p5,p3)
    p1=p1.replace(p4,p2)
    p1=p1.replace(str(F),svar)
    fres=eval(p1)
    fres=fres.subs(p3,dy2)
    fres=fres.subs(p2,dy)
    fres=fres.subs(F,Y)
    return fres 


def primediff2sympy(expr,var,var2):
    dp2=str(var2)+"''"
    dp1=str(var2)+"'"
    yvar=str(var2)
    xvar=str(var)
    f=Function(yvar)(var)
    df=diff(f,var)
    df2=diff(f,var,var)
    sf=yvar+'('+xvar+')'
    sdf=str(df)
    sdf2=str(df2)
    sexpr=str(expr)
    sexpr=sexpr.replace(dp2,'D2')
    sexpr=sexpr.replace(dp1,'D1')
    sexpr=sexpr.replace(yvar,sf)
    sexpr=sexpr.replace('D2',sdf2)
    sexpr=sexpr.replace('D1',sdf)
    return parse_expr(sexpr)


def changedview(expr,var1,var2):
    sexpr=str(expr)
    if "'" in sexpr:
         return primediff2sympy(expr,var1,var2)
    else:
         return sympydiff2prime(expr,var1,var2) 
        