 

from sympy import *
 
from IPython.display import Math  # ,display
from matplotlib.pyplot import ylabel, plot, show, xlabel, title
from libaldo_math2 import *
from libaldo_algorith import *
 
import copy
import numpy as np
import pandas as pd
 
import matplotlib.pyplot as plt
 
 
 
# from lib_Func import *
import dill  # pip install dill --user

filename = 'workspace.pkl'


def savework():
    dill.dump_session(filename='workspace.pkl')


def loadwork():
    dill.load_session(filename='workspace.pkl')


# and to load the session again:




def show_main():
    for i in dataQ:
        i.s()



def upBagSys(ksys, Bag, kope=''):
    for i in ksys:
        i.upBag(Bag, kope=kope)



###########################################
#               END MyIteger Class
###########################################


def get_real_value(ee):
    if type(ee) == MyEq:
        kres = ee.ksym
    else:
        kres = ee
    return kres



###############################################
#  Mass center Inertia

def pQ(mm, vv, kope=''):
    rho = symbols('rho')

    kres = mm / vv
    sE([rho, '=', kres])
    return kres


#################################################
#   Solve Algorithm

def solved(*args,**kwargs):
    '''
    solved (var,exp1,exp2,**kwargs)
        input 
            var : variable to find
            exp1 : math expre or MyEq class that is  equation equal= 0
            exp2 (optional): math expre or MyEq class  
            if exp2 is given then the Eq to evalue is expr1 -expr2
            kwargs : conditions to evalue, example x=1, t=0..etc
        return
            return a MyEq class with name str(var)
    '''        
    var=args[0]
    ee=args[1]
    if type(ee)==MyEq:
        ee=ee.ksym
    keq=ee    
    if len(args)==3:
        ee2=args[2]
        if type(ee2)==MyEq:
            ee2=ee2.ksym
        keq=ee-ee2
    if len(kwargs)>0:             
            keq=real_subs(keq,**kwargs) 
              
    kres=solve(keq,var)
    if type(kres)==list:
        kres=kres[0]
                
    kname=str(var)
    ee0=MyEq(kres,kname=kname)
    return ee0       
    
def solverSys(*args, Bag=''):
    Ke = []
    Kv = []
    Kn = []

    for i in args:
        if type(i) == MyEq:
            if Bag != '':
                i.upBag(Bag, kshow=False)
            Ke.append(i)
        if type(i) == Symbol:
            Kv.append(i)
            Kn.append(i.name)
    # return(Ke,Kv,Kn)

    return MyEqSolveLin(Ke, Kv, Kn, Bag=Bag)


def MyEqSolveLin(Ke, Kv, Kn, Bag=''):  # Solve n MyEq with n unknow variable
    '''
    Example
        Ke=[e2,e2,e0]  MyEqs Matrix
        Kv=[N1,T,a]    unKnow Vriables
        Kn=['N_1','T','a_c']  New Name Eq

        N11,T1,ac = MyEqSolveLin(Ke,Kv,Kn)
        returns resepective  answer
    '''
    vecs = []
    qq = len(Ke)
    kres = []
    for i in range(qq):
        ee = Ke[i]
        ksym = Kv[i]
        ks = ee.solve(ksym,kshow=False)
        if type(ks) == list:
            rr = max(ks)
            ks = rr

        vecs.append(ks)
        Ker = Ke[i + 1::]
        for e1 in Ker:
            e1.set(ksym, ks, kshow=False)
            e1.reduFac(kshow=False)
            e1.simplify(kshow=False)

    for i, kname in zip(vecs, Kn):
        ee = MyEq(i, kname, kshow=False)
        kres.append(ee)
    ueq = kres[-1]
    ksym = ueq()
    vsym = Kv[-1]
    for ee in kres[0:-1]:
        ee.set(vsym, ksym, kshow=False)
        ee.reduFac(kshow=False)
        ee.simplify(kshow=False)
    for i in kres:
        i.s()
    return kres


def Solve2Eq(ksym=[], kvar=[], knom=[], kope=''):
    e1, e2 = ksym
    v1, v2 = kvar
    t1, t2 = knom

    r1 = e1.solve(v1)
    e2.set(v1, r1, kshow=False)
    r2 = e2.solve(v2)
    r2 = opemat(r2, kope=kope)
    e1.set(v2, r2, kshow=False)
    r1 = e1.solve(v1)
    r1 = opemat(r1, kope=kope)
    aa = MyEq(r1, t1)
    bb = MyEq(r2, t2)
    return (aa, bb)


def Diff(ksym, kvar, kname=''):
    kres = ksym
    kres = kres.diff(kvar)
    if kname == '':
        return kres
    else:
        return MyEq(kres, kname)


def Diff2(ksym, kvar, kname=''):
    kres = ksym
    kres = kres.diff(kvar)
    kres = kres.diff(kvar)
    if kname == '':
        return kres
    else:
        return MyEq(kres, kname)



def upBag2sys(vecEq, kBag):
    for i in vecEq:
        i.upBag(kBag)




def eQSolver(*args):
    vec1 = []
    uk1 = []
    for i in args:
        if type(i) == list:
            for j in i:
                if type(j) == MyEq:
                    vec1.append(j())
                elif fpoly(j, 'n') > 1:
                    vec1.append(j)
                else:
                    uk1.append(j)
        else:
            if type(i) == MyEq:
                vec1.append(i())
            elif fpoly(i, 'n') > 1:
                vec1.append(i)
            else:
                uk1.append(i)

    vec2 = []
    kres = []
    for i in vec1:
        if type(i) == MyEq:
            vec2.append(i())
        else:
            vec2.append(i)

    mm = solve(vec2, uk1)
    if type(mm) == dict:
        kk, vv = kunpakDic(mm)

        for i, j in zip(kk, vv):
            kres.append(MyEq(j, i))
        return kres
    else:
        for i, j in zip(mm[0], uk1):
            j = MyEq(i, str(j))
            kres.append(j)
        return (kres)


def solvelin(*args, kope='', Eq=True):  # solveLinearSys(e1,e2,mu1,mu2)
    mS = []
    mV = []

    for i in args:
        if type(i) == MyEq:
            mS.append(i())
        elif type(i) == eQ:
            ee = MyEq(i.ksym, kname=i.name, kshow=False)
            mS.append(ee())
        elif type(i) == str:
            kope = i
        else:
            mV.append(i)
    kres = solve(mS, mV)

    kk, vv = kunpakDic(kres)
    if kope != '':
        vv = opemat_vec(vv, kope)
    if Eq:
        EqR = []
        for i, j in zip(kk, vv):
            EqR.append(MyEq(j, i))
        return EqR

    else:
        for i, j in zip(kk, vv):
            sE([i, '=', opemat(j, kope=kope)])
    if Eq:
        kres = []
        for i, j in zip(kk, vv):
            kres.append(MyEq(opemat(j, kope=kope), i, kshow=False))
        return kres

    return vv


def get_squareMono(ksym):
    if type(ksym) == MyEq:
        ksym = ksym.ksym
    kres = ksym
    mm = fpoly(ksym, 'list')
    mr = []
    ms = []
    rr = []
    centra = 0
    ksigno = 1
    for i in mm:
        mr.append(opemat(rpow(i, 2), 'r'))
        ms.append(str(opemat(rpow(i, 2), 'r')))
    for i, j, k in zip(ms, mr, mm):
        if 'sqrt' in i:
            central = k
            if '-' in str(central):
                ksigno = -1
        else:
            rr.append(j)
    if len(rr) == 2:
        kres = kpow(rr[1] + ksigno * rr[0], 2)
    return kres


#######
def expand2MyEq(ee):
    ktype = ee.type
    var2 = ee.var2
    mm = ee.list()
    cc = 1
    kname = ee.name
    kres = []
    for i in mm:
        nname = kname + str(cc)
        nname = MyEq(i, nname, var2=var2)
        kres.append(nname)
        cc += 1
    return kres


def upgrade(*args, kshow=True, andsolve=[]):
    if andsolve != []:
        vv = andsolve[0]
        ee = andsolve[1]
        vv = ee.solve(parse_expr(vv), vv, kope='s')
    eev = []
    evv = []
    for i in args:
        if type(i) == MyEq:
            if i.type == 'C':
                eev.append(i)
            else:
                evv.append(i)

    for i in eev:
        for j in evv:
            try:
                i.upgrade(j, kshow=False)
            except:
                pass
    for i in eev:
        if i.ksym != 0:
            i.s()
    for i in evv:
        if type(i) == MyEq:
            i.simplify(kshow=False)
            if i.ksym != 0:
                i.s()


def upgradeList(*args, kshow=True, andsolve=[], kope='s'):
    eev = []
    evv = []
    for i in args:
        if type(i) == MyEq:
            if i.type == 'C':
                if i != andsolve[1]:
                    eev.append(i)

    if andsolve != []:
        vv = andsolve[0]
        ee = andsolve[1]
        kres = ee.solve(vv)
        kres = opemat(kres, 's')
        vv = MyEq(kres, str(vv), ktype='C', kshow=False)
        ee.type = 'P'

    for i in eev:
        if i.type == 'C':
            i.upgrade(vv, kshow=False)
            i.simplify(kshow=False)

    for i in eev:
        if i.ksym != 0:
            i.s()
    vv.s()
    return vv



def func_sig(kf, x1, x2, var=x):
    ee = MyEq(kf, var2=var, kshow=False)
    xx = (x2 - x1) / 2
    return ee(xx)


def get_intersec_2func(y1, y2, var=x):  # y1(x), y2(x), return intersec y1 and y2
    ee = MyEq(y1 - y2, kshow=False)  # return vector
    return ee.solve(var)


def reset_ee(*args):
    eeFull = []
    for i in args:
        i.init = False


def Upgrade(*args, kope='', kshow=True):
    newa = []
    for i in args:
        if type(i) == MyEq:
             
            if i.ksym != 0:
                newa.append(i)
    args = newa
    antes = []
    for i in args:
        antes.append(str(i))
    qq = len(args)
    for i in range(qq):

        mm = []
        for j in range(qq):
            if j != i:
                mm.append(args[j])
        args[i].upgrade(mm, kshow=False, kope=kope)
    if kshow:
        for i, j in zip(args, antes):
            if str(i) != newa:
                if i.ksym != 0:
                    i.s()


def presolve(ksym, val):
    kres = solve(ksym, val)
    if kres == []:
        try:
            kres = solve(opemat(ksym, 'esf'), val)
            if kres != []:
                return kres
            else:
                ksym = factorSec(ksym, val)
                kres = solve(opemat(ksym, 'esf'), val)
                if kres != []:
                    return kres
        except:
            done = False
    return kres


def eQsolve(ksym, kname, kope=''):
    kval = parse_expr(kname)
    kres = csolve(ksym, kval)
    kres = opemat(kres, kope)
    kval = MyEq(kres, kname)
    return kval
    
def Qsolve(*args):
    '''
        N1,mu1=Qsolve(FxA,FyA,N1,mu1)
    '''    
    eqq=[]
    evv=[]
    for i in args:
        if type(i)==MyEq:
            eqq.append(i)
        else:
            evv.append(i)
    kres=solve(eqq,evv)
    try:
        vsym,vval=kunpakDic(kres)
    except:
        vsym=evv
        vval=list(kres[0])
         
    vres=[]
    for i ,j in zip(vsym,vval):
        kname=str(i)
        ee=MyEq(j,kname=i)
        vres.append(ee)
    return vres   



def Diff2flat(kres,kvar,var2): # ksym,kvar,var2
    
    for i in kvar:
        f=Function(str(i))(var2)
        df=diff(f)
        kname='d'+alphaname(i)
        nf=symbols(kname)
        kres=kres.subs(df,nf)
    ee=MyEq(kres,kshow=False)
    for i in kvar:
        ee.setdiff(i,i,kshow=False)
        
    return ee.ksym
    
    
    
#####################################
#           list
#####################################  

def solvelist(*args):
    '''
    input: [vector with all eq=0], variables to find ..
    output: MyEq of each variable
    example:
        a+2*b=0 and 3*a-b=0
        ee=[a+2*b,3*a-b]
        then :
        a,b=solvelist(ee,a,b)
        return a,b in MyEq ecuation class
    '''
    vecs=args[1::]
    kres= solve(*args)
    var,value=kunpakDic(kres)
    vres=[]
    for i ,j in zip(var,value):
        ee=MyEq(j,str(i))
        vres.append(ee)
    return vres  

def one_ksym2ksym(ksym):  
    r'''
    return denominator if numerator =1
    input 1/(a+b) return (a+b)
    input 3/(a+b) return  3/(a+b)
    input c/(a+b) return c/(a+b)
    ''' 
    
    if numer(ksym)==1:
        return denom(ksym)
    else:
        return ksym

def num_ksym2ksym(ksym):
    r'''
    return denominator if numerator = numeric
    input 1/(a+b) return (a+b)
    input 3/(a+b) return  (a+b)
    input c/(a+b) return c/(a+b)
    '''
    if Is_Number(numer(ksym)):
        return denom(ksym)
    else:
        return ksym      
        
#####################################
#           real subs
#####################################  
        
def subskwargs(QQ,**kwargs):
    mkey=[]
    vvalue=[]
    
    for key, value in kwargs.items():
        mkey.append(key)
        vvalue.append(value)
    kres=QQ 
    for i,j in zip(mkey,vvalue):
        valor=j
        if type(j)==MyEq:
            valor=j.ksym 
        kres=kres.subs(i,valor)
    return kres 

def real_subs(QQ,**kwargs):
    ''' 
    QQ= symbols function
    ** kwargs c7=6,g4=z..etc..
    RETURN real substitucion when variable have underscore name like 'c_7' 'g_4'    
    '''
    
    vvar=fpoly(QQ,'free')
    mvar=[]
    for i in vvar:
        nname=str(i)
        nname=nname.replace('_','')
        nvar=symbols(nname)
        mvar.append(nvar)
         
    kres=QQ 
    sres=str(kres)
    sres=sres.replace('_','')
    kres=parse_expr(sres)
    mkey=[]
    vvalue=[]
    for key, value in kwargs.items():
        mkey.append(key)
        vvalue.append(value)
        
        kres=kres.subs(parse_expr(key),value)
    for i,j in zip(mvar,vvar):
        kres=kres.subs(i,j)
    return (kres)  

#####################################
#           algebra
#####################################
def get_seudocofactor(e2,e3,var2):
    '''
    tyr to get polynomie complement and multipli complete degree secuence
    e1=x*x-2
    e2=x+1
        if e1=e2*k
        then k=a*x+b
        maybe x*x-2= (x+1)*(a*x+b)
        
    return (a*x+b),[a,b]
    '''
    vecvar=[a,b,c,d]
    vecfvar=[]
    ee2=e2
    if type(e2)==MyEq:
        ee2=e2.ksym
        
    ee3=e3
    if type(e3)==MyEq:
        ee3=e3.ksym    
    qq=degree(ee2,gen=var2)-degree(ee3,gen=var2)
    vres=0
    cc=0
    for i in range(qq+1):
        vres+=vecvar[i]*var2**(qq-cc)
        vecfvar.append(vecvar[i])
        cc+=1
    return (vres,vecfvar)
            
      
def passdoitI(kstr):
    klist=['>','<','=','True','∧',',']
    done=True
    for i in klist:
        if i in kstr:
            done=False
    return done        



def inversefunc(*args):
    '''return inverse function of expr 
       resoect var
        input(function math, var
        inversefunc(x**3+4,x)
        return (x - 4)**(1/3)
    '''   
    if len(args)==0:
        helplib('inversefunc')
        return
    expr=args[0]
    var=args[1]
    yy=symbols('yy')    
     
    expr=subsubs(expr,var,yy) 
     
      
    try:
        kres=ksolve(var-expr,yy)
    except:
        kres=0
        
    return kres

## Wolfram solve Tools

def wolfram_solve_integral(Id,expr,var=x):
     
    sres='Integrate['+str(expr)+','+str(var)+']'
    sres=wolfralan(sres,Id)
    
    p1=sres.find('=')
    sres=sres[p1+1::]
    sres=sres.replace('constant','C1')
    return wolf2sympy(sres)
         
    
 
def wolfralan(sres,Id):
    client = wolframalpha.Client(Id)
    q = sres
    res = client.query(q)
    answer = next(res.results).text
    return answer


def wolf2sympy(sres):
    try:
        kres=parse_mathematica(sres)
        return kres
    except:
        try:
            kres=parse_expr(sres)
            return kres
        except:
            return sres
            

def Is_Diff_expr(expr):
    done=False
    if "'" in str(expr):
        done= True
    return done
    
    
def dothis(expr,func):
    ddone=False
    if Is_Diff_expr(expr):
        expr=amperdiff2normaldiff(expr)
        ddone=True
        
    kres=expr
    try:
        sexpr=func+'('+str(expr)+')'
        kres=parse_expr(sexpr)
    except:
        pass
    if ddone:
        kres=normaldiff2amperdiff(str(kres))
    return kres  

def factorinte(expr,kfac):
    if type(expr)==Add:
        kres=0
        for i in expr.args:
            kres=kres+factorinte(i,kfac)
        return kres    
    done=False
    if signo(kfac)==-1:
        done=True
        kfac=-1*kfac
    kres=expr 
    if type(expr)==Integral:
        p1,p2=expr.args
        p2=p2.args[0]
        if str(p2) not in str(kfac):
            if type(p1)==Mul or type(p1)==Pow:
                mm=fpoly(p1,'list')
                if kfac in mm:
                    if done:
                        kfac=-1*kfac
                    p1=p1/kfac
                    kres=kfac*Integral(p1, p2)
                    
    return kres
    
def datadiff(var,*args):
    '''
       var=main var t
       args=[x,y]= x(t),y(t)
       Vs=['x','y']
       Vn=[x,y]
       Vf=[x(t),y(t)]
       Vd=[dxt,dyt]
       Vd2=[dx2t,dy2t]
       Pd=[x',y']
       Pd2=[x'',y'']
    '''
    
    Vs,Vn,Vf,Vd,Vd2,Pd,Pd2 =[],[],[],[],[],[],[]
    P1,P2=[],[]
    for i in args:
        Vs.append(str(i))
        Vn.append(i)
    for i in Vs:
        Vf.append(Function(i)(var))
    for i in Vf:
        Vd.append(i.diff(var))
        Vd2.append(i.diff(var,var))
    for i in Vn:
        Pd.append(symboldiff(i))
        Pd2.append(symboldiff2(i))
        
    return Vs,Vn,Vf,Vd,Vd2,Pd,Pd2 

def hightfunc(expr,var,*args):
    '''
    expr=pi*h*r*r*t/3
    var=t
    args=h,r 
    return pi*t*h(t)*r(t)**2/3 
    '''
    Vs,Vn,Vf,Vd,Vd2,Pd,Pd2=datadiff(var,*args)
    for i,j in zip(Vn,Vf):
        expr=subsubs(expr,i,j)
    return expr 

def hightdiff(expr,var,*args):
    '''
    expr=pi*h*r*r*t/3
    var=t
    args=h,r 
    return 2*pi*t*h(t)*r(t)*Derivative(r(t), t)/3 + pi*t*r(t)**2*Derivative(h(t), t)/3 
           + pi*h(t)*r(t)**2/3
    '''

    Vs,Vn,Vf,Vd,Vd2,Pd,Pd2=datadiff(var,*args)
    for i,j in zip(Vn,Vf):
        expr=subsubs(expr,i,j)
    return diff(expr,t)    
    
def sympyEq2prime(expr,var,*args):
    '''
    return expr diff sympy expr whith diff prime symbols
    '''

    Vs,Vn,Vf,Vd,Vd2,Pd,Pd2=datadiff(var,*args)
    for i,j in zip(Vd2,Pd2):
        expr=expr.subs(i,j)
    for i,j in zip(Vd,Pd):
        expr=expr.subs(i,j)
    for i,j in zip(Vf,Vn):
        expr=expr.subs(i,j)   
    return expr
  
  
def prime2sympy(expr,var,*args):
    '''
    return expr diff prime symbols 2  diff sympy expr  
    '''
    Vs,Vn,Vf,Vd,Vd2,Pd,Pd2=datadiff(var,*args)
    for i,j in zip(Pd2,Vd2):
        expr=expr.subs(i,j)
    for i,j in zip(Pd,Vd):
        expr=expr.subs(i,j)
    for i,j in zip(Vn,Vf):
        expr=expr.subs(i,j)   
    return expr

def normaldiff(expr,var,*args):
    ff=prime2sympy(expr,var,*args)
    dff=diff(ff,var)
    return dff    