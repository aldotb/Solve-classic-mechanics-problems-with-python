from sympy import *
from latex2sympy2 import latex2sympy
from functools import reduce
import matplotlib.pyplot as plt
import sys
sys.path.append('../Libaldo')
from IPython.display import Math  # ,display 
 
import re
 
# VARIABLE 
t,a,b,c,d,x,y,z,w,i,j,k,alpha=symbols('t a b c d x y z w i j k alpha')
AA=symbols('AA')
########################################## 
#   Variable Crea Transforma
##########################################
def cleanvar(*args):
    mm=[]
    for i in args:
        mm.append(symbols(str(i)))
    if len(mm)==1:
        return mm[0]
    else:    
        return mm 

def clearvar(*args):
    mm=[]
    for i in args:
        mm.append(symbols(str(i)))
    if len(mm)==1:
        return mm[0]
    else:    
        return mm
def diffvar(*args): #crea variables dieferenciables
        mm=[]
        for i in args:
            sres='d_'+alphaname(str(i))
            mm.append(symbols(sres))
        if len(mm)==1:
            return mm[0]
        else:
            return mm
def c2c(kval):
    if Is_Number(kval) and type(kval)==float and (kval-int(kval))==0:
        return int(kval)
    return kval 
        

       
def get_varfunc(ksym): # input symbols 'v(t)'  return symbols 't'  
    sres=str(ksym)     # inputs and output are symbols variables
    p1=sres.find('(')
    p2=sres.find(')')
    sres=sres[p1+1:p2]
    return parse_expr(sres)
    
def antiprimitiva(ksym): # input symbols 'v(t)'  return symbols 'v'
                          # input symbols 'v'  return symbols 'v' 
    
    sres=str(ksym)     # inputs and output are symbols variables
    p1=sres.find('(')
    if p1!=-1 and p1>0:
        sres=sres[0:p1]
        return parse_expr(sres)
    else:
        return ksym
        
        

    
    
def unisymbols(ksym):   
    '''
    unisymbols() :  this function homegenize diferent variables whit 
                    the same symbolic name in omly one in all symbols expresion
    '''                         
    if type(ksym)==list:
        return [unisymbols(i) for i in ksym]
    else:    
    
        try:
            kres=parse_expr(str(ksym))
        except:
            kres=ksym
        return(kres)
    
def sydem(ksym):
    
    '''
    sydem() :   symbol-idem 
                try to return the original function 
                with out auto math transformathis  
    '''
    kres=UnevaluatedExpr(ksym)
    return unisymbols(kres)
    
def sym2func(nval,ksym):
    '''
    sym2func()
    -----------------------    
    input  : nval = symbolic name function
             ksym = posible dependient variable    
    return : nval function
    '''    
    return convFunc(nval,ksym)    
    
def convFunc(nval,ksym): # symbol to Function args=(symbols, var2)
    newF=Function(nval)(ksym) 
    return newF
    
def sym2func(nval,ksym):
    return convFunc(nval,ksym)
    
def primitivename(ksym):
    kres=ksym
    sres=str(ksym)
    if '(' in sres:
        nsres=sres[0:sres.find('(')]
        kres=symbols(nsres)
    return kres

def get_symbols(*args):
    '''
    return vector insode all expr in args
    '''
    mm=[]
    for i in args:
        if type(i)==list:
            for j in i:
                vecm=fpoly(j,'free')
                for k in vecm:
                    if k not in mm:
                        mm.append(k)
        else:
            vecm=fpoly(i,'free')
            for k in vecm:
                if k not in mm:
                    mm.append(k)
    return mm    
def sym2Function(*args):  # symbols vector to Function vector *args=(symbols,symbols,......,var2)
    
    mm=[]
    vt=args[-1]
    for i in args[:-1]:
        mm.append(convFunc(i,vt))
    return mm

def symboldiff(*args):
    if len(args)==1:
        svar=str(args[0])
        dsvar=svar+"'"
        return symbols(dsvar)
    else:
        mm=[]
        for i in args:
            mm.append(symboldiff(i))
        return mm
    
def symboldiff2(*args):
    if len(args)==1:
        svar=str(args[0])
        dsvar=svar+"''"
        return symbols(dsvar)
    else:
        mm=[]
        for i in args:
            mm.append(symboldiff2(i))
        return mm 
#################################
#   create list variables
#################################

def replacestrinvec(vec,val1,val2):
    '''
    vec=['1','2','3','4','5']
    we need change '3' by 'x'
    replacestrinvec(vec,'3','x')
    return ['1','2','x','4','5']
    '''
    kres=[val2 if i==val1 else i for i in vec]
    return kres

def mzero(kval,kop=0):  # return list with [kop,kop,kop,.....] kval times
    mm=[]
    for i in range(kval):
        mm.append(kop)
    return mm    


def subtrac_list(lst1, lst2): # return list2 without similares in list2
    lst3 = [value for value in lst1 if value not in lst2]
    return lst3
    
def listmul(list1):
    kres=1
    for i in list1:
        kres=kres*i
    return kres
    

        
#################################
#             Math
#################################
def nformat(ksym,qn):
    return N(ksym,qn)


def try2float(kk):
    if type(kk)==list:
        mm=[]
        for i in kk:
            mm.append(try2float(i))
        return mm
    else:
        try:
            kk2=float(kk)
            return kk2
        except:
            return kk
    

#  cfrac() *****************
def cfrac(*args,evaluate=False):
    if len(args)==0:
        helplib('cfrac')
        return
    if len(args)==1 and type(args[0])!=str:
        p1=1
        p2=args[0]  
    elif len(args)==1:
        p1=numer(args[0])
        p2=denom(args[0])
    else: 
        p1=args[0]
        p2=args[1]
      
    if not evaluate:
        return UnevaluatedExpr(UnevaluatedExpr(p1)/(p2))
    else:
        try:
            return frs(k1=p1,k2=p2)
        except:
            return p1/p2
    
def frs(k1,k2): # return   S('k1/k2')
    
    if Is_Number(k1) and Is_Number(k2):
        return Rational(k1, k2)
    else:    
        return sydem(k1/k2)

#  sqrs() *****************

def rpow(*args):
    kval=args[0]
    if len(args)==1:
        kres= sqrt(kval)
    else:
        rr=args[1]
        kres= ppow(kval,cfrac(1,rr))
    return kres
  
    
def sqrs(k1,k2='',k3=1): # return pow(k1,S('k2'))
    '''
    sqrs()
    input:  sqrs(a)  return sqrt(a)
    input:  sqrs(a,b)  return sqrt(b) of (a)
    input:  sqrs(a,b,c)  return sqrt(b) of (a)**c

    '''
    k1=1*k1
        
    if k2=='':
        k2=2
    else:
        k2=1*k2
    v1=str(k2)    
    if k3!=1:
        k3=1*k3
        v2='('+str(k3)+')/('+v1+')'
    else:
        v2='1/('+v1+')'
    kres=pow(k1,S(v2))
    
    return(kres)


#  ppow() ***************** 
def ppow(*args):
    if len(args)==0:
        helplib('ppow')
        return
 
    return kpow(*args)
    
def kpow(*args): 
    '''
    kpow(a) return a**2
    kpow(a,b) return a** b 
    kpow(a,b,c) return a** (b/c)
    
    '''
    if len(args)==0:
        helplib('kpow')
        return
    ksym=args[0]
    if len( args)==1:
        return ksym**2
    elif len( args)==2:
        ee=args[1]
        return ksym**ee
    else:
        rr=args[2]
        ee=args[1]
        if rr==2:
            return sqrt(ksym**ee)
        else:
            return ksym**cfrac(ee,rr)
             
     
        
    
    
#  sex2rad() *****************
def sex2rad(k):  # convert segsadesimal to radial
    if k=='':
        helplib('sex2rad')
        return
    k=simplify(pi*k/180)
    return(k)

#  rad2sex() *****************
def rad2sex(k):  # convert radial to sexages.
    if k=='':
        helplib('sex2rad')
        return
    k=simplify(cfrac(180*k,pi))
    return(k)

#  sex2rad_i() *****************
def sex2rad_i(kang,s='r'):
    if s=='s':
        kang=sex2rad(kang)
    return(kang)

#  killabs() ***************** 
def kilabs(ksym):
    if ksym=='':
        helplib('killabs')
        return
    
    '''
        input: abs(x)
        return : x
    '''
    return killAbs(ksym=ksym)
    
def killAbs(ksym):
    kres=ksym 
    try:
        mm=str(kres)
        mm=mm.replace('Abs','')
        return parse_expr(mm)
    except:
        return kres

#  signo() ***************** 
def makepos(kres):
    return kres*signo(kres)
    
    
def signo(ksym):
    if ksym=='':
        helplib('signo')
        return
    '''
        input: x
        return : -1 if x<0
        return : 0 if x=0
        return : 1 if x>0
        
    '''
    kres=1*ksym 
    mm=fpoly(ksym,'free')
    for i in mm:
        kres=kres.subs(i,1)
     
    try:
        if kres<0:
            return -1
        elif kres>0:
            return 1
        else:
            return 1
    except:
        return 1
 
#  unfloatdiv() *****************  
def unfloatdiv(p1,p2): # return p1/p2 whitout avaluate
    P1=str(p1)
    P2=str(p2)
    P3='('+P1+')/('+P2+')'
    return parse_expr(P3,evaluate=False) 
 
#####################
#   get info
#####################

#  insideroot() ***************** 
def insideroot(ksym=''):
    '''
        input: x**(a/b)
        return : x**a

        
    '''
    if ksym=='':
        helplib('insideroot')
        return
    
    return get_inside_root(ksym=1*ksym)
    
def get_inside_root(ksym):
    kres=ksym 
    if Is_Root(kres): 
        mm=fpoly(kres,'list')
        return mm[0]
    return kres 
    
#  insidepow() ***************** 
def insidepow(ksym=''):
    if ksym=='':
        helplib('insidepow')
        return
    '''
        input: x**(a/b)
        return : x**(1/b)
  
    '''
    return get_inside_Pow(ksym=1*ksym)    
def get_inside_Pow(ksym):
    kres=ksym 
    if Is_Pow2(kres): 
        mm=fpoly(kres,'list')
        return mm[0]
    return kres

#  getexponent() ***************** 


#  getbase() ***************** 

     
   
def get_killexpo(ksym):  
    if Is_Mono(ksym):
        mm=fpoly(ksym,'list')
        return mm[0]        


#####################
#   get info Geometry
#####################

#  gethipo() ***************** 
def gethipo(a,b,kope=''):
    if a=='' or b=='':
        helplib('gethipo')
        return
        
    a=a*1
    b=b*1
    return get_hipo(a=a,b=b,kope=kope) # return raiz(a*a+b*b)
    
def get_hipo(a,b,kope=''):
  
    kres=rpow(kpow(a,2)+kpow(b,2),2)
    kres=opemat(kres,kope=kope)
    return kres

#  gethipo() ***************** 
def getcateto(a,b,kope=''):
    if a=='' or b=='':
        helplib('getcateto')
        return
    a=a*1
    b=b*1
    return get_cateto(a=a,b=b,kope=kope) # return raiz(a*a-b*b)
    
def get_cateto(a,b,kope=''): 
    kres=rpow(kpow(a,2)-kpow(b,2),2)
    kres=opemat(kres,kope=kope)
    return kres 

#  componentx()  *****************
def componentx(kr,kang):
    kr=kr*1
    kang=kang*1
    return x_pol(kr=kr,kang=kang)
    
    
def x_pol(kr,kang):
    return kr*cos(kang)
    
#  componenty()  *****************    
def componenty(kr,kang):
    kr=kr*1
    kang=kang*1
    return y_pol(kr=kr,kang=kang)    
def y_pol(kr,kang):
    return kr*sin(kang)    
 
def lenght2point(x1,y1,x2,y2,kope=''):
    x1,y1,x2,y2=1*x1,1*y1,1*x2,1*y2
    
    l1=x2-x1
    l2=y2-y1
    kres=rpow(l1*l1-l2*l2)
    kres=opemat(kres,kope=kope)
    return kres
    
#  algebra1() ***************** 
     
def get_sum2(a,b,kope=''): # return  a*a+b*b
    a=unisymbols(a)
    b=unisymbols(b)
    kres=  kpow(a,2)+kpow(b,2) 
    kres=opemat(kres,kope=kope)
    return kres

def get_dif2(a,b,kope=''):
    a=unisymbols(a)
    b=unisymbols(b)
    kres= kpow(a,2)-kpow(b,2)
    kres=opemat(kres,kope=kope)
    return kres
 

#####################
#   SOLVE
#####################

def ksolve(kres,var,*args,**kwargs):
        try:
            Vk=var.name
        except:
            Vk=str(var)
         
        solu=solve(kres,Vk)
        if type(solu)==list:
            if 'all' in args:
                return solu
            else:
                try:
                    solu=solu[0]
                except:
                    pass
        if len(args)>0:
            if 'positive' in args:
                solu=solu*signo(solu)
        return solu 
        
def csolve(Eq1,ksym,kd='',korden='',kpositive=False,kope='',kdisp=False,unifique=False):
    Eq1=unisymbols(1*Eq1)
    ksym=unisymbols(ksym)

    
    if not(kd==''):
        kdisp=True
    
    if unifique:
        kres=solve(Eq1,sympy_unifiq(Eq1,ksym))
    else:    
        kres=solve(Eq1,ksym)
    qq=len(kres)
    if qq==1:
        kres=kres[0]
    if qq>1:
        kres=list(kres)
    if kpositive:
        kres=[i for i in kres if i>0]
        
        if len(kres)==1:
            kres=kres[0]
    if not(korden==''):
        if qq>1:
            kres=kres[korden]
            
    kres=opemat(kres,kope)
    #qq=len(kres)
    #if qq==1:
        #kres=kres[0]
    if kd!='':
        sR=kd+' ='
        display(Math(sR+latex(kres)))        
    #if kreturn:    
    return(kres)


            
def csolveR(Eq1,ksym,kd='',kope=''): # solve 2do grade Eq return positive part
    nksym=ppow(1*ksym,2)
    kres=csolve(1*Eq1,nksym,kope=kope)
    kres=opemat(kres,kope)
    kres=rpow(kres,2)
    if kd!='':
        sR=kd+' ='
        display(Math(sR+latex(kres)))        
         
    return(kres)



#####################
#   Transformation
#####################
def tsimplify(kres):
    if kres=='':
        helplib('tsimpify')
        return
    return trigsimp(kres)


def subsubs(*args,force=False):
    if len(args)==0:
        helplib('subsubs')
        return
    expr,expr1,expr2=args    
    
    sres=str(expr)
    sres1=str(expr1)
    sres2=str(expr2)
    sres=sres.replace(sres1,sres2)
    if force:
        try:
            kres=parse_expr(sres,evaluate=False)
            return kres
        except:
            return expr
    else:
        try:
            kres=parse_expr(sres)
            return kres
        except:
            return expr    
        
def texpand(kres):
    if kres=='':
        helplib('texpand')
        return
    if type(kres)==Add:
        rtt=0
        mm=fpoly(kres,'list')
        for i in mm:
            rtt=rtt+expand_trig(i)
        return rtt
    else:    
        return expand_trig(kres)
    
    
def opemat(ksym,kope=''):
    if kope=='':
        return ksym
    '''
    opemat(Equation,kope=opt)
    opt 
        'f'= factor(Equation),
        'e'= expand(Equation),
        's'= simplify(Equation),
        't'= trigsimp(Equation)
        'x'= expand trig(Equation)
        'v'= Equation.evalf()
        'a'= apart(Equation)
        'c'= cancel(Equation)
        'E'= Equation.expand(force=True)
        '2' = kill sqrt( x^2 )
    '''
    if '(' in kope:
        newEq=creaQ(ksym,kope)
        return newEq
        
    kres=unisymbols(ksym)
    for i in kope:
    
        if i=='e' :
            try:
                kres=expand(kres)
            except:
                done=True    
        if i=='f' :
            try:
                kres=factor(kres)
            except:
                done=True    
        if i=='t' :
            try:
                kres=trigsimp(kres)
            except:
                done=True
        if i=='x' :
            try:
                kres=expand_trig(kres)
            except:
                done=True    
        if i=='s' :
            try:
                kres=simplify(kres)
            except:
                done=True
        if i=='v' :
             
            try:
                try:
                    kres=float(evalf(kres))
                    float(kres)
                    return(kres)
                except:
                    try:
                        kres=float(evalf(kres))
                    except:
                        kres=evalf(kres)
                    return(kres)
            except:    
                if len(fpoly(kres,'free'))==0:
                    try:
                        kres=kres.evalf()
                        return(kres)
                    except:  
                        try: 
                            if not(type(kres)==float or type(kres)==int):
                                kres=kres.evalf()
                                return(kres)
                        except:
                            kres=kres
                            return(kres)
                    
                        
                        
        if i=='a' :
            kres=apart(kres)
        if i=='c' :
            kres=cancel(kres)
        if i=='r' :
            try:
                kres=kill_root_poly(kres)
            except:    
                if Is_Mono(ksym):
                    kres=kres+1
                    kres=kill_root_poly(kres)
                    kres=kres-1
                else:
                    kres=kill_root_poly(ksym)
            try:
                mm=fpoly(kres,'list')
                skres=str(kres)
                smm=[]
                vres=[]
                for i in mm:
                    vres.append(str(opemat(i,'r')))
                    smm.append(str(i))
                 
                for i,j in zip(smm,vres):
                    skres=skres.replace(i,j)
                kres=parse_expr(kres)    
            except:
                done=true
        if i=='-' :
            kres2=kres
            try:
                kres2=kpow(kres,-1)
                kres2=opemat(kres2,'r')
                kres=kpow(kres2,-1)
            except:    
                done=False
                
        if i=='2' :
            kres=kill_RootPow(kres)
             
        if i=='E' :
            kres=kres.expand(force=True)
         
        if i=='F' :
             
            kres=nsimplify(kres)
        if i=='R' :    
            kres=reduFac2(kres)
        
        if i=='C' :
            mm=fpoly(kres,'free')
            for i in mm:
                try:
                    kres=cut_root2(kres,i)
                except:
                    done=true
        if i=='K'   or i=='k' :
            kres=expand(kres)
            if Is_Mono(kres):
                kres=signed_sqrt(kres*kres)
            if Is_Poly(kres):
                kres=kill_root_poly(kres)
        if i=='N':
            if Is_Number(kres) and not Is_Symbol(kres):
                kres=float(kres)
            kres=N(kres)
            
    return(kres)
    
def reduFac2(ksym):
    kres=ksym 
    try:
        mm=fpoly(factor(kres),'list')
        try:
            mkres=1
            for i in mm:
                if Is_Poly(i):
                    mkres=mkres*i 
            return mkres
        except:
            return kres    
             
    except:
        return kres

def reduMonoFac(ksym):
    kres=ksym
    nkres=1
    dfactor=False
    if Is_Poly(kres):
        kres=factor(kres)
        if Is_Mono(kres):
            dfactor=True
            return denom(kres)
    if Is_Mono(kres):
        mm=fpoly(kres,'list')
         
        done=False
        for i in mm:
            done=False
            eqs=str(i)
            if '+ ' in eqs or '- ' in eqs:
                done=True
            if done:
                nkres=nkres*i
        if dfactor:
            return expand(nkres)
        else:
            return nkres
    else:
        return kres        
        
def opemat_vec(ksym,kope):
    mm=[]
    for i in ksym:
        kres=opemat(i,kope)
        mm.append(kres)
    return mm
    
def opematsec(ksym,kope=''): # equal to opemat but secuential
    kres=unisymbols(ksym)
    for op in kope:
        kres=opemat(kres,kope=op)
    
    return kres    
        
def opematPolySec(ksym,kope=''):
    kres=ksym
    if  Is_Add(kres):
            mm=0
            for i in fpoly(kres,'list'):
                mm+=opemat(i,kope=kope)
            kres=mm    
        
    else:
        kres=opemat(kres,kope)
    
    return kres
    
def opemat_deno(ksym,kope=''):
    knume=numer(ksym)
    kdeno=denom(ksym)
    kdeno=opemat(kdeno,kope=kope)
    kres=knume/kdeno
    return kres

def opemat_nume(ksym,kope=''):
    knume=numer(ksym)
    kdeno=denom(ksym)
    knume=opemat(knume,kope=kope)
    kres=knume/kdeno
    return kres


#####################
#   get info ALgebra
#####################


def fpoly(ksym,kopt='',op2='',op3=''):
    
    '''
    'n': return number of args   
    'list': return list of args
    'get': return args No =op2 in list args
    'get_inv': return inverse of args No=op2 in,list
    'gets': return sum of term(a)+term(b)+... op2= 'ab...'
    'getp': return multi  of term(a)+term(b)+... op2= 'ab...'
    'free': return list of all variable symbols in ksym
    'filt':  
    'unfilt':  
    'frac2sum':  
    'lam':  return lambdify(op2,ksym)
    'simpow': 
    'subs': return(ksym.subs(op2,op3)
    'facmono': factorize ksym whit op2
    'in_allexp':
    'simp_fac':
    'simp_facs':
    '''
    if kopt!='free' and kopt!='n' and kopt!='forze_subs':
        if kopt!='forze_subs':
            karg=ksym.args
            klist=list(karg)
            knum=len(klist)
    kres=ksym
    done=False
    if kopt=='h':
        kres=kres.subs(op2,op3)
        print(op2,op3,kres.subs(op2,op3))
         
        
    if kopt=='n':
        if type(ksym)==int or type(ksym)==float:
            kres=0
            done=True
        elif type(ksym)==Symbol:
            kres=1
            done=True
            
        else:
            kres=len(ksym.args)
            done=True
       
    if kopt=='list':
        if type(ksym)==int or type(ksym)==float:
            kres=[]
            done=True
        elif type(ksym)==Symbol:
             
            mm=[]
            mm.append(ksym)
             
            kres=mm
            done=True
        elif op2=='poly':
            if Is_Mono(ksym):
                    if not Is_Poly(ksym):
                        kres=[ksym]
                        done=True
            
        else:
            karg=ksym.args
            kres=list(karg)
            done=True
    
    if kopt=='get':
        kres=(klist[op2])
        done=True
    if kopt=='get_inv':
        kres=(klist[op2])
        kres=kres**-1
        done=True    
    if kopt=='gets':
        done=True
        for i in op2:
            nsym=klist[int(i)]
            if done:
                mm=nsym
                done=False
            else:
                mm=mm+nsym
        kres=(mm)        
        done=True 
    if kopt=='getp':
        mm=1
        for i in op2:
            nsym=klist[int(i)]
            mm=mm*nsym
        kres=(mm)        
        done=True     
    if kopt=='free':
        try:
            kres=(list(ksym.free_symbols ))
            done=True
        except:
            kres=[]
            done=True
            
         
    if kopt=='frees':
       vsym=fpoly(ksym,'free')
       vsyms=[]
       for i in vsym:
            vsyms.append(str(i))
       kres=vsyms     
    if kopt=='filt' and op2!='':
        kres=0
        for i in fpoly(ksym,'list'):
            if op2 in fpoly(i,'free'):
                kres=kres+i
    
    if kopt=='filtp' and op2!='':
        kres=0
        for i in fpoly(ksym,'list'):
            if op2 in fpoly(i,'free'):
                kres=kres+i    
                
    if kopt=='unfilt' and op2!='':
        kres=0
        for i in fpoly(ksym,'list'):
            if op2 not in fpoly(i,'free'):
                kres=kres+i             
         
         
    if kopt=='frac2sum': 
        kres=klist
         
        kq0=klist[0]
        kq1=kq0.args[0]
        kq2=klist[1]
        kres=(kq2,kq1)
         
        done=True
    if kopt=='rqt':
        kstr='1/'+str(op2)
        kres=Pow(ksym,S(kstr))
    if kopt=='lam':
        kres=lambdify(op2,ksym)
    if kopt=='simpow': 
        kres=pow(fpoly(ksym,'get',0).args[0],fpoly(ksym,'list')[0].args[1]*fpoly(ksym,'get',1))    
    if kopt=='zubs':
        kres=ksym.subs(op2,op3)
    if kopt=='zubsV':
        for i, j in zip(op2,op3):
            kres=kres.subs(i,j)
             
    if kopt=='facmono':
        klist=poly(ksym,'list')
        mm=0
        for i in klist:
            if op2 in fpoly(i,'free'):
                mm=mm+i
        kres=mm    
    if kopt=='in_allexp':
        bres=True
        klist=fpoly(ksym,'list')
        for i in klist:
            if op2 not in fpoly(i,'free'):
                bres=False
        kres=bres
    
    if kopt=='simp_fac':
        newm=0
        oldm=0
        klist=fpoly(ksym,'list')
        kvvar=fpoly(ksym,'free')
        kvar=fpoly(op2,'free')
        kvar=kvar[0]
        for i in klist:
            kres1=fpoly(i/op2,'free')
            if kvar not in kres1:
                newm=newm+i/op2
            else:
                oldm=oldm+i
        if op3==1:
            kres=(newm)
        elif op3==2:
            kres=(oldm)
        elif op3==0:
            kres=(op2*newm)
        else:    
            kres=op2*(newm)+oldm
    
    if kopt=='simp_facs':
         
        kres=ksym
        skres=0
        veck=op2
        for op2 in veck:

            klist=fpoly(kres,'list')
            km1=0
            km2=0
            for i in klist:
                try: 
                    mm= fpoly(i,'list')
                    if len(mm)>0:
                        if op2 in mm:
                            #print(mm)
                            km1=km1+i/op2
                            km2=km2+i
                except:
                    done=False
            kres=kres-km2
            skres=skres+op2*km1
        skres=skres+(kres) 
        kres=skres
        
    if kopt=='list_tree':
        mm=[]
        kres=mm
        for i in fpoly(ksym,'list'):
            mm.append(short_type_name(i))                                  #mm.append(short_type_name(i))
        kres=mm     
        
    if kopt=='get_type':
         kres=short_type_name(ksym)
        
       
    if kopt=='forze_subs':
         
        kexp=parse_expr(str(ksym))
        ksym=parse_expr(str(op2))
        kval=parse_expr(str(op3))
        kres=kexp.subs(ksym,kval)
    
    if kopt=='if_have':
         
        klist=fpoly(ksym,'list')
        ksym=unisymbols(op2)
        kres=0
        for i in klist:
            vsym=fpoly(i,'free')
            if unisymbols(ksym) in vsym:
                kres+=i
    if kopt=='simplypol':
        klist=fpoly(ksym,'list')
        kres=0
        for i in klist:
            kres+=opemat(i,kope=op2)
            
        
    return(kres)
   
# # factorizar Polinomios lineales
def getargs(kres,*args):
    if kres=='':
        helplib('getargs')
        return
    return get_args(kres=kres,*args)
    
def get_args(kres,*args):
         
        for i in args:
            kres=kres.args[i]
        return kres 
def get_factor_with(eqq,kx,kcomple=True):
    mm=fpoly(eqq,'list') # gcf(2*x*x+3*L*x+2*x+4,x) return (3*L+2)
    kres=0
    for i in mm:
        try :
            newm=fpoly(i,'list')
            if kx in newm:
                if kcomple:
                    kres+=i/kx
                else:
                    kres+=i 
        except:
            pass
    return kres 
    
def get_rest(eqq,kx):
    mm=fpoly(eqq,'list') # gcf(2*x*x+3*L*x+2*x+4,x) return (3*L+2)
    kres=0
    for i in mm:
        try :
            vhay=fpoly(i,'free')
            if kx not in vhay:
                kres+=i 
        except:
            pass
    return kres 
def numerexpand(expr):
    p1,p2=fraction(expr)
    p1=expand(p1)
    return cfrac(p1,p2)
def numerfactor(expr):
    p1,p2=fraction(expr)
    p1=factor(p1)
    return cfrac(p1,p2)    
def numersimplify(expr):
    p1,p2=fraction(expr)
    p1=simplify(p1)
    return cfrac(p1,p2)    
def denoexpand(expr):
    p1,p2=fraction(expr)
    p2=expand(p2)
    return cfrac(p1,p2)
def denofactor(expr):
    p1,p2=fraction(expr)
    p2=factor(p2)
    return cfrac(p1,p2)    
def denosimplify(expr):
    p1,p2=fraction(expr)
    p2=simplify(p2)
    return cfrac(p1,p2) 

    
def factorize(expr,fexpr):
    return factors(expr,fexpr) 

def factors(expr,kfac):
    if Is_Add(expr):
        p1=0
        p2=0
        for i in expr.args:
            kres=simplify(i/kfac)
            if denom(kres)==1:
                p2=p2+simplify(i/kfac)
            else:
                p1=p1+i
        if type(p2)==Add:        
            return p1+Mul(kfac,p2,evaluate=False)
        else:
            return expr
    elif Is_Mul(expr):
        if not Is_Poly(expr):
            kres=1
            for i in expr.args:
                if '+' in str(i) or '-' in str(i):
                    kres=kres*factors(i,kfac)
            return kres
        else:
            return expr
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        return factors(p1,kfac)/factors(p2,kfac)
    elif Is_Pow(expr):
        bb=getbase(expr)
        ee=getexpo(expr)
        bb=factors(bb,kfac)
        return bb**ee
    elif Is_Root(expr):
        bb=getbase(expr)
        ee=getroot(expr)
        bb=factors(bb,kfac)
        return rpow(bb,ee)    
    else:
        return expr      
def factoriza(expr,fexpr):
    if Is_Add(expr):
        return factorize(expr,fexpr)
    elif Is_Div(expr):
        knum,kden=fraction(expr)
        p1=factorize(knum,fexpr)
        p2=factorize(kden,fexpr)
        return p1/p2
    elif Is_Pow(expr):
        bb=getbase(expr)
        ee=getexpo(expr)
        bb2=factorize(bb,fexpr)
        return bb2**ee
    elif Is_Mul(expr):
        if  '+' in str(expr):
            kres=1
            for i in expr.args:
                kres=kres*factorize(i,fexpr)
            return kres
        else:
            return expr
    else:
        return expr     
def dfactor(expr,var,var1,op=''):
    '''
    input diff expr( dx/dt then var=t, var1=x
    op= like factorize 
        dD2=diff(x,t,t)
        dD1=diff(x,t)
        dD=x(t)
    
    '''
    Da2,Da1,Da=simplediff(var,var1)
    sexpr=str(expr)
    sD2=str(Da2)
    sD1=str(Da1)
    sD=str(Da)
    sexp=str(expr)
    sexp=sexp.replace(sD2,'dD2')
    sexp=sexp.replace(sD1,'dD1')
    sexp=sexp.replace(sD,'dD')
    dD2,dD1,dD=symbols('dD2 dD1 dD')
    exp2=parse_expr(sexp)
    if op!='':
        mode=parse_expr(op)
        exp2=factorize(exp2,mode)
    else:
        exp2=factor(exp2)
        
    exp2=exp2.subs(dD2,Da2)
    exp2=exp2.subs(dD1,Da1)
    exp2=exp2.subs(dD,Da)
    return(exp2)

def dsimplify(expr,var,var1):
    '''
    input diff expr( dx/dt then var=t, var1=x
    op= like factorize 
        dD2=diff(x,t,t)
        dD1=diff(x,t)
        dD=x(t)
    
    '''
    Da2,Da1,Da=simplediff(var,var1)
    sexpr=str(expr)
    sD2=str(Da2)
    sD1=str(Da1)
    sD=str(Da)
    sexp=str(expr)
    sexp=sexp.replace(sD2,'dD2')
    sexp=sexp.replace(sD1,'dD1')
    sexp=sexp.replace(sD,'dD')
    dD2,dD1,dD=symbols('dD2 dD1 dD')
    exp2=parse_expr(sexp)
     
         
    exp2=simplify(exp2)
     
        
    exp2=exp2.subs(dD2,Da2)
    exp2=exp2.subs(dD1,Da1)
    exp2=exp2.subs(dD,D)
    return(exp2)    
def factorSec(kEq,ksym,kfiltro='.'):
    if type(ksym)==list:
        return MgrupFac(kEq=kEq,ksym=ksym,kfiltro=kfiltro)
    else:
        return grupFac(kEq=kEq,ksym=ksym,kfiltro=kfiltro)

def grupFac(kEq,ksym,kfiltro='.'):
    return My_factor(kEq=kEq,ksym=ksym,kfiltro=kfiltro)


def part(expr,address):
    r"""
    Returns part of an expression
    
    Arguments
    ---------
        expr : sympy expression
        address : (list of integers) indexes of the part
           of the expression tree to be recovered
    Returns
    -------
        requested part of the expression tree
    """
    for num in address:
        expr = expr.args[num]
    return expr

def inpart(expr,repl,address):
    r"""
    Replaces a part of the tree of an expression (and returns
    the copy)
    
    Arguments
    ---------
        expr: (sympy expression) expression to be intervened
        repl: (sympy expression) modified part of the expression
        address: (list of integers) indexes of the part of the
           expression tree to be replaced (see 'part()')
    Returns
    -------
        new expression with the replacement done
    """
    if len(address) == 1:
        largs = list(expr.args)
        largs[address[0]] = repl
        return expr.func(*largs)
    else:
        largs = list(expr.args)
        largs[address[0]] = inpart(expr.args[address[0]],repl,address[1:])
        new = expr.func(*largs)
    return new
    
def cpart(expr,address):
    r"""
    makes easier to visualize walking the tree. It returns a set of two expressions:
    the original expression with the part located by 'address' substituted
    by the symbol 'PIECE' and the part requested.
    """
    PART = Symbol(r'{\color{red}{PART}}')
    return Set(inpart(expr,PART,address),part(expr,address))  

def kreturn(ksym):
    unisymbols(ksym)
    

#####################
#   Differential
#####################


# Fix sympy functions to not have problem whit difernet value to same symbols
def kdiff(ksym,kvar,kope=''): # Force differential

    ksym=unisymbols(ksym)
    kvar=unisymbols(kvar)
    kres=diff(ksym,kvar)
    kres=opemat(kres,kope)
    return kres
    
def kintegrate(ksym,kvar,kope=''):  # force integration
    ksym=unisymbols(ksym)
    kvar=unisymbols(kvar)
    
    kres=integrate(ksym,kvar)
    kres=opemat(kres,kope)
    return kres


#####################
#   Substitution
#####################


def ksubs (ksym,kval,nkval,kope=''):
    ksym=unisymbols(ksym)
    kval=unisymbols(kval)
    nkval=unisymbols(nkval)
    kres=ksym.subs(kval,nkval)
    kres=opemat(kres,kope)
    return kres
    
def psimplify(ksym,op1='',op2='',kope=''):
    ksym=unisymbols(ksym)
    kres=operacion(ksym,op1=op1,op2=op2,kope=kope)
    return kres 



#####################
#   Verifica (Is)
#####################
def typedata(ksym):  # typedata((x*x*5)/8) return
    kres=' '
    if Is_Symbol(ksym):
        kres+='Symbol, '
    if Is_Number(ksym):
        kres+='Number, '
    if Is_Add(ksym):
        kres+='Add, '
    if Is_Mul(ksym):
        kres+='Mul, '
    if Is_Pow(ksym):
        kres+='Pow, '
    if Is_Div(ksym):
        kres+='Div, '
    if Is_Mono(ksym):
        kres+='Mono, '
    if Is_Poly(ksym):
        kres+='Poly, '
    if Is_Pow2(ksym):
        kres+='Pow2, '
    if Is_Root(ksym):
        kres+='Root, '
    if Is_Real(ksym):
        kres+='Real, '
    if Is_Integer(ksym):
        kres+='Integer, '
    if Is_Even(ksym):
        kres+='Even, '
    sE([kres])
    
    
def allType(ksym,kop='list'):
    if kop=='list':
        sE([ksym])
        sE(['Is Polynomie = ',Is_Poly(ksym)]);
        sE(['Is Symbols= ',Is_Symbol(ksym)]);
        sE(['Is Number= ',Is_Number(ksym)]);
        sE(['Is Real= ',Is_Real(ksym)]);
        sE(['Is Integer= ',Is_Integer(ksym)]);
        sE(['Is Even= ',Is_Even(ksym)]);
        sE(['Is Odd= ',Is_Odd(ksym)]);
        sE(['Is Monomie= ',Is_Mono(ksym)]);
        sE(['Is Add= ',Is_Add(ksym)]);
        sE(['Is Mul= ',Is_Mul(ksym)]);
        sE(['Is Pow=',Is_Pow(ksym)]);
        sE(['Is Pow2= ',Is_Pow2(ksym)]);
        sE(['Is Root= ',Is_Root(ksym)])
    else:
        sE([ksym])
        sE(['Is Polynomie = ',Is_Poly(ksym),' Is Symbols= ',Is_Symbol(ksym),' Is Number= ',Is_Number(ksym)]);
        sE(['Is Real= ',Is_Real(ksym),'  Is Integer= ',Is_Integer(ksym),' Is Even= ',Is_Even(ksym),' Is Odd= ',Is_Odd(ksym)]);
        sE(['Is Monomie= ',Is_Mono(ksym),'  Is Add= ',Is_Add(ksym),' Is Mul= ',Is_Mul(ksym)]);
        sE(['Is Pow=',Is_Pow(ksym),'  Is Pow2= ',Is_Pow2(ksym),' Is Root= ',Is_Root(ksym)])
        sE(['---------------------------------------------------------'])

def Is_Poly(ksym):
    done=True
    if Is_Mono(ksym):
        done=False
    
    return done  
    
def Is_Symbol(ksym):
    done=False
    if type(ksym)==Symbol:
        done=True
    if Is_Mono(ksym):
        vmm=fpoly(ksym,'free')
        for i in vmm:
            if type(i)==Symbol:
                done=True

    return done
    
def Is_notSymbol(ksym):
    if Is_Symbol(ksym):
        return False 
    return True    
    
def Is_Number(ksym):
    try:
        if ksym.is_number:
            return True
        else:
            return False
    except:
        return False
def Is_Real(ksym):
    return TrFa(sympify(ksym).is_real) 
   
def Is_Integer(ksym):
    return  TrFa(sympify(ksym).is_integer)
    
def Is_Even(ksym):
    return  (sympify(ksym).is_even ) 
    
def Is_Odd(ksym):
    return  (sympify(ksym).is_odd ) 

def TrFa(kval): # is True False
    if kval==True or kval==False:
        return(kval)
    else:
        return False
    
def Is_Mono(ksym):
    ksym=expand(ksym)
    if type(ksym)==Mul or type(ksym)==Pow or type(ksym)==Symbol:
        return True
    try:
        kn=len(fpoly(ksym,'list0'))
        if kn==1:
            return True
        else:
            return False 
    except:
        return False
                
def Is_Add(ksym):
    kres=ksym 
    if type(kres)==Add:
        return True
    else:
        return False
    
def Is_Mul(ksym):
    kres=ksym 
    if type(kres)==Mul:
        return True
    else:
        return False

def Is_NMul(ksym):
    kres=ksym 
    if Is_Mul(kres):
        if ksym.args[0]==-1:
            return True    
    return False 
        
def Is_PowPow(expr):
    if Is_Pow(expr):
        base=getbase(expr)
        if Is_Pow(base):
            return True
    return False
    
def Is_MulPow(expr):
    if Is_Mul(expr):
        if 'Pow' in srepr(denom(expr)):
            return  True
        if 'Pow' in srepr(numer(expr)):
            return  True
    return False
    
def Is_NPow(ksym):
    if type(ksym)==Mul:
        if ksym.args[0]==-1:
            if Is_Pow(ksym.args[1]):
                return True    
    return False
    
def Is_Pow(ksym):
    kres=ksym
 
    
    if type(kres)==Mul:
        if kres.args[0]==-1:
            if Is_Pow(kres.args[1]):
                return True
    
    elif type(kres)==Pow :
        return True
    else:
        return False
        
def Is_Pow2(ksym):
    try:
        mm=fpoly(ksym,'list')
        p1=mm[0]
        p2=mm[1]
    except:
        return False
 
    if type(ksym)==Pow and p2==2:
        return True
    else:
        return False        
def Is_NumberPow(expr):
    done=False
    if Is_Pow(expr):
        bb=getbase(expr)
        if Is_Number(bb):
            done=True
    return done  

def Is_Log(ksym):
    if type(ksym)==log:
        return True
    else:
        return False

def Is_Root(ksym):
    kres=False
    try:
        if Is_Pow(ksym):
            p1,p2=ksym.args
            if Is_Div(p2):
                kres=True
                return True
    except:
        return False
    return False

    
def Is_MonoMul(ksym):
    if type(ksym)==type(expand(ksym)):
        return True
    else:
        return False

def Is_NMono(expr):
    if not Is_Add(expr):
        try:
            mm=expr.args
            if Is_Number(mm[0]) and signo(mm[0])==-1:
                return True
        except:
            return False
    return False 
    
def Is_MonoPoly(ksym):
    return not Is_MonoMul(ksym)
    
def Is_Integral(ksym):
    mm=str(ksym.args)
    if (mm[-3::])==',))':
        return True
    else:
        return False 
        
def Is_Div(ksym):
    pi,p2=fraction(ksym)
    if p2!=1:
        return True
    else:
        return False
        
def Is_Inversa(expr):
    if Is_Div(expr):
        p1,p2=fraction(expr)
        if p1==1 or p1==-1:
            return True
    return False
def Is_Inverse(expr):
    if Is_Div(expr):
        p1,p2=fraction(expr)
        if p1==1 or p1==-1:
            return True
    return False
    
#######################################
###  Fix Trigomo,etrisc subs functions

def fixremptrig(ksym,alpha,kkval):
    valsin=sin(kkval)
    valcos=cos(kkval)
    msym=[kpow(sin(alpha),x) for x in range(4)]
    vsym=[kpow(valsin,x) for x in range(4)]
    kres=ksym
    for i,j in zip(msym,vsym):
        kres=kres.subs(i,j)
    msym=[kpow(cos(alpha),x) for x in range(4)]
    vsym=[kpow(valcos,x) for x in range(4)] 
    for i,j in zip(msym,vsym):
        kres=kres.subs(i,j)
    return kres  
    valsin=opemat(sin(alpha1),'v')
    msym=[kpow(sin(alpha1),x) for x in range(4)]
    vsym=[opemat(kpow(valsin,x),'v') for x in range(4)]

def simplifac(p1,p2):
    p1=factor(p1)
    p2=factor(p2)
    try:
        if (Is_Mul(p1) and  Is_Symbol(p2)) or (Is_Mul(p2) and Is_Symbol(p2)) or (Is_Mul(p2) and Is_Mul(p2)):
            mm1=fpoly(p1,'list')
            mm2=fpoly(p2,'list')
            P1=p1
            P2=p2
            TF1=False
            TF2=False
            Nu1=0
            Nu2=0
            for i in mm1:
                if Is_Number(i):
                    TF1=True
                    Nu1=i
            for i in mm2:
                if Is_Number(i):
                    TF2=True
                    Nu2=i 
            if TF1 and TF2:
                Nu3=min(Nu1,Nu2)
                p1=p1/Nu3
                p2=p2/Nu3
            P1=p1
            P2=p2
            mm1=fpoly(P2,'free')
            mm2=fpoly(P1,'free')
            for i in mm1:
                if i in mm2:
                    grade1=degree(P1,gen=i)
                    grade2=degree(P2,gen=i)

                    grade3=min(grade1,grade2)
                    nfac=i**grade3
                    p1=p1/nfac
                    p2=p2/nfac
            return p1,p2
        else:
            return p1,p2
    except:
        return p1,p2
def simpliadd(p1,p2):
    P1=p1
    P2=p2

    if Is_Add(p1) or Is_Add(p2):
        mm1=fpoly(p1,'list')
        mm2=fpoly(p2,'list')
        mm3=fpoly(p1,'list')
         
         
        for i in mm3:
            if i in mm1:
                p1=p1-i
                p2=p2-i
         
                 
    return  p2,p1
def simplifyr(ksym):
    kres=ksym
    kres=simplify(kres*kres)
    kres=rpow(kres,2)
    return kres




    
def lexp_simplify(ksym):
    sksym=str(ksym)
    if 'log(exp(' in sksym:
        kint=insidepar(insidepar(kk,'log(exp('),'exp(')
        skint=str(kint)
        tkill='log(exp('+skint+'))'
        sksym=sksym.replace(tkill,skint)
        return parse_expr(sksym)
    else:
        return ksym
        
      
        
        
def fix_sqrt2pow(ksym):
    try:
        kres=ksym
        mm=str(ksym)
        mm1= mm.replace('**2','**1')
        mm2=mm1.replace('sqrt(','(')
        mm3=parse_expr(mm2)
        return mm3
    except:
        return ksym
    

def tintegral_def(keq,alpha,a1,a2,kope=''):
    kfun=kintegrate(keq,alpha)
    val1=fixremptrig(kfun,alpha,a1)
    val2=fixremptrig(kfun,alpha,a2)
    kres=val2-val1
    kres=opemat(kres,kope=kope)
    return kres

 
def change_diff(ksym,y,x,newQ=''): # ksym =Integral func, y old v,x =mew v, new Func
    if Is_Integral(ksym):
         
        Isol=ksym.doit()
        if newQ!='':
            Isol=Isol.subs(y,newQ)
        else:
            Isol.subs(y,x)
        return Integral(Isol,x)
    else:
        return ksym       
            
def cut_fac(ksym,kval):
    if type(ksym)==Mul:
        return(simplify(unisymbols(ksym/kval)))
    elif type(ksym)==Add:
        mlist=fpoly(ksym,'list')
        mm=0
        for i in mlist:
            mm+=cut_fac(i,kval)
        return mm   
         
def cut_root2(ksym,kval):
    kk2='sqrt('+str(kval)+'**2)'
    if ksym!=0:
        try:
            if type(ksym)==Mul:
                kk= str(ksym)
                 
                kk3=kk.replace(kk2,str(kval))
                ksol=parse_expr(kk3)
                return ksol
            if type(ksym)==Add:
                nksym=0
                mm=fpoly(ksym,'list')
                for i in mm:
                    nksym+=cut_root2(i,kval)
                    
                return nksym    
        except:
            return ksym
    else:
        return ksym



        
class MyTriang:
    def __init__(self, hipo='',cat1='',cat2='',kope=''):

        self.khipo=hipo
        self.kcat1=cat1
        self.kcat2=cat2
        
        if hipo=='':
            kres=get_hipo(cat1,cat2)
            self.khipo=opemat(kres,kope=kope)
        if cat2=='':
            kres=get_cateto(hipo,cat1)
            self.kcat2=opemat(kres,kope=kope)
        if cat1=='':
            kres=get_cateto(hipo,cat2)
            self.kcat1=opemat(kres,kope=kope)   
    
    def sin(self,kope=''):
        hipo=self.khipo
        cat1=self.kcat1
        kres=cat1/hipo
        kres=opemat(kres,kope=kope)
        
        return kres
    
    def cos(self,kope=''):
        hipo=self.khipo
        cat2=self.kcat2
        kres=cat2/hipo
        kres=opemat(kres,kope=kope)
        
        return kres
    
    def tan(self,kope=''):
        cat1=self.kcat1
        cat2=self.kcat2
        kres=cat1/cat2
        kres=opemat(kres,kope=kope)
        
        return kres
    
    def hipo(self,kope=''):
        kres=self.khipo
        
        return kres 
        
    def cat1(self,kope=''):
        kres=self.kcat1
        
        return kres
        
    def cat2(self,kope=''):
        kres=self.kcat2
        
        return kres             
    
    def s(self):
        sE(['sin()=',self.sin(),'  cos()=',self.cos(),'tan()=',self.tan()])
        
        
def sqrt2fracpow(expr):
    return(signed_sqrt(expr))
    
def signed_sqrt(expr):  # This function from WenyinWei founded in GitHub
    """Signed sqrt operator
    Args:
        expr (sympy.expr): sympy expression
    Returns:
        sympy.expr: A simplified expression
    """


    expr = expr.factor()
    # recurse the function on each arg if the top function is a multiplication 
    # e.g. signed_sqrt( 4 * b^2 ) == 2 * b 
    if expr.func == Mul: 
        args_signed_sqrt = [signed_sqrt(arg) for arg in expr.args]
        return reduce(Mul, args_signed_sqrt)
    elif expr.func == Pow:
        base, exponent = expr.args
        if exponent.is_even:
            return base**(exponent/2)
    return sqrt(expr) 

 
def fixrootpow(expr,k=''):
    if k=='':
        expr=fixrootpow(expr,2)
        expr=fixrootpow(expr,3)
        return expr
    else:    
        sres=str(expr)
        snum=str(k)
        srot='**'+snum+')**(1/'+snum+')'
        if srot in sres:
            sres=sres.replace(srot,')')
            nexpr=parse_expr(sres)
            return nexpr
        else:
            return expr 
    
def KrP(ksym,kope=''):
    return kill_root_mono(ksym,kope=kope)
    
def kill_root_mono(ksym,kope=''): # kill root(pow(ksym))   if  Is_Mono(ksym) = True
    
    if Is_Root(ksym):
        kres= signed_sqrt(ksym*ksym)
    else:
        kres= ksym

    return kres

def kill_root_poly(ksym,kope=''):  # kill root(pow(ksym1)) + root(pow(ksym1))   if  Is_Poly(ksym) = True
    if Is_Poly(ksym):                #  ksym=ksym1+ksym2+ ....
        mm=0
        vksym=fpoly(ksym,'list')
        for i in vksym:
            mm+=kill_root_poly(i,kope='')
        mm=opemat(mm,kope=kope)
        return mm
    elif Is_Mul(ksym):                #  ksym=ksym1+ksym2+ ....
        kt=1
        vksym=fpoly(ksym,'list')
        for i in vksym:
            kt*=kill_root_poly(i)
         
        return kt   
    else:
        kres=kill_root_mono(ksym,kope='')
        return kres   
        
        
def sin2cos(ksym,angu,korden=2,kope=''):
    e1=unisymbols(ksym)
    e1=e1.subs(unisymbols(kpow(sin(angu),4)),(unisymbols(kpow(sin(angu),2)*kpow(sin(angu),2))))
    e1=e1.subs(unisymbols(kpow(sin(angu),3)),(unisymbols(kpow(sin(angu),2)*sin(angu))))
    e1=e1.subs(unisymbols(kpow(sin(angu),2)),(unisymbols(1-kpow(cos(angu),2))))
    if korden==1:
        e1=e1.subs(unisymbols(sin(angu)),(unisymbolsrpow(1-kpow(cos(angu),2))) )
    
    kres=e1 
    kres=opemat(kres,kope=kope)
    return kres
    
def cos2sin(ksym,angu,korden=2,kope=''):
    e1=unisymbols(ksym)
    e1=e1.subs(unisymbols(kpow(cos(angu),4)),unisymbols((kpow(cos(angu),2)*kpow(cos(angu),2))))
    e1=e1.subs(unisymbols(kpow(cos(angu),3)),unisymbols((kpow(cos(angu),2)*cos(angu))))
    e1=e1.subs(unisymbols(kpow(cos(angu),2)),unisymbols((1-kpow(sin(angu),2))))
    if korden==1:
        e1=e1.sub(unisymbols(cos(angu)),unisymbols(rpow(1-kpow(sin(angu),2))))
    
    kres=e1 
    kres=opemat(kres,kope=kope)
    return kres    
    
def MaT(x,y=''):
    if y=='':
        y=x[1]
        x=x[0]
    return(Matrix([[x,y]]))
    
def moduloMat(kMat):
    xx1=kMat[0]
    yy1=kMat[1]

    return get_hipo(xx1,yy1)    
        

#  algebra
def sortPoly(ksym,kvar,kindi):
    kres=ksym
    for i in range(kindi):
        kres=factorSec(kres,kpow(kvar,i+1))
    
    return kres    
    

    #  Diccionario
def unpack(mm):
    return kunpakDic(mm=mm)
    
def kunpakDic(mm):
     
    kkey=list(mm.keys())
    kvalu=list(mm.values())

    return( kkey,kvalu)  

#######################
##  used by MyEq
#######################

def multiSet(ksym, kval, vecEq=[]):
    if type(ksym) == list:
        for i, j in zip(ksym, kval):
            for kQ in vecEq:
                kQ.setValue(i, j, kshow=False)
    else:
        for kQ in vecEq:
            kQ.setValue(ksym, kval, kshow=False)

    for i in vecEq:
        i.s()
        
def multiRedFac(kvec=[]):
    for i in kvec:
        i.reduFac()
        i.s()
        
# convinaciones de soluciones
def solveFrom(eeV,ssV,kope=''):
    vec1=[x() for x in eeV]
    kres= solve(vec1,ssV)
    kres1=[opemat(x,kope=kope) for x in kres]
    kres2=kres1[0]
    mainres=[]
    for i ,j in zip(kres2,ssV):
        mainres.append(i)
        ee=MyEq(i,j.name,kshow=False,kope=kope)
        ee.s()
    return mainres
#def solve2sysEq(kvar=[], keq=[], kname=[]):

def killPwise(sksym): # Kill otherwise answer in simple str
    kres=sksym
    if 'Piecewise' in kres:
        x1=sksym.find('Piecewise(')
        x2=x1+len('Piecewise(')
        x3=sksym.find(', Ne(')
        x4=sksym.find('True)')
        x5=x4+len('True)')
        kres=sksym[0:x1]+sksym[x2:x3]+sksym[x5::]
    return kres
    
def fix_otherwise(ksym,kop=''): # Kill otherwise answer in answer
    kres2=ksym
    if kop=='odb':
        ksym=ksym.rhs
        
    mm=fpoly(ksym,'list')
    kres=''
    done=True
    kk='+'
    for i in mm:
        ss=str(i)
        if done:
            kres+=killPwise(ss)
            done=False
        else:
            kres+='+'+killPwise(ss)
    try: 
        return unisymbols(parse_expr (kres))
    except:
        return kres2
    
       
#######################
##  used by MyInteger
#######################

    
def Ope2Inte(e1,e2,kope='Add'):
    ope1=e1.kinte
    ope2=e2.kinte
    if kope=='Mul':
        return ope1*ope2
    elif kope=='Div':
        return ope1/ope2

    else:
        return ope1+ope2
        
def opeInteSolu(val1,val2,kope='Add'):
    if kope=='Mul':
        return val1*val2
    elif kope=='Div':
        return val1/val2 

    else:
        return val1+val2         

def miniopI(val1,val2,ktype='Add'):
    if ktype=='Mul':
        return val1.kinte*val2.kinte
    elif ktype=='Div':
        return val1.kinte/val2.kinte
    else:
        return val1.kinte+val2.kinte

        

##############  LAtex
def lxprint(*args):
    vec=''
    for i in args:
        if type(i)==str:
            vec+= i+'\;'
        else:
            vec+= latex(i)+'\;'
    display(Math(vec))    
  
def symb_diff(*args):
    kres=''
    for i in args:
        kres=kres+' d'+alphaname(i)
    return kres
    
def diff_name(ksym):
    kres='d_'+ alphaname(ksym)
    return kres
        
def diff_name_prima(ksym):
    kres=alphaname(ksym)+"'"
    return kres 
 
def diffname(k1,k2):
    if type(k1)!=str:
        k1=k1.name
    k1=alphaname(primitivename(k1))
    k2=alphaname(primitivename(k2))
    xx='d_'+k1
    tt='d_'+k2
    xkname='\\frac{'+xx+'}{'+tt+'}'
    return xkname 
    
def difffuncname(kfunc,ksym):
    kres='d'+alphaname(kfunc)+'('+alphaname(ksym)+')'
    return kres
    
def funcname(kfunc,ksym):
    kres=alphaname(kfunc)+'('+alphaname(ksym)+')'
    return kres

# def diffname(ksym):
    # kres='d'+ alphaname(ksym)    
    # return kres
def clean_underline(ksym):
        sres=str(ksym)
        sres=sres.replace('_','')
        return sres

def alphasubname(ksym,op=1):
    aaname=alphaname(ksym)
    aaname=aaname+'_'+str(op)
    return aaname



def eQrec(x1=0,y1=0,x2=0,y2=0,var2=''):
    mm=frs((y2-y1),(x2-x1))
    bb=y2-x2*mm
    kres=opemat(var2*mm+bb,'s')
    return kres

def diffvariable(k1,k2):
    return get_diff_name(k1,k2)
 
def get_diff_name(k1,k2):
    k1=alphaname(k1)
    k2=alphaname(k2)
    xx='d_'+k1
    tt='d_'+k2
    xt='\\frac{'+xx+'}{'+tt+'}'
    dxt=symbols(xt)
    return dxt

def difvar(*args): #crea variables dieferenciables

        mm=[]
        for i in args:
            sres='d'+diffsymbol(i)
        mm.append(symbols(sres))
        return mm

 
        
def Cg2func(f1,f2,x,x1,x2):
    Area=integrate(f1,(x,x1,x2))-integrate(f2,(x,x1,x2))
    X=integrate((f1-f2)*x,(x,x1,x2))
    X=X/Area
     
    Y=integrate((f1-f2)*(f1+f2)/2,(x,x1,x2))
    Y=Y/Area
    return X,Y 
    
def findSubFunc(ksym,sval,inside=''):
         
        kini=0
        kini2=0
        sroot=[]
        done=0
        while kini<len(str(ksym)) and kini2!=-1:
            kini2,sword=in_ope_string(ksym,sval,kini)
            if kini2!=-1:
                if inside!='':
                    if inside in sword:
                        sroot.append(sword)
                else:
                    sroot.append(sword)
                        
                 
            kini=kini2+len(sval)
        return sroot 

def balance(ssym):  # usada x findSubFunc
    cc=0
    for i in ssym:
        if i=='(':
            cc+=1
        if i==')':
            cc-=1
    return cc 
def mirror_parse_expr(sexpr):
    try:
        return parse_expr(sexpr,evaluate=False)  
    except:
        return parse_expr(sexpr)
def in_ope_string(ksym,sval,kini=0):  # usada x findSubFunc
    ssym=str(ksym)
    qq=len(ssym)
    cc=ssym.find(sval,kini,qq)
    inip=cc+len(sval)
    
    for i in range(inip+1,len(ssym)):
        sward=ssym[inip:i]
        veri=balance(sward)
        if veri==0:
            return cc,sval+sward 
            
def get_midle_str(ssym,p1,p2):  #get_midle_str('123456789','123','89') return '4567' 
     
    qq=len(ssym)
    q1=len(p1)
    q2=len(p2)
     
    ssym=ssym.replace(p1,'')
    ssym=ssym.replace(p2,'')
    return ssym
            
def get_vecposstrfind(kstr,sval):
    indices_object = re.finditer(pattern=sval, string=kstr)
    indices = [index.start() for index in indices_object]
    return indices 

# DIFERENCIAL 

def derivarespect(expr,var1,Dd):
    try:
        p1=Math(diff(expr,var1)*Dd)
    except:
        p1=diff(expr,var1)*Dd
        
    return p1
##  CREA PRIMITIVA    crea_primitiva(e1),crea_primitiva(ksym,[alpha,r],crea_primitiva(ksym,t,[alpha.r])



def func2primi(*args): 
    if len(args)==1:
        ee=args[0]
        ksym=ee.ksym
        varf=ee.varf
        var2=ee.var2
    elif len(args)==2 :
        ksym=args[0]
        varf=args[1]
        var2=t
    else:
        ksym=args[0]
        var2=args[1]
        varf=args[2]     
  
    if varf==[]:
        return ksym
    else:


        # variable
        nname=[]
        for i in varf:
            nname.append(str(i))
        nF=[]
        for i in nname:
            nF.append(Function(i)(var2))
           
        for i,j in zip(varf,nF):
            ksym=ksym.subs(i,j) 
         
        return ksym   


#  return primitiva to normal  short_primitiva(e1),short_primitiva(ksym,[alpha,r],short_primitiva(ksym,t,[alpha.r])
def primi2func(*args):   
    if len(args)==1:
        ee=args[0]
        ksym=ee.primitiva
        varf=ee.varf
        var2=ee.var2
    elif len(args)==2 :
        ksym=args[0]
        varf=args[1]
        var2=t
    else:
        ksym=args[0]
        var2=args[1]
        varf=args[2]
    sres=str(ksym)
    svar=[str(x) for x in varf]
    svar2=str(var2)
    
    oldvar=[x+'('+svar2+')' for x in svar]
    
    for i,j in zip(oldvar,svar):
        sres=sres.replace(i,j)
    return parse_expr(sres) 
    
def diffrespect(expr,V,t):
    
    F=Function(str(V))(t)
    sF=str(F)
    sres=str(expr) 
    sres=sres.replace(str(V),sF)
    kres=parse_expr(sres)
    kres=diff(kres,t)
    return kres
    
    
def Diff2diff(kres,kvar,var2): # ksym,kvar,var2
    
    for i in kvar:
        f=Function(str(i))(var2)
        df=diff(f)
        kname='d'+alphaname(i)
        nf=symbols(kname)
        kres=kres.subs(df,nf)
    return kres  

#  str Functions



def flat_diff(ksym,fd):
     
    
    if 'Derivative' in str(ksym):
        sres=between_par(str(ksym))
        s1='d'+sres[0]
        s2='d'+between_par(sres)

        d1=symbols(s1)
        d2=symbols(s2)
        return(ksym.subs(fd,d1/d2))
    else:
        return ksym
    
def reduFac(ksym):  # retun Eq=0,= a*b+c*b.. = b(a+c)..=0  then return (a+c)

        kres = ksym
        kres = factor(kres)
        kres2 = 1
        if Is_Mono(kres):
            kres = numer(kres)

            mm = fpoly(kres, 'list')
            for i in mm:
                if Is_Poly(i):
                    kres2 = kres2 * i

        if kres2 != 0 and kres2 != 1:
            return kres2
        else:
            return ksym
   
   
# algoritmos 2022

def mono_sin_numeros(ksym): # retorna un monomios sin coeficientes numericos
    if Is_Mono(ksym):
        if Is_Mul(ksym):
            kres=ksym
            vsym=fpoly(kres,'free')
            for i in vsym:
                kres=kres.subs(i,1)
        ksym=simplify(ksym/kres)
    return ksym   
    
def killrpow(expr):
    if Is_Add(expr):
        kres=0
        for i in fpoly(expr,'list'):
            kres=kres+killrpow(i)
        return kres
    elif denom(expr)!=1:
        return cfrac(killrpow(numer(expr)),kilrpow(denom((expr))))
    elif Is_Mul(expr):
        kres=1
        for i in fpoly(expr,'list'):
            kres=kres*killrpow(i)
        return kres
    elif Is_Pow(expr):
        try:
            kres=powdenest(sqrt(expr, force=True))
            return kres2

        except:    
            bb=getbase(expr)
            ee=getexpo(expr)
            return kpow(kilrpow(bb),killrpow(ee))
    else:
        return powdenest(sqrt(expr, force=True))
        
        
def kill_RootPow(ksym): # kill rootPow in Polynomie one pass
    return killrpow(ksym)   

def simplifyexpand(kres):
    firsti=True
    if type(kres)==Add:
        val=0
        sres=''
        mm=fpoly(kres,'list')
        
        for i in mm:
            val2=simplifyexpand(i)
            sval=str(val2)
            if firsti:
                sres=sval
                firsti=False
            elif sval[0]=='-':
                sres=sres+sval
            else:
                sres=sres+'+'+sval
        kres=parse_expr(sres,evaluate=False)       
        return kres
    else:
        return apart(expand(simplify(kres)))
        
def base_exponent(ksym):
    if Is_Pow(ksym):
        mm=fpoly(ksym,'list')
        return mm[0],mm[1]
    else:
        return ksym,1

def reducecero(ksym):
    kres= reducecero2(ksym)
    kres2=reducecero2(kres)
    return kres2
    
    
def reducecero2(ksym): 
    if type(ksym)==Add:
        done=false
        ksym=factor(ksym)
        if type(ksym)==Mul:
            ksym= reducecero2(ksym)
            return ksym
        else:
            return ksym
            
    elif Is_Div(ksym):
        ksym=numer(ksym)
        ksym=reducecero2(ksym)
        return ksym
    elif Is_Mul(ksym):
        mm=ksym.args
        kres1=1
        kres2=1
        for i in mm:
            if Is_Add(i):
                kres1=kres1*i
            else:
                kres2=kres2*i
        if kres1!=1:
            return kres1
        else:
            return kres2
            
    else:
        return ksym

def tfunc53(ksym,angle):
    kres=ksym
    kres=kres.subs(sin(alpha),cfrac(4,5))
    kres=kres.subs(cos(alpha),cfrac(3,5))
    kres=kres.subs(tan(alpha),cfrac(4,3))
    return kres

def tfunc37(ksym,angle):
    kres=ksym
    kres=kres.subs(sin(alpha),cfrac(3,5))
    kres=kres.subs(cos(alpha),cfrac(4,5))
    kres=kres.subs(tan(alpha),cfrac(3,4))
    return kres    
    
def tfunc16(ksym,angle):
    kres=ksym
    kres=kres.subs(sin(alpha),cfrac(7,25))
    kres=kres.subs(cos(alpha),cfrac(24,25))
    kres=kres.subs(tan(alpha),cfrac(7,24))
    return kres

def tfunc74(ksym,angle):
    kres=ksym
    kres=kres.subs(sin(alpha),cfrac(24,25))
    kres=kres.subs(cos(alpha),cfrac(7,25))
    kres=kres.subs(tan(alpha),cfrac(24,7))
    return kres    
    
def cos_if_tan(cat1,cat2=1):
    hipo=gethipo(cat1,cat2)
    return cfrac(cat2,hipo)
def cos_if_sin(cat1,hipo=1):
    cat2=getcateto(hipo,cat1)
    return cfrac(cat2,hipo)


def sin_if_tan(cat1,cat2=1):
    hipo=gethipo(cat1,cat2)
    return cfrac(cat1,hipo)
def sin_if_cos(cat1,hipo=1):
    cat2=getcateto(hipo,cat1)
    return cfrac(cat2,hipo)

def tan_if_sin(cat1,cat2=1):
    return cfrac(cat1,cat2)
def tan_if_cos(cat1,cat2=1):
    return cfrac(cat2,cat1)



def diffS(p1,p2,op=''):
    # input variables example V,t*gg
    # return symbols with name dV/dt
    sp1='d'+alphaname(p1)
    sp2='d'+alphaname(p2)
    name='\\frac{'+sp2+'}{'+sp1+'}'
    dp1=symbols(name)
    if op==2:
        dp2=diffS2(p1,p2)
        return dp1,dp2
    else:
        return dp1
def diffS2(p1,p2):
    # input variables example V,t*gg
    # return symbols with name dV2/dt2
    sp1='d'+alphaname(p1)+'_2'
    sp2='d'+alphaname(p2)+'_2'
    name='\\frac{'+sp2+'}{'+sp1+'}'
    return symbols(name)

def viewnicediff(expr,var,var1,var2=''):
      
    nd1x,nd2x=diffS(var,var1,2)
    od2x,od1x,odx=simplediff(var,var1)
    sod2x,sod1x,sodx=str(od2x),str(od1x),str(odx)
    sexpr=str(expr)
     
    sexpr=sexpr.replace(sod2x,'D2x')
    sexpr=sexpr.replace(sod1x,'D1x')
    sexpr=sexpr.replace(sodx,'Dx')
    if var2!='':
         
        nd2y,nd1y=diffS(var,var2,2)
        od2y,od1y,ody=simplediff(var,var2)
        sod2y,sod1y,soy=str(od2y),str(od1y),str(ody)
        sexpr=sexpr.replace(sod2y,'D2y')
        sexpr=sexpr.replace(sod1y,'D1y')
        sexpr=sexpr.replace(sody,'Dy')
    D2x,D1x,Dx,D2y,D1y,Dy=symbols('D2x D1x Dx D2y D1y Dy')
    kres=parse_expr(sexpr)
    kres=kres.subs(D2x,nd2x)
    kres=kres.subs(D1x,nd1x)
    kres=kres.subs(Dx,var1)
    if var2!='':
        
        kres=kres.subs(D2y,nd2y)
        kres=kres.subs(D1y,nd1y)
        kres=kres.subs(Dy,var2)
    
    return  kres

    
    

def diffuntion(expr,var,var1,var2=''):
    fname1=var1
    f1=Function(fname1)(var)
    if var2=='':
        expr=expr.subs(var1,f1)
        kres=expr.diff(var)
    else:
        fname2=var2
        f2=Function(fname2)(var)
        expr=expr.subs(var1,f1)
        expr=expr.subs(var2,f2)
        kres=expr.diff(var)
    return kres 



def simplediff(t,x):
    F=Function(str(x))(t)
    return  F.diff(t,t),F.diff(t),F  



def pack(**kwargs):
    return kwargs

def killpar(expr):
    sexpr=str(expr)
    if sexpr[0]=='[' and sexpr[-1]==']':
        sexpr=sexpr[1:-1]
    return parse_expr(sexpr)
    
def ganswer(expr,*args):
    if type(expr)==dict:
        exprt = zip(expr.keys(), expr.values())
        sname,valor=unpack(expr)
    if type(expr)==list:
        sname=[]
        valor=[]
        for i in expr:
            sname.append(i.lhs)
            valor.append(i.rhs)
         
         
   
    if 'value' in args and not 'varname' in args:
        return valor
    if not 'value' in args and   'varname' in args:
        return sname
    
    if 'varname' and   'varname' in args:
        return sname,valor
 
#  diff procediments 
def traducediff(*args):
    expr=args[0]
    var=args[1]
    
    vecv=args[2:len(args)]
    for i in vecv:
        f=Function(str(i))(var)
        sD2=str(i)+"''"
        eD2=f.diff(var,var)
        expr=expr.replace(sD2,str(eD2))
        sD1=str(i)+"'"
        eD1=f.diff(var)
        expr=expr.replace(sD1,str(eD1))
    return  parse_expr(expr) 
    
def diff2mark(expr,var,var1):
    '''
        expr= Eq() equation equallity
        vx= independ variable x
        vy= dependent var y(x)

        diff2mark(diff1=diff2,y,x)

    '''
       
    vx=var
    vy=var1
    dy=symboldiff(y)
    dy2=symboldiff2(y)
    Y=Function(str(vy))(vx)

    expr=expr.subs(Y.diff(vx,vx),dy2)
    expr=expr.subs(Y.diff(),dy)
    expr=expr.subs(Y,vy)
    return expr    
    
def diff2marko(expr,var,var1):
    o=symbols(str(var1))
    vx=var
    vy=var1
    dy=symboldiff(o)
    dy2=symboldiff2(o)
    O=Function(str(vy))(vx)

    expr=expr.subs(O.diff(vx,vx),dy2)
    expr=expr.subs(O.diff(),dy)
    expr=expr.subs(O,vy)
    return expr    
    
    
class Myhistory:
    def __init__(self,ginput=[]):
        self.ginput=ginput
           
           
    def  Add(self,gadd):
        self.ginput.append(gadd)
        
        
def Div(self,dval):
    if Is_Div(self):
        p1,p2=fraction(self)
        P1=p1
        P2=p2*dval
        return cfrac(P1,P2)
    else:
        return cfrac(self,dval)       
        
 

    
    
def getbase(expr): # return base in expr
    '''
    expr = x**y,x*y
    return x ,x*y
    '''
    
    if Is_Pow(expr):
        ee=expr.args[0]
        return ee
    else:
        return expr  
        
def getexpo(expr,op=''): #   return exponente in expr
    if Is_Pow(expr):
        mm=expr.args
        return mm[1]
        
    else:
        return 1   
 
def alphaname(ksym,op1=''):
    kk=str(ksym)
    if type(ksym)!=str:
        kk=str(ksym)
        
    if kk=='alpha':
        return 'α'
    if kk=='alpha1':
        return 'α1'
    if kk=='alpha2':
        return 'α2'
    if kk=='alpha3':
        return 'α3' 
        
    if kk=='betha':
        return 'ß'
    if kk=='betha1':
        return 'ß1'
    if kk=='betha2':
        return 'ß2'    
    if kk=='betha3':
        return 'ß3'

        

    if kk=='tetha':
        return 'θ' 
    if kk=='tetha1':
        return 'θ1'
    if kk=='tetha2':
        return 'θ2' 
    if kk=='tetha3':
        return 'θ3'
            
    if kk=='calpha':
        return 'cos('+'α'+')'        
    if kk=='salpha':
        return 'sin('+'α'+')'    
    if kk=='talpha':
        return 'tan('+'α'+')'
        
    if kk=='cbetha':
        return 'cos('+'ß'+')'
    if kk=='sbetha':
        return 'sin('+'ß'+')'    
    if kk=='tbetha':
        return 'tan('+'ß'+')'
        
    if kk=='ctetha':
        return 'cos('+'θ'+')'
    if kk=='stetha':
        return 'sin('+'θ'+')'    
    if kk=='ttetha':
        return 'tan('+'θ'+')'
        
    if kk=='dalpha':
        return 'dα'
    if kk=='dbetha':
        return 'dß'
    if kk=='dtetha':
        return 'dθ' 
        
        
    if kk=='mu':
        return '𝜇'
    if kk=='mu1':
        return '𝜇1'
    if kk=='mu2':
        return '𝜇2'
    if kk=='mu3':
        return '𝜇3'    
    else:
        return kk  
vecreatr=["<class 'sympy.core.symbol.Symbol'>","<class 'int'>","<class 'float'>","<class 'sympy.core.numbers.Pi'>","<class 'sympy.core.numbers.Rational'>"]
def arglist(expr,deep=3):
     
    infoexpr=[]
    infopos=[]
    cc=''
    A,B,C,D=ruta(expr,infoexpr,infopos,cc)
    BC=[[i,j] for i,j in zip(B,C) if  not Is_NMono(i)]
    B1,C1=[],[]
    for i in BC:
        B1.append(i[0])
        C1.append(i[1])
    BC2=[[i,j] for i,j in zip(B1,C1) if  len(j)<deep+1]
    B2,C2=[],[]
    for i in BC2:
        B2.append(i[0])
        C2.append(i[1])
    BC3=[[i,j] for i,j in zip(B2,C2) if not Is_Inverse(i)]
    B3,C3=[],[]
    for i in BC3:
        B3.append(i[0])
        C3.append(i[1])
    return B3,C3

def showarglist(expr,*args,deep=3,format=None,side=None):
    '''
    expr= Math epxr like (2*a**3-b**3-c**3)/(a*b+b*c+a*c)
    deep=2 -- args(0,1,2,3,... deep-1)
    format=list..return list    
           defaul vector
    '''
    
    vecv,vecp=arglist(expr=expr,deep=deep)
    if len(args)==0:
        svecargs=[]
        for i in vecp:
            if side=='L':
                sres='(0,'
            elif side=='R':
                sres='(1,'
            else:
                sres='('
            for k in i:
                sres=sres+k+','
            sres=sres[0:-1]
            sres=sres+')='
            svecargs.append(sres)
        mm=''
        if format=='list':
            for i,j in zip(svecargs,vecv):
                mm= 'args'+i+latex(j)  
                display(Math(mm))
        else:        
            for i,j in zip(svecargs,vecv):
                mm=mm+ '[ '+i+latex(j)+'] ,' 
            display(Math(mm))  
        
def arglist(expr,deep=3):
     
    infoexpr=[]
    infopos=[]
    cc=''
    A,B,C,D=ruta(expr,infoexpr,infopos,cc)
    BC=[[i,j] for i,j in zip(B,C) if  not Is_NMono(i)]
    B1,C1=[],[]
    for i in BC:
        B1.append(i[0])
        C1.append(i[1])
    BC2=[[i,j] for i,j in zip(B1,C1) if  len(j)<deep+1]
    B2,C2=[],[]
    for i in BC2:
        B2.append(i[0])
        C2.append(i[1])
    BC3=[[i,j] for i,j in zip(B2,C2) if not Is_Inverse(i)]
    B3,C3=[],[]
    for i in BC3:
        B3.append(i[0])
        C3.append(i[1])
    return B3,C3
def ruta(expr,infoexpr,infopos,cc):
    mm=expr.args
    if len(mm)>0:
        for i in range(len(mm)):
            nexp=mm[i]
            npos=cc+str(i)
             
            if nexp not in infoexpr:
                if str(type(nexp)) not in vecreatr :
                    if nexp not in infoexpr:
                        if not Is_Number(nexp):
                             
                            infoexpr.append(nexp)
                            infopos.append(npos)
                            try:
                                nexp,ninfo,ncc=ruta(nexp,infoexpr,infopos,npos)
                                return nexp,ninfo,ncc
                            except:
                                pass
        return  expr,infoexpr,infopos,cc  
    else:
        return  expr,infoexpr,infopos,cc

def dothis(*args):
    if 'numer' in args or 'denom' in args:
        args2=args
        expr=args[0]
        
        done=None
        for i in args:
            if i=='numer':  
                done='numer'
            if i=='denom':
                done='denom'
        if Is_Div(args[0]) and  done!=None:
            p1,p2=fraction(expr)
            args2=[i for i in args if i!='numer']
            args3=[i for i in args2 if i!='denom']
            if done=='numer':
                args3[0]=p1
                kres=dothis(*args3)
                kres=unisymbols(Div(kres,p2))
            else:
                args3[0]=p2
                kres= dothis(*args3) 
                kres=unisymbols(Div(p1,kres))
                 
            return kres      
    else:           

        expr=args[0]
        func=args[1]
        if len(args)==2:
            try:
                sexpr=func+'('+str(expr)+')'
                kres=eval(sexpr)
                return unisymbols(kres)
            except:
                 
                return expr

        if len(args)==3:
            try :
                expr2=args[2]
                sexpr=func+'('+str(expr)+','+str(expr2)+')'
                return unisymbols(eval(sexpr))
            except:
                 
                return expr
    
    
def helplib(kimg=''):
    if kimg=='':
        return get_filesh()
    else:
        try:
            kfile=kimg+'.png'
            camino='imghelp/'+kfile
            baseA = np.asarray(Image.open(camino))
            plt.imshow(baseA)
            plt.axis('off')
            plt.gcf().set_facecolor("black")
            plt.rcParams['figure.dpi'] = 800
            plt.show()
        except:
            print( "don't find function, please type helplib() for more options")
    
import os

def helplib(kimg=''):
    
    if kimg=='':
        return get_filesh()
    else:
        try:
            kfile=kimg+'.png' 
            camino='/imghelp/'+kfile
             
            baseA = np.asarray(Image.open(camino))
            plt.imshow(baseA)
            plt.axis('off')
            plt.gcf().set_facecolor("black")
            plt.rcParams['figure.dpi'] = 800
            plt.show()
        except:
            print( "don't find function, please type helplib() for more options")
 

#  CONJUNTOS
def UnionL(A,B):
    C=[x for x in A if x not in B]
    return C+B
def SubstracL(A,B):
    C=[x for x in A if x not in B]
    return C     
def IntersecL(A,B):
    C=[x for x in A if x  in B]
    return C

def maxint(sexpr):
    sval='ǁ'+sexpr+'ǁ'
    return symbols(sval)

def supersubs(expr,v1,v2):
    done_exp=False
    sexpr=str(expr)
    if 'exp' in sexpr:
        done_exp=True
        sexpr=sexpr.replace('exp','M')
    sv1=str(v1)
    sv2='('+str(v2)+')'
    sexpr=sexpr.replace(sv1,sv2)
    if done_exp:
        sexpr=sexpr.replace('M','exp')
    return parse_expr(sexpr)

def trinom2binom(expr,sexpr):
    fexpr=factor(sexpr)
    expr=expr.subs(sexpr,AA)
    expr=expr.subs(AA,fexpr)
    return expr    