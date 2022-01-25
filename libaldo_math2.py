from sympy import *


from functools import reduce
#from IPython.display import  Math ,display 
import matplotlib.pyplot as plt


from  libaldo_show import *
 



# VARIABLE 

def unisymbols(ksym):
    try:
        kres=parse_expr(str(ksym))
    except:
        kres=ksym
    return(kres) 


 
        
# Math Functions

def sqrs(k1,k2): # return pow(k1,S('k2'))
    v1=str(k2)
    v2='1/'+v1
    kres=pow(k1,S(v2))
    return(kres)

def frs(k1,k2): # return   S('k1/k2')
    v1=str(k1)
    v2=str(k2) 
    v3=v1+'/'+v2
    return(S(v3))

def rpow(ksym,op1=2,op2='',kope=''): # return  rpow(a,b) return a**(1/b)
    if type(op2)==int:
        ke=frs(op1,op2)
        kope=kope
    else:
        ke=frs(1,op1)
        kope=op2
    kres=(pow(ksym,ke))
    kres=opemat(kres,kope)
    return(kres)
    
def kpow(ksym,op1='',op2='',kope=''): # return  kpow(a,b) return a** b 
    if type(op2)==int:
        ke=frs(op1,op2)
        kope=kope
    else:
        ke=op1
        kope=op2
    kres=(pow(ksym,ke))
    kres=opemat(kres,kope)
    return(kres)

def sex2rad(k):  # convert segsadesimal to radial
    k=simplify(pi*k/180)
    return(k)

def rad2sex(k):  # convert radial to sexages.
    k=simplify(180*k/pi)
    return(k.evalf())

def sex2rad_i(kang,s='r'):
    if s=='s':
        kang=sex2rad(kang)
    return(kang)
   
def killAbs(ksym):
    kres=ksym 
    try:
        mm=str(kres)
        mm=mm.replace('Abs','')
        return parse_expr(mm)
    except:
        return kres
 
def get_inside_root(ksym):
    kres=ksym 
    if Is_Root(kres): 
        mm=fpoly(kres,'list')
        return mm[0]
    return kres 
    
def get_inside_Pow(ksym):
    kres=ksym 
    if Is_Pow2(kres): 
        mm=fpoly(kres,'list')
        return mm[0]
    return kres

def get_expo(ksym): # return exponente from monomie expresion
    if Is_Mono(ksym):
        mm=fpoly(ksym,'list')
        return mm[1]

def get_killexpo(ksym): # return exponente from monomie expresion
    if Is_Mono(ksym):
        mm=fpoly(ksym,'list')
        return mm[0]        

#  GEOMETYRY Util Formules

def get_hipo(a,b,kope=''): # return raiz(a*a+b*b)
    kres=rpow(kpow(a,2)+kpow(b,2),2)
    kres=opemat(kres,kope=kope)
    return kres

def get_cateto(a,b,kope=''):# return raiz(a*a-b*b)
    kres=rpow(kpow(a,2)-kpow(b,2),2)
    kres=opemat(kres,kope=kope)
    return kres 
    
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
    

    
    
def x_pol(kr,kang):
    return kr*cos(kang)
    
def y_pol(kr,kang):
    return kr*sin(kang)    
    
#  SOLVE
def lenght2point(x1,y1,x2,y2,kope=''):
    l1=x2-x1
    l2=y2-y1
    kres=rpow(l1*l1-l2*l2)
    kres=opemat(kres,kope=kope)
    return kres
    
def signo(ksym):
        kres=ksym 
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
def csolveR(Eq1,ksym,kd='',kope=''): # solve 2do grade Eq return positive part
    nksym=kpow(ksym,2)
    kres=csolve(Eq1,nksym,kope=kope)
    kres=opemat(kres,kope)
    kres=rpow(kres,2)
    if kd!='':
        sR=kd+' ='
        display(Math(sR+latex(kres)))        
         
    return(kres)
    
def csolve(Eq1,ksym,kd='',korden='',kpositive=False,kope='',kdisp=False,unifique=False):
    Eq1=unisymbols(Eq1)
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

#  SIMPLIFICATION



def opemat(ksym,kope=''):
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
        if i=='-' :
            kres2=kres
            try:
                kres2=kpow(kres,-1)
                kres2=opemat(kres2,'r')
                kres=kpow(kres2,-1)
            except:    
                done=False
                
        if i=='2' :
            kres=kpow(kres,2)
            kres=rpow(kres,2)
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

    
def factorSec(kEq,ksym,kfiltro='.'):
    if type(ksym)==list:
        return MgrupFac(kEq=kEq,ksym=ksym,kfiltro=kfiltro)
    else:
        return grupFac(kEq=kEq,ksym=ksym,kfiltro=kfiltro)

def grupFac(kEq,ksym,kfiltro='.'):
    return My_factor(kEq=kEq,ksym=ksym,kfiltro=kfiltro)
    
def MgrupFac(kEq,ksym,kfiltro='.'):
    Kres=0
    for i in ksym:
        kEq,kfin=My_factor2(kEq=kEq,ksym=i,kfiltro=kfiltro)
        Kres+=kfin
        
    return Kres+kEq    
    
def My_factor(kEq,ksym,kfiltro='.'):
    # My_factor(x*x*y+x*x*z+x*y+x*z,x*x) = return x*x*(y+z)+x*y+x*z

    # My_factor(x*x*(y+z)+x*y+x*z,x) = return x*(x*(y+z)+y+z)
    #    but....
    # My_factor(x*x*(y+z)+x*y+x*z,x,'+') = return x*x*(y+z) +x*(y+z)
    
    
    mm=fpoly(kEq,'list')
    mm1=[]
    mm2=[]
    vstr1=str(ksym)
    kres1=0
    kres2=0
    for i in mm:
        if Is_Mono(i) :
            vstr2=str(i)
            if vstr1 in vstr2 and kfiltro not in vstr2 :
                mm1.append(simplify(i/ksym))
                kres1+=simplify(i/ksym)
            else:
                mm2.append(i)
                kres2+=i
        else:
            mm2.append(i)
            kres2+=i
    
    kres=kres2+ksym*kres1 
    return  kres

    
def My_factor2(kEq,ksym,kfiltro='.'):
    # My_factor(x*x*y+x*x*z+x*y+x*z,x*x) = return x*x*(y+z)+x*y+x*z

    # My_factor(x*x*(y+z)+x*y+x*z,x) = return x*(x*(y+z)+y+z)
    #    but....
    # My_factor(x*x*(y+z)+x*y+x*z,x,'+') = return x*x*(y+z) +x*(y+z)
    
    
    mm=fpoly(kEq,'list')
    mm1=[]
    mm2=[]
    vstr1=str(ksym)
    kres1=0
    kres2=0
    for i in mm:
        if Is_Mono(i) :
            vstr2=str(i)
            if vstr1 in vstr2 and kfiltro not in vstr2 :
                mm1.append(simplify(i/ksym))
                kres1+=simplify(i/ksym)
            else:
                mm2.append(i)
                kres2+=i
        else:
            mm2.append(i)
            kres2+=i
    
    kres=(kres2,ksym*kres1) 
    return  kres

def linear_fac(ksym,kvec):
    kfinal=0
    kres=ksym
    for i in kvec:
        poly0=kres
        poly1=0
        poly2=0
        for j in fpoly(poly0,'list'):
            if i in fpoly(j,'free'):
                poly1+=j
            else:
                poly2+=j
        kfinal+=factor(poly1)
        kres=poly2
    kfinal+=kres
    return(kfinal) 
      
def polytrans(op1='',k1='',k2='',k3='',kdis=''):
    x,y,z=symbols('x y z')
    
    if op1=='h':
        
        show_res(ae.Equation((x+y)*(x**2-x*y+y**2),(x**3+y**3)),'op9')
        show_res(ae.Equation((x**3+y**3),(x+y)*(x**2-x*y+y**2)),'op10')
        show_res(ae.Equation(x**3+y**3+z**3-3*x*y*z,(x+y+z)*(x**2+y**2+z**2-x*y-y*z-z*x)),'op11')
        
    else:
        x=k1
        y=k2
        z=k3
        if op1=='9':
            kres=(x**3+y**3)
        if op1=='10':
            kres=(x+y)*(x**2-x*y+y**2) 
        if op1=='11':
            kres=(x+y+z)*(x**2+y**2+z**2-x*y-y*z-z*x) 

        return(kres)
def make_grp(ksym,kvar):
    kpor,kpoc=0,0
    kres=ksym
    if type(kres)==Add:
        kpor,kpoc=0,0
        klist=fpoly(kres,'list')
        # print(klist)
        for i in klist:
            kcos=simplify(i/kvar)
            if (denom(kcos))==1:
                kpoc+=kcos
            else:
                kpor+=i
        kres=kvar*(kpoc)+kpor
        return(factor(kvar*kpoc)+kpor)
    else:
        return(kres)

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


# Get Info

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
        sE(['Is Polynomie = ',Is_Poly(ksym),' Is Symbols= ',Is_Symbol(ksym),' Is Number= ',Is_Number(ksym)]);
        sE(['Is Real= ',Is_Real(ksym),'  Is Integer= ',Is_Integer(ksym),' Is Even= ',Is_Even(ksym),' Is Odd= ',Is_Odd(ksym)]);
        sE(['Is Monomie= ',Is_Mono(ksym),'  Is Add= ',Is_Add(ksym),' Is Mul= ',Is_Mul(ksym)]);
        sE(['Is Pow=',Is_Pow(ksym),'  Is Pow2= ',Is_Pow2(ksym),' Is Root= ',Is_Root(ksym)])
        
def Is_Poly(ksym):
    done=False
    if type(ksym)==Mul or type(ksym)==Pow:
        done= False
    if Is_Add(ksym):
        kk=fpoly(ksym,'list')
        xx=[Is_Mono(x) for x in kk]
        if True in xx:
            done= True
    else:
        try:
            kn=len(fpoly(ksym,'list0'))
            if kn>1:
                done= True
            else:
                done= False
        except:
            done= False
    if done:
        if Is_Mono(ksym):
            mm=fpoly(ksym,'list')
            kres2=1
            for i in mm:
                kres2=kres2*i
            if kres2==ksym:
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
    return TrFa(sympify(ksym).is_number)

def Is_Real(ksym):
    return TrFa(sympify(ksym).is_real) 
   
def Is_Integer(ksym):
    return  TrFa(sympify(ksym).is_integer)
    
def Is_Even(ksym):
    return TrFa(sympify(ksym).is_even ) 
    
def Is_Odd(ksym):
    return TrFa(sympify(ksym).is_odd ) 

def TrFa(kval): # is True False
    if kval==True or kval==False:
        return(kval)
    else:
        return False
    
def Is_Mono(ksym):
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
        
def Is_Pow(ksym):
    kres=ksym 
    if type(kres)==Pow:
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

def Is_Root(ksym):
    try:
        mm=fpoly(ksym,'list')
        p1=mm[0]
        p2=mm[1]
        if type(ksym)==Pow and p2==1/2:
            return True
        else:
            return False
    except:
        return False
        
def Is_MonoMul(ksym):
    if type(ksym)==type(expand(ksym)):
        return True
    else:
        return False
        
def Is_MonoPoly(ksym):
    return not Is_MonoMul(ksym)
    
def Is_Integral(ksym):
    mm=str(ksym.args)
    if (mm[-3::])==',))':
        return True
    else:
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

def alphaname(ksym):
    kk=str(ksym)
    if kk=='alpha':
        return 'α'
    if kk=='beta':
        return 'ß'
    else:
        return kk
        
###  geometria Analitica

def eQrec(x1=0,y1=0,x2=0,y2=0,var2=''):
    mm=frs((y2-y1),(x2-x1))
    bb=y2-x2*mm
    kres=opemat(var2*mm+bb,'s')
    return kres
 
def get_diff_name(k1,k2):
    xx='d_'+k1
    tt='d_'+k2
    xt='\\frac{'+xx+'}{'+tt+'}'
    dxt=symbols(xt)
    return dxt  

def get_ddiff_name(k1,k2):
    xx='d^{2}'+k1
    tt='d'+k2+'^{2}'
    xt='\\frac{'+xx+'}{'+tt+'}'
    dxt=symbols(xt)
    return dxt    


def put_dot_in(ssym):  
    kres=ssym
    kname='\dot{'+kres+'}'
    return symbols(kname)
    
def put_ddot_in(ssym):
    kres=ssym 
    kname='\ddot{'+kres+'}'
    return symbols(kname)
    
def kprint(val):
    sE([val])    