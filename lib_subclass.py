from sympy import *
from IPython.display import display, Math


Fz,Agv,V,H,W,Ec,Ep,M,G= symbols('Fz Agv V H W Ec Ep M G')

def sex2rad(k):
    k=simplify(pi*k/180)
    return(k)

def rad2sex(k):
    k=simplify(180*k/pi)
    return(k.evalf())

def sex2rad_i(kang,s='r'):
    if s=='s':
        kang=sex2rad(kang)
    return(kang)
class cvector:
    def __init__(self,fz=Fz,agv=Agv,typea='r'):
                 
        self.fz=fz
        self.agv=sex2rad_i(agv,typea)
        
    
def sum_fx(mm,s='s'):
    qq=len(mm)
    kres=0
    for i in range(qq):
        vec=mm[i]
        aa=vec.agv
        ff=vec.fz
        fx=ff*cos(aa)
        kres+=fx
    return(kres)

def sum_fy(mm,s='s'):
    qq=len(mm)
    kres=0
    for i in range(qq):
        vec=mm[i]
        aa=vec.agv
        ff=vec.fz
        fy=ff*sin(aa)
        kres+=fy
    return(kres)

def Ekinetic(v=V,m=M,g=G,s='s'):
    Ek=m*v**2
    Ek=Ek/2
    if s=='e':
            Ek=Ek.evalf()
    return(Ek)
    
def Epotencial(h=H,m=M,g=G,s='s'):
    pp=m*g*h
    if s=='e':
            pp=pp.evalf()
    return(pp)
def filter_aswer(kans,kparameter='knum'):
    qq=len(kans)
    mm=[]
    for i in range(qq):
        kval=kans[i]
        if kparameter=='trig':
            if kval>=0 and kval <= pi/2:
               mm.append(kval) 
            
        else:
            if(kval.is_Float or kval.is_Integer):
                mm.append(kval)
    return(mm)        
        
def clean_res(k,s='s'):
    kk=k[0]
    if s=='e':
        kk=float(kk)
    return(kk)    
    
def csolve(Eq,kv,kd=' ',kabs=False,s='s',noncero=True,kpositive=False,kdisp=False):
    sR=kd+' ='
    res=solve(Eq,kv)
    if noncero and len(res)==2:
        if res[0]==0:
            kres=res[1]
        else:
            kres=clean_res(res,s=s)
    else:
       kres=clean_res(res,s=s) 
    if kabs:
        kres=abs(kres)
    if kpositive:
        if kres.compare(0)==-1:
            kres=-1*kres
    if s=='e':
           kres=kres.evalf()
    if s=='d':
           kres=float(kres)
    if kdisp:       
        display(Math(sR+latex(kres)))     
    return(kres)
        
def csimplify(rv1,rv2):
    l1=list(rv1.args)
    l2=list(rv2.args)
    l3=list(set(l1) & set(l2))   
    kres=1
    for i in range(len(l1)):
        k=l1[i]
        if k  not in l3:
            if kres==1:
                kres=k
            else:
                kres=kres*k
    rvv1=kres
                
    kres=1
    for i in range(len(l2)):
        k=l2[i]
        if k  not in l3:
            if kres==1:
                kres=k
            else:
                kres=kres*k                    
    rvv2=kres
    return(rvv1,rvv2) 
       
def kclean(kfunc,kv):
    kres=kfunc.subs(kv,1)
    return(kres)

def ssolve(Eq,kv,kd=' ',kabs=False,s='s',keval=' ',noncero=True):
    sR=kd+' ='
    res=solve(Eq,kv)
    if noncero and len(res)==2:
        if res[0]==0:
            kres=res[1]
        else:
            kres=clean_res(res,s=s)
    else:
       kres=clean_res(res,s=s)
    if keval!=' ':
        p1=keval[0]
        p2=keval[1]
        kres=super_subs(kres,p1,p2)
    if kabs:
        kres=abs(kres)
    if s=='e':
           kres=kres.evalf()
    if s=='d':
           kres=float(kres)
    display(Math(sR+latex(kres)))     
    return(kres)
    
    
 
def ctrigremp(Eq1,a1,a2='37'):
    v35=Rational('3/5')
    v45=Rational('4/5')
    v34=Rational('3/4')
    v43=Rational('4/3')
    if a2=='53':
        Eq1=Eq1.xreplace({cos(a1):Rational('3/5')})
        Eq1=Eq1.xreplace({sin(a1):Rational('4/5')})
        Eq1=Eq1.xreplace({tan(a1):Rational('4/3')})
    
    elif a2=='45':
        Eq1=Eq1.xreplace({cos(a1):sqrt(2)/2})
        Eq1=Eq1.xreplace({sin(a1):sqrt(2)/2})
        Eq1=Eq1.xreplace({tan(a1):1})    
    
    elif a2=='60':
        Eq1=Eq1.xreplace({cos(a1):Rational('1/2')})
        Eq1=Eq1.xreplace({sin(a1):Rational('sqrt(3)/2')})
        Eq1=Eq1.xreplace({tan(a1):Rational('sqrt(3)/3')})
        
    elif a2=='30':
        Eq1=Eq1.xreplace({cos(a1):sqrt(3)/2})
        Eq1=Eq1.xreplace({sin(a1):1/2})
        Eq1=Eq1.xreplace({tan(a1):sqrt(3)/3})
    
    
    else:
        Eq1=Eq1.xreplace({cos(a1):Rational('4/5')})
        Eq1=Eq1.xreplace({sin(a1):Rational('3/5')})
        Eq1=Eq1.xreplace({tan(a1):Rational('3/4')})    
    return(Eq1)
    
def catan(k):
    k=atan(k)
    kres=rad2sex(k)
    return(kres)
 
def casin(k):
    k=asin(k)
    kres=rad2sex(k)
    return(kres) 
    
def cacos(k):
    k=acos(k)
    kres=rad2sex(k)
    return(kres)

def csubs(Eq1,ks,kv,kd=' ',s='s'):
    kres=Eq1.subs(ks,kv)
    sR=kd+' ='
    if s=='e':
        kres=kres.evalf()
    if s=='t':
        kres=trigsimp(kres)
    display(Math(sR+latex(kres)))    
    return(kres)

def show_res(Eq1,kd=' ',s='s'):
    kres=Eq1
    sR=kd+' ='
    if s=='e':
        kres=kres.evalf()
    if s=='t':
        kres=trigsimp(kres)
    display(Math(sR+latex(kres)))    
    return(kres)

def super_subs(Eq1,tup1,tup2):
    qq=len(tup2)
    for i in range(qq):
        ss=tup1[i]
        sv=tup2[i]
        Eq2=Eq1.subs(ss,sv)
        Eq1=Eq2
    return(Eq1)

def s_subs(Eq1,keval):
    tup1=keval[0]
    tup2=keval[1]
    kres=super_subs(Eq1,tup1,tup2)

    return(kres) 
    
def changue_val(nval,ksym):
    ksym=nval
