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
    
    
         
#R#
def clean_res(k,s='s'):
    kk=k[0]
    if s=='e':
        kk=float(kk)
    return(kk)    
 

    
def csolve(Eq,kv,kd=' ',kabs=False,s='s',noncero=True,kpositive=False,kdisp=False,ksolveset=False):
    sR=kd+' ='
    if ksolveset :
        kres=solveset(Eq,kv)
        if len(kres)==2:
            [kres1,kres2]=kres
            if kpositive:
                if kres1>0:
                    return(kres1)
                if kres2>0:
                    return(kres2)
            else:
                return(kres1,kres2)
    
    else:    
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
    else:    
        return(kres)
    
def kdisplay(skres):
        display(Math(latex(skres)))
        

def supersotore(kclass,kvar,kval):
    mm=kclass.kvalue
    

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
    #return(kres)

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

def opemat(kres,kope=''):
    if 'e' in kope:
        kres=expand(kres)
    if 'f' in kope:
        kres=factor(kres)    
    if 't' in kope:
        kres=trigsimp(kres)
    if 'x' in kope:
        kres=expand_trig(kres)    
    if 's' in kope:
        kres=simplify(kres)
    if 'v' in kope:
        kres=float(kres)     
    return(kres)


def kgetsol(ksym,eq1,eq2,kope=''):
    kres=csolve(eq1-eq2,ksym,kdisp=False)
    kres=opemat(kres,kope=kope)
    return(kres)


#R#
def designe_val(Eq1,kd=' ',kope=''):
    kres=Eq1
    sR=kd+' ='
    kres=opemat(kres,kope)
    display(Math(sR+latex(kres)))
    return(kres)    