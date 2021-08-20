from sympy import * 
from IPython.display import display, Math  
 

def show_res(Eq1,kd='',kr=''):
    kres=Eq1
    sR=kd+' ='
    display(Math(sR+latex(kres)))
    if kr!='':
        return(kres)

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

def rpow(ksym,op1=2,op2='',kope=''):
    if type(op2)==int:
        ke=frs(op1,op2)
        kope=kope
    else:
        ke=frs(1,op1)
        kope=op2
    kres=(pow(ksym,ke))
    kres=opemat(kres,kope)
    return(kres)
    
def kpow(ksym,op1='',op2='',kope=''):
    if type(op2)==int:
        ke=frs(op1,op2)
        kope=kope
    else:
        ke=op1
        kope=op2
    kres=(pow(ksym,ke))
    kres=opemat(kres,kope)
    return(kres)

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

def sympy_unifiq(keq,ksym):
    mm=fpoly(keq,'free')
    qq=len(mm)
    for i in range(qq):
        if str(ksym)==str(mm[i]):
            return(mm[i])
        
    return(ksym)

def str2sym(kstr):
    return (parse_expr(kstr))
    
def kclean(keq):  # Limpia la equation keq de factores repetidos
    ee1=factor(keq)
    ee2=ee1
    vec1=list(ee1.args)
    vec2=fpoly(ee1,'free')
    qq=fpoly(ee1,'n')
    if type(ee1)==Mul:
        if qq>0:
            for i in vec2:
                ee2=simplify(ee2/i)
         
    return ee2    

def ksolve(Eq1,ksym):
    Eq1=unisymbols(Eq1)
    ksym=unisymbols(ksym)
    try:
        Eq2=unisymbols(kclean(Eq1))
    except:
        Eq2=Eq1
        
    kres=csolve(Eq2,ksym)
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

def opemat(ksym,kope=''):
    '''
    opemat(Equation,kope=opt)
    opt 
        'f'= factor(Equation),    
        't'= trigsimp(Equation)
        'x'= simplify(Equation)
        'v'= Equation.evalf()
        'a'= apart(Equation)
        'c'= cancel(Equation)
        'E'= Equation.expand(force=True)
        '2' = kill sqrt( x^2 )
    '''
    kres=ksym
    
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
         
        try:
            try:
                kres=float(evalf(kres))
                return(kres)
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
                
                    
                    
    if 'a' in kope:
        kres=apart(kres)
    if 'c' in kope:
        kres=cancel(kres)
    if '2' in kope:
        kres=powsimp(kres**2)
        kres=powsimp(sqrt(kres))
    if 'E' in kope:
        kres=kres.expand(force=True)
     
    if 'F' in kope:
         
        kres=nsimplify(kres)        
    return(kres)
 


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
         
    return(kres)
 
def simple_fac(ksym,vecv):
    kres=ksym
    sexp=0
    for i in vecv:
        exp1=fpoly(kres,'facmono',i) 
        exp2=fpoly(exp1,'simp_fac',i) 
        sexp=sexp+exp2
        kres=kres-exp1 
    return(sexp)    
     
       
    
    i
def polyargs(ksym,kargs):
    kres=ksym
    sexp=0
    for i in kargs:
        numar=int(i)
        kres=kres.args[numar]
    return(kres)

    
def pao(ksym,kargs):
    return(polyargs(ksym,kargs)) 


def paolist(ksym,kop=''):
    done=True
    cc=0
    while done:
        try:
            newb=kop+str(cc)
            kres=pao(ksym,newb)
            show_res(kres,str(cc))
            cc+=1
        except:
            done=False
        
    
def argolist(ksym,kv=''): 
    try:    
        mm=fpoly(ksym,'list')
        nn=0
        for i in mm:
            nkv=kv+str(nn)
            get_args_list(i,nkv)
            nn+=1
            show_res(i,nkv)  
    except:
        show_res(ksym,kv)      
      
# TRANSFORMACIONE CON lambdify


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

def perfect_sqrt_exp(ksym,kop=''):
    klist=fpoly(ksym,'list')
    korden=([0,1,2],[0,2,1],[1,2,0])
    for kk in korden:
        e1=ksymroot(klist[kk[0]]) 
        e2=klist[kk[1]] 
        e3=ksymroot(klist[kk[2]])
    
        if 2*e1*e3==e2:
            kres=e1+e3
            if '1' in kop:
                kres=(e3+e1)
            if '2' in kop:
                kres=kpow(kres,2)
            return(kres)    
        elif -2*e1*e3==e2:
            kres=e1-e3
            if '1' in kop:
                kres=(e3-e1)
            if '2' in kop:
                kres=kpow(kres,2)
            return(kres)    
                
        else:
            kres=ksym 
    return(kres)        

def ksymroot(ksym): #perfect_sqrt_exp
    kres=1
    mm=gek(ksym)
    for i in mm:
        kv=i[0]
        ke=i[1]
        kres=kres*kpow(kv,ke,2)
    return(kres)

def gek(ksym): # used by ksymroot
    veriv=get_var_mono(ksym,1)
    verin=get_num_mono(ksym,1)
    mm=[]
    if verin:
        knum=get_num_mono(ksym)
        mm.append([knum,1])
    if veriv:
        kvar=get_var_mono(ksym)
        if type(kvar)==Symbol:
            mm.append([kvar,1])
        else:    
            if type(kvar)==Pow:
                mm.append([kvar.args[0],kvar.args[1]])
            else:
                klist=fpoly(kvar,'list')
                for i in klist:
                    if type(i)==Pow:
                        mm.append([i.args[0],i.args[1]])
                    else:
                        mm.append([i,1])
                    
    return(mm)                
    
    
def simpify_monoexp(ksym):
    try:  # metodo 1
        mm=fpoly(ksym,'list')
        if len(mm)==2:
            m1=mm[0]
            e1=mm[1]
            m1=opemat(m1,'s')
            mm1=m1.args[0]
            mm2=m1.args[1]
            mm3=opemat(e1*mm2,'s')
            return(kpow(mm1,mm3))
    except:
        return(ksym)
        
        
def simp_exp(ksym=1,op1='h'):
    if op1=='h':
        show_exp(['op=1,',kpow(kpow(x,a),b),'=',kpow(x,a*b)])
        show_exp(['op=2,',kpow(x,a)*kpow(x,b),'=',kpow(x,a+b)])
    elif op1==1:
            if type(ksym)==Pow:
                try:
                    mm=fpoly(ksym,'list')
                    ke1=mm[1] 
                    e2=mm[0] 
                    ke2=e2.args[1] 
                    m1=e2.args[0] 
                    kres=kpow(m1,opemat(ke1*ke2,'s'))
                    return(kres)
                except:
                    return(ksym)
            else:
                return(ksym)
    elif op1==2:
        if type(ksym)==Mul:
            try:
                mm=fpoly(ksym,'list')
                v1=mm[0] 
                v2=mm[1] 
                vv1=v1.args[0] 
                ve1=v1.args[1] 
                vv2=v2.args[0] 
                ve2=v2.args[1] 
                ke=opemat(ve1+ve2,'s') 
                kres=kpow(vv1,ke)
                return(kres)
            except:
                return(ksym)
        else:
            return(ksym)
    else:        
        return(ksym)
    

    
# factorizar Polinomios lineales
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


        
# Operaciones con Monomios    
def monomath(ksym='',op1='',kope=''):
    x,a,b,c=symbols('x a b c')
    kres=ksym
    if op1=='h' or ksym=='h' or ksym=='':
        show_exp(['all letters includinn in op1 will be ejecuted'])
        show_exp(['monomath(f(x),SeOf...) or Eq.polymath(SeOf..) '])
        show_exp(['O','=',kpow(kpow(x,a),b),'to',kpow(x,a*b)])
        show_exp(['S','=',kpow(x,a),'to',kpow(x,(b+c))])
        show_exp(['F','=',kpow(x,a),'to',kpow(x,(b*c))])
        show_exp(['e expand, f factor,  s implify, 1 only left Eq, 2 only right Eq'])
    else:
        for i in op1:
        
            if i=='O':
                try:

                    if type(kres)==Pow:
     
                        kres2,exp1=kres.args[0],kres.args[1]
                        try:
                            if type(kres2)==Pow:
                                
                                 
                                kres3,exp2=kres2.args[0],kres2.args[1]
                                newexp=exp1*exp2
                                newexp=opemat(newexp,kope)
                                kres=kpow(kres3,newexp)
                                 
                                 
                            
                        except:
                            done=False
                except:
                    done=False        
            if i=='S':
                try: 
                    if type(kres)==Pow:
                         
                        kres2,exp1=kres.args[0],kres.args[1]
                        if type(exp1)==Mul:
                             
                            exp2=expand(exp1)
                             
                            if type(exp2==Add):
                                kres=(kpow(kres2,exp2))
                    else:
                        done=False  
                except:
                    done=False                
            if i=='F':
                try :
                    kres2,exp1=kres.args[0],kres.args[1]
                    exp2=factor(exp1)
                    kres=(kpow(kres2,exp2))
                except:
                    done=False
                            
            kres=opemat(kres,i) 
        
        return(kres) 


def simple_m(ksym,op=''):
    '''
    simple_m(ksym,'N') , return numerator from ksyn
    simple_m(ksym,'D') , return denominator from ksyn
    simple_m(ksym,'O') , simplify exponent y pow2 pow , (x**y)**z return x**(y*z)
    simple_m(ksym,'E') , expand exponent   , (x)**(y*z) return x**y**z)
    '''
    
     
    
    
    kres=ksym
    vec='fse'
    
    # factor expand, simplify
    if op in vec:
        kres=opemat(kres,op)
        return(kres)
    # Numerator
    elif op=='N':
        try:
            return(numer(kres))
        except:
            return(kres)
    
    # Denominator
    elif op=='D':
        try:
            return(denom(kres))
        except:
            return(kres)   
    # simplify Pow2Pow
    elif op=='O': #redice the expo if it possible
        try:
            return(Mpowpow2pow(kres))
        except:
            return(kres)
    
    elif op=='E':  # Expand the exponeneteif is possible
        try:
            return(Mpow2powpow(kres))
        except:
            return(kres) 
 
            
    # Sin Accion 
    
    else:
        return(kres) 

def redu_pow2pow(ksym):
    kres=ksym
    
    try:
        if type(kres)==Pow:

            kres2,exp1=kres.args[0],kres.args[1]
            try:
                if type(kres2)==Pow:
                    
                     
                    kres3,exp2=kres2.args[0],kres2.args[1]
                    newexp=exp1*exp2
                    newexp=opemat(newexp,kope)
                    kres=kpow(kres3,newexp)
                     
                     
                
            except:
                done=False
    except:
        done=False
    
# especial functions to fix symbols crazy identifications

def cdiff(kfunc,ksym): # return diff kfunc whit respect to ksym , but unifique all symbols to same
    return diff(unisymbols(kfunc),unisymbols(ksym)) 
    
    

def get_expo_mono1(ksym): # mono 1 = monomio de una sola variable
    kres=ksym
    if type(ksym)==Pow:
            kvar=ksym.args[1]
            return(kvar)
    else:
        return(kres)
        
def get_exp_in_mono(ksym,ksvar=''):
    if type(ksym)==Pow:
        kvar=ksym.args[0]
        kexp=ksym.args[1]
        if kvar==ksvar:
            return(kexp)
        else:
            return(ksym)
    elif type(ksym)==Mul:
        for i in fpoly(ksym,'list'):
            kres=get_exp_in_mono(i,ksvar)
            if kres!=ksym:
                return(kres)
    else:
        return ksysm
                
#   #################################
#   #      New Monomial Functions   #
#   #################################  
      
def get_expo_mono2(ksym,op=''): # mono 2 = mono monomio de varias variables
    
    done=True
    
    if op=='' and type(ksym)==Pow:
        kres=get_expo_mono1(ksym)
        return(kres)
    elif type(ksym)==Mul and op!='':
        klist=fpoly(ksym,'list')
        for i in klist:
            if type(i)==Pow:
                if i.args[0]==op:
                    try:
                        kres1=get_expo_mono1(i)
                        return(kres1)
                    except:
                        done=False
    else:
        return(ksym)

def get_M_expon(ksym,op=''): # get_M_expon(ksym)-> return ksym expont if ksym is one variable
                             # get_M_expon(ksym,x)-> return expont of x in ksym if ksym is multvar
    return(get_expo_mono2(ksym,op=op))
    
def is_pow2pow(ksym): # return True if ksym have (x**y)**z formate
    done=False
    if type(ksym)==Pow:
        if type(ksym.args[0])==Pow:
            done=True
    return(done)
    
def getexpmono(ksym,evar):
    if type(ksym)==Mul:
        for i in fpoly(ksym,'list'):
            if type(i)==Pow:
                vv,ee=pow2list(i)
                if evar==vv:
                    return(ee)
                
    elif type(ksym)==Pow:
        vv,ee=pow2list(ksym)
        if evar==vv:
            return(ee)
    
    else:
        return(ksym)
        
        
def pow2list(ksym):
    kres=ksym
    if type(kres)==Pow:
        kvar,kexp = kres.args[0],kres.args[1]
        return(kvar,kexp)
    else:
        return(kres)
        
 
def Mpowpow2pow(ksym):
    kres=ksym
    if is_pow2pow(ksym):
        var1,exp1=pow2list(ksym)
        if type(var1)==Pow:
            newvar,exp2=pow2list(var1)
            newexp=exp2*exp1
            kres=kpow(newvar,newexp)
            return(kres)
        else:
            return(kres)
    else:
        return(kres)        
        
def Mpow2powpow(ksym):
    kres=ksym
    if type(kres)==Pow:
        kvar,kexp=pow2list(ksym)
        cc=1
        if type(kexp)==Add:
            klist=fpoly(kexp,'list')
            for i in klist:
                cc=(cc)*(kpow(kvar,i))
            return(cc)
        elif type(kexp)==Mul:
            cc=kvar
            klist=fpoly(kexp,'list')
            for i in klist:
                cc=(kpow(cc,i))
            return(cc)    
                
        else:
            return(kres)
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

def jack(ksym,vec):
    vecr=[]
    restp=0
    for i in vec:
        try:
            op1=i[0]
            op2=[]
            for j in i[1:len(i)]:
                op2.append(int(j))
            if op1=='n':
                restp=numer(part(ksym,op2))
            elif op1=='d':
                restp=denom(part(ksym,op2))
            else:
                restp=part(ksym,op2)
        except:
            restp=0
        vecr.append(restp)
    return(vecr)   
    
        
#   #################################
#   #      multiproces Functions    #
#   #      applied polyclass        #
#   #################################  
def multiproc(ksym,vecp,show=False):
     
    i=0
    sEE=[]
    qq=len(vecp)
    sE(['Math procedure in ',ksym.Q])
    while i<qq:
        kpro=vecp[i]
         

        if kpro=='s':
            a=ksym.opemat('s',1)
            sEE.append('simplify->')
            if show:
                sE(['simplify->',ksym.Q])
            i+=1
        elif kpro=='f':
            a=ksym.opemat('f',1)
            sEE.append('factor->')
            if show:
                sE(['factor->',ksym.Q])
            i+=1
        elif kpro=='e':
            a=ksym.opemat('e',1)
            sEE.append('expand->')
            if show:
                sE(['expand->',ksym.Q])
            i+=1
        elif kpro=='S':
            i+=1
            a=ksym.ope4('S',vecp[i],1)
            sEE.append('Add(')
            sEE.append(vecp[i])
            sEE.append(')->')
            if show:
                sE(['Add ',vecp[i],',',ksym.Q])
            i+=1
        elif kpro=='M':
            i+=1
            a=ksym.ope4('M',vecp[i],1)
            sEE.append('Mul(')
            sEE.append(vecp[i])
            sEE.append(')->')
            if show:
                sE(['Mul ',vecp[i],',',ksym.Q])
            i+=1
        elif kpro=='P':
            i+=1
            a=kksym.ope4('P',vecp[i],1)
            sEE.append('Pow(')
            sEE.append(vecp[i])
            sEE.append(')->')
            if show:
                sE(['Pow ',vecp[i],',',ksym.Q])
            i+=1     
    sE(sEE)
    return(ksym)

def unisymbols(ksym):
    try:
        kres=parse_expr(str(ksym))
    except:
        kres=ksym
    return(kres) 

def kreturn(ksym):
    unisymbols(ksym)
    
def dis2point(P1,P2):
    xx1,yy1=P1
    xx2,yy2=P2
    L1=kpow(xx2-xx1,2)
    L2=kpow(yy2-yy1,2)
    LL=rpow(L2+L1)
    return(LL)
    
def CentMass2func(x,f1,f2,x1,x2,kshow=False):
    y1= lambdify(x,f1)
    if f2==0:
        y3= lambdify(x,y1(x))
        y4= lambdify(x, (y1(x))*frs(1,2))
    else:    
        y2= lambdify(x,f2) 
        y3= lambdify(x,y1(x)-y2(x))  
        y4= lambdify(x, (y1(x)+y2(x))*frs(1,2))
    
    At= integrate(y3(x),(x,x1,x2)) 
    Cx= opemat((integrate(x*y3(x),(x,x1,x2)))/At,'sfF')
     
    Cy= opemat(integrate(y3(x)*y4(x),(x,x1,x2))/At,'sfF')
    if kshow:
        show_res(Cx,'Cx')
        show_res(Cy,'Cy')
    return([Cx,Cy])

def CentMass2func2(x,f1,f2,x1,x2,kshow=False):
    y1= f1
    y2= f2 
    y3= y1-y1  
    y4=  (y1+y2)*frs(1,2) 
    
    At= integrate(y3,(x,x1,x2)) 
    Cx= opemat((integrate(x*y3,(x,x1,x2)))/At,'sfF')
     
    Cy= opemat(integrate(y3*y4,(x,x1,x2))/At,'sfF')
    if kshow:
        show_res(Cx,'Cx')
        show_res(Cy,'Cy')
    return([Cx,Cy])

def MassCe(ksym,y1,y2,x1,x2,kope='',kshow=False):

    At=integrate((y1-y2),(ksym,x1,x2)) # Total Area
    y3= y1-y2 # Area Dif
    y4=(y1+y2)*frs(1,2) # Midel point
    Cx=integrate(ksym*y3,(ksym,x1,x2))/At
    Cx=opemat(Cx,kope)
    Cy=integrate(y3*y4,(ksym,x1,x2))/At
    Cy=opemat(Cy,kope)
    if kshow:
        show_res(Cx,'Cx')
        show_res(Cy,'Cy')
    return [Cx, Cy]

x,fx=symbols('x fx')

class dataFunc: # tools for get slope,angle etc on function 
    def __init__(self,x=x,fx=fx):
                 
        self.x=x  
        self.f=fx
        self.fx=fx

    
    def evalue(self,kval='',kope=''):
        kres=self.fx
        ksym=self.x
        if kval=='':
            return kres
        else:
            kres= kres.subs(ksym,kval) 
        kres=opemat(kres,kope)
        return kres
 
    
    def slope(self,kval='',kope=''):
        xExp=self.fx
        ksym=self.x
        kres=diff(xExp,ksym)
        if kval!='':
            kres=kres.subs(ksym,kval)
        kres=opemat(kres,kope)    
        return kres
    
    def slopeO(self,kval='',kope=''):
        kres=self.slope(kval=kval,kope='')
        kres=-1/kres
        kres=opemat(kres,kope)    
        return kres 
    
    def angTan(self,kval='',kope='',ksex=False):
        kres=self.slope(kval=kval,kope=kope)
        
        kres=opemat(kres,kope)
        kres=atan(kres)
        if ksex:
            kres=rad2sex(kres)
        return kres 
    
    def angOrto(self,kval='',kope='',ksex=False):
        kres=self.slopeO(kval=kval,kope='')
        
        kres=opemat(kres,kope)
        kres=atan(kres)
        if ksex:
            kres=rad2sex(kres)
        return kres