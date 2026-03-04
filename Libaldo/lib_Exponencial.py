
from sympy import *
from lib_Mathbasic import *
from lib_Algorith import *

from sympy import symbols, Add, Mul, Pow, sympify
'''
x, y = symbols('x y')

def fungtp(expr):
    # Si es lista o tupla
    if isinstance(expr, (list, tuple)):
        return type(expr)(fungtp(a) for a in expr)

    # Si es expresión simbólica con hijos
    if hasattr(expr, "args") and expr.args:
        new_args = [fungtp(a) for a in expr.args]
        return expr.func(*new_args)

    # Caso base: hoja (número, símbolo)
    try:
        return expr % 5
    except Exception:
        return expr


'''
 

 
 
###  CORE
def sfunc2func(expr,kfunc):
    '''
    args:
        expr=math exprssion
        kfunc= str,function name
    return:
        kfunc(expr)
    '''    
    sexpr=kfunc+'('+str(expr)+')'
    expr2=parse_expr(sexpr)
    return expr2

def corefunc(expr,kfunc):
    ''' 
    simplify loop of:
        Is_Add
        Is_Div
        Is_Mul
    '''    
        
    expr=unisymbols(expr)
    if Is_Add(expr):
        kres=0
        for i in expr.args:
            kres=kres+sfunc2func(i,kfunc)
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        kres=Div(sfunc2func(p1,kfunc),sfunc2func(p2,kfunc))
    elif Is_Mul(expr):
        kres=1
        for i in expr.args:
            kres=kres*sfunc2func(i,kfunc)
    else:
        kres=expr
    return kres

### End Core
 
def simplesimplifypow(expr):
    return getbase(expr)**getexpo(expr)
  
def expandexpoexpand(expr): # simplify each exponent  in expr
    kres=corefunc(expr,'expandexpo')
    if Is_Pow(kres):
        ee=getexpo(kres)
        base=getbase(kres)
        if Is_Pow(ee):
            ee=expandexpo(ee)
            kres= bb**ee
    return kres


def Is_MulPow(expr):
    if type(expr)==Mul:
        for data in expr.args:
            if Is_Pow(data):
                return True
    return False 

def mulexpo(expr,force=True): # simplify each ecponet in expr
    expr=unisymbols(expr)
    if type(expr)==Add:
        kres=0
        for data in expr.args:
            kres+=mulexpo(data,force=force)
        return kres
    elif type(expr)==Mul:
        kres=1
        for data in expr.args:
            kres=kres*mulexpo(data,force=force)
        return kres         
    elif Is_Pow(expr) and Is_Pow(getbase(expr)):
        ee=getexpo(expr)
        bb=getbase(expr)
        bbb=getbase(bb)
        eee=getexpo(bb)
        return bbb**(ee*eee)

    elif Is_Pow(expr) and Is_MulPow(getbase(expr)):
        bb=getbase(expr)
        ee=getexpo(expr)
        kres=1
        for data in bb.args:
            if Is_Pow(data):
                bb2=getbase(data)
                ee2=getexpo(data)
                kres=kres*bb2**(ee2*ee)
            else:
                kres=kres*data**ee
        return kres        
    elif Is_Pow(expr) and Is_Mul(getbase(expr)):
        ee=getexpo(expr)
        bb=getbase(expr)
        kres=1
        for data in bb.args:
            kres=kres*mulexpo(data**ee)
 
        return kres    
    else:
        return powdenest(expr,force=force)       
        

 




 

def partPow(expr): # if ksym = sqrt(x**3) x,3,2 ,  if ksym= x**3, return x,3,1 ,if is x ret x,1,1
    if Is_Pow(expr):
        if Is_Root(expr):
            p1,p2=expr.args
            nP=numer(p2)
            dP=denom(p2)
            bb=p1
            ee=nP
            rr=dP
            if Is_Pow(bb):
                ff=getexpo(p1)
                bb=getbase(p1)
                ee=ee*ff
        else:
            ee=getexpo(expr)
            bb=getbase(expr)
            rr=1
    else:
        bb=expr
        ee=1
        rr=1
    return bb,ee,rr

def partpow(expr): # if ksym = sqrt(x**3) x,3,2 ,  if ksym= x**3, return x,3,1 ,if is x ret x,1,1
    if Is_Pow(expr):
        return getbase(expr),getexpo(expr)
    else:
        return expr,1    


def pow2powpow(*args):
    expr=args[0]
    ktype='out'
    exp1=''
    if len(args)==2:
        exp1=args[1]
    if len(args)==3:
        ktype=args[2]
         

           
    if Is_Add(expr):
        mm=0
        for i in fpoly(expr,'list'):
            mm=mm+pow2powpow(i,exp1)
        return mm
    elif Is_Div(expr):
        return cfrac(pow2powpow(numer(expr),exp1),pow2powpow(denom(expr),exp1))
    elif Is_Mul(expr):
         
        mm=1
        for i in fpoly(expr,'list'):
            mm=mm*pow2powpow(i,exp1)
        return mm    
    elif Is_Pow(expr):
        ee=getexpo(expr)
        bb=getbase(expr)
        
        if Is_Mul(ee):
        
            mm2=fpoly(ee,'list')
            if exp1!='' and exp1 in mm2:
                ee1=exp1
                ee2=simplify(ee/ee1)
                if ktype!='out':
                    return((bb**ee2)**ee1)
                else:
                    return((bb**ee1)**ee2)
            else:
                 
                ee1=mm2[0]
                ee2=simplify(ee/ee1)
                return((bb**ee1)**ee2)
        else:
            return expr
         
    else:
        return expr 

        
def powexpand(ksym,op=''):
        '''
        input (x**(a*b))   ---->   return(x**a)**b
        input (x**(a*b),b)   ---->   return(x**b)**a
        '''
        kres=ksym 
        if Is_Number(ksym):
            
            return ksym
        
        
        if Is_Div(ksym):
             
            p1=powexpand(numer(ksym))
            p2=powexpand(denom(ksym))
             
            return p1/p2 
             
        
        elif Is_Mul(ksym):
            mm=fpoly(ksym,'list')
            kres=1
            for i in mm:
                kres=kres*powexpand(i,op=op)
             
            return kres
             
        elif type(ksym)==Pow:
            mm=fpoly(ksym,'list')
            val=mm[0]
            vexp=mm[1]
            if type(vexp)==Pow:

                return ksym
            elif type(vexp)==Mul:
                mm2=fpoly(vexp,'list')
                if len(mm2)==2:
                    p1=mm2[0]
                    p2=mm2[1]
                    if op!='':
                        kres=(val**p2)
                        kres=kres**p1
                    else:    
                        kres=(val**p1)
                        kres=kres**p2
                    return kres
                
            else:
                return ksym
             

             
        elif Is_Add(ksym):
            mm=fpoly(ksym,'list')
            mmr=0
            for i in mm:
                mmr+=powexpand(i,op=op)
            return mmr
        else:
            return ksym
            

    
def powsimplifybase(ksym,combine='base'):
        '''
        input (x**y*x**z*y**z   ---->   returnx**(x+y)*y**z
         
        '''
        expr=ksym 
        if Is_Add(expr):
            mm=0
            for i in fpoly(expr,'list'):
                mm=mm+powsimp(i, combine='base', force=True)
            return mm
        if Is_Number(expr):
            
            return ksym
        
        
        if fraction(expr)[0]==1:
             
            p1,p2==fraction(expr)
             
             
            return powsimp(p1, combine='base', force=True)/powsimp(p2, combine='base', force=True)
             
        
        elif Is_Mul(ksym):
            mm=fpoly(ksym,'list')
            kres=1
            for i in mm:
                kres=kres*powsimp(p1, combine='base', force=True)
             
            return kres
             
        elif type(ksym)==Pow:
            mm=fpoly(ksym,'list')
            val=mm[0]
            vexp=mm[1]
            if type(vexp)==Pow:

                return ksym
            elif type(vexp)==Mul:
                mm2=fpoly(vexp,'list')
                if len(mm2)==2:
                    p1=mm2[0]
                    p2=mm2[1]
                    if op!='':
                        kres=(val**p2)
                        kres=kres**p1
                    else:    
                        kres=(val**p1)
                        kres=kres**p2
                    return kres
                
            else:
                return ksym
             

             
        elif Is_Add(ksym):
            mm=fpoly(ksym,'list')
            mmr=0
            for i in mm:
                mmr+=powexpand(i,op=op)
            return mmr
        else:
            return ksym

def div2mulexp(ksym):
    '''
        input ((a/b)**x   ---->   return(a**x)*(b**(-x))
         
    '''
    if Is_Div(ksym):
        p1=numer(ksym)
        p2=denom(ksym)
        kres=p1*simplify(p2**(-1))
        return kres
    if Is_Pow(ksym):
        mm=fpoly(ksym,'list')
        vald=mm[0]
        vale=mm[1]
        if denom(vald)!='1':
            p1=numer(vald)
            p2=denom(vald)
            kres=(p1**vale)*(p2**(-1*vale))
            return kres
        else:
            return ksym
    if Is_Add(ksym):
        kres=0
        mm=fpoly(ksym,'list')
        for i in mm:
            kres+=div2mulexp(i)
        return kres
    if Is_Mul(ksym):
        kres=1
        mm=fpoly(ksym,'list')
        for i in mm:
            kres=kres*div2mulexp(i)
        return kres
    else:
        return ksym
def powsimpnumber(expr):
    if Is_MulPow(expr):
        nvar=list(expr.atoms(Integer))
        svar=[str(data) for data in nvar]
        sexpr=str(expr)
        v1,v2,v3,v4=symbols('v1 v2 v3 v4')
        rvar=[v1,v2,v3,v4]
        lvar=['v1','v2','v3','v4']
        cc=0
        qq=len(svar)
        for i in range(qq):
            sexpr=sexpr.replace(svar[cc],lvar[cc])
        nexpr=parse_expr(sexpr,evaluate=False)
        kres=powsimp(nexpr,force=True)
        sexpr=str(kres)
        for i in range(qq):
            sexpr=sexpr.replace(lvar[cc],svar[cc])
        return parse_expr(sexpr,evaluate=False)
    else:
        return expr   
def joinbase(expr):
    if Is_symbols(expr):
        return(expr)
    elif Is_Number(expr):
        return expr
    elif Is_Add(expr):
        mm=expr.args
        kres=joinbase(mm[0])
        for data in mm[1::]:  
            kres=kres+joinbase(data)
        return kres
    elif Is_Div(expr):
 
        p1,p2=fraction(expr)
 
        if Is_Pow(p1) and Is_Pow(p2) and getexpo(p1)==getexpo(p2):
 
            P1=getbase(p1)
            P2=getbase(p2)
            return spow(P1/P2,getexpo(p1))
        else:
 
            P1=joinbase(p1)
 
            P2=joinbase(p2)
            return sdiv(P1,P2)
    elif Is_Mul(expr):
        try:
            return powsimpnumber(expr)
        except:     
            return powsimp( expr ,force=True)
  
    elif Is_Pow(expr):
        bb=getbase(expr)
        ee=getexpo(expr)
        bb=powsimp(factor(bb),force=True)     
        return spow(bb,ee)    
    else:
        return expr

def disjoinbase(expr):
    if Is_Add(expr):
        mm=expr.args
        kres=0
        for data in mm:  
            kres=kres+disjoinbase(data)
        return kres
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        P1=disjoinbase(p1)
        P2=disjoinbase(p2)
        return Sdiv(P1,P2)
    elif Is_Mul(expr):
        mm=expr.args
        kres=disjoinbase(mm[0])
        for data in mm[1::]:  
            kres=kres*disjoinbase(data)
        return kres 
    elif Is_Pow(expr):
        bb=getbase(expr)
        ee=getexpo(expr)
        if Is_Mul(bb):
            mm=bb.args
            kres=mm[0]**ee
            for data in mm[1::]:
                kres=kres*data**ee
            return kres
        if Is_Add(ee):
            mm=ee.args
            kres=bb**mm[0]
            for data in mm[1::]:
                kres=kres*bb**data
            return kres       
        return expr    
    else:
        return expr

   
def joinbasemul(expr):
    if type(expr)==Mul:
        factores = expr.as_ordered_factors()  # Obtenemos los factores del producto
        potencias = {}  # Diccionario para agrupar bases por exponentes
        
        for factor in factores:
                base, exp =  powdenest(factor).as_base_exp()  # Obtenemos base y exponente
                
                # Si el exponente es negativo, tomamos el inverso de la base y el exponente positivo
                if signo(exp)==-1:
                    base = 1 / base
                    exp = -exp
                
                if exp in potencias:
                    potencias[exp].append(base)  # Si el exponente ya existe, añadimos la base
                else:
                    potencias[exp] = [base] 
        vece=[]
        vecb=[]
        for exp, bases in potencias.items():
            vece.append(exp)
            kres=1
            for data in bases:
                if kres==1:
                    kres=data
                else:
                    kres=kres*data
            vecb.append(kres) 
        kres=1
        for bb,ee in zip(vecb,vece):
            if ee==1:
                kresa=bb
            else:    
                kresa=spow(bb,ee)
            if kres==1:
                kres=kresa
            else:
                kres= kres*kresa  
        return kres
    else:
        return expr


def expopos(expr):
    if Is_e(expr):
        mm=expr.args
        sexp=mm[0]
        if signo(sexp)==-1:
            kres=parse_expr(str(1)+'/'+str(exp(-1*sexp)),evaluate=False)
            obj=kres.args
            return obj[1]
            
    bb,ee,rr=dataPow(expr)
    if Is_Add(ee):
        p1=1
        p2=1
        for i in ee.args:
            if signo(i)==-1:
                p2=p2*Sqrt(bb,rr,-1*i)
            else:
                p1=p1*Sqrt(bb,rr,i)
        if p2!=1:
            return parse_expr(str(p1)+'/('+str(p2)+')',evaluate=False)
        else:
            return p1
    elif Is_Mul(ee):
        if signo(ee)==-1:
            kres=parse_expr(str(1)+'/'+str(Sqrt(bb,rr,-1*ee)),evaluate=False)
            obj=kres.args
            
            return obj[1]
             
        else:
            return Sqrt(bb,rr,ee)
                            
    else:
        return expr
def sfactor(expr):
    if Is_Root(expr):
        bb,ee,rr=dataPow(expr)
        if rr==2 and ee==1:
            return parse_expr('sqrt('+str(sfactor(bb))+')',evaluate=False)
        elif rr==2 and ee!=1:
            return parse_expr('(sqrt('+str(sfactor(bb))+'))**'+str(ee),evaluate=False)
    if Is_Add(expr):
        return cleansumvec(expr)
   
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        p1=sfactor(p1)
        p2=sfactor(p2)
        return Div(p1/p2)
    elif Is_Mul(expr):
     
        cc=1
        for i in expr.args:
            if cc==1: 
                kres= str(sfactor(i))
                cc=cc+1
            else:
                kres= kres+'*'+str(sfactor(i))
            
        return   parse_expr(kres,evaluate=False)   
    elif Is_Root(expr):
        BB=insideroot(expr)
        bb=getbase(BB)
        ee=getexpo(BB)
        rr=getroot(expr)
        nbb=sfactor(bb) 
        try:
            return Sqrt(nbb,rr,ee)
        except:
            return expr
    else:
        return expr 


        
        
def positivexpo(kres):
    return joinexpofrac(kres)
     
def cleanMul(expr):
    mm=expr.args
    if len(mm)==2 and mm[0]==1:
        return mm[1]
    else:
        return expr
def cleansumvec(expr):
    vec=expr.args
    k=1
    for i in vec:
        k=k*denom(i)
    kk=sum(Array([simplify(i*k)  for i in vec]))
    return parse_expr('('+str(kk)+')'+'/'+str(k),evaluate=False)
            
def Is_e(ksym):
    if type(ksym)==exp:
        return True
    else:
        return False

        
def separebase(expr): 
    '''
    if expr= x**(a*z) expandbasepow(expr)= x**a*x**z
    '''    
    if Is_Add(expr):
        kres=0
        for i in fpoly(expr,'list'):
            kres=kres+separebase(i)
        return kres    
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        P1=separebase(p1)
        P2=separebase(p2)
        return cfrac(P1,P2)
         
    elif Is_Mul(expr):
        kres=1
        for i in expr.args:
            if Is_Pow(i):
                kres=kres*separebase(i)
            else:
                kres=kres*i
        return kres
    elif Is_Pow(expr):
        ee=getexpo(expr)
        bb=getbase(expr)
        kk=1
        mm=list(ee.args)
        if mm==[]:
            mm=[ee]
        for i in mm:
            kk=kk*ppow(bb,i)
            #kk=Mul(kk,ppow(bb,i),evaluate=False)
         
        return kk
    else:
        return expr
def sum2mulexpo(expr):
    if Is_Add(expr):
        mm=0
        for i in fpoly(expr,'list'):
            mm=mm+sum2mulexpo(i)
        return mm    
     
    elif Is_Mul(expr):
         
        mm=1
        for i in fpoly(expr,'list'):
             
            mm=mm*sum2mulexpo(i)
        return mm
    elif Is_Pow(expr):
        ee=getexpo(expr)
        bb=getbase(expr)
        if Is_Add(ee):
            kres=1
            for i in fpoly(ee,'list'):
                kres=kres*bb**i
            return kres    
        else:
            return expr
    else:
        return expr 

def nrsimply(expr):
    if 'sqrt' in str(expr):
        k=mathinsidepar(str(expr),'sqrt')
        if Is_Pow(k):
            bb=getbase(k)
            ee=getexpo(k)
            if ee%2==0:
                oldsexpr=str(sqrt(k))
                sexpr=str(expr)
                nsexpr=sexpr.replace(oldsexpr,str(bb**(ee/2)))
                return parse_expr(nsexpr)
            else:
                return expr
        else:
            return expr
    else:
        return expr
        
 
    


def krsimplify(expr):
    if Is_Root(expr):
        r=symbols('r')
        rr=getroot(expr)
        exp2=insideroot(expr)
        nexpr=ppow(exp2,cfrac(1,r))
        nexpr=mulexpo(nexpr)
        nexpr=simplifyexpo(nexpr)
        nexpr=subsubs(nexpr,r,rr)
        return nexpr
    else:
        return expr
        

		
	
def swapPow2Pow(expr):
	 
    if Is_Pow(expr):
        bb=getbase(expr)
        ee=getexpo(expr)
         
        if Is_Pow(bb):
             
            bb2=UnevaluatedExpr(getbase(bb))
            ee2=UnevaluatedExpr(getexpo(bb))
            return ppow(bb2**ee,ee2)

def checkrpow(expr,rval):
    if Is_Pow(expr):
        bb=getbase(expr)
        ee=getexpo(expr)
        if rval==ee:
            return bb
        else:
            nee=simplify(frs(ee,rval))
            return ppow(bb,nee)
    elif Is_Number(expr):
        expr=simplify(ppow(expr,1,rval))
        return expr
    else:
        return ppow(expr,frs(1,rval))
        
        

'''
def joinexpo(expr,op=''):
    if Is_Add(expr):
        kres=0
        for i in fpoly(expr,'list'):
            kres=kres+joinexpoA(i,op=op)
        return kres
    elif Is_Mul(expr) or Is_Div(expr):
        return joinexpoA(expr,op=op)
    else:
        return expr
'''        
        
def joinexpoA(expr,op=''):
    '''
    input x**a*z**a/y**a
    return (x*z/y)**a
    '''
    vecexpo=[]
    vecbases=[]
    for i in expr.args:
        if Is_Pow(i):
            bb=getbase(i)
            ee=getexpo(i)
            if ee in vecexpo:
                kpos=getposinvec(vecexpo,ee)
                pbase=vecbases[kpos]
                pbase.append(bb)
                vecbases[kpos]=pbase
            elif -1*ee in vecexpo:
                kpos=getposinvec(vecexpo,-1*ee)
                pbase=vecbases[kpos]
                pbase.append(1/bb)
                vecbases[kpos]=pbase

            else:
                vecexpo.append(ee)
                vecbases.append([bb])
    mres=1
    for i,j in zip(vecexpo,vecbases):
        kres=1
        for ii in j:
            if op=='simplify':
                kres=simplify(expand(kres*ii))
            else:
                kres=kres*ii
        kres=kres**i
        mres=mres*kres
    return mres

def getposinvec(vec,expr):
    pp=-1
    cc=0
    qq=len(vec)
    for i in range(qq):
        expr1=vec[i]
        if str(expr1)==str(expr):
            return cc
        cc=cc+1
    return -1  



def getrealexpo(expr):
     
    rr=getroot(expr)
    expr=mulexpo(expr)
    if rr==1:
        return getexpo(expr)
    else:
        ee=factor(getexpo(mulexpo(expr)))         
        if Is_Div(factor(ee)):
            p1,p2=fraction(factor(ee))
            if rr==p2:
                return p1
            else:
                return ee
        else:
            return ee

def realdataexpo(expr):
    if Is_PowRoot(expr):
        p1,p2=expr.args
        if Is_Pow(p1):
            bb,ee=p1.args
            r1,r2=fraction(p2)
            rr=r2
             
        else:
            bb=p1
            ee=p2
            rr=1
        
    else :
        bb=getbase(expr)
        ee=getexpo(expr)
        rr=1
    return bb,ee,rr

def fight_r_e(ee,rr):
    if rr==ee:
        return 1,1
    else:
        dd=simplify(factor(ee)/factor(rr))
        p1,p2=fraction(dd)
        return p1,p2

def rfactor(expr):
    if Is_Add(expr):
        kres=0
        for i in expr.args:
            kres=kres+rfactor(i)
        return kres   
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        p1=rfactor(p1)
        p2=rfactor(p2)
        return Div(p1,p2)
    elif Is_Mul(expr):
        kres=1
        for i in expr.args:
             
            kres=kres*rfactor(i)
        return kres    
    elif Is_Root(expr):
        BB=insideroot(expr)
        rr=getroot(expr) 
        BB=factor(BB)
        if denom(simplify(getexpo(BB)/rr))==1:
            return spow(getbase(BB),simplify(getexpo(BB)/rr)) 
        else:
            return rpow(BB,rr)     
 
    else:
        return expr
def rexpand(expr):
    if Is_Add(expr):
        kres=0
        for i in expr.args:
            kres=kres+rexpand(i)
        return kres   
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        p1=rexpand(p1)
        p2=rexpand(p2)
        return Div(p1,p2)
    elif Is_Mul(expr):
        kres=1
        for i in expr.args:
             
            kres=kres*rexpand(i)
        return kres    
    elif Is_Root(expr):
        BB=insideroot(expr)
        rr=getroot(expr) 
        BB=expand(BB)
        return rpow(BB,rr)     
 
    else:
        return expr
def killpow(expr):
    if Is_PowRoot(expr):
        try:
            bb,ee,rr=realdataexpo(expr)
            if ee==rr:
                return bb
                
            ee2=simplify(ee/rr)
            return bb**ee2

        except:
            return expr
    else:
        return expr

def rexpsimplify(expr):
    if Is_ExpRoot(expr):
        r=getroot(expr)
        b=insideroot(expr)
        ee=b.args[0]
        nexpr=simplify(cfrac(ee,r))
        if denom(nexpr)==1:
            return exp(nexpr)
        else:
            nn=numer(nexpr)
            dd=denom(nexpr)
            return rpow(exp(nn),dd)
         
    else:
        return expr
        
def simplersimplify(expr):
    bb,ee,rr=dataPow(expr)
    if Is_Pow(bb):
        return rsimplify(expr)
    elif Is_Div(bb):
         
        p1=simplersimplify(Sqrt(numer(bb),rr))
        p2=simplersimplify(Sqrt(denom(bb),rr))
        return sDiv(p1,p2)    
    elif  Is_Mul(bb):
        mm=bb.args
        vec=[simplersimplify(Sqrt(i,rr)) for i in mm]
        kres=list2mul(vec)
        return kres
    
    else:
        return expr        

        
def viewroot(expr):
    return rootview(expr)
def rootview(expr):
    if Is_Add(expr):
        kres=0
        for i in expr.args:
            kres=kres+rootview(i)
        return kres   
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        p1=rootview(p1)
        p2=rootview(p2)
        return Div(p1,p2)
    elif Is_Mul(expr):
        kres=1
        for i in expr.args:
             
            kres=kres*rootview(i)
        return kres    
    elif Is_Number(getexpo(expr)) and Is_Div(getexpo(expr)) and Is_Integer(denom(getexpo(expr))):
        bb=getbase(expr)
        ee=getexpo(expr)
        nee=numer(ee)
        nrr=denom(ee)
        B1=spow(bb,nee)
        return (B1)**cfrac(1,nrr)
    else:
        return expr        
      


        



#######################
########  used to rsimplify
################

def Is_RootMonoPow(expr):
    done=False
    if Is_Root(expr):
        if Is_Pow(insideroot(expr)):
            return True
    return False
def preSqrt(expr):
    bb,ee,rr= dataPow(expr)
    if type(ee)==Float:
        nee=Rational(ee).limit_denominator(10**3)
        knumer=numer(nee)
        kdenom=denom(nee)
        rr=rr*kdenom
        ee=knumer
    if Is_Pow(bb):
        bb2,ee2,rr2=preSqrt(bb)
        bb=bb2**cfrac(ee2,rr2)
    return bb,ee,rr   
def prePow(expr):
    bb,ee,rr=preSqrt(expr)
    if rr==1 and ee==1:
        return bb
    elif rr==1:
        return bb**ee
    elif ee==1:
        return bb**sDiv(1,rr)
    else:
        return bb**sDiv(ee,rr)

def Sqrt(*args):
    kres=Sqrt2(*args)
    return prePow(kres)
        
def Sqrt2(*args,evaluate=True):
    if len(args)==1 and Is_Add(args[0]):
        kres=0
        mm=args[0].args
        for i in mm:
            kres=kres+Sqrt2(i)
        return kres
    
    if len(args)==1 and Is_Div(args[0]):
        P=args[0]
        p1=prePow(numer(P))
        p2=prePow(denom(P))
        P1=Sqrt2(p1)
        P2=Sqrt2(p2)
        if Is_Root(P1) and Is_Root(P2) and (getroot(P1)==getroot(P2)):
            bb1,ee1,rr1=preSqrt(P1)
            bb2,ee2,rr2=preSqrt(p2)
            return Sqrt2(sDiv(bb1**ee1,bb2**ee2),rr1)
        else:
            return sDiv(P1,P2)

    if len(args)==1 and Is_Mul(args[0]):
        cc=1
        mm=args[0].args
        for i in mm:
            if cc==1:
                kres=Sqrt2(i)
                cc=cc+1
            else:
                kres=kres*Sqrt2(i)
        return kres
        
    if len(args)==1 and Is_Pow(args[0]):
        bb,ee,rr=preSqrt(args[0])
        bb=Sqrt2(bb) 
        return Sqrt2(bb,ee,rr)
    
        
    if len(args)==2:
        bb=args[0]
        rr=args[1]
        ee=1
    elif len(args)==3:
        bb=args[0]
        ee=args[1]
        rr=args[2]
    else:
        bb,ee,rr=preSqrt(args[0])
        bb=prePow(bb)
        
    if not evaluate:
        bb=prePow(bb)
        bb=parse_expr(str(bb),evaluate=False)
    if rr==1:
        bb=prePow(bb)
        return bb**ee
    elif rr==2:
        if ee==1:
            return parse_expr('sqrt('+str(bb)+')',evaluate=False)
        else:
            return parse_expr('sqrt(('+str(bb)+')**'+str(ee)+')',evaluate=False)
    else:
        bb=prePow(bb)
        if Is_Number(rr):
             
            sexpr='\\sqrt['+str(latex(rr))+']'
            if ee==1:
                sexpr=sexpr+'{'+str(latex(bb))+'}'
            else:    
                sexpr=sexpr+'{'+str(latex(bb))+'^{'+str(latex(ee))+'}}'
            return latex2sympy(sexpr)
        else:
            return  bb**sDiv(ee,rr) 

            
def Sqrt1(val,sq,ee=1,evaluate=True):
    if not evaluate:
        val=parse_expr(str(val),evaluate=False)
    if sq==1:
        return val**ee
    elif sq==2:
        if ee==1:
            return parse_expr('sqrt('+str(val)+')',evaluate=False)
        else:
            return sparse_expr('sqrt('+str(val)+'**'+str(ee)+')',evaluate=False)
    else:
        if Is_Number(sq):
            sexpr='\\sqrt['+str(latex(sq))+']'
            if ee==1:
                sexpr=sexpr+'{'+str(latex(val))+'}'
            else:    
                sexpr=sexpr+'{'+str(latex(val))+'^{'+str(latex(ee))+'}}'
            return latex2sympy(sexpr)
        else:
            return  val**cfrac(ee,sq)

def fracmixta(expr):
    expr=simplify(factor(expr))
    if denom(expr)!=1:
        p1=numer(expr)
        p2=denom(expr)
        P1=1
        K1=1
        if Is_Number(p1):
            P1=p1
            
        elif Is_Mul(p1):
            for i in p1.args:
                if Is_Number(i):
                    P1=i
                else:
                    K1=K1*i
        else:
            return (0,0,expr)
        P2=1
        K2=1
        if Is_Number(p2):
            P2=p2
        elif Is_Mul(p2):
            for i in p2.args:
                if Is_Number(i):
                    P2=i
                else:
                    K2=K2*i
        else:
            return (0,0,expr)
            
        kres=P1%P2
        return (P1-kres)/P2,kres*K1/K2,P2



    else:
        return 0,expr,1
def rsimplifyRPow(expr,*args):
    if Is_RootMonoPow(expr):
        base,ee,rr=partPow(expr)
        p1,p2,p3=fracmixta(ee/rr)
        if p1!=0:
            P3=Sqrt(base,p3,ee=p2)
            if 'parts' in args:
                return base,P3
            else:    
                return  Mul(base**p1,P3,evaluate=False) 
        else:
            if 'parts' in args:
                return 1,expr
            else:    
                return  expr
 
    if 'parts' in args:
        return 1,expr
    else:    
        return  expr

def mono_rsimplify(expr): # version 1
    bb,ee,rr=dataPow(expr)
    nexpr=bb**(simplify(Div(ee,rr)))
    bb,ee,rr=dataPow(nexpr)
    if  rr==1:
        return nexpr
    elif ee==rr:
        return bb
    elif ee==1 and rr!=1:
        if Is_Div(bb):
            p1=numer(bb)
            p1=mono_rsimplify(Sqrt(p1,rr))
            p2=denom(bb)
            p2=mono_rsimplify(Sqrt(p2,rr))
            return Div(p1,p2)
        elif Is_Mul(bb):
            kres=1
            for i in bb.args:
                kres=kres*mono_rsimplify(rpow(i,rr))
            return kres
        else:
            try:
                eres=ee%rr
                esal=ee-eres
                nee=esal/rr
                kres=Mul(bb**nee,Sqrt(bb,rr),evaluate=False)
                return kres
            except:
                return Sqrt(bb,rr,ee)
    elif ee==1 and rr==1:
        return expr 
        
   
    else:
        try:
            eres=ee%rr
            esal=ee-eres
            nee=esal/rr
            kres=Mul(bb**nee,Sqrt(bb,rr),evaluate=False)
            return kres
        except:
            return Sqrt(bb,rr,ee)
def mono_rsimplify2(expr):# version 1
    bb,ee,rr=dataPow(expr)
    nexpr=bb**(simplify(Div(ee,rr)))
    bb,ee,rr=dataPow(nexpr)
    if rr==1:
        if ee==1:
            return bb
        else:
            return bb**ee
    elif ee==1:
        if rr==1:
            return bb
        else:
            try:
                return Sqrt(bb,rr)
            except:
                try:
                    return mulexpo(disjoinexpo(bb**ee))
                except:
                    return expr
                    
    else:        
        try:
            eres=ee%rr
            esal=ee-eres
            nee=esal/rr
            kres=Mul(bb**nee,Sqrt(bb,rr),evaluate=False)
            return kres
        except:
            return Sqrt(bb,rr,ee) 


def dataPow(expr):
    if Is_e(expr):
        rr=getroot(expr)
        bb=exp(1)
        ee=expr.args[0]
    else:
        bb,ee,rr=partPow(expr)
        
    bb,ee,rr=unisymbols(bb),unisymbols(ee),unisymbols(rr)
    if denom(ee)!=1:
        rr=rr*denom(ee)
        ee=numer(ee)
    if denom(rr)!=1:
        ee=ee*denom(rr)
        rr=numer(rr)
    return bb,ee,rr
def monofactor(expr):
    if Is_Add(expr):
        return 1
    elif Is_Div(expr):
        return sdiv(monofactor(numer(expr)),monofactor(denom(expr)))
    elif Is_Mul(expr):
        kres=1
        for data in expr.args:
            kres=kres*monofactor(data)
        return kres
    elif Is_Pow(expr):
        ee=getexpo(expr)
        bb=getbase(expr)
        if Is_Add(bb):
            return 1
        else:
            return expr
    else:
        return expr

def simplifymonomies(expr):
    mm=monofactor(expr)
    kk=simplify(disjoinexpo(mm))
    mm2=simplify(expr/mm)
    return mm2*kk
def factorwithexp(expr,ee):
    if Is_Mul(expr) or Is_Pow(expr):
        if Is_Pow(expr) and getexpo(expr)==ee:
            return getbase(expr),1
        elif Is_Mul(expr):
            for data in expr.args:
                fact,resto= factorwithexp(data,ee)
                if fact!=1:
                    return fact,expr/data
            return 1,expr
        else:
            return 1,expr
    else:
        return 1,expr
        
def joinexpofrac(kres):
    if type(kres)==log:
        (kman,)=kres.args
        return log(joinexpofrac(kman))
        
         
    
    if type(kres)==Pow and type(getbase(kres))==Add:
        return joinexpofrac(getbase(kres))**getexpo(kres)
    
    if type(kres)==Mul and type(numer(kres))==Pow and type(denom(kres))==Pow and getexpo(numer(kres))==getexpo(denom(kres)):
        return spow(sdiv(getbase(numer(kres)),getbase(denom(kres))),simplify(getexpo(numer(kres))))
    if type(kres)==Mul and kres.args[0]==-1 and numer(kres)!=-1:
        kres=-1*kres
        kres=joinexpofrac(kres)
        return -kres
    
    if Is_Number(kres):
        return kres
    elif numer(kres)==-1 and signo(getexpo(denom(kres))) ==1:
        return sdiv(-1,denom(kres))
    if type(kres)==Pow:
         
        bb,ee=getbase(kres),getexpo(kres)
        if signo(ee)==-1:
            mm=bb**simplify(-ee)
            return sinversa(mm)
        else:
            return kres
 
    elif numer(kres)==1 and signo(getexpo(denom(kres))) ==1:
        return sinversa(1,denom(expr))
    elif numer(kres)==1 and signo(getexpo(denom(kres))) ==-1:
        return getebase(denom(kres))**-geteexpo(denom(expr)) 
 
    elif Is_Add(kres):
         
        kres2=0
        done=True
        for data in kres.args:
            if done:
                kres2=joinexpofrac(data)
                done=False
            else:
                kres2=kres2+joinexpofrac(data)
            
        return kres2 
    
    elif Is_Div(kres)and (signo(numer(kres)))**2!=1:
        done=False
        if kres.args[0]==-1:
            kres=kres*-1
            done=True
        pdd=denom(kres)
        pnn=numer(kres)
        veri1=Is_Mul(pnn) or Is_Pow(pnn)
        veri2=Is_Mul(pdd) or Is_Pow(pdd)
        if veri1 and veri2:
            vecb1=[]
            vecb2=[]
            vece1=[]
            if Is_Mul(pnn):
                vecp=pnn.args
            if Is_Pow(pnn):
                vecp=[pnn]
            for data in vecp:
                if Is_Pow(data):
                    bb,ee=getbase(data),getexpo(data)
                    kfac,kee=factorwithexp(pdd,ee)
                    if kfac!=1:
                        vecb1.append(bb)
                        vece1.append(ee)
                        vecb2.append(getbase(kfac))
            if len(vecb1)>0:            
                for bb1,bb2,ee1 in zip(vecb1,vecb2,vece1):
                    kres=kres/(bb1**simplify(ee1))
                for bb1,bb2,ee1 in zip(vecb1,vecb2,vece1):
                    kres=kres*(bb2**simplify(ee1))
                for bb1,bb2,ee1 in zip(vecb1,vecb2,vece1):
                    kres=kres*(sdiv(bb1,bb2)**simplify(ee1))
            if done:
                return -kres
            else:
                return kres
        else:
            if thereis_Div(pnn):
                pnn=joinexpofrac(pnn)
            if thereis_Div(pdd):
                pdd=joinexpofrac(pdd)
             
            if done: 
                return -Div(pnn,pdd)
            else:
                return Div(pnn,pdd)    
                    
    elif type(kres)==Mul and kres.args[0]==-1:
       
        kres=-1*kres
        return -joinexpofrac(kres)
        
    elif type(kres)==Mul:
        expr=kres
        mm=expr.args
        pnn=1
        pdd=1
        for data in mm:
            if type(data)==Pow:
                bb=getbase(data)
                ee=getexpo(data)
                if signo(ee)==-1:
     
                    pdd=pdd*Pow(bb,simplify(-ee))
                else:
                    pnn=smul(pnn,data)
            else:
                pnn=pnn*data
        if pdd!=1:
            kres=sdiv(pnn,pdd)
        else:
            kres=expr
        return kres
 
    else:
        return kres
        
from sympy import Expr, Add, Mul, Pow, fraction, factor
from typing import Callable

def expression_walker(
    expr: Expr,
    handle_add: Callable,
    handle_mul: Callable,
    handle_pow: Callable,
    handle_div: Callable,
    handle_base: Callable
) -> Expr:
    """
    Recorre recursivamente una expresión y aplica operaciones personalizadas a cada tipo de nodo.
    
    Parámetros:
        expr: Expresión SymPy a procesar.
        handle_add: Función que procesa nodos Add.
        handle_mul: Función que procesa nodos Mul.
        handle_pow: Función que procesa nodos Pow.
        handle_div: Función que procesa divisiones.
        handle_base: Función que procesa hojas (números, símbolos, etc.).
    """
    if expr.is_Add:
        return handle_add([expression_walker(arg, handle_add, handle_mul, handle_pow, handle_div, handle_base) for arg in expr.args])
    elif expr.is_Mul:
        return handle_mul([expression_walker(arg, handle_add, handle_mul, handle_pow, handle_div, handle_base) for arg in expr.args])
    elif expr.is_Pow:
        base = expression_walker(expr.base, handle_add, handle_mul, handle_pow, handle_div, handle_base)
        expo = expression_walker(expr.exp, handle_add, handle_mul, handle_pow, handle_div, handle_base)
        return handle_pow(base, expo)
    elif expr.is_Rational and not expr.is_Integer:
        num, den = fraction(expr)
        return handle_div(
            expression_walker(num, handle_add, handle_mul, handle_pow, handle_div, handle_base),
            expression_walker(den, handle_add, handle_mul, handle_pow, handle_div, handle_base)
        )
    else:
        return handle_base(expr)

def expofactor(expr: Expr) -> Expr:
    def handle_add(args): return Add(*args)
    def handle_mul(args): return Mul(*args)
    def handle_pow(base, expo): return Pow(base, factor(expo))  # <-- Factoriza el exponente
    def handle_div(num, den): return num / den
    def handle_base(x): return x    
    return expression_walker(expr, handle_add, handle_mul, handle_pow, handle_div, handle_base)
    
def basefactor(expr: Expr) -> Expr:
    def handle_add(args): return Add(*args)
    def handle_mul(args): return Mul(*args)
    def handle_pow(base, expo): return Pow(factor(base), expo)  # <-- Factoriza el exponente
    def handle_div(num, den): return num / den
    def handle_base(x): return x    
    return expression_walker(expr, handle_add, handle_mul, handle_pow, handle_div, handle_base)    

def expoexpand(expr: Expr) -> Expr:
    def handle_add(args): return Add(*args)
    def handle_mul(args): return Mul(*args)
    def handle_pow(base, expo): return Pow(base, expand(expo))  # <-- Factoriza el exponente
    def handle_div(num, den): return num / den
    def handle_base(x): return x    
    return expression_walker(expr, handle_add, handle_mul, handle_pow, handle_div, handle_base)
    
def expandexpo(expr: Expr) -> Expr:
    def handle_add(args): return Add(*args)
    def handle_mul(args): return Mul(*args)
    def handle_pow(base, expo): return Pow(base, expand(expo))  # <-- Factoriza el exponente
    def handle_div(num, den): return num / den
    def handle_base(x): return x    
    return expression_walker(expr, handle_add, handle_mul, handle_pow, handle_div, handle_base) 
    
def baseexpand(expr: Expr) -> Expr:
    def handle_add(args): return Add(*args)
    def handle_mul(args): return Mul(*args)
    def handle_pow(base, expo): return Pow(expand(base), expo)  # <-- Factoriza el exponente
    def handle_div(num, den): return num / den
    def handle_base(x): return x    
    return expression_walker(expr, handle_add, handle_mul, handle_pow, handle_div, handle_base) 
    
def expandbase(expr: Expr) -> Expr:
    def handle_add(args): return Add(*args)
    def handle_mul(args): return Mul(*args)
    def handle_pow(base, expo): return Pow(expand(base), expo)  # <-- Factoriza el exponente
    def handle_div(num, den): return num / den
    def handle_base(x): return x    
    return expression_walker(expr, handle_add, handle_mul, handle_pow, handle_div, handle_base) 
    
def exposimplify(expr: Expr) -> Expr:
    def handle_add(args): return Add(*args)
    def handle_mul(args): return Mul(*args)
    def handle_pow(base, expo): return Pow(base, simplify(expo))  # <-- Factoriza el exponente
    def handle_div(num, den): return num / den
    def handle_base(x): return x    
    return expression_walker(expr, handle_add, handle_mul, handle_pow, handle_div, handle_base)
    
def basesimplify(expr: Expr) -> Expr:
    def handle_add(args): return Add(*args)
    def handle_mul(args): return Mul(*args)
    def handle_pow(base, expo): return Pow(simplify(base), expo)  # <-- Factoriza el exponente
    def handle_div(num, den): return num / den
    def handle_base(x): return x    
    return expression_walker(expr, handle_add, handle_mul, handle_pow, handle_div, handle_base) 

def sinverse(expr):
    p1=numer(expr)
    p2=denom(expr)
    if p2==1:
        kres = sdiv(1,p1)
        return kres.args[1]
    else:
        return sdiv(p2,p1)    