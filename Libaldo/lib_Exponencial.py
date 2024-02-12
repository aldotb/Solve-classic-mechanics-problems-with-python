


from sympy import *
from lib_Mathematica import *
from lib_Mathbasic import *
from lib_Algorith import *
## exponenciales functions...

x,s,t,y=symbols('x s  t y')
# ussefull variables to use in partial fracction algorith
y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16=symbols('y1 y2  y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 y13 y14 y15 y16')
fabr=[y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16]

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
    expr2=eval(sexpr)
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


def getexpo(expr,op=''): #   return exponente in expr
    if Is_Pow(expr):
        mm=expr.args
        return mm[1]
        
    else:
        return 1
    
def getpow(expr): #   return exponente in expr
    '''
    expr = x**y,  x*y
    return y ,1
    '''
    if Is_Root(expr):
        kres=insideroot(expr)
        return getexpo(kres)
    elif Is_Pow(expr):
        ee=expr.args[1]
        p1,p2=fraction(ee)
        
        return p1
    elif Is_Root(expr):
        kres=insideroot(expr)
        return getexpo(kres)
    else:
        return 1        


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
        
# ***Multifunc to Add  Div and Mul****
# ------------------------------
def getroot(expr):
    if Is_Root(expr):
        mm=expr.args
        kres=mm[1]
        p1,p2=fraction(kres)
        return p2
    
    else:
        return 1


'''
def simplifyexpo(expr): # simplify each exponent  in expr
    if Is_Add(expr):
        mm=0
        for i in fpoly(expr,'list'):
            mm=mm+simplifyexpo(i)
        return mm
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        P1=simplifyexpo(p1)
        P2=simplifyexpo(p2)
        sP='('+str(P1)+')/('+str(P2)+')'
        return mirror_parse_expr(sP)   
    elif Is_Mul(expr):
        mm=1
        for i in fpoly(expr,'list'):
            mm=mm*simplifyexpo(i)
        return mm
    elif Is_Pow(expr):
        ee=getexpo(expr)
        bb=getbase(expr)
        ee=simplify(ee)
        return bb**ee
    else:
        return expr
'''
def simplifyexpo(expr): # simplify each exponent  in expr
    kres=corefunc(expr,'simplifyexpo')    
    if Is_Pow(kres):
        kres=getbase(kres)**simplify(getexpo(kres))
    return kres
'''
def expandexpo(expr): # simplify each exponent  in expr
    if Is_Add(expr):
        mm=0
        for i in fpoly(expr,'list'):
            mm=mm+expandexpo(i)
        return mm
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        P1=expandexpo(p1)
        P2=expandexpo(p2)
        sP='('+str(P1)+')/('+str(P2)+')'
         
        return mirror_parse_expr(sP)   
    elif Is_Mul(expr):
        mm=1
        for i in fpoly(expr,'list'):
            mm=mm*expandexpo(i)
        return mm
    elif Is_Pow(expr):
         
        bb=getbase(expr)
        ee=getexpo(expr)
        ee=expand(ee)
        return ppow(bb,ee)
        
    else:
        return expr
'''
def expandexpo(expr): # simplify each exponent  in expr
    kres=corefunc(expr,'expandexpo')
    if Is_Pow(kres):        
        kres=getbase(kres)**expand(getexpo(kres))
    return kres
def expandexpoexpand(expr): # simplify each exponent  in expr
    kres=corefunc(expr,'expandexpo')
    if Is_Pow(kres):
        ee=getexpo(kres)
        base=getbase(kres)
        if Is_Pow(ee):
            ee=expandexpo(ee)
            kres= bb**ee
    return kres
'''
def factorexpo(expr):
    if Is_Add(expr):
        mm=0
        for i in fpoly(expr,'list'):
            mm=mm+factorexpo(i)
        return mm    
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        P1=factorexpo(p1)
        P2=factorexpo(p2)
        sP='('+str(P1)+')/('+str(P2)+')'
        return mirror_parse_expr(sP)
    elif Is_Mul(expr):
        mm=1
        for i in fpoly(expr,'list'):
            mm=mm*factorexpo(i)
        return mm 
    elif Is_Pow(expr):
        ee=getexpo(expr)
        bb=getbase(expr)
        ee=factor(ee)
        return bb**ee
    else:
        return expr
'''
def factorexpo(expr): # simplify each exponent  in expr
    kres=corefunc(expr,'factorexpo')
    if Is_Pow(kres):
        kres=getbase(kres)**factor(getexpo(kres))    
    return kres
        
'''
def simplifybase(expr): # simplify each ecponet in expr
    if Is_Add(expr):
        mm=0
        for i in fpoly(expr,'list'):
            mm=mm+simplifybase(i)
        return mm    
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        P1=simplifybase(p1)
        P2=simplifybase(p2)
        sP='('+str(P1)+')/('+str(P2)+')'
        return mirror_parse_expr(sP)
    elif Is_Mul(expr):
        mm=1
        for i in fpoly(expr,'list'):
            mm=mm*simplifybase(i)
        return mm
    elif Is_Pow(expr):
        ee=getexpo(expr)
        bb=getbase(expr)
        bb=simplify(bb)
        return bb**ee
    else:
        return expr
'''
def simplifybase(expr): # simplify each exponent  in expr
    kres=corefunc(expr,'simplifybase')
    if Is_Pow(kres):        
        kres=getbase(kres)**simplify(getexpo(kres))
    return kres


'''
def expandbase(expr): # simplify each ecponet in expr
    if Is_Add(expr):
        mm=0
        for i in fpoly(expr,'list'):
            mm=mm+expandbase(i)
        return mm    
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        P1=expandbase(p1)
        P2=expandbase(p2)
        sP='('+str(P1)+')/('+str(P2)+')'
        return mirror_parse_expr(sP)
    elif Is_Mul(expr):
        mm=1
        for i in fpoly(expr,'list'):
            mm=mm*expandbase(i)
        return mm
    elif Is_Pow(expr):
        return expand_power_base(expr)
    else:
        return expr
'''
def expandbase(expr): # simplify each exponent  in expr
    kres=corefunc(expr,'expandbase')
    if Is_Pow(kres):        
        kres=expand(getbase(kres))**getexpo(kres)
    return kres        
'''
def factorbase(expr): # simplify each ecponet in expr
    if Is_Add(expr):
        mm=0
        for i in fpoly(expr,'list'):
            mm=mm+factorbase(i)
        return mm    
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        P1=factorbase(p1)
        P2=factorbase(p2)
        sP='('+str(P1)+')/('+str(P2)+')'
        return mirror_parse_expr(sP)
    elif Is_Mul(expr):
        mm=1
        for i in fpoly(expr,'list'):
            mm=mm*factorbase(i)
        return mm
    elif Is_Pow(expr):
        ee=getexpo(expr)
        bb=getbase(expr)
        bb=factor(bb)
        return bb**ee
    else:
        return expr
'''
def factorbase(expr): # simplify each exponent  in expr
    kres=corefunc(expr,'factorbase')
    if Is_Pow(kres):        
        kres=factor(getbase(kres))**getexpo(kres)
    return kres
 
def mulexpo(expr,force=False): # simplify each ecponet in expr
    if Is_Symbol(expr):
        return expr
    if Is_Number(expr):
        return expr
    if Is_Root(expr):
        rr=getroot(expr)
        bb=insideroot(expr)
        ee=getexpo(expr)
        if Is_e(bb):
            nee=bb.args[0]
            if rr==1:
                return  exp(simplify(ee*nee))
            elif rr==2:
                return  sqrt(exp(simplify(ee*nee)))
            else:
                return  Sqrt(exp(simplify(ee*nee)),rr)
                
            
            
    elif Is_Add(expr):
        mm=0
        for i in expr.args:
            mm=mm+mulexpo(i)
        return mm    
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        P1=mulexpo(p1)
        P2=mulexpo(p2)
        return unisymbols(cfrac(P1,P2))
    elif Is_Mul(expr):
        mm=1
        for i in expr.args:
            mm=mm*mulexpo(i)
        return mm
    elif Is_Pow(expr):
        if Is_PowPow(expr):
            bb,ee=getbase(expr),getexpo(expr)
            bb2,ee2=getbase(bb),getexpo(bb)
            return bb2**(ee*ee2)
        else:
            bb=getbase(expr)
            ee=getexpo(expr)
            if Is_Mul(bb):
                mm=bb.args
                kres=1
                for i in mm:
                    if Is_Pow(i):
                        bi=getbase(i)
                        ei=getexpo(i)
                        kres=kres*bi**(ei*ee)
                    else:
                        kres=kres*i**ee
                return kres       
            else:
                bb2=mulexpo(bb) 
                ee2=mulexpo(ee)
                return bb2**ee2
    else:
        return expr
  

 
def lexpand(expr,force=True): # simplify each exponent  in expr
     
    if Is_Add(expr):
        mm=0
        for i in expr.args:
            mm=mm+lexpand(i)
        expr=mm    
    expr=corefunc(expr,'lexpand')
    if Is_Log(expr):
        return expand_log(expr, force=True)     
           
    else:
        return expr
        
def lmul2lpow(expr):
    if Is_Add(expr):
        kres=0
        for i in fpoly(expr,'list'):
            kres=kres+lmul2pow(expr)
        return kres
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        p1=lmul2pow(p1)
        p2=lmul2pow(p2)
        return cfrac(p1,p2)
    elif Is_Mul(expr):
        
        if 'log' in str(expr):
            ee=1
            vlog=1
            for i in fpoly(expr,'list'):
            
                if type(i)==log:
                    vlog=i.args[0]
                else:
                    ee=ee*i
            vlog=vlog**ee        
            return log(vlog)       
    else:
        return expr
        
# def lfactor(expr): # expand logaritm  
    # if Is_Add(expr):
        # mm=0
        # for i in fpoly(expr,'list'):
            # mm=mm+lfactor(i)
        # return mm    
    # elif Is_Div(expr):
        # p1,p2=fraction(expr)
        # P1=lfactor(p1)
        # P2=lfactor(p2)
        # sP='('+str(P1)+')/('+str(P2)+')'
        # return mirror_parse_expr(sP)
    # elif Is_Mul(expr):
        # mm=1
        # for i in fpoly(expr,'list'):
            # mm=mm*lfactor(i)
        # return mm
       
    # elif Is_Log(expr):
        # return logcombine(expr, force=True) 
           
    # else:
        # return expr

def lfactor(expr):
    kres=corefunc(expr,'lexpand')
    if Is_Log(kres):
        return expand_log(expr, force=True)     
           
    else:
        return expr        
        
def lexponent(expr):
    if Is_Add(expr):
        mm=0
        for i in fpoly(expr,'list'):
            mm=mm+lexponent(i)
        return mm
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        P1=lexponent(p1)
        P2=lexponent(p2)
        sP='('+str(P1)+')/('+str(P2)+')'
        return mirror_parse_expr(sP)
    elif Is_Mul(expr):
        mm=1
        for i in fpoly(expr,'list'):
            mm=mm*lexponent(i)
        return mm    
    elif Is_Log(expr):
        expr2=expr.args[0]
        ee=getexpo(expr2)
        bb=getbase(expr2)
         
        return ee*log(bb)
    else:
        return expr        

def lcombine(expr):
    expr=logcombine(expr,force=True)
    return expr



 

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
        
def joinbase(expr): 
    '''
    if expr= x**(a+b)*x(z) joinbase(expr)= x**(a+b+z)
    '''    
    if Is_Add(expr):
        kres=0
        for i in fpoly(expr,'list'):
            kres=kres+joinbase(i)
        return kres    
    elif Is_Mul(expr) and denom(expr)!=1:
         
        return simplify(expr)
    elif Is_Pow(expr):
            bb=getbase(expr)
            ee=getexpo(expr)
            bb=joinbase(bb)
            ee=joinbase(ee)
            return ppow(bb,ee)
    elif Is_Mul(expr):
        kres=1
        mm=fpoly(expr,'list')
        bb=[]
        ee=[]
        for i in mm:
            bb1=getbase(i)
            ee1=getexpo(i)
            if bb1 not in bb:
                bb.append(bb1)
                ee.append(ee1)
            else:
                p1=bb.index(bb1)
                nexp=ee[p1]
                nexp=nexp+ee1
                ee[p1]=nexp
        kres=1
        for i,j in zip(bb,ee):
            kres=kres*ppow(i,j)
        if kres==expr:
            return powsimp(kres,force=True)
        else:    
            return kres
    else:
        return expr
        
def disjoinbase(expr):
    return separebase(expr)
    
def disjoinexpo(expr):
    kres=corefunc(expr,'disjoinexpo')
    if Is_Pow(expr):
        base=getbase(expr)
        ee=getexpo(expr)
        ee2=disjoinbase(ee)
        kres=base**ee2
    return kres 
    
    
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
def positivexpo(expr):
    if type(expr)==Symbol:
        return expr
    elif Is_Number(expr):
        return expr
        
    if Is_Root(expr):
         
        bb,ee,rr=dataPow(expr)
        if Is_Mul(bb):
            bb=cleanMul(positivexpo(bb))
        return Sqrt(positivexpo(bb),rr,ee)
    elif Is_e(expr):
         return cleanMul(expopos(expr))
         
    elif Is_Add(expr):
         
        mm=0
        for i in expr.args:
            mm=mm+cleanMul(positivexpo(i))
        return mm    
   
    elif Is_Div(expr):
         
        nn=numer(expr)
        dd=denom(expr)
         
        return parse_expr(str(nn)+'/('+str(dd)+')',evaluate=False)
         
    elif Is_Mul(expr):
 
        kres=1 
        for i in expr.args:
             kres=kres*cleanMul(positivexpo(i))
        return cleanMul(kres)      
 
    elif Is_Pow(expr):
        return expopos(expr)
   
    else:
        return expr 
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


def rsimplify(expr):
    if type(expr)==list:
        return [rsimplify9(i) for i in expr]
    elif type(expr)==Array:
        vec=list(expr)
        kres=[rsimplify9(i) for i in vec]
        return Array(kres)
    else:
        return rsimplify9(expr) 


def rsimplify3(expr): # simplify each ecponet in expr
    if Is_Add(expr):
        mm=0
        for i in fpoly(expr,'list'):
            mm=mm+rsimplify3(i)
        return mm    
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        P1=rsimplify3(p1)
        P2=rsimplify3(p2)
        sP='('+str(P1)+')/('+str(P2)+')'
        return mirror_parse_expr(sP)
    elif Is_Mul(expr):
        mm=1
        for i in fpoly(expr,'list'):
            mm=mm*rsimplify3(i)
        return mm
    elif Is_Pow(expr) and not Is_Root(expr):
        bb=getbase(expr)
        ee=getexpo(expr)
        bb=rsimplify3(bb)
        return bb**ee
         
    elif Is_Root(expr):
        rr=getroot(expr)
        mm=insideroot(expr)
        if Is_Pow(mm):
            bb=getbase(mm)
            ee=getexpo(mm)
            return(bb**(simplify(cfrac(ee,rr))))
        elif Is_Mul(mm):
            kres=1
            kl=mm.args
            for i in kl:
                bb=getbase(i)
                ee=getexpo(i)
                kres=kres*bb**(simplify(cfrac(ee,rr)))
            return kres
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
        
def rsimplify2(expr): # simplify each ecponet in expr
    if Is_Add(expr):
        mm=0
        for i in fpoly(expr,'list'):
            mm=mm+rsimplify2(i)
        return mm    
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        P1=rsimplify2(p1)
        P2=rsimplify2(p2)
        sP='('+str(P1)+')/('+str(P2)+')'
        return mirror_parse_expr(sP)
    elif Is_Mul(expr):
        mm=1
        for i in fpoly(expr,'list'):
            mm=mm*rsimplify2(i)
        return mm
    elif Is_Pow(expr):
        return krsimplify(expr)
         
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
        
        
def joinexpo(expr):
    if Is_Add(expr):
        kres=0
        for i in expr.args:
            if Is_MulPow(i):
                kres=kres+joinexpo(i)
            else:
                kres=kres+i
        return kres
    elif Is_Div(unisymbols(expr)):
        p1,p2=unisymbols(fraction(expr))
        if Is_MulPow(p1):
            p1=joinexpo(p1)
        if Is_MulPow(p2):
            p2=joinexpo(p2)
        
        kres=Div(p1,p2)
        return kres
    elif Is_MulPow(unisymbols(expr)):
        vecnor=[]
        vecpow=[]
        for i in expr.args:
            if type(i)==Pow:
                vecpow.append(i)
            else:
                vecnor.append(i)
        vecee=[]
        for i in vecpow:
            ee=getexpo(i)
            if ee not in vecee:
                vecee.append(ee)
        vecbases=[]
        for i in vecee:
            vecbase=[]
            for j in vecpow:
                base=getbase(j)
                ee=getexpo(j)
                if ee==i:
                    vecbase.append(base)
            vecbases.append(vecbase)
        kres=1
        vecpro=[]
        for i  in vecbases:
            kres=1
            for j in i:
                kres=kres*j
            vecpro.append(kres)
        vecexp=[]
        for i in vecpro:
            try:
                kres=expand(i)
            except:
                kres=i
            vecexp.append(kres)
        kres=1
        for i,j in zip( vecexp,vecee):
            kres=kres*(i**j)
        for i in vecnor:
            kres=kres*i 
        return kres
    else:
        return expr
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
        return factor(kres)   
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        p1=rfactor(p1)
        p2=rfactor(p2)
        return factor(p1/p2)
    elif Is_Mul(expr):
        kres=1
        for i in expr.args:
             
            kres=kres*rfactor(i)
        return factor(kres)    
    elif Is_Root(expr):
        BB=insideroot(expr)
        BB3=factor(BB) 
        if Is_Div(BB3):
             
            p1=factor(mulexpo(numer(factor(BB3))))
            p2=factor(mulexpo(denom(factor(BB3))))
            BB2=factor(p1/p2)
            return factor(expr.subs(BB,BB2))
        elif Is_Add(BB):
            BB2=factor(unisymbols(BB))
            return factor(expr.subs(BB,BB2))
        else:
            return expr
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
def rsimplify9(expr):
    if Is_Add(expr):
        kres=0
        for i in expr.args:
            kres=kres+rsimplify9(i)
        return kres 
    elif Is_Div(expr):
        p1=numer(expr)
        p1=rsimplify9(p1)
        p2=denom(expr)
        p2=rsimplify9(p2)
        return Div(p1,p2)
    elif Is_Mul(expr):
        mm=expr.args
        vec=[rsimplify9(i) for i in mm]
        kres=list2mul(vec)
        return kres
    elif Is_Root(expr):
        if (Is_Root(expr) and (Is_Div(insideroot(expr)) or Is_Mul(insideroot(expr)))):
            return simplersimplify(expr) 
        else:
            return mono_rsimplify(Sqrt(expr))
    else:
        return expr
def rsimplify8(expr):
    if Is_Add(expr):
        kres=0
        for i in expr.args:
            kres=kres+rsimplify8(i)
        return kres 
    elif Is_Div(expr):
        p1=numer(expr)
        p1=rsimplify8(p1)
        p2=denom(expr)
        p2=rsimplify8(p2)
        return Div(p1,p2)
    elif Is_Mul(expr):
        kres=1
        for i in expr.args:
            kres=kres*rsimplify8(i)
        return kres
    elif Is_Root(expr):
        if (Is_Root(expr) and (Is_Add(insideroot(expr)) or Is_Mul(insideroot(expr)))):
            kres=insideroot(expr)
            kres=rsimplify8(kres)
            rr=getroot(expr)
            return Sqrt(kres**sDiv(1,rr))
        else:
            return mono_rsimplify(Sqrt(expr))
    else:
        return expr        
def rsimplify7(expr):
    if Is_Add(expr):
        kres=0
        for i in expr.args:
            kres=kres+rsimplify7(i)
        return kres 
    elif Is_Div(expr):
        p1=numer(expr)
        p1=rsimplify7(p1)
        p2=denom(expr)
        p2=rsimplify7(p2)
        return Div(p1,p2)
    elif Is_Mul(expr):
        kres=1
        for i in expr.args:
            kres=kres*rsimplify7(i)
        return kres
    elif Is_Root(expr):
        return mono_rsimplify(expr)
    else:
        return expr

        
def rsimplify6(expr):
    if Is_Add(expr):
        kres=0
        for i in expr.args:
            kres=kres+rsimplify6(i)
        return kres 
    elif Is_Div(expr):
        p1=numer(expr)
        p1=rsimplify6(p1)
        p2=denom(expr)
        p2=rsimplify6(p2)
        return Div(p1,p2)
    elif Is_Mul(expr):
        kres=1
        for i in expr.args:
            kres=kres*rsimplify6(i)
        return kres
    elif Is_Root(expr):
        return mono_rsimplify(expr)
    else:
        return expr


def rsimplify5(expr):

    
    if Is_Root(expr) and not Is_Add(insideroot(expr)):
        rr=getroot(expr)
        bb=getbase(expr)
        P1=1
        P2=1
        mm=bb.args
        
        for i in mm:
            kres=expr
            for j in mm:
                if j!=i:
                    kres=kres.subs(j,1)
            exprr=kres
            res1,res2=rsimplifyRPow(exprr,'parts')
 
            if res1!=1:
                if P1==1:
                    P1=res1
                else:    
                    P1=Mul(P1,res1,evaluate=False) 
            if res2!=1:
                if P2==1:
                    P2=res2
                else:    
                    P2=Mul(P2,res2,evaluate=False)            
        if P1!=1:
            return Mul(rsimplify(P1),joinbase(P2),evaluate=False) 
        else:
            return rsimplify(P2)
        
    elif Is_ExpRoot(expr):
            return rexpsimplify(expr)
            
    elif Is_Add(expr):
        kres=0
        for i in expr.args:
            kres=kres+unisymbols(rsimplify5(i))
        return kres

            
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        p1=unisymbols(rsimplify5(p1))
        p2=unisymbols(rsimplify5(p2))
        return p1/p2
    elif Is_Mul(expr):
        kres=1
        for i in expr.args:
            kres=kres*unisymbols(rsimplify5(i))
        return kres
    elif Is_Root(expr):
        BB=insideroot(expr)
        rr=getroot(expr)
        BB3=factor(BB)
        if Is_Div(BB3):
            p1=unisymbols(factor(mulexpo(numer(factor(BB3)))))
            p2=unisymbols(factor(mulexpo(denom(factor(BB3)))))
            p1=simplify(rsimplify(rpow(p1,rr)))
            p2=simplify(rsimplify(rpow(p2,rr)))
            p1=killpow(p1)
            p2=killpow(p2)
            return p1/p2
        else:
            return killpow(expr)

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
            
