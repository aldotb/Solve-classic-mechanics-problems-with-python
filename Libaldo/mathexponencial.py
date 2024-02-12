

from  libaldo_math2 import *
from  libaldo_algorith import *
from  mathbasic import *
from sympy import *
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
    if Is_Add(expr):
        mm=0
        for i in fpoly(expr,'list'):
            mm=mm+mulexpo(i)
        return mm    
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        P1=mulexpo(p1)
        P2=mulexpo(p2)
        sP='('+str(P1)+')/('+str(P2)+')'
        return mirror_parse_expr(sP)
    elif Is_Mul(expr):
        mm=1
        for i in fpoly(expr,'list'):
            mm=mm*mulexpo(i)
        return mm
    elif Is_Pow(expr):
        if Is_PowPow(expr):
            bb=getbase(expr)
            ee=getexpo(expr)
            bb2=getbase(bb)
            ee2=getexpo(bb)
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
    

def positiveexpo(expr):
    return positivexpo(disjoinbase(expr))
    
def positivexpo(expr):
    if Is_Mul(expr):
        return disjoinbase(expr)
    elif Is_Add(expr):
        mm=0
        for i in fpoly(expr,'list'):
            mm=mm+positivexpo(i)
        return mm
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        P1=positivexpo(p1)
        P2=positivexpo(p2)
        sP='('+str(P1)+')/('+str(P2)+')'
        return mirror_parse_expr(sP)
    elif Is_Mul(expr) and denom(expr)!=1:
         
        mm=1
        for i in fpoly(expr,'list'):
            try:
                mm=mm*positivexpo(i)
            except:
                mm=mm*i
        return mm 
    elif Is_Pow(expr):
            
            ee=getexpo(expr)
            bb=getbase(expr)
            if signo(ee)==-1:
                if Is_Div(bb):
                    p1,p2=feaction(bb)
                    p3=frs(p2,p1)
                    nee=-1*ee
                    return p3**nee
                else:    
                    return mirror_parse_expr('1/'+str(bb**(-1*ee)))
                 
            else:
                return expr
   
    else:
        return expr    
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
    return rsimplify3(expr)


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