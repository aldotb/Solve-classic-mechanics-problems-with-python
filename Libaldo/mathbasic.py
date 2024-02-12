

from  libaldo_math2 import *
from  libaldo_algorith import *
from  mathexponencial  import *
from  mathexponencial  import getbase,getexpo
from sympy import *
## PARTIAL FRACTION...

x,s,t,y=symbols('x s  t y')
# ussefull variables to use in partial fracction algorith
y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16=symbols('y1 y2  y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 y13 y14 y15 y16')
fabr=[y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16]
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
def str2mathstr(sexpr):
    '''
        input str math expr
        return str math expr
        this f fix space between + and - symbols
        example:
        str2mathstr('x**2+z*(x+1)')
        return 'x**2 + z*(x + 1)'
    '''
    
    
    qq=len(sexpr)
    sres=''
    
    for i in range(qq):
        sval=sexpr[i]
        if i==0:
            sres=sres+sval
        elif sval=='+' or sval=='-':
            if sexpr[i-1] != ' ':
                sres=sres+' '
            sres=sres+sval    
            if i<qq:
                if sexpr[i+1] != ' ':
                    sres=sres+' '
        else :
            sres=sres+sval
    return sres



def complete2partfrac(vec):
    # complete missing monomies to run partial fracction algorith
    vec2=[]
    vecmono=fpoly(vec,'list')
    for i in vecmono:
        if Is_Pow(i):
            ee=getexpo(i)
            bb=getbase(i)
 
            for j in range(ee):
                vec2.append(bb**(j+1))
        else:
            vec2.append(i) 
    return vec2
    
    
def mdaughter(expr,vecf,gen=x): # used in partialfracction
    # create one degree bellow monomio formate to part fracct algorith
    # input monimie, vec with available variables, main variab default x
    # return daughter monimie, ejem 3*x+1..y1, 3*x*x+2*x+1..x*y2+y3...,
    dg=degree(expr)
    expr2=0
    for i in range(dg):
         
        expr2=expr2+gen**i*vecf[i]
    vecf=vecf[dg:-1]
    return expr2,vecf 

def vec_daughter(vecm,vecf,gen=x): # used in partialfracction
    #    retun vector whit respevctive nuerator formato por partf
    vec2=vecm
    kres=[]
    for j in range(len(vec2)):
        exp1,vecf=mdaughter(vec2[j],vecf,gen=gen)
        kres.append(exp1)
    return kres 
    
def partialfraction(*args):
    ff=args[0]
    if len(args)==2:
        return apart(ff,args[1])
    else:
        return apart(ff)
     
# def partialfraction(*args): # Return divition monimies in partial fraction 
 
    # y1,y2,y3,y4,y5,y6,y7,y8,y9,y10=symbols('y1 y2  y3 y4 y5 y6 y7 y8 y9 y10')
    # y11,y12,y13,y14,y15,y16=symbols('y11 y12 y13 y14 y15 y16')
    # fabr=[y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16]
    # expr=args[0]
    # if len(args)==2:
        # var=args[1]
    # else:
        # var=x
    # expr=factor(expr)
    # p1=numer(expr)
    # p2=denom(expr)
    # vec2=complete2partfrac(p2)
    # vec1=vec_daughter(vec2,fabr,gen=var)
    # P=0
    # for i,j in zip(vec1,vec2):
        # P=P+cfrac(i,j)
    # P2=factor(P)
    # P3=numer(P2)
    # P4=expand(P3)
    # dg1=degree(p1,gen=var)
    # dg2=degree(P4,gen=var)
    # if dg1<dg2:
        # vec1=coef_list(p1,size=dg2,var2=var)
        # vec2=coef_list(P4,var2=var)
    # elif dg1>dg2:
        # vec1=coef_list(p1,var2=var)
        # vec2=coef_list(P4,size=dg1,var2=var)
    # else: 
        # vec1=coef_list(p1,var2=var)
        # vec2=coef_list(P4,var2=var)
        
    # vec3=[vec2[i]-vec1[i] for i in range(len(vec2))]
    # kres=solve(vec3)
    # svar,vval=ganswer(kres,'varname','value')
    # for i,j in zip(svar,vval):
        # P=P.subs(i,j)
    # return P
    
#  getinsidepar

def vecfind(sexp,sval):
    sexp=str(sexp)
    mv=[]  
    qq=len(sexp)
    if sexp[-1]==sval:
        mv=[qq]
    sexp=str(sexp)
    mm=[]
    done=True
    while done:
        p=sexp.find(sval)
        if p==-1:
            done=False
        else:    
            mm.append(p)
            sexp=sexp[p+1:-1]
    mm=mm+mv 
    return mm 
    
def insidepar(expr,svar=''):
    if svar=='':
        sexpr=str(expr)
        p1=sexpr.find('(')
        vp2=vecfind(sexpr,')')
        qq=len(vp2)
        sres=sexpr[p1+1:vp2[qq-1]+1]
        return parse_expr(sres)
    else:
        sexpr=str(expr)
        pp=sexpr.find(svar)
        sexpr2=sexpr[pp::]
        if sexpr2[0]=='(':
            sres= sexpr2
        else:
             
            p3=sexpr2.find('(')
            sres=sexpr2[p3::]
                
        vec1=vecfind(sres,'(') 
        vec2=vecfind(sres,')') 
        mm=[]
        for i in vec1:
            mm.append([i,1])
        for i in vec2:
            mm.append([i,-1]) 
        mm.sort()
        qq=len(mm)
        p1=0
        p2=mm[qq-1][0]
        cc=1
        for i in range(1,qq):
            cc=cc+mm[i][1]
            if cc==0:
                p2=mm[i][0]
                break
        kres=sres[p1:p2]
        return parse_expr(kres)
            
## Pprimefactor...


def factorp(expr):
    kres=''
    if Is_Add(expr):
        done=True
        for i in expr.args:
            if done:
                kres=factorp(i)
                done=False
            else:
                kres=kres+'+'+factorp(i)
         
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        P1=factorp(p1)
        P2=factorp(p2)
        kres= '('+P1+')/('+P2+')'
    elif Is_Mul(expr):
        done=True
        for i in expr.args:
            if done:
                kres=factorp(i)
                done=False
            else:
                kres=kres+'*('+factorp(i)+')'
         
    elif Is_Root(expr):
        bb=insideroot(expr)
        rr=getroot(expr)
        kres=factorp(bb)
        kres='('+kres+')**('+str(frs(1,rr))+')'
            
     
    elif Is_NumberPow(expr):
        bb=getbase(expr)
        ee=getexpo(expr)
        vecbase,vecexp=factoresprimosvec(bb)
        done=True
        for i,j in zip(vecbase,vecexp):
            newe=simplify(j*ee)
            if done:
                kres=str(i)+'**('+str(newe)+')'
                done=False
            else:
                kres=kres+'*'+str(i)+'**('+str(newe)+')'
               
    elif Is_Number(expr):
        mm=factorint(expr)
        vecbase,vecexp=unpack(mm)
        qq=len(vecbase)
        done=True
        for i,j in zip(vecbase,vecexp):
            if j!=1:
                if done:
                    sres=str(i)+'**'+str(j) 
                    done=False
                else:
                    sres=sres+'*'+str(i)+'**'+str(j) 
            else:
                if done:
                    sres=str(i) 
                    done=False
                else:
                    sres=sres+'*'+str(i) 
        kres=sres
    
                
         
    else:
         kres=str(expr)
    return kres
    
def factorp(expr):
    mm=factorint(expr)
    vecbase,vecexp=unpack(mm)
    return vecbase,vecexp 


def getfaclist(expr):
    vecb=[]
    vece=[]
    mm=expr.args
     
    for i in mm:
        try:
            base,eee=unpack(factorint(i))
            for k,j in zip(base,eee):
                vecb.append(k)
                vece.append(j)
        except:
            vecb.append(i)
            vece.append(1)
    kres=''
    done=True
     
    for i,j in zip(vecb,vece):
        if j==1:
            subres=str(i)
        else:
            subres=str(i)+'**'+str(j)
        if done:
            kres=subres
            done=False
        else:
            kres=kres+'*'+subres
    return sympify( kres,evaluate=False)                
            

def primefactor(expr):
    if Is_Add(expr):
         
        kres=''
        done=True
        for i in expr.args:
            if done:
                kres=str(primefactor(i))
                done=False
            else:
                kres=kres+"+"+str(primefactor(i))
          
    elif Is_Div(expr):
         
        p1,p2=fraction(expr)
        p1=str(primefactor(p1))
        p2=str(primefactor(p2))
        kres= p1+'/'+p2
   
    elif Is_Mul(expr):
         
        
        kres=''
        done=True
        mm=expr.args
        for i in mm:
            if done:
                kres= str(primefactor(i))
                done=False
            else:
                kres= kres+'*'+str(primefactor(i))
    elif Is_Log(expr):
         
        mm=expr.args
        expr2=mm[0]
        kres=str(primefactor(expr2))
        kres= 'log('+kres+')'
    
    

    elif type(expr)==Rational:
         
        knume=numer(expr)
        kdeno=denom(expr)
        knum=primefactor(knume)
        kden=primefactor(kdeno)
        return  cfrac(knume,kdeno)
        
    
    elif type(expr)==int or type(expr)==Integer:
         
        if isprime(expr):
            kres=str(expr)
        else:    
            val=expr
            vecn,vece=kunpakDic(factorint(val))
            newf=''
            for i,j in zip(vecn,vece):
                if j==1:
                    newf=str(i)+'*'+newf
                else:    
                    newf=str(i)+'**'+str(j)+'*'+newf    
            newf=newf[0:-1]
            newf=str(newf)
            kres=sympify(newf,evaluate=False)
         
   
    

    elif type(expr)==Pow:
         
        base=getbase(expr)
        ee=getexpo(expr)
        if Is_Symbol(base) or isprime(base):
            kres=str(expr)
         
        else: 
            if type(base)==Rational:
                return (primefactor(base))**ee
                
            elif not isprime(base):  
                vec1,vec2=unpack(factorint(base))
                done=True
                kres=''
                for i,j in zip(vec1,vec2):
                    if done:
                        kres= str(i)+'**'+'('+str(j*ee)+')'
                        done=False
                    else:
                        kres=kres+'*'+str(i)+'**'+'('+str(j*ee)+')'
                     
            else:    
                kres= getfaclist(expr)
    elif Is_Root(expr):
         
        rr=getroot(expr)
        B=insideroot(expr)
        try:
            kres=primefactor(B)
            ee=getexpo(kres)
            bb=getbase(kres)
            return bb**cfrac(ee,rr)
        except:
            return expr
    else:
         
        kres=expr
    return  (sympify(kres,evaluate=False))
    
    
    
 
        
def superfactor(expr,kdiv,*args):
    ''' 
    superfactor(exp2factorize, factor expr, op)
        op = 
            'factor',factorize before doit
            'simplify' simplify ...
            'expand'..
            'expo', applied suerfactor in exponenets
    '''        
    if 'factor' in args:
        expr=factor(expr)
    if 'expand' in args:
        expr=expand(expr)
    if 'simplify' in args:
        expr=simplify(expr)
    if Is_Pow(expr):
        bb=superfactor(getbase(expr),kdiv)
        ee=getexpo(expr)
        if 'expo' in args:
            ee=superfactor(expo,kdiv)
        return ppow(bb,ee)
    if Is_Div(expr):
        p1,p2=fraction(expr)
        return cfrac(superfactor(p1,kdiv),superfactor(p2,kdiv))
    
    if  Is_Add(expr):
        P1=0
        P2=0
        
        for i in fpoly(expr,'list'):
            i=superfactor(i,kdiv)
            kres=i/kdiv
            p1,p2=fraction(kres)
            if p2==1:
                P1+=p1
            else:
                P2+=i
        return P2+P1*kdiv
    elif '+' in str(expr) or '-' in str(expr):
        kres=1
        for i in fpoly(expr,'list'):
            if Is_Add(i):
                kres=kres*superfactor(i,kdiv)
            else:
                kres=kres*i
        return kres
    else:
        return expr  


      
# basic ope 


        
def Div(p1,p2):
    kres= p1/p2
    try:    
        if int(kres)==kres:
            return int(kres)
    except:
        return kres
#  math whit list

def SubstracList(vec1,vec2):
    '''
        vec1=[a1,b1,..]
        vec2=[a2,b2,..]
        return [a1-a2,b1-b2,...]
    '''
    
    kres=[]
    for i,j in zip(vec1,vec2):
        kres.append(i-j)
    return kres
        
def AddList(vec1,vec2):
    '''
        vec1=[a1,b1,..]
        vec2=[a2,b2,..]
        return [a1+a2,b1+b2,...]
    '''
    kres=[]
    for i,j in zip(vec1,vec2):
        kres.append(i+j)
    return kres

def MulList(vec1,vec2):
    '''
        vec1=[a1,b1,..]
        vec2=[a2,b2,..]
        return [a1*a2,b1*b2,...]
    '''
    kres=[]
    for i,j in zip(vec1,vec2):
        kres.append(i*j)
    return kres 
def DivList(vec1,vec2):
    '''
        vec1=[a1,b1,..]
        vec2=[a2,b2,..]
        return [a1/a2,b1/b2,...]
    '''
    kres=[]
    for i,j in zip(vec1,vec2):
        kres.append(cfrac(i,j))
    return kres    
    
    
def symbolslist(vec,var=''):
    '''
        symbolslist(x+a+b+4)=[x,a,b]
        symbolslist(x+a+b+4,x)=[a,b]
        symbolslist([x+a+b+4,a,z])=[a,b,x,z]
        symbolslist([x+a+b+4,a,z],x)=[a,b,xz]

    '''
    if type(vec)==list:
        qq=len(vec)
        kres=0
        cc=10
        for i in vec:
            kres=expand(kres+i/cc)
            cc=cc*10
    else:
        kres=vec
        
    slist=kres.free_symbols
    vecs=list(slist)
    if var!='':
        kres=[]
        for i in vecs:
            if i!=var:
                kres.append(i)
    else:
        kres=slist
    return kres

# DEGREE LIST ALGORTHMICS

def pre_vec2compare(vec1,vec2):
    '''
        filter (0=2) in coefflist

    '''
    nv1=[]
    nv2=[]
    for i,j in zip(vec1,vec2):
        if not Is_Number(i) or not Is_Number(j):
            nv1.append(i)
            nv2.append(j)
    return nv1,nv2

def pre_coeff2list(expr1,expr2,var=x):
    '''
    prepare 2 coefflist two not generate error
    '''
    
    d1=degree(expr1,var)
    d2=degree(expr2,var)
    d3=max(d1,d2)
    vec1=coef_list(expr1,var,d3)
    vec2=coef_list(expr2,var,d3)
    vec1,vec2=pre_vec2compare(vec1,vec2)
    return vec1,vec2


# imaginary number
def module_img(expr):
    p1=re(expr)
    p2=im(expr)
    return get_hipo(p1,p2)

def imodule(expr):
    p1=re(expr)
    p2=im(expr)
    return get_hipo(p1,p2)    
    
def powi(expr,ee):
    return expand(expr**ee)
   
   
def rpart(expr):
    if Is_Div(expr):
        p1,p2=fraction
    sres=str(expr)
    sres=sres.replace('I','0')
    kres=parse_expr(sres)
    return kres

def repart(expr):
    sres=str(expr)
    sres=sres.replace('I','0')
    kres=parse_expr(sres)
    return kres


def imgpart(expr):
    p1=expr-repart(expr)
    p2=str(p1)
    p2=p2.replace('I','1')
    p3=parse_expr(p2)
    return p3
    
def ipart(expr):
    p1=expr-repart(expr)
    p2=str(p1)
    p2=p2.replace('I','1')
    p3=parse_expr(p2)
    return p3    

def simplify_img(expr):
    if Is_Add(expr):
        kres=(repart(expr))+imgpart(expr)*I
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        p1=simplify_img(p1)
        p2=simplify_img(p2)
        return cfrac(p1,p2)
    elif Is_Mul(expr):
        return simplify(expr)
    else:
        return simplify(expr)
    return kres 

def polar2img(*args):
    '''
    polar2img(sqrt(2),pi/4)
    return  1+i
    '''
    p1=args[0]
    p2=args[1]
    angle=p2
    radio=p1
    if 'pi' in str(p1):
        angle=p1
        radio=p2
    pre=radio*cos(angle)
    pim=radio*sin(angle)
    sres=str(pre)+'+'+str(pim)+'*I'
    return parse_expr(sres)    

def img2polar(expr):
    '''
    img2polar(1+I)
    return sqrt(2),pi/4
    '''
    p1=re(expr)
    p2=im(expr)
    p3=module_img(expr)
    alpha=atan(cfrac(p2/p1))
    alpha=simplify(alpha)
    return p3,alpha

def conjugada(z):
    if Is_Pow(z):
        bb=getbase(z)
        ee=getexponent(z)
        bb2=conjugada(bb)
        ee2=conjugada(ee)
        return bb**ee2
    else:
        sz=str(z)
        sz=sz.replace('I','-I')
        return parse_expr(sz)

def premath(expr,*args):
    '''
    pre trasforms answer in simplesolve
    '''
    if len(args)>0:
        if 'factor' in args:
            try:
                expr=factor(expr)
            except:
                pass
        if 'simplify' in args or 'simpli' in args:
            try:
                expr=simplify(expr)
            except:
                pass    
        if 'expand' in args:
            try:
                expr=expand(expr)
            except:
                pass    
        if 'reduce' in args:
            try:
                expr=reduce(expr)
            except:
                pass    

        if 'unisymbols' in args or 'unisys' in args:
            try:
                expr=unisymbols(expr)
            except:
                pass    
        if 'mulexpo' in args:
            try:
                expr=mulexpo(expr)
            except:
                pass    
        if 'rsimplify' in args:
            try:
                expr=rsimplify(expr)
            except:
                pass 
    return expr  

   
def doindenom(*args):
    expr=args[0]
    func=args[1]
    p1,p2=fraction(expr)
    if len(args)==3:
        expr2=args[2]
        p2=dothis(p2,func,expr2)
    else:
        p2=dothis(p2,func,expr2)
    
    return unisymbols(Div(p1,p2))
    
    
def doinnumer(*args):
    expr=args[0]
    func=args[1]
    p1,p2=fraction(expr)
    if len(args)==3:
        expr2=args[2]
        p1=dothis(p1,func,expr2)
    else:
        p1=dothis(p1,func,expr2)
    
    return unisymbols(Div(p1,p2))
 
