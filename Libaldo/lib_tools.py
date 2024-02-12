from sympy import *
from lib_Mathematica import *
from lib_Mathbasic import *
from lib_Algorith import *

from lib_MyEq import *
from lib_MyEqEq import *
 
 
def getexpr(expr):
    if type(expr)==MyEq:
        return expr.ksym
    elif type(expr)==MyEqEq:
        return expr.e1.ksym-expr.e2.ksym
    else:
        return expr
 
def compareandsolve(*args,var=x,eQ=True):
    vec,var=presolvecoeff(*args,var=var)
    VAR=[]
    VSOLU=[]
    done=True
    while len(vec)>0 and done:
        for i in vec:
            vS=symbolslist(i)
            if len(vS)==1:
                var=list(vS)
                var=var[0]
                kres=solve(i,var)
                kres=kres[0]

                VAR.append(var)
                VSOLU.append(kres)
                break
        vec2=[]
        for i in vec:
            vec2.append(i.subs(var,kres))
        vec3=[]
        for i in vec2:
            vS=symbolslist(i)
            if len(vS)>0:
                vec3.append(i)
        vec=vec3
    if eQ:
        eeres=[]
        for i,j in zip(VSOLU,VAR):
            ee=MyEq(i,str(j),var=var)
            eeres.append(ee)
        return eeres
    else:
        return VAR, VSOLU
def simplecompareandsolve(*args):
    VAR=[]
    VSOLU=[]
    done=True
    vec=args[0] 
     
    while len(vec)>0 and done: 
        for i in vec:
            vS=symbolslist(i)
            if len(vS)==1:
                var=list(vS)
                var=var[0]
                kres=solve(i,var)
                kres=kres[0]

                VAR.append(var)
                VSOLU.append(kres)
                break
        vec2=[]
        for i in vec:
            vec2.append(i.subs(var,kres))
        vec3=[]
        for i in vec2:
            vS=symbolslist(i)
            if len(vS)>0:
                vec3.append(i)
        vec=vec3
     
    eeres=[]
    for i,j in zip(VSOLU,VAR):
        ee=MyEq(i,str(j),var=var)
        eeres.append(ee)
    return eeres
def presolvecoeff(*args,var=x):
    vec,vecv=presolvedata(*args,var=var)
    vec=vec[0]
    kres=expand(vec)
    kres=sortdegree(kres,var)
    kres=coef_list(kres,var)
    
    #vec1=coef_list(vec,var)
    return kres,vecv
    
def presolvedata(*args,var=x):
    '''
    return ( vec expr to solve, vec variables)...
    '''
    vecf=[]
    vecv=[]
    for i in args:
        if type(i)==MyEqEq:
            vecf.append(i.L-i.R)
        elif type(i)==MyEq:
            vecf.append(i.ksym)
        elif Is_Symbol(i):
            vecv.append(i)
        elif type(i)==list:
            if Is_Symbol(i[0]):
                vecv=vecv+i
        else:
             vecf.append(i)
    if vecv==[]:
        p3=sum(vecf)
        vecv=symbolslist(p3,var)
    return vecf,vecv 


    
def solveSys(*args,var=x,eQ=True):
    '''
        input args[0]= [31𝑎+6𝑏+339, −30𝑎−5𝑏−320]
              args[1]=[a,b]
        return
            a=-9
            b=10
    '''
    done=False    
    if type(args[0])!=list:
        args2=[]
        for i in args:
            if type(i)!=str:
                args2.append(i)
            if i=='Eq':
                done=True
        args=args2        
        args=presolvedata(*args,var=var)
        
    vec=args[0]
    vec2=[]
    for i in vec:
        if type(i)==MyEq:
            vec2.append(i.ksym)
        elif type(i)==MyEqEq:
            vec2.append(i.L-i.R)
        else:
            vec2.append(i)
            
    var=args[1]
    kres=solve_poly_system(vec2,var)
    kres=kres[0]
    mm=[]
    if eQ:
        for i,j in zip(var,kres):
            display(Math(str(i)+'= '+latex(j)))
        return kres 
    
    
    else: 
        return kres
      
        
def pre_solvecoef(v1,v2):
    V1=[]
    V2=[]
    for i,j in zip(v1,v2):
        if i==0 and j==0:
            pass
        else:
            V1.append(i)
            V2.append(j)
    return (V1,V2)
    
    
def MySolve(*args,eQ=True,variables=False):
    veceq=[]
    vecvar=[]
    if len(args)==1 and type(args[0])==list:
        args=args[0]
    for i in args:
        if type(i)==MyEqEq:
            veceq.append(i.L-i.R)
        elif type(i)==MyEq:
            veceq.append(i.ksym)
            
        elif Is_Symbol(i):
            vecvar.append(i)
        else:
            veceq.append(i)
            
    if len(vecvar)==0:
        for i in veceq:
            lsym=symbolslist(i)
            for j in lsym:
                if j not in vecvar:
                    vecvar.append(j)
         
        vecss=[str(k) for k in vecvar]
        vecss.sort()
        vecvar=[symbols(k) for k in vecss]
     
    if variables:
        return solveSys(veceq,vecvar,eQ=False), vecvar
    else:
        return solveSys(veceq,vecvar,eQ=eQ)
        
# SOLVE OEFFIECIENT LIST ATB
   
def pre_vec2compare(vec1,vec2): #atb

    '''
    pre_vec2compare([a,b,c,5,d],[z,1,3,2,x]) =[a,b,c,d],[z,1,3,x])
                    preveent that 5=2 
     
    '''
    nv1=[]
    nv2=[]
    for i,j in zip(vec1,vec2):
        if not Is_Number(i) or not Is_Number(j):
            nv1.append(i)
            nv2.append(j)
    return nv1,nv2    

def pre_coeff2list(expr1,expr2,var=x): #atb
    '''
    match two list of coeff with the same haight degree
    input (math exp, math exp)
    output( coeff_vec,coeff_vec])
        
        pre_coeff2list(x**4+1,a*x**4+b*x+c])=
        ([1,0,0,0,1],[a,0,0,b,c])
    
    '''
    d1=degree(expr1,var)
    d2=degree(expr2,var)
    d3=max(d1,d2)
    vec1=coef_list(expr1,var,d3)
    vec2=coef_list(expr2,var,d3)
    vec1,vec2=pre_vec2compare(vec1,vec2)
    return vec1,vec2 

def pre_coeff0(expr,var=x): #atb
    '''
    input: mathexp,main var
    return: coeffvec
    
    pre_coeff0((a+1)*x**3+(b-2)*x*x+4,x)= [a+1,b-2]
        a+1=0,b+a=0, 4=0???
        
    '''    
    expv=coef_list(expr,var)
    kres=[]
    for i in expv:
        vres=symbolslist(i,var)
        if len(vres)>0:
            kres.append(i)
    return kres 


def solvecoefficients(*args,var=x): #atb
    '''
    return solutions when compare coeffircient in expr 
    or between two expr,
    the algorith try to get var to solve or
    set who var will be solve 
    
    input args:
        * MyEqEq
        * MyEq
        * matexpr
        * MyEq,matexpr
        * MyEq,MyEq
        * matexpr,matexpr
        * MyEqEq,var1,var2..
        * MyEq,var1,var2..
        * matexpr,var1,var2..
        * MyEq,matexpr,var1,var2..
        * MyEq,MyEq,var1,var2..
        * matexpr,matexpr ,var1,var2..  
        
    '''
    vecss=[]
    for i in args:
        if Is_Symbol(i):
            vecss.append(i)
    done=False
    if len(vecss)>0:
        done=True
    
    if len(args)==1:
        data=args[0]
        if type(data)==MyEqEq:
            p1=data.L
            p2=data.R
            vec1,vec2=pre_coeff2list(p1,p2,var=var)
            vec3=SubstracList(vec1,vec2)
            p3=p1+p2
        elif type(data)==MyEq:
            p1=data.ksym
            #vec3=coef_list(p1,var)
            vec3=pre_coeff0(p1,var)
            p3=p1
        else:
            p1=data 
            #vec3=coef_list(data,var)
            vec3=pre_coeff0(p1,var)
             
            p3=p1
    elif len(args)>1:
        p1=args[0]
        if type(p1)==MyEqEq:
            P1=p1.L
            P2=p1.R
            vec1,vec2=pre_coeff2list(P1,P2,var=var)
            vec3=SubstracList(vec1,vec2)
            p3=P1+P2
        elif type(p1)==MyEq:
            p1=p1.ksym
            p3=p1
            vec3=pre_coeff0(p1,var)
            if type(args[1])!=Symbol:
                p2=args[1]
                if type(p2)==MyEq:
                    p2=p2.ksym
                p3=p1-p2
                vec1,vec2=pre_coeff2list(p1,p2,var=var)
                vec3=SubstracList(vec1,vec2)
        else:
            vec3=coef_list(p1,var)
            p3=p1

             
    if done:
        vecv=vecss
    else:
        vecv=symbolslist(p3,var)

    vecargs=[]
    for i in vec3:
        vecargs
    return solveSys(vec3,vecv,'Eq')
    
    
def solvecoefflist(*args,var=x):
    vecss=[]
    for i in args:
        if Is_Symbol(i):
            vecss.append(i)
    done=False
    if len(vecss)>0:
        done=True
    
    if len(args)==1:
        data=args[0]
        if type(data)==MyEqEq:
            p1=data.L
            p2=data.R
            vec1,vec2=pre_coeff2list(p1,p2,var=var)
            vec3=SubstracList(vec1,vec2)
            p3=p1+p2
        elif type(data)==MyEq:
            p1=data.ksym
            #vec3=coef_list(p1,var)
            vec3=pre_coeff0(p1,var)
            p3=p1
        else:
            p1=data 
            #vec3=coef_list(data,var)
            vec3=pre_coeff0(p1,var)
             
            p3=p1
    elif len(args)>1:
        p1=args[0]
        if type(p1)==MyEqEq:
            P1=p1.L
            P2=p1.R
            vec1,vec2=pre_coeff2list(P1,P2,var=var)
            vec3=SubstracList(vec1,vec2)
            p3=P1+P2
        elif type(p1)==MyEq:
            p1=p1.ksym
            p3=p1
            vec3=pre_coeff0(p1,var)
            if type(args[1])!=Symbol:
                p2=args[1]
                if type(p2)==MyEq:
                    p2=p2.ksym
                p3=p1-p2
                vec1,vec2=pre_coeff2list(p1,p2,var=var)
                vec3=SubstracList(vec1,vec2)
        else:
            vec3=coef_list(p1,var)
            p3=p1

             
    if done:
        vecv=vecss
    else:
        vecv=symbolslist(p3,var)
        
    #return vec3,vecv
    return solveSys(vec3,vecv,'Eq')      
def solvecoeflist(*args):
    if len(args)==1 and type(args[0])==MyEqEq:
        expr=args[0]
        p1=expr.L
        p2=expr.R 
        ce1=coef_list(p1)
        ce2=coef_list(p2)
        ce1,ce2=pre_solvecoef(ce1,ce2)
        vecf=[]
        for i, j in zip(ce1,ce2):
            vecf.append(i-j)
        kres=solve(vecf)
        return kres
    elif len(args)==2 and type(args[0])==MyEqEq:
        expr=args[0]
        p1=expr.L
        p2=expr.R 
        ce1=coef_list(p1)
        ce2=coef_list(p2)
        ce1,ce2=pre_solvecoef(ce1,ce2)
        vecf=[]
        for i, j in zip(ce1,ce2):
            vecf.append(Eq(i,j))
        kres=solve(vecf,args[1])
        if type(kres)==list:
            if type(kres[0])==tuple:
                kres=kres[0]
                for i,j in zip(args[1],kres):
                    display(Math(str(i)+'= '+latex(str(j))))
        return kres
        
def preanswer(*args):
    kres=args[0]
    rres=kres
    if 'real' in args:
        rres=[G for G in kres if not 'I' in str(G)]
    if type(rres)==list and len(rres)==1:
        rres=rres[0]
    if type(rres)==dict:
        vecs,vval= kunpakDic(rres)
        if 'Eq' in args:
            mm=[]
            for i,j in zip(vecs,vval):
                ee=MyEq(j,str(i))
                mm.append(ee)
            return mm
        else:
            return vval    

##########################################
#   ALgorith FUnctions

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
    
    ee=MyEq(var-expr,'ee',var=yy,kshow=False)
      
    try:
        kres=ee.solve(yy,kshow=False)
    except:
        kres=0
    return kres
            
##########################################
#  PARTIAL FRACTION RE ULTRA NERD ALGORITH

def mdaughter(expr,fabr,var=x):
    '''
    retunr correspondients momonies that
    will be used to create partial fraction..
    exaple:
        mdaughter((x+1)) = y1
        mdaughter((x+1)*(x-2)) = y1,y2
        mdaughter((x+1)*(x-2)**3) = y1,y2,y3,y4,y5
        mdaughter((x+1)*(x*x-2)**2) = y1,(x*y2+y3),(x*y3+y4)
    '''    
    
    bb=getbase(expr)
    ee=getexpo(expr)
    dd=degree(bb,gen=var)
    mm2=[]
    cc=0
    for i in range(ee):
        mm3=0
        for j in range(dd):
            mm3=mm3+var**j*fabr[cc]
            cc=cc+1
        mm2.append(mm3)
    frb2=fabr[cc:len(fabr)]
    return mm2,frb2
    
    
def vec_daughter(expr2,var=x): # used in partialfracction
    #    retun vector whit respevctive nuerator formato por mdauther
    y1,y2,y3,y4,y5,y6,y7,y8,y9,y10=symbols('y1 y2  y3 y4 y5 y6 y7 y8 y9 y10')
    y11,y12,y13,y14,y15,y16=symbols('y11 y12 y13 y14 y15 y16')
    fabr=[y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16]
    if Is_Mul(expr2):
        mmm=fpoly(expr2,'list')
    else:
        mmm=[expr2]
    kres=[]
    for j in mmm:
        exp1,fabr=mdaughter(j,fabr,var=var)
        for i in exp1:
            kres.append(i)
    return kres 

def complete2partfrac(vec):
    # complete missing monomies to run partial fracction algorith
    vec2=[]
    if Is_Mul(vec):
        vecmono=fpoly(vec,'list')
    else:
        vecmono=[vec]
    for i in vecmono:
        if Is_Pow(i):
            ee=getexpo(i)
            bb=getbase(i)

            for j in range(1,ee+1):
                vec2.append(ppow(bb,j))
        else:
            vec2.append(i)
    return vec2

def partialfraction(expr,var=x):
    '''
    input : x/((x+1)*(x-1))
    return :-1/(x - 1) + 2/(x - 2)
    '''
    if type(expr)==MyEq:
            expr=expr.ksym
    Q=0    
    kp1,kp2=fraction(expr)
    gg1=degree(kp1,gen=var)
    gg2=degree(kp2,gen=var)
     
    if gg1>gg2 or gg1==gg2:
        Q=quo(kp1,kp2)
        R=rem(kp1,kp2)
         
        expr=cfrac(R,kp2)

     
 
     
     
    p1=numer(expr)
    p2=denom(expr)
    vec2=complete2partfrac(p2)
    vec1=vec_daughter(p2,var=var)
     
    P=0
    for i,j in zip(vec1,vec2):
        P=P+cfrac(i,j)

    P2=factor(P) 
    pp1,pp2=fraction(P2)

    pp1=sortdegree(pp1,var)
    pp2=sortdegree(expand(p1),var)
    gg1=degree(pp1,var)
    gg2=degree(pp2,var)
    ksize=max(gg1,gg2)           
    PP1=coef_list(pp1,var2=var,size=ksize)
    PP2=coef_list(pp2,var2=var,size=ksize)

    masterq=[]
    for i,j in zip(PP1,PP2):
        masterq.append(i-j)
     
    vecsolu,vecvar= MySolve(masterq,variables=True)
    for i,j in zip(vecvar,vecsolu):
        P=P.subs(i,j)
    return P+Q 


# diferential

def deltasym(var):
    #sdx= +str(var)
    sdx= +str(var)
    dx=symbols(sdx)
    return dx

def solveQ(*args):
    vecq=[]
    for i in args:
        if type(i)!=str:
            if type(i)==MyEqEq:
                vecq.append(i.L-i.R)
            elif type(i)==MyEq:
                vecq.append(i.ksym)

            else:
                vecq.append(i)
    kres=solve(vecq)
     
    if type(kres)==list:
        Kres=[]
        fval=kres[0]
        kvar,kval=unpack(fval)
        for i in kres:
            kvar2,kval2=unpack(i)
            Kres.append(kval2) 
        qq=len(kres)
         
        qq2=len(kvar)
        sk=[]
        for j in range(qq2):
            sk2=[]
            for i in range(qq):
                sk2.append(Kres[i][j])
            sk.append(sk2)    
        for i,j in zip(kvar,sk):
            display(Math(str(i)+" = "+latex(j)))
        return sk    
     
         
            
    if type(kres)==dict:
        kvar,kval=unpack(kres)
    
        if 'value' in args:
            return kval
        if 'eQ' in args:
            ee=[]
            for i,j in zip(kvar,kval):
                ee.append(MyEq(j,str(i)))
            return ee              
        
    return kres

def expr2var(expr,svar):
    var=symbols(svar)
    expr=subsubs(expr,svar,var)
    return expr


def getdataexpr(expr):
    if type(expr)==MyEqEq:
        return expr.L-expr.R 
    elif type(expr)==MyEq:
        return expr.ksym
    else:
        return expr 

def getvecvar(*args):
    vecvar=[]
    obj=obj2expr(*args)
    for i in obj:
        if type(i)!=str:
            kres=list(i.free_symbols )
            for j in kres:
                if type(j)==Symbol:
                    if j not in vecvar:
                        vecvar.append(j)

    return vecvar       

def getvecexpr(*args):
    vecexpr=[]
    for i in args:
        if type(i)!=str:
            if type(i)==list or type(i)==tuple:
                for j in i:
                    if type(j)!=Symbol:
                        kres=getdataexpr(j)
                        if j not in vecexpr:
                            vecexpr.append(float2int(kres))
            else:
                if type(i)!=Symbol:
                    kres=getdataexpr(i)
                    if i not in vecexpr:
                        vecexpr.append(kres)
    return vecexpr 

######################
## SIMPLE SOLVE

def filtermathexpr(expr):
    '''
    return significative math equation in expr
    else return ''
    '''
    if type(expr)==MyEqEq:
         
        kres= float2int(expr.L-expr.R)
        return kres
    elif type(expr)==MyEq:
        return float2int(expr.ksym)
    elif '+' in str(expr) or '-' in str(expr):
        return float2int(expr)
    else:
        return '' 

def getmatexpr(*args):
    '''
    return all posible mathh equation
    to will be evaluates in solve
    input,list,tuple,MyEq MyEqEq, math expr
    '''
    vecexpr=[]
    prevec=[]
    for i in args:
         
        if type(i)==Matrix:
            H,W=i.shape
            for ii in range(H):
                for jj in range(W):
                    kres=filtermathexpr(i[ii,jj])
                    if kres!='':
                        if type(kres)==list:
                            for j in kres:
                                vecexpr.append(kres)
                        else:
                            vecexpr.append(kres)
                    
                    
        elif type(i)==list or type(i)==tuple or type(i)==Array:
            for j in i:
                kres=filtermathexpr(j)
                if kres!='':
                    if type(kres)==list:
                        for j in kres:
                            vecexpr.append(float2int(kres))
                    else:
                        vecexpr.append(float2int(kres))
        else:
            kres=filtermathexpr(i)
            if kres!='':
                if type(kres)==list:
                    for j in kres:
                        vecexpr.append(float2int(kres))
                else:
                    vecexpr.append(float2int(kres))
    vecexpr= unisymbols(vecexpr)               
    return vecexpr

def presolve(*args):
    vecvar=[]

    for i in list(args):
        if type(i)==Symbol:
            vecvar.append(i)
    vecmath=getmatexpr(*args)
    vecvar2=[]
    for i in list(vecmath):
        vtup=i.free_symbols
        vlist=list(vtup)
        for j in vlist:
            if j not in vecvar2:
                vecvar2.append(j)
    return vecvar,vecvar2,vecmath      

def premath(expr,*args):
    '''
    pre trasforms aswer in simplesolve
    '''
    if 'factor' in args:
        expr=factor(expr)
    if 'simplify' in args:
        expr=simplify(expr)
    if 'expand' in args:
        expr=expand(expr)
    if 'reduce' in args:
        expr=reduce(expr)
    if 'unisymbols' in args:
        expr=unisymbols(expr)
    if 'float' in args:
        try:
           expr=float(expr)
        except:
            pass
    return expr
def spreanswer(expr):
    mm2=str(expr)
    smm=mm2.replace(',), (',',')
    smm=smm.replace('[(','[')
    smm=smm.replace(')]',']')
    try:
        mm3=parse_expr(smm)
        return mm3
    except:
        return expr
def prepresolve(*args):
    args2=[]
    for i in args:
        if type(i)==MyEq:
            if type(i.ksym)==list or type(i)==Array:
                if type(i)==Array:
                    L=list(i.ksym)
                else:
                    L=i.ksym
                for j in L:
                    args2.append(j)
            else:
                arg2.append(i)
        elif type(i)==MyEqEq:
            if type(i.L)==list or type(i.L)==Array:
                if type(i.L)==Array:
                    L=list(i.L-i.R)
                else:
                    L=i.L-i.R
                for j in L:
                    args2.append(j)    
            else:
                arg2.append(i)
        elif type(i)==Array:
            L=list(i)
            for j in L:
                args2.append(j)
                
        elif type(i)==list:
            for j in L:
                args2.append(j)       
        else:
            args2.append(i)
    return args2

def prepresolve(*args):
    args2=[]
    for i in args:
        if type(i)==MyEq:
            if type(i.ksym)==list or type(i.ksym)==Array:
                if type(i)==Array:
                    L=list(i.ksym)
                else:
                    L=i.ksym
                for j in L:
                    args2.append(j)
            else:
                args2.append(i)
        elif type(i)==MyEqEq:
            if type(i.L)==list or type(i.L)==Array:
                if type(i.L)==Array:
                    L=list(i.L-i.R)
                else:
                    L=[p2-p1 for p2,p1 in zip(i.L,i.R)] 
                for j in L:
                    args2.append(j)    
            else:
                args2.append(i)
        elif type(i)==Array:
            L=list(i)
            for j in L:
                args2.append(j)
                
        elif type(i)==list:
            for j in L:
                args2.append(j)       
        else:
            args2.append(i)
            
    return args2
    
def simplesolve(*args,kshow=True,**kwargs):
    '''
    input mathexpr eq, MyE,;yEqEq, tuple or list of math eq 
    and variables to solve
    if not was input variables 
    the system assum deafult there are in all expr..
    
    example:
        pp=(x+y-4,3*x-5*y+1)
        A=MyEq(z*3+y,'A')
        Q=MQ(3*sin(alpha),27)
        
        x,y,z=simplesolve(pp,A,Q,x,y,z, op= 'factor','simplify','expand','reduce'
                                    'unisymbols',
                                    'show'= nice display)
        return x,y,x values in solve the eqss
        

    '''
    args2=prepresolve(*args)
    args=args2
    vecvar,vecvar2,vecexpr=presolve(*args)
    vecmath=[]
    if vecvar==[]:
        vecvar=vecvar2
    qq=len(vecvar)
    if len(kwargs)>0:
        vecexpr=[real_subs(i,**kwargs) for i in vecexpr]
    #return vecmath,vecvar 
    kres=solve(vecexpr,vecvar)
    kres=spreanswer(kres) 
     
    if type(kres)==dict:
        svar,valor=unpack(kres)
        kres=valor 
    if 'all' in args:
        if 'noimg' in args or 'noimaginary' in args:
            kres2=[]
            for i in kres:
                if not 'I' in str(i):
                    kres2.append(i)
            if kshow:
                display(Math(latex(kres2)))
            return kres2
        else:    
            if kshow:
                display(Math(latex(kres)))
            return kres     
    elif type(kres)==tuple:
        kres=list(kres)
        svar=[str(i) for i in vecvar]

    elif type(kres)==list:
        if type(kres[0])==tuple:
            kres=list(kres[0])
        elif type(kres[0])==list:
            kres=kres[0]
        svar=[str(i) for i in vecvar]

    else:
        kres=[kres]
        svar=[str(i) for i in vecvar]
    
    kres2=[]
    for i in kres:
        kres2.append(premath(i,*args))
       
    if 'all' in args:
        return kres2
    
    if 'noshow' not in args:
         
        vecee=[] 
        for i,j in zip(kres2,svar):
            if 'pos' in args or 'positive' in args:
                ee=MyEq(i*signo(i),kname=str(j),kshow=kshow)
                 
            else:
                ee=MyEq(i,kname=str(j),kshow=kshow)
            vecee.append(ee)
        kres2=vecee     
    if 'value' in args:
        kres2=[i.ksym for i in kres2]
    if len(kres2)==1:
        kres2=kres2[0]
        
    return kres2
def simplesolves(*args, **kwargs):
     
    vecexpr, vecvar, vecstr2 = obj2vec3(*args)
    ops=['noshow','float']
    ops2=[]
    vecstr=[]
    for i in vecstr2:
        if i not in ops:
            vecstr.append(i)
        else:
            ops2.append(i)
    if vecvar == []:
        vecvar = getvecvar(*vecexpr)
    if len(kwargs)>0:
        vecexpr=[real_subs(i,**kwargs) for i in vecexpr]
    if vecvar == []:
        kres = solve(*vecexpr)
    else:
        kres = solve(vecexpr, vecvar)
     
    if type(kres) == dict:
        svar, valor = unpack(kres)       
        kres = valor
        if 'float' in args:
            kres2=[]
            for i in valor:
                try:
                    kres2.append(float(i))
                except:
                    kres2.append(i)
            kres=kres2
            valor=kres2
    elif type(kres)==list:
        if type(kres[0])==tuple:
            cc=1
            valor=[]
            svar=[]
            for i in kres:
                for j,k in zip(vecvar,i):
                    svar.append(sympify(str(j))+sympify(str(cc)))
                    if 'float' in args:
                        print('a')
                        try:
                            valor.append(float(k))
                        except:
                            valor.append(k)
                    else:
                        valor.append(k)
                cc=cc+1
            kres=valor
         
            
                
    
        
    if 'noshow' not in ops2:
        for i,j in zip(svar,valor):
            kvar=i
            svar=str(i)
            svar=svar.replace('_','')
            display(Math(svar+'='+latex(j)))
    if len(kres)==1:
        return kres[0]
    return kres 
'''
def simplesolve(*args,positive=''):
   
    input mathexpr eq, MyE,;yEqEq, tuple or list of math eq 
    and variables to solve
    if not was input variables 
    the system assum deafult there are in all expr..
    
    example:
        pp=(x+y-4,3*x-5*y+1)
        A=MyEq(z*3+y,'A')
        Q=MQ(3*sin(alpha),27)
        
        x,y,z=simplesolve(pp,A,Q,x,y,z, op= 'factor','simplify','expand','reduce'
                                    'unisymbols',
                                    'show'= nice display)
        return x,y,x values in solve the eqss
        
 
    
    
    vecvar,vecvar2,vecexpr=presolve(*args)
    
    vecmath=[]
    if vecvar==[]:
        vecvar=vecvar2
    qq=len(vecvar)
    
    #return vecmath,vecvar 
    if positive!='':
        positive=Symbol(positive, positive=True)
    kres=solve(float2int(vecexpr),vecvar)
    if 'all' in args:
        return solve(vecexpr)
     
    if 'float' in args:
        kres2=[]
        for i in kres:
            try:
                kres2.append(float(i))
            except:
                kres2.append(i)
        kres=kres2
         
    if type(kres)==dict:
        svar,valor=unpack(kres)
         
        kres=valor
        
         
    elif type(kres)==tuple:
        kres=list(kres)
        svar=[str(i) for i in vecvar]

    elif type(kres)==list:
        if type(kres[0])==tuple:
            kres=list(kres[0])
            if 'pos' in args or 'positive' in args:
                kres=kres[1]
 
        elif type(kres[0])==dict:
            kres=kres[0]
            if 'pos' in args or 'positive' in args:
                kres=kres[1]
            
        svar=[str(i) for i in vecvar]

    else:
        kres=[kres]
        svar=[str(i) for i in vecvar]

    kres2=[]
    for i in kres:
        kres2.append(premath(i,*args))
        
    if 'eq' in args or 'Eq' in args:
        kres3=[]
        for i,j in zip(kres2,svar):
            ee= MyEq(i,j)
            kres3.append(ee)
        if type(kres3)==list and len(kres3)==1:
            return kres3[0]
        else:        
            return kres3
        
         
 
        
    if 'noshow' not in args:
        vecee=[] 
        for i,j in zip(kres2,svar):
            mkres=i
            if 'positive' in args:
                mkres=signo(i)*i
                
            ee=MyEq(mkres,kname=str(j))
                
            vecee.append(ee)

    if type(kres2)==list and len(kres2)==1:
        kres=kres2[0]
        if 'positive' in args:
            kres=signo(kres)*kres
        return kres    
    return kres2        
            
'''
def showtuple(svar,vvar):
    for i,j in zip(svar,vvar):
        sp1=str(i)+'='
        display(Math(sp1+latex(j)))
        
####  tools Argssss

vecreatr=["<class 'sympy.core.symbol.Symbol'>","<class 'int'>","<class 'float'>","<class 'sympy.core.numbers.Pi'>","<class 'sympy.core.numbers.Rational'>"]
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
        
def str2vec(sexpr):
    kvec=[]
    for i in sexpr:
        kvec.append(int(i))
    return kvec

def argsmap (expr,kname='A',deep=2,data=False):
    infoexpr=[]
    infopos=[]
    cc=''
    A,B,C,D=ruta(expr,infoexpr,infopos,cc)
    mapval=[]
    mappos=[]
    for i,j in zip(B,C):
        if Is_Div(i):
            if numer(i)!=1:
                mapval.append(i)
                mappos.append(j)

        else:
            mapval.append(i)
            mappos.append(j)
    
    mapval,mappos=filterNegPos(mapval,mappos)
    if len(mapval)==0:
        return
    if data:
        return mapval,mappos
    svecargs=[]
    sres=''
    superres=''
    if kname!='':
         
        for i in mappos:
            sres=kname+'.args('
            for k in i:
                sres=sres+k+','
            sres=sres[0:-1]
            sres=sres+')='
            svecargs.append(sres) 
        mm=''
        for i,j in zip(svecargs,mapval):
            mm=mm+ "  "+'['+i+latex(j)+'],'
        display(Math(mm))
            
def realsub2(expr,**kwargs):
    vecs,vecv=unpack(kwargs)
    for i,j in zip(vecs,vecv):
        if i in nombresdiferen:
            expr=expr.subs(eval(i),j)
        elif type(j)==MyEq:
            expr=expr.subs(i,j.ksym)
        else:
            expr=expr.subs(i,j)
    return expr
   
def transform2(expr1,expr2,var=x):
    A,B,C,D=symbols('A B C D')
    vecvar=[A,B,C,D]
    svecvar=['A','B','C','D']
    sexpr2=str(expr2)
    vecres=[]
    for i,j in zip(svecvar,vecvar):
        if i in sexpr2:
            vecres.append(j)
    Q=expand(expr1)-expand(parse_expr(sexpr2))
    nvec= coef_list(Q,var) 
    vecres2=simplesolve(nvec,'noshow')
    if type(vecres2)!=list:
        vecres2=[vecres2]
    expr3=parse_expr(expr2)
    for i,j in zip(vecres,vecres2):
        expr3=expr3.subs(i,j)
    return expr3       
    
import math
import numpy as np

# Main Function takes in the coefficient of the Cubic Polynomial
# as parameters and it returns the roots in form of numpy array.
# Polynomial Structure -> ax^3 + bx^2 + cx + d = 0

def solve3(a, b, c, d):

    if (a == 0 and b == 0):                     # Case for handling Liner Equation
        return np.array([(-d * 1.0) / c])                 # Returning linear root as numpy array.

    elif (a == 0):                              # Case for handling Quadratic Equations

        D = c * c - 4.0 * b * d                       # Helper Temporary Variable
        if D >= 0:
            D = math.sqrt(D)
            x1 = (-c + D) / (2.0 * b)
            x2 = (-c - D) / (2.0 * b)
        else:
            D = math.sqrt(-D)
            x1 = (-c + D * 1j) / (2.0 * b)
            x2 = (-c - D * 1j) / (2.0 * b)
            
        return np.array([x1, x2])               # Returning Quadratic Roots as numpy array.

    f = findF(a, b, c)                          # Helper Temporary Variable
    g = findG(a, b, c, d)                       # Helper Temporary Variable
    h = findH(g, f)                             # Helper Temporary Variable

    if f == 0 and g == 0 and h == 0:            # All 3 Roots are Real and Equal
        if (d / a) >= 0:
            x = (d / (1.0 * a)) ** (1 / 3.0) * -1
        else:
            x = (-d / (1.0 * a)) ** (1 / 3.0)
        return np.array([x, x, x])              # Returning Equal Roots as numpy array.

    elif h <= 0:                                # All 3 roots are Real

        i = math.sqrt(((g ** 2.0) / 4.0) - h)   # Helper Temporary Variable
        j = i ** (1 / 3.0)                      # Helper Temporary Variable
        k = math.acos(-(g / (2 * i)))           # Helper Temporary Variable
        L = j * -1                              # Helper Temporary Variable
        M = math.cos(k / 3.0)                   # Helper Temporary Variable
        N = math.sqrt(3) * math.sin(k / 3.0)    # Helper Temporary Variable
        P = (b / (3.0 * a)) * -1                # Helper Temporary Variable

        x1 = 2 * j * math.cos(k / 3.0) - (b / (3.0 * a))
        x2 = L * (M + N) + P
        x3 = L * (M - N) + P

        return np.array([x1, x2, x3])           # Returning Real Roots as numpy array.

    elif h > 0:                                 # One Real Root and two Complex Roots
        R = -(g / 2.0) + math.sqrt(h)           # Helper Temporary Variable
        if R >= 0:
            S = R ** (1 / 3.0)                  # Helper Temporary Variable
        else:
            S = (-R) ** (1 / 3.0) * -1          # Helper Temporary Variable
        T = -(g / 2.0) - math.sqrt(h)
        if T >= 0:
            U = (T ** (1 / 3.0))                # Helper Temporary Variable
        else:
            U = ((-T) ** (1 / 3.0)) * -1        # Helper Temporary Variable

        x1 = (S + U) - (b / (3.0 * a))
        x2 = -(S + U) / 2 - (b / (3.0 * a)) + (S - U) * math.sqrt(3) * 0.5j
        x3 = -(S + U) / 2 - (b / (3.0 * a)) - (S - U) * math.sqrt(3) * 0.5j

        return np.array([x1, x2, x3])           # Returning One Real Root and two Complex Roots as numpy array.


# Helper function to return float value of f.
def findF(a, b, c):
    return ((3.0 * c / a) - ((b ** 2.0) / (a ** 2.0))) / 3.0


# Helper function to return float value of g.
def findG(a, b, c, d):
    return (((2.0 * (b ** 3.0)) / (a ** 3.0)) - ((9.0 * b * c) / (a **2.0)) + (27.0 * d / a)) /27.0


# Helper function to return float value of h.
def findH(g, f):
    return ((g ** 2.0) / 4.0 + (f ** 3.0) / 27.0)

import math
import numpy as np

# Main Function takes in the coefficient of the Cubic Polynomial
# as parameters and it returns the roots in form of numpy array.
# Polynomial Structure -> ax^3 + bx^2 + cx + d = 0

def solve3(a, b, c, d):

    if (a == 0 and b == 0):                     # Case for handling Liner Equation
        return np.array([(-d * 1.0) / c])                 # Returning linear root as numpy array.

    elif (a == 0):                              # Case for handling Quadratic Equations

        D = c * c - 4.0 * b * d                       # Helper Temporary Variable
        if D >= 0:
            D = math.sqrt(D)
            x1 = (-c + D) / (2.0 * b)
            x2 = (-c - D) / (2.0 * b)
        else:
            D = math.sqrt(-D)
            x1 = (-c + D * 1j) / (2.0 * b)
            x2 = (-c - D * 1j) / (2.0 * b)
            
        return np.array([x1, x2])               # Returning Quadratic Roots as numpy array.

    f = findF(a, b, c)                          # Helper Temporary Variable
    g = findG(a, b, c, d)                       # Helper Temporary Variable
    h = findH(g, f)                             # Helper Temporary Variable

    if f == 0 and g == 0 and h == 0:            # All 3 Roots are Real and Equal
        if (d / a) >= 0:
            x = (d / (1.0 * a)) ** (1 / 3.0) * -1
        else:
            x = (-d / (1.0 * a)) ** (1 / 3.0)
        return np.array([x, x, x])              # Returning Equal Roots as numpy array.

    elif h <= 0:                                # All 3 roots are Real

        i = math.sqrt(((g ** 2.0) / 4.0) - h)   # Helper Temporary Variable
        j = i ** (1 / 3.0)                      # Helper Temporary Variable
        k = math.acos(-(g / (2 * i)))           # Helper Temporary Variable
        L = j * -1                              # Helper Temporary Variable
        M = math.cos(k / 3.0)                   # Helper Temporary Variable
        N = math.sqrt(3) * math.sin(k / 3.0)    # Helper Temporary Variable
        P = (b / (3.0 * a)) * -1                # Helper Temporary Variable

        x1 = 2 * j * math.cos(k / 3.0) - (b / (3.0 * a))
        x2 = L * (M + N) + P
        x3 = L * (M - N) + P

        return np.array([x1, x2, x3])           # Returning Real Roots as numpy array.

    elif h > 0:                                 # One Real Root and two Complex Roots
        R = -(g / 2.0) + math.sqrt(h)           # Helper Temporary Variable
        if R >= 0:
            S = R ** (1 / 3.0)                  # Helper Temporary Variable
        else:
            S = (-R) ** (1 / 3.0) * -1          # Helper Temporary Variable
        T = -(g / 2.0) - math.sqrt(h)
        if T >= 0:
            U = (T ** (1 / 3.0))                # Helper Temporary Variable
        else:
            U = ((-T) ** (1 / 3.0)) * -1        # Helper Temporary Variable

        x1 = (S + U) - (b / (3.0 * a))
        x2 = -(S + U) / 2 - (b / (3.0 * a)) + (S - U) * math.sqrt(3) * 0.5j
        x3 = -(S + U) / 2 - (b / (3.0 * a)) - (S - U) * math.sqrt(3) * 0.5j

        return np.array([x1, x2, x3])           # Returning One Real Root and two Complex Roots as numpy array.


# Helper function to return float value of f.
def findF(a, b, c):
    return ((3.0 * c / a) - ((b ** 2.0) / (a ** 2.0))) / 3.0


# Helper function to return float value of g.
def findG(a, b, c, d):
    return (((2.0 * (b ** 3.0)) / (a ** 3.0)) - ((9.0 * b * c) / (a **2.0)) + (27.0 * d / a)) /27.0


# Helper function to return float value of h.
def findH(g, f):
    return ((g ** 2.0) / 4.0 + (f ** 3.0) / 27.0)  

def fmaxmin(obj,var=x):
    ksym=obj2func(obj)
    kres=diff(ksym,var)
    solu=solve(kres,var)
    if type(solu)==list and len(solu)==1:
        return solu[0]
    return solu
def finflextion(obj,var=x):
    ksym=obj2func(obj)
    kres=diff(ksym,var,var)
    solu=solve(kres,var)
    if type(solu)==list and len(solu)==1:
        return solu[0]
    return solu



###  SUPERSOLVE

def presolvedata(*args):
    '''
    input: all diferents exprs..MyEq,MQ, func var..all..
    retur
        vecf=[math expre finding]
        vecv=[symbols required]
        vecs=[str expretions]
    '''
    vecf=[]
    vecv=[]
    vecs=[]
    args2=[]
    for i in args:
        if type(i)==list:
            for j in i:
                args2.append(j)
        else:
            args2.append(i)
    args=args2

    for i in args:
        if type(i)==MyEqEq:
            vecf.append(i.e1.ksym-i.e2.ksym)
        elif type(i)==MyEq:
            vecf.append(i.ksym)
        elif type(i)==Symbol:
            vecv.append(i)
        elif type(i)==str:
            vecs.append(i)
        else:
            vecf.append(i)
    return vecf,vecv,vecs 

def getvecsymbols(*args):
    '''
    return symbols finding in whatwvwer expresions 
    '''
    
    vecs=[]
    for i in args:
        try:
            kres=i.free_symbols
            kres=list(kres)
            for j in kres:
                if j not in vecs:
                    vecs.append(j)
        except:
            pass
    return vecs

def supersol(*args):
    '''
    find answer to eQ
    '''
    vecv=None
    vecf,vecv,vecs=presolvedata(*args)
    if vecv==None or vecv==[]:
        vecv=getvecsymbols(*vecf)
         
    kres=solve(vecf,vecv)
    return kres,vecv    

def presolanswer(kres,vecv):
    '''
    fix, sort, and interprete the most answer format
    '''
    
    vecval=[]
    vecnam=[]
    vecsubname=[]
    if len(vecv)==1:
        if type(kres)==dict:
            vecnam,vecval=unpack(kres)
        elif type(kres)==list:
            cc=1 
            for i in kres:
                 
                il=list(i)
                for j,k in zip(il,vecv):
                    vecval.append(j)
                    vecnam.append(str(k)+str(cc))
                     
                cc=cc+1    
        return vecval,vecnam            
    else:
        if type(kres)==dict:
            vecnam,vecval=unpack(kres)
            
        elif type(kres)==list:
            cc=1 
            for i in kres:
                 
                il=list(i)
                for j,k in zip(il,vecv):
                    vecval.append(j)
                    vecnam.append(str(k)+str(cc))
                     
                cc=cc+1    
        return vecval,vecnam

def supersolve(*args,kshow=True):

    '''
    supersolve(fumc, Myeq, equations,..etc,,,var1,var2,var3...)
    '''
    
    kres,vecv=supersol(*args)
    valores,nombres=presolanswer(kres,vecv)
     
    ee=[]
    for i,j in zip(valores,nombres):
        ee.append(MyEq(i,kname=str(j),kshow=kshow))
    if 'eq' in args or 'Eq' in args:
        if type(ee)==list:
            return ee 
        else:
            return ee[0]
    else:
        if type(valores)==list:
            return valores 
        else:
            return valores[0]     
            
def Lsubs(vec,*args,**kwargs):
    kres=[]
    if type(vec[0])==list:
        for i in vec:
            kres.append(lsubs(i,*args,**kwargs))
        return kres
    else:
        return lsubs(vec,*args,**kwargs)
                        
def lsubs(vec,*args,**kwargs):
    vecres=[]
    if len(kwargs)>0:
        for i in vec:
            kres=i
            try:
                kres=real_subs(kres,**kwargs)
            except:
                pass
            vecres.append(kres)
    else:
        vecres=vec
    for i in args:
        if type(i)==MyEq: 
            sexpr=i.name
            expr=unisymbols(i.ksym)
            vecres2=[]
            for j in vecres:
                try:
                    kres=j.subs(sexpr,expr)
                    vecres2.append(kres)
                except:
                    vecres2.append(j)
            vecres=vecres2
    return vecres


def smartanswer(expr, *args,**kwargs):
    ops = ['simplify', 'factor', 'float', 'expand','rsimplify']
    kname = ''
    kops = []
    kvar = ''
    expr=real_subs(expr,**kwargs)
         
    for i in args:
        if type(i) == str:
            if i in ops:
                kops.append(i)
            else:
                kname = i
        if type(i) == Symbol:
            kvar = i
    if len(kops)>0:
        if 'expand' in kops:
            try:
                expr = expand(expr)
            except:
                pass
        if 'simplify' in kops:
            try:
                expr = simplify(expr)
            except:
                pass
        if 'factor' in kops:
            try:
                expr = factor(expr)
            except:
                pass    
        if 'float' in kops:
            try:
                expr = float(expr)
            except:
                pass
        if 'rsimplify' in kops:
            try:
                expr = rsimplify(expr)
            except:
                pass        
    if kname != '':
        if kvar !='':
            ee=MyEq(float2int(expr),kname=kname,var=kvar)
        else:
            ee=MyEq(float2int(expr),kname=kname)
        return ee
    else:
        return float2int(expr)

def obj2expr(*args):
    vecres=[]
    for i in args:
        if type(i)==MyEq:
            vecres.append(i.ksym)
        elif type(i)==MyEqEq:
            vecres.append(i.L-i.R)
        else:
            vecres.append(i)
    if len(vecres)==1:
        return vecres[0]
    else:
        return vecres
    
    
def mathobj2expr(obj):
    if type(obj)==str:
        return None
    elif type(obj)==MyEq:
        return  obj.ksym 
    elif type(obj)==MyEqEq:
        return  simplify(expand(obj.L-obj.R)) 
    else:
        return obj 

def realangle(*args):
    vx=args[0]
    vy=args[1]
    vxx=vx*signo(vx)
    vyy=vy*signo(vy)
    if vxx!=0:
        angle=atan(vyy/vxx)
    if vx==0:
        if vy>=0:
            kres= pi/2
        else:
            kres= -pi/2
    elif vy==0:
        if vx>=0:
            kres= 0
        else:
            return pi
    elif vy>0:
        if vx>0:
            kres= angle
        else:
            kres= pi/2+angle
    else:
        if vx>0:
            kres= -angle
        else:
            kres= pi+angle
    if 's' in args:
        kres= simplify(rad2sex(kres))
    if 'float' in args:
        kres= float(kres)
    return kres

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

def pivec(pi1,pi2,pix):
    '''
    return vector = [ pi1,pi1+pix,pi1+2*pix,..... pi2]
    '''
    done=True
    vecpi=[]
    cc=0
    while done:
        vecpi.append(pi1+pix*cc)
        cc=cc+1
        if pix*cc>pi2:
            done=False
    return vecpi 

def discretesolvefunc(expr,x1,x2,x3,var=x):
    '''
    try to retur vector ehit posible solutions of expr for 
    evaluate range( x1, x2) in x3 spaces...
    
    '''
    vec=pivec(x1,x2,x3)
    ee=obj2MyEq(expr,var=var)
    vecres=[]
    for i in vec:
        try:
            solu=simplify(ee(i))
            if solu==0:
                vecres.append(i)
        except:
            pass
    if len(vecres)==1:
        return vecres[0]
    return vecres

def sympyintgral2expr(expr):
    '''
    return sympy integral expreion whitout integral symbols, include dx
    '''
    slatex=latex(expr)
    dstr=slatex[-2::]
    dvar=symbols(dstr)
    slatex=slatex.replace('\\int ','')
    p1=slatex.find('\\')
    slatex=slatex[0:p1]
    nexpr=latex2sympy(slatex)
    if nexpr==1:
        return dvar
    else:
        return nexpr*dvar
        

def vecsubs(expr,vec1,vec2):
    '''
    return subtitutio follow vec1,vec2 indication
    '''
    expr=unisymbols(expr)
   
    for i,j in zip(vec1,vec2):
        try:
            expr=expr.subs(i,j)
        except:
            pass
         
    return expr
def completesquare(expr,var=x):
    '''
    expr=x*x+2*x
    if x*x+2*x+a=(x+b)**2
    return a,b
    ''' 
    a, b=symbols('a b')
    Q=MQ(expr,expand((var+b)**2),kshow=False)
    vec1=Q.coef_list(var)
    vec2=[]
    for i in vec1:
        if 'a' in str(i) or 'b' in str(i):
            vec2.append(i)
    a,b=simplesolve(*vec2,a,b,'noshow')
    return -b,a
def coef0(obj):
    '''
    return coefficient zero in f(x,y)
    '''
    
    expr=obj2expr(obj)
    w=list(expr.free_symbols)
    kres=0
    for i in expr.args:
        done=True
        for j in w:
            if str(j) in str(i):
                done=False
        if done:
            kres=kres+i        
    return kres     
def squarecompletexy(*args,var1=x,var2=y):
    var1=symbols(str(var1))
    var2=symbols(str(var2))
    exprx=args[0]
    var=var1
    Mx=exprx.subs(var2,0)
    Ry=exprx-Mx
    Mx2=squarecomplete(Mx,x)
    R=coef0(Mx2)
    P1=Mx2-R
    My=Ry+R
    var=y
    My2=squarecomplete(My,var)
    R=coef0(My2)
    P2=My2-R
    if 'parts' in args:
        return P1,P2,R
    else:
        return P1+P2+R

def squarefill(*args):
    '''
    squarefill(x*x+2*x+1,x)  return (x+1)**2
    squarefill(2*x*x+4*x+1,x,'factor') return 2*(x+1)**2-2    factor
    squarefill(2*x*x+4*x+1,x,'factor',parts) return (x+1)**2,-2,2     
 
    '''
    var=x
    vecop=['factor','parts']
    fdone=False
    cf=1
    nexpr=args[0]
    for i in args:
        if type(i)==MyEq:
            expr=i.ksym
        elif type(i)==Symbol:
            var=i
        elif i=='factor':
            fdone=True 
        elif type(i)!=str:
            expr=i
        else:
           pass  
    sresto=0
    w=list(expr.free_symbols)
    if len(w)>1:
        expr8=expr
        for i in w:
            if i!=var:
                expr8=expr8.subs(i,0)
        nexpr=expr8        
        sresto=expr-nexpr
        
    x1,x2,x3=symbols('x1,x2,x3')
    vecx=[x1,x2,x3]
    
    if fdone:
         
        mm=nexpr.args
        mm2=[i.subs(var,j) for i,j in zip(mm,vecx)]
        expr2=factor(sum(mm2))
        expr3=reducecero(expr2)
        cf=abs(simplify(expr2/expr3))
        nexpr=sum(mm)/cf
    clist=coef_list(nexpr,var)
    a=sqrt(clist[0])
    b=(clist[1])/(2*a)
    c=b*b
    d=clist[2]-c
    
    if 'parts' in args:
        return (a*var+c),d*cf+sresto,cf
    else:    
        if fdone:
            return  cf*(a*var+c)**2+d*cf +sresto
        else:
            return (a*var+c)**2+d +sresto


def square2polar(obj,var=t,kshow=True):
    expr=obj2expr(obj)
    mm=expr.args
    mx=sum([i for i in mm if 'x' in str(i)])
    my=sum([i for i in mm if 'y' in str(i)])
    mz=[mx,my]
    Nm=[numer(i) for i in mz]
    Dn=[denom(i) for i in mz]
    R=[rsimplify(sqrt(i)) for i in Nm]
    ss=[Nm[0].subs(x,0),Nm[1].subs(y,0)]
    kres=[unisymbols(ss[0]+sqrt(Dn[0])*cos(var)),unisymbols(ss[1]+sqrt(Dn[1])*sin(var))]
    if kshow:
        display(Math(latex(kres)))
    return kres


#######################
########  used to rsimplify
################

def Is_RootMonoPow(expr):
    done=False
    if Is_Root(expr):
        if Is_Pow(insideroot(expr)):
            return True
    return False
def Sqrt(val,sq,ee=1):
    sexpr='\\sqrt['+str(sq)+']'
    if ee==1:
        sexpr=sexpr+'{'+str(val)+'}'
    else:    
        sexpr=sexpr+'{'+str(val)+'^{'+str(ee)+'}}'
    return latex2sympy(sexpr) 

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

def makethis(expr,*args):
    vecveri=['unisymbols','factor','expand','simplify','texpand','rsimplify','Sqrt','reducecero']
    vecop=[i for i in vecveri if i in args]
    for i in vecop:
        try:
            sexpr=i+'('+str(expr)+')'
            expr=eval(sexpr)
        except:
            pass
    return expr


def solution(obj,var):
    kres=obj2expr(obj)
    skres=str(kres)
    sdone=False
    if type(var)==str:
        sdone=True
        svar=var
        var=parse_expr(var) 
    else:
        svar=str(var)
    mm=simplesolve(kres,var,'all',kshow=False)
    qq=len(mm)
    mms=[svar+str(cc) for cc in range(1,qq+1)]
    if sdone:
        evec=[]
        for eval,ename in zip(mm,mms):
            evec.append(MyEq(eval,kname=ename))
        return evec    
            
    else:
        if qq==1:
            display(Math(svar+'='+latex(mm[0])))
            return mm[0]
        else:
            return mm    