from sympy import *
from lib_Mathbasic import *
from lib_Algorith import *

from lib_MyEq import *
from lib_MyEqEq import *
 
import copy

def getobjdata(*args):
    vecobj=[]
    for data in args:
        if type(data)==list:
            for data2 in data:
                vecobj.append(data2)
        else:
                vecobj.append(data)
                
    vecexpr=[]
    vecsymb=[]
    for data in vecobj:
        
        if type(data)==MyEqEq:
            vecexpr.append(data.L-data.R)
        elif type(data)==MyEq:
            vecexpr.append(data.ksym)
        elif type(data)==symbols:
            vecsymb.append(data)
        else:
            vecexpr.append(data)
    return vecexpr+vecsymb
    
        
def obj2mathexpr(*args):
    vecres=[]
    for i in args:
        if type(i)==MyEq:
            vecres.append(i.ksym)
        elif type(i)==MyEqEq:
            vecres.append(i.e1.ksym-i.e2.ksym) 
        else:
            vecres.append(i)
    if len(vecres)==1:
        return vecres[0]
    else:
        return vecres
        
        
def getdataQ(*args):
    ops=['float','value','eQ','equation']
    vecexp=[]
    vecvar=[]
    vecnam=[]
    vecslist=[]
    args2=[]
    for data in args:
        if type(data)==list:
            for data2 in data:
                args2.append(data2)
        else:
            args2.append(data)
    args=args2        
    for data in args:
        if Is_String(data):
            if not data in ops:
                vecnam.append(data)
        elif Is_symbols(data):
            vecvar.append(data)
            vecnam.append(str(data))
        elif type(data)==list:
            for  item in data:
                vecexp.append(item)        
        else:
            vecp=list(data.free_symbols)
            if vecp!=[]:
                if not vecp in vecslist:
                    vecslist.append(vecp)
                    vecexp.append(obj2mathexpr(data))
    return vecexp,vecvar,vecnam
    
def predataQ(*args):
    vecexp,vecvar,vecnam=getdataQ(*args)
    if vecvar==[]:
        vecvar2=[]
        for data in vecexp:
            vres=list(data.free_symbols)
            for qvar in vres:
                if not qvar in vecvar2:
                    vecvar2.append(unisymbols(qvar))
        vecvar=vecvar2
    for data in vecvar:
        vecnam.append(str(data))
    return vecexp,vecvar,vecnam

def ndegree(vecexp,vecvar):
    etype=1
    if len(vecvar)>1:
        try:
            var=vecvar[0]
            lvar=vecvar[1::]
            qq=len(vecvar)
            for nv in lvar:
                nvec=[]
                for data in vecexp:
                    nvec.append(data.subs(nv,var))
                vecexp=nvec
        except:
            pass
    for dataq in vecexp:
        for datav in vecvar:
            kres=degree(dataq,datav)
            if kres>etype:
                etype=kres
    return etype 

    

def simplesolve2(*args,show=True,**kwargs):
    vecvar=[]
    vecexp=[]
    vecnam=[]
    qq2=0
    ops=['float','value','eq','equation','noimg']
    for data in args:
        if not data in ops:
            if Is_symbols(data):
                vecvar.append(data)
            elif type(data)==list:
                for  item in data:
                    vecexp.append(item)        
            else:
                vecexp.append(obj2mathexpr(data))
    for data in vecvar:
         vecnam.append(str(data))
    if len(kwargs)>0:
        vecexp2=[]
        for data in vecexp:
            vecexp2.append(real_subs(data,**kwargs))
        vecexp=vecexp2    
    if vecvar==[]:
        for data in vecexp:
            vecv=list(data.free_symbols)
            for dvar in vecv:
                if not dvar in vecvar:
                    vecvar.append(dvar)
        vecname=[str(data) for data in vecvar]            

    kres=solve(vecexp,vecvar)
     
    if type(kres)==dict: 
        sval,vvar=list(kres.keys()),list(kres.values())
        sname=[str(data) for data in sval]
        
    elif type(kres)==list:    
        if type(kres[0])==tuple:
            if len(kres)>1:
                qq=len(kres)
                qq2=len(list(kres[0]))
            else:
                vvar=kres[0]
                sval=vecvar[0]
                sname=str(sval)
                qq2=0
        else:
            if len(kres)==1:
                vvar=kres[0]
                sval=vecvar[0]
                sname=str(sval)
                qq2=0
        
    if qq2>0:
        nvecval=[]
        nvecname=[]
        for i in range(0,qq2):
            for cc in range(qq):
                lvalor=list(kres[i])
                nvecval.append(lvalor[i])
                nvecname.append(vecname[i]+str(i+1))
        sname=nvecname
        vvar=nvecval
 
    if 'float' in args:
        vecflo=[]
        for data in vvar:
            try:
                vecflo.append(float(data))
            except:
                vecflo.append(data)
        vvar=vecflo        
    veceq=[]
    for kname,kval in zip(sname,vvar):
        e1=MyEq(kval,kname)
        veceq.append(e1)
    if  'eq' in args or 'equation' in args:
        return veceq
    else:
        if type(vvar)==list and len(vvar)==1:
            return  vvar[0]
        else:
            return vvar     


def kunpack(**kwargs):
    sval=[]
    vval=[]
    for key, value in kwargs.items():
        sval.append(key)
        vval.append(value)
    return sval,vval    
        
    
 
def get_independent_term(expr,*args):
    if len(args)>0:
        lvar=[data for data in args]
    else:    
        lvar=list(expr.free_symbols)
    for data in lvar:
        expr=expr.subs(data,0)
    return expr

def completesquare(expr,var=x):
    '''
    expr=x*x+2*x
    if x*x+2*x+a=(x+b)**2
    return a,b
    ''' 
 
    vec=coef_list(expr,var)
    vec2=[cfrac(data,vec[0]) for data in vec] 
    nb=cfrac(vec2[1],2)
    nc=nb**2
    kres=spow((var+nb),2)-nc+vec[2]
    return kres



def SquareComplete(expr,*args):
    kres=0
    if len(args)>0:
        lvar=[data for data in args]
    else:    
        lvar=list(expr.free_symbols)
    tindep=get_independent_term(expr,*args)
    onlyv=expr-tindep
    for data in lvar:
        onlyv=expr-tindep
        lvar2=list(set(lvar)-{data})
        for data2 in lvar2:
            onlyv=onlyv.subs(data2,0)
        kres=kres+completesquare(onlyv,var=data)
    return kres +  tindep 
           
            
            
    
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
        elif type(i)==symbols:
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


def square2polar(obj,var=t,show=True):
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
    if show:
        display(Math(latex(kres)))
    return kres  
    
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
    
def completesquarexyz(QQQ,*args,var=t,show=True):
    '''
    input MyEq=A*x*x+B*y*y+C*z*z+D*x+E*y+F*z*G
    return Cx,Cy,Cx,R2
    '''

    QQ=copy.deepcopy(QQQ)
    if type(QQ)!=MyEq:
        QQ=MyEq(QQ,'QQ',show=False)
    R=QQ(x=0,y=0,z=0)
    QQ.Substrac(R,show=False)
    ex2=x*x*vecsubs(QQ.ksym,[x**2,y**2,z**2,x,y,z],[1,0,0,0,0,0])
    ey2=y*y*vecsubs(QQ.ksym,[x**2,y**2,z**2,x,y,z],[0,1,0,0,0,0])
    ez2=z*z*vecsubs(QQ.ksym,[x**2,y**2,z**2,x,y,z],[0,0,1,0,0,0])
    ex=x*vecsubs(QQ.ksym,[x**2,y**2,z**2,x,y,z],[0,0,0,1,0,0])
    ey=y*vecsubs(QQ.ksym,[x**2,y**2,z**2,x,y,z],[0,0,0,0,1,0])
    ez=z*vecsubs(QQ.ksym,[x**2,y**2,z**2,x,y,z],[0,0,0,0,0,1])
    X,ax=completesquare(ex2+ex+a,var=x)
    Y,ay=completesquare(ey2+ey+a,var=y)
    Z,az=completesquare(ez2+ez+a,var=z)
    P=factor(simplify(expand(ax+ay+az-R)))
     
    if  'center' in  args:
        if show:
            display(Math('CenterPoint=('+latex(X)+','+latex(Y)+','+latex(Z)+')'))
            display(Math('R2='+latex(P)))
        return Point3D(X,Y,Z),P
    elif 'parametric' in args:
        vecee=[]
        ee=MQ((x-ax)/X,cos(var)**2)
        vecee.append(ee)
        ee=MQ((y-ay)/Y,sin(var)**2)
        vecee.append(ee)
        return (vecee)
            
    else:
        if show:
            display(Math('Cx='+latex(X)))
            display(Math('Cy='+latex(Y)))
            display(Math('Cz='+latex(Z)))
            display(Math('R2='+latex(P)))
        return (X,Y,Z,P)   
        
def obj2float(expr):
    if type(expr)==Add:
        kres=0
        for data in expr.args:
            kres+=expr2float(data)
        return kres
    elif type(expr)==Mul:
        kres=1
        for data in expr.args:
            kres=kres*expr2float(data)
        return kres    
    else:
        try:
            kres=float(expr)
        except:
            kres=expr
        return kres

# matematic expr to sympy srags map symbols to analize equation structure

# input: 3*x+3 outout : 'A(M(I,S),I)'
# input: x**3+x+3*log(y)  outout : 'A(P(S,I),S,M(I,L(S)))'
def math2map(expr):
    mm=srepr(expr)
    done=True
    cc=0
    while done:
        try:
                
            p1,p2=find_simetric_end_par(mm,cc)
            if not'(' in mm[p1::]:
                done=False
            elif p2-p1<6 or onlynumber(mm[p1:p2]):
                mm=mm.replace(mm[p1:p2+1],'*')
                cc=cc+1
            else:
                cc=p1+1
        
        except:
            done=False
    vecs=['symbols*','Integer*','Mul','Add','Pow','Rational*',' ','log']
    vece=['S','I','M','A','P','R','','L']
    for v1,v2 in zip(vecs,vece):
        qq=mm.count(v1)
        for i in range(qq):
            mm=mm.replace(v1,v2)
    return mm

def onlynumber(sexpr):
    vec='0123456789, ()'
    done=True
    for data in sexpr:
        if not data in  vec:
            return False
    return True

def find_simetric_end_par(expr,p1):
    sexpr=str(expr)
    p2=find_first_par_from(str(expr),p1)
    done=True
    qbal=1
    p3=p2+1
    while done:
        if sexpr[p3]=='(':
            qbal=qbal+1
        if sexpr[p3]==')':
            qbal=qbal-1
        if qbal==0:
            return p2,p3
        else:
            p3=p3+1
def find_first_par_from(sexpr,p1):
    sres='('
    done=True
    while done:
        svar=sexpr[p1]
        if svar=='(':
            return p1
        else:
            p1=p1+1  


def solve_coefcompare(obj,var):
    p1=obj.e1.ksym
    p2=obj.e2.ksym
    L1=coef_list(p1,var)
    L2=coef_list(p2,var)
    vece=[]
    for data1,data2 in zip(L1,L2):
        vece.append(data1-data2)  
    return simplesolve(*vece)

 
 
def eQandsymbols(*args):
    vece1=[]
    vece2=[]
    vece3=[]
    vece4=[]
    vecsymbols=[]
    for data in args:
        if Is_symbols(data):
            vecsymbols.append(data)
            
    for data in args:
        if type(data)!=str:
            if type(data)!=symbols:
                vece1.append(data)
    for data in vece1:
        if type(data)==list:
            for data2 in data:
                vece2.append(data2)
        else:
            vece2.append(data)
    vecgrup=[]        
    for data in vece2:
        if type(data)==MyEqEq:
            expr=data.e1.ksym-data.e2.ksym
            if 'I' in str(expr):
                vece3.append(im(expr))
                if re(expr)!=0:
                    vece3.append(re(expr))
            else:
                vece3.append(expr)
        elif type(data)==MyEq:
            vece3.append(data.ksym)
        else:
            vece3.append(data)
    for data in vece3:
        try:
            veri=list(data.free_symbols)
            if len(veri)>0:
                if not data in vece4:
                    vecgrup.append(veri)
                    vece4.append(data)
        except:
            pass
    if vecsymbols==[]:
        vecvar=[]
        for data in vece4:
            try:
                vecv=list(data.free_symbols)
                for sdata in vecv:
                    if not sdata in vecvar:  
                        vecvar.append(sdata)
            except:
                pass
        vecvar.sort(key=str)
        vecsymbols=vecvar
    vec1=[]
    vec2=[]
    for data in vece4:
        if Is_Img(data):
            vec1.append(im(data))
            vec2.append(re(data))
        else:
            vec2.append(data)
    vece5=vec2+vec1        
    return vece5,vecsymbols 
    
def eqpresolve(vecexpr,vecvar,**kwargs):
    if len(kwargs)>0:
        vecval2=[]
        for data in vecexpr:
            kval=real_subs(data,**kwargs)
            vecval2.append(kval)
        vecexpr=vecval2    
    kres=solve(vecexpr,vecvar,dict=True)
    numvar=len(vecvar)
    if type(kres)==dict:
        svar,ssol=list(kres.keys()),list(kres.values())
        snam=[str(data) for data in svar]
        knames=snam
        kvalue=ssol
        qq=1         
    elif type(kres)==list:
        ssvar,sssol,ssnam=[],[],[]
        qq=len(kres)
        cc=1
        for data in kres:
            svar,ssol=list(data.keys()),list(data.values())
            snam=[str(data) for data in svar]
            ssvar=ssvar+svar
            sssol=sssol+ssol
            ssnam=ssnam+snam
            cc+=1
        knames=ssnam
        kvalue=sssol
 
    else:
        knames=[str(data) for data in vecvar]
        kvalue=kres
        qq=1
    return knames,kvalue,qq,numvar 

    
def simplesolve(*args,positive='',**kwargs):
    """
    Solve systems of equations with multiple processing options.
    
    This function provides an advanced interface for solving equations that combines
    symbolic resolution with filtering and formatting options for results.
    
    Parameters
    ----------
    *args : equations and options
        - Equations to solve (Eq objects or expressions)
        - Processing options (see below)
    
    positive : str or list, optional
        Name(s) of variable(s) that must be positive.
        If specified, filters solutions to keep only those where the
        specified variables have positive values.
    
    **kwargs : additional arguments
        Additional parameters passed to underlying solving functions.
    
    Processing Options (included in *args):
    --------------------------------------
    'positive' or 'pos' : Filter solutions to keep only positive values
                          (for single-variable equations)
    
    'float' : Convert solutions to floating-point numbers
    
    'factor' : Factorize the solutions
    
    'simplify' : Simplify solutions (also applies factorization)
    
    'noimg' : Filter out solutions containing imaginary numbers
    
    'fulldata' : Return both variable names and values as (vecname, vecval)
    
    'eQ' : Return results as Eq objects (equalities)
    
    'noshow' or 'noimg' : Suppress automatic display of results
    
    Returns
    -------
    Depending on the options used:
    
    - By default: Returns solution values
      * Single value if only one solution
      * List of values if multiple solutions
    
    - With 'fulldata': Returns tuple (vecname, vecval)
      vecname: List of variable names as strings
      vecval: List of corresponding values
    
    - With 'eQ': Returns list of Eq objects representing the solutions
    
    Raises
    ------
    Various errors from underlying functions:
    - If equations cannot be solved
    - If inconsistent constraints
    - If specified variable not found
    
    Notes
    -----
    1. The function first extracts equations and variables using eQandsymbols()
    2. Then solves them using eqpresolve() with provided kwargs
    3. Applies post-processing filters based on options in *args
    4. Optionally displays formatted results using showresolve2Eq()
    
    Processing order:
    1. Equation extraction (eQandsymbols)
    2. Initial solving (eqpresolve)
    3. Positive variable filtering (if positive parameter provided)
    4. Single-variable positive filtering (if 'positive'/'pos' in args)
    5. Float conversion (if 'float' in args)
    6. Factorization (if 'factor' in args)
    7. Simplification (if 'simplify' in args)
    8. Imaginary number filtering (if 'noimg' in args)
    
    Examples
    --------
    >>> # Basic equation solving
    >>> x, y = symbols('x y')
    >>> simplesolve(Eq(x + y, 5), Eq(x - y, 1))
    [x = 3, y = 2]
    
    >>> # Get solutions as floats
    >>> simplesolve(Eq(x**2, 2), 'float')
    [-1.4142135623730951, 1.4142135623730951]
    
    >>> # Only positive solutions
    >>> simplesolve(Eq(x**2, 4), 'positive')
    [2]
    
    >>> # Get full data (names and values)
    >>> simplesolve(Eq(x + y, 3), Eq(x - y, 1), 'fulldata')
    (['x', 'y'], [2, 1])
    
    >>> # Filter variables that must be positive
    >>> simplesolve(Eq(x + y, 3), Eq(x*y, 2), positive='x')
    Filters solutions where x > 0
    
    See Also
    --------
    eQandsymbols : Extracts equations and variables
    eqpresolve : Core equation solving function
    onlypositives : Filters positive solutions
    showresolve2Eq : Displays formatted solutions
    resolve2Eq : Converts to Eq objects
    
    Dependencies
    ------------
    Requires sympy for symbolic mathematics and several helper functions:
    - eQandsymbols, eqpresolve, onlypositives, resolve2Eq, showresolve2Eq
    - signo (for sign detection)
    - factor, simplify (from sympy)
    """
    vecexpr,vecvar=eQandsymbols(*args)
     
    vecname,vecval,qq,numvar=eqpresolve(vecexpr,vecvar,**kwargs)
    if positive!='':
        vecname,vecval,qq,numvar=onlypositives(positive,vecname,vecval,qq,numvar)
    if ('positive' in args or 'pos' in args) and numvar==1:
        nsol=[]
        for data in vecval:
            if signo(data)==1 or float(data)>=0:
                nsol.append(data)
        vecname=vecname[0:len(nsol)]
        qq=len(nsol)
        vecval=nsol   
    if 'float' in args:
        vecval2=[]
        for data in vecval:
            try:
                vecval2.append(float(data))
            except:
                vecval2.append(data)
        vecval=vecval2
    if 'factor' in args:
        vecval2=[]
        for data in vecval:
            try:
                vecval2.append(factor(data))
            except:
                vecval2.append(data)
        vecval=vecval2
    if 'simplify' in args:
        vecval2=[]
        for data in vecval:
            try:
                vecval2.append(simplify(data))
            except:
                vecval2.append(data)
        vecval=vecval2        
        vecval=[factor(data) for data in vecval]
        
    if 'noimg' in args:
        vecval2=[]
        for data in vecval:
            if not 'I' in str(data):
                vecval2.append(simplify(data))

        vecval=vecval2
    if 'fulldata' in args:
        return vecname,vecval
    elif 'eQ' in args:
        return resolve2Eq(vecname,vecval,qq,numvar)
    else:
        if not 'noshow' in args and not 'noimg' in args:
            showresolve2Eq(vecname,vecval,qq,numvar)
        if len(vecval)==1:
            return vecval[0]
        else:    
            return vecval
def onlypositives(vposi,knames,kvalue,qq,numvar):
    nvecn=[]
    nvecv=[]
    nqq=qq
    nknames=knames[0:numvar]
    nkvalue=[]
    for i in range(qq):
        nvecval=kvalue[numvar*i:numvar*(i+1)]
        done=True
        for tvec,vvec in zip(nknames,nvecval):
            if str(tvec)==vposi and signo(vvec)==-1:
                done=False
        if done:
            nvecn=nvecn+nknames
            nvecv=nvecv+nvecval
        else:
            nqq=nqq-1
    return(nvecn,nvecv,nqq,numvar)
    

        
def resolve2Eq(vecname,vecval,qq,numvar):
    veceq=[]
    if qq==1:
        for kname,kval in zip(vecname,vecval):
            veceq.append(MyEq(kval,kname))
        return veceq
    elif qq>1 and numvar==1:
        cc=1
        for data,sname in zip(vecval,vecname):
            veceq.append(MyEq(data,sname+str(cc)))
            cc+=1
        return veceq    
        
    else:
        cont=0
        for i in range(numvar):
            for ss in range(1,qq+1):
                veceq.append(MyEq(vecval[cont],vecname[cont]+str(i+1)))
                cont+=1
        return veceq

def showresolve2Eq(vecname,vecval,qq,numvar):
    veceq=[]
    if qq==1:
        for kname,kval in zip(vecname,vecval):
            veceq.append(MyEq(kval,kname))
         
    else:
        cont=0
        for i in range(qq):
            for ss in range(numvar):
                try:
                    veceq.append(MyEq(vecval[cont],vecname[ss]+str(i+1)))
                    cont+=1
                except:
                    pass    
                  
def newton_raphson(expr,x0, tol=1e-6, max_iter=100):
    x=symbols('x')
    x=unisymbols(x)
    dexpr=unisymbols(expr.diff(x))
    def f(w):
        return expr.subs(x,w)
    def f_prime(w):
        return dexpr.subs(x,w)
    x1 = x0
    for i in range(max_iter):
        fx = f(x1)
        fpx = f_prime(x1)
        if fpx == 0:
            raise ValueError("Derivada es cero. No se puede continuar con el método de Newton-Raphson.")
        
        x_new = x1 - fx / fpx
        
        # Verificar la convergencia
        if abs(x_new - x1) < tol:
            return x_new
        
        x1 = x_new
    
    raise ValueError("El método de Newton-Raphson no converge después de {} iteraciones.".format(max_iter))                 

def solveaprox(*args,tol=1e-6,niter=100):
    '''
    args[0]= math expr
    args[1,...] = point where will be evaluate
    tol= aproxximation to answer
    niter= number of max loops
    '''
    expr=args[0]
    if len(args)==1:
        vecx0=[0]
    elif len(args)==2:
        vecx0=[args[1]]
    else:
        vecx0=args[1::]
    kres=[newton_raphson(expr,x0,tol=tol, max_iter=niter) for x0 in vecx0]
    kres2=[]
    for data in kres:
        if not data in kres2:
            kres2.append(data)
    return kres2  
    
    
def find_intex1(expr,vecx,qq=100):
    x1,x2=vecx
    def fs(w):
        return expr.subs(x,w)
    X=list(np.arange(x1,x2,(x2-x1)/qq))
    Y=[fs(data) for data in X]
    S=[signo(data) for data in Y]
    So=S[0]
    Xo=X[0]
    vecc=[]
    vecm=[]
    for vs,vx in zip(S[1::],X[1::]):
        if vs!=So:
            vecc.append([Xo,vx])
            vecm.append((Xo+vx)/2)
            So=vs
        Xo=vx 
    return vecc,vecm
def find_intex2(expr,LP,qq=100):
    solu=[]
    for data in LP:
        sL,sP=find_intex1(expr,data,qq=100)
        solu.append(sP[0])
    return solu  
def find_intex(expr,vac,qq=100):
    LP,PP=find_intex1(expr,vac,qq=qq)
    vecsolu=find_intex2(expr,LP,qq=qq)
    return vecsolu   
    
def find_integ1(expr1,expr2,vac):
    x1,x2=vac
    def fs1(w):
        return expr1.subs(x,w)
    def fs2(w):
        return expr2.subs(x,w)
    
    qq=100
    X=list(np.arange(x1,x2,(x2-x1)/qq))
    Y=[fs1(data)-fs2(data) for data in X]
 
    S=[signo(data) for data in Y]
    So=S[0]
    Xo=X[0]
    vecc=[]
    vecm=[]
    for vs,vx in zip(S[1::],X[1::]):
        if vs!=So:
            vecc.append([Xo,vx])
            vecm.append((Xo+vx)/2)
            So=vs
        Xo=vx 

    return vecc,vecm

def find_integ2(expr1,expr2,LP):
    solu=[]
    for data in LP:
        sL,sP=find_integ1(expr1,expr2,data)
        solu.append(sP[0])
    return solu  

def find_integ(expr1,expr2,vac):
    LP,PP=find_integ1(expr1,expr2,vac)
    vecsolu=find_integ2(expr1,expr2,LP)
    return vecsolu

 # finding intersection of funcion whit x axis ,'roots'...
# input:
#    vec_xsolu(fx1,fx2,fx3,...,x1=-1,x2=1)
# output [[x11,x12,..],[x21,x22,,]]respective function:
#   
def vec_xsolu(*args,x1=-1,x2=1):
    vecsolu=[]
    vac=[x1,x2]
    for data in args:
        vecsolu.append(find_roots(data,vac))
    return vecsolu
    
def find_root1(expr,vecx,qq=100):
    x1,x2=vecx
    def fs(w):
        return expr.subs(x,w)
    X=list(np.arange(x1,x2,(x2-x1)/qq))
    Y=[fs(data) for data in X]
    S=[signo(data) for data in Y]
    So=S[0]
    Xo=X[0]
    vecc=[]
    vecm=[]
    for vs,vx in zip(S[1::],X[1::]):
        if vs!=So:
            vecc.append([Xo,vx])
            vecm.append((Xo+vx)/2)
            So=vs
        Xo=vx 
    return vecc,vecm
def find_r2(expr,LP,qq=100):
    solu=[]
    for data in LP:
        lim1,lim2=data
        ylim1=expr.subs(x,lim1)
        ylim2=expr.subs(x,lim1)
        if abs(ylim1)<0.0000001:
            solu.append(vlim1)
        elif abs(ylim2)<0.0000001:
            solu.append(vlim2)
        else:    
            sL,sP=find_r1(expr,data,qq=100)
          
            solu.append(sP[0])
    return solu
def newton_raphson(expr,x0, tol=1e-6, max_iter=100):

    dexpr=unisymbols(expr.diff(x))
    def fu(w):
        return expr.subs(x,w)
    def f_prime(w):
        try:
            kres=dexpr.subs(x,w)
        except:
            kres=dexpr
        return kres    
    x1 = x0
    for i in range(max_iter):
        fx = fu(x1)
        fpx = f_prime(x1)
        if fpx==0:
            return x1
        x_new = x1 - fx / fpx
        
        # Verificar la convergencia
        if abs(x_new - x1) < tol:
            return x_new
        
        x1 = x_new
    
    return x1   
    
def find_roots(expr,vac,*args):
    LP,PP=find_root1(expr,vac)
    vecsolu=[]
    for data in PP:
        vecsolu.append(newton_raphson(expr,data, tol=1e-6, max_iter=100))
     
    return vecsolu
def find_groots1(expr1,expr2,vac,qq=100):
    x1,x2=vac
    expr=expr1-expr2
    def fs(w):
        return expr.subs(x,w)
    X=list(np.arange(x1,x2,(x2-x1)/qq))
    Y=[fs(data) for data in X]
    S=[signo(data) for data in Y]
    So=S[0]
    Xo=X[0]
    vecc=[]
    vecm=[]
    for vs,vx in zip(S[1::],X[1::]):
        if vs!=So:
            vecc.append([Xo,vx])
            vecm.append((Xo+vx)/2)
            So=vs
        Xo=vx 
    return vecc,vecm

def newton_graphson(expr1,expr2,x0, tol=1e-6, max_iter=100):
    expr=expr1-expr2
    dexpr=unisymbols(expr.diff(x))
    def fu(w):
        return expr.subs(x,w)
    def f_prime(w):
        try:
            kres=dexpr.subs(x,w)
        except:
            kres=dexpr
        return kres    
    x1 = x0
    for i in range(max_iter):
        fx = fu(x1)
        fpx = f_prime(x1)
        if fpx==0:
            return x1
        x_new = x1 - fx / fpx
        
        # Verificar la convergencia
        if abs(x_new - x1) < tol:
            return x_new
        
        x1 = x_new
    
    return x1   
    
def find_groots(expr1,expr2,vac,*args):
    LP,PP=find_groots1(expr1,expr2,vac)
    vecsolu=[]
    for data in PP:
        vecsolu.append(newton_graphson(expr2,expr2,data, tol=1e-6, max_iter=100))
     
    return vecsolu 
    
def solvecomplex(*args):
    ''' 
    return simplesolve of expresion including Img numbers and complex numbers
    '''
    eqs=[]
    var=[]
    sexpr=[]
    for data in args:
        if type(data)==MyEqEq:
            expr=data.e1.ksym-data.e2.ksym
            if 'I' in str(expr):
                eqs.append(ipart(expr))
                eqs.append(rpart(expr))
            else:
                eqs.append(expr)
        elif type(data)==MyEq:
            epxr=data.ksym
            if 'I' in str(expr):
                eqs.append(ipart(expr))
                eqs.append(rpart(expr))
            else:
                eqs.append(expr)
        elif Is_symbols(data):
            var.append(data)
        elif type(data)==str:
            sexpr.append(data)
        else:
            eqs.append(data)
    if len(var)>0:
        eqs=eqs+var
    if len(sexpr)>0:
        eqs=eqs+sexpr
    return simplesolve(*eqs)

#####  Geometry
def eQ2points(*args):
    '''
    eQ2points(x1,y1,x2,y2)
    eQ2points((x1,y1),(x2,y2))
    eQ2points(x1,y1,m), m=slope
    eQ2points(x1.....,'name') return MyEq name=name
    '''
    vec=[]
    name=''
    for data in args:
        if type(data)==str:
            name=data
        else:
            vec.append(data)
    args=vec
    if len(args)==4:
        x1=args[0]
        y1=args[1]
        x2=args[2]
        y2=args[3]
        m=cfrac((y2-y1),(x2-x1))
    elif len(args)==3:
        x1=args[0]
        y1=args[1]
        m=args[2]
    else:
        x1,y1=args[0]
        x2,y2=args[1]
        m=cfrac((y2-y1),(x2-x1))
    b=y1-m*x1
    if name!='':
        ee=MyEq(m*x+b,name)
        return ee
    return m*x+b    
    
    
import time

def timestart():
    # see finishtime
    return time.time()
def starttime(): 
    return time.time()
 
    
    
    
def endtime(t):
    return time.time()-t  
def timeend(t):
    return time.time()-t    
def finishtime(t):
    return time.time()-t 

###  FACORIAL REDUCCION
def factorfac(expr):
    '''
    factorfac(x*(x-1)*(x-2)*(x-3)*(x-4)+3) return factorial(x)/factorial(x-5)+3
    '''
    
    if type(expr)==factorial:
       return expr 
    if Is_Add(expr):
        m=expr.args 
        kres=0
        for data in m:
            kres=kres+factorfac(data)
        return kres

    elif Is_Div(expr):
        p1=numer(expr)
        p2=denom(expr)
        return sdiv(factorfac(p1),factorfac(p2))
    
    elif Is_Mul(expr):
        return factorfactorial(expr)
    else:
        return expr

def factorfactorial(expr):
    '''
    return true if expr=(x-a)*(x-a-1)*(x-a-2)*..(x-a-n)
    '''
    var=list(expr.free_symbols)[0]
    kres=[Is_MiniMionomio(data) for data in expr.args]
    if prod(kres)==1:
        vecti=[data.subs(var,0) for data in expr.args]
    else:
        return expr
        
    vecti=sorted(vecti)
    done=True
    qq=len(vecti)
    for i in range(0,qq-1):
        if vecti[i]+1!=vecti[i+1]:
            done=False
    t0=vecti[-1]
    tn=vecti[0]
    kres= factorial(var+t0)/factorial(var+tn-1) 
    return kres

def Is_MiniMionomio(expr): 
    '''
    return True if expr= x,x+1,x+a
    '''
    if type(expr)==Symbol:
        return True
    if Is_Add(expr):
        mm=expr.args
        if len(mm)==2:
            donevar=False
            donenum=False
            for data in mm:
                if type(data)==Symbol:
                    donevar=True   
                if Is_Number(data):
                    donenum=True
            if donevar and donenum:
                return True
            else:
                return False
        return False
    return False


def Line2Eq(*args):
    '''
    from 2 Point (x1,y1),(x2,y2) or Point2D (P1,P2) return m*x+b
    '''
    if type(args[0])==Line2D:
        L=args[0]
 
        Q=MQ(unisymbols(L.equation()),0,show=False)
        Q.alone(y,show=False)
        return Q.R
    else:
        L=Line2D(args[0],args[1]) 
        return Line2Eq(L)   

def despeja(expr,var,**kwargs):
    expr=real_subs(unisymbols(expr),**kwargs)
    Q=MQ(expr,0,show=False)
    Q.alone(var,show=False)
    return Q.R
    
def getterminosxy(Obj):
    ''' 
    Obj=ay+bx+c, MyEq(ay+bx+c), MQ(ay+bx+c,0)
    return (a,b,c) 
    '''
    obj=obj2expr(Obj)
    iexpr=obj.subs(x,0)
    iexpr=iexpr.subs(y,0)
    xexpr=obj-iexpr
    yexpr=obj-iexpr
    xexpr=xexpr.subs(y,0)
    yexpr=yexpr.subs(x,0)
    if yexpr!=0:
        ksig=signo(yexpr)
        if ksig==-1:
            return -1*yexpr,-1*xexpr ,-1*iexpr
        else:
            return yexpr,xexpr,iexpr
    else:
        ksig=signo(xexpr)
        if ksig==-1:
            return 0,-1*xexpr,-1*iexpr  
        else:
            return 0,xexpr,iexpr 
            
        
def find2square(expr,var=x,*args):
    '''
    loop x1=1,2,3,......x2=100
    to check if expr(var=x1) is perfect Square
    expr=x*x+3+4, var=x defauly
    expr=y*Y*7+... var=y
    x1=ini value, x2=end value, 
    x1=1,x2=100 default value
    '''
    done=True
    x1=1
    x2=100
    if len(args)>0:
        x1=args[0]
    if len(args)>1:
        x2=args[1]

    while done:
        kres=expr.subs(var,x1)
        if Is_Square(kres):
            return x1
        elif x1>=x2:
            return False
        else:
            x1+=1        
            
def sumatoria(*args):
    '''
    nice display and get answer to summation(*args)
    '''
    ksum=Sum(*args)
    ssum=ksum.doit()
     
    if type(ssum)==Piecewise:
        kres=ssum.args[1][0]
        mdisplay(ksum,' = ',kres, ssum)
        return kres
    else:
        mdisplay(ksum,' = ',ssum)
        return ssum 

        
        
def simplify_minusone(expr):
    ''' 
    simplifica las expresiones (-1)**x a 1 si pued detectar si x es par o impar
    return 1 si es parreturn -1 si es impar
    '''
    
    if Is_Add(expr):
        mm=0
        for i in expr.args:
            mm=mm+simplify_minusone(i)
        return mm    
    elif Is_Div(expr):
        p1=simplify_minusone(numer(expr))
        p2=simplify_minusone(denom(expr))
        return cf(p1,p2)
    elif type(expr)==Mul:
        mm=expr.args
        kres=simplify_minusone(mm[0])
        for i in mm[1::]:
            kres=kres*simplify_minusone(i)
        return kres
    elif Is_Pow(expr):
        if Is_PowMinusOne(expr):
            ee=getexpo(expr)
            dkres=False
            try:
                if ee%2==0:
                    return 1
            except:
                pass
            if Is_Mul(ee):
                for data in ee.args:
                    try:
                        if data%2==0:
                            return 1
                    except:
                        pass
            return expr            
        return expr
         
    elif Is_Root(expr):
        rr=getroot(expr)
        mm=insideroot(expr)
        return rpow( simplify_minusone(mm),rr)
        
    else:
        return expr           
        
        
def mmodulo(expr, n):
    # Caso base: expr es un entero
    if isinstance(expr, Integer):
        return expr % n
    
    # Caso suma: expr = a + b + ...
    elif isinstance(expr, Add):
        total = 0
        for arg in expr.args:
            total += mmodulo(arg, n)
        return total % n if isinstance(total, Integer) else Add(*[mmodulo(arg, n) for arg in expr.args])
    
    # Caso multiplicación: expr = a * b * ...
    elif isinstance(expr, Mul):
        # Si algún factor es múltiplo de n, el producto es 0 módulo n
        if any(isinstance(arg, Integer) and arg % n == 0 for arg in expr.args):
            return 0
        # Si no, procesamos cada factor recursivamente
        product = 1
        for arg in expr.args:
            arg_mod = mmodulo(arg, n)
            if isinstance(arg_mod, Integer):
                product *= arg_mod
            else:
                return Mul(*[mmodulo(arg, n) for arg in expr.args])
        return product % n
    
    # Caso general (símbolos, potencias, etc.)
    else:
        return expr % n if Is_Number(expr) else expr 

class Msolve:
 
    def __init__(self,*args):
        vecq=[]
        vecv=[]
        for data in args:
            if type(data)==MyEqEq:
                vecq.append(data)
            else:
                vecv.append(data)
 
        
        self.Qs=vecq
        self.Vs=vecv
        self.Qs2=[]
        self.Vsolved=[]
        self.Vvalue=[]
        self.setgrado()
        self.s()
    def setgrado(self):
        for Qx in self.Qs:
            p1,p2=getsymbolsMQ(Qx)
            vecs=list(set(p1+p2))
            Qx.grado=len(vecs) 
                     
    def upsolu(self):
        for Qx in self.Qs:
            if Qx.grado==1:
                self.Qs.remove(Qx)
                self.Qs2.append(Qx)
                
    def presolve(self):
        for Qy in self.Qs2:
            for Qx in self.Qs:
                Qx.set(Qy,show=False)
        self.s()
    def alone(self):
        for Qx in self.Qs2:
            p1,p2=getsymbolsMQ(Qx)
            pp=p1+p2
            Qx.alone(pp[0],show=False)
             
    def cleard(self):
        self.setgrado()
        vecok=[]
        v1=[]
        v2=[]
        for Qx in self.Qs2:
            if Qx.grado>0:
                kres=Qx.L-Qx.R
                if not kres in vecok:
                    vecok.append(kres)
                    v1.append(Qx)
        self.Qs2=v1
        for Qx in self.Qs:
            if Qx.grado>0:
                kres=Qx.L-Qx.R
                if not kres in vecok:
                    vecok.append(kres)
                    v2.append(Qx)
        self.Qs=v2        
         
    def solve(self,Qx,varx):
        Qx.alone(varx,show=False)
        for Qy in self.Qs:
            if Qx!=Qy:
                Qy.set(Qx,show=False)
        self.setgrado()
        self.upsolu()
        self.alone()
        self.presolve()
        self.cleard()
        self.revisa()
         
 
 
    def s(self):
        for Qx in self.Qs2:
            name=Qx.name
            p1=Qx.L
            p2=Qx.R
            display(Math(name+': '+latex(p1)+' = '+latex(p2)+'\ grado: '+latex(Qx.grado)+',\Qv'))
        for Qx in self.Qs:
            name=Qx.name
            p1=Qx.L
            p2=Qx.R
            display(Math(name+': '+latex(p1)+' = '+latex(p2)+',\  grado: '+latex(Qx.grado)+',\Qs'))

    def revisa(self):
        done=True 
        while done:
            done=False
            self.setgrado()
            for Qx in self.Qs:
                if Qx.grado==1:
                    done=True
            if done:
                self.setgrado()
                self.upsolu()
                self.alone()
                self.presolve()
                self.cleard()
            clearcell()
            self.s()
            time.sleep(3)

        clearcell()
        self.s()
            
def getsymbolsMQ(Qq):
    try:
        p1=list(Qq.L.free_symbols)
    except:
        p1=[]
    try:    
        p2=list(Qq.R.free_symbols)
    except: 
        p2=[]
    return p1,p2            