from sympy import *
from lib_Mathbasic  import * 
from lib_Exponencial import * 

 

k,x,t=symbols('k x t')
# Trandformacion de monimios and informacion Funciones
s,v,u,w=symbols('s v u w')


def diffsymbolprime(var):
    svar=str(var)
    dsvar1=svar+"'"
    dsvar2=svar+"''"
    kres=  dsvar1+' '+dsvar2
     
    return symbols(kres)

def taylor(expr, var=x, x0=0, order=5):
    """
    Devuelve el polinomio de Taylor de una expresión simbólica, ordenado en forma creciente por grado.

    Parámetros:
    - expr : Expr
        La expresión simbólica a expandir.
    - var : Symbol
        La variable respecto a la cual expandir.
    - x0 : número o símbolo (opcional)
        El punto alrededor del cual hacer la expansión (por defecto, 0).
    - order : int (opcional)
        Orden del desarrollo (por defecto, 6).

    Retorna:
    - Expr
        El polinomio de Taylor truncado, ordenado por grados crecientes.
    """
    expansion = expr.series(var, x0, order).removeO()
    terms = expansion.as_ordered_terms(order='lex')  # ordena por potencia creciente
    kres=Add(*terms)
    m=list(kres.args)
    return order_srender(m)


def exp_in(ksym): # retorna True si hay una funcion exp en ksym 
                  #  exp_in(4*exp(m)) return True
                  #  exp_in(4*m) return False
    sres=str(ksym)
    if 'exp' in sres:
        return True
    return False
    
def Is_Exp(ksym):
    return exp_in(ksym=ksym)

################################################

    
def get_exp_Mono(ksym):       # retorna la parte exp de un monomio
    if exp_in(ksym):          # get_exp_Mono(4*exp(x)) return exp(x)
        if Is_Mul(ksym):       
            mm=fpoly(ksym,'list')
            for i in mm:
                if exp_in(i):
                    return i
    return ksym  

def get_Euler (ksym):
    return get_cofac_exp(ksym=ksym) 

##########################################    
    
def get_cofac_exp(ksym):      # retorna la parte  no exp de un monomio
    kexp=get_exp_Mono(ksym)   # get_cofac_exp (4*exp(x)) return 4
    if kexp!=ksym:
        return simplify(ksym/kexp)
    return 1 

def get_facEuler(ksym):
    return get_cofac_exp(ksym=ksym)

####################################


def get_expoinexp(ksym):
    kres=ksym
    if exp_in(ksym):
        mm=get_exp_Mono(ksym)
        mm2=fpoly(mm,'list')
        kres=mm2[0]
        return kres
    return 1    

def get_expnEuler(ksym):
    return get_expoinexp(ksym=ksym)


#  LOG

def log_in(ksym): # retorna True si hay una funcion log en ksym 
                  #  log_in(4*exp(m)) return False
    sres=str(ksym)              #  log_in(log(4)*m) return True
    if 'log' in sres:
        return True
    return False 

def get_log_Mono(ksym):       # retorna la parte exp de un monomio
    if log_in(ksym):          # get_exp_Mono(4*exp(x)) return exp(x)
        if Is_Mul(ksym):       
            mm=fpoly(ksym,'list')
            for i in mm:
                if log_in(i):
                    return i
    return ksym

def get_cofac_log(ksym):      # retorna la parte  no exp de un monomio
    kexp=get_log_Mono(ksym)   # get_cofac_exp (4*exp(x)) return 4
    if kexp!=ksym:
        return simplify(ksym/kexp)
    return ksym 

def get_insidelog(ksym):      #  
    kexp=get_log_Mono(ksym)   #  
    return simplify(exp(kexp)) 

#####################################
#           Monomie Tools 
#####################################
    
def monodata(*args):
    r'''
       data_mono(monomie, op1, op2)
        inputs args:
           op1=='isexp' ,return True if ksym is monomie whit exp factor
           op1=='getexp' ,return only exp factor 
           op1=='getexpfac' ,return all factors less exp fac

    '''
    ksym=args[0]
    if len(args)==2:
        op1=args[1]
    if len(args)==3:
        op2=args[1]
     
    # Exp
    if op1=='isexp':
        return exp_in(ksym)

    elif op1=='getexp':
        return get_exp_Mono(ksym)
    elif op1=='getexpfac':
        return get_cofac_exp(ksym)
    elif op1=='getexpoinexp':
        return get_expoinexp(ksym)

    # Log
    if op1=='islog':
        return log_in(ksym)

    elif op1=='getlog':
        return get_log_Mono(ksym)
    elif op1=='getlogfac':
        return get_cofac_log(ksym)
    elif op1=='getinlog':
        return get_insidelog(ksym)
        
        
def clear_exp_QQ(ksym1,ksym2,vmain=k,var2=t):
    p2=simplify(ksym2)
    p1=simplify(ksym1)
    if Is_Add(p2) and Is_Poly(p2):
        mm=0
        mmexp=0
        for i in fpoly(p2,'list'):
            if exp_in(i):
                mmexp+=i
            else:
                mm+=i
        p2=mmexp
        p1=p1-mm
        if Is_Mono(p2) and denom(p2)!=1:
            kdeno=denom(p2)
            p1=p1*kdeno
            p2=simplify(p2*kdeno)
            
        mm=get_facEuler(p2)
        if mm!=1:
            p1=p1/mm
            p2=simplify(p2/mm)
            
            
        return  (p1,p2) 

def getexponent(ksym): # return only exponente , not sqrt  
  
    kres=ksym
    ssym=str(ksym)
    if 'sqrt' in ssym:
        if 'sqrt'==ssym[0:4]:
            nssym=between_par(ssym,'sqrt')
            kres=parse_expr(nssym)
    if Is_Pow(kres):
        mm=fpoly(kres,'list')
        return mm[1]
    else:
        return 1
        
'''
def getexpo(ksym):
    return getexponent(ksym) 
     
def getbase(ksym):  # return base from Pow monomie expresion
     
        input: x**a
        return : x  
     
    ksym=ksym*1
    if Is_Pow(ksym):
        kres=ksym
        ssym=str(ksym)
        if 'sqrt' in ssym:
            if 'sqrt'==ssym[0:4]:
                nssym=between_par(ssym,'sqrt')
                kres=parse_expr(nssym)
                return getbase(kres)
         
        mm=fpoly(ksym,'list')
        return mm[0]
         
            
    else:
        return ksym
'''
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

def sqrt2pow(ksym): # convert sqrt(x**y) in x**(y/2)
    
    b1,e1,r1=partPow(ksym)
    return kpow(b1,cfrac(e1,r1))
    
def root2pow(ksym): # convert sqrt(x**y) in x**(y/2)
    
    b1,e1,r1=partPow(ksym)
    return kpow(b1,cfrac(e1,r1))    
    
def pow2sqrt(ksym): # convert (x**(y/2)) in sqrt(x**y)
    b1,e1,r1=partPow(ksym)
    if r1==2:
        return sqrt(kpow(b1,e1))
    else:
        return ksym

    
    b1,e1,r1=partPow(ksym)
    return kpow(e1,cfrac(e1,1)) 
def modulores(ksym,rr=2): # function use by reduceroot
    b,ee,r1=partPow(ksym)
     
    if denom(ee)!=1:
        nn=numer(ee)
        dd=denom(ee)
        p2=nn%rr
        p1=int(nn/rr)
        p3=cfrac(p2,dd)
        return(b,p1,p3)
    else:
        nn=ee
        dd=r1
        p2=nn%rr
        p1=int(nn/rr)
        p3=cfrac(p2,r1)
        return(b,p1,p3)

def reduceroot(ksym2,rr=2):
      
    P1=1
    P2=1
    ksym=kpow(ksym2,rr) 
    if Is_Mul(ksym):
        mm=fpoly(ksym,'list')
        for i in mm:
            if Is_Pow(i):
                b,p1,p2=modulores(i,rr)
                
                P1=P1*kpow(b,p1)
                P2=P2*kpow(b,p2)
            else:
                P2=P2*i
        P3=rpow(P2,rr)
        return Mul(P1, P3, evaluate=False)     
    else:
        return ksym
def simplifyroot(ksym):
     
    ksym=1*ksym 
    if Is_Pow(ksym):
         
        kres= root2pow(ksym)
        return kres
        
    elif Is_Add(ksym):
        kres=0
        for i in fpoly(ksym,'list'):
            kres+=simplifyroot(i)
        return kres
    elif Is_Mul(ksym):
        kres=1
        for i in fpoly(ksym,'list'):
            kres=kres*simplifyroot(i)
        return kres
    else:
        return ksym

    
#####################################
#             MyEqEq 
##################################### 


        


#####################################
#             str
#####################################
from sympy import parse_expr

def insiderfunc(expr, sexpr):
    """
    Extrae una subexpresión dentro de una cadena que representa una expresión matemática y la convierte en una expresión simbólica de SymPy.

    Parámetros:
    -----------
    expr : str
        Expresión principal de donde se desea extraer una subexpresión.
    sexpr : str
        Subexpresión inicial a buscar dentro de `expr`.

    Retorno:
    --------
    sympy.Expr
        La subexpresión completa encontrada en `expr`, convertida en una expresión simbólica de SymPy si es posible.

    Descripción:
    ------------
    - Convierte la expresión `expr` en una cadena (`smain`).
    - Busca la posición de `sexpr` en `smain` y determina el índice final.
    - Si `sexpr` no contiene paréntesis, se extiende hasta encontrar un espacio, paréntesis de cierre o el final de la cadena.
    - Si `sexpr` contiene paréntesis, se extiende hasta que los paréntesis estén equilibrados.
    - Finalmente, se usa `parse_expr(kres)` para convertir la cadena resultante en una expresión de SymPy.

    Ejemplo de uso:
    ---------------
    ```python
    from sympy import symbols

    expr = "(x + y) * z"
    sexpr = "(x + y)"
    result = insiderfunction(expr, sexpr)
    print(result)  # Convertirá "(x + y)" en una expresión de SymPy
    ```
    """
    smain = str(expr)
    p1 = smain.find(sexpr)
    p2 = p1 + len(sexpr)
    kres = sexpr
    balance = 0

    if "(" not in sexpr:
        done = True
        while done:
            if p2 == len(smain) or smain[p2] == " " or smain[p2] == ")":
                return parse_expr(kres)
            else:
                kres += smain[p2]
                p2 += 1
    else:
        balance = 1
        done = True
        while done:
            if p2 < len(smain):
                kres += smain[p2]
                if smain[p2] == "(":
                    balance += 1
                elif smain[p2] == ")":
                    balance -= 1
                if balance == 0:
                    return parse_expr(kres)
                else:
                    p2 += 1
            else:
                return parse_expr(kres)

def fix_mulpar(sexp):
    vec1='0123456789'
    vec2='abcdefghijklmnopqrstuvwxyz'
    vec3='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    p1=sexp.find('(')
    if p1==-1 or p1==0:
        return sexp
    elif sexp[p1-1] in vec1 or sexp[p1-1] in vec2 or sexp[p1-1] in vec3:
        sexp=sexp[0:p1]+'*'+sexp[p1:len(sexp)]
        return sexp
    else:
        return sexp        
def between_par(sall,ssym=''):  # between_par('***(12(**)3)***') return '12(**)3'
    if ssym!='':
        ssym=fix_mulpar(ssym)
    mm=str(sall)
    if ssym=='':
        ssym='('
    q1=mm.find(ssym)
    if q1==-1:
        return ''
    mm2=mm[q1::]
    
    q2=mm2.find('(')
    mm3=mm2[q2::]
    
    cc=1
    sres=''
    p=1
    i=1
    while cc>0:
        v=mm3[i]
        if v=='(':
            cc+=1
        if v==')':
            cc-=1
        if cc>0:    
            sres+=v
        i+=1
    return parse_expr(sres) 

def between_func(ksym,sfunc):
    '''
        input a function like 
            ksym=sqrt(y + 4)*sqrt((x - 2)**2) type: symbols o function
            sfunc='sqrt(' type=str
        return:
            vec ['y + 4', '(x - 2)**2'] ,all inside sqrt(*****) in tihs case
    
    '''
    svar=str(ksym)
    qq=svar.count(sfunc)
    vecf=[]
    posx=0
    for i in range(qq):
        ff=between_par(svar,sfunc)
        vecf.append(ff)
        posx=svar.find(sfunc)
        svar=svar[posx+1::]
    return vecf 

def fix_str_math_space(ssym):
    ''' 
        input:'3*x*(y+1)'
        return 3*x*(y + 1)'
    '''
    
    kres=''
    cc=0
    for i in ssym:
        if cc==0:
            kres=kres+i
        else:
            if i=='+' or i=='-':
                 
                if ssym[cc-1]!=' ':
                    kres+=' '
                kres+=i
                if ssym[cc+1]!=' ':
                    kres+=' '
            else:
                kres+=i
        cc+=1
    return kres

def sgetmonofactor(ksym,s1,op=''):
    ''' 
       return valid math expresion when convert math expr in str(expr)
       example:
            A=x*y*(1+a)*sqrt(z+1)+20*z
            sgetmonofactor(A,'x*y') return  'x*y*(1 + a)'
            sgetmonofactor(A,'x*y','f') return  'x*y*(1 + a)*sqrt(z + 1)'
    '''
    
    if type(ksym)!=str:
        smm=str(ksym)
    else:
        smm=ksym
    smm=fix_str_math_space(smm) 
    done=True
    pstar=False
    if type(op)==int:
        nsmm=smm[op::]
        p1=nsmm.find(s1)
        if p1==-1:
            return ksym
        p1=p1+op
    else:    
        p1=smm.find(s1)
    if p1==-1:
        return ksym
    sres=s1
    qq=len(s1)
    cc=0
    i=p1+qq
    while done:
        vf=smm[i]
        sres+=vf
        if vf=='(':
            cc+=1
            pstar=True
        if vf==')':
            cc+=-1
        if pstar and cc==0:
            done=False

        i+=1
    if op!='' and type(op)==str:
        if smm[i+1]!='':
            sres=sgetmonofactor(ksym,sres)
        
    return sres 
#####################################
#     polynomie Tools
#####################################

def getmonoexpr(ksym,sval,op=''):
    ''' 
       return valid math expresion that include sval expr
       example:
            A=x*y*(1+a)*sqrt(z+1)+20*z
            getmonoexpr(A,'x*y') return  x*y*(1 + a)
            getmonoexpr(A,'x*y','f') return  x*y*(1+a)*sqrt(z+1)
    '''
    skres=sgetmonofactor(ksym,sval,op=op)
    try:
        kres=parse_expr(skres)
        return kres
    except:
        return ksym
        
        
#   squaresum()
#   -----------
def squaresum(ksym,p1=0,p2=0):
    '''
    sqauresum(ksym,p1=x,p2=1)
        input :
            ksym = initial function to trasform
            p1,p2  = (p1+p2)**2 new formate to get convert
        return:
            (p1+p2)**2 + K ..such that (p1+p2)**2 + K = ksym
        example:
            the polynomie 4*x*x+8*x-6 tray transfor to (2*x+2)**2 +K
            
                squaresum(4*x*x+8*x-6,2*x,2) return (2*x + 2)**2 - 10
    '''
 
    vec1=expand((p1+p2)**2)
    vec2=ksym-vec1
    kres=kpow((p1+p2),2)+vec2
    return kres


def parcial_fraction(x,L,a,b,c,d):
    r'''
    input L/((a*x+b)*(c*x+d))
    
    return  A/(a*x+b) + B/(c*x+d)
    '''
    A=L*a/(d*(a - c))
    B=-L*c/(d*(a - c))
    p1=(a*x+b)
    p2=(a*x+b)
    
    kres=A/p1+B/p2
    return kres
    
def par_frac(ksym,var=''):  
    r'''
    FRACCION PARCIAL
    input C/[A+B]
    return if it possible C1/A + C2/A que es igula a  C/[A+B]
    '''
   
    
    kres=ksym
    
    if Is_Mono(kres) and Is_Mul(kres):
         
        p1=numer(kres) 
        done=False
        if var!='':
            X=symbols('X')
            kres=ksym.subs(var,X)
            done=True
             
         
        p2=denom(kres)
         

        mm= fpoly(p2,'list')
        
        qq=len(mm)
        A1,A2=symbols('A1 A2')
        f1=A1*mm[0]
        f2=A2*mm[1]
        kres2=simplify(expand(f1+f2))
        kres3=factorSec(kres2,X)
         
        A22=solve(parse_expr(between_par(str(kres3),'X*(')),A2)[0]
        kres4=kres3.subs(A2,A22)
        r1=solve(1-kres4,A1)[0]
        r2=A22.subs(A1,r1)
        dd1=mm[0]
        dd2=mm[1]
         
        if done:
            dd1=dd1.subs(X,var)
            dd2=dd2.subs(X,var)
         
        fac1=p1*r1/dd1
        fac2=p1*r2/dd2
        return fac1+fac2
 
    else: ksym 

 
def real_subs(QQ,**kwargs):
    ''' 
    QQ= symbols function
    ** kwargs c7=6,g4=z..etc..
    RETURN real substitucion when variable have underscore name like 'c_7' 'g_4'    
    '''
    
    vvar=fpoly(QQ,'free')
    mvar=[]
    for i in vvar:
        nname=str(i)
        nname=nname.replace('_','')
        nvar=symbols(nname)
        mvar.append(nvar)
         
    kres=QQ 
    sres=str(kres)
    sres=sres.replace('_','')
    kres=parse_expr(sres)
    mkey=[]
    vvalue=[]
    for key, value in kwargs.items():
        mkey.append(key)
        vvalue.append(value)
        
        kres=kres.subs(parse_expr(key),value)
    for i,j in zip(mvar,vvar):
        kres=kres.subs(i,j)
    return (kres)

def coef_list(sexpr,var=x):  # return coef list complete according max degree
        '''
        ee=x*x*a+x*b+c
        return [a,b,c]

        '''
        
        p = Poly(sexpr, var)
        return  p.all_coeffs()
   
        

def sortdegree(ksym,var2=x):
    '''
      group polynomies respect main varible
    '''
    if str(var2) not in str(ksym):
        return ksym
    mm=coef_list(ksym,var2)
    kres=0
    qq=len(mm)
    cc=qq-1
    for i in mm:
        kres+=i*var2**cc
        cc+=-1
    return kres

def solve_compare_coef(*args):
    '''
        input (Lista de coeficientes del poly1,Lista de coeficientes del poly2)
        return solucion al comparar vlores
        ejm:
            vec1=[x1+3*x2,x2-2*x2]
            vec2=[2,3]
            devuelve la sol de x1 y x2 en x1+3*x2=1 y x2-2*x2=3
            opcion:
                solve_compare_coef([x1+3*x2,x2-2*x2],[2,3],[x2,x1])
    '''            
    Vec1=args[0]
    Vec2=args[1]
    Lv=''
    if len(args)==3:
        Lv=args[2]
    vecQ=[]
    Sv=0
    for i,j in zip(Vec1,Vec2):
        vecQ.append(Eq(i,j))
        Sv+=i
    if Lv=='':
        Lv=get_symbols(Vec1)
    solu=solve(vecQ,Lv)
    try:
        key,value=unpack(solu)
        return value
    except:
        return solu
def max_degree(*args,gen=x):
    '''
    input (x*x+x+1,x*8+,5)
    output 2
    '''
    dmax=0
    if len(args)==1 and type(args[0])==list:
        args=args[0]
    for i in args:
        k=degree(i,gen=gen)
        if k>dmax:
            dmax=k
    return dmax 
def algebra_compare_coef(*args,var2=x):
    '''
    function that resuelve variables when compare 2 polinomies wits same grade and same coef    
    input= args:
            function,or,symbols
            gen= function with defaul gen=x
    output: solution all variables possibles
        
            
    '''
    p1=''
    p2=''
    vsym=[]
    for i in args:
        if Is_Symbol(i):
                vsym.append(i)
        else:
            if p1=='':
                p1=i 
            else:
                p2=i        
                
    
     
    ql1=coef_list(p1,var2=var2)
    ql2=coef_list(p2,var2=var2)
    qq=[]
    for i,j in zip(ql1,ql2):
        qq.append(i-j)
    kres=solve(qq,vsym)    
    if type(kres)==dict:
        kres= kunpakDic(kres)
        nvar=kres[0]
        nval=kres[1]
    if type(kres)==list:
        kres= kunpakDic(kres)
        nvar=kres[0]
        nval =kres[1]
 
    return nval



    
#####################################
#           functions 
#####################################    
def subsnumber(ksym,val1,val2): # replace an integer number by symbol 
    kres=ksym
    sres=str(kres)
    sv1=str(val1)
    sv2=str(val2)
    sres=sres.replace(sv1,sv2)
    kres2=parse_expr(sres)
    return kres2
    
def func2symbol(ksym,y,x):
    veri='('+str(x)+')'
    sf=str(y)
    if veri in sf:
        nf=sf.replace(veri,'')
    else:
        nf=sf
    oldf=nf+veri
    sres=str(ksym)
    sres=sres.replace(oldf,nf)
    return parse_expr(sres)
    
def unfog(expr,*args):
    if type(expr)==MyEq:
        expr=expr.ksym
    svec=[]
    vvec=[]
    for i in args:
        if type(i)==str:
            svec.append(i)
        else:
            vvec.append(i)
    kres=findFunc(expr,svec,vvec)
    return kres   
#def findFunc(f1,vecf,vecv):
def findFunc(f1,vecf,vecv):
    '''
     f((x+1)**2)= 2*x*x+4*x+5 , f1= 2*x+3
     procces, x--> Add(1) -->['A'],[1]  -->x+1  -->Pow(2)  -->['A','P'],[1,2]  -->  
     vecf = ['A','P'] ,Posibles...>add,'A',mul  'M',divide 'D',sqrt 'R',pow 'P'
     vecv = ['A','P']
     findFunc(2*x*x+4*x+5,['A','P'],[1,2])
         
     return f(x)= 2*x+3   
    '''
    y,z=symbols('y z')
    #Q=MyEqEq(y,f1)
    xx=solve(y-f1,x)
    if type(xx)==list:
        xx=xx[0]
    for i,j in zip(vecf,vecv):
        if i=='R':
            xx=sqrs(xx,j)
        elif i=='M':
            xx=xx*j
        elif i=='D':
            xx=xx/j
        elif i=='P':
            xx=xx**j
        else:
            xx+=j
    yy=solve(z-xx,y)
    if type(yy)==list:
        yy=yy[0]
    y=yy.subs(z,x)

    
    return y    




#####################################
#           monomies tools
#####################################

def baselist(ksym): # return list of all monomies in Monomie in Mul expr
    '''
        input (x+a)**a*(y-2)**b*z**c
        return [(x+a),(y-2),z]
    '''    
    if Is_Mul(ksym):
        mm=fpoly(ksym,'list')
        vres=[]
        for i in mm:
            if Is_Pow(i):
                vres.append(getbase(i))
            else:
                vres.append(i)
        return vres        
    else:
        return ksym
def expolist(ksym): # return list of each exponents  in Monomie from Mul expr
    '''
        input (x+a)**a*(y-2)**b*z**c
        return [(a,b,c]
    ''' 
    if Is_Mul(ksym):
        mm=fpoly(ksym,'list')
        vres=[]
        for i in mm:
            if Is_Pow(i):
                vres.append(getexponent(i))
            else:
                vres.append(1)
        return vres        
    else:
        return ksym

def expandfracpart(expr,gen=x):
    P1=numer(expr)
    P2=denom(expr)
    Ng=max_degree(P1,gen=gen)
    Dg=max_degree(P2,gen=gen)
    if Ng>Dg:
        Q=quo(P1,P2,gen)
        remain=rem(P1,P2,gen)
        mm=[Q,remain/P2]
        return mm
    else:
        return [0,expr]
        


def powlist(expr):
    if Is_Pow(expr):
        b,ee,rr=partPow(expr)
        mm=[]
        for i in range(ee):
            mm.append(b**(i+1))
        return mm
    else:
        return [expr]
        
def factorlist3(expr):
    '''
    input((x+1)(x-1))
    output[x+1,x-1]
    '''
    expr=factor(expr)
 
    if Is_Pow(expr):
        return powlist(expr)
    else:
        mm=fpoly(expr,'list')
        kres=[]    
        for i in mm:
            if Is_Mul(i):
                mm1= fpoly(i,'list')
                kres=kres+mm1
            elif Is_Pow(i):
                kres=kres+ powlist(i)
            else:
                kres.append(i)
        return kres        
'''
def partialfraction(*args): # Return divition monimies in partial fraction 
    y1,y2,y3,y4,y5,y6,y7,y8,y9,y10=symbols('y1 y2  y3 y4 y5 y6 y7 y8 y9 y10')
    y11,y12,y13,y14,y15,y16=symbols('y11 y12 y13 y14 y15 y16')
    fabr=[y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16]
    expr=args[0]
    if len(args)==2:
        var=args[1]
    else:
        var=x
    expr=factor(expr)
    p1=numer(expr)
    p2=denom(expr)
    vec2=complete2partfrac(p2)
    vec1=vec_daughter(vec2,fabr,gen=var)
    P=0
    for i,j in zip(vec1,vec2):
        P=P+cfrac(i,j)
    P2=factor(P)
    P3=numer(P2)
    P4=expand(P3)
    dg1=degree(p1,gen=var)
    dg2=degree(P4,gen=var)
    if dg1<dg2:
        vec1=coef_list(p1,size=dg2,var2=var)
        vec2=coef_list(P4,var2=var)
    elif dg1>dg2:
        vec1=coef_list(p1,var2=var)
        vec2=coef_list(P4,size=dg1,var2=var)
    else: 
        vec1=coef_list(p1,var2=var)
        vec2=coef_list(P4,var2=var)
        
    vec3=[vec2[i]-vec1[i] for i in range(3)]
    kres=solve(vec3)
    svar,vval=ganswer(kres,'varname','value')
    for i,j in zip(svar,vval):
        P=P.subs(i,j)
    return P
'''
def swaplist(klist):
    qq=len(klist)
    mm=[]
    for i in range(qq):
        mm.append(klist[qq-1-i])
    return mm  
    
#####################################
#           polynomies tools
#####################################

#   simplifysum()
#   -----------        
'''
def simplifysum(pval):
  
    return sum of each monomie of pval but uniformatin each one to get real sum
    

   
    if Is_Add(pval):
        mm=fpoly(pval,'list')
     
        kk=0
        for i in mm:
            k1=simplify(i)
            kk+=mulexp(k1)
        return kk
    else:
        return pval
'''
#   linfactor()
#   -----------        
def linfactor(ksym,kvar=''):
    '''
    factorize() return personality factorization of Polynomie summ 
    example:
        if kfunc = a*x*x+b*c*x*x+3*x
        factorize(kfunc,x*x) 
        return :
        x*x*(x+b*x)+3*x    

    '''

    p1=ksym
    if kvar!='' and type(p1)==Add:
        oldv=0
        newv=0
        kfac=kvar
        mm=fpoly(p1,'list')
        for i in mm:
            vec=simplify(i/kfac)
            if denom(vec)==1:
                newv+=vec
            else:
                oldv+=i
        kres=oldv+newv*kfac
        return kres
    else:
        return ksym



 
def Area2funcXY(ee1, ee2, var=x, cota='', kope=''):
    if type(ee1)==MyEq:
        kres1=ee1.ksym
    else:
        kres1=ee1    
    if type(ee2)==MyEq:
        kres2=ee2.ksym
    else:
        kres2=ee2    
    kres=solve(kres1-kres2,var)
    if cota=='' and len(kres)==2:
        segx=kres
    else:
        x1,x2=cota
        segx=[x1]
        for i in kres:
            if i>x1 and i<x2:
                segx.append(i)
        segx.append(x2)
        
    vecf=[]
    vecps=0
    qq=len(segx)
    for i in range(qq-1):
        ee=MyEq(ee1-ee2,'',x1=segx[i],x2=segx[i+1],var2=var,ktype='I',show=False)
        vecf.append(ee)
        vecps+=ee 
     
    ee2=MyEq(vecps,'Area')
    kres=0
    
    NoSYMB=False
    for i in vecf:
        if Is_Symbol(i.ksym):
            NoSYMB=False
            
    
    for i in vecf:
        value=i.ksym
        if NoSYMB:
            kres+=abs(value.doit())
        else:
            kres+=value.doit()
    ee3=MyEq(kres,'Area')     
    
def get_rem(kres1,kres2,kname='',var2=''): # return rest of div e1/e2
        if type(kres1)==MyEq:
            kres1=kres1.ksym
        if type(kres2)==MyEq:
            kres2=kres2.ksym
        kres=  rem(kres1,kres2)  
        if kname!='':
            
            return MyEq(kres,kname=kname,var2=var2) 
        else:    
            return kres

def get_quo(kres1,kres2,kname='',var2=''): # return quotient  of div e1/e2
        if type(kres1)==MyEq:
            kres1=kres1.ksym
        if type(kres2)==MyEq:
            kres2=kres2.ksym
        kres= quo(kres1,kres2) 
        if kname!='':
            return MyEq(kres,kname=kname,var2=var2) 
        else:    
            return kres
def get_cofactors(kres1,kname='',var2=''):  # return cofactors  of div e1/e2
        if type(kres1)==MyEq :
            kres1=kres1.ksym
        if type(kres2)==MyEq:
            kres2=kres2.ksym
        kres= cofactors(kres1,kres2) 
         
        if kname!='':
            return MyEq(kres,kname=kname,var2=var2) 
        else:    
            return kres
            
            
def get_seudocofactor(e2,e3,var2):
    '''
    tyr to get polynomie complement and multipli complete degree secuence
    e1=x*x-2
    e2=x+1
        if e1=e2*k
        then k=a*x+b
        maybe x*x-2= (x+1)*(a*x+b)
        
    return (a*x+b),[a,b]
    '''
    vecvar=[a,b,c,d]
    vecfvar=[]
    ee2=e2
    if type(e2)==MyEq:
        ee2=e2.ksym
        
    ee3=e3
    if type(e3)==MyEq:
        ee3=e3.ksym    
    qq=degree(ee2,gen=var2)-degree(ee3,gen=var2)
    vres=0
    cc=0
    for i in range(qq+1):
        vres+=vecvar[i]*var2**(qq-cc)
        vecfvar.append(vecvar[i])
        cc+=1
    return (vres,vecfvar)
    
# ---------  transformada 
   
def transformada(Q1,expr0,ssym,kope=''): 
    '''
    # return Q1 in format expr0 ,
    example:
         we want transform 2*x-y to 2(x+y) -3*y ,
        Q1= 2*x-y
        expr0= 'A+2(x+y)'
        ssym ='A'
        kope = 'e', help to simplify ..
        
        transformada(Q1,expr0,ssym,kope='e') return 2(x+y) -3*y
    '''
    S=symbols(ssym)
    
    expr0=expr0.replace(ssym,str(S))
    expr1=simplify(parse_expr(expr0))
    expr2=Q1-expr1
    
    if kope!='':
        for i in kope:
            expr2=opemat(expr2,i)
    kres=csolve(expr2,S)        
    expr3=expr1.subs(S,kres)
    return expr3    
 
def transformada2(expr1,epxr2,sexpr1,sA):
    p1=expr1-epxr2
    tr2=transformada(epxr2,sexpr1,sA)
    return p1+tr2

#####################################
#           geometry tools
#####################################    
def acircle(R,kname=''):
    if type(R)==MyEq:
        R=R.ksym
    kres=pi*R*R    
    if kname!='':
        return MyEq(kres,kname)
    return kres

def dvolume(L1,L2, kname=''):
    if type(L1)==MyEq:
            L1=L1.ksym 
    if type(L2)==MyEq:
            L2=L2.ksym
    kres=L1*L2    
    if kname!='':
        return MyEq(kres,kname)
    return kres 

def mpoint(L1,L2, kname=''):
    if type(L1)==MyEq:
            L1=L1.ksym 
    if type(L2)==MyEq:
            L2=L2.ksym
    kres=L1/L2    
    if kname!='':
        return MyEq(kres,kname)
    return kres     
    
#####################################
#           solver
#####################################     
    
  
#####################################
#           vectores
#####################################

def checvec(vv,qq): # if vv=symbol o number  then return  [vv,vv,vv.... ] qq times ,is not ret vv
    if type(vv)!=list:
        return mzero(qq,vv)
    else:
        return vv
        
def klen(ksym): # if type(ksym)=list then return len(ksym) if not return 1
    if type(ksym)!=list:
        return 1
    else:
        return len(ksym)
        
def addvec(v1,v2): # input v1=[a,b],,v2=[c,d]  return [a+c,b+d], v1=[a,b],,v2=3 return [a+3,b+3]
    v1=checvec(v1,klen(v2))
    v2=checvec(v2,klen(v1))
 
    kres=[]
    for i,j in zip(v1,v2):
        kres.append(i+j)
    return kres

def restvec(v1,v2):  # input v1=[a,b],,v2=[c,d]  return [a-c,b-d], v1=[a,b],,v2=3 return [a-3,b-3]
    v1=checvec(v1,klen(v2))
    v2=checvec(v2,klen(v1))
 
    kres=[]
    for i,j in zip(v1,v2):
        kres.append(i-j)
    return kres 

def mulvec(v1,v2): # input v1=[a,b],,v2=[c,d]  return [a*c,b*d], v1=[a,b],,v2=3 return [a*3,b*3]
    v1=checvec(v1,klen(v2))
    v2=checvec(v2,klen(v1))
 
    kres=[]
    for i,j in zip(v1,v2):
        kres.append(i*j)
    return kres   
    
def unionvec(*args):
    kres=[]
    for i in args:
        for j in i:
            if j not in kres:
                kres.append(j)
    return kres 
    
def getindex(ksym,vec): # getindex(c,[a,b,c,d]) return 2
    return vec.index(ksym) 

#####################################
#           list
#####################################  

def mulitem(ksym):
    '''
    input[a,b,c...]
    return a*b*c*...
    '''
    kres=1
    for i in ksym:
        kres=kres*i
    return kres 
def sumitem(ksym):
    '''
    input[a,b,c...]
    return a+b+c+...
    '''
    kres=1
    for i in ksym:
        kres=kres*i
    return kres
 

def one_ksym2ksym(ksym):  
    r'''
    return denominator if numerator =1
    input 1/(a+b) return (a+b)
    input 3/(a+b) return  3/(a+b)
    input c/(a+b) return c/(a+b)
    ''' 
    
    if numer(ksym)==1:
        return denom(ksym)
    else:
        return ksym

def num_ksym2ksym(ksym):
    r'''
    return denominator if numerator = numeric
    input 1/(a+b) return (a+b)
    input 3/(a+b) return  (a+b)
    input c/(a+b) return c/(a+b)
    '''
    if Is_Number(numer(ksym)):
        return denom(ksym)
    else:
        return ksym    
        
#####################################
#           simplify  Sqrt(Pow(2))
#####################################          
def simplifyrpow(ksym):
    r'''
     
    input sqrt(y + 4)*sqrt((x - 2)**2)
    return  (x - 2)*sqrt(y + 4)
     
    ''' 
    
    svec= between_func(ksym,'sqrt(')
    svar=str(ksym)
    kres=ksym
    if svec!=[]:
        qq=len(svec)
        svar=str(ksym)
        for i in svec:
            if i[-3::]=='**2':
                p1='sqrt('+i+')'
                p2='('+i[0:-3]+')'
                svar=svar.replace(p1,p2)
                kres=parse_expr(svar)
    return kres 

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

def factoriza(expr,p1,*args):

    ''' 
    factoriza(exp2factorize, factor expr, op)
        op = 
            'factor',factorize before doit
            'simplify' simplify ...
            'expand'..
            'expo', applied suerfactor in exponenets
    '''        
    return superfactor(expr,p1,*args)
    
        
        
def factorlist(expr):
    '''
    input((x+1)(x-1))
    output[x+1,x-1]
    '''
    if Is_Add(expand(expr)):
        if Is_Mul(expr):
            mm=fpoly(expr,'list')
        elif Is_Pow(expr):
            b,ee,rr=partPow(expr)
            mm=[b]*ee    
        else:
            mm=factor(expr)
            mm=fpoly(mm,'list')
        return mm
    else:
        return [expr]   
#####################################
#           conjuntos
##################################### 
def someinvec(vec1,vec2): 
    '''
    return True if some items in vec 2 are in vec1
    return false if nothing items in vec 2 are in vec1
    '''
    for i in vec2:
        if i in vec1:
            return True
    return False    
    




#####################################
#           exponents
##################################### 
''' 
def factorexp(ksym):
    if Is_Pow(ksym):
        v1,e1=base_exponent(ksym)
        e2=opemat(e1,'f')
        kres=v1**e2
        return kres
    elif Is_Mul(ksym):
        mm= fpoly(ksym,'list')
        kres=1
        for i in mm:
            kres=kres*factorexp(i)
        return kres    
    elif Is_Add(ksym):
        mm= fpoly(ksym,'list')
        kres=0
        for i in mm:
            kres=kres+factorexp(i)
        return kres    
        
    elif Is_Div(ksym):
        p1=factorexp(numer(ksym))
        p2=factorexp(denom(ksym))
        return p1/p2
            
    
    else:
         
        return ksym
'''       
 
    
    
def base2frac(ksym): # return nmonomie whit expand exponents
     
    if Is_Pow(ksym):
        kres=ksym
        mm=fpoly(ksym,'list')
        vval=mm[0]
        vexp=mm[1]
        
        if denom(vval)==1:
            vexp=-1*vexp
            kres=(1/vval)**vexp
            return kres
        else:
            return ksym
         
    elif Is_Div(ksym):
        p1=packexp(base2frac(ksym))
        p2=packexp(base2frac(ksym))
        return p1/p2
            
    elif Is_Mul(ksym):
        mm= fpoly(ksym,'list')
        kres=1
        for i in mm:
            kres=kres*base2frac(i)
        return kres    
    elif Is_Add(ksym):
        mm= fpoly(ksym,'list')
        kres=0
        for i in mm:
            kres=kres+base2frac(i)
        return kres    
    else:
        kres=ksym
        return kres
        
def packexp(ksym): # return nmonomie whit expand exponents
     
    if Is_Pow(ksym):
        kres=ksym
        mm=fpoly(ksym,'list')
        vval=mm[0]
        vexp=mm[1]
        if type(vexp)==Add:
            kres=1
            mm1=fpoly(vexp,'list')
            for i in mm1:
                kres=kres*(vval**i)
            return kres
        else:
            return ksym
         
    elif Is_Div(ksym):
        p1=packexp(numer(ksym))
        p2=packexp(denom(ksym))
        return p1/p2
            
    elif Is_Mul(ksym):
        mm= fpoly(ksym,'list')
        kres=1
        for i in mm:
            kres=kres*packexp(i)
        return kres    
    elif Is_Add(ksym):
        mm= fpoly(ksym,'list')
        kres=0
        for i in mm:
            kres=kres+packexp(i)
        return kres    
    else:
        kres=ksym
        return kres        
        
'''        
def expandexp(ksym,op=''): # return nmonomie whit expand exponents
     
    if Is_Pow(ksym):
        if op=='e':
            mm=fpoly(ksym,'list')
            p1=mm[0]
            p2=mm[1]
            kres=expandexp(p2,op='')
            kres=p1**kres
            return kres
        else:    
            kres=ksym
            mm=fpoly(ksym,'list')
            vval=mm[0]
            vexp=mm[1]
            vexp2=expand(vexp)
            kres=kpow(vval,vexp2)
            return kres
    elif Is_Div(ksym):
        p1=expandexp(numer(ksym))
        p2=expandexp(denom(ksym))
        return p1/p2
            
    elif Is_Mul(ksym):
        mm= fpoly(ksym,'list')
        kres=1
        for i in mm:
            kres=kres*expandexp(i)
        return kres    
    elif Is_Add(ksym):
        mm= fpoly(ksym,'list')
        kres=0
        for i in mm:
            kres=kres+expandexp(i)
        return kres    
    else:
        kres=ksym
        return kres
    
def simplifybase(ksym,kope=''):
         
        if type(ksym)==Pow:
            mm=fpoly(ksym,'list')
            p1=mm[0]
            p2=mm[1]
            p1=simplify(p1)
            kres=p1**p2
            if kope!='':
                kres=opemat(kres,kope=kope)
            return kres
        elif type(ksym)==Add:
            kres=0
            mm=fpoly(ksym,'list')
            for i in mm:
                kres+=simplifybase(i)
            if kope!='':
                kres=opemat(kres,kope=kope)
            return kres
        elif type(ksym)==Mul:
            kres=1
            mm=fpoly(ksym,'list')
            for i in mm:
                kres*=simplifybase(i)
            if kope!='':
                kres=opemat(kres,kope=kope)
            return kres    
        else:
            return ksym
                
 
def simplifyexp(ksym,kope=''):
    if '**' not in str(ksym):
        return ksym
    if Is_Add(ksym):
        mm=fpoly(ksym,'list')
        kres=0
        for i in mm:
            kres+=simplifyexp(i,kope=kope)
        return kres
    elif  Is_Pow(ksym):
        mm=fpoly(ksym,'list')
        val=mm[0]
        vexp=mm[1]
        if Is_Pow(val):
            mm2=fpoly(val,'list')
            val2=mm2[0]
            vexp2=mm2[1]
            nexp=vexp2*vexp
            nexp=opemat(nexp,kope=kope)
            return kpow(val2,nexp)
        else:
            nexp=simplify(vexp)
            nexp=opemat(nexp,kope=kope)
            return kpow(val,nexp)
            
        
    elif Is_Mul(ksym):
        try:
            vv=fpoly(ksym,'free')
            mm0=fpoly(ksym,'list')
            qq=len(vv)
            sexp=mzero(qq)
            kres=1
            mm=[]
            for i in mm0:
                if Is_Number(i):
                    kres=kres*i
                else:
                    mm.append(i)

            for i in mm:
                if Is_Pow(i):
                    mm2=fpoly(i,'list')
                    var=mm2[0]
                    nvexp=mm2[1]
                    knum=vv.index(var)
                    sexp[knum]=sexp[knum]+nvexp
                else:
                    knum=vv.index(i)
                    sexp[knum]=1

            for i,j in zip(vv,sexp):
                je=j
                if kope!='':
                    je=opemat(j,kope=kope)
                kres=kres*i**je
            return kres
        except:
            return ksym
    else:
        return ksym 
'''
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


        

def reducePow(ksym):
        kres= simplifyexp(mulexpo(div2mulexp(ksym)))
        return simplifyexp(kres)
#####################################
#           GDC MCM 
#####################################   

def GCD(ee1,ee2): # returm GDC between ee1,ee2 Maximo Comun Divisor
    b1=baselist(1*ee1)
    b2=baselist(1*ee2)
    e1=expolist(1*ee1)
    e2=expolist(1*ee2)
    kres=1
    for v1,ve1 in zip(b1,e1):
        for v2,ve2 in zip(b2,e2):
            if v2==v1:
                if Is_Number(ve1) and Is_Number(ve2):
                    nexp=min(ve1,ve2)
                    kres=kres*v1**nexp
    return kres
    
def LCM(v1,v2): # returm GDC between ee1,ee2 Maximo Comun Divisor
    ee1=1*v1
    ee2=1*v2
    b1=baselist(1*ee1)
    b2=baselist(1*ee2)
    e1=expolist(1*ee1)
    e2=expolist(1*ee2)
    bb=unionvec(b1,b2)
    
    kres=1
    for i in bb:
        ex1=0
        ex2=0
        ee=[]
        if i in b1 :
            p1=getindex(i,b1)
            ex1=e1[p1]
        if i in b2 :
            p2=getindex(i,b2)
            ex2=e2[p2]
        if ex1!=0:
            ee.append(ex1)
        if ex2!=0:
            ee.append(ex2)
         
        ef=min(ee)
        kres=kres*i**ef
    return kres  

    
def addexpand(expr): # converte in factor inside sum polinomies factor
    if Is_Div(expr):
        p1=addexpand(numer(expr))
        p2=addexpand(denom(expr))
         
        return cfrac(p1,p2)
    
    elif Is_Pow(expr):
        B,ee,rr=partPow(expr)
        B=addexpand(B)
         
        return kpow(B,ee,rr)

    elif Is_Mul(expr):
        mm=fpoly(expr,'list')
        kres=1
        for i in mm:
            kres=kres*addexpand(i)
            
        return kres    
    elif Is_Add(expr):
        kres=factor(expr)
         
        return kres
    else:
        return expr    
    
##########################################
#####  simplify Root !!!!!!!!!

def lsimplify(expr):
    if Is_Add(expr):
        expr2=logcombine(expr)
        if str(expr2)!=str(expr):
            return expr2
        else:
            mm=0
            for i in fpoly(expr,'list'):
                mm=mm+lsimplify(i)
            return mm    
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        P1=lsimplify(p1)
        P2=lsimplify(p2)
         
        return Div(P1,P2)
    elif Is_Mul(expr):
        mm=1
        for i in fpoly(expr,'list'):
            mm=mm*lsimplify(i)
        return mm 
    
    else:
        return expr

def lcombine(expr):
    if Is_Mul(expr):
        kres=1
        for i in expr.args:
            kres=kres*lcombine(i)
        return kres
    if Is_Div(expr):
        p1=numer(expr)
        p2=denom(expr) 
        p1=unisymbols(lcombine(p1))
        p2=unisymbols(lcombine(p2))
        return unisymbols(Div(p1,p2))
    else:
        return logcombine(unisymbols(expr),force=True)       
def simplifyroot(expr):
    if Is_Pow(expr):
        B,ee,rr=partPow(expr)
        if Is_Div(B):
            p1=numer(B)
            p2=denom(B)
            P1=kpow(p1,ee,rr)
            P2=kpow(p2,ee,rr)

            expr=cfrac(P1,P2)
            expr=killroot(expr)
        if Is_Mul(B) and not Is_Div(B):
            kres=1
            mm=fpoly(B,'list')
            for i in mm:
                kres=kres*rpow(i**ee,rr)
            expr=kres    
            expr=killroot(expr) 

              
    
    A=0
    B,ee,rr=partPow(expr)
    
    if Is_Div(ee):
        p2=denom(ee)
        rr=rr*p2
        ee=numer(ee)
    if Is_Pow(B):
        B2,ee2,rr2=partPow(B)
            
        ee=ee*ee2
        rr=rr*rr2
        B=B2
         
            
        
    if Is_Div(ee):
        p2=denom(ee)
        rr=rr*p2
        ee=numer(ee) 
    nM,nE=compdiv(ee,rr)
    P1=kpow(B,nM)
    P2=kpow(B,nE,rr)
    if P1==1:
        P3=P2
    elif P2==1:
        P3=P1
    else:
        P3=Mul(P1,P2,evaluate=False)
 
    return P3

      
def compdiv(p1,p2):
    R=p1%p2
    C=cfrac((p1-R),p2)
    return(C,R)


def simplify_root_exp(expr):
    if Is_Pow(expr):
        B,ee,rr=partPow(expr)
        if Is_Div(B):
              
            p1=numer(B)
            p2=denom(B)
            P1=rpow(p1**ee,rr)
            P2=rpow(p2**ee,rr)
            expr=cfrac(P1,P2)
    if Is_Add(expr):
        mm=0
        for i in fpoly(expr,'list'):
            mm=mm+simplify_root_exp(i)
        return mm
    elif Is_Div(expr):
        p1=simplify_root_exp(numer(expr))
        p2=simplify_root_exp(denom(expr))
        return cfrac(p1,p2)
    elif Is_Mul(expr):
        mm=1
        for i in fpoly(expr,'list'):
            mm=mm*simplify_root_exp(i)
        return mm
    else:
        return simplifyroot(expr)        
        
        

    
    
def truetable (n):
    if n < 1:
        return [[]]
    subtable = truetable(n-1)
    return [ row + [v] for row in subtable for v in [0,1] ] 


##   PLOt    

from matplotlib.patches import Polygon
def Area2(self,*args,typos=-1,tsize=10):
    var=self.var
    done=False
    ksolu=solve(self.ksym,var)
    if len(ksolu)>0:
        done=True
    if len(args)==2:
        x1=args[0]
        x2=args[1]
        x3=100
        X1=args[0]
        X2=args[1]
    elif len(args)==3:
        x1=args[0]
        x2=args[1]
        x3=args[2]
        X1=args[0]
        X2=args[1]
    else:
         
        ksolu=solve(self.ksym,var)
        x1=min(ksolu)
        X1=min(ksolu)
        x2=max(ksolu)
        X2=max(ksolu)
        x3=100
        dx=(x2-x1)/10
        x1=x1-dx
        x2=x2+dx
        done=True
    fac=(x2-x1)/x3
    X=[]
    Y=[]
    yvalor=self('value',var=fac)
     
    #for i in range(x3):
        #valor=float(x1+fac*i)
        #X.append(float(valor))
        #Yvalor=yvalor.subs(var,valor)
        
        #Y.append(float(Yvalor))
    allvar=[]
    if done:
        if len(ksolu)>0:
            allvar=ksolu
    if X1 not in allvar:
        allvar.append(X1)
    if X2 not in allvar:
        allvar.append(X2)
    ksolu=allvar
    ksolu.sort()
    X1=min(ksolu)
    X2=max(ksolu)
    fig, ax = plt.subplots()
    ix,iy =funclist(self,X1,X2,100)
     
     
    verts = [(X1,0),*zip(ix, iy),(X2,0)]
     
    poly = Polygon(verts, facecolor='0.9', edgecolor='0.5')
    ax.add_patch(poly)
    ax.grid(True, which='both')
    ax.plot(ix, iy)
    Lx=[x1,x2]
    Ly=[0,0]
    ymin,ymax=ax.get_ylim()
    ax.plot(Lx, Ly,color='black')
    if done:
         
        ysolu=funclist(self,ksolu)
         
        for i ,j in zip(ksolu,ysolu):
            ax.plot([i,i],[j,0],linestyle='-.',color='black')
            ax.plot([i],[j],'o',color='red')
        tvec=textpos(ksolu,typos,tsize)
        tex=f'${latex(Integral(self.primi,self.var))}$'
        Yp=1
        for i in tvec:
            tx,ty,tf=i
            ty=ty*Yp
            Yp=-1*Yp
            ax.text(tx, ty, tex,fontsize=tf)

            
    ax.grid(True, linestyle='-.')
    ax.tick_params(labelcolor='r', labelsize='medium', width=3)

    plt.show() 
    return ksolu
    
    
def funclist(self,*args):
    
    if len(args)>1:
        x1=args[0]
        x2=args[1]
        qq=50
        if len(args)>2:
            qq=args[2]
        xr=list(range(qq+1))    
        fac=(x2-x1)/qq
        xx=[x1+fac*i for i in xr]
    else:
        xx=args[0]

    kfunc=self.ksym
    yy=[kfunc.subs(self.var,i) for i in xx]
    if type(args[0])==list:
        return yy
    else:
        return xx,yy
        
def squarecomplete(expr, var):
    """
    Complete the square for a quadratic expression.
    
    Parameters:
    expr: sympy expression - the quadratic expression
    var: sympy symbol - the variable to complete the square for
    
    Returns:
    sympy expression - the completed square form
    """
     
    expr = simplify(expr)
    
    # Extraer coeficientes usando poly si es posible
    try:
        poly = Poly(expr, var)
        x2 = poly.coeff_monomial(var**2)
        x1 = poly.coeff_monomial(var)
        x0 = poly.coeff_monomial(1)
    except:
        # Fallback: método por coeficientes
        x2 = expr.coeff(var**2)
        x1 = expr.coeff(var)
        x0 = expr.subs(var, 0) - x2*0  # Más robusto
    
    # Si x2 es 0, no es cuadrática
    if x2 == 0:
        return expr
    
    p2 = sqrt(x2)
    p1 = x1 / (2 * p2)
    p3 = x0 - p1**2
    
    # Construir la forma
    forma = (p2 * var + p1)**2 + p3
    
    return forma 

def setpack(expr,**kwargs):
    '''
    replace value of kwaargs in exprend return expr
    '''
    doned=False
    if "'" in str(expr):
        doned=True
    if not Is_Number(expr):
        if doned:
            expr=amperdiff2normaldiff(expr)
            
        svec,vvec=unpack(kwargs)
        for i,j in zip(svec,vvec):
            if doned:
                expr=expr.replace(i,str(j))
            else:    
                expr=expr.subs(i,j)
            
        if doned:
            expr=normaldiff2amperdiff(expr)
        return expr 
    else:
        return expr

def amperdiff2normaldiff(expr):
    sQ=str(expr)
    vecvaldiff=['d2x','d2y','d2z','d2u','d2v','d2w','dx','dy','dz','du','dv','dw']
    vecvaldiffp=["x''","y''","z''","u''","v''","w''","x'","y'","z'","u'","v'","w'"]
    for i,j in zip(vecvaldiff,vecvaldiffp):
        sQ=sQ.replace(j,i)
     
    return sQ

def normaldiff2amperdiff(expr):
    sQ=parse_expr(expr)
    vecvaldiff=['d2x','d2y','d2z','d2w','d2u','d2v','dx','dy','dz','dw','du','dv']
    vecvaldiffp=[d2x,d2y,d2z,d2w,d2u,d2v,dx,dy,dz,dw,du,dv]
    for i,j in zip(vecvaldiff,vecvaldiffp):
        sQ=sQ.subs(i,j)
     
    return sQ    
    
def normalizediff(expr,sf,var):
    ''' 
    traduce sympy diff expressions to
    MyEqEq diff mat expressions
    '''
    B,C=symbols('B C')
    t=symbols(var)
    f=Function(sf)(t)
    expr=expr.subs(f.diff(var),B)
    expr=expr.subs(f,C)
    f=symbols(sf)
    df,d2f=diffsymbolprime(f)
    dt,d2t=diffsymbolprime(t)
    expr=expr.subs(C,f)
    expr=expr.subs(B,df/dt)
    x,y,z,u,v,w=symbols('x y z u v w')
    mm=normaldiff2amperdiff(amperdiff2normaldiff(expr)) 
    return mm
    
def conjugate(expr):
    '''
    if expr=a+b return a-b
    if expr=a-b return a+b
    '''
    mm=expr.args
    return (mm[1]-mm[0])   

def squareadd2squaremul(expr):
    ''''
    input (x*x-y*y)
    return (x+y)*(x-y)
    '''
    var1,var2=expr.free_symbols 
    svar1=str(var1)
    svar2=str(var2)
    k1=symbols(svar1,positive=True)
    k2=symbols(svar2,positive=True)
    expr=expr.subs(var1,k1)
    expr=expr.subs(var2,k2)
    p1,p2= expr.args
    p11=sqrt(p1)
    p111=str(p11)
    p111=p111.replace('I','1')
    p11=parse_expr(p111)
    done=False
    if signo(p2)==-1:
        done=True
        p2=-1*p2
    p22=sqrt(p2)
    p222=str(p22)
    p222=p222.replace('I','1')
    p22=parse_expr(p222)
    sp1=str(p11)
    sp2=str(p22)
    sexpr="("+sp2+"+"+sp1+")*("+sp2+"-"+sp1+")"
    kres=sympify(sexpr,evaluate=False)
    return UnevaluatedExpr(kres)


def Divunevaluate(expr):
    '''
        expr is a unevaluate expr
        expr= a*b/b
        Divunevaluate(expr) return a*
    '''
    
    expr2=expr.args[0]
    kn=numer(expr2)
    kn=sympify(str(kn))
    kd=denom(expr2)
    kd=sympify(str(kd))
    nkn=1
    nkd=1
    if type(kn)==Mul:
        nlist=kn.args
    else:
        nlist=[kn]
    if type(kd)==Mul:
        dlist=kd.args
    else:
        dlist=[kd]    
    for i in nlist:
        if i not in dlist:
            nkn=nkn*i
    for i in dlist:
        if i not in nlist:
            nkd=nkd*i
    return  nkn/nkd

 

def simplifyMul(p1,p2):
    P1=p1
    P2=p2
    if Is_Mul(p1) and Is_Mul(p2):
        P1=1
        P2=1
        mm1=p1.args
        mm2=p2.args
        for i in mm1:
            if i not in mm2:
                P1=P1*i
        for j in mm2:
            if j not in mm1:
                P2=P2*j
                
    return P1,P2
    
    
def simplifySum(p1,p2):
    P1=p1
    P2=p2
    if Is_Add(p1) and Is_Add(p2):
        P1=0
        P2=0
        mm1=p1.args
        mm2=p2.args
        for i in mm1:
            if i not in mm2:
                P1=P1+i
        for j in mm2:
            if j not in mm1:
                P2=P2+j
                
    return P1,P2
    
def precross(expr):
    if Is_Add(expr):
        expr=factor(expr)
        p1,p2=fraction(expr)
        
    elif Is_Div(expr):
        p1,p2=fraction(expr)
    elif Is_Pow(expr):
        
        bb=getbase(expr)
        ee=getexpo(expr)
        if Is_Div(bb):
            p1,p2=fraction(bb)
            p1=p1**ee
            p2=p2**ee
    else:
        p1=expr
        p2=1
    return p1,p2

def monodegree(obj,var=x): 
    expr= obj 
    if not str(var) in str(expr):
        return 0
    else:    
        p=expr.subs(var,1)
        p2=simplify(expr/p)
        p3=getexpo(p2)
        return p3
def monocoef(obj,var=x): 
    expr= obj 
    p=expr.subs(var,1)
    return p
    
from operator import itemgetter
def infodegreemat(obj,var=x):
    kres=[]
    expr= obj 
    mm=expr.args
    for i in mm:
        coe=monocoef(i,var=var)
        gra=monodegree(i,var=var)
        kres.append([coe,gra,i])
    kres=sorted(kres, key=itemgetter(1), reverse=True)
    vecc=[]
    vecg=[]
    vece=[]
    for i in range(len(kres)):
        vec=kres[i]
        vecc.append(vec[0])
        vecg.append(vec[1])
        vece.append(vec[2])
    return vecc,vecg,vece

def trinomio2squarebinomio(obj,var=x):
    expr=obj
    mm=expr.args
    if len(mm)==3:
        vecc,vecg,vece=infodegreemat(expr,var=var)
        #if (vecg[0]+vecg[2])/2==vecg[1]:
        ksigno=signo(vece[1])
        a=rsimplify(sqrt(vece[0]))             
        b=vece[1]/2/a
        
        ebin=(a+b)**2
        eebin=expand(ebin)
        resto=simplify(eebin-vece[0]-vece[1]-vece[2])
        return  ebin+resto
             
    else:
        return expr
def cos2halfcos(expr,alpha):
    return eval(str((expr.subs(cos(alpha),2*scos(alpha/2)**2-1,evaluate=False))))
def sin2halfsin(expr,alpha):
    return eval(str((expr.subs(sin(alpha),2*ssin(alpha/2)*scos(alpha/2)))))
def alpha2halfalpha(expr,alpha):
    return sin2halfsin(cos2halfcos(expr,alpha),alpha/2)
def ang2halfang(expr,alpha):
    return sin2halfsin(cos2halfcos(expr,alpha),alpha)    
def angle2halfangle(expr,alpha):
    return sin2halfsin(cos2halfcos(expr,alpha),alpha)    
def ang2halfang(expr,alpha):
    return sin2halfsin(cos2halfcos(expr,alpha),alpha)
def tri2squ(obj,var=x):
    expr=obj
    mm=expr.args
    if len(mm)==3:
        vecc,vecg,vece=infodegreemat(expr,var=var)
        #if (vecg[0]+vecg[2])/2==vecg[1]:
        ksigno=signo(vece[1])
        a=sqrt(vece[0])             
        b=vece[1]/2/a
        
        ebin=(a+b)**2
        eebin=expand(ebin)
        resto=simplify(eebin-vece[0]-vece[1]-vece[2])
        return  ebin+resto
             
    else:
        return expr
def squaresin2cos(expr,alpha):
    b = Wild("b")
    expr.match(expr.subs(alpha,b))
    kres= str((expr.subs(sin(alpha)**2,1-cos(b)**2)).xreplace(expr.match(expr.subs(alpha,b))))
    return parse_expr(kres,evaluate=False)
def tri2bin(expr,sexpr,var=''):
    if var=='':
        var=self.var
    nexpr=parse_expr(sexpr)
    expr2=tri2squ(nexpr,var=var)
    msexpr=str(expr)
    kres=subsubs(expr,nexpr,expr2)
    return kres

def functiondiff(*args):
    '''
    if Area= b*h/2 in function of time...
    then A(t)=b(t)*h(t)/2...
    input :  functiondiff(expr, var, varF1, varF2,varF3...etc)
    in this case:
     A=functiondiff(A,t,A) = dA
     R=functiondiff(b*h/2,t,b,h) = b*dh/2 + db*h/2 
 
    '''
    expr=args[0]
    var=args[1]
    if len(args)==2:
        return diff(expr,var)
    else:    
        v1=[]
        for i in range(2,len(args)):
            v1.append(args[i])
        v2=[]
        cc=1
        for i in v1:
            F='f'+str(i)
            F=Function(str(i))(var)
            v2.append(F)
            cc=cc+1
        v3=[]
        for i in v2:
            v3.append(diff(i,var))
        v4=[]
        for i in v3:
            v4.append(str(i))
        v5=[]
        for i in v1:
            svar='d'+str(i)
            svar=symbols(svar)
            v5.append(svar)
            
        
        for i,j in zip(v1,v2):
            expr=expr.subs(i,j)
        dexpr=diff(expr,var)
        for i,j in zip(v3,v5):
            dexpr=dexpr.subs(i,j)
        for i,j in zip(v2,v1):
            dexpr=dexpr.subs(i,j)    
        return dexpr    
        
        
def get_fact_inside(obj):
    ss=str(obj)
    ss=ss.replace('factorial(','')
    ss=ss.replace(')','')
    return parse_expr(ss)  
   
def expandfact(obj,val=1):
    expr=get_fact_inside(obj)
    kres=Product(expr-x, (x, 0, val-1)).doit()
    kres=kres*factorial(expr-val)
    return kres

def imagip(expr,op=''):
    cdone=False
    ddenomm=1
    if denom(expr)!=1 and not 'I' in str(expr):
        cdone=True
        ddenom=denom(expr)
        expr=numer(expr)
        
    vecr= 0
    veci= 0
    expr=expand(expr)
    for data in expr.args:
        if 'I' in str(data):
            sdata=str(data)
            sdata=sdata.replace('I','1')
            veci+=parse_expr(sdata)            
        else:
            vecr+=data
    vecr=factor(simplify(vecr))
    if cdone:
        vecr=sdiv(vecr,ddenom)
    veci=factor(simplify(veci))
    if cdone:
        veci=sdiv(veci,ddenom)    
    sre='('+str(vecr)+')'
    sim='('+str(veci)+')*I'
    
    if op=='re':
        return vecr
    elif op=='im':
        return veci
    else:
        return parse_expr(sre+'+'+sim,evaluate=False)
        
        
def frac2onesum(expr,*args):
    '''
    if expr is a+b,frac2onesum(expr) return a*(1+b/a)
    if expr is a+b,frac2onesum(expr,'swap') return b*(1+a/b)   
        
    if expr is (a+b)**c,frac2onesum(expr) return a**c*(1+b/a)**c
    if expr is (a+b)**c,frac2onesum(expr,'swap') return b**c*(1+a/b)**c 
    '''
    if Is_Add(expr):
        if len(expr.args)==2:
            mm=expr.args
            p1=mm[0]
            p2=mm[1]
            if 'swap' in args:
                p1=mm[1]
                p2=mm[0]
            if p2==1:
                ladiv=sinversa(p1)
            else:    
                ladiv=sdiv(p2,p1)    
            return p1*(1+ ladiv) 
        else:
            return expr
    elif Is_Div(expr):
        data1=numer(expr)
        data2=denom(expr)
        return parse_expr('('+str(frac2onesum(data1,*args))+')/('+str(frac2onesum(data2,*args))+')',evaluate=False)
    elif Is_Mul(expr):
        kres=1
        for data in expr.args:
            kres=smul(kres,frac2onesum(simplify(data),*args))
        return kres
    elif Is_Pow(expr):
        bb=getbase(expr)
        ee=getexpo(expr)
        if Is_Add(bb):
            if len(bb.args)==2:
                mm=bb.args
                p1=mm[0]
                p2=mm[1]
                if 'swap' in args:
                    p1=mm[1]
                    p2=mm[0]
                if p2==1:
                    ladiv=sinversa(p1)
                else:    
                    ladiv=sdiv(p2,p1)    
                return smul(spow(p1,ee),spow(1+ladiv,ee))
            else:
                return spow(frac2onesum(bb,*args),ee)
        else:
            return spow(frac2onesum(bb,*args),ee)
            
    else:
        return expr
def mysortlist(L, index, *args):
    """
    Ordena una lista de listas según el elemento en la posición `index`.
    
    Parámetros:
    L (list of list): Lista de listas a ordenar.
    index (int): Índice del elemento por el cual ordenar.
    reverse (bool): Si es True, ordena en orden descendente.
    
    Retorna:
    list of list: Lista ordenada según el índice dado.
    """
    if 'reverse' in arga:
        return sorted(L, key=lambda x: x[index], reverse=reverse)
    else:
        return sorted(L, key=lambda x: x[index])        



def primefactorcomplemento(p1,p2,*args):
    '''
    primefactorcomplemento(p1,p2) = p3 then primefactor(p3)*primefactor(p1) is part of primefactor(p2)
    args:
        'base',  return only bases as list
        'expo',  return exponents
        'data',  return bases and expo
        'value',  return value bases ** expo
    '''
    bb1,ee1=primefactor(p1,'data')
    bb2,ee2=primefactor(p2,'data')
    kres=1
    bb3,ee3=[],[]
    for kb,ke in zip(bb1,ee1):
        if not kb in bb2:
            bb3.append(kb)
            ee3.append(ke)
            kres=kres*kb**ke
        else:
            ni=bb2.index(kb)
            nbb=kb
            nee=ee2[ni]
            if ke>nee:
                bb3.append(kb)
                ee3.append(ke-nee)
                kres=kres*kb**(ke-nee)
            
    if len(args)>0:        
        if 'base' in args:
            return bb3
        elif 'expo' in args:
            return ee3
        elif 'data' in args:
            return bb3,ee3
        elif 'value' in args:
            return float2int(kres)
        else:    
            return primefactor(kres)
    return  primefactor(kres) 

def primefactorlist(expr):
    k=primefactor(expr)
    mm=k.args 
    vec=[]
    for data in mm:
        bb=getbase(data)
        ee=getexpo(data)
        vec.append([bb,ee])
    return vec

def primefactornumber(numb,*args):
    '''
    return primefactors data in numb
    op base: return base
    op expo return expo
    '''
    kres= factorint(numb) 
    if 'base' in args:
        return primefactors(numb)
    if 'expo' in args:
        k=list(kres.keys())
        return [kres[i] for i in k]
        
    return kres     
from sympy import (
    Integer, Symbol, Mul, Pow, Add, log, exp,
    factorint, numer, denom, sympify
)

# --- Prime Factorization ---
def primefactor(expr):
    """
    Factoriza números primos dentro de una expresión simbólica,
    respetando la estructura del árbol.
    """
    # Caso número entero
  
    
    if isinstance(expr, (int, Integer)):
        if expr == 1:
            return Integer(1)
        fac = factorint(expr)  # dict {p: exp}
        claves = list(fac.keys())  # [2, 3, 5]
        valores = list(fac.values())  # [3, 2, 1]
        if valores[0]==1:
            kres=claves[0]
        else:
            kres=spow(claves[0],valores[0])
        if len(claves)>1:    
            for p,e in zip(claves[1::],valores[1::]):
                if e==1:
                    kres=smul(kres,p)
                else:
                    kres=smul(kres,spow(p,e))
        return kres
    # Caso suma
    if expr.is_Add:
        return Add(*[primefactor(a) for a in expr.args], evaluate=False)

    # Caso multiplicación
    if expr.is_Mul:
        return Mul(*[primefactor(a) for a in expr.args], evaluate=False)

    # Caso división (Mul con potencia negativa, pero explícito)
    if expr.is_Rational:
        return Mul(
            primefactor(numer(expr)),
            Pow(primefactor(denom(expr)), -1, evaluate=False),
            evaluate=False
        )

    # Caso potencia
    if expr.is_Pow:
        base, expn = expr.args
        return Pow(primefactor(base), expn, evaluate=False)

    # Caso logaritmo
    if expr.func == log:
        (arg,) = expr.args
        return log(primefactor(arg), evaluate=False)

    # Caso exponencial
    if expr.func == exp:
        (arg,) = expr.args
        return exp(primefactor(arg))

    # Caso símbolo u otra cosa: se deja igual
    return expr


# --- Helpers ---
def primefactorbase(numb,full=False):
    if full:
        bb=primefactorbase(numb)
        ee=primefactorexpo(numb)
        kres=[]
        for e1,b1 in zip(ee,bb):
            kres=kres+e1*[b1]
        return kres    

    else:    
        d=primefactor(numb) 
        if Is_Pow(d):
            return [getbase(d)]
        else:
            vec=[]
            for data in d.args:
                vec.append(getbase(data))
            return vec    
def primefactorexpo(numb):
    """
    Retorna lista de bases en la factorización.
    """
    d=primefactor(numb) 
    if Is_Pow(d):
        return [getexpo(d)]
    else:
        vec=[]
        for data in d.args:
            vec.append(getexpo(data))
        return vec 

def primefactorcountbase(numb):
    '''
    return count list of bases in factorint(numb)
    '''
    return len(list(factorint(numb).keys()))
def primefactorcountexpo(numb):
    '''
    return count list of exponentes  in factorint(numb)
    '''
    return len(list(factorint(numb).values()))
    
def primefactormod(expr):
    """
    Factoriza números primos dentro de una expresión simbólica,
    respetando la estructura del árbol.
    """
    # Caso número entero
    if isinstance(expr, (int, Integer)):
        if expr == 1:
            return Integer(1)
        fac = factorint(expr)  # dict {p: exp}
        factors = []
        for p, e in fac.items():
            if e == 1:
                factors.append(p)
            else:
                factors.append(Pow(p, e, evaluate=False))
        return Mul(*factors, evaluate=False)

    # Caso suma
    if expr.is_Add:
        return Add(*[primefactormod(a) for a in expr.args], evaluate=False)

    # Caso multiplicación
    if expr.is_Mul:
        return Mul(*[primefactormod(a) for a in expr.args], evaluate=False)

    # Caso división (Rational)
    if expr.is_Rational:
        return Mul(
            primefactormod(numer(expr)),
            Pow(primefactormod(denom(expr)), -1, evaluate=False),
            evaluate=False
        )

    # Caso potencia
    if expr.is_Pow:
        base, expn = expr.args
        return Pow(primefactormod(base), expn, evaluate=False)

    # Caso logaritmo
    if expr.func == log:
        (arg,) = expr.args
        return log(primefactormod(arg), evaluate=False)

    # Caso exponencial
    if expr.func == exp:
        (arg,) = expr.args
        return exp(primefactormod(arg))

    # Caso símbolo u otra cosa
    return expr
    
from math import gcd

def IsCoprime(x, y):
    ''' 
    return True if are comprime x and y 
    '''
    return gcd(x, y) == 1     
def testfunction(sfunc,examples):
    kres=[]
    for data in examples:
        kres.append(eval(sfunc+'('+str(data)+')'))
    for exam,kres in zip(examples,kres):
        display(Math(sfunc+'('+latex(exam)+')='+latex(kres)))

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
        
     
def secdiv(expr,factor):
    if Is_Add(expr):
        kres=0
        for data in expr.args:
            kres=kres+secdiv(data,factor)
        return kres
    elif Is_Div(expr):
        nnumer=numer(expr)
        dd=denom(expr)
        return Div(secdiv(nnumer,factor),secdiv(dd,factor))
    elif Is_Mul(expr):
        kres=1
        for data in expr.args:
            kres=kres*secdiv(data,factor)
        return kres
    elif Is_Root(expr):
        rr=getroot(expr)
        bb=insideroot(expr)
        return rpow(secdiv(bb,factor**rr),rr)
    else:
        return simplify(sdiv(expr,factor)) 
        
def getfullexpo(expr):
    if isinstance(expr,Pow):
        bb=expr.args[0]
        ee=expr.args[1]
        if isinstance(bb,Pow):
            bb2=bb.args[0]
            ee2=bb.args[1]*ee
            return bb2,numer(ee2),denom(ee2)
        return bb,numer(ee),denom(ee)
    else:
        bb=expr
        ee=1
        rr=1
        return bb,ee,rr  
        
def rsimplify(expr):
    if Is_Number(expr):
        return expr
    elif Is_Symbols(expr):
        return expr
    elif Is_Add(expr):
        vec=[ rsimplify(data) for data in expr.args]
        return sum(vec)
    elif Is_Div(expr):
        return cf( rsimplify(numer(expr)), rsimplify(denom(expr)))
    elif Is_Mul(expr):
        return prod([ rsimplify(data) for data in expr.args])


        
    elif Is_Root(expr):
        if Is_Mul(insideroot(expr)) and not Is_Div(insideroot(expr)):
            bb,ee,rr=getfullexpo(expr)
            bb=joinbase(bb)
            if Is_Pow(bb):
                bb2,ee2,rr2=getfullexpo(bb)
                bb=bb2
                ee=ee*ee2
                rr=rr*rr2
            if rr==1:
                return bb**ee
            else:
                return bb**cf(ee,rr)    
            
        if Is_Pow(insideroot(expr)):
            bb,ee,rr=getfullexpo(expr)
            if rr==1:
                return bb**ee
            else:
                return bb**cf(ee,rr)
  
        elif Is_Div(insideroot(expr)):
            bb,ee,rr=getfullexpo(expr)
            return cf( rsimplify(numer(bb)**cf(ee,rr)), rsimplify(denom(bb)**cf(ee,rr)))
        else:
            bb,ee,rr=getfullexpo(expr)
            return bb**cf(ee,rr)
    else:
        return expr 

def uproot(expr,rr):
    ee=getexpo(expr)
    if Is_Number(ee) and Is_Number(rr):
        if ee%rr>=0:
            ee1=ee%rr
            ee0=int((ee-ee1)/rr)
            return ee0,ee1
        else:
            return 0,getexpo(expr) 
             
    else:
        return 0,getexpo(expr) 
def rmulroot(expr):
    BB=getonlybase(expr)
    EE=getpow(expr)
    RR=getroot(expr)
    if Is_Mul(BB):
        fac1=1
        fac2=1
        for data in BB.args:
            bb=getbase(data)
            ee=getexpo(data)
            e0,e1=uproot(data,RR)
            if e0>0:
                fac1=fac1*bb**e0
            if e1>0:
                fac2=fac2*bb**e1
        if fac2==1:
            return fac1
        elif fac1==1:
            return Pow(fac2,cf(1,RR))
        else:    
            return Mul(fac1,Pow(fac2,cf(1,RR)),evaluate=False)
    else:
        return expr
def rrsimplify(expr):
    if Is_Number(expr):
        return expr
    elif Is_Symbols(expr):
        return expr
    elif Is_Add(expr):
        vec=[rrsimplify(data) for data in expr.args]
        return sum(vec)
    elif Is_Div(expr):
        return sdiv(rrsimplify(numer(expr)),rrsimplify(denom(expr)))      
    elif Is_Mul(expr):
        vec=[rrsimplify(data) for data in expr.args]
        return prod(vec)
    elif Is_Root(expr):
        base=getonlybase(expr)
        if Is_Mul(base):
            return rmulroot(expr)
        else:    
            ee=getpow(expr)
            rr=getroot(expr)
            fac1=''
            fac2=''
            if Is_Mul(base):
                vec0=[]
                vec1=[]
                for data in base.args:
                    bb=getbase(data)
                    efuera,edentro=uproot(bb,ee,rr)
                    if efuera>0:
                        if fac1=='':
                            fac1=bb**efuera
                        else:
                            fac1=fac1*bb**efuera
                    if edentro>0:
                        if fac2=='':
                            fac2=bb**edentro
                        else:
                            fac1=fac1*bb**efuera                        
                    
            bb=getbase(insideroot(expr))
            if Is_Mul(bb):
                lvar=list(bb.free_symbols)
                for nvar in lvar:
                    svar=str(nvar)
                    svar=symbols(str(nvar),positive=True)
                expr2=simplify(parse_expr(str(expr)))
                return unisymbols(expr2)
            else:       
                bb=rrsimplify(getbase(insideroot(expr)))
                ee=getpow(expr)
                rr=getroot(expr)
                newexpo=simplify(cf(ee,rr))
                ee=numer(newexpo)
                rr=denom(newexpo)
                if ee==1 and rr==1:
                    return bb
                elif rr==1:
                    return Pow(bb,ee,evaluate=False)
                else:
                    return Pow(bb,cf(1,rr),evaluate=False)
                
    elif Is_Pow(expr):
        bb=rrsimplify(getbase(expr))
        ee=getpow(expr)
        return bb**ee
    else:
        return expr

    
def rootsimplify(expr): # simplify each ecponet in expr
    if Is_Add(expr):
        kres = sum([rootsimplify(data) for data in expr.args])    
    elif Is_Div(expr):
        p1,p2=fraction(expr)
        P1=rootsimplify(p1)
        P2=rootsimplify(p2)
        kres = sdiv(P1,P2)
    elif Is_Mul(expr):
        kres = prod([rootsimplify(data) for data in expr.args])
         
    elif Is_Root(expr):
         
        rr=getroot(expr)
        base=insideroot(expr)
        bb = rootsimplify(getbase(base))
        ee=getexpo(base)
        if denom(simplify(cfrac(ee,rr)))==1:
            kres = bb**simplify(cfrac(ee,rr))
     
        elif Is_Div(base):
            p1=rpow(numer(base)**ee,rr)
            kres1=rootsimplify(p1)
             
            p2=rpow(denom(base)**ee,rr)
            kres2=rootsimplify(p2)

            if not Is_Root(kres1) or not Is_Root(kres2):
                kres =  sdiv(kres1,kres2) 
            else:
                kres =  expr
        else:    
            try:
                p1=ee%rr
                p2=int(ee/rr)
                if p1==0:
                    kres=bb**p2 
                else:
                     
                    kres=Mul(bb**p2,rpow(bb**p1,rr))
                kres = kres
            except:
                kres = expr
 
        
    elif Is_Pow(expr):
        bb = getbase(expr)
        ee = getexpo(expr)
        kres = spow(rootsimplify(bb),ee)
    else:
        kres = expr 

    if Is_Mul(kres) and kres.args[0]==1:
        kres=prod(kres.args[1::])
    return kres       
from sympy import Integral, Symbol, latex        
class CustomIntegral(Integral):
    def _repr_latex_(self):
        # Integral indefinida ∫ 1 dx
        if self.function == 1 and len(self.limits) == 1:
            # Compruebo si límite es tipo (var,) o (var,a,b)
            limit = self.limits[0]
            if len(limit) == 1:
                var = limit[0]
                return r"$\int\, d%s$" % latex(var)
            elif len(limit) == 3:
                var, a, b = limit
                return r"$\int_{%s}^{%s} d%s$" % (latex(a), latex(b), latex(var))
        # Por defecto
        return super()._repr_latex_()

    def _latex(self, printer=None):
        if self.function == 1 and len(self.limits) == 1:
            limit = self.limits[0]
            if len(limit) == 1:
                var = limit[0]
                return r"\int d%s" % latex(var)
            elif len(limit) == 3:
                var, a, b = limit
                return r"\int_{%s}^{%s} d%s" % (latex(a), latex(b), latex(var))
        return super()._latex(printer)     
        
from IPython.display import display, Math 

def creteIntegral(expr,var,x1='',x2=''):
    return createintegral(expr,var,x1=x1,x2=x) 
def createIntegral(expr,var,x1='',x2=''):
    return createintegral(expr,var,x1=x1,x2=x)    
def createintegral(expr,var,x1='',x2=''):
    if 'd'+str(var) in str(expr):
        expr= expr.subs('d'+str(var),1)
    if str(expr)=='1':
        if x1=='':
            return CustomIntegral(expr,var)
        else:
            return CustomIntegral(expr,(var,x1,x2))
    else:
        if x1=='':
            return Integral(expr,var)
        else:
            return Integral(expr,(var,x1,x2))
def integra(expr,var,vecv='',x1='',x2=''):
    vecdt=['d_x','d_y','d_z','d_v','d_w','d_u']
    sexpr=str(expr)
    for data in vecdt:
        if data in sexpr:
            sexpr=sexpr.replace(data,'1')
            
        expr=parse_expr(sexpr)  
     
    if vecv!='':
        expr=renewvar(expr,vecv)
    if x1 != "":
        kres=createintegral(expr, var, x1=x1, x2=x2)
    else:
        kres=createintegral(expr, var)
    return kres            
def creatediff(expr,var):
    if expr==1:
        final_latex='d'+str(var)
        return Math(final_latex )
    else:    
        lexpr = latex(expr)
        final_latex = rf"{lexpr} \, "+"d"+str(x)
    
        return Math(final_latex )
        
def calculateFraction(a, b):
  
    # If the numerator is zero, answer is "0"
    if a == 0:
        return "0"

    # If exactly one of the numerator or denominator
    # is negative, then result will be negative
    res = "-" if (a < 0) ^ (b < 0) else ""

    a = abs(a)
    b = abs(b)

    # Calculate and Append the part before decimal point
    res += str(a // b)

    rem = a % b

    # If completely divisible, return res
    if rem == 0:
        return res

    res += "."
    mp = {}

    while rem > 0:
      
        # If this remainder is already seen,
        # then there exists a repeating fraction.
        if rem in mp:
            res = res[:mp[rem]] + "(" + res[mp[rem]:] + ")"
            break
        
        # If the remainder is seen for the first time,
        # store its index
        mp[rem] = len(res)

        rem = rem * 10

        # Calculate quotient, append it to result and
        # calculate next remainder
        res += str(rem // b)
        rem = rem % b

    return res

def decimalparts(a,b):
    '''
    decimal=1/7 = 0.14285714285714285 , decimalparts(1,7) = ('0', '0', '142857')
    decimal=1/3 = 0.33333.. , decimalparts(1,3) = ('0', '0', '3')
    decimal=1/6= 0.1666666.. , decimalparts(1,6) = ('0', '1', '6')
    '''
    
    sexpr=calculateFraction(a, b)
    if '(' in sexpr:
        p1=sexpr.find('(')
        p2=sexpr.find(')')
        pperiodic=sexpr[p1+1:p2]
    else:
        pperiodic='0'
    p1=sexpr.find('.')
    if '.(' in sexpr:
        pnoperiodic='0'
        
    elif not '(' in sexpr:
        pnoperiodic=sexpr[p1+1::]
    else:    
        p2=sexpr.find('(')
        pnoperiodic=sexpr[p1+1:p2]
 
    pentera= sexpr[0:p1]
    return (pentera,pnoperiodic,pperiodic)
     
     
def findconsecutiveprimelist(L,*args):
    '''
    retorna una lista con todas las listas de numeros primos consecutivo
    findconsecutiveprimelist([p1,p2,p3,p7,p9,p10,p11,p12,p13,p14,p15,p20])=  ( [p1,p2,p3] , [p7] , [p9,p10,p11,p12,p13,p14,p15] , [p20] )
    findconsecutiveprimelist([p1,p2,p3,p7,p9,p10,p11,p12,p13,p14,p15,p20],'long')=  [p9,p10,p11,p12,p13,p14,p15]
    '''
    cc=primepi(L[0])
    vec=[]
    p1=L[0]
    vec2=[p1]
    for data in L[1::]:
        if nextprime(vec2[-1])==data:
            vec2.append(data)
        else:
            vec.append(vec2)
            vec2=[data]
    vec.append(vec2)
    if 'long' in args:
        vecm=[len(data) for data in vec]
        kres=vec[vecm.index(max(vecm))]
        return kres
    else:
        return vec     
        
def Is_Consecutiveprimelist(L):
    '''
    return True if L is a consecutive prime list
    '''
    qq=len(L) 
    pnumb=L[0]
    done=True
    for i in range(1,qq):
        if nextprime(pnumb)!=L[i]:
            return False
        else:
            pnumb=L[i]
    return done         
    
def binomialfact(x,y):
    return cfrac(fac(x),fac(y)*fac(x-y))

def reducefact(expr,nn):
    if Is_Add(expr):
        kres=0
        for data in expr.args:
            kres=kres+reducefact(data,nn)
        return kres
    elif Is_Div(expr):
        return cf(reducefact(numer(expr),nn),reducefact(denom(expr),nn))
    elif Is_Mul(expr):
            mm=expr.args 
            kres= reducefact(mm[0],nn)
            for data in mm[1::]:
                kres=kres*reducefact(data,nn)
            return kres
    elif type(expr)==factorial:
        (data,)=expr.args
        obj1=0
        for obj in data.args:
            if Is_Number(obj):
                obj1=obj
                break
                
        if obj1>nn:
            lim=obj1-nn
            return expandfact(expr,lim+1)
        else:
            return expr
    else:
        return expr    
    
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


def _doit(expr):
    if Is_Add(expr):
        mm=expr.args
        kres=mm[0].doit()
        for data in mm[1::]:
            kres=kres+data.doit()
        return kres
        
            
    
def evalif(expr,**kwargs):
    '''
    evalif( x+y,x=1,y=2) return 3
    '''
    expr=real_subs(expr,**kwargs)
    return expr   
    
def absexpand(expr):
    # Si ya no hay absolutos, devolver la expresión
    if not expr.has(Abs):
        return (expr,)
    
    # Tomar un absoluto cualquiera (el primero encontrado)
    abs_sub = next(iter(expr.atoms(Abs)))
    arg = abs_sub.args[0]
    
    # Dos casos: reemplazar Abs(arg) por +arg o -arg
    expr_pos = expr.xreplace({abs_sub: arg})
    expr_neg = expr.xreplace({abs_sub: -arg})
    
    # Recursión: seguir descomponiendo dentro de cada caso
    return descompone_asb(expr_pos) + descompone_asb(expr_neg)    
    
def get_facpart(expr):
    '''
    Analiza un factor y devuelve (coeficiente, exponente)
    '''
    if Is_Number(expr):
        ff = expr
        nn = 0
    elif Is_Mul(expr):
        ff = expr.subs(x, 1)
        opart = simplify(expr / ff)
        if Is_Div(opart) and 'x' in str(denom(opart)):
            bb = denom(opart)
            nn = -getexpo(bb)
        else:
            nn = getexpo(opart)
    elif Is_Pow(expr):
        ff = 1
        nn = getexpo(expr)
    else:
        ff = 1
        nn = 0
    return ff, nn

def binomialcoef(expr,target_exp,*args, factor=''):
    if 'sumterm' in args:
        return primefactor(expr.subs(x,1))
    else:    
        return float2int(termbinomial(expr, target_exp, factor=factor))
    
def termbinomial(expr, target_exp, factor=''):
    from sympy import symbols, binomial, Mul, Add, Pow, factorial, Symbol, Integer
    
    x = symbols('x')
    
    # CASO 1: Si hay factor externo
    if factor != '':
        if isinstance(factor, Add):
            kres = 0
            for data in factor.args:
                kres = kres + termbinomial(expr, target_exp, factor=data)
            return 
        else:
            ff, nn = get_facpart(factor)
            return ff * termbinomial(expr, target_exp - nn)
    
    # CASO 2: Si expr es una multiplicación
    elif isinstance(expr, Mul):
        factors = expr.args
        binomial_part = None
        other_factors = []
        
        for f in factors:
            if isinstance(f, Pow) and isinstance(f.base, Add):
                binomial_part = f
            else:
                other_factors.append(f)
        
        if binomial_part is None:
            raise ValueError("No se encontró parte binomial en la expresión")
        
        if other_factors:
            other_expr = Mul(*other_factors)
            return termbinomial(binomial_part, target_exp, factor=other_expr)
        else:
            return termbinomial(binomial_part, target_exp)
    
    # CASO 3: expr es una potencia binomial (a*x^m + b)^n
    elif isinstance(expr, Pow):
        n = getexpo(expr)  # exponente de expr
        bb = getbase(expr)  # base de expr
        
        # Si n o target_exp son simbólicos (no numéricos enteros), usar fórmula general
        if not isinstance(n, Integer) or not isinstance(target_exp, Integer):
            b = bb.subs(x, 0)   # término independiente
            expr2 = bb - b      # término dependiente de x
            
            # Si expr2 = 0 (base constante)
            if expr2 == 0:
                # Solo hay término si target_exp == 0
                from sympy import KroneckerDelta
                return bb**n * KroneckerDelta(target_exp, 0)
            
            a = expr2.subs(x, 1)  # coeficiente de x^m
            
            # Si a = 0 (base constante)
            if a == 0:
                return bb**n * KroneckerDelta(target_exp, 0)
            
            expr3 = expr2 / a   # expresión de x en el término lineal
            
            if expr3 == x:
                m = 1
            else:
                m = getexpo(expr3)
            
            # Si m = 0 (base constante)
            if m == 0:
                return bb**n * KroneckerDelta(target_exp, 0)
            
            # Fórmula simbólica general
            # k = n - (target_exp / m)
            k = n - target_exp / m
            
            from sympy import factorial
            coef = factorial(n) / (factorial(k) * factorial(n - k))
            return coef * (a**(n - k)) * (b**k)
        
        else:
            # n y target_exp numéricos (código original)
            b = bb.subs(x, 0)
            expr2 = bb - b
            
            if expr2 == 0:
                if target_exp == 0:
                    return bb**n
                else:
                    return 0
            
            a = expr2.subs(x, 1)
            
            if a == 0:
                if target_exp == 0:
                    return bb**n
                else:
                    return 0
            
            expr3 = expr2 / a
            
            if expr3 == x:
                m = 1
            else:
                m = getexpo(expr3)
            
            if m == 0:
                if target_exp == 0:
                    return bb**n
                else:
                    return 0
            
            if target_exp % m != 0:
                return 0
            
            k = n - (target_exp // m)
            
            if k < 0 or k > n:
                return 0
            
            coeficiente = binomial(n, k) * (a**(n - k)) * (b**k)
            return coeficiente
    
    else:
        if target_exp == 0:
            return expr
        else:
            return 0 
            
def sbinomial(n,m):
    if m==0:
        return 1
    a,b=symbols('a b')
    '''retorna el binomial sin evaluara '''
    k2=Mul(factorial(m,evaluate=False), factorial(n-m,evaluate=False)) 
    k1=factorial(n,evaluate=False)
    mm=cfrac(a,b)
    return mm.subs(a,k1).subs(b,k2)

def sfac(expr):
    return factorial(expr,evaluate=False)   

#####  imaginary Numbers    
def rpart(expr):
    expr=expand(expr)
    vec=0
    if Is_Add(expr):
        for data in expr.args:
            if not 'I' in str(data):
                vec+=data
        return vec
    else:
        if 'I' in str(data):
            return 0
        else:
            return data    
def ipart(expr):
    expr=expand(expr)
    vec=0
    if Is_Add(expr):
        for data in expr.args:
            if 'I' in str(data):
                sdata=str(data)
                sdata=sdata.replace('I','1')
                data = parse_expr(sdata)
                vec+=data
                
        return vec
    else:
        if 'I' in str(data):
            sdata=str(data)
            sdata=sdata.replace('I','1')
            return parse_expr(sdata)
        else:
            return 0 
def imodule(expr):
    rp=rpart(expr)
    ip=ipart(expr)
    kres=sqrt(rp*rp+ip*ip)
    return simplify(kres) 

def simplifyimg(expr):
    if Is_Div(expr) and Is_Img(denom(expr)):
        knum=numer(expr)
        kden=denom(expr)
        p2=ipart(kden)
        p1=rpart(kden)
        newd=p1-p2*I
        knum=knum*(newd)
        kden=simplify(expand(kden*newd))
        return cf(knum,kden)
    else:
        return simplify(expr)    