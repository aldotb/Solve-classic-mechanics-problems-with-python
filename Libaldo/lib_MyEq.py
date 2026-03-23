from sympy import *

from lib_Variables import *
from lib_Mathbasic import *
from lib_Mathematica import *
from lib_Algorith import *
from lib_Exponencial import *
 
import copy
import numpy as np


from sympy import *
from IPython.display import display, Math
 
def get_diff_var(var):
    # find therespetive differential varible for x,y,z....
    return dLvar[Lvar.index(var)] 
    
from sympy import Add, Mul, Pow

def _map_expr(expr, f, *types):

    # --- atomic nodes ---
    if expr.is_Atom:
        return expr

    # --- default types ---
    if not types:
        types = ('Add', 'Pow')

    # --- atomic transform attempt ---
    new = f(expr)
    if new != expr:
        return new

    for t in types:

        # --- Add ---
        if t == 'Add' and expr.is_Add:
            return Add(*(_map_expr(a, f, *types) for a in expr.args))

        # --- Mul ---
        if t == 'Mul' and expr.is_Mul and not Is_Div(expr):
            return Mul(*(_map_expr(a, f, *types) for a in expr.args))

        # --- Div ---
        if t == 'Div' and Is_Div(expr):
            num, den = expr.args
            num2 = _map_expr(num, f, *types)
            den2 = _map_expr(den, f, *types)
            return cf(num2, den2)

        # --- Pow ---
        if t == 'Pow' and Is_Pow(expr):
            base, exp = expr.args
            base2 = _map_expr(base, f, *types)
            exp2  = _map_expr(exp,  f, *types)
            return Pow(base2, exp2)

        # --- Root ---
        if t == 'Root' and Is_Root(expr):
            base = insideroot(expr)
            base2 = _map_expr(base, f, *types)
            exp = expr.args[1]
            return Pow(base2, exp)

    return expr   
def setBag(**kwargs):
    return {Symbol(k): v for k, v in kwargs.items()} 

    
class MyEq:
    def __init__(self, ksym, name='', var=x, var2=y, var3=z, type='P',vfunc=[], show=True,render=True):
          
        
        self.type = type
        self.vfunc = vfunc
        # Convertir a expresión sympy si es string
        if isinstance(ksym, str):
            self.ksym = parse_expr(ksym)
        elif isinstance(ksym, MyEq):
            self.ksym = ksym.ksym   
        else:
            self.ksym = ksym
            
        self.name = name
        if name=='alpha':
            name='α'
        if name=='beta':
            name='β'    
        self.var = var
        try:
            self.dvar=get_diff_var(var)
        except:
            self.dvar='d'+str(var)    
        self.var2 = var2
        self.var3 = var3    
        self.showroot=False 
        self.render=render
        if show:
            self.s()
    
    def s(self):
        if self.render:
            latex_str = latex_freeze(self.ksym)
            if self.showroot:
                latex_str = niceroot(self.ksym)
            sR = self.name + ' = '  
            display(Math(sR + latex_str))
 
 
    def _repr_latex_(self):
        """Se ejecuta automáticamente cuando escribes P solo en Jupyter"""
        sR = self.name + ' = ' if self.name else ''
        return f"${sR}{latex(self.ksym)}$"
    def _sympy_(self):
        return self.ksym
    # MÉTODOS DE OPERACIONES QUE MODIFICAN LA INSTANCIA ACTUAL
    def Add(self, other,show=True):
        """Convierte P a P + other"""
        if isinstance(other, MyEq):
            self.ksym = self.ksym + other.ksym
        else:
            self.ksym = self.ksym + other
        if show:
            self.s()
         
    
    def Sub(self, other,show=True):
        """Convierte P a P - other"""
        if isinstance(other, MyEq):
            self.ksym = self.ksym - other.ksym
        else:
            self.ksym = self.ksym - other
        if show:
            self.s()
    def Substrac(self, other,show=True):
        """Convierte P a P - other"""
        if isinstance(other, MyEq):
            self.ksym = self.ksym - other.ksym
        else:
            self.ksym = self.ksym - other
        if show:
            self.s()         
    
    def Mul(self, other,show=True):
        """Convierte P a P * other"""
        if isinstance(other, MyEq):
            self.ksym = self.ksym * other.ksym
        else:
            self.ksym = self.ksym * other
        if show:
            self.s()
         
    def mulupdown(self,expr,*args,show=True):
        p1=self.numer
        p2=self.denom
        P1=p1*expr
        P2=p2*expr
        if 'expand' in args:
            P1=expand(P1)
            P2=expand(P2)      
        if 'simplify' in args:
            P1=simplify(P1)
            P2=simplify(P2)
        if 'factor' in args:
            P1=factor(P1)
            P2=factor(P2)    
        self.ksym=sdiv(P1,p2)
        if show:
            self.s()

        
    
    def Div(self, other,show=True):
        """Convierte P a P / other"""
        if isinstance(other, MyEq):
            self.ksym = cf(self.ksym,other.ksym)
        else:
            self.ksym = cf(self.ksym,other)
        if show:
            self.s()
         
    
    def Pow(self, other,show=True):
        """Convierte P a P ** other"""
        if isinstance(other, MyEq):
            self.ksym = self.ksym ** other.ksym
        else:
            self.ksym = self.ksym ** other
        if show:
            self.s()
        
    def Rpow(self, ee,show=True):
        nee=cf(1,ee)
        self.ksym=Pow(self.ksym,nee)
        if show:
            self.s()    
         
    def fullhelp(self):
        """
        🚀 MYEQ - CLASE COMPACTA PARA EXPRESIONES MATEMÁTICAS 🚀
    
        📥 CREACIÓN:
            P = MyEq(3*x + 2*y, 'P')           # Con nombre
            Q = MyEq('sin(x) + cos(y)')         # Desde string
            R = MyEq(x**2 - 1, show=False)      # Sin mostrar
    
        👀 VISUALIZACIÓN:
            P                  # Muestra automáticamente
            P.s()              # Muestra explícitamente
    
        🔢 OPERACIONES QUE MODIFICAN:
            P.Add(5)           # P = P + 5
            P.Sub(2*x)         # P = P - 2*x  
            P.Mul(3)           # P = P * 3
            P.Div(2)           # P = P / 2
            P.Pow(2)           # P = P ** 2
    
        ➕ OPERACIONES QUE CREAN NUEVAS:
            Q = P + 10         # Nueva ecuación Q
            R = 3 * P          # Multiplicación
            S = P / 2          # División
    
        📊 EVALUACIÓN:
            P(4)               # P.subs(x, 4)
            P(2, 3)            # P.subs({x:2, y:3})
            P(x=1, y=2)        # Sustitución por nombre
    
        🧮 CÁLCULO:
            P.diff(x)          # Derivada respecto a x
            P.diff()           # Derivada con variable por defecto
            P.integrate(x)     # Integral respecto a x
            P.expand()         # Expansión
            P.simplify()       # Simplificación
    
        🔍 ANÁLISIS DE EXPRESIÓN:
            P.args()           # Muestra estructura completa
            P.args(0, 1)       # Accede a sub-argumento [0][1]
            P.numer            # Numerador
            P.denom            # Denominador
    
        ✅ VERIFICACIÓN (is_*):
            P.is_Add           # ¿Es suma?
            P.is_Mul           # ¿Es multiplicación?
            P.is_Symbol        # ¿Es símbolo?
            P.is_real          # ¿Es real?
            P.is_positive      # ¿Es positivo?
    
        📈 MÉTODOS SYMPY AUTOMÁTICOS:
            CUALQUIER método de SymPy funciona:
            P.factor(), P.series(), P.limit(), 
            P.solve(), P.evalf(), etc.
    
        🎯 EJEMPLOS RÁPIDOS:
            P = MyEq(x**2 + 2*x + 1).factor()     # (x + 1)**2
            deriv = MyEq(sin(x)).diff(x)           # cos(x)
            valor = MyEq(x + y)(1, 2)              # 3
    
        ⚡ ACCESO DIRECTO A SYMPY:
            P.ksym             # Expresión SymPy pura
        """
        print(self.fullhelp.__doc__)

    def help(self=None):
        """
        🚀 MYEQ - REFERENCIA RÁPIDA 🚀
        
        Creación:    MyEq(expr, nombre)
        Visualizar:  P, P.s()
        Operaciones: P.Add(val), P.Sub(val), P.Mul(val), P.Div(val), P.Pow(val)
        Evaluación:  P(4), P(x=1, y=2)
        Cálculo:     P.diff(), P.integrate(), P.expand(), P.simplify()
        Análisis:    P.args(), P.numer, P.denom
        Verificación: P.is_Add, P.is_real, etc.
        SymPy:       CUALQUIER método de SymPy funciona automáticamente!
        """
        print(self.help.__doc__ if self else MyEq.help.__doc__)
     
    # MÉTODOS PARA OPERACIONES QUE RETORNAN NUEVAS INSTANCIAS
    def __mul__(self, other):
        if isinstance(other, MyEq):
            return MyEq(self.ksym * other.ksym, show=False)
        return MyEq(self.ksym * other, show=False)
    
    def __rmul__(self, other):
        return MyEq(other * self.ksym, show=False)
    
    def __add__(self, other):
        if isinstance(other, MyEq):
            return MyEq(self.ksym + other.ksym, show=False)
        return MyEq(self.ksym + other, show=False)
    
    def __radd__(self, other):
        return MyEq(other + self.ksym, show=False)
    
    def __sub__(self, other):
        if isinstance(other, MyEq):
            return MyEq(self.ksym - other.ksym, show=False)
        return MyEq(self.ksym - other, show=False)
    
    def __rsub__(self, other):
        return MyEq(other - self.ksym, show=False)
    
    def __truediv__(self, other):
        if isinstance(other, MyEq):
            return MyEq(self.ksym / other.ksym, show=False)
        return MyEq(self.ksym / other, show=False)
    
    def __rtruediv__(self, other):
        return MyEq(other / self.ksym, show=False)
    
    def __pow__(self, other):
        if isinstance(other, MyEq):
            return MyEq(self.ksym ** other.ksym, show=False)
        return MyEq(self.ksym ** other, show=False)
 
    def __radd__(self, other): return MyEq(other + self.ksym, show=False)
    def __rmul__(self, other): return MyEq(other * self.ksym, show=False)
    def __rsub__(self, other): return MyEq(other - self.ksym, show=False)
    def __rtruediv__(self, other): return MyEq(other / self.ksym, show=False)

        
    # EVALUACIÓN CON PARÉNTESIS
    def __call__(self, *args, **kwargs):
        """Permite P(4), P(x=4), P('float',x=4)"""

        # --- modo float ---
        float_mode = 'float' in args
        args = tuple(a for a in args if a != 'float')

        if args and kwargs:
            raise ValueError("Usa solo args o kwargs, no ambos")

        # --- sin argumentos ---
        if not args and not kwargs:
            return self.ksym

        # --- args posicionales ---
        if args:
            variables = [self.var, self.var2, self.var3]
            subs_dict = {}

            for var, val in zip(variables, args):
                if var is not None:
                    subs_dict[var] = val

            kres = self.ksym.subs(subs_dict)

        # --- kwargs ---
        else:
            kres = self.ksym.subs(kwargs)

        # --- float opcional ---
        if float_mode:
            try:
                kres = float(N(kres))
            except:
                kres = N(kres)

        return kres
    
    # MÉTODOS SYMPY DELEGADOS
    def Rationalize(self,expr):
        expr=sqrt(expr)
        p1,p2=numer(self.ksym),denom(self.ksym)
        p1=p1*expr
        p2=p2*expr
        p2=rsimplify(sqrt((simplify(p2*p2))))
        self.ksym=cf(p1,p2)
        self.s()
    
    def diff(self, *args,var=''):
        """Derivada con variable por defecto o variables especificadas"""
        if var=='':
            var=self.var
        expr=self.ksym
        dexpr=diff(expr,var)
        
        if len(args) == 0:
            return  dexpr
        else:
            return MyEq(dexpr,args[0])
    
    def integrate(self, *args):
        """Integral con variable por defecto o variables especificadas"""
        if len(args) == 0:
            return integrate(self.ksym, self.var)
        else:
            return integrate(self.ksym, *args)
    
    def expand(self,show=True):
        """Expande la expresión"""
        self.ksym=expand(self.ksym)
        if show:
            self.s()
    
    def simplify(self,show=True):
        """SImplify  la expresión"""
        self.ksym=simplify(self.ksym)
        if show:
            self.s()
    def factor(self,show=True):
        """SImplify  la expresión"""
        self.ksym=factor(self.ksym)
        if show:
            self.s()        
        
    def rsimplify(self,show=True):
        self.ksym=rsimplify(self.ksym)
        if show:
            self.s()
            
    def reducecero(self,show=True):
        self.ksym=reducecero(self.ksym)
        if show:
            self.s()
        
    def tsimplify(self,show=True):
        self.ksym=trigsimp(self.ksym)
        if show:
            self.s()
    def texpand(self,show=True):
        self.ksym=expand_trig(self.ksym)
        if show:
            self.s()
  
    def sin2cos(self, ang,show=True):
        self.ksym = sin2cos(self.ksym,ang)
        if show:
            self.s() 
    def cos2sin(self, ang,show=True):
        self.ksym = cos2sin(self.ksym,ang)
        if show:
            self.s()
        
    # REPRESENTACIONES
    def __str__(self):
        return str(self.ksym)
    
    def __repr__(self):
        return f"MyEq({self.ksym})"
    
    # Delegar automáticamente todos los métodos de sympy no implementados
    def __getattr__(self, name):
        if hasattr(self.ksym, name):
            return getattr(self.ksym, name)
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

    @property
    def numer(self):
        return self.ksym.as_numer_denom()[0]
    
    @property
    def denom(self):
        return self.ksym.as_numer_denom()[1]

    @property
    def nu(self):
        return numer(self.ksym)
    
    @property
    def de(self):
        return denom(self.ksym)      

    ### arguments  ###
    
    def args(self, *args, deep=4, format='list'):
        """
        Navega por los argumentos de la expresión
        - Sin argumentos: muestra la estructura
        - Con índices: retorna el sub-argumento especificado
        """
        if len(args) == 0:
            self._showarglist(self.ksym, deep=deep, format=format)
        else: 
            kres = self.ksym
            for i in args:
                kres = kres.args[i]
            return kres
    
    def _showarglist(self, expr, deep=4, format='list', current_depth=0):
        """Muestra la estructura de argumentos de la expresión"""
        if current_depth >= deep:
            return
        
        indent = "  " * current_depth
        if format == 'list':
            print(f"{indent}{type(expr).__name__}: {expr}")
            if hasattr(expr, 'args') and expr.args:
                for i, arg in enumerate(expr.args):
                    print(f"{indent}  [{i}]: ", end="")
                    self._showarglist(arg, deep, format, current_depth + 1)
 
    def loadBag(self, Bag):
        expr=self.ksym
        used = {k:v for k,v in Bag.items() if k in expr.free_symbols}
        return expr.subs(used)
        
    def set(self, *args, evaluate=True, **kwargs):
        """
        Sustituye variables y MODIFICA permanentemente la expresión.
        
        Ejemplos:
            P.set(x, z)          # Reemplaza x por z permanentemente
            P.set(MyEq)  # Reemplaza value MyEq in name in MyEq
            P.set(x=z, y=w)      # Reemplaza con argumentos nombrados
        """
        kres=self.ksym
        if len(kwargs)>0:
            self.ksym=real_subs(self.ksym,**kwargs)
        if len(args)==2:     
            try:
                self.ksym=self.ksym.subs(args[0],args[1])
                 
            except:
                pass
        if len(args)==1:
            data=args[0]
            if isinstance(data, MyEq):
                valor=data.ksym
                name=data.name
                self.ksym=self.ksym.subs(str(name),valor)
            from lib_MyEqEq import MyEqEq    
            if isinstance(data, MyEqEq):
                p1=unisymbols(data.e1.ksym)
                p2=data.e2.ksym
                try:
                    self.ksym=self.ksym.subs(p1,p2)
                except:
                    self.ksym=self.ksym.subs(str(p1),p2)    
        self.s()  # Mostrar el nuevo resultado

    def limit(self,*args):
        ksym=self.ksym
        var=self.var
        op=""
        valr=args[0]
        if len(args)==2:
            op=args[1]
            return limit(ksym, var, valr, dir=op)
        else:
            return limit(ksym, var, valr)
 
    def primefactor(self,show=True):
            kres=self.ksym
            kres=primefactor(kres)
             
            self.ksym=kres
            if show:
                self.s() 
    def getbase(self):
        return getbase(self.ksym)
    def getexpo(self):
        return getexpo(self.ksym)        
    def insideroot(self):
        return insideroot(self.ksym)  

    ### Solve ### 
    
    def alone(self,var):
        var2=symbols(self.name)
        expr=var2-self.ksym
        e1=MyEq(expr,'e1',show=False)
        kres=e1.solve(var)
        return kres
        
    def solve(self,ivar,*args,**kwargs):
        self.unisymbols()
        sdone=False
        svmain=False
        expr=self.ksym
        if type(ivar)==str:
            sdone=True
        if len(kwargs)>0:
            vecs=[]
            vecv=[]
            for key, value in kwargs.items():
                vecs.append(key)
                vecv.append(value)
            if self.name in vecs:
                for svar,sval in zip(vecs,vecv):
                    if svar==self.name:
                        expr=expr-sval
            expr=real_subs(expr,**kwargs)
        svar=str(ivar)
        kres=solve(expr,svar)
            
        if 'float' in args:
            kres2=[]
            for data in kres:
                try:
                    kres2.append(float(data))
                except:
                    kres2.append(data)
            kres=kres2
            
        if 'positive' in args:
            kres2=[]
            for data in kres:
                if signo(data)>=0:
                    kres2.append(simplify2(data))
            kres=kres2
        if 'noimg' in args:
            kres2=[data for data in kres if not 'I' in str(data)]
            kres=kres2
        
        if len(kres)==1:
            if sdone:
                return MyEq(kres[0],svar)
            else:
                return kres[0]
        else:
            if not sdone:
                return kres
            else:
                evec=[]
                cc=1
                for i in range(len(kres)):
                    evec.append(MyEq(kres[i],svar+str(i+1)))
                return evec

    ### Exponenciales Funcion ###
    def mulexpo(self,show=True):
        """Aplica mulexpo a la expresión"""
        self.ksym = mulexpo(self.ksym)
        if show:
            self.s()
 
    
    def lexpand(self,show=True):
        """Expande logaritmos"""
        self.ksym = _map_expr(self.ksym, lexpand,'Add','Mul','Pow','Root')
        if show:
            self.s()
 
    
    def lfactor(self,joinBase=False,show=True):
        """Factoriza logaritmos"""
        self.ksym = lfactor(self.ksym,joinBase=joinBase)
        if show:
            self.s()
 
    
    def lexponent(self,show=True):
        """Trabaja con exponentes logarítmicos"""
        self.ksym = lexponent(self.ksym)
        if show:
            self.s()
 
    
    def lcombine(self,show=True):
        """Expande logaritmos"""
        self.ksym = _map_expr(self.ksym, lcombine,'Add','Mul','Pow','Root')
        if show:
            self.s()
 
    
    def pow2powpow(self,show=True):
        """Convierte potencias a formato powpow"""
        self.ksym = pow2powpow(self.ksym)
        if show:
            self.s()
 
    
    def powexpand(self,show=True):
 
        self.ksym = powexpand(self.ksym)
        if show:
            self.s()
 
    
    def div2mulexp(self,show=True):
        """Convierte divisiones a multiplicaciones con exponentes negativos"""
        self.ksym = div2mulexp(self.ksym)
        if show:
            self.s()
 
    
    def joinbase(self,show=True):
        """Une bases comunes"""
        self.ksym = joinbase(self.ksym)
        if show:
            self.s()
 
    
    def disjoinbase(self,show=True):
        """Separa bases comunes"""
        self.ksym = disjoinbase(self.ksym)
        if show:
            self.s()
 
    
    def positivexpo(self,show=True):
        """Convierte exponentes a formato positivo"""
        self.ksym = positivexpo(self.ksym)
        if show:
            self.s()
 
    
    def rexpand(self,show=True):
        """Expansión radical"""
        self.ksym = rexpand(self.ksym)
        if show:
            self.s()
 
    def expofactor(self ,show=True):
        """Factoriza exponentes"""
        self.ksym = expofactor(self.ksym)
        if show:
            self.s()
 
    
    def basefactor(self,show=True):
        """Factoriza bases"""
        self.ksym = basefactor(self.ksym)
        if show:
            self.s()
 
    
    def expoexpand(self,show=True):
        """Expande exponentes"""
        self.ksym = expoexpand(self.ksym)
        if show:
            self.s()
 
    
    def baseexpand(self,show=True):
        """Expande bases"""
        self.ksym = baseexpand(self.ksym)
        if show:
            self.s()
 
    
    def exposimplify(self,show=True):
        """Simplifica exponentes"""
        self.ksym = exposimplify(self.ksym)
        if show:
            self.s()
 
    
    def basesimplify(self,show=True):
        """Simplifica bases"""
        self.ksym = basesimplify(self.ksym)
        if show:
            self.s()
 
    
    def add2numer(self,fac,show=True):
        p1=numer(self.ksym)
        p2=denom(self.ksym)
        p1=p1+fac
        self.ksym = p1/p2
        if show:
            self.s()		

    def add2denom(self,fac,show=True):
        p1=numer(self.ksym)
        p2=denom(self.ksym)
        p2=p2+fac
        self.ksym = p1/p2
        if show:
            self.s()

    def mul2numer(self,fac,show=True):
        p1=numer(self.ksym)
        p2=denom(self.ksym)
        p1=p1*fac
        self.ksym = p1/p2
        if show:
            self.s() 
           
    def mul2denom(self,fac,show=True):
        p1=numer(self.ksym)
        p2=denom(self.ksym)
        p2=p2*fac
        self.ksym = p1/p2
        if show:
            self.s()

    def pow2numer(self,fac,show=True):
        p1=numer(self.ksym)
        p2=denom(self.ksym)
        p1=p1**fac
        self.ksym = p1/p2
        if show:
            self.s() 
    def exponencial(self,show=True):
        self.ksym=exp(self.ksym)
        if show:
            self.s()      
    def pow2denom(self,fac,show=True):
        p1=numer(self.ksym)
        p2=denom(self.ksym)
        p2=p2**fac
        self.ksym = p1/p2
        if show:
            self.s()

    def xfactor(self,fac,show=True):
        kres=factorsec(self.ksym,fac)
        self.ksym=kres
        if show:
            self.s()

    def factorsec(self,fac,show=True):
        kres=factorsec(self.ksym,fac)
        self.ksym=kres
        if show:
            self.s() 

            
        ### make all ###            
    def dothis(self, action, *args, **kwargs):
        """
        Ejecuta CUALQUIER método que ya esté cargado en memoria.
        SIN verificar, SIN diccionarios, DIRECTAMENTE.
        """
        # ¡DIRECTO AL GRANO!
        method = getattr(self, action)
        return method(*args, **kwargs)   
    
        
    def domain(self,var=''):
        from sympy import S
        from sympy.calculus.util import continuous_domain
        if var=='':
            var=self.var
        kres=continuous_domain(self.ksym,var,S.Reals)
        return kres    
    def range(self,var='',x1=-1,x2=1):
        from sympy import Interval, Symbol, S, exp, log, pi, sqrt, sin, tan
        from sympy.calculus.util import function_range
        if var=='':
            var=self.var
        kres=function_range(self.ksym, var, Interval(x1, x2))
        strx1=str(x1)
        if x1==-oo:
            strx1='-∞'
        strx2=str(x2)
        if x2==oo:
            strx2='∞'
        mdisplay('for interval <'+ strx1 +','+strx2+'>')
        return kres
    
    
    def numersimplify(self, show=True):
        self.ksym = _map_expr(self.ksym, numersimplify)
        if show:
            self.s()
    def numerfactor(self, show=True):
        self.ksym = _map_expr(self.ksym, numerfactor)
        if show:
            self.s()
    def numerexpand(self, show=True):
        self.ksym = _map_expr(self.ksym, numerexpand)
        if show:
            self.s()            
    def denosimplify(self, show=True):
        self.ksym = _map_expr(self.ksym, denosimplify)
        if show:
            self.s()
    def denofactor(self, show=True):
        self.ksym = _map_expr(self.ksym, denofactor)
        if show:
            self.s()
    def denoexpand(self, show=True):
        self.ksym = _map_expr(self.ksym, denoexpand)
        if show:
            self.s()            
     
             
     
    def float(self):
        try:
            return float(self.ksym)
        except:
            return self.ksym 
         
     
    def factoriza(self,xfac,*args,show=True):
        expr=self.ksym
        kres=factoriza(expr,xfac,*args)
        self.ksym=kres
        if show:
         self.s()     
        
    def unisymbols(self):
        kres=self.ksym
        kres=unisymbols(kres)
        self.ksym=kres
        
    def rpart(self):
        expr=self.ksym
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
    def ipart(self):
        expr=self.ksym
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
    def imodule(self):
        expr=self.ksym
        rp=rpart(expr)
        ip=ipart(expr)
        kres=sqrt(rp*rp+ip*ip)
        return simplify(kres)

        
        
    def simplifyimg(self,show=True):
        expr=self.ksym
        expr=simplifyimg(expr)
        self.ksym= expr 
        if show:
            self.s()

        
def niceroot(expr):
    """
    Versión mejorada que detecta patrones de exponentes simbólicos
    """
    from sympy import Pow, latex, Add, Mul, Integer
    
    # Función auxiliar para detectar exponentes fraccionarios simbólicos
    def is_symbolic_fraction(power_expr):
        """Detecta si es x^(a/b) donde a y b son simbólicos"""
        if not isinstance(power_expr, Pow):
            return False, None, None
        
        expo = power_expr.exp
        
        # Patrón 1: x^(a/b) donde a/b es Rational (ya cubierto)
        if expo.is_Rational:
            num, den = expo.as_numer_denom()
            return (den != 1, num, den)
        
        # Patrón 2: x^(a * 1/b) que es x^(a/b)
        if expo.is_Mul:
            numerator = 1
            denominator = 1
            for arg in expo.args:
                if arg.is_Pow and arg.exp == -1:  # 1/b
                    denominator = arg.base
                elif arg.is_Number:
                    numerator *= arg
                else:
                    numerator = arg  # asumimos que es el numerador
            
            if denominator != 1:
                return True, numerator, denominator
        
        # Patrón 3: x^(1/b)
        if expo.is_Pow and expo.exp == -1:
            return True, 1, expo.base
        
        return False, None, None
    
    # Procesar la expresión
    is_frac, num, den = is_symbolic_fraction(expr)
    if is_frac:
        base = expr.base
        if num == 1:
            return r'\sqrt[' + latex(den) + r']{' + latex(base) + r'}'
        else:
            return r'\sqrt[' + latex(den) + r']{' + latex(base**num) + r'}'
    
    # Resto de la función igual...
    elif isinstance(expr, Add):
        terms = [niceroot(term) for term in expr.args]
        return ' + '.join(terms)
    
    elif Is_Div(expr):
        num_str = niceroot(numer(expr))
        den_str = niceroot(denom(expr))
        return r'\frac{' + num_str + r'}{' + den_str + r'}'
    
    elif isinstance(expr, Mul):
        factors = [niceroot(factor) for factor in expr.args]
        return ' '.join(factors)
    
    else:
        return latex(expr)
        
        
    def unisymbols(self):
        kres=self.ksym
        kres=unisymbols(kres)
        self.ksym=kres    
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
        
def sendEqshow(data,name):
    if name!='':
        return MyEq(data,name)
    else:
        return data        
        
def out2data(data,name):
    if name!='':
        return MyEq(data,name)
    else:
        return data        