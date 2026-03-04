from sympy import *
from lib_Variables import *
from lib_Mathematica import *
from lib_Mathbasic import *
from lib_Algorith import *
 
import copy
import numpy as np


from sympy import *
from IPython.display import display, Math
 
def get_diff_var(var):
    # find therespetive differential varible for x,y,z....
    return dLvar[Lvar.index(var)] 
    
    
class MyEq:
    def __init__(self, ksym, name='', var=x, var2=y, var3=z, type='P',vfunc=[], show=True,render=True):
          
        
        self.type = type
        self.vfunc = vfunc
        # Convertir a expresión sympy si es string
        if isinstance(ksym, str):
            self.ksym = parse_expr(ksym)
        else:
            self.ksym = ksym
            
        self.name = name
        self.var = var
        self.dvar=get_diff_var(var)
        self.var2 = var2
        self.var3 = var3    
        self.showroot=False 
        self.render=render
        if show:
            self.s()
    
    def s(self):
        if self.render:
            latex_str = latex(self.ksym)
            if self.showroot:
                latex_str = niceroot(self.ksym)
            sR = self.name + ' = '  
            display(Math(sR + latex_str))
 
 
    def _repr_latex_(self):
        """Se ejecuta automáticamente cuando escribes P solo en Jupyter"""
        sR = self.name + ' = ' if self.name else ''
        return f"${sR}{latex(self.ksym)}$"
    
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
        """Permite P(4) o P(x=4, y=5) para evaluar la expresión"""
        if args and kwargs:
            raise ValueError("Usa solo argumentos posicionales O argumentos con nombre, no ambos")
        
        if args:
            # P(4) - sustituye la variable por defecto
            if len(args) == 1:
                return self.ksym.subs(self.var, args[0])
            else:
                # P(4, 5) - sustituye var, var2, var3 en orden
                substitutions = {}
                variables = [self.var, self.var2, self.var3]
                for i, value in enumerate(args):
                    if i < len(variables) and variables[i] is not None:
                        substitutions[variables[i]] = value
                return self.ksym.subs(substitutions)
        
        elif kwargs:
            # P(x=4, y=5) - sustituye variables específicas
            return self.ksym.subs(kwargs)
        
        else:
            # P() - retorna la expresión original
            return self.ksym
    
    # MÉTODOS SYMPY DELEGADOS
    def diff(self, *args):
        """Derivada con variable por defecto o variables especificadas"""
        if len(args) == 0:
            return  diff(self.ksym,self.var)
        else:
            return diff(self.ksym,self.var,*args)
    
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
            
    def tsimplify(self,show=True):
        self.ksym=trigsimp(self.ksym)
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
 
    
    def set(self, *args, evaluate=True, **kwargs):
        """
        Sustituye variables y MODIFICA permanentemente la expresión.
        
        Ejemplos:
            P.set(x, z)          # Reemplaza x por z permanentemente
            P.set({x: z, y: w})  # Reemplaza múltiples variables
            P.set(x=z, y=w)      # Reemplaza con argumentos nombrados
        """
        if args and kwargs:
            raise ValueError("Usa solo argumentos posicionales O argumentos con nombre")
        
        if len(args) == 2:
            # P.set(x, z) - sustituye variable por valor
            old_var, new_var = args
            if evaluate==False:
                sexpr=str(self.ksym)
                sexpr=sexpr.replace(str(old_var),str(new_var))
                kres=parse_expr(sexpr,evaluate=False)
                if type(kres)==Mul:
                    if kres.args[0]==1:
                        kres=kres.args[1]
                self.ksym=kres
            else:    
                self.ksym = self.ksym.subs(old_var, new_var)
        elif len(args) == 1 and isinstance(args[0], dict):
            # P.set({x: z, y: w}) - sustituye con diccionario
            self.ksym = self.ksym.subs(args[0], evaluate=evaluate)
        elif kwargs:
            # P.set(x=z, y=w) - sustituye con argumentos nombrados
            self.ksym = self.ksym.subs(kwargs, evaluate=evaluate)
        else:
            raise ValueError("Formato inválido. Usa: set(var, valor) o set({var: valor}) o set(var=valor)")
        
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
                    kres2.append(data)
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
        self.ksym = lexpand(self.ksym)
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
        """Combina logaritmos"""
        self.ksym = lcombine(self.ksym)
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
        
    def numesimplify(self,show=True):
        self.simplifynumer(show=show)
    def simplifynumer(self,show=True):
        expr=self.ksym
        kres=simplifynumer(expr)
        self.ksym=kres
        if show:
            self.s()
             
    def numefactor(self,show=True):
        self.factornumer(show=show)         
    def factornumer(self,show=True):
        expr=self.ksym
        kres=factornumer(expr)
        self.ksym=kres
        if show:
            self.s()
    def float(self):
        try:
            return float(self.ksym)
        except:
            return self.ksym 
    
    def numeexpand(self,show=True):
        self.expandnumer(show=show)    
    def expandnumer(self,show=True):
        expr=self.ksym
        kres=expandnumer(expr)
        self.ksym=kres
        if show:
            self.s()

    def denosimplify(self,show=True):
        self.simplifydenom(show=show)             
    def simplifydenom(self,show=True):
        expr=self.ksym
        kres=simplifydenom(expr)
        self.ksym=kres
        if show:
            self.s()
            
    def denofactor(self,show=True):
        self.factordenom(show=show)        
    def factordenom(self,show=True):
        expr=self.ksym
        kres=factordenom(expr)
        self.ksym=kres
        if show:
            self.s()
            
    def denoexpand(self,show=True):
        self.expanddenom(show=show)        
    def expanddenom(self,show=True):
        expr=self.ksym
        kres=expanddenom(expr)
        self.ksym=kres
        if show:
            self.s()    
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