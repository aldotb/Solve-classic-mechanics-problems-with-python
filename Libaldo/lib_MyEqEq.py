from sympy import *
from lib_Variables import *
from lib_Mathematica import *
from lib_Mathbasic import *
from lib_Algorith import *
 
import copy
import numpy as np


from sympy import *
from IPython.display import display, Math

from lib_MyEq import *

inesym=['<','<=','=','>=','>','<>']
inesym2=['<','≤','=','≥','>','≠']
def getdata(expr):
    kres=expr
    if type(expr)==MyEq:
        kres=expr.ksym
    return kres
def get_diff_var(var):
    # find therespetive differential varible for x,y,z....
    return dLvar[Lvar.index(var)]

def obj2expr(obj):
    """Convierte cualquier objeto a MyEq"""
    if isinstance(obj, MyEq):
        return obj.ksym
    elif isinstance(obj, str):
        return  parse_expr(obj) 
    else:
        return obj
    
class MyEqEq:
    def __init__(self, *args, var=x, show=True):
        """
        Representa la ecuación: obj1 = obj2 con todos los métodos de MyEq
        """
        # Convertir obj1 y obj2 a MyEq si no lo son
        self.type='EQ'
        self.ssym='='
        self.e1 = MyEq(obj2expr(args[0]),'e1',var=var,show=False)
        self.var=var
        if len(args)==2:
            self.e2 = MyEq(obj2expr(args[1]),'e2',var=var,show=False)
        else:
            
            self.type='IQ'
            self.ssym =  get_inesym(args[1])   
            self.e2 = MyEq(obj2expr(args[2]),'e2',var=var,show=False)
            
        self.showroot = False  # ← Inicializar showroot como False
        
        # La ecuación como expresión SymPy
        self.ksym = self.e1.ksym - self.e2.ksym
        self.eQ=Eq(self.e1.ksym,self.e2.ksym)
 
        if show:
            self.s()
    
    def _to_myeq(self, obj):
        """Convierte cualquier objeto a MyEq"""
        if isinstance(obj, MyEq):
            return obj.ksym
        elif isinstance(obj, str):
            return  parse_expr(obj) 
        else:
            return obj 
            
            
    def __call__(self,  **kwargs):
        """Permite P(4) o P(x=4, y=5) para evaluar la expresión"""
        if len(kwargs)>0:
            p1=self.e1.ksym
            P1=real_subs(p1,**kwargs)
            p2=self.e2.ksym
            P2=real_subs(p2,**kwargs)
            Q2=MyEqEq(P1,P2,show=false)
            Q2.s()
    def s(self,*args,kshow=True):
        self.e1.ksym=self.e1.ksym
        self.e2.ksym=self.e2.ksym
        p1=self.e1.ksym
        p2=self.e2.ksym
        ps=self.ssym
                 
        if kshow:
            if len(args)>0:
                sexpr2=str(args[0])
                display(Math(latex(p1)+' '+ps+' '+latex(p2)+', '+sexpr2))
            else:    
                display(Math(latex(p1)+' '+ps+' '+latex(p2)))
             
    @property     
    def L(self):
        return  self.e1.ksym 
    @property    
    def R(self):
        return  self.e2.ksym  
        
    def display_root(self,show=True):
        """Activa modo raíces bonitas y muestra"""
        self.showroot = True
        return self.s()
    
    def display_normal(self,show=True):
        """Desactiva modo raíces bonitas y muestra"""
        self.showroot = False
        return self.s()
    
    # 🔥 MÉTODOS QUE APLICAN A AMBOS LADOS
    
    def Add(self, other,show=True):
        """Suma a ambos lados: (left + other) = (right + other)"""
        p1,p2=obj2exprMQ(other) 
        self.e1.ksym=self.e1.ksym+p1
        self.e2.ksym=self.e2.ksym+p2
        self.ksym = self.e1.ksym - self.e2.ksym
        if show:
            self.s()
   
    
    def Sub(self, other,show=True):
        """Resta a ambos lados: (left - other) = (right - other)"""
        p1,p2=obj2exprMQ(other)
        self.e1.ksym=self.e1.ksym-p1
        self.e2.ksym=self.e2.ksym-p2
        self.ksym = self.e1.ksym - self.e2.ksym
        if show:
            self.s()
      
    
    def Mul(self, other,show=True):
        """Multiplica ambos lados: (left * other) = (right * other)"""
        p1,p2=obj2exprMQ(other)
        self.e1.ksym=self.e1.ksym*p1
        self.e2.ksym=self.e2.ksym*p2
        self.ksym = self.e1.ksym - self.e2.ksym
        if show:
            self.s()
    
    def Div(self, other,show=True):
        """Divide ambos lados: (left / other) = (right / other)"""
        p1,p2=obj2exprMQ(other)
        self.e1.ksym=float2int(self.e1.ksym/p1)
        self.e2.ksym=float2int(self.e2.ksym/p2)
        self.ksym =  self.e1.ksym - self.e2.ksym 
        if show:
            self.s()
       
    def Pow(self, *args,show=True):
        if len(args)==0:
            other=2
        else:
            other=args[0]          
        """Divide ambos lados: (left / other) = (right / other)"""
        self.e1.ksym=self.e1.ksym**other
        self.e2.ksym=self.e2.ksym**other
        self.ksym = self.e1.ksym - self.e2.ksym
        if show:
            self.s()

    def Rpow(self, other='',show=True):
        """Divide ambos lados: (left / other) = (right / other)"""
        if other=='':
            self.e1.ksym=sqrt(self.e1.ksym)
            self.e2.ksym=sqrt(self.e2.ksym)
        else:    
            self.e1.ksym=self.e1.ksym**cf(1,other)
            self.e2.ksym=self.e2.ksym**cf(1,other)
        if show:
            self.s()
        
    def set(self, *args,show=True, **kwargs):
        P1=unisymbols(self.e1.ksym)
        P2=unisymbols(self.e2.ksym)
        if len(kwargs)>0:
            P1=real_subs(P1,**kwargs)
            P2=real_subs(P2,**kwargs)
            self.e1.ksym=P1
            self.e2.ksym=P2

        if len(args)==0:
            pass 
        elif len(args)==1 and type(args[0])==MyEqEq:
            Q2=args[0]
            sexpr1=str(Q2.L)
            sexpr2=str(Q2.R)
            sP1=str(P1)
            sP2=str(P2)
            sP1=sP1.replace(sexpr1,'('+sexpr2+')')
            sP2=sP2.replace(sexpr1,'('+sexpr2+')')
            self.e1.ksym=parse_expr(sP1)
            self.e2.ksym=parse_expr(sP2)
            
     
        elif len(args)==1 and type(args[0])==MyEq:
            var=args[0].name
            value=args[0].ksym
            
            self.e1.ksym=P1.subs(var,value)
            self.e2.ksym=P2.subs(var,value)
 
        else:
            var=args[0]
            value=args[1]
            
            self.e1.ksym=P1.subs(var,value)
            self.e2.ksym=P2.subs(var,value)    
        if show:
            self.s()
            
    def setL(self,*args,show=True,**kwargs):
        if len(args)==1:
            kres=args[0]
        else:    
            kres=self.e1.ksym
            kres=real_subs(kres,**kwargs)
        self.e1.ksym=kres        
        if show:
            self.s()
            
    def setR(self,*args,show=True,**kwargs):
        if len(args)==1:
            kres=args[0]
        else:    
            kres=self.e2.ksym
            kres=real_subs(kres,**kwargs)
        self.e2.ksym=kres        
        if show:
            self.s() 

            
    def subs(self, p1,p2,show=True):
        """Aplica sustituciones a ambos lados"""
        P1=self.e1.ksym
        P2=self.e2.ksym
        P1=P1.subs(p1,p2)
        P2=P2.subs(p1,p2)
        self.e1.ksym=P1
        self.e2.ksym = P2
        if show:
            self.s()
 
    def expand(self,OPS='RL',show=True):
        """Expande ambos lados"""
        if 'L' in OPS:
            self.e1.ksym=expand(self.e1.ksym)
        if 'R' in OPS:    
            self.e2.ksym=expand(self.e2.ksym)
        self.ksym = self.e1.ksym - self.e2.ksym
        if show:
            self.s()
         
    def baseexpand(self,OPS='RL',show=True):
        """Expande bases en ambos lados"""
        self.e1.ksym=baseexpand(self.e1.ksym)
        self.e2.ksym=baseexpand(self.e2.ksym)
        self.ksym = self.e1.ksym - self.e2.ksym
        if show:
            self.s()  
         
    def expandbase(self,OPS='RL',show=True):
        """Expande bases en ambos lados"""
        if 'L' in OPS:
            self.e1.ksym=baseexpand(self.e1.ksym)
        if 'R' in OPS:
            self.e2.ksym=baseexpand(self.e2.ksym)
        self.ksym = self.e1.ksym - self.e2.ksym
        if show:
            self.s()
        
    def factor(self,OPS='RL',show=True):
        """Factoriza ambos lados"""
        if 'L' in OPS:
            self.e1.ksym=factor(self.e1.ksym)
        self.e2.ksym=factor(self.e2.ksym)
        self.ksym = self.e1.ksym - self.e2.ksym
        if show:
            self.s()
            
    def lsimplify(self,show=True):
        self.e1.ksym =lsimplify(self.e1.ksym)
        self.e2.ksym =lsimplify(self.e2.ksym)
        if show:
            self.s() 
    def lcombine(self,show=True):
        self.e1.ksym =lcombine(self.e1.ksym)
        self.e2.ksym =lcombine(self.e2.ksym)
        if show:
            self.s()
    def lexpand(self,show=True):
        self.e1.ksym =lexpand(self.e1.ksym)
        self.e2.ksym =lexpand(self.e2.ksym)
        if show:
            self.s()        
    def lfactor(self,show=True):
        self.e1.ksym =lfactor(self.e1.ksym)
        self.e2.ksym                                                                                        =lfactor(self.e2.ksym)
        if show:
            self.s()   
    def basefactor(self,OPS='RL',show=True):
        """basefactor ambos lados"""
        if 'L' in OPS:
            self.e1.ksym=basefactor(self.e1.ksym)
        if 'R' in OPS:
            self.e2.ksym=basefactor(self.e2.ksym)
        
        self.ksym = self.e1.ksym - self.e2.ksym
        if show:
            self.s()
            
    def factorbase(self,OPS='RL',show=True):
        """basefactor ambos lados"""
        if 'L' in OPS:
            self.e1.ksym=basefactor(self.e1.ksym)
        if 'R' in OPS:
            self.e2.ksym=basefactor(self.e2.ksym)
        self.ksym = self.e1.ksym - self.e2.ksym
        if show:
            self.s()  
            
    def rsimplify(self,OPS='RL',show=True):
        if 'L' in OPS:
            self.e1.ksym=rsimplify(self.e1.ksym)
        if 'R' in OPS:
            self.e2.ksym=rsimplify(self.e2.ksym)    
        if show:
            self.s()
    def simplify(self,OPS='RL',show=True):
        
        """Expande bases en ambos lados"""
        if 'L' in OPS:
            self.e1.ksym=simplify(self.e1.ksym)
        if 'R' in OPS:
            self.e2.ksym=simplify(self.e2.ksym)
        self.ksym = self.e1.ksym - self.e2.ksym

        if isinstance(self.L,Pow) and isinstance(self.R,Pow):
            if getbase(self.L)==getbase(self.R):
                self.e1.ksym=getexpo(self.L)
                self.e2.ksym=getexpo(self.R)
            if getexpo(self.L)==getexpo(self.R):
                self.e1.ksym=getbase(self.L)
                self.e2.ksym=getbase(self.R)    
        self.addsimplify(show=False)
        self.mulsimplify(show=False)
        
        if show:
            self.s()
         
    def addsimplify(self,show=True):
        p1=self.L
        p2=self.R
        sdone=False
        if Is_Add(p1) or Is_Add(p2):
            if Is_Add(p1):
                V1=p1.args
            else:
                V1=[p1]
        
            if Is_Add(p2):
                V2=p2.args
            else:
                V2=[p2]

            kr1=0
            kr2=0
            for data in V1:
                if not data in V2:
                    kr1=kr1+data
            for data in V2:            
                if not data in V1:
                    kr2=kr2+data
            self.e1.ksym=kr1
            self.e2.ksym=kr2
            
        if show:
            self.s()

    def mulsimplify(self,show=True):
        p1=self.L
        p2=self.R
         
        if Is_Mul(p1) or Is_Mul(p2):
            if Is_Mul(p1):
                V1=p1.args
            else:
                V1=[p1]
        
            if Is_Mul(p2):
                V2=p2.args
            else:
                V2=[p2]
            V1n=[]
            V1d=[]
            for data in V1:
                if Is_Div(data):
                    V1d.append(denom(data))
                else:
                    V1n.append(data)
            V2n=[]
            V2d=[]
            for data in V2:
                if Is_Div(data):
                    V2d.append(denom(data))
                else:
                    V2n.append(data)
            P1n=1
            P1d=1
            P2n=1
            P2d=1
            if len(V1n)>0:
                if len(V2n)>0:
                    for data in V1n:
                        if not data in V2n:
                            P1n=P1n*data
                else:
                    P1n=prod(V1n)
            if len(V2n)>0:
                if len(V1n)>0:
                    for data in V2n:
                        if not data in V1n:
                            P2n=P2n*data
                else:
                    P2n=prod(V2n)
            if len(V1d)>0:
                if len(V2d)>0:
                    for data in V1d:
                        if not data in V2d:
                            P1d=P1d*data
                else:
                    P1d=prod(V1d)
            if len(V2d)>0:
                if len(V1d)>0:
                    for data in V2d:
                        if not data in V1d:
                            P2d=P2d*data
                else:
                    P2d=prod(V2d)
            if P1d==1:
                Pa=P1n
            else:
                Pa=cf(P1n,P1d)

            if P2d==1:
                Pb=P2n
            else:
                Pb=cf(P2n,P2d)
            self.e1.ksym=Pa
            self.e2.ksym=Pb
        if show:
            self.s()    
    def mulexpo(self,OPS='RL',show=True):
        """Aplica mulexpo a ambos lados"""
        if 'L' in OPS:
            self.e1.ksym=mulexpo(self.e1.ksym)
        if 'R' in OPS:
            self.e2.ksym=mulexpo(self.e2.ksym)
        self.ksym = self.e1.ksym - self.e2.ksym
        if show:
            self.s()
            
    def joinbase(self,OPS='RL',show=True):
        """Aplica joinbase"""
        if 'L' in OPS: 
            self.e1.ksym=joinbase(self.e1.ksym)
        if 'R' in OPS:
            self.e2.ksym=joinbase(self.e2.ksym)
        if show:
            self.s()    
        
    def disjoinbase(self,OPS='RL',show=True):
        """Aplica disjoinbase"""
        if 'L' in OPS:
            self.e1.ksym=disjoinbase(self.e1.ksym)
        if 'R' in OPS:
            self.e2.ksym=disjoinbase(self.e2.ksym)
        if show:
            self.s()
            
    def joinexpo(self,OPS='RL',show=True):
        """Aplica joinexpo"""
        if 'L' in OPS:
            self.e1.ksym=joinexpo(self.e1.ksym)
        if 'R' in OPS:
            self.e2.ksym=joinexpo(self.e2.ksym)
        if show:
            self.s()    
        
    def disjoinexpo(self,OPS='RL',show=True):
        """Aplica disjoinexpo"""
        if 'L' in OPS:
            self.e1.ksym=disjoinexpo(self.e1.ksym)
        if 'R' in OPS:
            self.e2.ksym=disjoinexpo(self.e2.ksym)
        if show:
            self.s()
            
    def joinbasemul(self,OPS='RL',show=True):
        """Aplica  joinbasemul"""
        if 'L' in OPS:
            self.e1.ksym= joinbasemul(self.e1.ksym)
        if 'R' in OPS:
            self.e2.ksym= joinbasemul(self.e2.ksym)
        if show:
            self.s()
 
            
    # ... AGREGAR TODOS LOS MÉTODOS QUE QUIERAS de la misma forma ...
    
    def solve(self, *args, **kwargs):
        """Resuelve la ecuación"""

        QQ=MyEq(self.L-self.R,'QQ',show=False)
        return QQ.solve(*args,**kwargs)
        
  
    def solveandset(self,expr,show=True,**kwargs):
        p1 =real_subs(self.e1.ksym,**kwargs)
        p2 =real_subs(self.e2.ksym,**kwargs)
    
        P=MyEq(p1-p2,'ee',show=False)
        kres=P.solve(expr,show=False)
        ee=MyEq(kres,str(expr))
        self.set(expr,kres,show=False)
        if show:
            self.s()
            
    def alone(self,expr,show=True,**kwargs):
        self.e1.ksym =real_subs(self.e1.ksym,**kwargs)
        self.e2.ksym =real_subs(self.e2.ksym,**kwargs)
        P=MyEq(self.e1.ksym-self.e2.ksym,'ee',show=False)
        kres=P.solve(expr,show=False)
        self.e1.ksym=expr
        self.e2.ksym=kres
        
        if show:
            self.s()
        
    def _repr_latex_(self):
        """Representación en LaTeX"""
        left_str = self.e1.name if self.e1.name else latex(self.e1.ksym)
        right_str = self.e2.name if self.e2.name else latex(self.e2.ksym)
        return f"${left_str} = {right_str}$"
    
    # 🔥 MÉTODO MÁGICO: delega automáticamente métodos no implementados
    def __getattr__(self, name):
        """Aplica métodos de MyEq a ambos lados automáticamente"""
        if hasattr(MyEq, name):
            def method_wrapper(*args, **kwargs):
                # Aplicar al lado izquierdo
                left_method = getattr(self.e1, name)
                left_method(*args, **kwargs)
                
                # Aplicar al lado derecho  
                right_method = getattr(self.e2, name)
                right_method(*args, **kwargs)
                
                # Actualizar ecuación
                self.ksym = self.e1.ksym - self.e2.ksym
                self.s()
                
            return method_wrapper
        raise AttributeError(f"'{type(self).__name__}' no tiene atributo '{name}'")
 
    def crosmul(self,show=True):
        p1=numer(self.e1.ksym)*denom(self.e2.ksym)
        p2=numer(self.e2.ksym)*denom(self.e1.ksym)
        self.e1.ksym=p1
        self.e2.ksym=p2
        if show:
            self.s()

    def all2R(self,show=True):
        self.e2.ksym=self.e2.ksym-self.e1.ksym
        self.e1.ksym=0
        if show:
            self.s()
    def all2L(self,show=True):
        self.e1.ksym=self.e1.ksym-self.e2.ksym
        self.e2.ksym=0       
        if show:
            self.s()
        
    def flip(self,show=True):

        k=self.e1.ksym
        self.e1.ksym=self.e2.ksym
        self.e2.ksym=k  
        if show:
            self.s()
        
    def switch(self,show=True):

        k=self.e1.ksym
        self.e1.ksym=-self.e2.ksym
        self.e2.ksym=-k 
        if show:
            self.s()
    def cross(self,show=True):
        p1=numer(self.e1.ksym)
        p2=denom(self.e1.ksym)
        p3=numer(self.e2.ksym)
        p4=denom(self.e2.ksym)
        if p4==1:
            P1=p1
        else:
            P1=p1*p4
        if p2==1:
            P2=p3
        else:
            P2=p3*p2
        self.e1.ksym=P1
        self.e2.ksym=P2 
        if show:
            self.s() 

        
    def args(self,*args,deep=2,format='list'):
        if len(args)==0:
            showarglist(self.e1.ksym,deep=deep,format=format,side='L')
            showarglist(self.e2.ksym,deep=deep,format=format,side='R')
        else: 
            kres=Eq(self.e1.ksym,self.e2.ksym)
            for i in args:
                kres=kres.args[i]
            return kres 
    
    def xfactor(self,factor,show=True):
        self.secfactor(factor,show=show)
    def factorsec(self,factor,show=True):
        self.secfactor(factor,show=show)        
    def secfactor(self,factor,show=True):
        p1=self.e1.ksym
        p2=self.e2.ksym
        try:
            self.e1.ksym=secfactor(p1,factor)
        except:
            pass
            
        try:
            self.e2.ksym=secfactor(p2,factor)
        except:
            pass    
        if show:
            self.s()
            
    def exponencial(self,show=True):
        self.e1.ksym=exp(self.e1.ksym)
        self.e2.ksym=exp(self.e2.ksym)
        if show:
            self.s()   
            
    def Integral(self,show=True):
        from lib_MyDiff import findvardiff, Integral2
        p1=self.L
        p2=self.R
        var1=findvardiff(p1)
        var2=findvardiff(p2)
        P1=Integral2(p1,var1)
        P2=Integral2(p2,var2)
        self.e1.ksym=P1
        self.e2.ksym=P2
        if show:
            self.s()
    def doitI(self,show=True):
        self.e1.ksym=self.e1.ksym.doit()
        self.e2.ksym=self.e2.ksym.doit()+C
        if show:
            self.s()  

    def diff(self,var1,var2,show=True):
        dvar1=get_diff_var(var1)
        dvar2=get_diff_var(var2)
        dp1=diff(self.e1.ksym,var1)*dvar1
        dp2=diff(self.e2.ksym,var2)*dvar2
        return MyEqEq(dp1,dp2,show=show)
        
    def changediff(self,var1,var2,expr,show=True):
        '''
           changediff(var to change, new var, var in func newvar)
           Q:dy=3*dx
           Q.changediff(x,w,2*w)
           
        '''
        dv1=get_diff_var(var1)
        dv2=get_diff_var(var2)
        dexpr=diff(expr,var2)
        self.set(dv1,dexpr*dv2,show=False)
        self.set(var1,expr,show=false)
        self.rsimplify(show=False)
        self.simplify(show=False)
        if show:
            self.s()        
def get_inesym(ssym):
    return inesym2[inesym.index(ssym)]
    
    
def MQ(expr1,expr2):
    if type(expr2)!=str:
        return MyEqEq(expr1,expr2)
    else:
        return MyEq(expr1,expr2)    
 
def obj2exprMQ(expr):
    if type(expr)==MyEq:
        return expr.ksym,expr.ksym
    elif type(expr)==MyEqEq:
        return expr.e1.ksym, expr.e2.ksym
    else:
        return expr,expr    