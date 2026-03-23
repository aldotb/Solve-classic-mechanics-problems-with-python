from sympy import *
 
from lib_Variables import *
from lib_Mathbasic import *
from lib_MyEq import *
from lib_MyEqEq import * 

 

 
import copy
import numpy as np
import pandas as pd
 
 
# from lib_Func import *
import dill  # pip install dill --user
import pickle
 
 
 
def savework():
    dill.dump_session(filename='workspace.pkl')


def loadwork():
    dill.load_session(filename='workspace.pkl')

Lvar=[x, y, z, w, v, u, t, f, h, V, A,m,r,alpha,beta,theta]
dLvar=[dx, dy, dz, dw, dv, du, dt, df, dh, dV, dA,dm,dr,dalpha,dbeta,dtheta] 
    
# and to load the session again:
def creadiffvar(var):
    if var==alpha:
        return dalpha
        
    svar=str(var)
    sdvar='d'+svar
    return symbols(sdvar)
    
def creteintegral(expr,var=x,x1='',x2=''):
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
            

#direx,direy=symbols('\overrightarrow{x} \overrightarrow{y}')
 
def findvardiff(expr):
    lvar=list(expr.free_symbols)
    for data in lvar:
        if data in dLvar:
            return Lvar[dLvar.index(data)]
    for data in lvar:
        if len(str(data))==2 and str(data)[0]=='d':
            return  symbols(str(data)[1]) 
 
class MyIntg(MyEq):

    def __init__(self, expr, name, var=x, x1='', x2='', show=True, done=False):

        expr = obj2expr(expr)
        self.type = 'I'

        if var == x or var == y or var == z:
            if var == x:
                dvar = dx
            if var == y:
                dvar = dy
            if var == z:
                dvar = dz
        else:
            dvar = creadiffvar(var)

        self.dvar = dvar
        self.sdvar = str(dvar)

        expr = unisymbols(expr)

        try:
            expr = expr.subs(dvar, 1)
            expr = expr.subs(str(dvar), 1)
        except:
            pass

        self.ksym = expr
        self.name = name
        self.x1 = getdata(x1)
        self.x2 = getdata(x2)
        self.var = var

        self.Iksym = creteintegral(expr, var, x1=x1, x2=x2)

        self.varc = False
        self.done = False
        self.changv = []

        if show:
            self._display_integral()

    # ✅ DISPLAY CENTRALIZADO
    def _display_integral(self):

        from sympy import Symbol, latex
        from IPython.display import display, Math

        name_tex = latex(Symbol(str(self.name)))

        if self.ksym == 1:
            if self.x1 != '':
                expr_tex = rf"\int_{{{self.x1}}}^{{{self.x2}}} {latex(self.dvar)}"
            else:
                expr_tex = rf"\int {latex(self.dvar)}"
        else:
            expr_tex = latex(self.Iksym)

        display(Math(f"{name_tex} = {expr_tex}"))

    def __call__(self, done=False, **kwargs):
        if self.type == 'P':
            kres = unisymbols(self.ksym)
        if self.type == 'I':
            kres = unisymbols(self.Iksym)
        if len(kwargs) > 0:
            kres = real_subs(kres, **kwargs)
        return kres

    def __repr__(self):
        if self.type == 'P':
            kres = str(self.ksym)
        else:
            kres = str(self.Iksym)
        return kres

    def _latex(self, obj):
        return latex(self.ksym)

    def __str__(self):
        return self.__repr__()

    def rebuild(self):
        self.Iksym = creteintegral(
            unisymbols(self.ksym),
            var=unisymbols(self.var),
            x1=self.x1,
            x2=self.x2
        )

    def s(self):
        if self.type=='I':
            expr_tex = latex(self.Iksym)
            name_tex = latex(Symbol(self.name))
            display(Math(f"{name_tex} = {expr_tex}"))
        else:
            expr_tex = latex(self.ksym)
            name_tex = latex(Symbol(self.name))
            display(Math(f"{name_tex} = {expr_tex}"))

    def Add(self,obj):
        p1=self.Iksym
        if obj.type=='I':
            p2=obj.Iksym
        else:
            p2=obj2expr
         
        self.Iksym=p1+p1
        self.s()    
            
        
    def value(self):
        kres=self.Iksym.doit()
        return kres
        
    def doit(self,*args,show=True):
        if self.type=='P':
            pass
        else:
            Ires = self.Iksym 
            kres = Ires.doit()

            if 'float' in args:
                try:
                    kres = float(kres)
                except:
                    pass

            if 'positive' in args:
                kres = signo(kres)*kres
                    
            if 'C' in args:
                kres = kres + C1
            if 'C1' in args:
                kres = kres + C1
            if 'C2' in args:
                kres = kres + C2    

            # CONVERTIR OBJETO
            self.__class__ = MyEq
            if len(args)==1 and isinstance(args[0],str):
                self.__init__(kres, args[0], show=False)
            else:
                self.__init__(kres, self.name, show=False)

        if show:
            self.s()
  
     
    def insideIntg(self):
        return self.ksym
    
            
    def getdiffvar2(self):
        var=self.var
        svar=str(var)
        sdvar='d'+svar
        return symbols(sdvar) 
        
    def getintegrlaexpretion(self):
        kres=self.ksym
        diffvar=self.getdiffvar2()
        return kres*diffvar
        
    def transfordiff(self,var2,show=True,**kwargs):
        '''
        var2: new variablefunction diff
        **kwarg: dx=y*y*z, z=3...
        '''
        var1=self.var
        expr=self.ksym*self.dvar
        expr=real_subs(expr,**kwargs)
        self.ksym=expr.subs(var1,1)
        self.var=var2
        self.dvar=symbols('d'+str(var2))
        self.sdvar='d'+str(var2)
        if show:
            self.s()  
            
 
    def changevar(self,var2,expr,x1='',x2='',relimit=False,show=True):
        '''
        var2=newvar
        expr= self.var in funcionof var2
        '''
        dvar=diffvar(var2) 
        var1=self.var
        newdiff=diff(expr,var2)
        ksym=self.ksym.subs(var1,expr)*newdiff
        if x1=='' and relimit:
            oldx1 = self.x1
            oldx2 = self.x2
        
            # resolver expr = oldx1 y expr = oldx2
            y1 = solve(expr-oldx1, var2)[0]
            y2 = solve(expr-oldx2, var2)[0]
            
            x1 = y1
            x2 = y2
            
        Iksym=creteintegral(ksym,var=var2,x1=x1,x2=x2)
        self.ksym=ksym
        self.Iksym=Iksym
        self.var=var2
        if show:
            self.s()
    def changediff(self,var2,dexpr,x1='',x2='',relimit=False,show=True):
        '''
        changediff(newvar, diff expr...
        '''
        if x1=='':
            x1=self.x1
            x2=self.x2
        expr=self.ksym
        if isinstance(dexpr,MyEq):
            dexpr=dexpr.ksym
            
        E1=createintegral(dexpr,var2)
        expr=E1.doit()
        return self.changevar(var2,expr,x1=x1,x2=x2,relimit=relimit,show=show)        
            
    def changedifferential(self,var2,expr1,x1='',x2=''):
        '''
        var2= new variable 
        expr = if var = x and new var = z and x=3*z*z then  expr1= 3*z*z
        x1= new limit1 when var is z
        x2 = new limit2 when var is z
        '''
        
        var1=self.var
        dvar1=creadiffvar(var1)
        dvar2=creadiffvar(var2)
        expr2=self.getintegrlaexpretion()
        dexpr2=diff(expr1,var2)
        expr3=expr2.subs(var1,expr1)
        kres=dvar2*expr3.subs(dvar1,dexpr2)
        kres=kres.subs(dvar2,1)
        self.ksym=kres
  
        self.var=var2
        if x1!='':
            self.x1=x1
            self.x2=x2
        self.s()

    
    def getdiffvar(var):
        if var==x:
            return dx
        if var==y:
            return dy
        if var==z:
            return dz
        if var==w:
            return dw
        if var==v:
            return dv
        if var==u:
            return du
        if var==t:
            return dt 
            
    def info(self):        
        mdisplay('Iksym:',self.Iksym )
        mdisplay('var:',self.var,',     dvar:',self.dvar,',     sdvar:',self.sdvar)
        mdisplay('ksym:',self.ksym,',      name:',self.name) 
        mdisplay('Limits     ','x1:',self.x1,',      x2:',self.x2)    

        
    def restorevar(self):
        if len(self.changv)>0:
            ksym=self.ksym
            nvar,nexpr=self.changv[-1]
            ovar=self.var
            ksym=ksym.subs(ovar,nexpr)
            self.var=nvar
            self.ksym=ksym
            self.changv=self.changv[0:-1]
            self.rebuild()
        self.s()
 
        
    def changevar2(self,*args):
        if len(args)==1:
            expr1=self.var
            expr2=obj2expr(args[0])
            v1=self.var
            v2=list(symbolslist(expr2))[0]
        elif len(args)==2:
            expr1=args[0]
            expr2=args[1]
            v1=self.var
            v2=list(symbolslist(expr2))[0]
        if len(args)==3:
            expr1=args[0]
            expr2=args[1]
            v1=self.var
            v2=args[2]
         
        kres=rulerchain(self.ksym,v1,v2,expr1,expr2)   
        x1=self.x1
        x2=self.x2
        
         
        if x1!='':
            x1=self.x1
            x2=self.x2
             
            Q=MQ(expr1,expr2,render=False)
 
             
            V=Q.solve(v2,show=False)
            V=obj2expr(V)
            xx1=V.subs(v1,x1,show=False)
            xx2=V.subs(v1,x2,show=False)
             
            self.x1=unisymbols(xx1)
            self.x2=unisymbols(xx2)
        self.var=v2
        self.ksym=unisymbols(kres)
        self.s()
         
    def Area(self,var=''):
        if var=='':
            var=self.var
        f=self.ksym
        vecx=xlimitplot(f,var,self.x1,self.x2) 
        vecp=vecxzonaplot(f,var,self.x1,self.x2)
        At=0
        qq=len(vecx)
        for i in range(qq-1):
            
            At=At+vecp[i]*creteInt(f,self.var,x1=vecx[i],x2=vecx[i+1])
        self.type='P'
        self.ksym=At
        self.s()
        
    def set(self,*args,**kwargs):
        expr=self.ksym
        if len(kwargs)>0:
            expr=real_subs(expr,**kwargs) 
        if len(args)==1 and type(args[0])==MyEq:
            svar=args[0].name
            vvar=args[0].ksym
            expr=expr.subs(svar,vvar)    
        if len(args)==2 and Is_Math(args[0]) and Is_Math(args[1]): 
            expr=expr.replace(args[0],args[1])
     
        self.Iksym = creteintegral(expr, self.var, x1=self.x1, x2=self.x2)
        self.ksym=expr
        self.s()
      
    def positivexpo(self):
        kres=self.ksym
        self.ksym=parse_expr(str(positivexpo(kres)),evaluate=False)
   
        if self.type=='I':
            self.rebuild()        
        self.s()
    
    def setdiff(self,var,expr,x1='',x2=''):
        var1=self.var
        var2=var
        
        dvar1=self.dvar
        dvar2=creadiffvar(var)

        expr1=self.ksym
        try:
            nexpr=expr.subs(dvar2,1)
        except:
            nexpr=expr
            
        if str(dvar2) in str(expr):
            expr2=integrate(nexpr,var)
        else:
            expr2=nexpr
        try:    
            expr3=expr1.subs(var1,expr2)
        except:
            expr3=expr1

        nexppd=diff(expr2,var2)
        expr4=expr3*nexppd
        self.ksym=expr4
        if x1!='':
            self.x1=x1
            self.x2=x2
        self.Iksym=creteintegral(self.ksym,var2,x1=x1,x2=x2)
        self.var=var2
        
        self.s()  
    def tan2sincos(self,ang=alpha):
        kres=self.ksym
        kres2=tan2sincos(kres,ang=ang)
        self.ksym=kres2
        self.rebuild()
        self.s()
  
 
    def diff(self):
        self.type='P'
        return self.ksym
    def setlimits(self,x1,x2):
        self.x1=x1
        self.x2=x2
        self.rebuild()
        self.s()
    def rsimplify(self):
        self.ksym=rsimplify(self.ksym)
        self.s()
        
    def rootsimplify(expr): # simplify each ecponet in expr
        if Is_Add(expr):
            return sum([rootsimplify(data) for data in expr.args])    
        elif Is_Div(expr):
            p1,p2=fraction(expr)
            P1=rootsimplify(p1)
            P2=rootsimplify(p2)
            return sdiv(P1,P2)
        elif Is_Mul(expr):
            return prod([rootsimplify(data) for data in expr.args])
             
        elif Is_Root(expr):
             
            rr=getroot(expr)
            base=insideroot(expr)
            
            bb = getbase(base)
            ee = getexpo(base)
            if not Is_Div(bb):
                try:
                    p1=ee%rr
                    p2=int(ee/rr)
                    if p1==0:
                        kres=bb**p2 
                    else:
                        kres=Mul(bb**p2,rpow(bb**p1,rr))
                    return kres
                except:
                    return expr
            else:
                bb1,ee1=getbase(numer(bb)),getexpo(numer(bb))
                bb2,ee2=getbase(denom(bb)),getexpo(denom(bb))
                if ee1/rr==1 and ee2/rr==1:
                    return spow(sdiv(bb1,bb2),ee)
                elif ee1/rr==1 and ee2/rr!=1:
                    return sdiv(bb1**ee,rpow(spow(bb2,ee2*ee),rr))
                elif ee2/rr==1 and ee1/rr!=1:
                    return sdiv(rpow(spow(bb1,ee1*ee),rr),bb2**ee)
                else:
                    return expr
                    

        elif Is_Pow(expr):
            bb = getbase(expr)
            ee = getexpo(expr)
            return spow(rootsimplify(bb),ee)
        else:
            return expr
        def rootsimplify(self):
            self.ksym=rootsimplify(self.ksym)
            self.s()    
def rulerchain(expr,v1,v2,expr1,expr2):
    ee=MyEq(expr1-expr2,'ee',show=False)
    V1=ee.solve(v1,show=False)
    nv1=V1.ksym
    newd=diff(nv1,v2)
    expr=expr.subs(v1,nv1)
    expr=expr*newd
    return expr 

        

class My2Intg (MyEq):
    '''
        My2Intg(x*y,'P',(x,0,4-4*x/3),(y,0,3))
    
    '''
    
    def __init__(self,expr,name,tuple1,tuple2):
      
        self.ksym=expr
        self.name=name
        self.tuple1=tuple1
        self.tuple2=tuple2
        var1=tuple1[0]
        var2=tuple2[0]
        dvar1=diffvar(var1)
        dvar2=diffvar(var2)
        expr=expr.subs(dvar1,1)
        expr=expr.subs(dvar2,1)
        self.Iksym=Integral(Integral(expr,tuple1),tuple2)
         
        self.type='2I'
        self.s()
    '''    
    def rebuild(self,*args):
        kres= Integral(Integral(self.Ikres,self.tuple1),self.tuple2) 
        self.Iksym=kres 
    ''' 

    def s(self):
        if self.type=='2I':
            expr_tex = latex(self.Iksym)
            name_tex = latex(Symbol(self.name))
            display(Math(f"{name_tex} = {expr_tex}"))
        else:
            expr_tex = latex(self.ksym)
            name_tex = latex(Symbol(self.name))
            display(Math(f"{name_tex} = {expr_tex}"))
        
    def doit(self,kname=''):
        kres=self.Iksym
        kres=kres.doit()
        self.type='P'
        self.ksym=kres
        self.s()
            
    
    def firstIntg(self,name=''):
        kres=self.ksym 
        if name!='':
            ee=MyIntg(kres,name,self.varx, x1=self.x1,x2=self.x2)             
            return ee
        else:
            return self.integral1
    def double2simple(self,expr,name=''):
        if name=='':
            name=self.name+'_1'
        ee=MyIntg(getdataexpr(expr),name,self.vary,x1=self.y1,x2=self.y2)
        return ee
     
    
    
    
    def swapIntegral(self):
        x1,x2,y1,y2=self.x1,self.x2,self.y1,self.y2
        var1=self.varx
        var2=self.vary
        
        try:
            yy2=x1.subs(var2,y2)
        except:
            yy2=0
        try:
            xx2=simplesolve(x1-var1,var2)
            if type(xx2)==list:
                xx2=xx2[-1]
            if xx2==[]:
                xx2=0
        except:
            xx2=0
        
        try:
            yy1=x2.subs(var2,y1)
        except:
            yy1=0
        try:
            xx1=simplesolve(x2-var1,var2)
            if type(xx1)==list:
                xx1=xx1[-1]
            if xx1==[]:
                xx1=0
        except:
            xx1=0
            
        self.x1,self.x2,self.y1,self.y2 =xx1,xx2,yy1,yy2
        self.varx=var2
        self.vary=var1
        self.s()
 
class My3Intg (MyEq):
    '''
        My2Intg(x*y*z,'P',(x,0,4-4*x/3),(y,0,3),(z,0,100))
    
    '''
    def __init__(self,expr,name,tuple1,tuple2,tuple3):
      
        self.ksym=expr
        self.name=name
        self.tuple1=tuple1
        self.tuple2=tuple2
        self.tuple3=tuple3
        kres= Integral(Integral(Integral(expr,tuple1),tuple2),tuple3)
        self.Iksym=kres
        display(Math(name +'='+ latex(kres)))
        self.type='3I'
        
    def rebuild(self):
        kres= Integral(Integral(Integral(self.kres,self.tuple1),self.tuple2),self.tuple3)
        self.Iksym=kres
         
    def s(self):
        
        expr_tex = latex(self.Iksym)
        name_tex = latex(Symbol(self.name))
        display(Math(f"{name_tex} = {expr_tex}"))
         
        
    def value(self):
        kres=self.Iksym.doit()
        return kres
        
    def doit(self,*args,show=True):
        if self.type=='P':
            pass
        else:
            Ires = self.Iksym 
            kres = Ires.doit()

            if 'float' in args:
                try:
                    kres = float(kres)
                except:
                    pass

            if 'positive' in args:
                kres = signo(kres)*kres
                    
            if 'C' in args:
                kres = kres + C1
            if 'C1' in args:
                kres = kres + C1
            if 'C2' in args:
                kres = kres + C2    

            # CONVERTIR OBJETO
            self.__class__ = MyEq
            if len(args)==1 and isinstance(args[0],str):
                self.__init__(kres, args[0], show=False)
            else:
                self.__init__(kres, self.name, show=False)

        if show:
            self.s()   
    
    def firstIntg(self,name=''):
        kres=self.ksym
        if name!='':
            ee=MyIntg(self.ksym,name,self.varx ,x1=self.x1,x2=self.x2)
            return ee
        else:
            return self.integral1
        
    def secondIntg(self,name=''):
        kres=self.ksym
        if name!='':
            ee=My2Intg(self.ksym,name,self.varx,self.vary,x1=self.x1,x2=self.x2,y1=self.y1,y2=self.y2)
            return ee
        else:
            return self.integral2
     
    def triple2double(self,expr,name=''):
        if name=='':
            name=self.name
        ee=My2Intg(getdataexpr(expr),name,self.vary,self.varz,x1=self.y1,x2=self.y2,y1=self.z1,y2=self.z2)
        return ee
    def triple2simple(self,expr,name=''):
        if name=='':
            name=self.name
        ee=MyIntg(getdataexpr(expr),name,self.varz,x1=self.z1,x2=self.z2)
        return ee
    
    
    def swapIntegral(self):
        x1,x2,y1,y2=self.x1,self.x2,self.y1,self.y2
        var1=self.varx
        var2=self.vary
        
        try:
            yy2=x1.subs(var2,y2)
        except:
            yy2=0
        try:
            xx2=simplesolve(x1-var1,var2)
            if type(xx2)==list:
                xx2=xx2[-1]
            if xx2==[]:
                xx2=0
        except:
            xx2=0
        
        try:
            yy1=x2.subs(var2,y1)
        except:
            yy1=0
        try:
            xx1=simplesolve(x2-var1,var2)
            if type(xx1)==list:
                xx1=xx1[-1]
            if xx1==[]:
                xx1=0
        except:
            xx1=0
            
        self.x1,self.x2,self.y1,self.y2 =xx1,xx2,yy1,yy2
        self.varx=var2
        self.vary=var1
        self.s()
    def changevar(self, var_old, expr_new, var_new):
        """
        Cambio de variable en integral definida o indefinida.
        var_old : variable original (x)
        expr_new : expresión en nueva variable (2*u+1)
        var_new : nueva variable (u)
        """

        expr = self.ksym
        x1, x2 = self.x1, self.x2

        # Jacobiano
        J = diff(expr_new, var_new)

        # sustituir integrando
        newexpr = expr.subs(var_old, expr_new) * J

        # intentar invertir transformación
        try:
            inv = solve(var_old - expr_new, var_new)
            if isinstance(inv, list):
                inv = inv[0]
        except:
            inv = None

        # transformar límites si existen
        if x1 is not None and inv is not None:
            newx1 = inv.subs(var_old, x1)
            newx2 = inv.subs(var_old, x2)
        else:
            newx1, newx2 = None, None

        # actualizar integral
        self.ksym = newexpr
        self.varx = var_new
        self.x1 = newx1
        self.x2 = newx2
        self.type = 'I'
        self.s()

def obj2str(obj):
    if type(obj)==str:
        return obj
    elif type(obj)==MyEq:
        return str(obj.ksym)
    elif type(obj)==MyEqEq:
        return str(simplify(expand(obj.L-obj.R)))
    else:
        return str(obj)
def obj2func(obj):
    if type(obj)==str:
        return parse_expr(obj)
    elif type(obj)==MyEq:
        return  obj.ksym 
    elif type(obj)==MyEqEq:
        return  simplify(expand(obj.L-obj.R)) 
    else:
        return obj

def obj2expr(obj):

    # MyEqEq → L-R
    if isinstance(obj, MyEqEq):
        return simplify(expand(obj.L - obj.R))

    # MyEq → ksym
    if isinstance(obj, MyEq):
        return obj.ksym

    # string → parse
    if isinstance(obj, str):
        return parse_expr(obj)

    # sympy expr → directo
    if isinstance(obj, Basic):
        return obj

    # fallback
    return obj        
        
def obj2MyEq(obj,var=x):
    obj=obj2func(obj)
    ee=MyEq(obj,'ee',var=var,show=False)
    return ee
    
def getinterx(obj,var=x):
    ee=obj2MyEq(obj,var=var)
    vx=ee.roots(show=False)
    if 'I' in str(vx) and ee.degree()==3:
        vec3=ee.coef_list()
        vx= list(simplesolve(obj,var))
    vx2=[]
    for i in vx:
        if not 'I' in str(i):
            vx2.append(i)
            
    return vx2

def zonegraph(f,x1,x2):
    expr=obj2expr(f)
    x3=(x2+x1)/2
    ee=obj2MyEq(f)
    kres=float(simplify(ee(x3)))
 
    if kres<0:
        return -1
    elif kres>0:
        return 1
    else:
        return 0

def xlimitplot(obj,var,x1,x2):
    kres=simplesolve(obj,var,'noimg','noshow')
    vecl=[x1,x2]
    for data in kres:
        if data>x1 and data<x2:
            vecl.append(data)
    vecl.sort()
    return vecl 

def vecxzonaplot(obj,var,x1,x2):
    f=obj2expr(obj)
    vecx=xlimitplot(f,var,x1,x2)
    vecp=[]
    qq=len(vecx)
    for i in range(qq-1):
        vecp.append(zonegraph(f,vecx[i],vecx[i+1]))
    return vecp 


 
def MQI(p1,p2):
    P1=p1
    if type(p1)==MyIntg:
        P1=p1.integral1
    P2=p2
    if type(p2)==MyIntg:
        P2=p2.integral1    
    return MQ(P1,P2)

def arclenght(*args,x1=0,x2=1):
    '''
    return Integral to find Arc Lenght fron Fx from x1 to x2
    input:                            output
        arclenght(x)..                   L=Integrate(arcf(x),(x,1,0))
        arclenght(x,'L2')..              L2=Integrate(arcf(x),(x,1,0))
        arclenght(x,'L3',x1=a,x2=b)..    L3=Integrate(arcf(x),(x,a,b))
    '''    
    expr=obj2expr(args[0])
    kname='L'
    var=x
    for i in args[1::]:
        if type(i)==str:
            kname=i
        elif type(i)==symbols:
            var=i
        else:
            pass        
    dx=diff(expr,var)
    dx2=dx*dx+1
    dx3=sqrt(dx2)
    return MyIntg(dx3,kname,var=var,x1=x1,x2=x2)

def surfacerev(*args,x1=0,x2=1):
    '''
    return Integral to find revolution surface fron Fx from x1 to x2
    input:                            output
        surfacerev(x)..                   L=Integrate(surff(x),(x,1,0))
        surfacerev(x,'L2')..              L2=Integrate(surff(x),(x,1,0))
        surfacerev(x,'L3',x1=a,x2=b)..    L3=Integrate(surff(x),(x,a,b))
    '''    
    expr=obj2expr(args[0])
    kname='L'
    var=x
    for i in args[1::]:
        if type(i)==str:
            kname=i
        elif type(i)==symbols:
            var=i
        else:
            pass        
    dx=diff(expr,var)
    dx2=dx*dx+1
    dx3=sqrt(dx2)
    return MyIntg(2*pi*expr*dx3,kname,var=var,x1=x1,x2=x2)


def get_diffvardata(expr):
    '''
    if expr= x*x*dx
    get_diffvardata(expr)  return x,dx
    '''
    vec=list(expr.free_symbols)
    svec=[str(i) for i in vec]
    var=''
    dvar=''
    for i in svec:
        if len(i)>1:
            if i[0]=='d':
                var=unisymbols(parse_expr(i[1::]))
                dvar= creadiffvar(var)
                nexpr=expr.subs(dvar,1)
                 
    return var,dvar,nexpr

def exprdiffparts(expr):
    '''
    if expr= x*x*dx
    get_diffvardata(expr)  return x,dx
    '''
    vec=list(expr.free_symbols)
    svec=[str(i) for i in vec]
    var=''
    dvar=''
    for i in svec:
        if len(i)>1:
            if i[0]=='d':
                var=unisymbols(parse_expr(i[1::]))
                dvar= creadiffvar(var)
                nexpr=expr.subs(dvar,1)
                 
    return var,dvar,nexpr

def get_jacobian(*args):
    '''
    example:
        main variable are x,y
        and cahge var:
            u=x+y
            v=y-exp(x)
        get_jacobian(x,y,x+y,y-exp(x)
        return 1/(expr(x)+1)
    '''        
    vars=[]
    vvec=[]
    svec=['u','v']
    for data in args:
        if type(data)==symbols:
            vars.append(data) 
        else:
            vvec.append(data)
    vece=[MyEq(vdata,sdata,show=False) for vdata,sdata in zip(vvec,svec)]        
    rdata=[]
    for edata in vece:
        for dvar in vars:
            rdata.append(edata.diff(dvar))
    try:        
        M=MyMat(*rdata,2,2,show=False)
        return cfrac(M.det())
    except:
        print('Please load My_Matrix library before')