
from sympy import *
from lib_MyEqEq import MyEqEq,MQ
from lib_MyEq import *
from libaldo_math2 import *
import copy
 







from sympy import *
 
from IPython.display import Math  # ,display
from matplotlib.pyplot import ylabel, plot, show, xlabel, title
from libaldo_math2 import *
from libaldo_algorith import *
 
from lib_toolsEq import *
from lib_tools import * 
from mathbasic import *
 
import copy
import numpy as np
import pandas as pd
 
import matplotlib.pyplot as plt
# from lib_Func import *
import dill  # pip install dill --user
import pickle
filename = 'workspace.pkl'
x1=symbols('x1')
 
def savework():
    dill.dump_session(filename='workspace.pkl')


def loadwork():
    dill.load_session(filename='workspace.pkl')


# and to load the session again:
def creteInt(expr,var,x1='',x2=''):
    if x1=='': 
        kres=Integral(expr,var)
    else:
        kres=Integral(expr, (var,x1,x2))
    return kres  

#direx,direy=symbols('\overrightarrow{x} \overrightarrow{y}')
class MyIntg (MyEq):
    def __init__(self,expr,name,var=x,x1='',x2='',kshow=True,doit=False):
        self.type='I'
        self.ksym=obj2expr(expr) 
        self.name=name
        self.x1=getdata(x1)
        self.x2=getdata(x2)
        self.var=var 
        self.integral1=creteInt(expr,var,x1=x1,x2=x2)
        self.varc=False
        self.doit=doit
         
        if kshow:
            display(Math(self.name +'='+ latex(self.integral1)))
    def __call__(self,doit=False):
        if doit:
            kres=self.doit()
        else:
            kres=self.integral1
            
        return kres
       

    def __repr__(self):
        kres = str(self.ksym)

        return kres

    def _latex(self, obj):
        return latex(self.ksym)

    def __str__(self):
        return self.__repr__()    
    def rebuild(self):
        self.integral1=creteInt(self.ksym,var=self.var,x1=self.x1,x2=self.x2)
        
    def s(self):
        if self.type=='I':
            self.rebuild()
            display(Math(self.name +'='+ latex(self.integral1)))
        else:
            display(Math(self.name +'='+ latex(self.ksym)))
    def doit(self,*args):
        if self.type=='P':
            kres=self.ksym
            kres=kres.doit()
             
        else:
            Ires=self.integral1
            kres=Ires.doit()
        if 'C' in args:
            print(1)
            kres=kres+C1
        else:    
            if 'float' in args:
                try:
                    kres=float(kres)
                except:
                    pass
            if 'positive' in args:
                kres=signo(kres)*kres
        if 'up' in args or 'update' in args:
            self.ksym=kres
            self.type='P'
            self.s()
            
        else:
            return kres
     
    def insideIntg(self):
        return self.ksym
    
            
       
    
 
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
            
    def changevar(self,*args):
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
 
             
            V=Q.solve(v2,kshow=False)
            V=obj2expr(V)
            xx1=V.subs(v1,x1,kshow=False)
            xx2=V.subs(v1,x2,kshow=False)
             
            self.x1=unisymbols(xx1)
            self.x2=unisymbols(xx2)
        self.var=v2
        self.ksym=unisymbols(kres)
        self.s()
         
    def Area(self):
        f=self.ksym
        vecx=xlimitplot(f,self.x1,self.x2) 
        vecp=vecxzonaplot(f,self.x1,self.x2)
        At=0
        qq=len(vecx)
        for i in range(qq-1):
            
            At=At+vecp[i]*creteInt(f,self.var,x1=vecx[i],x2=vecx[i+1])
        self.type='P'
        self.ksym=At
        self.s()
    
    
    
def rulerchain(expr,v1,v2,expr1,expr2):
    ee=MyEq(expr1-expr2,'ee',kshow=False)
    V1=ee.solve(v1,kshow=False)
    nv1=V1.ksym
    newd=diff(nv1,v2)
    expr=expr.subs(v1,nv1)
    expr=expr*newd
    return expr            

class My2Intg (MyEq):
    def __init__(self,expr,name,varx,vary,x1='',x2='',y1='',y2='',kshow=True):
        self.type='II'
        self.ksym=expr
        self.name=name
        self.x1=x1
        self.x2=x2
        self.varx=varx
        self.integral1=creteInt(expr,varx,x1=x1,x2=x2)
       
        self.vary=vary
        self.y1=y1
        self.y2=y2
        self.integral2=creteInt(self.integral1,vary,x1=y1,x2=y2)
         
        if kshow:
            display(Math(name +'='+ latex(self.integral2)))
        
    def rebuild(self):
        self.integral1=creteInt(self.ksym,var=self.varx,x1=self.x1,x2=self.x2)
        self.integral2=creteInt(self.integral1,var=self.vary,x1=self.y1,x2=self.y2)
    def s(self):
        self.rebuild()
        display(Math(self.name +'='+ latex(self.integral2)))
    def doit(self,kname=''):
        kres=self.integral2
        kres=kres.doit()
        if kname!='':
            ee=MyEq(kres,kname)
            return ee
        else:
            return kres    
    
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
    def __init__(self,expr,name,varx,vary,varz,x1='',x2='',y1='',y2='',z1='',z2='',kshow=True):
        self.ksym=expr
        self.name=name
        self.x1=x1
        self.x2=x2
        self.varx=varx
        self.integral1=creteInt(expr,varx,x1=x1,x2=x2)
       
        self.vary=vary
        self.y1=y1
        self.y2=y2
        self.integral2=creteInt(self.integral1,vary,x1=y1,x2=y2)
        
        self.varz=varz
        self.z1=z1
        self.z2=z2
        self.integral3=creteInt(self.integral2,varz,x1=z1,x2=z2)
         
        if kshow:
            display(Math(name +'='+ latex(self.integral3)))
        
    def rebuild(self):
        self.integral1=creteInt(self.ksym,var=self.varx,x1=self.x1,x2=self.x2)
        self.integral2=creteInt(self.integral1,var=self.vary,x1=self.y1,x2=self.y2)
        self.integral3=creteInt(self.integral2,var=self.varz,x1=self.z1,x2=self.z2)
    def s(self):
        self.rebuild()
        display(Math(self.name +'='+ latex(self.integral3)))
    def doit(self,kname=''):
        kres=self.integral3
        kres=kres.doit()
        if kname!='':
            ee=MyEq(kres,kname)
            return ee
        else:
            return kres    
    
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
    if type(obj)==str:
        return parse_expr(obj)
    elif type(obj)==MyEq:
        return  obj.ksym 
    elif type(obj)==MyEqEq:
        return  simplify(expand(obj.L-obj.R)) 
    else:
        return obj        
        
def obj2MyEq(obj,var=x):
    obj=obj2func(obj)
    ee=MyEq(obj,'ee',var=var,kshow=False)
    return ee
    
def getinterx(obj,var=x):
    ee=obj2MyEq(obj,var=var)
    vx=ee.roots(kshow=False)
    if 'I' in str(vx) and ee.degree()==3:
        vec3=ee.coef_list()
        vx= list(solve3(*vec3))
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

def xlimitplot(obj,x1,x2):
    obj=obj2expr(obj)
    vecx=getinterx(obj)
    vecp=[x1,x2]
    for i in vecx:
        if i>x1 and i<x2:
            vecp.append(i)
    vecp.sort()
    return vecp 

def vecxzonaplot(obj,x1,x2):
    f=obj2expr(obj)
    vecx=xlimitplot(f,x1,x2)
    vecp=[]
    qq=len(vecx)
    for i in range(qq-1):
        vecp.append(zonegraph(f,vecx[i],vecx[i+1]))
    return vecp 


 
        