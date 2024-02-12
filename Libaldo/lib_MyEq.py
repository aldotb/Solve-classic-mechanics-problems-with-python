

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
from lib_funcdiff import *
x1=symbols('x1')

 
nombresdiferen=['dx','dy','dz','du','dv','dw','d2x','d2y','d2z','d2u','d2v','d2w'] 

def savework():
    dill.dump_session(filename='workspace.pkl')


def loadwork():
    dill.load_session(filename='workspace.pkl')


# and to load the session again:

25
C1, C2, C3, C4, t, x, y, z,t1,t2 = symbols('C1 C2 C3 C4 t x y z t1 t2')
dataQ = []
e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12 = symbols('e1 e2 e3 e4 e5  e6 e7 e8 e9 e10 e11 e12')
alpha=symbols('alpha')
#direx,direy=symbols('\overrightarrow{x} \overrightarrow{y}')
class MyEq:
    def __init__(self, ksym, kname='', var=x,var1=y,var2=t,  varx=x, vary=y, varz=z, kp=False, kope='', kshow=True, ktype='P',
                 xx='', dtype=1, depen=False, x1='', x2='', y1='', y2=0, z1='', z2=0, Pobj='', ee='', ssym='', init=True,
                 kfull=True, andsolve='',varf=[],diffEq='',eQshow='',vfunc=[],vmain=x,mark=False,link=''):
        self.mark=mark
        if type(ksym)==str:
            self.ksym=traducediff(ksym,var,var1)
             
        else:    
            self.ksym=ksym
        self.t=ktype
        self.pname=kname
        self.type=ktype
        self.kinte=''
        self.ee=ee
        self.depen=depen
        self.primi=ksym
        self.var1=var1
        self.var2=var
        self.var=var        
        self.link=link 
         
        self.x1=x1
        self.x2=x2
                
        self.varL=''        
          
         
        self.y1=y1
        self.y2=y2
        self.z1=z1
        self.z2=z2
        self.varx=varx
        self.vary=vary
        self.varz=varz
        self.eeq=ksym
        self.vmain=vmain
        self.diffvar=''   
        if kname!='' and len(str(kname))==1:
            self.vmain=symbols(kname)
            
         
         
         
         
        self.name = alphaname(kname)
        self.sname=''
        self.kinte=''
        self.xx=xx
        self.EqDiff=''
        self.Adiff=''
        self.Bdiff=''
        self.EqInte=''
        self.iniG=0
        self.odeQ=''
        self.ode=''
        self.oldprimi=''
        self.Pobj=Pobj
        self.ssym=ssym
        self.father=''
        self.backup=[] 
        self.kshow=kshow
        self.init=init 
        self.histo=ksym
        self.primitiva=func2primi(ksym,var2,varf)
        self.tipoView=True
        try:
            self.pdiff=diff(func2primi(ksym,var2,varf),var2)
        except:
            self.pdiff=0
        try:
            self.flatdiff=Diff2diff(diff(func2primi(ksym,var2,varf),var2),var2,varf)
        except:
            self.flatdiff=0
        
        self.varf=varf    
        self.primi_eq=ksym
        self.diff_eq=''
        self.primi_diff_eq='' 
        self.origen=1
        self.var1=var1
        self.var2=var2
        self.diffname='d'+str(var)
        if type(ksym)==MyEq:
            seq=str(ksym())
            self.ksym = unisymbols(opemat(ksym(), kope=kope))
            self.v = unisymbols(opemat(ksym(), kope=kope))
            self.eeq=ksym
        #elif type(ksym)==MyEqMat:
        #    self.ksym = unisymbols(opemat(ksym(), kope=kope))
        #elif type(ksym)==Equality:
        #    mm=ksym.e2.ksym
        #    self.ksym=mm=kQ.e2.ksym
        #     
        #    try:
        #        self.var2=get_varfunc(ksym.lhs)
        #    except:
        #        done=False
        else:    
            seq=str(ksym)
            #self.ksym = unisymbols( ksym )
            self.v = unisymbols( ksym )

        if ktype=='F' or ktype=='Ff':
            if var2!='':
                s1=kname
                if type(var2)==list:
                    s2=alphaname(var2[0])+','+alphaname(var2[1])
                else:
                    s2=alphaname(var2)
                kname=s1+'_{('+s2+')}'
                self.name = kname
         
        
        
        
            
        if 'I' in ktype:
            ksym=self.ksym
            if Is_Number(ksym):
                self.ksym=ksym*var+C1
            else:    
                ddx='d' + str(self.var)
                ssym=str(ksym)
                
                if ddx in ktype:
                    ssym=str(self.ksym) 
                    ssym=ssym.replace(ddx,'1')
                    ksym=parse_expr(ssym)
                    
                    self.ksym=ksym
                    self.primi=self.ksym
                if self.x1=='':
                    ksym=self.ksym
                    if Is_Number(ksym):
                        self.ksym=ksym*var+C1
                    else:    
                        self.ksym=Integral(ksym,var)
                else:
                    self.ksym=Integral(self.ksym, (var,self.x1,self.x2))
                self.histo = unisymbols(opemat(self.ksym, kope=kope))
        if '2I' in ktype or '3I' in ktype:
             
            if self.y1=='': 
                self.ksym=Integral(self.ksym, self.vary)
            else:
                self.ksym=Integral(self.ksym, (self.vary,self.y1,self.y2))
            self.histo = unisymbols(opemat(self.ksym, kope=kope))

        if '3I' in ktype:
             
            if self.z1=='': 
                self.ksym=Integral(self.ksym, self.varz)
            else:
                self.ksym=Integral(self.ksym, (self.varz,self.z1,self.z2))
            self.histo = unisymbols(opemat(self.ksym, kope=kope))            
        
        if ktype=='Diff':
        
            svard='d'+str(var)
            sres=str(ksym)
            if svard in sres:
                sres=sres.replace(svard,'1')
                self.ksym=parse_expr(sres)
            self.Adiff=Derivative(self.ksym, var)
            self.sname='d_'+kname 
            self.name=kname
            self.var2=var 
                
         
            
        if ktype=='diff' or ktype=='diff2':
            f=Function(kname)(var2)
            kres=self.ksym
            kres=kres.subs(var1,f)
             
            self.Bdiff=kres
                
            if ktype=='diff':
                self.Adiff=Derivative(f,var1)
                 
            else:
                self.Adiff=Derivative(f, var1,var1)
            self.EqDiff=Eq(self.Adiff,self.Bdiff)
                

        if ktype=='fdiff':
            self.name='d'+kname+'_'+str(var2)
            nksym=get_diff_name(str(ksym),str(var2))
             
            
            self.ksym=nksym
            self.var1=ksym
            
        if self.Pobj!='':
            self.type='Ph'
            
        if kshow:
            if ktype=='diff' or ktype=='diff2':
                display(Math(latex(self.EqDiff)))
            elif ktype=='D':
                sres=str(ksym)
                pdd=['d'+str(self.var),'d_'+str(self.var), 'd'+alphaname(self.var),'d_'+alphaname(self.var)]
                for i in pdd:
                    if i in sres:
                        sres=sres.replace(i,'1')
                fres=parse_expr(sres)
                self.ksym=fres    
 
                sR = self.name + ' ='
                kres=self.ksym
                self.diffname=diff_name(self.var)
                display(Math(sR + latex(kres)+' '+self.diffname))
                  
                
            elif self.mark:
                ksym=self.ksym
                dexpr=diff2mark2(ksym,var,var1)
                ee=MyEq(dexpr,kname=self.kname)
                
                
            else:
                kres=self.ksym 
            
                if self.name == '':
                    display(Math(latex(kres)))
                else:
                    sR = self.name + ' ='
                    display(Math(sR + latex(kres)))
        
        
        if andsolve!='':
            kval=andsolve
            kname=str(andsolve)
             
            kres=csolve(self.ksym,kval)
            kres=opemat(kres,kope)
            andsolve=MyEq(kres,kname)
            return (andsolve) 
        if dtype==1:
            if self not in dataQ and self.name!='': 
                dataQ.append(self)
        #self.ksym=sympify(str(self.ksym), locals={str(var2): var2})     
  
     

        self.primi_eq=ksym 
        self.diff_eq=''
        self.primi_diff_eq=''
        self.origen=1
        self.eQshow=eQshow
        self.vfunc=vfunc
        self.vdf=[]
        self.vd2f=[]
        self.vdp=[]
        self.vd2p=[]
        if "'" in str(self.ksym):
            self.type='DI'
        
     
    def __call__(self,*args,kshow=True, **kwargs):
        op=[]
        nargs=[]
        vecargs=['value','expand','eval']
        kname=''
        if len(args)==0 and len(kwargs)==0:
            return self.ksym
        if len(args)==1 and len(kwargs)==0: 
            ksym=self.ksym
            var=self.var
            valor=args[0]
            kres=ksym.subs(var,valor)
            
            return kres
        elif len(kwargs)>0 and len(nargs)==0:
            kres=self.ksym
            kres=real_subs(kres,**kwargs)
            return kres
        else:
            if len(kwargs)>0 and len(nargs)==0:
                kres=self.ksym
                kres=real_subs(kres,**kwargs)
                return kres
                
                
            for i in args:
                if type(i)==str:
                    if i not in vecargs:
                        kname=i 
                    else:    
                        op.append(i)
                else:
                    nargs.append[i]

             
            if len(kwargs)==0 and len(args)==1:
                kvalue=args[0]
                expr=self.ksym
                expr=expr.subs(str(self.var),kvalue)

                return expr


            kres=self.ksym
            if len(kwargs)==0 and len(args)==0:
                return self.ksym
            if len(kwargs)>0:
                if self.type=='dP':
                    kres= self.primi_diff_eq
                    kres=real_subs(kres,**kwargs) 
                    return kres
                kres=real_subs(kres,**kwargs)

            if len(args)==1 and len(kwargs)==0 and Is_Number(args[0]):
                kres=self.ksym
                var=self.var
                val=args[0]
                return kres.subs(var,val)
            if len(args)>0:
                if 'float' in args:
                    try:
                        kres=float(kres)
                    except:
                        pass
                if 'value' in args:
                    return kres

            nname=''
            if len(args)==1 and len(kwargs)==0 and type(args[0])!=str:
                kres=kres.subs(self.var1,args[0])



            if 'update' in args:
                self.ksym=kres
                self.s()
                return
            for i in args:
                if i!='float' and i!='update' and type(i)==str:
                    nname=i
            if nname=='':
                ee=self.xcopy(self.name,kshow=False)
                ee.ksym=kres
                ee.s()
                return
            else:
                ee=MyEq(kres,kname=nname,kshow=kshow)
                return ee
       

    def __repr__(self):
        kres = str(self.ksym)

        return kres

    def _latex(self, obj):
        return latex(self.ksym)

    def __str__(self):
        return self.__repr__()

    ###########################################
    #               variables                 #
    ###########################################

    def norm_variable(self):
        clist = [C1, C2]
        m1 = self.free()
        s1 = [str(x) for x in m1]
        for kvar in clist:
            svar = str(kvar)
            for i in range(len(s1)):
                if s1[i] == svar:
                    kvar = m1[i]
    def validatesymbols(self,svar):
        self.ksym=expr2var(self.ksym,svar)
         
        self.s()

    ###########################################
    #               Update                    #
    ###########################################

    def __add__(self, other):
        """ Returns the vector addition of self and other """
        if type(other) == MyEq :
            kres = self.ksym + other.ksym
        else:
            kres = self.ksym + other
        if type(kres) == MyEq:
            kres.s()
        return kres

    def __radd__(self, other):
        if type(other) == MyEq :
            kres = self.ksym + other.ksym
        else:
            kres = self.ksym + other
        return kres

    def __sub__(self, other):
        if type(other) == MyEq :
            kres = self.ksym - other.ksym
        else:
            kres = self.ksym - other
        return kres

    def __rsub__(self, other):
        if type(other) == MyEq :
            kres = self.ksym - other.ksym
        else:
            kres = self.ksym - other
        return kres

    def __mul__(self, other):
        """ Returns the vector addition of self and other """
        if type(other) == MyEq :
            kres = self.ksym * other.ksym
        else:
            kres = self.ksym * other
        return kres

    def __rmul__(self, other):
        """ Returns the vector addition of self and other """
        if type(other) == MyEq :
            kres = self.ksym * other.ksym
        else:
            kres = self.ksym * other
        return kres
        
    
    def __truediv__(self, other):
        if type(other)==MyEq:
            kres=other.ksym
            return  self.ksym / kres
        else:    
            return  self.ksym / (1*other) 

    def __rtruediv__(self, other):
        if type(other)==MyEq:
            kres=other.ksym
            return  kres /self.ksym
        else:    
            return  (1*other)/self.ksym 
    
        
    def Di(self):
        if self.varf!=[]:
            return self.primi_eq
            
    ######################################## 
    #       show                           #
    ########################################   
    
    def s2(self,op=''):

        kres = self.ksym
        self.ksym=unisymbols(kres)
        if op=='2':
            kres=self.eQshow
        if self.type == 'Diff':

            ps1 = self.sname + '='
            ps2 = diff_name(self.var)

            display(Math(ps1 + latex(kres) + ps2))

        
        elif self.type == 'D':
            sR = self.name + ' ='
            kres=self.ksym
            display(Math(sR + latex(kres)+' '+self.diffname))
            
        elif self.type == 'I':
            sR = self.name + ' ='
            kres=self.ksym
            display(Math(sR + latex(kres)))
            
            

        elif self.type == 'diff' or self.type == 'diff2':
            display(Math(latex(self.EqDiff)))
            
        elif self.type=='LA':
            try:
                kres=self(L=latex2sympy('\\mathcal{L}'))
            except:
                kres=kres.subs('L',latex2sympy('\\mathcal{L}'))
            sR = self.name + ' ='
            display(Math(sR + latex(kres)))
            
        else:
            if self.name == '':
                display(Math(latex(kres)))
            else:
                sR = self.name + ' ='
                display(Math(sR + latex(kres)))
    def update(self,kres):
        self.ksym=kres
        
    def ss(self):
        sR = self.name + ' ='
        display(Math(sR + latex(self.ksym)))
        
    def s(self,op='', **kwargs):
        if type(self.ksym)==UnevaluatedExpr:
            self.ksym=unisymbols(self.ksym)
     
        if self.tipoView:
            self.ss()
        else:    
            kres = self.ksym
            self.ksym=unisymbols(kres)
            if self.type==('dP'):
                kres=self.primi_diff_eq 
                ps1 = self.name + '='
                ps2 = kres

                display(Math(ps1 + latex(ps2)))
            else:
                #hvar = pickle.load(open("trial.p","rb")) 
                if self.father=='':
                    qk = len(kwargs)
                    if qk > 0:
                        ee = MyEq(self.v, kshow=False)
                        for key, value in kwargs.items():
                            ee.set(parse_expr(key), value, kshow=False, ktype=self.type)
                        ee.s2(op=op)
                         
                        #hvar.append([type(ee),type(ee.ksym),ee.ksym])
                         
                    else:
                        self.s2(op=op)
                        #hvar.append([type(self),type(self.ksym),self.ksym])
                    
                    #pickle.dump(hvar, open("trial.p", "wb"))
         
             
                else:
                    mQ.s()
            if self.link!='':
                self.link.s()
    def showd(self):
        kres=self.ksym
        for i in self.vfunc:
            kres=show_modefunc(kres,i,ktype=2)
        ee=MyEq(kres,kname=self.name)
        
    def showp(self):
        kres=self.ksym
        for i in self.vfunc:
            kres=show_modefunc(kres,i,ktype=3)
        ee=MyEq(kres,kname=self.name)
      

    def unisymbols(self):
        kres=self.ksym
        kres=unisymbols(kres)
        self.ksym=kres
    ######################################## 
    #       set & update
    ########################################    
    def MsetValue(self, vknom, vkval, kope=''):
        for knom, kval in zip(vknom, vkval):
            knom = unisymbols(knom)
            kval = unisymbols(kval)
            kres = self.ksym
            kres = kres.subs(knom, kval)
            kres = opemat(kres, kope=kope)

            kres = self.simplify()
            #self.ksym=kres
        self.s()

    def multiSet(self, vknom, vkval, kope=''):
        for knom, kval in zip(vknom, vkval):
            knom = unisymbols(knom)
            kval = unisymbols(kval)
            kres = self.ksym
            kres = kres.subs(knom, kval)
            kres = opemat(kres, kope=kope)

            kres = self.simplify()
            self.ksym=kres
        self.s()




    def sset(self,kname='',vmain='',var2='',kshow=True,**kwargs):
        if vmain=='':
            vmain=self.vmain
        if var2=='':
            var2=self.var2
        if len(kwargs)>0:
             
            kres=self.ksym
            for key, value in kwargs.items():
                nsym=parse_expr(key)
                if nsym not in fpoly(kres,'free'):
                    kname=str(nsym)
                    var2=self.var2
                    kres=kres.subs(Function(kname)(var2),value)
                    
                else:
                    kres=kres.subs(parse_expr(key),value)
                                   
        if kname!='':
             
            nQQ=MyEq(kres,kname=kname,vmain=vmain,var2=var,kshow=kshow)
            return nQQ
        else:    
            self.ksym=kres    
        if kshow:    
            self.s()    
    
    def evalif(self, kname='', **kwargs):
        ee = self.xcopy(kname=kname, kshow=False)

        qk = len(kwargs)
        if qk > 0:
            for key, value in kwargs.items():
                if type(value) == MyEq:
                    value = value.ksym
                ee.set(parse_expr(key), value, kshow=False)
                try:
                    ee.primi = ee.primi.subs(parse_expr(key), value)
                except:
                    done = False

            if kname != '':
                ee.s()
                return ee
            else:
                ee.s()

        else:
            self.s()

    def set(self,*args, kshow=True,**kwargs):
        ksym=self.ksym
        if len(kwargs) > 0:
            ksym=realsub2(ksym,**kwargs)
        if len(args)==2:
            if type(args[0])!=MyEq:
                p1=args[0]
                p2=getdata(args[1])
                try:
                    ksym=ksym.subs(p1,p2)
                except:
                    ksym=subsubs(ksym,p1,p2)
                    
                 
        if len(args)>0 and type(args[0])==MyEq:
            for i in args:
                p1=kname(i)
                p2=i.ksym
                try:
                    ksym=ksym.subs(p1,p2)
                except:
                    ksym=subsubs(ksym,p1,p2)
                            
                
        self.ksym=ksym     
        if kshow:     
            self.s()
        

    def setValue(self, knom, kval, kshow=True, kope='', kret=False):
        if type(knom) != list:
            knom = [knom]
            kval = [kval]
        for i, j in zip(knom, kval):
            i = unisymbols(i)
            j = unisymbols(j)
            kres = self.ksym
            try:
                kres = kres.subs(i, j)
            except:
                pass
            kres = opemat(kres, kope)

            self.ksym=kres

        if kret:
            return self.ksym

        if kshow:
            self.s()
 
    def setdiff(self, kvar, kval, kshow=True, kope='', kret=False):
        try :
            kval=alphaname(kval)
        except:
            done=False
        self.set(knom=kvar.diff(), kval=kval, kshow=kshow, kope=kope, kret=kret)


    def eval(self,*args, **kwargs):
         
        ksym = self.ksym
        ksym=real_subs(ksym,**kwargs)
        if 'float' in args:
            ksym=float(ksym)
        return ksym
        
 
    def set_solve(self, kset, kvset, kvars, knames):
        ee = MyEq(self.ksym, kshow=False)
        ee.set(kset, kvset, kshow=False)
        kres = ee.solve(kvars)
        ee2 = MyEq(kres, knames)
        return ee2

    def set_solveR(self, kset, kvset, kvars, knames):
        ee = MyEq(self.ksym, kshow=False)
        ee.set(kset, kvset, kshow=False)
        kres = ee.solveR(kvars)
        ee2 = MyEq(kres, knames)
        return ee2



    def v(self):
        return self.v
    def value(self):
        return self.ksym
        
    def update(self, *args):
        kres=self.ksym
        for i in args:
            val=i.ksym
            name=i.name
            kres=kres.subs(name,val)
        self.ksym=kres
        self.s()
            
         

    def upgradeinteger(self,*args,**kwargs):
        kres=self.ksym
        p1,p2=kres.args
        ee=MyEq(p1,'ee',kshow=False)
        for i in args:
            ee.upgrade(*args,kshow=false)
        eres=ee.ksym
        eres=real_subs(eres,**kwargs)
        return Integral(eres, p2)
          
    def upgradeintegerdiff(self,newdiff,x1='',x2=''):
        
        newv=newdiff.var
        kres=self.ksym
        p1,p2=kres.args
        dsym=newdiff.ksym
        p1=p1*dsym
        p2=(newv,p2[1],p2[2]) 
        if x1!='':
            p2=(newv,x1,x2)

        kres= Integral(p1,p2)
        self.var=newv
        self.diffname=newdiff.diffname
        self.ksym=kres
            
            

    def up(self,*args,x1='',x2=''):
        for i in args:
            self.upgrade(i,x1=x1,x2=x2,kshow=False)
        self.s()
        
    def change_integral(self,newval='',nvar='',x1='',x2=''):
        kres=self.ksym
        p1,p2=kres.args
        if newval!='':
            p1=newval
        var=self.var
        if nvar!='':
            var=nvar
        xx1=self.x1
        xx2=self.x2
        if x1!='':
            xx1=x1
            xx2=x2
        p2=(var,xx1,xx2)
        self.primi=p1
        self.ksym= Integral(p1,p2)
    def get_intedata(self):
        kres=self.ksym
        p1,p2=kres.args
        kres=p1
        var,x1,x2=p2
        return(kres,var,x1,x2)
    
    ##############################################
    ##  Transformacion
    
    def replaceexpr(self,val1,val2):
        ssym=str(self.ksym)
        ssym=ssym.replace(str(val1),str(val2))
        try:
            self.ksym=parse_expr(ssym)
        except:
            pass
        self.s()

        
    def replace(self,*args,kshow=True,**kwargs): #@tif
        p1=self.ksym
        p1=real_subs(p1,**kwargs)
        
        if len(args)>0:
            for i in args:
                
                if type(i)==MyEq:
                    sname=i.name
                    value=i.ksym
                else:
                    sname=str(i)
                    value=i
                p1=p1.subs(sname,value)
                 
            self.ksym=p1
             
        if kshow:
            self.s()

            

    def upgrade(self, *args,x1='',x2='',kope='', kshow=True, **kwargs):
           
        if self.type=='I':
             
            kdiffe=args[0]
            ktypo=kdiffe.type
            if ktypo=='D':
                self.upgradeintegerdiff(args[0],x1=x1,x2=x2)
            if ktypo=='P':
                kres,var,x1,x2=self.get_intedata()
                ee=MyEq(kres,'ee',kshow=False)
                ee.upgrade(kdiffe,kshow=False)
                kres=ee.ksym
                self.change_integral(newval=kres)
              
            
        else:
            if len(args) == 1:
                if type(args[0]) == list:
                    args = args[0]

            for i in args:
                kname = i.name
                if len(kname)==2:
                    sres=str(self.ksym)
                    oname=kname[0]+'_'+kname[1]
                    if oname in sres:
                        self.set(oname, i, kshow=False, kope=kope)
                    else:
                        self.set(kname, i, kshow=False, kope=kope)
                else:
                    self.set(kname, i, kshow=False, kope=kope)
                if 'I' in self.type:
                    ee = MyEq(self.primi, kshow=False)
                    ee.set(kname, i, kshow=False, kope=kope)
                    self.primi = ee.ksym
                try:
                    self.primi = self.primi.subs(kname, i)
                except:
                    done = False
            if len(kwargs) > 0:
                self.set(**kwargs, kshow=False)
        if kshow:
            self.s()

    def kret(self):
        return self.ksym
        sR = self.name + ' ='
        
    ###########################################
    #            edit copy
    ###########################################

    def xcopy(self, kname, kshow=True):
        ee = copy.deepcopy(self)
        if kname != '':
            ee.name = kname
        if kshow:
            ee.s()
        return ee
    def sExp(self, kval):
        if self.name == '':
            display(Math(latex(self.ksym)))
        else:
            sR = self.name + ' ='
            display(Math(sR + latex(kval)))

    def undo(self,kshow=True):
        self.ksym = self.histo
        self.v = self.histo
        if kshow:
            self.s()

    def cut_denom(self, kope=''):
            kres = numer(self.ksym)
            kres = opemat(kres, kope)
            self.ksym=kres
            return self.ksym


    

    def reduF(self, kshow=True, kupdate=True):

        kres = self.ksym
        kres = factor(kres)
        kres2 = 1
        if Is_Mono(kres):
            kres = numer(kres)

            mm = fpoly(kres, 'list')
            for i in mm:
                if Is_Poly(i):
                    kres2 = kres2 * i

        if kres2 != 0 and kres2 != 1:
            self.ksym = kres2
        self.s()
    ###########################################
    #            Math operation
    #  Add, SUbs, Mul, Div, Pow, Rpow
    ###########################################
    
    ##   Add    
    def add(self,*args,kshow=True,**kwargs): #@tif
        return self.Add(*args,kshow=True,**kwargs)
        
    def Add(self,avar,kname='', kshow=True,**kwargs): #@tif
        kres=self.ksym
        kres=real_subs(kres,**kwargs)
        kres=kres+avar
        if kname!='':
             
            ee=MyEq(kres,kname,var=self.var,kshow=kshow)
            return ee
        else:
            self.ksym=kres
            if kshow:
                self.s()
    def Substrac(self,avar,kname='', kshow=True,**kwargs): #@tif
        kres=self.ksym
        kres=real_subs(kres,**kwargs)
        kres=kres-avar
        if kname!='':
             
            ee=MyEq(kres,kname,var=self.var,kshow=kshow)
            return ee
        else:
            self.ksym=kres
            if kshow:
                self.s()
    def subsubs(self,*args,force=False,text=False):
        expr=self.ksym
        if type(args[0])==list:
            vec1=args[0]
            vec2=args[1]
            for i,j in zip(vec1,vec2):
                expr=subsubs(expr,i,j,force=force)    
        else:
            expr=subsubs(expr,args[0],args[1],force=force)
         
        self.ksym=expr
        self.s()
        
    def subs(self,val1,val2,kshow=True):
        if type(val1)==MyEq:
            val1=val1.ksym
        if type(val2)==MyEq:
            val2=val2.ksym   
        kres=self.ksym
        kres=kres.subs(val1,val2)
        self.ksym=kres
        if kshow:
            self.s()
        
    
    def Subs(self,*args,kshow=True,**kwargs): #@tif
        kres=self.ksym
        kres=real_subs(kres,**kwargs)
        kname=''
        
        for i in args:
            if type(i)==str:
                kname=i 
            else:
                kres=kres-gval(i)
        if kname=='':
            self.ksym=kres
            if kshow:
                self.s()
        else:
            ee=MyEq(kres,kname=kname,kshow=kshow)
            return ee        
    ##   Mul   
    def mul(self,*args,kshow=True,**kwargs): #@tif
        return self.Mul(*args,kshow=True,**kwargs)    
    def Mul(self,*args,ope='',kshow=True,**kwargs): #@tif
        kres=self.ksym
        kres=real_subs(kres,**kwargs)
        kname=''
        
        for i in args:
            if type(i)==str:
                kname=i 
            else:
                kres=kres*gval(i)
        if ope!='':
            kres=opemat(kres,ope)
            
        if kname=='':
            self.ksym=kres
            if kshow:
                self.s()
        else:
            ee=MyEq(kres,kname=kname,kshow=kshow)
            return ee    

    

    ##   Mul   
    def div(self,dval,kshow=True):  
        return self.Div(dval=dval,kname=kname,kshow=kshow)
        
    def Div(self,*args ,kshow=True): 
        vecop=['simplify','factor','expand']
        dval=args[0]
        kres=self.ksym
        kres=Div(kres,dval)
        if 'simplify' in args:
            kres=simplify(kres)
        if 'factor' in args:
            kres=factor(kres)    
        if 'expand' in args:
            kres=expand(kres)
        self.ksym=kres
        if kshow:
            self.s()
        
        

    def Pow(self,kpo,kname='',ope='',kshow=True,**kwargs): #@tif
        kres=self.ksym
        kres=real_subs(kres,**kwargs)
        kres=kpow(kres,gval(kpo))
        if ope!='':
            kres=opemat(kres,ope)
        if kname!='':
            ee=MyEq(kres,kname=kname,kshow=kshow)
            return ee
        else:
            self.ksym=kres
            if kshow:
                self.s()
    
    def sqrs(self,*args,ope='',kupdate=True, kshow=True): # #@tif  return pow(k1,S('k2'))
        k1=self.ksym
        k2=''
        k3=''
        kname=''
        qq=len(args)
        if qq==0:
            kres=sqrs(k1=k1,k2=k2,k3=k3)
             
        elif qq==1:
            if type(args[0])==str:
                kname=args[0]
            else:
                k2=args[0]
        elif qq==2:
            if type(args[0])==str:
                kname=args[0]
                k2=args[1]
                
            else:
                k2=args[0]
                k3=args[1]
        else:
            kname=args[0]
            k2=args[1]
            k3=args[2]
            
        kres=sqrs(k1=k1,k2=k2,k3=k3)
        if ope!='':
            kres = opemat(kres, kope=kope)
        
        if kname!='':
            ee=MyEq(kres,kname)
            return ee
            
        if kupdate:
            self.ksym=kres
            if kshow:
                self.s()
        else:
            return kres    
        
    
    def Rpow(self,kpo,kname='',kshow=True):   
        kres=self.ksym
        if kpo==2:
            kres=sqrt(kres)
        else:    
            kres=ppow(kres,1,kpo)
        self.ksym=kres
        self.s()

    ###########################################
    #            Math operation
    #  Add, SUbs, Mul, Div, Pow, Rpow
    ###########################################


    def cut_fac(self, kval):
        kres = cut_fac(self.ksym, kval)
        self.ksym=kres
        return kres

    def signo(self):
        kres = self.ksym
        return signo(kres)
        
    def partialfraction(self,gen=''):
        if gen!='':
            var=gen
        else:
            var=self.var
        kres=self.ksym
        if Is_Add(kres):
            total=0
            mm=fpoly(kres,'list')
            for i in mm:
                total+=partialfraction(i,var)
            return total    
        kres=partialfraction(kres,var)
        self.ksym=kres
        self.s()
        
        
    def par_frac(self,var='',kshow=True):
        try:
            kres=par_frac(self.ksym,var=var)
            self.ksym=kres    
        except:
            pass
        if kshow:
            self.s()
    
    @property
    def float(self):
        try:
            return float(self.ksym)
        except:
            return self.ksym
            
    ###########################################
    #            get info
    ###########################################
     
    def baselist(self): # if is Monomie return list of bases
        kres=self.ksym
        return baselist(kres)
    
    def expolist(self): # if is Monomie return list of exponenets
        kres=self.ksym
        return expolist(kres)    
          
    
    def get_dview(self):
        kres=self.ksym
        for i in self.vfunc:
            kres=show_modefunc(kres,i,ktype=2)
        return kres
        
    def get_pview(self):
        kres=self.ksym
        for i in self.vfunc:
            kres=show_modefunc(kres,i,ktype=3)
        return kres 
        
        
    def list(self, kopt=''): #@tif
        kres = self.ksym
        kres = fpoly(kres, 'list')
        if kopt != '':
            return kres[kopt]
        else:
            return kres
    
    def slist(self, kopt=''):  # symbols list
        kres = self.ksym
        kres = fpoly(kres, 'list')
        vres=[]
        for i in kres:
            vres.append(num_ksym2ksym(i))
        if kopt != '':
            return vres[kopt]
        else:
            return vres


    def free(self): #@tif
        kres = self.ksym
        kres = fpoly(kres, 'free')
        return kres
        
    def get_args(self,*args): #@tif
        kres=self.ksym
        for i in args:
            kres=kres.args[i]
        return kres
    
    
    def args(self,*args,deep=2,format='list'):
        if len(args)==0:
            showarglist(self.ksym,deep=deep,format=format)
        else: 
            kres=self.ksym
            for i in args:
                kres=kres.args[i]
            return kres

    def get_primitive(self,kname='',kshow=True):
        kres=self.ksym
        if self.varf!=[]:

            ee=self.xcopy(kname,kshow=False)
            nf=[Function(x)(ee.var2) for x in ee.varf]
            for i,j in zip(ee.varf,nf):
                ee.set(i,j,kshow=False)
            
            self.primi_eq=ee.ksym
            if kname!='':
                ee.s()
                return ee
            else:
                return ee.ksym
        else:
            self.s()

    

    
    def forceSetValue(self, knom, kval):

        sknom = str(knom)
        skval = str(kval)
        skvalue = str(self.ksym)
        skvalue.replace(sknom, skval)
        kfac = 1
        if skavlue[0] == '-':
            skvalue = skvalue[1:-1]
            kfac = -1
        kres = parse_expr(skvalue)
        kres = kres * kfac
        self.ksym=kres
        self.s()



    def kdiff(self, kval, kope='', kupdate=False):
        kres = kdiff(self.ksym, kval)
        kres = opemat(kres, kope)
        if kupdate:
            self.ksym=kres
            self.s()
        else:
            return kres
            
    
    def integrate(self,var1=''):
        P1=self.ksym
        if var1=='':
            var1=self.var
         
        v1,v2=varDiff(str(var1),'vv')
        P1=P1.subs(v1,1)
         
        sv1='d'+str(v1)
         
        P1=P1.subs(sv1,1)
         
        ee1=MyIntg(P1,'ee1',var=var1,kshow=False)
         
        self.ksym=ee1.ksym
         
        self.s()    
    def integral(self,kname='',x1='', x2='', var=''):
        from lib_MyIntegral import MyIntg 
        if kname=='':
            name2=self.name
            name2=name2.replace('d','')
            kname=name2
        var=self.var
        obj=obj2func(self)

        ee=MyIntg(obj,kname,var=var,x1=x1,x2=x2)
        ee.varc=True
        return ee   
                
        

    def modo_integral(self,*args):
        Iee=self.xcopy('iee',kshow=False)
        kres=Iee.ksym
        kname=''
        for i in args:
            if type(i)==str:
                kname=i
            else:
                ival=i[0]
                delI=symbols(diffname(ival))
                Iee.set(delI,1,kshow=False)
        kres=Iee.ksym
        for i in args:
            if type(i)!=str:
                if len(i)==1:
                    kres=Integral(kres,i[0])
                elif len(i)==2:
                    kres=Integral(kres,(i[0],0,i[1]))
                else:
                    kres=Integral(kres,(i[0],i[1],i[2]))
        if kname=='':
            self.ksym=kres
            self.s()
        else:
            newee=MyEq(kres,kname=kname)
            return newee

    def tintegral_def(self, alpha1, a1, a2, kupdate=False):
        keq = self.ksym
        kres = tintegral_def(keq, alpha1, a1, a2)
        if kupdate:
            self.ksym=kres
            self.s()
        else:
            return kres

    def only_nume(self, kope=''):
        kres = self.get_nume()
        kres = opemat(kres, kope)
        self.ksym=kres
        self.s()

    def only_deno(self, kope=''):
        kres = self.get_deno()
        kres = opemat(kres, kope)
        self.ksym=kres
        rself.s()

    def get_MonoExp(self):
        kres2 = self.ksym
        kres = kres2.fpoly('get', 1)
        return kres

    # def evalue(self, *args, kope='', kshow=False, **kwargs):
        # qk = len(kwargs)
        # if qk > 0:
            # for key, value in kwargs.items():
                # self.set(parse_expr(key), value, kshow=False)
                # try:
                    # self.primi = self.primi.subs(parse_expr(key), value)
                # except:
                    # done = False
        # if len(args) == 1:
            # var2 = self.var2

            # kres = self(var2=args[0])
        # return kres

    def evalueArray(self, knom, vkval):

        kres = [self.evalue(knom, xx, kshow=False) for xx in vkval]
        return kres

    ###########################################
    #            Functions
    ###########################################
    
    def rad2sex(self): #convert rad 2 sex
        kres=self.ksym
        try:
            kres=rad2sex(kres)
        except:
            pass
        self.ksym=kres
        self.s()

    def sex2rad(self): #convert sex 2 rad
        kres=self.ksym
        try:
            kres=sex2rad(kres)
        except:
            pass
        self.ksym=kres
        self.s()
            
    
    
    @property
    def Type(self):
        return type(self.ksym)
        
    def slope(self, xx, kval='', kope=''):
        xExp = self.ksym

        kres = self.kdiff(xx)
        if kval != '':
            kres = kres.subs(xx, kval)
        kres = opemat(kres, kope)
        return kres

    def slopeO(self, xx, kval='', kope=''):
        xExp = self.ksym

        kres = self.kdiff(xx)
        kres = -1 / kres
        if kval != '':
            kres = kres.subs(xx, kval)
        kres = opemat(kres, kope)
        return kres

    def get_aTan(self, kname='', **kwargs):
        if kname != '':
            return MyEq(self.angTan(**kwargs), kname=kname)

        return self.angTan(**kwargs)

    def angTan(self, kname='', **kwargs):
        kres = self.ksym
        var2 = self.var2
        kres = kres.diff(var2)
        if len(kwargs) > 0:
            ee = MyEq(kres, kshow=False)
            kres = ee(**kwargs)
            if len(fpoly(kres, 'list')) > 1:
                mm = fpoly(kres, 'list')
                qq = len(mm)
                kres = mm[qq - 1]
        if kname != '':
            return MyEq(atan(kres), kname)
        else:
            return atan(kres)

    def get_aOrto(self, kname='', **kwargs):
        if kname != '':
            return MyEq(self.angOrto(**kwargs), kname=kname)
        return self.angOrto(**kwargs)

    def angOrto(self, kname='', **kwargs):
        kres = self.angTan(**kwargs)
        if kname != '':
            return MyEq(kres, kname)
        else:
            return (kres + pi / 2)

    def ang_vecTan(self, **kwargs):
        kres = ee.diff(**kwargs)
        kres = atan(kres)
        return kres

    def ang_vecOrto(ee, **kwargs):
        kres = ee.diff(**kwargs)
        kres = atan(kres)
        return kres + pi / 2

    def Is_Poly(self):
        return Is_Poly(self.ksym)

    def Is_Mono(self):
        return (Is_Mono(self.ksym))

    def Is_Add(self):
        kres = self.ksym
        if type(kres) == Add:
            return True
        else:
            return False

    def Is_Mul(self):
        kres = self.ksym
        if type(kres) == Mul:
            return True
        else:
            return False

    def Is_Pow(self):
        kres = self.ksym
        if type(kres) == Pow:
            return True
        else:
            return False

    ###########################################
    #            Transformation
    ###########################################

    def transformada(self,expr0,ssym,kope=''):
        kres=self.ksym
        try:
            kres=transformada(kres,expr0=expr0,ssym=ssym,kope=kope)
            self.ksym=kres
        except:
            pass
        self.ksym=kres
        self.s()
    
    def simplifyroot(self,kshow=True):
     
        kres=self.ksym
        kres=simplifyroot(kres)
        self.ksym=kres
        if kshow:
            self.s()
            
        
    def subsnumber(self,val1,val2,kope='',kshow=True):
        kres=self.ksym
        kres=subsnumber(kres,val1,val2)
        if kope!='':
                kres=opemat(kres,kope=kope)
        self.ksym=kres
        if kshow:
                self.s()
    
    def all_type(self):
        self.s()
        sE(['Monomie= ', self.Is_Mono(), '  ', 'Polynomie = ', self.Is_Poly()])
        sE(['Is Add= ', self.Is_Add(), '  ', 'Is Mul= ', self.Is_Mul()])

    
    def addexpand(self):
        self.ksym=addexpand(self.ksym)
        self.s()
    
    
    def expand(self, kname='',kshow=True):  #@tif  
        kres = self.ksym
        kres=expand(kres)
        if kname!='':
            ee=MyEq(kres,kname=kname,kshow=kshow)
            return ee
        else:
            self.ksym=kres
            if kshow:
                self.s()
    def expanddenom(self, kname='',kshow=True):  #@tif  
        kres = self.ksym
        p1=numer(kres)
        p2=denom(kres)
        p2=expand(p2)
        
        kres=p1/p2
        if kname!='':
            ee=MyEq(kres,kname=kname,kshow=kshow)
            return ee
        else:
            self.ksym=kres
            if kshow:
                self.s()
    def expandnumer(self, kname='',kshow=True):  #@tif  
        kres = self.ksym
        p1=numer(kres)
        p2=denom(kres)
        p1=expand(p1)
        
        kres=p1/p2
        if kname!='':
            ee=MyEq(kres,kname=kname,kshow=kshow)
            return ee
        else:
            self.ksym=kres
            if kshow:
                self.s()
    def factor(self,*args, kshow=True ):  
        kres = self.ksym
        if 'sqrt' in args:
            kres=factor(kres,extension=sqrt(2))
        elif len(args)==1 and type(args[0])!=str:
            kres=collect(kres,args[0])
        else:
            kres=factor(kres)
        self.ksym=kres
        if kshow:
            self.s()        
    def apart(self,kshow=True):
         
        kres=self.ksym
        kres=apart(kres)
        self.ksym=kres
        if kshow:
            self.s()
            
    def inverse(self,kshow=True):
         
        kres=self.ksym
        kres=iinverse(kres)
        self.ksym=kres
        if kshow:
            self.s()    
    
    def simplifybase(self,kname='',kshow=True): # simplify each ecponet in expr
        ksym = self.ksym
        kres=simplifybase(ksym)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            self.ksym=kres
            self.s()
                
    def simplify(self):
        kres = self.ksym
        kres = simplify(kres)
        self.ksym=kres
        self.s()
            
    def simplify_img(self,kshow=True):
        kres = self.ksym
        kres = simplify_img(kres)
        self.ksym=kres
        if kshow:
            self.s()
    def lmul2lpow(self, kupdate=True, kshow=True, kope=''):
        kres = self.ksym
        kres = lmul2lpow(kres)
        kres = opemat(kres, kope=kope)
        if kupdate:
            self.ksym=kres
        if kshow:
            self.s()
    
    
    def lsimplify(self, kupdate=True, kshow=True, kope=''):
        kres = self.ksym
        kres = lsimplify(kres)
        kres = opemat(kres, kope=kope)
        if kupdate:
            self.ksym=kres
        if kshow:
            self.s()        
    
    def numerexpand(self):
        expr=self.ksym
        p1,p2=fraction(expr)
        p1=expand(p1)
        self.ksym = cfrac(p1,p2)
        self.s()
    def numerfactor(self):
        expr=self.ksym        
        p1,p2=fraction(expr)
        p1=factor(p1)
        self.ksym = cfrac(p1,p2)
        self.s()        
    def numersimplify(self):
        expr=self.ksym
        p1,p2=fraction(expr)
        p1=simplify(p1)
        self.ksym = cfrac(p1,p2)    
        self.s()
    def denoexpand(self):
        expr=self.ksym
        p1,p2=fraction(expr)
        p2=expand(p2)
        self.ksym = cfrac(p1,p2)
        self.s()
    def denofactor(self):
        expr=self.ksym
        p1,p2=fraction(expr)
        p2=factor(p2)
        self.ksym = cfrac(p1,p2)    
        self.s()
    def denosimplify(self):
        expr=self.ksym
        p1,p2=fraction(expr)
        p2=simplify(p2)
        self.ksym = cfrac(p1,p2)
        self.s()
    
    
    def dfactor(self, var,var1,op=''):
        kres=self.ksym
        kres = dfactor(kres,var=var,var1=var1,op=op)
         
        self.ksym=kres
        self.s()
         
    def dsimplify(self, var,var1):
        kres=self.ksym
        kres = dsimplify(kres,var=var,var1=var)
         
        self.ksym=kres
        self.s()
        
    def dothis(self,*args):
        ksym=self.ksym
        args2=[ksym]
        for i in args:
            args2.append(i)
        kres=dothis(*args2)
        self.ksym=kres
        self.s()
    def doindenom (self,*args):
        ksym=self.ksym
        func=args[0]
        if len(args)==1:
            kres=doindenom(ksym,func)
        else:
            expr2=args[1]
            kres=doindenom(ksym,func,expr2)
        
        self.ksym=kres
        self.s()   
    def doinnumer (self,*args):
        ksym=self.ksym
        func=args[0]
        if len(args)==1:
            kres=doinnumer(ksym,func)
        else:
            expr2=args[1]
            kres=doinnumer(ksym,func,expr2)
        
        self.ksym=kres
        self.s()     
    def simplifyexpand(self, kupdate=True, kshow=True, kope=''):
        kres = self.ksym
        kres = simplifyexpand(kres)
        kres = opemat(kres, kope=kope)
        if kupdate:
            self.ksym=kres
        if kshow:
            self.s()        
            
    def expandbase(self,kname='',kshow=True):
        ksym = self.ksym
        kres=expandbase(ksym)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            self.ksym=kres
            self.s()
    
    def disjoinbase(self,kname='',kshow=True):
        ksym = self.ksym
        kres=disjoinbase(ksym)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            self.ksym=kres
            self.s()
            
            
    def disjoinexpo(self,kname='',kshow=True):
        ksym = self.ksym
        kres=disjoinexpo(ksym)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            self.ksym=kres
            self.s()        
            
    def joinbase(self,kname='',kshow=True):
        ksym = self.ksym
        kres=joinbase(ksym)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            self.ksym=kres
            self.s()
            
    def joinexpo(self,kshow=True):
        ksym = self.ksym
        kres=joinexpo(ksym)
        if kshow:
            self.s()        
    def separebase(self,kname='',kshow=True):
        ksym = self.ksym
        kres=separebase(ksym)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            self.ksym=kres
            self.s()

            
    def sum2mulexpo(self,kname='',kshow=True):
        ksym = self.ksym
        kres=sum2mulexpo(ksym)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            self.ksym=kres
            self.s()        
    def expandexpo(self,kname='',kshow=True):
         
        ksym = self.ksym
        kres=expandexpo(ksym)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            self.ksym=kres
            self.s()
    def rsimplify(self,kname='',kshow=True):
         
        ksym = self.ksym
        kres=rsimplify(ksym)
        kres=simplify(kres)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            self.ksym=kres
            self.s()        
    
    def swappow2pow(self,kname='',kshow=True):
         
        ksym = self.ksym
        kres=swapPow2Pow(ksym)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            self.ksym=kres
            self.s()

    def swapPow2Pow(self,kname='',kshow=True):
         
        ksym = self.ksym
        kres=swapPow2Pow(ksym)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            self.ksym=kres
            self.s()

    def pow2powpow(self,*args,kshow=True):
        ksym = self.ksym
        exp1=''
        ktype='out'
        kname=''
        for i in args:
            if type(i)==str:
                ktype=i
            else:
                exp1=i 
        args=[ksym,exp1,ktype]        
        kres=pow2powpow(*args)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            self.ksym=kres
            self.s()      
    def packexp(self,kshow=True):
        kres=self.ksym
        kres=packexp(kres)
        self.ksym=kres
        if kshow:
            self.s()        
     
    def base2frac(self,kshow=True):
        kres=self.ksym
        kres=base2frac(kres)
        self.ksym=kres
        if kshow:
            self.s()
    
    def basefactor(self,kope='',kshow=True):
        return self.simplifyexp(kope=kope,kshow=kshow)
        
    
    def simplifyexpo(self,kname='',kshow=True):
        ksym = self.ksym
        kres=simplifyexpo(ksym)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            return kres
    def primefactors(self,kshow=True):
        kres=self.ksym
        kres=primefactor(kres)
         
        self.ksym=kres
        if kshow:
            self.s()
    def primefactor(self,kshow=True):
        self.ksym=UnevaluatedExpr(primefactor(self.ksym))
        self.tipoView=True
        if kshow:
            self.ss()
        
    
     


            
    def div2mulexp(self,kope='',kshow=True):
        kres=self.ksym
        kres=div2mulexp(kres)
        if kope!='':
            kres=opemat(kres,kope=kope)
        self.ksym=kres
        if kshow:
            self.s()
            
    def simplifyrpow(self, kupdate=True, kshow=True, kope=''):
        kres = self.ksym
        kres = simplifyrpow(kres)
        kres = opemat(kres, kope=kope)
        if kupdate:
            self.ksym=kres
        if kshow:
            self.s()
    def reducePow(self, kupdate=True, kshow=True, kope=''):
        kres = self.ksym
        kres = reducePow(kres)
        kres = opemat(kres, kope=kope)
        if kupdate:
            self.ksym=kres
        if kshow:
            self.s()
            
    def simplify_sec(self,kshow=True):
        kres = self.ksym
        if self.Is_Add():
            sres = 0
            mm=fpoly(kres, 'list')
            for i in mm:
                sres += simplify(i)
            self.ksym=sres    
             
        self.update(self.ksym)
        if kshow:
            self.s()
            
         
    def roots(self,kshow=True):
        return self.solve(self.var,'all',kshow=kshow) 
  
    def root(self):
        return self.solve(self.var,'all')
        
    def linfactor(self, kvar,kshow=True):
        p1=self.ksym
        p1=linfactor(p1,kvar)
        self.ksym=p1
        if kshow:
            self.s()

    
    def cancel(self, kvar='', kshow=True, kope=''):
        kres = self.ksym
        if kvar != '':
            g1 = self.fpoly('filt', kvar)
            g2 = kres - g1
            g3 = factor(g1)
            kres = g3 + g2
        else:
            kres = cancel(kres)
        kres = opemat(kres, kope=kope)
        self.ksym=kres
        if kshow:
            self.s()
            
    def texpand(self, kshow=True, kope=''):
        kres = self.ksym
         
                
        kres = texpand(kres)
        kres = opemat(kres, kope=kope)
        self.ksym=kres
        if kshow:
            self.s()

    def tsimplify(self, kshow=True, kope=''):
        kres = self.ksym
        kres = trigsimp(kres)
        kres = opemat(kres, kope=kope)
        self.ksym=kres
        if kshow:
            self.s()
    def tfunc53(self,alpha=alpha,kshow=True):
        kres=self.ksym
        kres=tfunc53(kres,alpha)
        self.ksym=kres
        if kshow:
            self.s()
            
    def tfunc37(self,alpha=alpha,kshow=True):
        kres=self.ksym
        kres=tfunc37(kres,alpha)
        self.ksym=kres
        if kshow:
            self.s()        
    def tfunc16(self,alpha=alpha,kshow=True):
        kres=self.ksym
        kres=tfunc16(kres,alpha)
        self.ksym=kres
        if kshow:
            self.s()
       
    def tfunc74(self,alpha=alpha,kshow=True):
        kres=self.ksym
        kres=tfunc74(kres,alpha)
        self.ksym=kres
        if kshow:
            self.s()
            
    def lexpand(self, kshow=True):
        kres = self.ksym
        kres = lexpand(kres)
        self.ksym=kres 
        if kshow:
            self.s() 
    def lexponent(self):
        expr=self.ksym
        expr=lexponent(expr)
        self.ksym=expr 
        self.s()
            
    def lfactor(self, kshow=True, kope=''):
        kres = self.ksym
        kres = logcombine(kres,force=True)
        self.ksym=kres
         
        if kshow:
            self.s()
    def positivexpo(self,kname='',kshow=True,force=False):
        ksym = self.ksym
        kres=positivexpo(ksym,force=force)
        if self.get_deno()!=1:
            kres=cfrac(positivexpo(self.get_nume(),force=force),positivexpo(self.get_deno(),force=force))
        
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            self.ksym=kres
            self.s()    
    def exp(self, kshow=True, kope=''):
        kres = self.ksym
        kres = exp(kres)
        kres = opemat(kres, kope=kope)
        self.ksym=kres
        if kshow:
            self.s() 
    def lexponent(self,kname='',kshow=True):
        ksym = self.ksym
        kres=lexponent(ksym)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            self.ksym=kres
            self.s()
    def log(self, kshow=True, kope=''):
        kres = self.ksym
        kres = log(kres)
        try:
            kres=lsimplify(kres)
        except:
            pass        
        kres = opemat(kres, kope=kope)
        self.ksym=kres
        if kshow:
            self.s()    

    def opemat(self, kope='', kshow=True):
        kres = self.ksym
        kres = opemat(kres, kope=kope)
        self.ksym=kres
        if kshow:
            self.s()

    def opematsec(self, kope='', kshow=True):  # equal to pemat but secuential
        kres = unisymbols(self.ksym)
        kres = opematsec(kres, kope=kope)
        self.ksym=kres
        if kshow:
            self.s()

    def simplify_Rpow(self, kshow=True):
        kres = self.ksym
        kres = cut_root_of_pow(kres)

        self.ksym=kres
        if kshow:
            self.s()

    def sort(self):
        kres = self.ksym
        klist = self.fpoly('list')
        mm = 0
        for i in klist:
            mm += i
        kres = mm
        self.ksym=kres
        self.s()

    def expandExp(self):
        kres2 = self.ksym
        kvar = fpoly(kres2, 'get', 0)
        kexp = fpoly(kres2, 'get', 1)
        keList = fpoly(kexp, 'list')
        mm = 0
        for i in keList:
            mm += kpow(kvar, i)
        return mm

    def get_inside_root(self, kname=''):
        mm = str(self.ksym)
        cc = 0
        kk = mm.find('sqrt')
        kk2 = kk + 4
        sres = ''
        qq = len(mm)
        for i in range(kk2, qq):
            val = mm[i]
            if val == '(':
                cc += 1
            if val == ')':
                cc -= 1
            sres = sres + val
            if cc == 0:
                kres = parse_expr(sres)
                if kname != '':
                    ee = MyEq(kres, kname)
                    return ee
                else:
                    return kres
        return self.ksym
    def get_type(self):
        return type(self.ksym)
    
    
    @property    
    def numer(self):
        kres=self.ksym
        return kreturn(numer(kres))
    
    @property
    def get_nume(self):
        return kreturn(numer(self.ksym))


    @property
    def denom(self):
        kres=self.ksym
        return kreturn(denom(kres))
        

    def get_deno(self):
        return  kreturn(denom(self.ksym))

    def rem(self,val2=''):
        if val2!='':
            kres=self.ksym
            return rem(self.ksym,val2)
        else:    
            p1=self.numer
            p2=self.denom
            return rem(p1,p2)

    def diffValue(self, kval=''):
        if kval == '':
            kval = self.var2
        kres = diff(self.ksym, kval)
        return kres

    def fpoly(self, kopt='', op2='', op3=''):
        kres = self.ksym
        return fpoly(kres, kopt=kopt, op2=op2, op3=op3)
        
    def getexpo(self,kname='',kshow=True):
        ksym = self.ksym
        kres=getexpo(ksym)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:  
            return kres
    
    def getexponent(self,kname='',kshow=True):
        return self.getexpo(kname=kname,kshow=kshow) 
        
    def get_expo(self,kname='',kshow=True):  # return exponente from monomie expresion
        return self.getexpo(kname=kname,kshow=kshow) 
    
    def getbase(self,kname='',kshow=True):
        ksym = self.ksym
        kres=getbase(ksym)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            return kres

    
    def get_killexpo(self):
        ksym = self.v
        if Is_Mono(ksym):
            mm = fpoly(ksym, 'list')
            return mm[1]

    def get_sitem(self, vv=[]):
        mm = self.list()
        qq = len(mm)
        ksum = 0
        for i in range(qq):
            if i in vv:
                ksum += mm[i]
        return ksum

    def part(self, vec):

        kres = self.ksym
        try:
            return part(kres, vec)
        except:
            return kres

    def solvediff(self, kvar, kd='', kremp=False, kope='', korden=''):
        return self.solve(kvar=kvar.diff(), kd=kd, kremp=kremp, kope=kope, korden=korden)

    def solve_if_and(self, svar, eqv=0, kope='',korden='',kshow=True, **kwargs):
        r'''
        solve variable from  MyEq
        parameters :
            svar :type str , variablesin side the Eq taht we will find
            eqv  :type nemeric or symbols , if the value of all Eq
                  defaul Eq=0
            kwargs: t=0,g=10... etc
        return MyEq of svar
        example:
        **********
        R(t)= C1 + C2*t + g*sin(t*w)/w**2
        C1= solve_if_and('C1',L,t=0)
        return:  C1=L

        R.upgrade(C1)
        return:  C2*t + L + g*sin(t*w)/w**2

        C2=solve_if_and(R,'C2',t=2)
        return:-L/2 - g*sin(2*w)/(2*w**2)

        R.upgrade(C2)
        return: L + g*sin(t*w)/w**2 + t*(-L/2 - g*sin(2*w)/(2*w**2))

        '''
        kres = self.ksym - eqv
        if len(kwargs) > 0:
            kres = self(**kwargs) - eqv
        nee = MyEq(kres, kname=svar, kshow=False, kope=kope,)
        kres=nee.ssolve(svar,korden=korden,kshow=kshow,kope=kope)
        self.set(svar,kres,kshow=kshow,kope=kope)
        return kres

    def solveif(self,var,**kwargs):
        kres=self.ksym-parse_expr(self.name)
        kres=real_subs(kres,**kwargs)
        ee=MyEq(kres,'ee',kshow=False)
        rvar=ee.solve(var,kname=str(var))
        return rvar
        
    def get_symbol_in_ee(self, svar=''):
        mm = self.free()
        for i in mm:
            if i.name == svar:
                return i
        return svar

    def left(self):
        return symbols(self.name)

    def right(self):
        return  self.ksym

    def ssolve(self, kname,korden='',kshow=True,kope=''):

        kvar = self.get_symbol_in_ee(svar=kname)
        return self.solve(kvar, kname=kname,korden=korden,kshow=kshow)

    ###########################################
    #            SOLVE
    ###########################################
    
    #   solve()
     
    def solve(self,*args,kshow=True,**kwargs):
        var=args[0]
        kres=self.ksym
         
            
        if len(kwargs)>0:
            vecs,vecv=unpack(kwargs)
            if self.name in vecs:
                Y=0
                for i,j in zip(vecs,vecv):
                    if i==self.name:
                        Y=j
                        kres=kres-Y
            kres=real_subs(kres,**kwargs)
        ksolu=solve(kres,var)
        newargs=[]
        if len(args)>1:
            newargs=args[1::]
             
        if 'float' in newargs:
            nksolu=[]
            for i in ksolu:
                try:
                    nksolu.append(float(i))
                except:
                    nksolu.append(i)
            ksolu=nksolu
        if 'signo' in newargs:
            nksolu=[]
            for i in ksolu:
                try:
                    nksolu.append(-1*i)
                except:
                    nksolu.append(i)
            ksolu=nksolu
        if 'all' in newargs:
            if kshow:
                display(Math(latex(ksolu)))
            return ksolu
         
        if 'noimg' in newargs or 'nonimaginary' in newargs:
            return [i for i in ksolu if not 'I' in str(i)]
        try:
            ksolu=ksolu[0]
        except:
            ksolu=ksolu
        if 'value' in newargs:
            return ksolu
           
         
        if 'hide' in newargs:
            kshow=False
        
        
        var=MyEq(ksolu,str(var),kshow=kshow)    
        if 'update' in newargs:
            self.set(var)
         
        else:
            return var
         
    
    
    
    # def solve(self,*args,kshow=True,**kwargs):
        # if len(args)==0:
            # helplib('solve')
            # return
        # var=args[0]
        # args=args[1:len(args)]        
        # self2=copy.deepcopy(self.ksym)
        # kres=self.ksym
        # eQ=copy.deepcopy(kres)
        # try:
            # kname=str(var) 
        # except:
            # kname=str(var)
        # if len(kwargs)>0:
             
            # mm=[]
            # for i in kwargs:
                # mm.append(i)
               
            # if self.name in mm:
                 
                # kres=self.ksym-symbols(self.name)
                 
        # if len(kwargs)>0:
            # kres=real_subs(kres,**kwargs)
            
             
        # solu=ksolve(kres,var)
         
        
        # if len(args)>0:
            # if 'all' in args:
                # solu=solve(kres,var)
                # if 'float' in args:
                    # try:
                        # solu=float(solu)

                    # except:
                        # solu=[float(i)  for i in solu]

                    # display(Math(latex(solu)))
                    # return solu    
                    
                # return solu
            # if 'float' in args:
                # try:
                    # solu=float(solu)
                    
                # except:
                    
                    # pass
                
            # if 'value' in args:
                # return solu
            # elif'update' in args:  
                # eQ=self2.subs(var,solu)
                # self.ksym=eQ
                # self.s()
                # ee=MyEq(solu,kname=str(var),kshow=kshow)
                # return ee
        
         
        # ee=MyEq(solu,kname=str(var),kshow=kshow)
        # return ee
         
    def solveset(self, *args,kshow=True,**kwargs):
        kres=self.solve(*args,kshow=True,**kwargs)
        self.upgrade(kres)
        return kres
        
        
    def solvediffsymbol(self,var2):
        var=self.var
        V2=symbols(str(var2))
        dy1,dy2=primesymboldiff(var2)
        yy=Function(str(var2))(var)
        sv2f=str(var2)+'('+str(self.var)+')'
        sv2=str(var2)
        kres=self.ksym
        dy=yy.diff(var)
        ee=self.solve(dy,kshow=False)
        ee.name=str(dy1)
         
        kres=ee.ksym
        kres2=subsubs(kres,sv2f,sv2)
        ee.ksym=kres 
        
        ee.s()
        return ee


    def solveB(self, kvar, kname='', kremp=False, kope='', korden='', ktype='P', var2='', Bag='',kpositive=False,noncero=False,kshow=True):
        mm = self.ksym
        self.ksym = opemat(mm, 's')
        if type(kvar) == 'str' and kname == '':
            kname = kvar
            kvar = parse_expr(kname)

        if str(kvar) == 'None':
            kname = symbols(kname)
        keq = self.ksym
        # kvar=unisymbols(kvar)
        if noncero and kname!='':
            nV=symbols(self.name)
            keq=keq-nV
            kname=nV.name
        kres = csolve(keq, kvar, kope=kope, korden=korden)
        if kres == []:
            kres = csolve(keq, kname, kope=kope, korden=korden)
            if kres == []:
                kres = presolve(self.ksym, kvar)

        if self.type == 'Ph':
            if not str(kvar) in 'alphaalpha1alpha2alpha3betha1betha2betha3':
                self.Pobj.store_val(kvar, kres)
                mm = self.Pobj.F
                qq = len(mm)
                for i in range(qq):
                    pkres = mm[i][0]
                    try:
                        pkres = pkres.subs(kvar, kres)
                    except:
                        done = False
                    mm[i][0] = pkres
                self.Pobj.F = mm
        if kname != '':
            if type(kres) == list:
                if len(kres) > 1:
                    if kpositive:
                        kres1 = MyEq(kres[1], kname, var2=self.var2)
                        return kres1

                    cc = 1
                    kkres = []
                    for i in kres:
                        kname1 = kname + str(cc)
                        kres1 = MyEq(i, kname1, kshow=False, var2=self.var2)
                        if Bag != '':
                            kres1.upBag(Bag, kope=kope)

                        kres1.s()
                        cc += 1
                        kkres.append(kres1)

                    return kkres

            else:

                kres = MyEq(kres, kname, kshow=False, var2=self.var2)
                if Bag != '':
                    kres.upBag(Bag, kope=kope)
                if kshow:
                    kres.s()

                return kres


        else:
            return kres

    def ssolveR(self, kname):
        kvar = self.get_symbol_in_ee(svar=kname)
        return self.solveR(kvar, kd=kname)

    def solveR(self, kvar, kd='', kremp=False, kope='', korden='', Bag=''):
        if Bag != '':
            self.upBag(Bag=Bag, kshow=False)
        keq = self.ksym

        try:
            kres = csolveR(keq, kvar, kope=kope)
        except:
            kres = csolve(keq, kname, kope=kope, korden=korden)
        if kd != '':
            return MyEq(opemat(kres, kope=kope), kd, var2=self.var2)
        else:
            return kres

    def toFoat(self):
        kres = self.ksym
        try:
            kres = float(self.ksym)

            self.ksym=kres
            self.s()
        except:
            self.s()

    def toEqual(self, items, kkname=''):
        mm = self.list()
        kres1 = 0
        kres2 = 0
        qq = len(mm)
        for i in range(qq):
            if i in items:
                kres1 += mm[i]
            else:
                kres2 += mm[i]
        kname = self.name
        kname1 = kname + '1'
        kname2 = kname + '2'
        kname1 = MyEq(kres1, kname1, kshow=False)
        kname2 = MyEq(kres2, kname2, kshow=False)
        return MyEqEq(kname1, kname2, kname=kkname)

    def opematsec(self, kope=''):  # equal to pemat but secuential
        kres = unisymbols(self.ksym)
        kres = opematsec(kres, kope=kope)
        self.ksym=kres
        return kres

    def opemat_deno(self, kope=''):
        kres = unisymbols(self.ksym)
        kres = opemat_deno(kres, kope=kope)
        self.ksym=kres
        return kres

    def opemat_nume(self, kope=''):
        kres = unisymbols(self.ksym)
        kres = opemat_nume(kres, kope=kope)
        self.ksym=kres
        return kres

    def solve_tan(self, alpha):
        c = tan(alpha)

        kres = self.ksym
        kres.subs(sin(alpha), c / rpow(c * c + 1, 2))
        kres.subs(cos(alpha), 1 / rpow(c * c + 1, 2))
        kres1 = csolve(kres, tan(alpha))
        return kres1
        
    def solve_angle(self,angle=alpha,kname='Ang'):

        if type(angle)==str:
            kname=angle
            angle=alpha

        valor= self.solve(angle,kshow=False)
        ee=MyEq(valor.ksym,kname)
        return ee    

    ###########################################
    #             Polinomial
    ###########################################
    def quo(self,val2=''):
        if val2!='':
            kres=self.ksym
            return quo(self.ksym,val2)
        else:    
            p1=self.numer
            p2=self.denom
            return quo(p1,p2)    
    def quotient(self,ksym,kname=''):
        
        P=self.ksym
        Q=ksym
        if type(Q)==MyEq:
            Q=ksym.ksym
        kres=quo(P,Q)
        if kname!='':
            kres=MyEq(kres,kname)
            return kres
        else: 
            return kres
        
    def rem(self,ksym,kname=''):
        
        P=self.ksym
        Q=ksym
        if type(Q)==MyEq:
            Q=ksym.ksym
        kres=rem(P,Q)
        if kname!='':
            kres=MyEq(kres,kname)
            return kres
        else: 
            return kres
        
    def coef_sum(self,var2=''):
        if var2=='':
            var2='x'
        kres=self.ksym
        kres=kres.subs(var2,1)
        return kres
        
    def get_coef(self,ksym,kname=''): # ee= x**2(a+b)+x*4*b+9
                           # get_coef(x) return 4*b
        mm=self.list()
        mmv=[]
        sres=str(ksym)
        for i in mm:
            mmv.append(i/ksym)
        for i,j in zip(mmv,mm):
            if not sres in str(i):
                kres=j
                kres=kres.subs(ksym,1)
        if kname!='':
            return MyEq(kres,kname=kname,var2=self.var2)
        else:
            return kres
            
    def get_ter_inde(self,ksym,kname=''):
        kres=0
        svar=str(ksym)
        mm=self.list()
        for i in mm:
            if not svar in str(i):
                kres+=i
        if kname!='':
            return MyEq(kres,kname=kname,var2=self.var2)
        else:
            return kres        
    
    def degree(self, kvar=0):
        kres = self.ksym
        return degree(kres, gen=kvar)

    def degree_list(self):
        kres = self.ksym
        return degree_list(kres, gen=kvar)
        

    def main_coef(self):
        kres = self.ksym
        return LC(kres)
        
    def main_monomio(self):
        kres = self.ksym
        return LM(kres)
    
    def main_term(self):
        kres = self.ksym
        return LT(kres)
        
    def coef_listK(self,vx=''):
        if vx=='':
            vx=self.var2
        mm=self.list()
        vlist=[]
        for i in mm:
            kres=i.subs(vx,1)
            vlist.append(kres)
        return vlist      
    

    def reduceroot(self,kshow=True):
        kres = self.ksym
        try:
            kres = reduceroot(kres)
            self.ksym=kres
             
        except:
            pass
        if kshow:
            self.s()
            
    def separable(self):
        knumer=numer(self.ksym)
        kdenom=denom(self.ksym)
        kres=rem(knumer,kdenom)
        nnumer=knumer-kres
        p1=simplify(nnumer/kdenom)
        p2=kres/kdenom
        self.ksym=p1+p2
        self.s()
        
    def reduce(self):
        kres = self.ksym
        kres = expand(kres)
        kres = simplify(kres)
        kres = factor(kres)
        kres = simplify(kres)
        try:
            kres = cut_root_of_pow(kres)
            self.ksym=kres
        except:
            self.ksym=kres
        return kres
        
    def reducecero(self,kshow=True):
        kres = self.ksym
        kres2 = reducecero(kres)
        if kres2!=kres:
            self.ksym=kres2
        if kshow:
            self.s()
         



    
    ###########################################
    #             Algebra
    ########################################### 
    def squaresum(self,p1=0,p2=0):
        kres=self.ksym
        if p1!=0:
            kres2=squaresum(kres,p1=p1,p2=p2)
            self.ksym=kres2
            self.s()
        else:
            self.s()
            
            
    # def expandexp(self):
        # kres=self.ksym
        # try:
            # kres=expandexp(kres)
            # self.ksym=kres
        # except:
            # pass
        # self.s()    
            
    def powexpand(self,op='',kshow=True):
        '''
        input (x**(a*b))   ---->   return(x**a)**b
        input (x**(a*b),b)   ---->   return(x**b)**a
        '''
        kres=self.ksym
        self.ksym=powexpand(kres,op=op)
        if kshow:
            self.s()
    
    def killsqrtpow(self,kshow=True):
      
        kres=self.ksym
        self.ksym=killsqrtpow(kres)
        if kshow:
            self.s()
    
    def opematexp(self,kope=''):
        '''
        apply opemat only in exponent monomie
        '''
        ksym=self.ksym
        if Is_Pow(ksym):
            mm=self.list()
            kres=mm[1]
            kres=opemat(kres,kope=kope)
            self.ksym=kpow(mm[0],kres)
            
        if Is_Add(ksym):
            mm=self.list()
            mkres=0
            for i in mm:
                if Is_Pow(i):
                    mm2=fpoly(i,'list')
                    kres=mm2[1]
                    kres=opemat(kres,kope=kope)
                    kres=kpow(mm2[0],kres)
                else:
                    kres=i
                mkres+=kres 
             
            self.ksym=mkres     
        self.s()
    
    def mulexpo(self,kname='',kshow=True,force=False):
        ksym = self.ksym
        kres=mulexpo(ksym,force=force)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            self.ksym=kres
            self.s()
        
    def expandexpA(self,kope=''):
        '''
        input (x**a)**b
        return(x**(a*b))
        '''
        ksym=self.ksym
        kres=expandexpA(ksym)
        self.ksym=kres
        if kope!='':
            self.opemat(kope=kope)
        self.s()
        
    def factorexpo(self,kname='',kshow=True):
        ksym = self.ksym
        kres=factorexpo(ksym)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            self.ksym=kres
            self.s()
    def factorbase(self,kname='',kshow=True):
        ksym = self.ksym
        kres=factorbase(ksym)
        if kname!='':
            return MyEq(kres,kname,var=self.var,kshow=kshow)
        else:
            self.ksym=kres
            self.s()
    def coef_list(self,var2=x):  # return coef list complete according max degree
        '''
        ee=x*x*a+x*b+c
        return [a,b,c]
        
        ee=a*x*x-b
        return [a,0,-b]
        
        '''
        
        kres=self.ksym
        return  coef_list(kres,var2)   
    
    def get_seudocofactor(self,ksym,kname='',var2=x):
        val1=self.ksym
         
        if type(ksym)==MyEq:
            ksym=ksym.ksym
            
        kres,veccc=get_seudocofactor(val1,ksym,var2)
         
         
        if kname!='':
            return MyEq(kres,kname=kname,var2=var2)
        else:
            return kres 
    
    def sortdegree(self,var2='',kope=''):
        '''
        input a*x*x+b*x*x+c*x+d
        return (a+b)x*x+c*x+d
        '''
        if var2=='':
            var2=self.var2
        ksym=self.ksym
        kres=sortdegree(ksym,var2=var2)
        self.ksym=kres
        if kope!='':
            self.opemat(kope=kope)
        self.s()
        
    
    def addexp(self,kope=''):  # input:  x**a * x**b  return x**(a+b)
        kres=self.ksym
        kres=addexp(kres)
        self.ksym=kres
        if kope!='':
            self.opemat(kope=kope)
        self.s()
    
    
    def insidepar(self,ssym=''):
        '''
        args:
            ssym dfault is '' return the first () founded 
            ssym = str( where search..., 'x(', 'sqrt(', 'log(' etc 
        return expresion inside  select function     
        '''
        ksym=self.ksym
        return insidepar(ksym,ssym=ssym)
        
        
    def GCD(self,val1,kname='',kshow=True,var=x):
        if type(val1)==MyEq:
            val1=val1.ksym
        kres=factor(gcd(self.ksym,val1))
        if kname!='':
            return MyEq(kres,kname=kname,var=var)
        else:
            if kshow :
                return kres
    def MCD(self,val1,kname='',kshow=True,var=x):
        if type(val1)==MyEq:
            val1=val1.ksym
        kres=factor(gcd(self.ksym,val1))
        if kname!='':
            return MyEq(kres,kname=kname,var=var)
        else:
            if kshow :
                return kres
    def LCM(self,val1,kname='',kshow=True,var=x):
        if type(val1)==MyEq:
            val1=val1.ksym
        kres=factor(lcm(self.ksym,val1))
        if kname!='':
            return MyEq(kres,kname=kname,var=var)
        else:
            if kshow :
                return kres    
    def MCM(self,val1,kname='',kshow=True,var=x):
        if type(val1)==MyEq:
            val1=val1.ksym
        kres=factor(lcm(self.ksym,val1))
        if kname!='':
            return MyEq(kres,kname=kname,var=var)
        else:
            if kshow :
                return kres 
    ###########################################
    #             Reduccion Algoritmo
    ###########################################
    
    def factorize(self,ksym,kshow=True):
     
       # ee=3*x**7+2*x**3
       #ee.factorize(x**3)
       #return x**3*(3*x**4+2)
       
     
        kres=self.ksym
        kres=factorize(kres,ksym)
        self.ksym=kres
        if kshow:
            self.s()
    def factoriza(self,ksym,kshow=True):
     
       # ee=3*x**7+2*x**3
       #ee.factorize(x**3)
       #return x**3*(3*x**4+2)
       
     
        kres=self.ksym
        kres=factoriza(kres,ksym)
        self.ksym=kres
        if kshow:
            self.s()        
        
    def nformat(self,nd):
        ksym=self.ksym
        try:
            ksym=nformat(ksym,nd)
            self.ksym=ksym
            self.s()
        except:
            self.s()
            
    def reduFac(self, kop='', kshow=True):  # retun Eq=0,= a*b+c*b.. = b(a+c)..=0  then return (a+c)

        kres = self.ksym
        kres = factor(kres)
        kres2 = 1
        if Is_Mono(kres):
            kres = numer(kres)

            mm = fpoly(kres, 'list')
            for i in mm:
                if Is_Poly(i):
                    kres2 = kres2 * i

        if kres2 != 0 and kres2 != 1:
            self.ksym = kres2
        else:
            self.ksym = kres
        if kshow:
            self.s()

    def factorSec(self, *args, kfiltro='.', kshow=True):
        kres = self.ksym
        for i in args:
            kvar=i
            kres = factorSec(kEq=kres, ksym=kvar, kfiltro=kfiltro)
        self.ksym=kres
        if kshow:
            self.s()

    def factorSecV(self, *args):
        for i in args:
            self.factorSec(i, kshow=False)
        self.s()

    def get_cofactor(self, kfactor):
        mm = self.list()
        for i in mm:
            j = fpoly(i, 'list')
            if kfactor in j:
                kres = i / kfactor
                return kres

    def grupFac(self, ksym, kfiltro='.', kshow=True):
        kres = self.ksym
        kres = grupFac(kres, ksym, kfiltro='.')
        self.ksym=kres
        if kshow:
            self.s()

    def get_factor_with(self, kx, kcomple=True):
        eqq = self.ksym
        return get_factor_with(eqq, kx=kx, kcomple=kcomple)

    def monoFactor(self, kval):
        kres = self.ksym
        if Is_Poly(kres):
            mm1 = 0
            mm2 = 0
            mm = fpoly(kres, 'list')
            for i in mm:
                if Is_Mono(i) and Is_Mul(i):

                    vmm = fpoly(i, 'list')
                    if kval in vmm:
                        newm = simplify(i / kval)
                        mm1 += newm
                    else:
                        mm2 += i
            kres2 = kval * mm1 + mm2
            self.update(kres2)
            return kres2

        return kres
    
    def tri2bin(self,sexpr,var=x):
        ksym=self.ksym
        kres=tri2bin(ksym,sexpr,var=var)
        self.ksym=kres
        self.s()
    
    def trinom2binom(self,sexpr,kshow=True):
        ksym=self.ksym
        kres=trinom2binom(ksym,sexpr)
        self.ksym=kres
        if kshow:
            self.s()
        
    
    
    def fixrootpow(self,op='',kshow=True):
        p1=self.ksym
        p1=fixrootpow(p1)
        self.ksym=p1    
        if kshow:
            self.s()

                
    def killroot(self):
        ksym=self.ksym
        nkasym=simplify_root_exp(ksym)
        self.ksym=nkasym
        self.s() 
        
    def fixRootPow(self, kksym):
        self.setValue(rpow(kpow(kksym, 2), 2), kksym)

    def fix_sqrtPow(self):
        ksym=self.ksym
        nkasym=fix_sqrt2pow(ksym)
        self.ksym=nkasym
        self.s()

    def find_root(self, ksym, x1, x2, xx):
        xx = np.linspace(x1, x2, xx)
        yy = self.evalueArray(ksym, xx)
        mm = []
        qq = len(yy)
        for i in range(qq - 1):
            if (yy[i] > 0 and yy[i + 1] < 0) or (yy[i] < 0 and yy[i + 1] > 0):
                mm.append((xx[i] + xx[i + 1]) / 2)
        qq = len(mm)
        for i in range(qq):
            e1 = MyEq(self.v, str(ksym), kshow=False)
            e1.setValue(ksym, mm[i])

    def upBag(self, Bag, kname='', kshow=True, kope=''):
        vs = Bag.vmain
        vv = Bag.vsolve

        for i, j in zip(vs, vv):
            self.set(i, j, kshow=False)
            self.opemat(kope, kshow=False)
        if kshow:
            self.s()

    def upTriang(self, angul, T3, kope=''):
        if type(angul) == MyTriang:
            angu2 = angul
            angul = T3
            T3 = angu

        v1 = [sin(angul), cos(angul), tan(angul), kpow(sin(angul), 2), kpow(cos(angul), 2), kpow(tan(angul), 2)]
        v2 = [T3.sin(), T3.cos(), T3.tan(), kpow(T3.sin(), 2), kpow(T3.cos(), 2), kpow(T3.tan(), 2)]

        for i in range(len(v2)):
            v2[i] = opemat(v2[i], kope=kope)

        self.set(v1, v2)

    def evalueBag(self, bag, kope=''):
        e1 = MyEq(self.v, self.name, kshow=False)
        e1.upBag(bag, kope=kope)

    def setVal_from_bag(self, bag, kshow=False, kope=''):
        for i, j in zip(bag.dataS, bag.dataV):
            self.setValue(i, j, kshow=kshow, kope=kope)
        if kshow:
            self.s()

    #  Algoritmos de reparacion

    def fix_reduc(self):
        kres = self.ksym
        try:
            kres = fix_reduc(kres)
            self.ksym=kres
            self.s()
        except:
                self.s()
    ###############################
    #  Diff
    ###############################
        
    def neginverse(self,kname=''):
        if kname!='':
            return  MyEq(-1/self.ksym,kname=kname)
        else:
            return -1/self.ksym
            
    def doDiff(self,kname='',var=''):
        if kname=='':
            kname= self.name 
        if var=='':
            var=self.var
            
        kres=self.ksym
        kres=diff(kres,var)
        ee=MyDiffe(kres,kname,var=var)
        return ee

        
        
            
    def diffEq(self, kname='', var2='', ktype='P', typeD=1):
        self.oldprimi = self.primi
        kres = self.ksym
        if var2 == '':
            var2 = self.var2
        if typeD == 2:
            return MyEq(Derivative(kres, var2), kname=kname, var2=var2, ktype='Diff')

    def changediffI(self, newvar2, newvalue='', x1='', x2=''):
        kres = self.primi
        if newvalue != '':
            kres = kres * newvalue
            self.primi = kres
        self.var2 = newvar2
        if x1 != '':
            self.x1 = x1
            self.x2 = x2
        if self.x1 == '':
            self.ksym = Integral(self.primi, self.var2)
        else:
            self.ksym = Integral(self.primi, (self.var2, self.x1, self.x2))
        self.s()
        
        
    def changeDiff(self,ee,x1='',x2=''):
        if self.type=='D' and ee.type=='D':
            oldv=self.var
            newv=ee.var
            kres1=self.ksym
            kres2=ee.ksym
            kres=kres1*kres2
            self.ksym=kres
            self.var=newv
            self.diffname=ee.diffname
        
        if self.type=='D' and ee.type=='P':
            oldv=self.var
            newv=ee.var
            kres1=self.ksym
             
            kres2=ee.ksym
            kres1=kres1.subs(oldv,kres2)
            kres3=diff(kres2,newv)
            kres1=kres1*kres3


            self.ksym=kres1
            self.var=newv
            self.diffname=diff_name(self.var)
            
            
        if self.type=='I' and ee.type=='P':
            self.upgrade(ee,kshow=False)
            ksym,oldv,x1,x2=self.get_intedata()
             
            kdiff=diff(ee.ksym,ee.var)
            kres=ksym*kdiff
             
            xx1=self.x1
            xx2=self.x2
            var=ee.var
             
           
            if x1!='':
                xx1=x1
                xx2=x2
            p2=(var,xx1,xx2)
            p1=kres
            self.primi=kres
            ksym= Integral(p1,p2)
            self.ksym=ksym
            self.x1=xx1
            self.x2=xx2
            self.var=var
    
        self.s() 
            
    def set_interval(self,x1,x2):
        self.change_integral(nvar=self.var,x1=x1,x2=x2)
        self.s()
        
    def changediff(self,newv,expr):
        r'''
        input, newv= new variable to changediff
        expr , is   thw value of old main variable 'var' 
        in function the new variable
        '''
        newd=expr.diff(newv)
        kres=self.ksym
         
        
        oldv=self.var
         
        
        if self.type=='I':
            val=kres.args
            eqq=val[0]
            eqq=eqq.subs(oldv,expr)
            eqq=eqq*newd
            # if self.x1!='':
                # x1=val[1][1] 
                # x2=val[1][2] 
                # qq=MyEq(oldv-expr,'qq',kshow=False)
                # valim=qq.solve(newv,kshow=False)
                # nx1=valim.ksym.subs(oldv,x1)
                # nx2=valim.ksym.subs(oldv,x2)

                # ee=MyEq(eqq,kname=self.name,x1=nx1,x2=nx2,var=newv,ktype='I',kshow=False)
                # self.x1=nx1
                # self.x2=nx2
                # self.var=newv
                # self.ksym=ee.ksym
                # self.s()
            # else:
            ee=MyEq(eqq,kname=self.name,var=newv,ktype='I',kshow=False,x1=self.x1,x2=self.x2)

            self.var=newv
            self.ksym=ee.ksym
            self.s()
        else:
             
            kres=kres.subs(oldv,expr)
            kres=kres*newd
            self.ksym=kres
            self.var=newv
            self.s()
        
    ###############################
    #  traasforma to function            
    
    def transform2(self,sexpr,var=x,kshow=True):
        ksym=self.ksym
        kres=transform2(ksym,sexpr,var=var)
        self.ksym=kres
        if kshow:
            self.s()
        



    
    ###############################
    #  inverse function
    def inversefunc(self,var=''):
        if var=='':
            helplib('inversefunc')
            return
        expr=self.ksym
        kres=inversefunc(expr,var)
        oldn=self.name
        nname=oldn+'^{-1}'
        nP=symbols(nname)
        Q=Eq(nP,kres)
        display(Math(latex(Q)))
        return kres
    
        
    
    ###############################
    #  Integral
    #######################    
    
    def solveI(self):
        kres=self.ksym

        kres=kres.doit()
        kres=fix_otherwise(kres)
        rr=ccode(kres)
        return rr

    def clean_ccode(self,kstr):
        kstr=kstr.replace('{','')
        kstr=kstr.replace('}','')
        mlist=kstr.split(',')
        return mlist
        
    def clean_solve(self,kstr):
        vecc=self.clean_ccode(kstr)
        vecr=[]
        for i in vecc:
            try:
                kres=parse_expr(i)
                if kres!=True and kres!=False:
                    vecr.append(kres)
            except:
                pass
        return vecr
    
    def Isolution(self):
        kstr=self.solveI()
        vecr=self.clean_solve(kstr)
        if len(vecr)==1:
            self.ksym=vecr[0]
            self.s()
        else:    
            return vecr
    
    
    def soveIntegral(self):
        self.ksym=self.ksym.doit()
        self.s()
        
        
                
    def doitI(self, kname='', kope='', kshow=True,C1=C1, **kwargs):
            
            kres=self.ksym
            
            if len(kwargs)>0:
                for key, value in kwargs.items():
                    kres=kres.subs(parse_expr(key),value)
                    
            kres2=kres.doit()   
             
            kres = opemat(kres2, kope=kope)
            kres=fix_otherwise(kres)
              
            rr=ccode(kres)
             
            if '\n' in rr:
                print('pason')
                 
                
                vecs=get_vecposstrfind(rr,'\n')
                exp1=rr[vecs[0]+1:vecs[1]]
                 
                if passdoitI(exp1):
                    kres=parse_expr(exp1)
                       
                else:
                    exp2=rr[vecs[len(vecs)-2]+3:vecs[len(vecs)-1]]
                    if passdoitI(exp2):kres=parse_expr(exp2)
                        

                 
                
           
                        
            if self.name == 'Vo':
                if sign(kres) == -1:
                    kres = -1 * kres
            try:
                kres=kres+C1 
            except:
                pass
            
            if kname != '':
                ee = MyEq(kres, kname=kname,kshow=kshow)
                return ee
            else:
                self.ksym=kres
                if kshow:    
                    self.s()

    def clearlimit(self,name=''):
        if self.type=='I':
            expr=self.primi
            var=self.var
            kres=Integral(expr,var)
            if name=='':
                return kres
            else:
                ee=MyIntg(expr,name,var=var)
                return ee
            
    def doit(self, *args,kshow=True, c1=0, C=0,C1=C1,keyAPI='',nolimit=False,**kwargs):
        vecpara=['area','positive']
        for i in args:
            if i not in vecpara:
                kname=i
                
               
        if self.type == 'diff':
            if len(kwargs) == 0:
                kres = self.ksym
                var2 = self.var2
                kres2 = diff(kres, var2)
                self.update(kres2)
                self.type = 'F'
                s1 = self.name
                s2 = alphaname(var2)
                kname = s1 + '_{(' + s2 + ')}'
                self.name = kname
                self.s()
                return

        if self.type == 'I':
            if keyAPI!='':
                sexpr=str(self.ksym)
                sexpr=sexpr.replace('Integral(','Integrate[')
                sexpr=sexpr.replace('), (','), {')
                sexpr=sexpr.replace(sexpr[-2::],'}]')
                import wolframalpha
                client = wolframalpha.Client(keyAPI)
                q = sexpr
                res = client.query(q)
                answer = next(res.results).text                
                p1=answer.find('=')+2
                p2=answer.find('≈')
                if p2!=-1:
                    sres=answer[p1:p2]
                else:
                    sres=answer[p1::]
                self.ksym=parse_expr(sres,evaluate=False)
                self.s()    
            
            else:
                
                kres=self.ksym
                if self.x1!='':
                    c1=0
                    C1=0
                    c1=0
                if len(kwargs) > 0:
                     
                    kres = real_subs(kres,**kwargs)
                      

                 
                kres = kres.doit()
                if 'positive' in args or 'area' in args:
                    if signo(kres)==-1:
                        kres=-1*kres
                if kname != '':
                    ee = MyEq(unisymbols(kres), kname=kname,var=self.var)
                    ee.type='I'
                    return ee
                else:
                    return unisymbols(kres)

                    
        elif self.type == 'Diff':
            kres = self.ksym
            kres = kres.diff(self.var2)
            self.ksym=kres
            self.type = 'P'
            self.name = 'd' + self.name
            self.s()


        else:
            try:
                ksym=self.ksym
                ksym=ksym.doit()
                self.ksym=ksym
                self.s()
            except:
                self.s()
         

    def get_diff(self, var2='', kname=''):
        if var2 == '' and kname == '':
            var2 = self.var2
        if var2 != '':
            if type(var2) == str:
                kname = var2
                var2 = self.var2
                if type(kname) != str:
                    var2 = kname

        kres = self.ksym
        kres = diff(kres, var2)
        if kname != '':
            ee = MyEq(kres, kname, var2=var2)
            return ee
        else:
            return kres
    def get_primitive(self):
        return self.primi_eq
    def change_primitive(self,expr):
        QQ=MyIntg(expr,self.name,var=self.var,x1=self.x1,x2=self.x2)
        self.ksym=QQ.ksym
        
    
    def get_dprimitive(self):
        return self.diff_eq
        
         
    def get_primitive(self,kname=''):
        if kname!='':
            ee=MyEq(self.primi_eq ,kname=kname,var=self.var)
            return ee
        else:    
            return self.primi_eq        
        
    def kdiff(self, kname='', var2='', kope='', kupdate=False, ktype=''):
        if ktype == '':
            ktype = self.type
        kres = self.ksym
        if var2 == '':
            var2 = self.var2
        kres = self.ksym
        kres2 = diff(kres, var2)
        return kres2
        if kname != '':
            ee = MyEq(kres2, kname=kname, var2=var2, ktype=ktype)
            return ee
        else:
            return kres2

    def update_inte(self):
        if self.type == 'I':
            kres = self.primi
            if self.x1 == '':
                self.ksym = Integral(kres, self.var2)
            else:
                self.ksym = Integral(kres, (self.var2, self.x1, self.x2))

    def fac_integral(self):
        ksym = self.primi
        var2 = self.var2

        if Is_Mono(ksym) and Is_Mul(ksym):
            kres = [x for x in ksym.args if str(var2) not in str(x)]
            kres2 = [x for x in ksym.args if str(var2) in str(x)]
            mono1 = 1
            mono2 = 1
            for i in kres:
                mono1 = mono1 * i
            for i in kres2:
                mono2 = mono2 * i
            self.primi = mono2
            kfac = mono1
            self.primi = mono2
            if self.x1 == '':
                self.ksym = Integral(mono2, self.var2)
            else:
                self.ksym = Integral(mono2, (self.var2, self.x1, self.x2))
            self.Mul(mono1)
        else:
            self.s()

    def maximun(self, ksym, ksave=True):
        '''
        Return maxumin value of ksym if self=0
        '''
        if ksave:
            savework()
        kres = self.ksym
        kres = diff(kres, ksym)
        kres = solve(kres, ksym)
        kname = ksym.name
        if len(kres) == 1:
            return MyEq(kres[0], kname)
        else:
            vres = []
            cc = 1
            for i in kres:
                nname = kname + str(cc)
                vres.append(MyEq(i, kname=nname))
                cc += 1
        return vres

    def plot(self,*args,ymax='',ymin=''):
        x1=-1
        x2=1 
        x3=100
        if len(args)==2:
            x1=args[0]
            x2=args[1]
        elif len(args)==3:
            x1=args[0]
            x2=args[1]
            x3=args[2]
        expr=self.ksym
        var =self.var       
        Pplot(expr,var,x1=x1,x2=x2,x3=x3,ymax=ymax,ymin=ymin)
        

        
    def kplot(self, ksym, x1, x2, x3):
        x = np.linspace(x1, x2, x3)
        y = np.array([float(self(x=x)) for x in xx])
        data_plot = pd.DataFrame({str(self.var2):x, self.name:y})
        sns.lineplot(x = str(self.var2), y = self.name, data=data_plot)
        plt.show()
        
        
    def length_arc(self, x1='', x2='', x='', ksolve=True):
        if x1 == '':
            x1 = self.x1
            x2 = self.x2
        if x == '':
            x = self.var2
        y = MyEq(self.ksym, 'y', x1=x1, x2=x2, varx=x, ktype='Diff')
        y.doit()
        if ksolve == 'dL':
            y.name = 'dL'
            y.s()
            return y
        y2 = y * y
        y2 = opemat(y2, 'ef')
        y3 = rpow(1 + y2)
        y3 = opemat(y3, 'r')
        L = MyEq(y3, 'L', x1=x1, x2=x2, varx=x, ktype='Diff')

        L.integral()
        if ksolve:
            return L.doitI()
        else:
            return L


    ################################################
    ##                Diferencial                ###
    ################################################
    def limit(self,*args):
        if len(args)==0:
            helplib('limit')
            return
        ksym=self.ksym
        var=self.var
        op=""
        valr=args[0]
        if len(args)==2:
            op=args[1]
            return limit(ksym, var, valr, dir=op)
        else:
            return limit(ksym, var, valr)
    
    def Limit(self,*args):
        if len(args)==0:
            helplib('limit')
            return
        ksym=self.ksym
        var=self.var
        op=""
        valr=args[0]
        if len(args)==2:
            op=args[1]
            return Limit(ksym, var, valr, dir=op)
        else:
            return Limit(ksym, var, valr)
        self.type='L'
        self.varL=valr

    
    def dsolve(self, kname='', C1='', C2=''):
        kres = self.EqDiff
        dres = dsolve(kres)
        mdres = dres.args
        kres = mdres[1]
        if C1 != '':
            kres = kres.subs('C1', C1)
        if C2 != '':
            kres = kres.subs('C2', C2)
        if kname != '':
            ee = MyEq(kres, kname=kname, var2=self.var2, ktype='F')
            return ee
        else:
            ee = MyEq(kres, kname=alphaname(self.var1), var2=self.var2, ktype='F')
            return ee
    def derivada(self,*args):
        kname='d'+self.name
        kres=diff(self.ksym,self.var)
        return MyDiffe(kres,kname,var=self.var)
        
    def doDiff(self):
    
        kname='d'+alphaname(self.name)
        var=self.var
         
         
        kres=self.ksym
        kres=diff(kres,var)
        ee=MyEq(kres,kname,var=var,ktype='D')
        return ee

        
    def Diff(self,*args,kshow=True):
         
        kres=self.ksym
        var=self.var
        kname=''
        for i in args:
            if type(i)==str:                 
                kname=i
            if Is_Symbol(i):                 
                var=i
        kres=diff(kres,var)
        if kname=='':
            kname='d_'+self.name
        return MyEq(kres,kname=kname,var=var,ktype='D',kshow=kshow)
        
        
        
        

    def applyDiff(self, var2=''):
        kres = self.ksym
        if var2 == '':
            var2 = self.var2
        kres = kres.diff(var2)
        self.ksym=kres
        self.s()

    def fab_diff(self):
        var2=self.var2
        kname= self.name
        f1,f2,f3=[Function('f1')(var2),Function('f2')(var2),Function('f3')(var2)]
        nvarf=['f1','f2','f3']
        vvarf=[f1,f2,f3]
        kres= self.ksym
        ee2= self.xcopy('ee2',kshow=False)
        qq=len( self.varf)
        vvarf=vvarf[0:qq]
        nvarf=nvarf[0:qq]
        varf= self.varf
        for i,j in zip(varf,vvarf):
            ee2.set(i,j,kshow=False)
        kres=ee2.ksym
        kres=kres.diff(var2)
        ee3=MyEq(kres,'ee3',kshow=False)
        ee3.primi_eq=self.primi_eq
         
        for i,j in zip(vvarf,varf):
            ee3.setdiff(i,1,kshow=False)
            ee3.set(i,j,kshow=False)
        ee3.name= diffname(self.var1,var2)
        ee3.primi_diff_eq=diff(ee3.primi_eq,ee3.var2)

        return ee3
        
    def Diff2diff(self,kshow=True):
        fksym=self.flatdiff
        for i in self.varf:
            kvar=str(i)
            kvar2='d'+kvar
            var2=self.var2
            kres=self.primitiva
            f=Function(i)(var2)
            news=symbols(alphaname(kvar2))
            fksym=fksym.subs(diff(f,var2),news)
        return(fksym)    

        
    def Diff2diff(self,kres): # ksym,kvar,var2
        var2=self.var2
        kvar=self.varf
        for i in kvar:
            f=Function(str(i))(var2)
            df=diff(f)
            kname='d'+alphaname(i)
            nf=symbols(kname)
            kres=kres.subs(df,nf)
        return kres
    def p(self):
        return self.primitiva
        
    def pf(self):
        try:
            return primi2func(self.ksym,self.var2,self.varf)
        except:
            return self.ksym
        
    def Fdiff(self,*args):
        kname=''
        kupdate=False
        for i in args:
            if type(i)==str:
                if i=='update':
                    kupdate=True
                else:
                    kname=i
        if self.varf==[]:
            kres=diff(self.ksym,self.var2)
            primi=0
        else:
            primi=diff(self.primitiva,self.var2)
            kres=Diff2diff(primi,self.varf,self.var2)
            
        if kname=='':
            if kupdate:
                self.ksym=kres
                self.primitiva=primi
                
                self.s()
                 
            else:
                
                return kres
        else:
            ee=self.xcopy(kname,kshow=False)
            ee.ksym=kres
            ee.primitiva=primi
            ee.dtype='diff'
            ee.s()
            return ee
            
    def diffQ(self,*args,niceview=True):
        ''' 
        
        
        input: diff(t,x,y,kname='')
        return diff(Q, in t and  x=func)     
        
        
        
        input: diff((t,x,y),kname='')
        return diff(Q, in t and  x,y =func in t) 
        
        args:
            if kname=='', self.update, 
            kname!='' return new Eq named kname
            
            niceview=True
            return example Diferential(x(t),t) = dx/dt but  value is original
        '''    
        exp1=self.ksym
         
        var=args[0]
        var1=args[1]
        kname=''
        var2=''
        if len(args)==3:
            kval=args[2]
            if type(kval)==str:
                kname=kval
                 
            else:
                var2=kval
            
        exp22=diffuntion(self.ksym,var,var1,var2)
        if kname!='':
            ee=MyEq(exp22,kname=kname,kshow=False)
            if niceview:
                ee2=MyEq(viewnicediff(exp22,var,var1,var2),kname)
                return ee    
        else:
            self.ksym=exp22
            if niceview:
                ee2=MyEq(viewnicediff(exp22,var,var1,var2),self.name)
                    
    def diff2(self,*args,kshow=True):
        kres=self.ksym
        var=self.var
        kres=diff(kres,var,var)
        if len(args)==0:
            return kres
        else:     
            if type(args[0])==str:
                return MyEq(kres,kname=args[0],var=var,kshow=kshow)
            else:
                return kres.subs(var,args[0])
    
    def limitdiff(self):
        var=self.var
        sx='d'+str(var)
        dx =symbols(sx) 
        kres=(self(var+dx)-self(var))/dx
        ee=MyEq(kres,kname='d'+self.name ,var=dx)
        return ee

    def simplifysum():
        kres=0
        ksym=self.ksym 
        if Is_Add(ksym):
            for i in fpoly(ksym,'list'):
                kres=kres+simplify(i)
            self.ksym=kres 
        else:
            self.ksym=simplify(ksym)    
                  
        self.s() 
    
    def factorsum():
        kres=0
        ksym=self.ksym 
        if Is_Add(ksym):
            for i in fpoly(ksym,'list'):
                kres=kres+factor(i)
            self.ksym=kres 
        else:
            self.ksym=factor(ksym)    
                  
        self.s()

                  
        self.s()
        
    def diff(self,*args,func='',kshow=True):
        kname=''
        var=self.var
        vecvar=[]
        for i in args:
            if type(i)==Symbol:
                vecvar.append(i)
            if type(i)==str:
                kname=i
                
        if func!='':
            dF=functiondiff(self,var,func)
            dF=sympydiff2prime(dF,var,func)
        elif len(vecvar)==0:
            dF=diff(self.ksym,var)
        else:
            dF=self.ksym
            for i in vecvar:
                dF=diff(dF,i)
                
        if kname!='':
            ee=MyEq(dF,kname=kname,var=var,kshow=kshow)
            return ee
        else:
            return dF
            
                    
        
            
    def diffFunction(self,var2):
        var=self.var
        expr=functiondiff(self,var,var2)
        kname='d'+self.name
        kres=sympydiff2prime(expr,var,var2)
        ee=MyEq(expr,kname=kname,var=var,ktype='dP',kshow=False)
        ee.primi_diff_eq=kres
        ee.s()
        return ee   
        
    def tanEq(self,px=''):
        var=self.var
        if px=='':
            svar=str(self.var)+'1'
            px=symbols(svar)
        df=MyEq(self.diff(),'ee',var=self.var,kshow=False)
         
        eT=df(px)*(self.var-px)+self(px)
        return eT

 
    def eQtangent(self,px=''):
        return self.tanEq(px)
    
    
    
    def eQnormal(self,px=''):
        var=self.var
        if px=='':
            svar=str(self.var)+'1'
            px=symbols(svar)
        df=MyEq(self.diff(),'ee',var=self.var,kshow=False)
         
        eT=-1*(self.var-px)/df(px)+self(px)
        return eT 
            
    def diffk(self, *args, respect='',kupdate=False,kshow=True,**kwargs):
               
        if self.varf==[]:
              self.primi_eq=self.ksym
              self.diff_eq=diff(self.ksym,self.var)
              self.primi_diff_eq=diff(self.ksym,self.var)
            
        if self.varf!=[]:


            if self.primi_eq=='':
             self.primi_eq=self.get_primitive()

            if self.primi_diff_eq=='':    
             ee=self.xcopy('ee',kshow=False)
             kres=ee.primi_eq
             kres=kres.diff(ee.var2)
             self.primi_diff_eq=kres
         

            
        
        if len(args)==0 and len(kwargs)==0:
            self.dtype='diff'
            kname=self.name
            ee3= self.fab_diff()
            self.ksym=ee3.ksym
            self.name=diffname(self.name,self.var2)
            if respect!='':
                self.name= diffname(self.name,respect)
            if self.varf!=[]:
                self.dtype='diffd'

            self.diff_eq=self.ksym

            self.s()

            return

        elif len(args)==1 and len(kwargs)==0 and type(args[0])==str and self.varf!=[]:
            kname=self.name
            ee3= self.fab_diff()
            ee3.dtype='diff'
            ee3.name= diffname(kname,ee3.var2)
            if respect!='':
                ee3.name= diffname(kname,respect)
                ee3.diffEq=self.diffEq
                 
            if ee3.varf!=[]:
                 ee3.dtype='diffd'
            
                 
            ee3.s()
            ee3.origen=self.ksym
            ee3.origen=self.ksym
            ee3.vfunc=self.vfunc
            return ee3
        else:
            kname = ''
            var = self.var2
            kres = self.ksym
            for i in args:
                if type(i) == str:
                    kname = i
                else:
                    var = i
            kres = diff(kres, var)
            if len(kwargs) > 0:
                ee = MyEq(kres, kshow=False)
                for key, value in kwargs.items():
                    ee.set(parse_expr(key), value, kshow=False)
                kres = ee.ksym

            if kname != '':
                ee = MyEq(kres, kname,var1=self.var1,kshow=kshow)
                ee.dtype='diff'
                ee.var2=self.var2
                ee.vfunc=self.vfunc
                return ee
            else:
                return kres

    def simple_form(self,kfunc):
        self.setdiff(kfunc,1,kshow=False)
        self.set(kfunc, self.var1)

    def flatdiff(self,kfunc,kshow=True):
        self.setdiff(kfunc,1,kshow=False)
        self.set(kfunc, self.var1,kshow=kshow)

    def diff_parcial(self,*args):
        kname=''
        vdiff=[]
        for i in args:
            if type(i)==str:
                kname=i
            else:
                vdiff.append(i)


        kres=self.ksym
        var2=self.var2
        mm=vdiff
        mm=list(mm)
        mm.append(var2)
        mm=tuple(mm)
        ddiff=sym2Function(*mm)
        for i,j in zip(vdiff,ddiff):
            kres=kres.subs(i,j)

        kres=diff(kres,var2)
        if kname!='':
            return MyEq(kres,kname,var2=var2)
        else:
            self.ksym=kres
            self.s()

    def killdiff(self,kfunc):

        self.setdiff(kfunc,1)

    def killfunc(self,kfunc):

        self.set(kfunc, self.var1)

    def killfunc(self,kfunc):

        self.set(kfunc, self.var1)

    def Length(self, kvar):
        kres = self.v
        kdiff = kres.diff(kvar)
        kres = rpow(1 + kpow(kdiff, 2))
        return kres

    #  Algoritmos Trigonometricos
    def sin2cos(self, angu, korden=2, kope='', kshow=True):
        self.set(kpow(sin(angu), 3), (1 - kpow(cos(angu), 2)) * cos(angu), kshow=False)
        kres = self.ksym
        kres = sin2cos(kres, angu=angu, korden=korden, kope=kope)
        self.ksym=kres
        if kshow:
            self.s()
        else:
            pass

    def cos2sin(self, angu, korden=2, kope='', kshow=True):
        self.set(kpow(cos(angu), 3), (1 - kpow(sin(angu), 2)) * sin(angu), kshow=False)
        kres = self.ksym
        kres = cos2sin(kres, angu=angu, korden=korden, kope=kope)
        self.ksym=kres
        if kshow:
            self.s()

    def get_diffEq(self, vx=''):

        kname = self.name + 'd'
        Vx = self.xx

        if Vx != '':

            kres = self.ksym

            kres = diff(kres, Vx)
            kname = kname + '_' + str(Vx)
            self.EqDiff = MyDiff(kres, Vx, kname)
            return self.EqDiff
        elif vx != '':
            self.xx = vx
            kres = self.ksym
            kname = kname + '_' + str(vx)
            kres = diff(kres, vx)
            self.EqDiff = MyDiff(kres, vx, kname)
            return self.EqDiff
        else:
            sE(['first set variable MyEq.xx= Valeue'])

    def Integral_eQ(self, kname='', var2='', x1='', x2='', ktype='P'):
        ksym = self.ksym

        if var2 == '':
            var2 == self.var2
        if x1 == '':
            x1 = 0
            x2 = var2
        else:
            x1 = x1
            x2 = x2
        kres = Integral(ksym, (var2, x1, x2))
        if kname == '':
            return kres.doit()
        else:
            ee = MyEq(kres, kname=kname, var2=var2, ktype='I')
            return ee

    def Area(self, kname='', var2='', x1='', x2='', kope=''):
        if var2 == '':
            var2 = self.var2
        if x1 == '':
            x1 = self.x1
            x2 = self.x2

        ee = MyEq(self.ksym, kname=kname, var2=var2, x1=x1, x2=x2, ktype='I', kshow=False)
        if kname != '':
            ee.s()
            ee.doitI()
            return ee
        else:
            ee.doitI(kshow=False)
            return ee.ksym
    def func(self,value):
        ksym=self.ksym
        var=self.var
        kres=ksym.subs(var,value)
        try:
            return float(kres)
        except:
            return kres
            
    def area_value(self):
        ee = MyEq(self.ksym, var2=self.var2, x1=self.x1, x2=self.x2, kshow=False)
        ee.area(kshow=False)
        ee.doit(kshow=False)
        return ee.ksym

    def area(self, var2='', x1='', x2='', kshow=True):
        if x1 != '':
            self.x1 = x1
            self.x2 = x2
        if var2 == '':
            var2 = self.var2
        else:
            self.var2 = var2
        self.primi = self.ksym

        if self.x1 == '':
            kres = Integral(self.ksym, self.var2)
        else:
            kres = Integral(self.ksym, (self.var2, self.x1, self.x2))
        self.type = 'I'
        self.ksym=kres
        if kshow:
            self.s()
 
    def Integral(self,x1='', x2=''):
        ktype=self.type
         
        if x1=='':
            if self.x1!='':
                x1=self.x1
                x2=self.x2
        else:
             
            self.x1=x1
            self.x2=x2
        
        if ktype=='I':
            self.s()
            return
        else:
            kres=self.ksym
            if x1=='':
                kres=Integral(kres,self.var)
            else:
                kres=Integral(kres,(self.var,x1,x2))
            self.ksym=kres
            self.type='I'
            self.s()
            return
         
        
            
            
        

        
        

    def solveIntegral(self, kname='', var2='', x1='', x2='', kupdate=True):
        if x1 == '':
            x1 = self.x1
        if x2 == '':
            x1 = self.x2
        if var2 == '':
            var2 = self.var2
        ksym = self.ksym
        kres = ksym.doit()
        if kname != '':
            ee = MyEq(kres, kname=kname, var2=var2)
            return ee
        else:
            if kupdate:
                self.ksym = kres
                self.s()
            else:
                return kres

    def get_InteEq(self, vx='', kname=''):

        if kname == '':
            kname = 's' + self.name
        else:
            kname = kname
        Vx = self.xx

        if Vx != '':

            kres = self.ksym

            self.EqInte = MyInteger(kres, Vx, kname)
            return self.EqInte
        elif vx != '':
            self.xx = vx
            kres = self.ksym

            self.EqInte = MyInteger(kres, vx, kname)
            return self.EqInte
        else:
            sE(['first set variable MyEq.xx= Value'])

    def convMyFunc(self, knom='', vV=''):
        ksym = self.v
        if vV == '':
            vV = self.free()
        elif type(vV) == tuple:
            vV = list(vV)
        ksym = self.v
        kk = MyFunc(knom, ksym, vV)
        return kk

    def newEq(self, kname='', kope='', **kwargs):
        ee = MyEq(self.ksym, kname=kname, kshow=False)
        if len(kwargs) > 0:
            for key, value in kwargs.items():
                ee.set(parse_expr(key), value, kshow=False, kope=kope)

        ee.s()
        return ee

    def findSubFunc(self, sval, inside=''):
        ksym = self.ksym
        kini = 0
        kini2 = 0
        sroot = []
        done = 0
        while kini < len(str(ksym)) and kini2 != -1:
            kini2, sword = in_ope_string(ksym, sval, kini)
            if kini2 != -1:
                if inside != '':
                    if inside in sword:
                        return sword
                else:
                    return sword
            kini = kini2 + len(sval)
        return sroot

    def killSimpleRoot(self):
        ss = self.findSubFunc('sqrt', '**2')
        sm = get_midle_str(ss, 'sqrt', '**2')
        exp = str(self.ksym)
        exp = exp.replace(ss, sm)
        self.ksym = parse_expr(exp)
        self.s()
    ################################################################################
    #               alone
    #################################################################################

    def alone(self,ksym):
        
        newv=symbols(self.name)
        ee=MyEq(newv-self.ksym)
        kres=ee.solve(str(ksym),kshow=False)
        kres.s()
        return kres

    def maximize(self,respp=''):
        if respp=='':
            respp=self.var2 
            
        dv=MyEq(self.diff(respp),kshow=False)
        dv.reducecero(kshow=False)
        kres=dv.ksym
        kres=numer(kres)
        dv.ksym=kres
        ksolu=dv.solve(respp)
        return ksolu
        
    #  Fourier Serie

    def fourierserie(self,n=''):
        s=fourier_series(self.ksym,(self.var,self.x1,self.x2))
        kres=s
        if n!='':
            kres=s.truncate(n=n)
        self.furs=s
        return kres
    def L(self):
         return self.x2-self.x1
    
    def an(self,n=''):
        if self.furs=='':
            s=fourier_series(self.ksym,(self.var,self.x1,self.x2))
            self.furs=s
            if n!='':
                return s.an[n]
            else:
                return s.an
    
    def bn(self,n=''):
        if self.furs=='':
            s=fourier_series(self.ksym,(self.var,self.x1,self.x2))
            self.furs=s
            if n!='':
                return s.bn[n]
            else:
                return s.bn  
    #               LAPLACE
    
    def laplace(self):
        from MyDiff import laplacetransform, func2laplace,laplace, Laplace
        return laplace(self)
     
            
    
    def Laplace(self):
        from MyDiff import laplacetransform, func2laplace,laplace, Laplace
        return Laplace(self) 

    def ilaplace(self,*args):
        kname='Ft'
        var=t
        for i in args:
            if type(i)==str:
                kname=i
            else:
                var=i
        from MyDiff import laplacetransform, func2laplace,laplace, Laplace,ilaplace
        kres= ilaplace(self.ksym)
        if var!=t:
            kres=kres.subs(t,var)
        return MyEq(kres,kname,var=var)

    
################################################################################
#               END MyEq Class
#################################################################################
# def MyIntg(ksym,kname,var=x,x1='',x2='',kshow=True,ktype='I',value=False):
    # if type(ksym)==MyEq:
        # ksym=ksym.ksym
    
    # ssym=str(ksym)
    # dsym='d'+str(var)
    # if dsym in ssym:
        # ssym=ssym.replace(dsym,'1')
        # ksym=parse_expr(ssym)
    # if value:
        # ee=MyEq(ksym,kname,var=var,x1=x1,x2=x2,ktype='I',kshow=False)
        # return ee.ksym
    # else:    
        # return MyEq(ksym,kname,var=var,x1=x1,x2=x2,ktype='I',kshow=kshow)
    
def MyDiffe(ksym,kname,var=t,ktype='D'):
    ssym=str(ksym)
    dsym='d'+str(var)
    if dsym in ssym:
        ssym=ssym.replace(dsym,'1')
        ksym=parse_expr(ssym)
     
    return MyEq(ksym,kname,var=var,ktype='D')
    
def MyDiff(ksym,kname,var=t,ktype='D'):
    ssym=str(ksym)
    dsym='d'+str(var)
    if dsym in ssym:
        ssym=ssym.replace(dsym,'1')
        ksym=parse_expr(ssym)
     
    return MyEq(ksym,kname,var=var,ktype='D')

def get_real_value(ee):
    if type(ee) == MyEq:
        kres = ee.ksym
    else:
        kres = ee
    return kres


def show_main():
    for i in dataQ:
        i.s()



def upBagSys(ksys, Bag, kope=''):
    for i in ksys:
        i.upBag(Bag, kope=kope)



###########################################
#               END MyIteger Class
###########################################






###############################################
#  Mass center Inertia

def pQ(mm, vv, kope=''):
    rho = symbols('rho')

    kres = mm / vv
    sE([rho, '=', kres])
    return kres


#################################################
#   Solve Algorithm

def solved(*args,**kwargs):
    '''
    solved (var,exp1,exp2,**kwargs)
        input 
            var : variable to find
            exp1 : math expre or MyEq class that is  equation equal= 0
            exp2 (optional): math expre or MyEq class  
            if exp2 is given then the Eq to evalue is expr1 -expr2
            kwargs : conditions to evalue, example x=1, t=0..etc
        return
            return a MyEq class with name str(var)
    '''        
    var=args[0]
    ee=args[1]
    if type(ee)==MyEq:
        ee=ee.ksym
    keq=ee    
    if len(args)==3:
        ee2=args[2]
        if type(ee2)==MyEq:
            ee2=ee2.ksym
        keq=ee-ee2
    if len(kwargs)>0:             
            keq=real_subs(keq,**kwargs) 
              
    kres=solve(keq,var)
    if type(kres)==list:
        kres=kres[0]
                
    kname=str(var)
    ee0=MyEq(kres,kname=kname)
    return ee0       
    
def solverSys(*args, Bag=''):
    Ke = []
    Kv = []
    Kn = []

    for i in args:
        if type(i) == MyEq:
            if Bag != '':
                i.upBag(Bag, kshow=False)
            Ke.append(i)
        if type(i) == Symbol:
            Kv.append(i)
            Kn.append(i.name)
    # return(Ke,Kv,Kn)

    return MyEqSolveLin(Ke, Kv, Kn, Bag=Bag)


def MyEqSolveLin(Ke, Kv, Kn, Bag=''):  # Solve n MyEq with n unknow variable
    '''
    Example
        Ke=[e2,e2,e0]  MyEqs Matrix
        Kv=[N1,T,a]    unKnow Vriables
        Kn=['N_1','T','a_c']  New Name Eq

        N11,T1,ac = MyEqSolveLin(Ke,Kv,Kn)
        returns resepective  answer
    '''
    vecs = []
    qq = len(Ke)
    kres = []
    for i in range(qq):
        ee = Ke[i]
        ksym = Kv[i]
        ks = ee.solve(ksym,kshow=False)
        if type(ks) == list:
            rr = max(ks)
            ks = rr

        vecs.append(ks)
        Ker = Ke[i + 1::]
        for e1 in Ker:
            e1.set(ksym, ks, kshow=False)
            e1.reduFac(kshow=False)
            e1.simplify(kshow=False)

    for i, kname in zip(vecs, Kn):
        ee = MyEq(i, kname, kshow=False)
        kres.append(ee)
    ueq = kres[-1]
    ksym = ueq()
    vsym = Kv[-1]
    for ee in kres[0:-1]:
        ee.set(vsym, ksym, kshow=False)
        ee.reduFac(kshow=False)
        ee.simplify(kshow=False)
    for i in kres:
        i.s()
    return kres


def Solve2Eq(ksym=[], kvar=[], knom=[], kope=''):
    e1, e2 = ksym
    v1, v2 = kvar
    t1, t2 = knom

    r1 = e1.solve(v1)
    e2.set(v1, r1, kshow=False)
    r2 = e2.solve(v2)
    r2 = opemat(r2, kope=kope)
    e1.set(v2, r2, kshow=False)
    r1 = e1.solve(v1)
    r1 = opemat(r1, kope=kope)
    aa = MyEq(r1, t1)
    bb = MyEq(r2, t2)
    return (aa, bb)


def Diff(ksym, kvar, kname=''):
    kres = ksym
    kres = kres.diff(kvar)
    if kname == '':
        return kres
    else:
        return MyEq(kres, kname)


def Diff2(ksym, kvar, kname=''):
    kres = ksym
    kres = kres.diff(kvar)
    kres = kres.diff(kvar)
    if kname == '':
        return kres
    else:
        return MyEq(kres, kname)



def upBag2sys(vecEq, kBag):
    for i in vecEq:
        i.upBag(kBag)




def eQSolver(*args):
    vec1 = []
    uk1 = []
    for i in args:
        if type(i) == list:
            for j in i:
                if type(j) == MyEq:
                    vec1.append(j())
                elif fpoly(j, 'n') > 1:
                    vec1.append(j)
                else:
                    uk1.append(j)
        else:
            if type(i) == MyEq:
                vec1.append(i())
            elif fpoly(i, 'n') > 1:
                vec1.append(i)
            else:
                uk1.append(i)

    vec2 = []
    kres = []
    for i in vec1:
        if type(i) == MyEq:
            vec2.append(i())
        else:
            vec2.append(i)

    mm = solve(vec2, uk1)
    if type(mm) == dict:
        kk, vv = kunpakDic(mm)

        for i, j in zip(kk, vv):
            kres.append(MyEq(j, i))
        return kres
    else:
        for i, j in zip(mm[0], uk1):
            j = MyEq(i, str(j))
            kres.append(j)
        return (kres)


def solvelin(*args, kope='', Eq=True):  # solveLinearSys(e1,e2,mu1,mu2)
    mS = []
    mV = []

    for i in args:
        if type(i) == MyEq:
            mS.append(i())
         
        elif type(i) == str:
            kope = i
        else:
            mV.append(i)
    solu = solve(mS, mV)
 
    display(Math(latex(solu)))
    return ganswer(solu,'value')
    
    
    
    


def get_squareMono(ksym):
    if type(ksym) == MyEq:
        ksym = ksym.ksym
    kres = ksym
    mm = fpoly(ksym, 'list')
    mr = []
    ms = []
    rr = []
    centra = 0
    ksigno = 1
    for i in mm:
        mr.append(opemat(rpow(i, 2), 'r'))
        ms.append(str(opemat(rpow(i, 2), 'r')))
    for i, j, k in zip(ms, mr, mm):
        if 'sqrt' in i:
            central = k
            if '-' in str(central):
                ksigno = -1
        else:
            rr.append(j)
    if len(rr) == 2:
        kres = kpow(rr[1] + ksigno * rr[0], 2)
    return kres


#######
def expand2MyEq(ee):
    ktype = ee.type
    var2 = ee.var2
    mm = ee.list()
    cc = 1
    kname = ee.name
    kres = []
    for i in mm:
        nname = kname + str(cc)
        nname = MyEq(i, nname, var2=var2)
        kres.append(nname)
        cc += 1
    return kres


def upgrade(*args, kshow=True, andsolve=[]):
    if andsolve != []:
        vv = andsolve[0]
        ee = andsolve[1]
        vv = ee.solve(parse_expr(vv), vv, kope='s')
    eev = []
    evv = []
    for i in args:
        if type(i) == MyEq:
            if i.type == 'C':
                eev.append(i)
            else:
                evv.append(i)

    for i in eev:
        for j in evv:
            try:
                i.upgrade(j, kshow=False)
            except:
                pass
    for i in eev:
        if i.ksym != 0:
            i.s()
    for i in evv:
        if type(i) == MyEq:
            i.simplify(kshow=False)
            if i.ksym != 0:
                i.s()


def upgradeList(*args, kshow=True, andsolve=[], kope='s'):
    eev = []
    evv = []
    for i in args:
        if type(i) == MyEq:
            if i.type == 'C':
                if i != andsolve[1]:
                    eev.append(i)

    if andsolve != []:
        vv = andsolve[0]
        ee = andsolve[1]
        kres = ee.solve(vv)
        kres = opemat(kres, 's')
        vv = MyEq(kres, str(vv), ktype='C', kshow=False)
        ee.type = 'P'

    for i in eev:
        if i.type == 'C':
            i.upgrade(vv, kshow=False)
            i.simplify(kshow=False)

    for i in eev:
        if i.ksym != 0:
            i.s()
    vv.s()
    return vv



def func_sig(kf, x1, x2, var=x):
    ee = MyEq(kf, var2=var, kshow=False)
    xx = (x2 - x1) / 2
    return ee(xx)


def get_intersec_2func(y1, y2, var=x):  # y1(x), y2(x), return intersec y1 and y2
    ee = MyEq(y1 - y2, kshow=False)  # return vector
    return ee.solve(var)


def reset_ee(*args):
    eeFull = []
    for i in args:
        i.init = False


def Upgrade(*args, kope='', kshow=True):
    newa = []
    for i in args:
        if type(i) == MyEq:
             
            if i.ksym != 0:
                newa.append(i)
    args = newa
    antes = []
    for i in args:
        antes.append(str(i))
    qq = len(args)
    for i in range(qq):

        mm = []
        for j in range(qq):
            if j != i:
                mm.append(args[j])
        args[i].upgrade(mm, kshow=False, kope=kope)
    if kshow:
        for i, j in zip(args, antes):
            if str(i) != newa:
                if i.ksym != 0:
                    i.s()


def presolve(ksym, val):
    kres = solve(ksym, val)
    if kres == []:
        try:
            kres = solve(opemat(ksym, 'esf'), val)
            if kres != []:
                return kres
            else:
                ksym = factorSec(ksym, val)
                kres = solve(opemat(ksym, 'esf'), val)
                if kres != []:
                    return kres
        except:
            done = False
    return kres


def eQsolve(ksym, kname, kope=''):
    kval = parse_expr(kname)
    kres = csolve(ksym, kval)
    kres = opemat(kres, kope)
    kval = MyEq(kres, kname)
    return kval
    
def Qsolve(*args):
    '''
        N1,mu1=Qsolve(FxA,FyA,N1,mu1)
    '''    
    eqq=[]
    evv=[]
    for i in args:
        if type(i)==MyEq:
            eqq.append(i)
        else:
            evv.append(i)
    kres=solve(eqq,evv)
    try:
        vsym,vval=kunpakDic(kres)
    except:
        vsym=evv
        vval=list(kres[0])
         
    vres=[]
    for i ,j in zip(vsym,vval):
        kname=str(i)
        ee=MyEq(j,kname=i)
        vres.append(ee)
    return vres   



def Diff2flat(kres,kvar,var2): # ksym,kvar,var2
    
    for i in kvar:
        f=Function(str(i))(var2)
        df=diff(f)
        kname='d'+alphaname(i)
        nf=symbols(kname)
        kres=kres.subs(df,nf)
    ee=MyEq(kres,kshow=False)
    for i in kvar:
        ee.setdiff(i,i,kshow=False)
        
    return ee.ksym
    
    
    
#####################################
#           list
#####################################  

def solvelist(*args):
    '''
    input: [vector with all eq=0], variables to find ..
    output: MyEq of each variable
    example:
        a+2*b=0 and 3*a-b=0
        ee=[a+2*b,3*a-b]
        then :
        a,b=solvelist(ee,a,b)
        return a,b in MyEq ecuation class
    '''
    vecs=args[1::]
    kres= solve(*args)
    var,value=kunpakDic(kres)
    vres=[]
    for i ,j in zip(var,value):
        ee=MyEq(j,str(i))
        vres.append(ee)
    return vres  

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
#           real subs
#####################################  
        
def subskwargs(QQ,**kwargs):
    mkey=[]
    vvalue=[]
    
    for key, value in kwargs.items():
        mkey.append(key)
        vvalue.append(value)
    kres=QQ 
    for i,j in zip(mkey,vvalue):
        valor=j
        if type(j)==MyEq:
            valor=j.ksym 
        kres=kres.subs(i,valor)
    return kres 

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

#####################################
#           algebra
#####################################
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
            
      
def passdoitI(kstr):
    klist=['>','<','=','True','∧',',']
    done=True
    for i in klist:
        if i in kstr:
            done=False
    return done        
def kreplace(ksym,var1,var2):
    if type(var2)==MyEq:
        var2=var2.ksym
    if type(ksym)==MyEq:
        ksym=ksym.ksym
       
    ksym=ksym.subs(var1.name,var2)
    return ksym
    
def kreplace(ksym,var1,var2):
    if type(var2)==MyEq:
        var2=var2.ksym
    if type(ksym)==MyEq:
        ksym=ksym.ksym
       
    ksym=ksym.subs(var1,var2)
    return ksym
def subskwargs2(ksym,**kwargs):
    var,val=unpack(kwargs)
     
    if type(ksym)==MyEq:
        kres=ksym.ksym
    for i,j in zip(var,val):
        
        kres=kreplace(ksym,i,j)
    return kres
 
def subskwargs(expr,**kwargs):
    ''' 
    QQ= symbols function
    ** kwargs c7=6,g4=z..etc..
    RETURN real substitucion when variable have underscore name like 'c_7' 'g_4'    
    '''
    key,value=unpack(kwargs)
    sres=str(expr)
    for i,j in zip(key,value):
        svalue=str(j)
        skey=i
        if len(skey)==2:
            nskey=skey[0]+'_'+skey[1]
            if nskey in sres:
                sres=sres.replace(nskey,svalue)
            else:
                sres=sres.replace(skey,svalue)
        else:
            sres=sres.replace(skey,svalue)
            
    return parse_expr(sres)

def subskwargs1(expr,**kwargs):
    ''' 
    QQ= symbols function
    ** kwargs c7=6,g4=z..etc..
    RETURN real substitucion when variable have underscore name like 'c_7' 'g_4'    
    '''
    if len(kwargs)>0:
        key,value=unpack(kwargs)
        sres=str(expr)
        for i,j in zip(key,value):
            svalue=str(j)
            skey=gval(i)
            if len(skey)==2:
                nskey=skey[0]+'_'+skey[1]
                if nskey in sres:
                    sres=sres.replace(nskey,svalue)
                else:
                    sres=sres.replace(skey,svalue)
            else:
                sres=sres.replace(skey,svalue)
                
        return parse_expr(sres)  
    else:
        return expr
# Tools expr value
def gval(kexpr):
     # return value of kexpr if expr is a MuEq type
     # ele return kexpr
    if type(kexpr)==MyEq:
        return kexpr.ksym
    else:
        return kexpr
        
def integerfactor(expr,var,name=''):
    # return integral( exp( expr),var)
    # if name!='' return MyEq 
     
    iexp=MyIntg(expr,'q',var=var,kshow=False)
     
    kres=exp(iexp.ksym)
    if name=='':
        return kres
    else:
        return MyEq(kres,kname=name)
            
def varDiff2(*args):
    mm=[]
    for i in args:
         
        k=symbols(i+"''")
        mm.append(k)
    return(mm) 
def varDiff(*args):
    mm=[]
    for i in args:
         
        k=symbols(i+"'")
        mm.append(k)
    return(mm) 
    
def textpos(ksolu,typos,tsize):
    qq=len(ksolu) 
    mm=[]
    for i in range(0,qq-1):
        
        p1=ksolu[i]
        p2=ksolu[i+1]
        L=p2-p1
        fac=L/10
        mm.append([p1+fac,typos,tsize,p1,p2])
    return mm

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
            ax.plot([i,i],[j,ymin],linestyle='-.',color='black')
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

def textpos(ksolu,typos,tsize):
    qq=len(ksolu) 
    mm=[]
    for i in range(0,qq-1):
        
        p1=ksolu[i]
        p2=ksolu[i+1]
        L=p2-p1
        fac=L/10
        mm.append([p1+fac,typos,tsize])
    return mm        
    
def laplaced(f,n):
    expr=laplace(f)
    expr=diff(expr,s)
    expr=(-1**n)*expr
    return expr
    
def laplaced(func):
    expr=laplace(func)
    expr2=integrate(expr,(s,s,oo))
    return expr2    

def Pplot(expr='',var=x,x1=-1,x2=1,x3=100,ymax='',ymin='',*args):
    '''
    expr = math function
    x1 = x min
    x2 = x max
    x3 = no point per line
    ymax= ymax
    ymin = ynim
    '''
    maa=getasint(expr,var)
    vecx=reparte(x1,x2,maa,x3)
    vecy=[]
    f = lambdify(var, expr, 'numpy')
    for i in vecx:
        yres = f(i)   
        vecy.append(yres)
    if ymax!='':
        if ymin=='':
            ymin=-1*ymax
        Vecx=[]
        Vecy=[]
        for i,j in zip(vecy,vecx):
            tempx=[]
            tempy=[]
            for yy,xx in zip(i,j):
                if yy<ymax and yy>ymin:
                    tempy.append(yy)
                    tempx.append(xx)
            Vecx.append(tempx)
            Vecy.append(tempy)
        vecx=Vecx
        vecy=Vecy
    fig, ax = plt.subplots()
    ax.grid(True, which='both')
    ax.axhline(0,color="black",linestyle='-.',alpha=0.5) #x-axis line
    ax.axvline(0,color="black",linestyle='-.',alpha=0.5)
    for X,Y in zip(vecx,vecy):
        ax.plot(X, Y)

    ax.grid(True, linestyle='-.')
    ax.tick_params(labelcolor='r', labelsize='medium', width=3)

    plt.show()

## fun used by Pplot()
def rrp(x1,x2,x3):
    fac=(x2-x1)/x3
    return np.arange(x1,x2+fac,fac)
def getasint(expr,var):
    ksolu=[]
    if Is_Div(expr):
        p1,p2=fraction(expr)
        ee=MyEq(p2,'ee',var=var,kshow=False)
        ksolu=ee.solve(var,'all')
        
    return ksolu
def reparte(x1,x2,vec,x3=100):
    mm=[]
    vec2=[]
    
    if len(vec)==0:
        kres=rrp(x1,x2,x3)
        mm.append(kres)
        return mm
    else:
        vec2=[]
        for i in vec:
            if i>x1  and i<x2 :
                vec2.append(i)
            vec=vec2    
        valores=[]
        p1=[x1]
        p2=[]
        for i in vec:
            p1.append(i+0.001)
            p2.append(i-0.001)
        p2.append(x2)
        for i,j in zip(p1,p2):
            mm.append(rrp(i,j,x3))
    return mm     
def deltasym(var):
    sdx='Δ'+str(var)
    dx=symbols(sdx)
    return dx


def realsub2(expr,**kwargs):
    vecs,vecv=unpack(kwargs)
    for i,j in zip(vecs,vecv):
        newj=j
        if type(j)==MyEq:
            newj=j.ksym
        
        if i in nombresdiferen:
            expr=expr.subs(eval(i),newj)
        
        else:
            expr=expr.subs(i,j)
    return expr 

def getdata(expr):
    kres=expr
    if type(expr)==MyEq:
        kres=expr.ksym
    return kres 

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
    mapval,mappos=filterdeep(mapval,mappos,deep)
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
        
def filterNegPos(vecd,vecp):
    NegMon=[]
    NegPos=[]
    OthMon=[]
    OthPos=[]
    for i,j in zip(vecd,vecp):
        if Is_NMono(i):
            NegMon.append(i)
            NegPos.append(j)
        else:
            OthMon.append(i)
            OthPos.append(j)
    for i,j in zip(OthMon,OthPos):
        if (-1*i) not in NegMon:
            NegMon.append(i)
            NegPos.append(j)
    return NegMon,NegPos
         
def filterdeep(vecd,vecp,deep):
 
    Othd=[]
    Othp=[]
    for i,j in zip(vecd,vecp):
        if len(j)<deep+1:
            Othd.append(i)
            Othp.append(j)

    return Othd,Othp 

def kreturn(obj):
    if type(obj)==UnevaluatedExpr:
        return unisymbols(obj)
    else:
        return obj
        
        
def area_triangulo(base, altura):
    return (base * altura) / 2

def area_rectangulo(base, altura):
    return base * altura

def area_rombo(diagonal_mayor, diagonal_menor):
    return (diagonal_mayor * diagonal_menor) / 2

def area_circulo(radio):
    return 3.1416 * radio ** 2  

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

def get_ksym(obj):
    if type(obj)==MyEqEq:
        return(obj.L-obj.R)
    elif type(obj)==MyEq:
        return obj.ksym
    elif type(obj)==str:
        return parse_expr(obj)
    else:
        return obj
        
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
def obj2MyEq(obj,var=x):
    obj=obj2func(obj)
    ee=MyEq(obj,'ee',var=var,kshow=False)
    return ee   
    
def get_diffvar(obj):
    var=obj.var
    if var==x:
        svar=symbols('dx')
    elif var==y:
        svar=symbols('dy')
    elif var==z:
        svar=symbols('dz')
    elif var==alpha:
        svar='d'+alphaname(alpha)
        svar=symbols(svar)
    else:
        svar='d'+str(var)
        svar=symbols(svar)
    return svar
def get_diffname(obj):
    var=obj.var
    if var==x:
        svar='dx'
    elif var==y:
        svar='dy'
    elif var==z:
        svar='dz' 
    
    if var==alpha:
        svar='da' 
    else:    
        svar='d'+str(var)
    return svar
def get_diffexpr(obj):
    dexpr=diff(obj.ksym,obj.var)
    return dexpr 

def get_fulldiff(obj):
    de=get_diffvar(obj)
    dd=get_diffexpr(obj)
    if dd==1:
        return de
    else:    
        return Mul(de,dd, evaluate=False) 
    
def MyFullDiff(obj):
    ksym=get_fulldiff(obj)
    var=obj.var
    kname=get_diffname(obj)
    ee=MyEq(ksym,'ee',var=var,kshow=False)
    ee.diffvar=get_diffvar(obj)
    return ee

def iinverse(expr):
    if Is_Div(expr):      
        p1,p2=fraction(expr)
        return cfrac(p2,p1)        
    elif expr==0:         
        return oo
    else:       
        return cfrac(1,expr)   
    return cfrac(1,expr)    
    
def eQrecta(*args):
    kname=''
    vecdata=[]
    for i in args:
        if type(i)==str:
            kname=i
        else:
            vecdata.append(i)
            
    if len(vecdata)==3:
        x1,y1,m=vecdata
        b=y1-x1*m
        kres=m*x+b
    else:
        x1,y1,x2,y2=vecdata 
        m=cfrac((y2-y1),(x2-x1))
        if type(m)==float:
            m=simplify(Rational(m))
            
        b=y1-x1*m
        kres=m*x+b
    if kname=='':
       return kres
    else:
       ee=MyEq(kres,kname=kname)
       return ee        