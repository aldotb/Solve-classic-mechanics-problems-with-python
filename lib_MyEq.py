from sympy import symbols
from functools import reduce
from IPython.display import  Math #,display 
import matplotlib.pyplot as plt
from libaldo_math2 import *
from libaldo_show import *
from lib_MyData import * 
import copy
#from lib_Func import * 
 
C1,C2,C3,C4,t,x=symbols('C1 C2 C3 C4 t x')
dataQ=[]
e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12=symbols('e1 e2 e3 e4 e5  e6 e7 e8 e9 e10 e11 e12')

        
            
class MyEq:
    def __init__(self, ksym, kname='',var2=t, var1=x, kp=False, kope='', kshow=True,ktype='P',xx='' ,dtype=1,depen=False,x1='',x2='',Pobj='',kfun=''):
        self.t=ktype
        self.type=ktype
        self.kinte=''
        self.depen=depen
        self.primi=''
        self.x1=x1
        self.x2=x2
        self.eeq=ksym
        self.var2=var2
        self.var1=var1
        self.name = alphaname(kname)
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
         
        if type(ksym)==MyEq:
            seq=str(ksym())
            self.ksym = unisymbols(opemat(ksym(), kope=kope))
            self.v = unisymbols(opemat(ksym(), kope=kope))
            self.eeq=ksym
        elif type(ksym)==MyEqMat:
            self.ksym = unisymbols(opemat(ksym(), kope=kope))
            
        else:    
            seq=str(ksym)
            self.ksym = unisymbols(opemat(ksym, kope=kope))
            self.v = unisymbols(opemat(ksym, kope=kope))

        if ktype=='F':
            if var2!='':
                s1=kname
                if type(var2)==list:
                    s2=alphaname(var2[0])+','+alphaname(var2[1])
                else:
                    s2=alphaname(var2)
                kname=s1+'_{('+s2+')}'
                self.name = kname

        if ktype=='I':
            self.primi=self.ksym
            if self.x1=='': 
                self.ksym=Integral(self.ksym, self.var2)
            else:
                self.ksym=Integral(self.ksym, (self.var2,self.x1,self.x2))
            self.histo = unisymbols(opemat(self.ksym, kope=kope))
            
        if ktype=='diff' or ktype=='diff2':
            f=Function(var1)(var2)
            kres=self.ksym
            kres=kres.subs(var1,f)
            var2=self.var2
            var1=self.var1
            self.Bdiff=kres
            if ktype=='diff':
                self.Adiff=Derivative(f, var2)
            else:
                self.Adiff=Derivative(f, var2,var2)
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
            else:
                kres=self.ksym
            
                if self.name == '':
                    display(Math(latex(kres)))
                else:
                    sR = self.name + ' ='
                    display(Math(sR + latex(kres)))
        
        

        # self.s()
        if dtype==1:
            if self not in dataQ and self.name!='': 
                dataQ.append(self)
            
    def __call__(self, kname='',kop='',kupdate=False,*args,**kwargs):
         
        if type(kname)!=str and len(kwargs)==0 and len(args)==0:
             
            var2=self.var2
            ksym=self.ksym
            return ksym.subs(var2,kname)
            
        if len(args)==1 and len(kwargs)==0 and self.iniG==0 and self.t=='Q':
            self.ksym=args[0]
            self.iniG=1
            self.s()
            return
        if len(args)==1 and len(kwargs)==0 and self.t!='Q':
            ksym=self.ksym
            var2=self.var2
            nvar=args[0]
            ee=MyEq(ksym,kname=self.name,kshow=False)
            ee.set(var2,nvar)
 
            return    
        if len(args)==1 and len(kwargs)==0 and self.iniG==1 and args[0]=='' and self.t=='Q':
            
            self.ksym=0
            self.iniG=0
            #self.s()
            return    
         
        if len(args)>0 and len(kwargs)==0:
            kvar=self.var2
            if kvar!='':
                kres=self.ksym
                kname=self.name
                ee=MyEq(kres,kname,kshow=False)
                if type(kvar)!=list:
                    ee.set(kvar,args[0])
                else:
                    for i in range(len(kvar)):
                       ee.set(kvar[i],args[i],kshow=False)
                    ee.s()    
                
                 
                
        if len(kwargs)>0:
            qk=len(kwargs)
            
            ee=MyEq(self.v,kshow=False)
            for key, value in kwargs.items():
                key=key.replace('_','')
                thesym=unisymbols(parse_expr(key))
                thevalue=unisymbols(value)
                ee.set(parse_expr(key),value,kshow=False)
            if kupdate:
                 
                self.update(ee.ksym)
                self.s()
                return
            if kname!='':
                ee.name=kname
                ee.s()
                return ee
            if len(args)==1:
                if type(args[0])==str:
                    ee.name=kname
                    ee.s()
                    return ee    
            return ee.ksym
        
        return self.ksym
     
    def __repr__(self):
        kres=str(self.ksym)

        return  kres
    def _latex(self, obj):
        return latex(self.ksym) 
         
    def __str__(self):
        return self.__repr__()    
        
    ###########################################
    #               Update
    ###########################################
    def __add__(self, other):
        """ Returns the vector addition of self and other """
        if type(other)==MyEq or type(other)==eQ:
            kres = self.ksym + other.ksym
        else:
            kres = self.ksym + other 
        return kres
    def __radd__(self, other):
        if type(other)==MyEq or type(other)==eQ:
            kres = self.ksym + other.ksym
        else:
            kres = self.ksym + other 
        return kres   
        
    def __sub__(self, other):
        if type(other)==MyEq or type(other)==eQ:
            kres = self.ksym - other.ksym
        else:
            kres = self.ksym - other 
        return kres
    def __rsub__(self, other):
        if type(other)==MyEq or type(other)==eQ:
            kres = self.ksym - other.ksym
        else:
            kres = self.ksym - other 
        return kres    
    def __mul__(self, other):
        """ Returns the vector addition of self and other """
        if type(other)==MyEq or type(other)==eQ:
            kres = self.ksym*other.ksym
        else:
            kres = self.ksym*other 
        return kres
    def __rmul__(self, other):
        """ Returns the vector addition of self and other """
        if type(other)==MyEq or type(other)==eQ:
            kres = self.ksym*other.ksym
        else:
            kres = self.ksym*other 
        return kres
       
 
    def __truediv__(self, other):
        if type(other)==MyEq or type(other)==eQ:
            kres = self.ksym/other.ksym 
        else:
            kres = self.ksym/other 
        return kres

    def __rtruediv__(self, other):
        if type(other)==MyEq or type(other)==eQ:
            kres = self.ksym/other.ksym 
        else:
            kres = self.ksym/other 
        return kres  
 
         

    
    
    def set_solve(self,kset,kvset,kvars,knames):
        ee=MyEq(self.ksym,kshow=False)
        ee.set(kset,kvset,kshow=False)
        kres=ee.solve(kvars)
        ee2=MyEq(kres,knames)
        return ee2
    def set_solveR(self,kset,kvset,kvars,knames):
        ee=MyEq(self.ksym,kshow=False)
        ee.set(kset,kvset,kshow=False)
        kres=ee.solveR(kvars)
        ee2=MyEq(kres,knames)
        return ee2    
        
    def xcopy(self,kname,kshow=True):
        ee=copy.deepcopy(self)
        if kname!='':
            ee.name=kname
        if kshow:
            ee.s()
        return ee
        
        
    def solveset(self,val1,val2,ksym,kshow=True):
         
        var2=unisymbols(self.var2)
        kres=unisymbols(self.ksym)
        kres1=kres.subs(var2,val1)
        if kshow:  
            kres2=csolve(kres1.subs(var2,val1)-val2,ksym,str(ksym)) 
        else:
            kres2=csolve(kres1.subs(var2,val1)-val2,ksym)
        self.set(ksym,kres2,kshow=kshow)
        
    def v(self):
        return self.v

    def update(self, kres):
         
        self.histo = self.ksym
        self.ksym = kres
        self.v = kres
        self.oldprimi=self.primi
        
    def upgrade(self,*args,kope='',kshow=True,**kwargs ):
        for i in args:
            kname=i.name
            self.set(kname,i,kshow=False,kope=kope)
        if len(kwargs)>0:
            self.set(**kwargs,kshow=False)
        if kshow:    
            self.s()
        
    def kret(self):
        return self.ksym

    def s2(self):
        kres=self.ksym
        if self.type=='diff' or self.type=='diff2':
            display(Math(latex(self.EqDiff)))
        else:
            if self.name == '':
                display(Math(latex(kres)))
            else:
                sR = self.name + ' ='
                display(Math(sR + latex(kres)))
            
    def s(self,**kwargs):
        
        qk=len(kwargs)
        if qk>0:
            ee=MyEq(self.v,kshow=False)
            for key, value in kwargs.items():
                ee.set(parse_expr(key),value,kshow=False,ktype=self.type)
            ee.s2()
   
        else:
            self.s2()       
            

    def sExp(self, kval):
        if self.name == '':
            display(Math(latex(self.ksym)))
        else:
            sR = self.name + ' ='
            display(Math(sR + latex(kval)))

    def undo(self):
        self.ksym = self.histo
        self.v = self.histo
        return self.ksym
        
    def norm_variable(self):
        clist=[C1,C2]
        m1=self.free()
        s1=[str(x) for x in m1]
        for kvar in clist:
            svar=str(kvar)
            for i in range(len(s1)):
                if s1[i]==svar:
                    kvar=m1[i]
        
    ###########################################
    #            Math operation
    ###########################################

    def reduF(self,kshow=True,kupdate=True):
          
        kres=self.ksym
        kres=factor(kres)
        kres2=1
        if Is_Mono(kres):
            kres=numer(kres)
             
            mm=fpoly(kres,'list')
            for i in mm:
                if Is_Poly(i):
                    kres2=kres2*i
                   
        if kres2!=0 and kres2!=1:
            self.ksym= kres2 
        self.s()

    def Add(self, kval, kope='', kupdate=True,kshow=True):
        kres = self.ksym
        if Is_Add(kres) and kope != '':
            mm = 0
            for i in fpoly(kres, 'list'):
                nval = i + kval
                nval = opemat(nval, kope=kope)
                mm += nval
            kres = mm

        else:
            kres = kres + kval
        kres = opemat(kres, kope=kope)
        if kupdate:
            self.update(kres)
            if kshow:
                return self.ksym
        else:
            return self.ksym
    

           
    def Mul(self, kval, kope='', kupdate=True,kshow=True):
         
        
            
        kres = self.ksym
        if Is_Add(kres):
            mm = 0
            for i in fpoly(kres, 'list'):
                nval = i * kval
                nval = opemat(nval, kope=kope)
                mm += nval
            kres = mm

        else:
            kres = kres * kval
        if '(' not in kope:
            kres = opemat(kres, kope=kope)
            if kupdate:
                self.update(kres)
                if kshow:
                    self.s()
            else:
                return self.ksym
        else:
            opemat(kres, kope=kope)
            
    def Div(self, kval, kope='', kupdate=True,kshow=True):
        kres = self.ksym
        if Is_Add(kres):
            mm = 0
            for i in fpoly(kres, 'list'):
                nval = i / kval
                nval = opemat(nval, kope=kope)
                mm += nval
            kres = mm

        else:
            kres = kres / kval

        kres = opemat(kres, kope=kope)
        if kupdate:
            self.update(kres)
            if kshow:
                self.s()
        else:
            return self.ksym

    def Pow(self, kval, kope='', ksec=False, kupdate=True,kshow=True):
        kres = kpow(self.ksym, kval)
        if ksec:
            kres = opematPolySec(kres, kope)
        else:
            kres = opemat(kres, kope)
        kres = opemat(kres, kope=kope)
        if kupdate:
            self.update(kres)
            if kshow:
                self.s()
        else:
            return self.ksym

    def Rpow(self, kval, kope='', ksec=False, kupdate=True,kshow=True):
        kres = rpow(self.ksym, kval)
        if ksec:
            kres = opematPolySec(kres, kope)
        else:
            kres = opemat(kres, kope)
        if kupdate:
            self.update(kres)
            if kshow:
                self.s()
        else:
            return self.ksym

    def cut_fac(self, kval):
        kres = cut_fac(self.ksym, kval)
        self.update(kres)
        return kres
    def signo(self):
        kres=self.ksym
        return signo(kres) 
            
    ###########################################
    #            get info
    ###########################################

    def list(self, kopt=''):
        kres = self.ksym
        kres = fpoly(kres, 'list')
        if kopt != '':
            return kres[kopt]
        else:
            return kres

    def free(self):
        kres = self.ksym
        kres = fpoly(kres, 'free')
        return kres

    def args(self, opt1=''):
        kres = self.ksym
        if opt1 == '':
            kres = kres.args
            return kres
        else:
            try:
                kres = fpoly(kres,'list')
                return  kres[opt1]
            except:
                return self.ksym

        ###########################################

    #            Set
    ###########################################
    def setdiff(self, kvar, kval, kshow=True, kope='', kret=False):
        self.set(knom=kvar.diff(), kval=kval, kshow=kshow, kope=kope, kret=kret)

    def evalif(self,*args):
        vsym=[]
        vval=[]
        done=True
        for i in args:
            if done:
                vsym.append(i)
                done=False
            else:
                vval.append(i)
                done=True
        kres=self.ksym
        for i,j in zip(vsym,vval):
            kres=kres.subs(i,j)
        ee=MyEq(kres,self.name)    
    
    def eval(self,**kwargs):
        kname=self.name
        ksym=self.ksym
         
        key, value = kunpakDic(kwargs)
        for i,j in zip(kkey,value):
            ksym=ksym.subs(i,j)
        MyEq(ksym,self.name)

    def float(self,  kupdate=True,kshow=True):
        kres=self.ksym
        try:
            kres=float(kres)
            
        except:
            self.s()
            return
        if kupdate:
            self.update(kres)
            if kshow:
                self.s()
        else:
            return self.ksym    
        
    def evalif(self, kname='',**kwargs):
        ee=self.xcopy(kname=kname,kshow=False)
        
        qk=len(kwargs)
        if qk>0:
            for key, value in kwargs.items():
                if type(value)==MyEq:
                    value=value.ksym
                ee.set(parse_expr(key),value,kshow=False)
                try:
                    ee.primi=ee.primi.subs(parse_expr(key),value)
                except:
                    done=False
        
            if kname!='':
                ee.s()
                return ee 
            else:
                return ee.ksym
            
        else:    
            self.s()
        
    def set(self, knom='', kval='', kshow=True, kope='', kret=False,andsolve=[],Bag='',**kwargs):
         
        if self.type=='Ph':
            self.Pobj.store_val(knom,kval)
            mm= self.Pobj.F
            qq=len(mm)
            for i in range (qq):
                pkres=mm[i][0]
                try:
                    pkres=pkres.subs(knom,kval)
                except:
                    done=False
                mm[i][0]=pkres
            self.Pobj.F=mm
         
        try:
            kres=self.primi
            self.primi=kres.subs(knom,kval)
        except:
            done=False
        qk=len(kwargs)
        if qk>0:
            for key, value in kwargs.items():
                self.set(parse_expr(key),value,kshow=False)
                try:
                    self.primi=self.primi.subs(parse_expr(key),value)
                except:
                    done=False
            if Bag!='':
                self.upBag(Bag,kshow=False)
            if kshow:
                self.s()    
            
        else:    
            if andsolve!=[]:

                if kval=='':
                    self.setValue(self.v,knom,kshow=False,kope=kope)
                    
                else:
                    self.setValue(knom=knom, kval=kval, kshow=False, kope=kope, kret=kret)
                    
                    ee=MyEq(self.solve(andsolve[0]),str(andsolve[1]))
                    return ee
                if Bag!='':
                    ee.upBag(Bag)
                return ee    
            
            else:
                if kval=='':
                    self.setValue(self.v,knom,kope=kope,kshow=False)
                    print(3)
                    if Bag!='':
                        self.upBag(Bag,kshow=False)
                else:
                    mm=self.setValue(knom=knom, kval=kval, kshow=False, kope=kope, kret=kret)
                     
                    if Bag!='':
                        self.upBag(Bag,kshow=False)
                if kshow:
                    if Bag=='':
                        self.s()
                    
    
                        
    def setValue(self, knom, kval, kshow=True, kope='', kret=False):
        if type(knom) != list:
            knom = [knom]
            kval = [kval]
        for i, j in zip(knom, kval):
            i = unisymbols(i)
            j = unisymbols(j)
            kres = self.ksym
            kres = kres.subs(i, j)
            kres = opemat(kres, kope)
            
            self.update(kres)
            # kres=self.simplify(kshow=False)
        if kret:
            return self.ksym

        if kshow:
            self.s()
        # return self.ksym

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
        self.update(kres)
        self.s()

    def MsetValue(self, vknom, vkval, kope=''):
        for knom, kval in zip(vknom, vkval):
            knom = unisymbols(knom)
            kval = unisymbols(kval)
            kres = self.ksym
            kres = kres.subs(knom, kval)
            kres = opemat(kres, kope=kope)

            kres = self.simplify()
            self.update(kres)
        self.s()
        # return self.ksym

    def multiSet(self, vknom, vkval, kope=''):
        for knom, kval in zip(vknom, vkval):
            knom = unisymbols(knom)
            kval = unisymbols(kval)
            kres = self.ksym
            kres = kres.subs(knom, kval)
            kres = opemat(kres, kope=kope)

            kres = self.simplify()
            self.update(kres)
        self.s()
        # return self.ksym

    def cut_denom(self, kope=''):
        kres = numer(self.ksym)
        kres = opemat(kres, kope)
        self.update(kres)
        return self.ksym

    def kdiff(self, kval, kope='', kupdate=False):
        kres = kdiff(self.ksym, kval)
        kres = opemat(kres, kope)
        if kupdate:
            self.update(kres)
            self.s()
        else:
            return kres

            
    def integral(self, kname='',x1='', x2='', kope='', kupdate=False,ktype='P',var2='',C1=C1):
        keq = self.ksym
        kres = keq
        if var2!='':
            kvar=var2
        elif self.var2!='':
            kvar=self.var2
        else:
            kvar=t
            
        
        if x1 == '':
            kres = integrate(keq, kvar)
            # self.update(kres)
        else:
            kres = integrate(keq, (kvar, x1, x2))
        kres = opemat(kres, kope)
        if kname!='':
            return MyEq(kres+C1,kname,var2=kvar,ktype=ktype)
        if kupdate:
            self.update(kres)
            self.s()
        else:
            return kres

    def tintegral_def(self, alpha1, a1, a2, kupdate=False):
        keq = self.ksym
        kres = tintegral_def(keq, alpha1, a1, a2)
        if kupdate:
            self.update(kres)
            self.s()
        else:
            return kres

    def only_nume(self, kope=''):
        kres = self.get_nume()
        kres = opemat(kres, kope)
        self.update(kres)
        self.s()

    def only_deno(self, kope=''):
        kres = self.get_deno()
        kres = opemat(kres, kope)
        self.update(kres)
        rself.s()

    def get_MonoExp(self):
        kres2 = self.ksym
        kres = kres2.fpoly('get', 1)
        return kres

    def evalue(self, *args, kope='', kshow=False,**kwargs):
        qk=len(kwargs)
        if qk>0:
            for key, value in kwargs.items():
                self.set(parse_expr(key),value,kshow=False)
                try:
                    self.primi=self.primi.subs(parse_expr(key),value)
                except:
                    done=False
        if len(args)==1:
            var2=self.var2
            
            kres=self(var2=args[0])
        return kres    
            
         

    def evalueArray(self, knom, vkval):

        kres = [self.evalue(knom, xx, kshow=False) for xx in vkval]
        return kres

    #   Functions Tools

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

    def angTan(self, **kwargs):
        kres = self.ksym
        var2 = self.var2
        kres = kres.diff(var2)
        if len(kwargs)>0:
            ee=MyEq(kres,kshow=False)
            kres=ee(**kwargs)
            if len(fpoly(kres,'list'))>1:
                mm=fpoly(kres,'list')
                qq=len(mm)
                kres=mm[qq-1]
         
        return atan(kres)


    def angOrto(self, **kwargs):
        kres =self.angTan(**kwargs)
        
        return  (kres+pi/2)

    def Is_Poly(self):
        return Is_Poly(self.ksym)

    def Is_Mono(self):
        return(Is_Mono(self.ksym)) 

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

        #   Transformacion
    def all_type(self):
        self.s()
        sE(['Monomie= ',self.Is_Mono(),'  ','Polynomie = ',self.Is_Poly()])
        sE(['Is Add= ',self.Is_Add(),'  ','Is Mul= ',self.Is_Mul()]) 
        
    
    def expand(self, kupdate=True, kshow=True, kope=''):
        kres = self.ksym
        kres = expand(kres)
        kres = opemat(kres, kope=kope)
        if kupdate:
            self.update(kres)
        if kshow:
            self.s()

    def simplify(self, kupdate=True, kshow=True, kope=''):
        kres = self.ksym
        kres = simplify(kres)
        kres = opemat(kres, kope=kope)
        if kupdate:
            self.update(kres)
        if kshow:
            self.s()

    def simplify_sec(self):
        kres = self.ksym
        if self.Is_Add():
            mm = 0
            for i in fpoly(kres, 'list'):
                mm += simplify(i)
            kres = mm
        kres = simplify(kres)
        self.update(kres)
        self.s()

    def factor(self, kvar='', kshow=True, kope=''):
        kres = self.ksym
        if kvar != '':
            g1 = self.fpoly('filt', kvar)
            g2 = kres - g1
            g3 = factor(g1)
            kres = g3 + g2
        else:
            kres = factor(kres)
        kres = opemat(kres, kope=kope)
        self.update(kres)
        if kshow:
            self.s()

    def texpand(self, kshow=True, kope=''):
        kres = self.ksym
        kres = expand_trig(kres)
        kres = opemat(kres, kope=kope)
        self.update(kres)
        if kshow:
            self.s()

    def tsimplify(self, kshow=True, kope=''):
        kres = self.ksym
        kres = trigsimp(kres)
        kres = opemat(kres, kope=kope)
        self.update(kres)
        if kshow:
            self.s()

    def opemat(self, kope='', kshow=True):
        kres = self.ksym
        kres = opemat(kres, kope=kope)
        self.update(kres)
        if kshow:
            self.s()

    def opematsec(self, kope='', kshow=True):  # equal to pemat but secuential
        kres = unisymbols(self.ksym)
        kres = opematsec(kres, kope=kope)
        self.update(kres)
        if kshow:
            self.s()

    def simplify_Rpow(self, kshow=True):
        kres = self.ksym
        kres = cut_root_of_pow(kres)

        self.update(kres)
        if kshow:
            self.s()

    def sort(self):
        kres = self.ksym
        klist = self.fpoly('list')
        mm = 0
        for i in klist:
            mm += i
        kres = mm
        self.update(kres)
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

    #  get Info
    def get_inside_root(self,kname=''):
        mm=str(self.ksym)
        cc=0
        kk=mm.find('sqrt')
        kk2=kk+4
        sres=''
        qq=len(mm)
        for i in range(kk2,qq):
            val=mm[i]
            if val=='(':
                cc+=1
            if val==')':
                cc-=1
            sres=sres+val
            if cc==0:
                kres= parse_expr(sres)
                if kname!='':
                    ee=MyEq(kres,kname)
                    return ee
                else:
                    return kres
        return self.ksym
        
    def get_nume(self):
        return numer(self.ksym)

    def get_deno(self):
        return denom(self.ksym)

    def get_diff(self, kval):
        kres = kdiff(self.ksym, kval)
        return kres

    def fpoly(self, kopt='', op2='', op3=''):
        kres = self.ksym
        return fpoly(kres, kopt=kopt, op2=op2, op3=op3)

    def get_expo(self): # return exponente from monomie expresion
        ksym=self.v
        if Is_Mono(ksym):
            mm=fpoly(ksym,'list')
            return mm[1]

    def get_killexpo(self):
        ksym=self.v
        if Is_Mono(ksym):
            mm=fpoly(ksym,'list')
            return mm[1]

    def part(self, vec):

        kres = self.ksym
        try:
            return part(kres, vec)
        except:
            return kres

    def solvediff(self, kvar, kd='', kremp=False, kope='', korden=''):
        return self.solve(kvar=kvar.diff(), kd=kd, kremp=kremp, kope=kope, korden=korden)

    def solve_if(self,ksym,kname='',andset='', **kwargs):
        kres=self.ksym
        sname=self.name
        ee=MyEq(kres,kname,kshow=False)
         
        Sval=0
        for key, value in kwargs.items():
            if  key==sname:
                Sval=value
            else:    
                ee.set(parse_expr(key),value,kshow=False)
        
        kres=ee.solve(ksym)+Sval
        if andset!='':
            self.set(andset,kres)
        else:    
            if kname=='':
                return kres
            else:
                ee1=MyEq(kres,kname)
                return ee1
        
        
    def solve(self, kvar, kname='', kremp=False, kope='', korden='',ktype='P',var2='',Bag=''):
        if str(kvar)=='None':
            kname=symbols(kname)
        keq = self.ksym
        # kvar=unisymbols(kvar)
        kres = csolve(keq, kvar, kope=kope, korden=korden)
        if kres==[]:
            kres = csolve(keq, kname, kope=kope, korden=korden)
        if self.type=='Ph':
            self.Pobj.store_val(kvar,kres)
            mm= self.Pobj.F
            qq=len(mm)
            for i in range (qq):
                pkres=mm[i][0]
                try:
                    pkres=pkres.subs(kvar,kres)
                except:
                    done=False
                mm[i][0]=pkres
            self.Pobj.F=mm
        if kname!='':
            if type(kres)==list:
                if len(kres)>1:
                    cc=1
                    kkres=[]
                    for i in kres:
                        kname1=kname+str(cc)
                        kres1=MyEq(i,kname1,kshow=False)
                        if Bag!='':
                            kres1.upBag(Bag,kope=kope)
                             
                        kres1.s()
                        cc+=1
                        kkres.append(kres1)
                    
                    return kkres
            else:
                 
                kres=MyEq(kres,kname,kshow=False)
                if Bag!='':
                    kres.upBag(Bag,kope=kope)
                     
                kres.s()
                 
                 
                return kres
                
                         
        else:
            return kres


    def solveR(self, kvar, kd='', kremp=False, kope='', korden='',Bag=''):
        if Bag!='':
            self.upBag(Bag=Bag,kshow=False)
        keq = self.ksym
        # kvar=unisymbols(kvar)
        kres = csolveR(keq, kvar, kope=kope)
        if kres==[]:
            kres = csolve(keq, kname, kope=kope, korden=korden)
        if kd != '':
            return MyEq(opemat(kres,kope=kope), kd)

        if kremp:
            nkres = kres.subs(kvar, kres)
        return opemat(kres,kope=kope)
    def toFoat(self):
        kres=self.ksym
        try:
            kres=float(self.ksym)
             
            self.update(kres)
            self.s()
        except:
            self.s()
    
    def opematsec(self, kope=''):  # equal to pemat but secuential
        kres = unisymbols(self.ksym)
        kres = opematsec(kres, kope=kope)
        self.update(kres)
        return kres

    def opemat_deno(self, kope=''):
        kres = unisymbols(self.ksym)
        kres = opemat_deno(kres, kope=kope)
        self.update(kres)
        return kres

    def opemat_nume(self, kope=''):
        kres = unisymbols(self.ksym)
        kres = opemat_nume(kres, kope=kope)
        self.update(kres)
        return kres

    def solve_tan(self, alpha):
        c = tan(alpha)

        kres = self.ksym
        kres.subs(sin(alpha), c / rpow(c * c + 1, 2))
        kres.subs(cos(alpha), 1 / rpow(c * c + 1, 2))
        kres1 = csolve(kres, tan(alpha))
        return kres1

    #   Polynomial
    def degree(self, kvar=0):
        kres = self.ksym
        return degree(kres, gen=kvar)

    def degree_list(self):
        kres = self.ksym
        return degree_list(kres)

    def reduce(self):
        kres = self.ksym
        kres = expand(kres)
        kres = simplify(kres)
        kres = factor(kres)
        kres = simplify(kres)
        try:
            kres = cut_root_of_pow(kres)
            self.update(kres)
        except:
            self.update(kres)
        return kres

    #   reductions Algorotmic

    def reduFac(self,kop='', kshow=True):  # retun Eq=0,= a*b+c*b.. = b(a+c)..=0  then return (a+c)
        
        
        kres=self.ksym
        kres=factor(kres)
        kres2=1
        if Is_Mono(kres):
            kres=numer(kres)
             
            mm=fpoly(kres,'list')
            for i in mm:
                if Is_Poly(i):
                    kres2=kres2*i
                   
        if kres2!=0 and kres2!=1:
            self.ksym= kres2 
        else:
            self.ksym= kres  
        self.s()   
        
        
        # kres = self.ksym
        # try:
            # mm = fpoly(factor(kres), 'list')
            # try:
                # mkres = 1
                # for i in mm:
                    # done=False
                    # if not Is_Poly(i):
                        # if Is_Root(i):
                            # mm=get_inside_root(i)
                            # if Is_Poly(mm):
                                # done=True
                    # else:
                        # done=True
                    # if done:                    
                        # mkres = mkres * i
                # if mkres != 0 and mkres != 1:
                    # self.update(mkres)
                # else:
                    # self.s(kopt='r')
                # if kshow:
                    # self.s()

            # except:
                # if kshow:
                    # self.s()
        # except:
            # if kshow:
                # self.s()

    def factorSec(self, kvar, kfiltro='.', kshow=True):
        kres = self.ksym
        kres = factorSec(kEq=kres, ksym=kvar, kfiltro=kfiltro)
        self.update(kres)
        if kshow:
            self.s()
    def factorSecV(self, *args):
        for i in args:
            self.factorSec(i,kshow=False)
        self.s()  
        
    def get_cofactor(self,kfactor):
        mm=self.list()
        for i in mm:
            j=fpoly(i,'list')
            if kfactor in j:
                kres=i/kfactor
                return kres
    def grupFac(self, ksym, kfiltro='.', kshow=True):
        kres = self.ksym
        kres = grupFac(kres, ksym, kfiltro='.')
        self.update(kres)
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

    def fixRootPow(self, kksym):
        self.setValue(rpow(kpow(kksym, 2), 2), kksym)

    #   Plot

    def kplot(self, ksym, x1, x2, x3):
        xx = np.linspace(x1, x2, x3)
        yy = [self(x) for x in xx]
        plt.plot(xx, yy)

        plt.ylabel(self.name)
        plt.xlabel(str(self.var2))
        ktitle = self.name + '(' + str(str(self.var2)) + ')'  
        plt.title(ktitle)
        plt.show()

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

    def upBag(self, Bag,kname='', kshow=True, kope=''):
        vs=Bag.vmain
        vv=Bag.vsolve
         
        for i,j in zip(vs,vv):
            self.set(i,j,kshow=False)
            self.opemat(kope,kshow=False)    
        if kshow:
            self.s()


    def upTriang(self, angul, T3,kope=''):
        if type(angul)==MyTriang:
            angu2=angul
            angul=T3
            T3=angu
        
            
        v1 = [sin(angul), cos(angul), tan(angul), kpow(sin(angul), 2), kpow(cos(angul), 2), kpow(tan(angul), 2)]
        v2 = [T3.sin(), T3.cos(), T3.tan(), kpow(T3.sin(), 2), kpow(T3.cos(), 2), kpow(T3.tan(), 2)]
        
        for i in range(len(v2)):
            v2[i]=opemat(v2[i],kope=kope)
         
         
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
            self.update(kres)
            self.s()
        except:
            self.s()
    
    def diffEq(self,kname='',var2='',ktype='P',typeD=1):
        self.oldprimi=self.primi
        kres=self.ksym
        if var2=='':
            var2=self.var2
        if typeD==2:
            return MyEq(Derivative(kres,var2),kname=kname,var2= var2,ktype='Diff') 
    def changediff(self,nkvar,nvalue='',x1='',x2=''):
        if self.type!='I':
            return
        if nvalue=='':
            xx=nkvar
        else:
            xx=nvalue
        if x1!='':
            self.x1=x1
            self.x2=x2
        else:
            if self.x1!='':
                x1=self.x1
                x2=self.x2
                var2=self.var2
    
                if str(var2) in str(x1):
                    self.x1=x1.subs(var2,xx)
                if str(var2) in str(x2):
                    self.x2=x2.subs(var2,xx)
        kprimi=self.primi
        
        kprimi=kprimi.subs(self.var2,xx)
        if str(self.var2) not in str(kprimi):
            kprimi=kprimi*nvalue
        else:
            kprimi=kprimi*xx.diff(nkvar)
        
        self.primi=kprimi
        self.var2=nkvar
        if self.x1=='': 
                self.ksym=Integral(kprimi, self.var2)
        else:
                self.ksym=Integral(kprimi, (self.var2,self.x1,self.x2))
             
        self.s()
                        
                    
                        
            
           
        
             
    def doit(self,kname='',kshow=True,c1=0,**kwargs):
        if len(kwargs)>0:
            ee=MyEq(self.ksym,kshow=False)
            kres=ee(**kwargs)
            ee=MyEq(kres,kshow=False)
            kres=ee.doit()+c1
             
        else:    
            kres=self.ksym
            kres=kres.doit()+c1
        if kname!='':
            ee=MyEq(kres,kname=kname)
            return ee
        else:    
            self.type='P'
            self.update(kres)
            if kshow:
                self.s()    

        
        
           
            

    def get_diff(self,kname='',var2=''):
        if kname=='' and var2=='':
            var2=self.var2
        if kname!='' and var2=='':
            if type(kname)==str:
                var2=self.var2      
            else:
                var2=kname
                kname=''
            
        kres =self.ksym
        kres2=diff(kres,var2)
        if kname!='':
            ee=MyEq(kres2,kname,var2=var2,ktype=self.type)
            return ee
        return kres2
        
    def kdiff(self, kname='',var2='',kope='', kupdate=False,ktype=''):
        if ktype=='':
            ktype=self.type
        kres = self.ksym
        if var2=='':
            var2=self.var2
        kres = self.ksym
        kres2=diff(kres,var2)    
        return kres2
        if  kname!='':
            ee=MyEq(kres2,kname=kname,var2=var2,ktype=ktype)
            return ee
        else:
            return kres2
    def update_inte(self):
        if self.type=='I':
            kres=self.primi
            if self.x1=='': 
                        self.ksym=Integral(kres, self.var2)
            else:
                self.ksym=Integral(kres, (self.var2,self.x1,self.x2))     
                
        
    def fac_integral(self):
        ksym=self.primi
        var2=self.var2
        
        if Is_Mono(ksym) and Is_Mul(ksym):
            kres=[x for x in ksym.args if str(var2) not in str(x)]
            kres2=[x for x in ksym.args if str(var2)   in str(x)]
            mono1=1
            mono2=1
            for i in kres:
                mono1=mono1*i
            for i in kres2:
                mono2=mono2*i    
            self.primi=mono2
            kfac=mono1
            self.primi=mono2
            if self.x1=='': 
                self.ksym=Integral(mono2, self.var2)
            else:
                self.ksym=Integral(mono2, (self.var2,self.x1,self.x2))     
            self.Mul(mono1)
        else:
            self.s()
    
    #  Algoritmos de diferencial
    # def diff2(self,  *args,kope='',var2='',ktype='P'):
        # kres = self.ksym
        # kope=''
        # kname='' 
        # if  len(args)>0:
             
            # for i in args:
                # if type(i)==str:
                    # kname=i
                # else:
                    # kres=kres.diff(i)
            # if kname!='':
                # ee=MyEq(kres,kname,var2=i,ktype='F')
                # return ee
            # else:
                # return kres
        # elif self.var2!='',and ktype=='F' and kname!='':
            # ee=MyEq
            
            
            
        
        # else:
            # kres=kres.diff(self.var2)
        # kres=opemat(kres,kope) 
        # return kres
    
    def dsolve(self,kname='',C1='',C2=''):
        kres=self.EqDiff
        dres=dsolve(kres)
        mdres=dres.args
        kres=mdres[1]
        if C1!='':
            kres=kres.subs('C1',C1)
        if C2!='':
            kres=kres.subs('C2',C2)    
        if kname!='':
            ee=MyEq(kres,kname=kname,var2=self.var2,ktype='F')
            return ee
        else:
            ee=MyEq(kres,kname=alphaname(self.var1),var2=self.var2,ktype='F')
            return ee
        
        
    
    #(self,ksym,kvar,kname='',kope='',kshow=True,Bag='',kcero=False,ktype='linear',ksolu='')    
    def Diff(self,kname='',var2='',ktype='P'):
        if var2=='':
            var2=self.var2
        kres=self.ksym
        kres=kres.diff(var2)
        if kname=='':
            return kres
        else:
            ee=MyEq(kres,kname=kname,var2=var2,ktype=ktype)
            return ee
    
    def applyDiff(self,var2=''):
        kres=self.ksym
        if var2=='':
            var2=self.var2
        kres=kres.diff(var2)
        self.update(kres)
        self.s()
    
    def diff(self,*args):
        kname=''
        mm=[]
        kres=self.ksym
        for i in args:
            if type(i)==str:
                kname=i 
            else:
                mm.append(i)
        for i in mm:
            kres=kres.diff(i)
        if kname=='':
            return kres
        else:
            ee=MyEq(kres,kname,var2=self.var2)
            return ee
                
                
        
        
        
           
         
             
            
            
    #
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
        self.update(kres)
        if kshow:
            self.s()
        else:
            pass
    def cos2sin(self, angu, korden=2, kope='', kshow=True):
        self.set(kpow(cos(angu), 3), (1 - kpow(sin(angu), 2)) * sin(angu), kshow=False)
        kres = self.ksym
        kres = cos2sin(kres, angu=angu, korden=korden, kope=kope)
        self.update(kres)
        if kshow:
            self.s()

    def get_diffEq(self,vx=''):

        
        kname=self.name+'d'
        Vx=self.xx
         
        if Vx!='':
             
            kres=self.ksym
            
            kres=diff(kres,Vx)
            kname=kname+'_'+str(Vx)
            self.EqDiff=MyDiff(kres,Vx,kname)
            return self.EqDiff 
        elif vx!='':
            self.xx=vx 
            kres=self.ksym
            kname=kname+'_'+str(vx)
            kres=diff(kres,vx)
            self.EqDiff=MyDiff(kres,vx,kname)
            return self.EqDiff       
        else:
            sE(['first set variable MyEq.xx= Valeue']) 
            
    def Integral_eQ(self,kname='',var2='',x1='',x2='',ktype='P'):
        ksym=self.ksym
        
        if var2=='':
            var2==self.var2
        if x1=='':
            x1=0
            x2=var2
        else:
            x1=x1
            x2=x2
        kres=Integral(ksym,(var2,x1,x2))
        if kname=='':
            return kres.doit()
        else:
            ee=MyEq(kres,kname=kname,var2=var2,ktype='I')
            return ee
    def area_value(self):
        ee=MyEq(self.ksym,var2=self.var2,x1=self.x1,x2=self.x2,kshow=False)
        ee.area(kshow=False)
        ee.doit(kshow=False)
        return ee.ksym
    
    def area(self,var2='',x1='',x2='',kshow=True):
        if x1!='':
            self.x1=x1
            self.x2=x2
        if var2=='':
            var2=self.var2
        else:
            self.var2=var2
        self.primi=self.ksym

        if self.x1=='': 
            kres=Integral(self.ksym, self.var2)
        else:
            kres=Integral(self.ksym, (self.var2,self.x1,self.x2))
        self.type='I'
        self.update(kres)
        if kshow:
            self.s()
        
            
        
    def Integral(self,kname='',var2='',x1='',x2='',ktype='P'):
        ksym=self.ksym
        
        if var2=='':
            var2==self.var2
        if x1=='':
            x1=self.x1
        if x2=='':
            x2=self.x2

        kres=Integral(ksym,(var2,x1,x2))
        if kname=='':
            self.ksym=kres
            self.s() 
        else:
            ee=MyEq(kres,kname=kname,var2=var2,ktype='I')
            return ee 

            
    def solveIntegral(self,kname='',var2='',x1='',x2='',kupdate=True):
        if x1=='':
            x1=self.x1
        if x2=='':
            x1=self.x2    
        if var2=='':
            var2=self.var2    
        ksym=self.ksym
        kres=ksym.doit()
        if kname!='':
            ee=MyEq(kres,kname=kname,var2=var2)
            return ee
        else:
            if kupdate:
                self.ksym=kres
                self.s() 
            else:
                return kres

        
    def get_InteEq(self,vx='',kname=''):

        if kname=='':
            kname='s'+self.name
        else:
            kname=kname    
        Vx=self.xx
         
        if Vx!='':
             
            kres=self.ksym
            
            
            
            self.EqInte=MyInteger(kres,Vx,kname)
            return self.EqInte
        elif vx!='':
            self.xx=vx 
            kres=self.ksym
            
            
            self.EqInte=MyInteger(kres,vx,kname)
            return self.EqInte       
        else:
            sE(['first set variable MyEq.xx= Value']) 
         

    def convMyFunc(self,knom='',vV=''):
        ksym=self.v 
        if vV=='':
            vV=self.free()
        elif type(vV)==tuple:
            vV=list(vV)
        ksym=self.v
        kk=MyFunc(knom,ksym,vV)
        return kk

    # def Integral(self,kvar,kname='',x1='',x2='',kope='',kshow=True,Bag='',kcero=False):
        # ksym=self.ksym
        # ee=MyInteger(ksym,kvar=kvar,kname=kname,x1=x1,x2=x2,kope=kope,kshow=False,Bag=Bag,kcero=kcero)
        # ee.s()
        # return ee
###########################################
#               END MyEq Class
###########################################        

def get_real_value(ee):
    if type(ee)==MyEq:
        kres=ee.ksym
    else:
        kres=ee
    return kres    
def show_main():
    for i in dataQ:
        i.s()
 
        
class MyEqMat:
    def __init__(self, ee1,ee2):
        if type(ee1)==MyEq:
            self.ksym1=ee1()
            k1=ee1()
            self.ee1=ee1
            self.Te1=1
        else:
            self.ksym1=ee1
            k1=ee1
            self.ee1=''
            self.Te1=0
             
        if type(ee2)==MyEq:
            self.ksym2=ee2()
            k2=ee2()
            self.ee2=ee2
            self.Te2=1
        else:
            self.ksym2=ee2 
            k2=ee2
            self.ee2=''
            self.Te2=0
            
        k1,k2=self.ksym1,self.ksym2
        
        self.ksym=MaT(k1,k2)
         
        
    def __call__(self,*args,kop='', **kwargs):
        if len(args)==0 and len(kwargs)==0:
            self.update()
            self.ksym=MaT(self.ksym1,self.ksym2)
            return self.ksym
            
 
    
    def update(self):
        if self.Te1==1:
            ee=self.ee1
            self.ksym1=ee()
        
        if self.Te2==1:
            ee=self.ee2
            self.ksym2=ee()    
         


 #
        


def upBagSys(ksys,Bag,kope=''):
    for i in ksys:
        i.upBag(Bag,kope=kope)
        
def My2Integer(*args,**kwargs):
        
        if len(args)==3:
            ksym,var1,kname=args
            var2=var1
        else:
            ksym,var1,var2,kname=args
        x1,x2,y1,y2,kope,kshow,Bag,kcero,ktype='','','','','',False,'',False,1
        Lksym=[x1,x2,y1,y2,kope,kshow,Bag,kcero,ktype]
        Lval=['','','','','',False,'',False,1]
             
        
        
        Lsval=['x1','x2','y1','y2','kope','kshow','Bag','kcero','ktype']
        kkey,kvalu=kunpakDic(kwargs)
        for i,j in zip(kkey,kvalu):
                Lksym[Lsval.index(i)]=j
            
        x1,x2,y1,y2,kope,kshow,Bag,kcero,ktype=Lksym
        ee1=MyInteger(ksym=ksym,kvar=var1,kname=kname,x1=x1,x2=x2,kope=kope,kshow=False,Bag=Bag,kcero=False,ktype=1) 
        kres=ee1.kinte
        ee2=MyInteger(ksym=kres,kvar=var2,kvar2=var1,kname=kname,x1=y1,x2=y2,y1=x1,y2=x2,kope=kope,kshow=True,Bag=Bag,kcero=False,ktype=2)
        
        return ee2
'''        
class MyInteger:
    def __init__(self,ksym,kvar,kname='',kvar2='' ,x1='',x2='',y1='',y2='',kope='',kshow=True,Bag='',kcero=False,ktype=1):
        self.ksym=unisymbols(ksym)
        self.kvar=unisymbols(kvar)
        self.kvar2=kvar2
        self.name=kname
        #if type(kvar2)==str:
            #self.name=kvar2
        self.x1=x1
        self.x2=x2
        self.y1=y1
        self.y2=y2
        self.type=ktype
        kres=ksym
        if self.type==2:
             try:
                kres=Integral(ksym,(kvar2,y1,y2))
             except:
                kres=Integral(ksym,kvar2)
        try:
            kres=Integral(kres,(kvar,x1,x2))
        except:
            kres=Integral(kres,kvar)        
           
        self.kinte=kres    
        if kshow:        
            sE([kname,'= ',self.kinte])
            
            
    def __call__(self):
        return self.kinte
         
        
    def set_limts(self,xx1=0,xx2=0):
        self.x1=xx1
        self.x2=xx2
        self.updateInte()
        self.s()
        
    def s(self):
        if self.ktype=='P':
            keq=self.ksym
        else:
            keq=self.kinte
        if self.name =='':
            display(Math(latex(keq)))
        else:
            sR=self.name+' ='
            display(Math(sR+latex(keq)))
    
    def updateInte(self):
        self.kback=[self.ksym,self.kvar,self.kinte]
        nksym=unisymbols(self.ksym)
        nkvar=unisymbols(self.kvar)
        if self.x2!='':
            self.kinte=Integral(nksym,(nkvar,self.x1,self.x2))
        else:
            self.kinte=Integral(nksym,nkvar )
    
    def undo(self):
        [self.ksym,self.kvar,self.kinte]=self.kback
        self.s()


    def opemat(self,kope):
        kres=unisymbols(self.ksym)
        e1=e1Q(kres,kope=kope,kshow=False)
         

        self.ksym=e1.v
        self.updateInte()
        self.s()
    def opematPolySec(self,kope=''):
        kres=self.ksym
        if  Is_Add(kres):
                mm=0
                for i in fpoly(kres,'list'):
                    e1=e1Q(i,kope=kope,kshow=False)
                    mm+=e1.v
                kres=mm    
        
        else:
            kres=opemat(kres,kope)
        self.ksym=kres 
        self.updateInte()
        self.s()    
            
    def setonlydiff(self,nvar):
        self.kvar=nvar
        self.updateInte() 
        self.s()

    def changediff(self,nvar,nksym='',kshow=True,kope=''):
        
        oksym=unisymbols(self.ksym)
         
        okvar=unisymbols(self.kvar)
             
        
        if nksym=='':
            kres=oksym.subs(okvar,nvar)
            kres=opemat(kres,kope=kope)
        else:
            kfac=integrate(nksym,nvar)
            kres=oksym
            try:
                kres=unisymbols(kres.subs(okvar,kfac))
                kres=unisymbols(opemat(kres,kope=kope))
            except:
                done=True
            kres=unisymbols(kres*nksym)
            kres=opemat(kres,kope=kope)
        
        
        self.ksym=kres
        self.kvar=nvar
        self.updateInte()       
        if kshow:
            self.s()
    def set_var_noncero(self,kv):
         
        vstr=str(kv) 
         
        v1=symbols(vstr,nonzero=True,extended_nonzero=True)
        kres=self.ksym 
        kres=kres.subs(unisymbols(kv),v1)
        self.ksym=kres
        self.updateInte()


    def set(self,knom,kval,kope='',kshow=True):
        if type(knom) != list:
            knom=[knom]
            kval=[kval]
        for i,j in zip(knom,kval):
            i=unisymbols(i)
            j=unisymbols(j)
            kres=unisymbols(self.ksym)
            try:
                kres=kres.subs(i,j)
            except:
                done=False    
             
            kres=opemat(kres,kope)
            #kres=strSubs(kres,knom,kval)
            self.ksym=unisymbols(kres)
         
        self.updateInte()       
        if kshow:
            self.s()
    
    def convertMyEq(self):
        ee=self.name
        eek=self.ksym
        ee=MyEq(self.ksym,self.name)
        
        
    def doit(self,x1='',x2='' ):
        print(1)
        if x1=='':
             
            kres=self.kinte
            kres=kres.doit()
            kname=self.name
            newF=parse_expr(kname)
            newF= MyEq(kres,kname=kname)
            newF.s()
            print(1)
            return newF
        

    def integer(self,x1='',x2=''):
        if x1=='':
            kres=self.kinte
            return kres.doit()
        else:
            nksym=self.ksym
            nkvar=self.kvar
            kres=integrate(nksym,(nkvar,x1,x2))
            return kres
    
    def upTriang(self,angul,T3):
        v1=[sin(angul),cos(angul),tan(angul),kpow(sin(angul),2),kpow(cos(angul),2),kpow(tan(angul),2)]
        v2=[T3.sin(),T3.cos(),T3.tan(),kpow(T3.sin(),2),kpow(T3.cos(),2),kpow(T3.tan(),2)]
        self.set(v1,v2)

    # Math simplify
    
    def expand(self,kupdate=True,kshow=True):
        kres=self.ksym
        kres=opemat(kres,'e')
        if kupdate:
            self.ksym=kres
            self.updateInte()    
        else:
            return kres
            
        if kshow:
            self.s()
    
    def simplify(self,kupdate=True,kshow=True):
        kres=self.ksym
        kres=opemat(kres,'s')
        if kupdate:
            self.ksym=kres
            self.updateInte()    
        else:
            return kres
            
        if kshow:
            self.s()
            
    
         
        
        
    
    def texpand(self,kupdate=True,kshow=True):
        kres=self.ksym
        kres=opemat(kres,'x')
        if kupdate:
            self.ksym=kres
            self.updateInte()    
        else:
            return kres
            
        if kshow:
            self.s() 
    
    def tsimplify(self,kupdate=True,kshow=True):
        kres=self.ksym
        kres=opemat(kres,'t')
        if kupdate:
            self.ksym=kres
            self.updateInte()    
        else:
            return kres
            
        if kshow:
            self.s()

    def Mul(self,kval,kupdate=True,kshow=True):
        kres=self.ksym
        kres=kres*kval
        if kupdate:
            self.ksym=kres
            self.updateInte()    
        else:
            return kres
            
        if kshow:
            self.s()

    def Add(self,kval,kupdate=True,kshow=True):
        kres=self.ksym
        kres=kres+kval
        if kupdate:
            self.ksym=kres
            self.updateInte()    
        else:
            return kres
            
        if kshow:
            self.s()                

    def kk(self):
        return MyEq(5,'kk')
'''
###########################################
#               END MyIteger Class
########################################### 


class MyDiff:
    def __init__(self,ksym,kvar,kname='',kope='',kshow=True,Bag='',kcero=False,ktype='linear',ksolu=''):
        f=symbols(kname, cls=Function)
        self.ksym=unisymbols(ksym)
        self.kvar=unisymbols(kvar)
        self.name=kname
        
        self.dvar='dot'
        self.ksymbols=symbols(str(kvar)+'dot')
        self.grade=1
        self.v=unisymbols(ksym)
        self.ktype=ktype
        x=unisymbols(kvar)
        if ktype=='square':
            self.eq1=f(x).diff(x,x)
        else:
            self.eq1=f(x).diff(x)	
        self.eq2=unisymbols(ksym)
        
        self.ode=Eq(self.eq1,self.eq2)
        mm=self.ode 
        if kshow:
            MyEq(mm)
    
    def __call__(self):
        kres=self.ode
        return kres.rhs
        
        
    def s(self):
        kk=self.ksymbols
        ksym=self.ksym
        self.upv()
        return self.ode


    def upv(self):
        self.v=self.ksym


    def diff(self,nkvar=''):
        kres=self.ksym
        kvar=self.kvar
        kgrade=self.grade
        if nkvar=='':
            self.ksym=kres.diff(kvar)
            self.grade+=1
            self.dvar='d'+self.dvar 
            self.ksymbols=symbols(str(kvar)+self.dvar)    
        self.s()

    def changediff(self,nkvar='',nksym=''):
        if nkvar!='':
            kres=self.ksym
            kvar=self.kvar
            kgrade=self.grade
            if nksym=='' and kgrade==1  :
                kres=kres.subs(kvar,nkvar)
            else:
                kres=kres.subs(kvar,nksym)
                kresd= diff(nksym,nkvar)
                kres=kres*kresd   
            self.ksym=kres
            self.kvar=nkvar 
            self.dvar='dot'  
            self.ksymbols=symbols(str(nkvar)+self.dvar)
        self.s()

    def updateDiff(self):
        f=  symbols(self.name, cls=Function)
        self.kback=[self.ksym,self.kvar]
        nksym=unisymbols(self.ksym)
        nkvar=unisymbols(self.kvar)
        #self.kinte=Derivative(nksym,nkvar)
        self.upv()
        x=nkvar
        self.eq2=self.ksym

        if self.ktype=='square':
            self.eq1=f(x).diff(x,x)
        else:
            self.eq1=f(x).diff(x)
        self.ode=Eq(self.eq1,self.eq2)

    def ss(self):
        mm=self.ode 
        MyEq(mm)
    
    def undo(self):
        [self.ksym,self.kvar]=self.kback
        self.s()


    def opemat(self,kope):
        kres=unisymbols(self.ksym)
        e1=e1Q(kres,kope=kope,kshow=False)
         

        self.ksym=e1.v
        self.updateDiff()
        self.s()
    def opematPolySec(self,kope=''):
        kres=self.ksym
        if  Is_Add(kres):
                mm=0
                for i in fpoly(kres,'list'):
                    e1=e1Q(i,kope=kope,kshow=False)
                    mm+=e1.v
                kres=mm    
        
        else:
            kres=opemat(kres,kope)
        self.ksym=kres 
        self.updateDiff()
        self.s()    
            
    def setonlydiff(self,nvar):
        self.kvar=nvar
        self.updateDiff() 
        self.s()

    def changediff2(self,nvar,nksym='',kshow=True,kope=''):
        
        oksym=unisymbols(self.ksym)
         
        okvar=unisymbols(self.kvar)
             
        
        if nksym=='':
            kres=oksym.subs(okvar,nvar)
            kres=opemat(kres,kope=kope)
        else:
            kfac=integrate(nksym,nvar)
            kres=oksym
            try:
                kres=unisymbols(kres.subs(okvar,kfac))
                kres=unisymbols(opemat(kres,kope=kope))
            except:
                done=True
            kres=unisymbols(kres*nksym)
            kres=opemat(kres,kope=kope)
        
        
        self.ksym=kres
        self.kvar=nvar
        self.updateDiff()       
        if kshow:
            self.s()

    def set(self,knom,kval,kope='',kshow=True):
        if type(knom) != list:
            knom=[knom]
            kval=[kval]
        for i,j in zip(knom,kval):
            i=unisymbols(i)
            j=unisymbols(j)
            kres=unisymbols(self.ksym)
            try:
                kres=kres.subs(i,j)
            except:
                done=False    
             
        kres=opemat(kres,kope)
            #kres=strSubs(kres,knom,kval)
        self.ksym=unisymbols(kres)
         
        self.updateDiff()       
        if kshow:
            self.s()
            return self.ode

    def doit(self,kname='',kop=''):
        kres=self.ksym 
        kvar=self.kvar
        if kop=='':
            kres=kres.diff(kvar)
        if kname!='':
            return MyEq(kres,kname)
        else:
            return kres
         

        
    def integer(self,x1='',x2=''):
        if x1=='':
            kres=self.kinte
            return kres.doit()
        else:
            nksym=self.ksym
            nkvar=self.kvar
            kres=integrate(nksym,(nkvar,x1,x2))
            return kres

    def upTriang(self,angul,T3):
        kres=self.ksym
        v1=[sin(angul),cos(angul),tan(angul),kpow(sin(angul),2),kpow(cos(angul),2),kpow(tan(angul),2)]
        v2=[T3.sin(),T3.cos(),T3.tan(),kpow(T3.sin(),2),kpow(T3.cos(),2),kpow(T3.tan(),2)]
        ee=e0Q(kres,kshow=False)
        ee.set(v1,v2,kshow=False)
        kres=ee.v
        self.ksym=kres
        self.updateDiff()       
        return self.ode
        
         

    # Math simplify

    def expand(self):
        kres=self.ksym
        kres=opemat(kres,'e')
        self.ksym=kres
        self.s()

    def simplify(self):
        kres=self.ksym
        kres=opemat(kres,'s')
        self.ksym=kres
        self.s()

    def tsimplify(self):
        kres=self.ksym
        kres=opemat(kres,'t')
        self.ksym=kres
        self.updateDiff()
        self.ss()

    def texpand(self):
        kres=self.ksym
        kres=opemat(kres,'x')
        self.ksym=kres
        self.s()

    def opemat(self,kope=''):
            kres=self.ksym
            kres=opemat(kres,kope=kope)
            self.ksym=kres
            self.s()

    def Mul(self,kval,kupdate=True,kshow=True):
        kres=self.ksym
        kres=kres*kval
        if kupdate:
            self.ksym=kres
            self.updateDiff()    
        else:
            return kres
            
        if kshow:
            self.s()

    def Add(self,kval,kupdate=True,kshow=True):
        kres=self.ksym
        kres=kres+kval
        if kupdate:
            self.ksym=kres
            self.updateDiff()    
        else:
            return kres
            
        if kshow:
            self.s()


    def Pow(self,kupdate=True,kshow=True):
        self.s()
         
        kres=self.ksym
        kres=diff(kres,self.kvar)
        kres=kpow(kres,2)
        if kupdate:
            self.ksym=kres
            self.ktype='square'
            self.updateDiff()
            return self.s()
        else:
            return kres
        if kshow:
            self.s()
    def set_var_noncero(self,kv):
         
        vstr=str(kv) 
         
        v1=symbols(vstr,nonzero=True,extended_nonzero=True)
        kres=self.ksym 
        kres=kres.subs(unisymbols(kv),v1)
        self.ksym=kres
         
        self.updateDiff()
        
    def dsolve(self,Vv='',var2=''):
        
            
        f=  symbols(self.name, cls=Function)
        x=unisymbols(self.kvar) 
        if self.ktype=='square':
            self.eq1=unisymbols(f(x).diff(x,x))
        else:
            self.eq1=unisymbols(f(x).diff(x))
        self.eq2=unisymbols(self.ksym)
        
        self.ode=Eq(self.eq1,self.eq2)
        ode=unisymbols(self.ode)
        sol = unisymbols(dsolve(ode,f(x)))
        kres= unisymbols(sol.rhs)
        skres=str(kres)
        if 'Piecewise(' in skres:
            kres=fix_otherwise(kres)
            kres=unisymbols(kres)
        if Vv=='':
            return kres
        else:
            e1=MyEq(kres,Vv,var2=var2,kshow=False)
            e1.norm_variable()
            e1.s()
            return e1
            
        
             
            
           
        
    def diffdiff(self,kname=''):
        vv=self.kvar 
        kres=self.ksym
        kres=diff(kres,vv)
        kres=MyDiff(kres,vv,kname)
        kres.s()
        return kres
  

        
class MyEqEq:
    def __init__(self, e1,e2,kname='',kshow=True):
        self.e1=e1
        self.e2=e2

        self.ksym=e1.v-e2.v
        self.eQ=Eq(e1.v,e2.v) 
        self.v=e1.v-e2.v
        self.name=kname  
        if kshow:
            if self.name == '':
                display(Math(latex(self.eQ)))
            else:
                sR = self.name + ')..  '
                display(Math(sR +latex(self.eQ)))    

    def update(self):
        e11=self.e1
        e22=self.e2
        self.ksym=e11.v-e22.v
        self.eQ=Eq(e11.v,e22.v)
        self.v=self.e1.v-self.e2.v    

    def s(self):
        self.update()
        if self.name == '':
            display(Math(latex(self.eQ)))
        else:
            sR = self.name + ')..'
            display(Math(sR   +latex(self.eQ)))
    
    def set(self,ksym,kval,kop='RL'):
        if 'L' in kop:
            self.e1.set(ksym,kval,kshow=False)
        if 'R' in kop:
            self.e2.set(ksym,kval,kshow=False)
        self.s()
        
    def expand(self,kop='RL'):
        if 'L' in kop:
            self.e1.expand(kshow=False)
        if 'R' in kop:
            self.e2.expand(kshow=False)
        self.s()

    def simplify(self,kop='RL'):
        if 'L' in kop:
            self.e1.simplify(kshow=False)
        if 'R' in kop:
            self.e2.simplify(kshow=False)
        self.s()

    def factor(self,kop='RL'):
        if 'L' in kop:
            self.e1.factor(kshow=False)
        if 'R' in kop:
            self.e2.factor(kshow=False)
        self.s()    

    def tsimplify(self,kop='RL'):
        if 'L' in kop:
            self.e1.tsimplify(kshow=False)
        if 'R' in kop:
            self.e2.tsimplify(kshow=False)
        self.s()

    def tfactor(self,kop='RL'):
        if 'L' in kop:
            self.e1.tfactor(kshow=False)
        if 'R' in kop:
            self.e2.tfactor(kshow=False)
        self.s()        

    def Add(self,kval,kop='RL'):
        if 'L' in kop:
            self.e1.Add(kval,kshow=False)
        if 'R' in kop:
            self.e2.Add(kval,kshow=False)
        self.s() 

    def Mul(self,kval,kop='RL'):
        if 'L' in kop:
            self.e1.Mul(kval,kshow=False)
        if 'R' in kop:
            self.e2.Mul(kval,kshow=False)
        self.s() 
        
    def Div(self,kval,kop='RL'):
        if 'L' in kop:
            self.e1.Div(kval,kshow=False)
        if 'R' in kop:
            self.e2.Div(kval,kshow=False)
        self.s() 
        
    def Pow(self,kval,kop='RL'):
        if 'L' in kop:
            self.e1.Pow(kval,kshow=False)
        if 'R' in kop:
            self.e2.Pow(kval,kshow=False)
        self.s()  
          
    def Rpow(self,kval,kop='RL'):
        if 'L' in kop:
            self.e1.Rpow(kval,kshow=False)
        if 'R' in kop:
            self.e2.Rpow(kval,kshow=False)
        self.s()
    
    def reduFac(self,kop='RL'):
        if 'L' in kop:
            self.e1.reduFac()
        if 'R' in kop:
            self.e2.reduFac()
        self.s()
    
    
    def solve(self,kval,kname='',kope='', korden=''): 
        kres=self.e1.v-self.e2.v
        e0=e0Q(kres,kshow=False)
        if kname=='':
            return e0.solve(kval,kope=kope,korden=korden)
        else:
            return e0.solve(kval,kd=kname, kope=kope, korden=korden)

#  New Eqes

class MyBag:
    def __init__(self,vmain=[],vsolve=[],vsecu=[],**kwargs):
        if len(kwargs)>0:
            vv,vs=kunpakDic(kwargs)
            mm=[]
            for i in  vv:
                mm.append(parse_expr(i)) 
            self.vmain=mm
            self.vsolve=vs
        else: 
            self.vmain=vmain
            self.vsolve=vsolve
        self.dataV=[]
        self.dataS=[]
        self.qq=0

        self.vsecu=vsecu
     
    def __call__(self,ksym='',kope=''):
        kres=ksym
        ee=MyEq(kres,kshow=False)
         
        ee.upBag(self,kshow=False)
        kres=ee.ksym
        if Is_Number(kres):
            ee.float(kshow=False)
            kres=ee.ksym 
        return kres   

def creaQ(ksym,kope):
    kres=ksym
    qq=len(kope)
    mm=kope[1:-1]
    ktype='P'
     
     
    vec=list(mm.split(","))
    qq=len(vec)
    if qq==1:
        return MyEq(ksym,vec[0])
    if qq==2:
        return MyEq(ksym,vec[0],ktype=vec[1])
    if qq==3:
        if vec[1]=='F':
         return MyEq(ksym,vec[0],ktype='F',var2=vec[2])
         
            
    
    

    


def FxQ(P, kname='Fx', kope=''):
    return MyEq(P.x_res(), kname, kope=kope,ktype='Ph',Pobj=P)
    
def FrxQ(P,ang=0, kname='Frx', kope=''):
    xx=P.x_res()
    yy=P.y_res()
    kres=yy*sin(ang)+xx*cos(ang)
     
    return MyEq(kres, kname, kope=kope,ktype='Ph',Pobj=P)    
    
def FyQ(P, kname='Fy', kope=''):
    return MyEq(P.y_res(), kname, kope=kope,ktype='Ph',Pobj=P)
    
def FryQ(P,ang=0, kname='Fry', kope=''):
    xx=P.x_res()
    yy=P.y_res()
    kres=yy*cos(ang)-xx*sin(ang)
     
    return MyEq(kres, kname, kope=kope,ktype='Ph',Pobj=P)    

def ToQ(P, kname='To', kope='', x1=0, y1=0,kshow=True,forza_pos=True):
    kres=P.torque(x1, y1)
    kres=kres*(signo(kres))
    return MyEq(kres, kname, kope=kope,ktype='Ph',Pobj=P,kshow=kshow)
def XmaxQ(P, kname='Xmax', kope='',kshow=True):
    kres=P.x_max()
    return MyEq(kres, kname, kope=kope,ktype='Ph',Pobj=P,kshow=kshow)
def YmaxQ(P, kname='Ymax', kope='',kshow=True):
    kres=P.y_max()
    return MyEq(kres, kname, kope=kope,ktype='Ph',Pobj=P,kshow=kshow)
  
def TflyQ(P, kname='Tfly', kope='',kshow=True):
    kres=P.t_fly()
    return MyEq(kres, kname, kope=kope,ktype='Ph',Pobj=P,kshow=kshow)

    
def EtotalQ(P, kname='E_t', kope='',kshow=True):
    return MyEq(P.Etotal(), kname, kope=kope,ktype='Ph',Pobj=P,kshow=kshow)
    
def MaQ (P,kname='Ma',kope='',kshow=True):
     
    kres=P.get_InerTotal()
    if kres==0:
        kres=P.In
    try:    
        kres=kres*aw
    except:
        kres=kres*P.aw
    return MyEq(kres, kname, kope=kope,ktype='Ph',Pobj=P,kshow=kshow)

def eQrot (P,kname='eqR',kope='s',andsolve='',set_dire='positive'):
    ''' 
        Data from : 
                P.add_forza( ) ... Torque
                P.add_Inercia() ...  Momen Angular
                
        Return :
                MyEq = To - Mo=I()*aw
        Options:
                eQrot (P, ,andsolve='aw') ..solve acc angular
                eQrot (P, ,andsolve='w') ..solve vel angular
                eQrot (P, ,andsolve='vt') ..solve tang velocity
                eQrot (P, ,andsolve='at') ..solve acc tangen
    '''    
    aw=symbols('aw')
    kshow=True
    if andsolve!='':
        kshow=False
    To=ToQ(P,kshow=kshow)
    
    Ma=MaQ(P,kshow=kshow)
    kres=unisymbols(To.ksym-Ma.ksym)
    if set_dire!='positive':
        kres=unisymbols(To.ksym+Ma.ksym)
        
    ee=MyEq(kres,kname=kname,kope=kope,ktype='Ph',Pobj=P,kshow=kshow)
    
    

    kaw= unisymbols(ee.solve(aw))
     
    if andsolve=='aw':
        kres3= unisymbols(ee.solve(aw))
         
        
        return MyEq(kres3, 'a_w', kope=kope,ktype='Ph',Pobj=P)
        
    
    
    elif andsolve=='w':
        kres3= unisymbols(ee.solve(aw))
        kres2=kres3*(P.t)
        P.w=kres2
        return MyEq(kres2, 'w', kope=kope,ktype='Ph',Pobj=P)
        
    elif andsolve=='at':
        kres3= ee.solve(aw)
        kres2=kres3*(P.r)
        return MyEq(kres2, 'a_t', kope=kope,ktype='Ph',Pobj=P) 
    elif andsolve=='vt':
        kres3= ee.solve(aw)
        kres2=kres3*(P.r)/t
        return MyEq(kres2, 'v_t', kope=kope,ktype='Ph',Pobj=P)
    elif andsolve=='an':
        kres3= unisymbols(ee.solve(aw))
        kres=kres3*kres3/(P.r*P.r)
         
        return MyEq(kres, 'a_n', kope=kope,ktype='Ph',Pobj=P)
         
    else :    
        return MyEq(To()-Ma(),kname=kname,kope=kope,ktype='Ph',Pobj=P,kshow=False) 
        
def FxyQ(P, kname='Fxy', kope=''):
    return MyEq(P.xy_res(), 'Fxy', kope=kope,ktype='Ph',Pobj=P)

def PoQ(P, kname='', kope=''):
    if kname == '':
        return MyEq(P.energia('P'), 'Po', kope=kope,ktype='Ph')
    else:
        return MyEq(P.energia('P'), kname, kope=kope,ktype='Ph')
def EkQ(P, kname='', kope=''):
    if kname == '':
        return MyEq(P.energia('K'), 'Ek', kope=kope,ktype='Ph')
    else:
        return MyEq(P.energia('K'), kname, kope=kope,ktype='Ph')
def axQ(P, kname='', kope=''):
    if kname == '':
        return MyEq(P.x_res() / P.m, 'a_x', kope=kope,ktype='Ph')
    else:
        return MyEq(P.x_res() / P.m, kname, kope=kope,ktype='Ph')
def ayQ(P, kname='', kope=''):
    if kname == '':
        return MyEq(P.y_res() / P.m, 'a_y', kope=kope,ktype='Ph')
    else:
        return MyEq(P.y_res() / P.m, kname, kope=kope,ktype='Ph')
def axyQ(P, kname='', kope=''):
    if kname == '':
        return MyEq(P.xy_res() / P.m, 'a_x', kope=kope,ktype='Ph')
    else:
        return MyEq(P.xy_res() / P.m, kname, kope=kope,ktype='Ph')
def EtQ(P, kname='', kope=''):
    if kname == '':
        return MyEq(P.energia(), 'E_t', kope=kope,ktype='Ph')
    else:
        return MyEq(P.energia(), kname, kope=kope,ktype='Ph')
        
        
def InerQ(P, kname='', kope=''):
    if kname == '':
        return MyEq(P.I_n, 'I_n', kope=kope,ktype='Ph')
    else:
        return MyEq(P.I_n , kname, kope=kope,ktype='Ph')
def YcQ(P, kname='', kope=''):
    if kname == '':
        return MyEq(P.cgravedad('y'),'Y_c', kope=kope,ktype='Ph') 
    else:
        return MyEq(P.cgravedad('y'),kname, kope=kope,ktype='Ph')
def XcQ(P, kname='', kope=''):
    if kname == '':
        return MyEq(P.cgravedad('x'),'X_c', kope=kope,ktype='Ph')  
    else:
        return MyEq(P.cgravedad('x'),kname, kope=kope)               
def e0Q(ksym, kope='', kshow=True):
    return MyEq(ksym, 'e0', kope=kope, kshow=kshow,ktype='C')
def e1Q(ksym, kope='', kshow=True):
    return MyEq(ksym, 'e1', kope=kope, kshow=kshow,ktype='C')
def e2Q(ksym, kope='', kshow=True):
    return MyEq(ksym, 'e2', kope=kope, kshow=kshow,ktype='C')
def e3Q(ksym, kope='', kshow=True):
    return MyEq(ksym, 'e3', kope=kope, kshow=kshow,ktype='C')
def e4Q(ksym, kope='', kshow=True):
    return MyEq(ksym, 'e4', kope=kope, kshow=kshow,ktype='C')
def e5Q(ksym, kope='', kshow=True):
    return MyEq(ksym, 'e5', kope=kope, kshow=kshow,ktype='C')
def e6Q(ksym, kope='', kshow=True):
    return MyEq(ksym, 'e6', kope=kope, kshow=kshow,ktype='C')
def e7Q(ksym, kope='', kshow=True):
    return MyEq(ksym, 'e7', kope=kope, kshow=kshow,ktype='C')
def e8Q(ksym, kope='', kshow=True):
    return MyEq(ksym, 'e8', kope=kope, kshow=kshow,ktype='C')

#################################################
#   Solve Algorithm
    

def solverSys (*args,Bag=''):
    Ke=[]   
    Kv=[]     
    Kn=[]
    
    for i  in args:
        if type(i)== MyEq:
            if Bag!='':
                i.upBag(Bag,kshow=False)
            Ke.append(i)
        if type(i)== Symbol:     
            Kv.append(i)
            Kn.append(i.name)
    #return(Ke,Kv,Kn)
    
    return MyEqSolveLin(Ke,Kv,Kn,Bag=Bag)
    
def MyEqSolveLin(Ke,Kv,Kn,Bag=''): # Solve n MyEq with n unknow variable
    ''' 
    Example 
        Ke=[e2,e2,e0]  MyEqs Matrix
        Kv=[N1,T,a]    unKnow Vriables
        Kn=['N_1','T','a_c']  New Name Eq
        
        N11,T1,ac = MyEqSolveLin(Ke,Kv,Kn)
        returns resepective  answer  
    '''    
    vecs=[]                 
    qq=len(Ke)
    kres=[]
    for i in range (qq):
        ee=Ke[i]
        ksym=Kv[i]
        ks=ee.solve(ksym)
        if type(ks)==list:
            rr=max(ks)
            ks=rr 
             
        vecs.append(ks)
        Ker=Ke[i+1::]
        for e1 in Ker:
            e1.set(ksym,ks,kshow=False)
            e1.reduFac(kshow=False)
            e1.simplify(kshow=False)
        
    for i ,kname in zip( vecs,Kn):
        ee=MyEq(i,kname,kshow=False)
        kres.append(ee)
    ueq=kres[-1]
    ksym=ueq()
    vsym=Kv[-1]
    for ee in kres[0:-1]:
         
        ee.set(vsym,ksym,kshow=False)
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

def Diff(ksym,kvar,kname=''):
    kres=ksym
    kres=kres.diff(kvar)
    if kname=='':
        return kres 
    else:
        return MyEq(kres,kname)
def Diff2(ksym,kvar,kname=''):
    kres=ksym
    kres=kres.diff(kvar)
    kres=kres.diff(kvar)
    if kname=='':
        return kres 
    else:
        return MyEq(kres,kname)        
class MyAnswEq:
    def __init__(self,e1,setkope=''):
        self.ee=e1
        self.kope=setkope


    def eA(self):
        kres=self.ee.v
        kope=self.kope
        try:
            return opemat(kres,kope=kope)
        except:
            return kres        
for i in range(1,9):
    kname='Q'+str(i)
    globals()[kname] = MyEq(0,kname,kshow=False,ktype='G',dtype=0)  




 

            

def upBag2sys(vecEq,kBag):
    for i in vecEq:
        i.upBag(kBag)
        
# Centro de Gravedad X Y  en funcion de dos MyEq y segmentada en x1,x2

def Cx_2func(ee1,ee2,kname,x1,x2,kvar2,kope=''):
    '''
        Centro de Gravedad X en funcion de dos MyEq y segmentada en x1,x2
        
        ee1=lower MyEq
        ee2=upper MyEq
        kname= name of new MyEq
        kname='' return Value
        x1,x2 = limits
        kope=simplify
    '''    
    A=MyEq(ee2()-ee1(),x1=x1,x2=x2,var2=kvar2,ktype='I',kshow=False)
    A.doit(kshow=False)
    Cx=MyEq(kvar2*(ee2()-ee1()),x1=x1,x2=x2,var2=kvar2,ktype='I',kshow=False)
    Cx.doit(kshow=False)
    kres=Cx()/A()
    kres=opemat(kres,kope=kope)
    if kname=='':
        return kres
    ee=MyEq(kres,kname=kname)
    return ee
    
def Cy_2func(ee1,ee2,kname='',x1='',x2='',kvar2='',kope=''):
    '''
        Centro de Gravedad Y en funcion de dos MyEq y segmentada en x1,x2
        
        ee1=lower MyEq
        ee2=upper MyEq
        kname= name of new MyEq
        kname='' return Value
        x1,x2 = limits
        kope=simplify
    '''    
    A=MyEq(ee2()-ee1(),x1=x1,x2=x2,var2=kvar2,ktype='I',kshow=False)
    A.doit(kshow=False)
    Cy=MyEq((ee2()+ee1())/2*(ee2()-ee1()),x1=x1,x2=x2,var2=kvar2,ktype='I',kshow=False)
    Cy.doit(kshow=False)
    kres=Cy()/A()
    kres=opemat(kres,kope=kope)
    if kname=='':
        return kres
     
    ee=MyEq(kres,kname=kname)
    return ee


def superset(*args):
    Kargs=[]
    for i in args:
        if type(i)==list:
            for j in i:
                Kargs.append(j)
        else:
            Kargs.append(i)

    
    mee=[]
    msym=[]
    mval=[]
    isSym=True
    for i in Kargs:
        if type(i)==MyEq or type(i)==eQ:
            mee.append(i)
        else:
            if isSym:
                msym.append(i)
                isSym=False
            else:
                mval.append(i)
                isSym=True
    kres=[]
    for j,k in zip(msym,mval):
        kres.append(MyEq(k,str(j),kshow=False))
    for i in mee:
        for j,k in zip(msym,mval):
            i.set(j,k,kshow=False)
    for i in mee:
        for j,k in zip(msym,mval):
            i.set(j,k,kshow=False)           
     
    for i in mee:
        if i.ksym!=0:
            kres.append(i) 
    for i in kres:
        i.s()
    return(kres)
def supersetVec(kVec):
    mee=[]
    msym=[]
    mval=[]
    isSym=True
 
    for i in kVec:
        if type(i)==MyEq or type(i)==eQ:
            mee.append(i)
        else:
            if isSym:
                msym.append(i)
                isSym=False
            else:
                mval.append(i)
                isSym=True
    for i in mee:
        for j,k in zip(msym,mval):
            i.set(j,k,kshow=False)
    kres=[]
    for i in mee:
        kres.append(i) 
        if i.ksym!=0:
            i.s()            
    return kres
    
    return(kres)
    
def set3(*args):
    vee=[]
    veee=[]
    for i in args:
        if type(i)==str:
            kname=i
        elif type(i)==MyEq:
            vee.append(i)

        else:
            ksym=i
    ee=vee[0]
    for j in range(1,len(vee)):
        veee.append(vee[j])
    kres=ee.solve(ksym)
    for i in veee:
        i.set(ksym,kres)
        
    return MyEq(kres,kname) 
    
def setset(*args):
    mee=[]
    msym=[]
    mval=[]
    isSym=True
 
    for i in args:
         
            
            
        if type(i)==MyEq or type(i)==eQ:
            mee.append(i)
        elif type(i)==list:
            for j in i:
                mee.append(j)
        else:
            if isSym:
                msym.append(i)
                isSym=False
            else:
                mval.append(i)
                isSym=True
    for i in mee:
        for j,k in zip(msym,mval):
            i.set(j,k,kshow=False)
     
    for i in mee:
        i.s()            
         

class eQ(MyEq):
    def __init__(self, ksym, kname='', kp=False, kope='', kshow=True,ktype='P',xx='',kvar2='',var2=t,dtype=1,depen=False,x1='',x2=''):
        self.t=ktype
        self.type=ktype
        self.depen=depen
        self.eeq=ksym
        self.primi=''
        self.xx=xx
        self.EqDiff=''
        self.EqInte=''
        self.var2=var2        
        self.iniG=0
        self.odeQ=''
        self.ode=''
        self.oldprimi=''
        if type(ksym)==MyEq:
            seq=str(ksym())
            self.ksym = unisymbols(opemat(ksym(), kope=kope))
            self.v = unisymbols(opemat(ksym(), kope=kope))
            self.eeq=ksym
        elif type(ksym)==MyEqMat:
            self.ksym = unisymbols(opemat(ksym(), kope=kope))
            
        else:    
            seq=str(ksym)
            self.ksym = unisymbols(opemat(ksym, kope=kope))
            self.v = unisymbols(opemat(ksym, kope=kope))
        self.name = alphaname(kname)
        
        
        

        self.histo = unisymbols(opemat(ksym, kope=kope))
        if kshow:
            if self.name == '':
                display(Math(latex(self.ksym)))
            else:
                sR = self.name + ' ='
                display(Math(sR + latex(self.ksym)))
        self.xx=xx
        self.EqDiff=''
        self.EqInte=''
        self.var2=kvar2
        if kvar2=='':
            self.var2=t
        self.iniG=0
        self.odeQ=''
        self.ode=''
        if x1!='':
            self.x1=x1
            self.x2=x2
        else:    
            self.x1=0
            self.x2=t
            if kvar2!='':
                self.x2=kvar2
        
        # self.s()
        if dtype==1:
            if self not in dataQ and self.name!='': 
                dataQ.append(self) 
    def __call__(self,*args,**kwargs):
        if len(args)==1:
            if self.ksym==args[0]:
                self.s()
                return
        if self.ksym=='' and len(args)==0:
            print ('si')
            return
              
        elif self.ksym=='' and len(args)==1:
            self.ksym=args[0]
            self.s()
            return        
        elif len(args)==1 and  args[0]==self.ksym:
            self.s()
        else:
            if len(args)==1 and args[0]=='cls':
                self.ksym=''
                return
            
            if len(args)==1 and len(kwargs)==0:
                var2=self.var2
                ksym=self.ksym
                return ksym.subs(var2,args[0])

            if len(args)==1 and len(kwargs)==0 and self.iniG==0 and self.t=='Q':
                self.ksym=args[0]
                self.iniG=1
                self.s()
                return
            if len(args)==1 and len(kwargs)==0 and self.t!='Q':
                ksym=self.ksym
                var2=self.var2
                nvar=args[0]
                ee=MyEq(ksym,kname=self.name,kshow=False)
                ee.set(var2,nvar)

                return    
            if len(args)==1 and len(kwargs)==0 and self.iniG==1 and args[0]=='' and self.t=='Q':

                self.ksym=0
                self.iniG=0
                #self.s()
                return    
            
            if len(args)>0 and len(kwargs)==0:
                kvar=self.var2
                if kvar!='':
                    kres=self.ksym
                    kname=self.name
                    ee=MyEq(kres,kname,kshow=False)
                    if type(kvar)!=list:
                        ee.set(kvar,args[0])
                    else:
                        for i in range(len(kvar)):
                            ee.set(kvar[i],args[i],kshow=False)
                        ee.s()    



            if len(kwargs)>0:
                vv,vs=kunpakDic(kwargs)
                mmv=[]
                mms=[]
                kope=''
                for i,j in zip(vv,vs):
                    if i=='kope':
                        kope=j 
                    else:    
                        mmv.append(parse_expr(i))
                        mms.append(j)
                 
                kres=self.ksym
                for i,j in zip(mmv,mms):
                    kres=kres.subs(i,j)
                return kres    
                

        return self.ksym

        
def reset_eq(*args):
    for i in args:
        i.ksym=''
        
        
e1=eQ(ksym='',kname='e1',kshow=False,ktype='C')   
e2=eQ(ksym='',kname='e2',kshow=False,ktype='C')  
e3=eQ(ksym='',kname='e3',kshow=False,ktype='C')  
e4=eQ(ksym='',kname='e4',kshow=False,ktype='C')  
e5=eQ(ksym='',kname='e5',kshow=False,ktype='C')   
e6=eQ(ksym='',kname='e6',kshow=False,ktype='C')   
e7=eQ(ksym='',kname='e7',kshow=False,ktype='C') 
e8=eQ(ksym='',kname='e8',kshow=False,ktype='C')
e9=eQ(ksym='',kname='e9',kshow=False,ktype='C')  
e10=eQ(ksym='',kname='e10',kshow=False,ktype='C')   
e11=eQ(ksym='',kname='e11',kshow=False,ktype='C')  
e12=eQ(ksym='',kname='e12',kshow=False,ktype='C')

def eQSolver(*args):
    vec1=[]
    uk1=[]
    for i in args:
        if type(i)==list:
            for j in i:
                if type(j)==MyEq:
                    vec1.append(j())
                elif fpoly(j,'n')>1:
                    vec1.append(j)
                else:
                    uk1.append(j)
        else:
            if type(i)==MyEq:
                vec1.append(i())
            elif fpoly(i,'n')>1:
                vec1.append(i)
            else:
                uk1.append(i)
            
    
    vec2=[]
    kres=[]
    for i in vec1:
        if type(i)==MyEq:
            vec2.append(i())
        else:
            vec2.append(i)
        
    mm=solve(vec2,uk1)
    if type(mm)==dict:
        kk,vv=kunpakDic(mm)
         
        for i,j in zip(kk,vv):
            kres.append(MyEq(j,i))
        return kres  
    else:  
        for i,j in zip(mm[0],uk1):
            j=MyEq(i,str(j))
            kres.append(j)
        return(kres)
    
def solvelin(*args,kope='',Eq=False):  #solveLinearSys(e1,e2,mu1,mu2)
        mS=[]
        mV=[]
        
        for i in args:
            if type(i)==MyEq:  
                mS.append(i())
            elif type(i)==eQ:
                ee=MyEq(i.ksym,kname=i.name,kshow=False)
                mS.append(ee())
            elif type(i)==str:
                kope=i
            else:
                mV.append(i)
        kres=solve(mS,mV)
        
        kk,vv=kunpakDic(kres)
        if kope!='':
            vv=opemat_vec(vv,kope)
        if Eq:
            EqR=[]
            for i,j in zip(kk,vv):
                EqR.append(MyEq(j,i))
            return EqR    
            
        else:
            for i,j in zip(kk,vv):
                sE([i,'=',opemat(j,kope=kope)])
        if Eq:
            kres=[]
            for i,j in zip(kk,vv):
                kres.append (MyEq(opemat(j,kope=kope),i,kshow=False))
            return kres
             
            
            
        return vv    
def get_squareMono(ksym):
    if type(ksym)==MyEq:
        ksym=ksym.ksym
    kres=ksym 
    mm=fpoly(ksym,'list')
    mr=[]
    ms=[]
    rr=[]
    centra=0
    ksigno =1
    for i in mm:
        mr.append(opemat(rpow(i,2),'r'))
        ms.append(str(opemat(rpow(i,2),'r')))
    for i,j,k  in zip(ms,mr,mm):
        if 'sqrt' in i:
            central=k
            if '-' in str(central):
                ksigno=-1
        else:
            rr.append(j)
    if len(rr)==2:
        kres=kpow(rr[1]+ksigno*rr[0],2)
    return kres        
#######
def expand2MyEq(ee):
    ktype=ee.type
    var2=ee.var2
    mm=ee.list()
    cc=1
    kname=ee.name
    kres=[]
    for i in mm:
        nname=kname+str(cc)
        nname=MyEq(i,nname,var2=var2)
        kres.append(nname)
        cc+=1
    return kres   

def upgrade(*args,kshow=True,andsolve=[]):
    if andsolve!=[]:
        vv=andsolve[0]
        ee=andsolve[1]
        vv=ee.solve(parse_expr(vv),vv)
    eev=[]
    evv=[]
    for i in args:
        if type(i)==MyEq: 
            if i.type=='C':
                eev.append(i)
            else:
                evv.append(i)
            
    for i in eev:
        for j in evv:
            i.upgrade(j,kshow=False)
    for i in eev:
        if i.ksym!=0:
            i.s()
    for i in evv:
        if type(i)==MyEq:
            i.s()
    