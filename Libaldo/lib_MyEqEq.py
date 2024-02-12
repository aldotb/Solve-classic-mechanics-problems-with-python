from sympy import *
 
from IPython.display import Math  # ,display
from matplotlib.pyplot import ylabel, plot, show, xlabel, title
from libaldo_math2 import *
from libaldo_algorith import *
 
from lib_toolsEq import *
from lib_tools import * 
from mathbasic import *
from lib_MyEq import * 
 
 
 
import copy 
x,t,k,y=symbols('x t k y')

inesymbol=['<','<=','=','>=','>']
inesymbol2=['<','≤','=','≥','>']
C1,C2,C3,C4,C5=symbols('C1 C2 C3 C4 C5')
def creteInt(expr,var,x1='',x2=''):
    if var!='':
        if x1=='': 
            kres=Integral(expr,var)
        else:
            kres=Integral(expr, (var,x1,x2))
        return kres
    else:
        return expr


def MQ(*args, var=x,var1=y,var2=z,kshow=True,ktype='eq',render=True):
    return MyEqEq(*args, var=var,var1=var1,var2=var2,kshow=kshow,ktype=ktype,render=render)
class MyEqEq:
    def __init__(self, *args, var=x,var1=y,var2=z,kshow=True,ktype='eq',render=True):
        '''
        self.insys=''
        self.clone=''
        self.cloneb=''
        self.var=var
        self.type=ktype
        self.symb="="
        P1=args[1]
        self.CC=0         
        self.vecCC=[C1,C2,C3,C4,C5] 
        self.e1=MyEq(0,'e1',var=var,kshow=False)
        self.e2=MyEq(0,'e2',var=var,kshow=False)
        ineqtype=notaclass=[GreaterThan,LessThan,StrictGreaterThan,StrictLessThan]
        ineqsynbo=['>=','<=','>','<']
        if  type(P1) in ineqtype:  # input one math expr tpe  x+3>8
            ksym1,ksym2=P1.args
            self.psymb=ineqsynbo[ineqtype.index(type(P1))]
            self.symb=nice_iq_symbol(self.psymb)
            self.e1.ksym=ksym1
            self.e2.ksym=ksym2
            self.type='IQ'            
        '''     
        
        
        self.insys=''
        self.clone=''
        self.cloneb=''
        self.var=var
        self.type=ktype
        self.symb="="
        P1=args[1]
         
        self.CC=0         
        self.vecCC=[C1,C2,C3,C4,C5] 
        self.e1=MyEq(0,'e1',var=x,kshow=False)
        self.e2=MyEq(0,'e2',var=x,kshow=False)
        ineqtype=notaclass=[GreaterThan,LessThan,StrictGreaterThan,StrictLessThan]
        ineqsynbo=['>=','<=','>','<']
        if  P1 in ineqsynbo:  # input one math expr tpe  x+3>8
             
            ksym1,ksym2=args[0],args[2]
            self.psymb=args[1]
            self.symb=nice_iq_symbol(self.psymb)
            self.e1.ksym=ksym1
            self.e2.ksym=ksym2
            self.type='IQ'     
        elif len(args)==3:  # input tree  args tpe  x+3,'>',8
            ksym1,ksym2=args[0],args[2]
            if type(ksym1)==MyEq:
                ksym1=ksym1.ksym
            if type(ksym2)==MyEq:
                ksym2=ksym2.ksym
            self.psymb=args[1]
            self.symb=nice_iq_symbol(self.psymb)
            self.e1.ksym=ksym1
            self.e2.ksym=ksym2
            self.type='IQ'
             
         
        elif len(args)==2:
        
            ksym1,ksym2=args[0],args[1]
            if type(ksym1)==MyEq:
                ksym1=ksym1.ksym
            if type(ksym2)==MyEq:
                ksym2=ksym2.ksym
                
                
            self.psymb='='
            self.symb='='
            self.e1.ksym=ksym1
            self.e2.ksym=ksym2
            self.type='EQ'
             
        elif len(args)==4:
        
            ksym1,ksym2=args[0],args[1]
            if type(ksym1)==MyEq:
                ksym1=ksym1.ksym
            if type(ksym2)==MyEq:
                ksym2=ksym2.ksym
                
                
            self.psymb='='
            self.symb='='
            self.e1.ksym=ksym1
            self.e1.var=args[2]
            self.e2.ksym=ksym2
            self.e2.var=args[3]
            self.type='EQ'    
        
    
        p1=self.e1.ksym
        p2=self.e2.ksym
        ps=self.symb
        self.render=render
        if "'" in str(p1) or "'" in str(p2):
            self.type='DI'
        if render:
            if kshow:
                display(Math(latex(p1)+' '+ps+' '+latex(p2)))
        
        self.lundo=[]
        self.lundo.append([self.e1.ksym,self.e2.ksym])    
    def __call__(self,*args, **kwargs):
        if len(args)==0 and len(kwargs)>0:
             
            p1=self.L
            p1=real_subs(p1,**kwargs)
             
            p2=self.R 
            p2=real_subs(p2,**kwargs)
            QQ=MQ(p1,p2) 
             
            return
        kfloat=False
        for i in args:
            if i=='float':
                kfloat=True
        QQ=self.xcopy()
         
        p1=QQ.left
        p2=QQ.right
        p11=subskwargs(p1,**kwargs)
        p22=subskwargs(p2,**kwargs)
        if kfloat:
            try:
                p11=float(p11)
            except:
                pass
            try:
                p22=float(p22)
            except:
                pass
            
        QQ.e1.ksym=p11
        QQ.e2.ksym=p22
        QQ.s()
    def __repr__(self):
        kres = self.eQ
        return kres
        
    def _latex(self, obj):
        return latex(self.eQ)    
    def __str__(self):
         
        return str(self.__repr__())
    
    def up_sides(self,p1='',p2=''):
        self.up_ee(p1=p1,p2=p2)
        
    def up_ee(self,p1='',p2=''):
        if p1!='':
            self.e1.ksym=p1
        if p2!='':
            self.e2.ksym=p2
    def unisymbols(self):
        self.e1.ksym=unisymbols(self.e1.ksym)
        self.e2.ksym=unisymbols(self.e2.ksym)
    ###########################################
    #               Update                    #
    ###########################################

    def __add__(self, other):
        p1=self.e1.ksym
        p2=self.e2.ksym
        if type(other)==MyEqEq:
            p1=p1+other.left
            p2=p2+other.right
        elif type(other)==MyEq:
            p1=p1+other.ksym
            p2=p2+other.ksym

        else:
            p1=p1+other 
            p2=p2+other 
         
        return Eq(p1,p2)

    def __radd__(self, other):
        p1=self.e1.ksym
        p2=self.e2.ksym
        if type(other)==MyEqEq:
            p1=p1+other.left
            p2=p2+other.right
        elif type(other)==MyEq:
            p1=p1+other.ksym
            p2=p2+other.ksym

        else:
            p1=p1+other 
            p2=p2+other 
         
        return Eq(p1,p2)

    def __sub__(self, other):
        p1=self.e1.ksym
        p2=self.e2.ksym
        if type(other)==MyEqEq:
            p1=p1-other.left
            p2=p2-other.right
        elif type(other)==MyEq:
            p1=p1-other.ksym
            p2=p2-other.ksym

        else:
            p1=p1-other 
            p2=p2-other 
         
        return Eq(p1,p2)

    def __rsub__(self, other):
        p1=self.e1.ksym
        p2=self.e2.ksym
        if type(other)==MyEqEq:
            p1=p1-other.left
            p2=p2-other.right
        elif type(other)==MyEq:
            p1=p1-other.ksym
            p2=p2-other.ksym

        else:
            p1=p1-other 
            p2=p2-other 
         
        return Eq(p1,p2)
    
    def __mul__(self, other):
        p1=self.e1.ksym
        p2=self.e2.ksym
        if type(other)==MyEqEq:
            p1=p1*other.left
            p2=p2*other.right
        elif type(other)==MyEq:
            p1=p1*other.ksym
            p2=p2*other.ksym

        else:
            p1=p1*other 
            p2=p2*other 
         
        return Eq(p1,p2)
    
    def __rmul__(self, other):
        p1=self.e1.ksym
        p2=self.e2.ksym
        if type(other)==MyEqEq:
            p1=p1*other.left
            p2=p2*other.right
        elif type(other)==MyEq:
            p1=p1*other.ksym
            p2=p2*other.ksym

        else:
            p1=p1*other 
            p2=p2*other 
         
        return Eq(p1,p2)
        
        
    def __truediv__(self, other):
        p1=self.e1.ksym
        p2=self.e2.ksym
        if type(other)==MyEqEq:
            p1=p1/other.left
            p2=p2/other.right
        elif type(other)==MyEq:
            p1=p1/other.ksym
            p2=p2/other.ksym

        else:
            p1=p1/other 
            p2=p2/other 
         
        return Eq(p1,p2)

    def __rtruediv__(self, other):
        p1=self.e1.ksym
        p2=self.e2.ksym
        if type(other)==MyEqEq:
            p1=p1/other.left
            p2=p2/other.right
        elif type(other)==MyEq:
            p1=p1/other.ksym
            p2=p2/other.ksym

        else:
            p1=p1/other 
            p2=p2/other 
         
        return Eq(p1,p2)        
    
    #               Update                    #
    ########################################### 
    
    def add_ics(self,var1,var2,var3):
        Vics=self.Vics
        
         
        if 'd2' in var1:
             
            arma=d2ficc(self.vmain,self.var2,var2,var3)
            Vics.append(arma)
        elif  'd' in var1 and 'd2' not in var1 :
             
            arma=dficc(self.vmain,self.var2,var2,var3)
            Vics.append(arma)
        else:
            v1=str(self.vmain)+'('+str(var2)+')'
            arma=  fics(v1,var3)
            Vics.append(arma)
        self.Vics=Vics
        self.arma_ics()
        sE(['ics=',self.ics])
    
    def arma_ics(self):
        kres='{'
        Vics=self.Vics
        for i in Vics:
            kres+=i+','
        kres=kres[0:-1]
        kres+='}'
        self.ics=parse_expr(kres)
    ################################## 
    #      Apariencia
    ################################## 
    def s2(self):
    
        self.up_eQshow()
        func=self.func
        dfunc=self.dfunc
        qq=len(func)
        kres=self.eQshow
        for i in range(qq):
            kfunc=func[i]
            kdfunc=dfunc[i]
            kres=kres.subs(kfunc[0],kdfunc[0])
            #res=kres.subs(kfunc[0],kdfunc[0])
        self.eQshow=kres
        self.s('2')
        
    def s3(self):
         
        func=self.func
        dfunc=self.dfunc
        qq=len(func)
        kres=self.eQshow
        for i in range(qq):
            kfunc=func[i]
            kdfunc=dfunc[i]
            kres=kres.subs(kfunc[1],kdfunc[1])
            #res=kres.subs(kfunc[0],kdfunc[0])
        self.eQshow=kres
        self.s('2')    
    
    def save_eQshow(self,keQ):
        self.eQshow=keQ
        
    def up_eQshow(self):
        self.eQshow=self.eQ
    
        
    def up2primitiva(self):
        eq=self.eQ
        self.e1.ksym=eq.lhs
        self.e2.ksym=eq.rhs
        self.modeNormal=False
        
        
    
    def up2normalMode(self):
        return self.primitiva2func() 
    
    def primitiva2func(self):
        try:
            self.e1.ksym=flat_diff(self.left,self.df)
        except:
            pass
        vmain=self.vmain
        vmain=antiprimitiva(vmain)
        nx=symbols(str(vmain))
        fmain=Function(str(vmain))(self.var2)
        if self.modeNormal:
            self.primieQ=Eq(self.e1.ksym,self.e2.ksym)
            self.modeNormal=False
        
        try:
            ee1=self.e1
            ee1.set(fmain,nx,kshow=False)
            self.e1=ee
        except:
            pass
            
        try:
            ee2=self.e2
            ee2.set(fmain,nx,kshow=False)
            self.e2=ee2
        except:
            pass
         
            
        self.s()
        
    def flat_diff(self):
        try:
            self.e1.ksym=self.e1.ksym.subs(Q.df,flat_diff(Q.df))
        except:
            pass
        self.s()
    def updateEq(self):
        p1=self.e1.ksym
        p2=self.e2.ksym
        self.eQ =Eq(p1,p2)


    def updatee(self,p1,p2):
        self.e1.ksym=p1
        self.e2.ksym=p2
        self.eQ =Eq(p1,p2)
        self.exp1=p1 
        self.exp2=p2


    def update(self):

        self.eQ = Eq(sydem(self.e1.ksym), sydem(self.e2.ksym))
         
    def update2  (self,nEq):
        kres=nEq
        self.e1.ksym=kres.lhs
        self.e2.ksym=kres.rhs
        self.update() 
    
    def reformate(self,p1,p2,kshow=True):
        self.e1.ksym=p1
        self.e2.ksym=e2
        self.s(kshow=kshow)
        
    ################################## 
    #         Show()
    ##################################
    
    def set_vista(self,ktype):
        self.vista=ktype
    def swap(self):
        ee=self.e2
        self.e2=self.e1
        self.e1=ee
        self.s()
    def sshow(self,kshow):
        if kshow:
            self.s()
            
    def s(self,kshow=True,kclone=True):
        if self.render: 
            if self.insys!='':
                self.insys.s()
            else:
                p1=self.e1.ksym
                p2=self.e2.ksym
                ps=self.symb
                 
                if kshow:
                    display(Math(latex(p1)+' '+ps+' '+latex(p2)))
        P1,P2=self.lundo[-1]
        if P1!=p1 or P2!=p2:
            self.lundo.append([self.L,self.R])
     
    def xcopy(self,op=''):
        kres=copy.deepcopy(self)
        if op!='':
            kres.s()
            
        return kres

    @property
    def rval(self,**args):  # return lhs from Eq
        return self.e2.ksym

    @property
    def right(self):  # return lhs from Eq
        return self.e2.ksym    
    
    @property
    def R(self):
        return unisymbols(self.right)
        
    @property
    def L(self,val=''):
        if val=='':    
            return unisymbols(self.left)
        else:
            ee=unisymbols(self.e1)
            ee.var=unisymbols(self.var)
            
            return ee(val)
            
         
    def setvarL(self,var):
        self.e1.var=var
    
        
    def setvarR(self,var):
        self.e2.var=var
    
    def all2L(self,kshow=True):
        p1=self.L
        p2=self.R
        self.e1.ksym=simplify(p1-p2)
        self.e2.ksym=0
        if kshow:
            self.s()
    def all2R(self,kshow=True):
        p1=self.L
        p2=self.R
        self.e2.ksym=simplify(p2-p1)
        self.e1.ksym=0
        if kshow:
            self.s()   
        
    @property
    def u(self):
        self.undo()
        
    @property
    def lval(self):  # return lhs from Eq
        return self.e1.ksym
    @property
    def left(self):  # return lhs from Eq
        return self.e1.ksym    

    def getdata(self):
        return self.e1.ksym,self.e2.ksym
    
    def dothis(self,*args):
        op='LR'
        args3=[]
        for i in args:
            if i!='L' and i!='R':
                args3.append(i)
            else:
                op=i
        if 'L' in op:
            ksym=self.L
            args2=[ksym]
            for i in args3:
                args2.append(i)
            kres=dothis(*args2)
            self.e1.ksym=kres
        
        if 'R' in op:
            ksym=self.R
            args2=[ksym]
            for i in args3:
                args2.append(i)
            kres=dothis(*args2)
            self.e2.ksym=kres       
    
        
        self.s()
    def doindenom(self,*args):
        op=''
        vecL=[self.L]
        vecR=[self.R]
        for i in args:
            if i=='L' or i=='R':
                op=i 
            else:
                vecL.append(i)
                vecR.append(i)
        if op=='L':
            self.e1.ksym=doindenom(*vecL)
        elif op=='R':
            self.e2.ksym=doindenom(*vecR)
        else:
            self.e1.ksym=doindenom(*vecL)
            self.e2.ksym=doindenom(*vecR)
        self.s()    

    def doinnumer(self,*args):
        op=''
        vecL=[self.L]
        vecR=[self.R]
        for i in args:
            if i=='L' or i=='R':
                op=i 
            else:
                vecL.append(i)
                vecR.append(i)
        if op=='L':
            self.e1.ksym=doinnumer(*vecL)
        elif op=='R':
            self.e2.ksym=doinnumer(*vecR)
        else:
            self.e1.ksym=doinnumer(*vecL)
            self.e2.ksym=doinnumer(*vecR)
        self.s() 
    
    def updata(self,p1,p2):
        self.e1.ksym=p1
        self.e2.ksym=p2
        
    
    def setL(self,ksym='',kshow=True):
        '''
            if ksym!='' then self.left will be ksym
            else try to replace kwars in left
        '''
        
        if ksym!='':
            self.e1.ksym=ksym 

        if kshow:
            self.s()
    def setR(self,ksym='',kshow=True):
        '''
            if ksym!='' then self.left will be ksym
            else try to replace kwars in left
        '''
        
        if ksym!='':
            self.e2.ksym=ksym 

        if kshow:
            self.s()   
    ################################## 
    #         set()
    ##################################
    def validatesymbols(self,svar):
        self.e1.ksym=expr2var(self.L,svar)
        self.e2.ksym=expr2var(self.R,svar)
        self.s()
        
    def set_eQshow(self,*args,kshow=True):
        qq=len(args)
        ks1=[]
        ks2=[]
        
        for i in range(qq):
            if i%2==0:
                ks1.append(args[i])
            else:
                ks2.append(args[i])
        
        eQshow=self.eQshow
        for i,j in zip(ks1,ks2):
            eQshow=eQshow.subs(i,j)
        self.eQshow=eQshow    
        
        if kshow:
            self.s('2')
        
        
    def changesigno(self,kshow=True):
        self.e1.ksym=-1*self.e1.ksym
        self.e2.ksym=-1*self.e2.ksym
        if kshow:
            self.s()
    
    def evalue_if(self,*args,kshow=True,kope='',**kwargs):
        QQ=self.xcopy('QQ')
        margs=args
        mkwargs=kwargs
        QQ.set(*margs,kshow=kshow,kope=kope,**mkwargs)
    
    def evalueif(self,*args,kshow=True,kope='',**kwargs):
        QQ=self.xcopy('QQ')
        margs=args
        mkwargs=kwargs
        QQ.set(*margs,kshow=kshow,kope=kope,**mkwargs)
    
    # def set(self,*args,kshow=True,**kwargs):
    
         
        # self.e1.set(*args,kshow=False,**kwargs)
        # self.e2.set(*args,kshow=False,**kwargs)
          
    
           
        # self.s()    
    # def set(self,**kwargs):
        # p1=self.L
        # p1=setpack(p1,**kwargs)
        # p2=self.R
        # p2=setpack(p2,**kwargs)
        # self.e1.ksym=p1
        # self.e2.ksym=p2
        # self.s()
    def set(self,*args,**kwargs):
        p1=self.L
        p2=self.R
        if len(kwargs)>0:
            p1=realsub2(p1,**kwargs)
            p2=realsub2(p2,**kwargs)
           
        if len(args)>0:
            for i in args:
                if type(i)==MyEq:
                    kres=i.ksym
                    svar=i.name
                    try:
                        sp1=str(p1)
                        p1=p1.subs(svar,kres)
                        if str(p1)==sp1:
                            p1=subsubs(p1,str(svar),kres)
                    except:
                        pass
                    try:
                        sp2=str(p2)
                        p2=p2.subs(svar,kres)
                        if str(p2)==sp2:
                            p2=subsubs(p2,str(svar),kres)
                    except:
                        pass
        self.e1.ksym=p1
        self.e2.ksym=p2
        self.s()
        
        
    def undo(self,kshow=True):
        qq=len(self.lundo)
        if qq>1:
            qu=qq-1
            kres=self.lundo[0:qu]
            p1,p2=kres[-1]
            self.e1.ksym=p1
            self.e2.ksym=p2
            ps=self.symb
            self.lundo=kres
            
            if kshow:
                display(Math(latex(p1)+' '+ps+' '+latex(p2)))
        else:
            self.s()
        
        
    ################################## 
    #         simplify
    ################################## 

    

    def addexpand(self):
        self.e1.ksym=addexpand(self.e1.ksym)
        self.e2.ksym=addexpand(self.e2.ksym)
        self.s()
        
        
    def expand(self, kop='RL',kshow=True):
        if 'L' in kop:
            self.e1.expand(kshow=False)
        if 'R' in kop:
            self.e2.expand(kshow=False)
        if kshow:
            self.s()
    
    
        
     
        
    def Expand(self,kshow=True,kope=''):
        p1=self.left
        p2=self.right
        p1=factor(p1)
        if Is_Mono(p1):
            p2=p2*denom(p1)
            p2=factor(p2)
            p1=numer(p1)
            
        if Is_Mono(p2):
            p1=p1*denom(p2)
            p2=numer(p2)
            
        p1=opemat(p1,kope=kope)
        p2=opemat(p2,kope=kope)
        self.e1.ksym=p1
        self.e2.ksym=p2
        self.s()
            
        
    def opemat(self,kope='',op='LR'):
        p1=self.left
        p2=self.right
        
        if 'L' in op:
            p1=opemat(p1,kope=kope)
        if 'R' in op:
            p2=opemat(p2,kope=kope)    
        self.e1.ksym=p1
        self.e2.ksym=p2
        self.s()    
    def getpart(self,*args):
        return self.get_args(*args)

    def args(self,*args,deep=2,format='list'):
        if len(args)==0:
            showarglist(self.L,deep=deep,format=format,side='L')
            showarglist(self.R,deep=deep,format=format,side='R')
        else: 
            kres=Eq(self.e1.ksym,self.e2.ksym)
            for i in args:
                kres=kres.args[i]
            return kres
        
    def get_args(self,*args):
        self.eQ=Eq(self.L,self.R)
        kres=self.eQ
        for i in args:
            kres=kres.args[i]
        return kres
    
    def kargs(self,*args):
        kres=self.eQ
        for i in args:
            kres=kres.args[i]
        return kres
        
    def exp_alone(self):
        eQs=str(self.eQ)
        sres=''
        try:
            smm=between_par(eQs,'exp(')
        except:
            self.s()
            return
        smm='exp('+smm+')'
        try:
            mm=parse_expr(smm)
        except:
            self.s()
            return
        self.alone(mm)
    
    def nformat(self,nd):
        self.e1.ksym=nformat(self.e1.ksym,nd)
        self.e2.ksym=nformat(self.e2.ksym,nd)
        self.s()
    
    def reducecero(self):
        p1=self.L 
        p2=self.R 
        if p1==0:
            p2=reducecero(p2)
        else:
            if p2==0:
                p1=reducecero(p1)
        self.e1.ksym=p1        
        self.e2.ksym=p2
        self.s()
        
        
        
    def reduce(self,*args,kshow=True):
    
        p1=self.e1.ksym
        p2=self.e2.ksym
        if Is_Pow(self.L) and Is_Pow(self.L):
            base1=getbase(p1)
            base2=getbase(p2)
            ee1=getexpo(p1)
            ee2=getexpo(p2)
            if base1==base2:
                self.e1.ksym=ee1
                self.e2.ksym=ee2
            elif ee1==ee2:
                self.e1.ksym=base1
                self.e2.ksym=base2
            elif ee1==base1 and ee2==base2:
                self.e1.ksym=base1
                self.e2.ksym=base2
            elif ee1==base1:  
                self.e1.ksym=base2
                self.e2.ksym=ee2    
            elif ee2==base2: 
                self.e1.ksym=base1
                self.e2.ksym=ee1
            else:
                self.s()
                return
        else:    

        
            if type(p1)==log and type(p2)==log:
                self.e1.ksym=p1.args[0]
                self.e2.ksym=p2.args[0]
            elif type(p1)==Pow and type(p2)==Pow:    
                    b1=getbase(p1)
                    b2=getbase(p2)
                    e1=getexpo(p1)
                    e2=getexpo(p2)
                    if b1==b2:
                        self.e1.ksym=e1
                        self.e2.ksym=e2
                    if e1==e2:
                        self.e1.ksym=b1
                        self.e2.ksym=b2
            else:            
                p1=self.e1.ksym
                p2=self.e2.ksym    
                p1=p1.factor()
                p1=unisymbols(p1)
                p2=self.e2.ksym
                p2=p2.factor()
                p2=unisymbols(p2)
                if len(args)>0:
                    for i in args:
                        p1=p1/i
                        p2=p2/i
                    self.e1.ksym=p1
                    self.e2.ksym=p2
                else: 
                    if Is_Mul(p1) and Is_Mul(p2):
                        m1=fpoly(p1,'list')
                        m2=fpoly(p2,'list')
                        L1=subtrac_list(m1,m2)
                        pp1=listmul(L1)
                        L2=subtrac_list(m2,m1)
                        pp2=listmul(L2)
                        self.e1.ksym=pp1
                        self.e2.ksym=pp2
                     
                 
        if kshow:
            self.s() 
            
    
    def simplifySum(self):
        p1=self.e1.ksym 
        p2=self.e2.ksym
        try:
            p1,p2=simplifySum(p1,p2)
        except:
            pass
        self.e1.ksym=p1
        self.e2.ksym=p2
        self.s()
    
    def simplify(self,*args, kop='RL' ,force=False,kshow=True,ssum=True,smul=True):
        

        p1=self.e1.ksym 
        p2=self.e2.ksym
            
        if self.type=='iq3':
            self.e1.ksym=simplify(p1)
            self.e2.ksym=simplify(p2)
            if kshow:
                self.s()
            return    
        
        if Is_Div(p1) and Is_Div(p2):
            if denom(p1)==denom(p2):
                self.e1.ksym=numer(p1)
                self.e2.ksym=numer(p2)
                if kshow:
                    self.s()
                return
        
        if Is_Pow(p1) and Is_Pow(p2):
            if getbase(p1)==getbase(p2):
                self.e1.ksym=getexpo(p1)
                self.e2.ksym=getexpo(p2)
                if kshow:
                    self.s()
                return
            elif getexpo(p1)==getexpo(p2):
                self.e1.ksym=getbase(p1)
                self.e2.ksym=getbase(p2)
                if kshow:
                    self.s()
                return
           
        self.e1.ksym=simplify(p1)
        self.e2.ksym=simplify(p2) 

        if kshow:
            self.s()
    def lsimplify(self,kop='RL' ,kshow=True):
        p1=self.L
        p2=self.R
         
        p1=lsimplify(p1)
        p2=lsimplify(p2)
        self.e1.ksym=p1
        self.e2.ksym=p2

        if kshow:
            self.s()
            
    def lcombine(self,kop='RL' ,kshow=True):
        p1=self.L
        p2=self.R
         
        p1=logcombine(p1,force=True)
        p2=logcombine(p2,force=True)
        self.e1.ksym=p1
        self.e2.ksym=p2

        if kshow:
            self.s()            
            
    def expandexpo(self,kop='LR',kshow=True):
        p1=self.e1.ksym
        if 'L' in kop:
            p1=expandexpo(p1)
            self.e1.ksym=p1
        p2=self.e2.ksym
        if 'R' in kop:
            p2=expandexpo(p2)
            self.e2.ksym=p2
        if kshow:
            self.s()
            
            
    def rsimplify(self,kop='LR',kshow=True):
        p1=self.e1.ksym
        if 'L' in kop:
            if 'sqrt' in str(p1):
                p1=rsimplify(p1) 
                self.e1.ksym=p1
        p2=self.e2.ksym
        if 'R' in kop:
            if 'sqrt' in str(p2):
                p2=rsimplify(p2)
                self.e2.ksym=p2
        if kshow:
            self.s()        
    def simplify_img(self,kop='LR',kshow=True):
        p1=self.e1.ksym
        if 'L' in kop:
            p1=simplify_img(p1) 
            self.e1.ksym=p1
        p2=self.e2.ksym
        if 'R' in kop:
            p2=simplify_img(p2)
            self.e2.ksym=p2
        if kshow:
            self.s()        
    def expandbase(self,kop='LR',kshow=True):
        p1=self.e1.ksym
        if 'L' in kop:
            p1=expandbase(p1)
            self.e1.ksym=p1
        p2=self.e2.ksym
        if 'R' in kop:
            p2=expandbase(p2)
            self.e2.ksym=p2
        if kshow:
            self.s()        
    def factorbase(self,kop='LR',kshow=True):
        p1=self.e1.ksym
        if 'L' in kop:
            p1=factorbase(p1)
            self.e1.ksym=p1
        p2=self.e2.ksym
        if 'R' in kop:
            p2=factorbase(p2)
            self.e2.ksym=p2
        if kshow:
            self.s()
    def simplifybase(self,kop='LR',kshow=True):

        p1=self.left
        if 'L' in kop:
            p1=simplifybase(p1)
            self.e1.ksym=p1
        p2=self.right
        if 'R' in kop:
            p2=simplifybase(p2)
            self.e2.ksym=p2
         
        if kshow:
            self.s()
    def powexpand(self,kop='LR',kshow=True):
        '''
        input (x**(a*b))   ---->   return(x**a)**b
        input (x**(a*b),b)   ---->   return(x**b)**a
        '''
        op=''
        if 'i' in kop:
            op='i'
        p1=self.left
        if 'L' in kop:
            p1=powexpand(p1,op=op)
            self.e1.ksym=p1
        p2=self.right
        if 'R' in kop:
            p2=powexpand(p2,op=op)
            self.e2.ksym=p2
         
        if kshow:
            self.s()

    def mulexpo(self,kop='LR',kshow=True,force=False):
        p1=self.e1.ksym
        if 'L' in kop:
            p1=mulexpo(p1,force=force)
            self.e1.ksym=p1
        p2=self.e2.ksym
        if 'R' in kop:
            p2=mulexpo(p2,force=force)
            self.e2.ksym=p2
        if kshow:
            self.s()
    def pow2powpow(self,exp1='',kop='LR',kshow=True):
        if exp1!='' and type(exp1)==str:
            kop=exp1	
        p1=self.e1.ksym
        if 'L' in kop:
            p1=pow2powpow(p1,exp1)
            self.e1.ksym=p1
        p2=self.e2.ksym
        if 'R' in kop:
            p2=pow2powpow(p2,exp1)
            self.e2.ksym=p2
        if kshow:
            self.s()

    def positivexpo(self,kop='LR',kshow=True):
        p1=self.e1.ksym
        if 'L' in kop:
            p1=positivexpo(p1)
            self.e1.ksym=p1
        p2=self.e2.ksym
        if 'R' in kop:
            p2=positivexpo(p2)
            self.e2.ksym=p2
        if kshow:
            self.s()        
    def simplifyexpo(self,op='LR',kope='',kshow=True):
        p1=self.e1.ksym
        if 'L' in op:
            p1=simplifyexpo(p1)
            self.e1.ksym=p1
        p2=self.e2.ksym
        if 'R' in op:
            p2=simplifyexpo(p2)
            self.e2.ksym=p2

        if kshow:
            self.s()
            
            
    def basefactor(self,op='LR',kope='',kshow=True):
        return self.simplifyexp(op=op,kope=kope,kshow=kshow)
        
        
    def simplifyexp(self,op='LR',kope='',kshow=True):
        p1=self.left
        if 'L' in op:
            p1=simplifyexp(p1,kope=kope)
            self.e1.ksym=p1
        p2=self.right
        if 'R' in op:
            p2=simplifyexp(p2,kope=kope)
            self.e2.ksym=p2

        if kshow:
            self.s()

    def div2mul(self,ops='LR'):
        if 'L' in ops and Is_Div(self.L):
            p1=numer(self.L)
            p2=denom(self.L)
            self.e1.ksym=p1*p2
        if 'R' in ops and Is_Div(self.R):
            p1=numer(self.R)
            p2=denom(self.R)
            self.e2.ksym=p1*p2 
        self.s() 
        
    def div2mulexpo(self,op='LR',kope='',kshow=True):
        p1=self.left
        if 'L' in op:
            p1=div2mulexpo(p1)
            if kope!='':
                p1=opemat(p1,kope=kope)
            self.e1.ksym=p1
            
        p2=self.right
        if 'R' in op:
            p2=div2mulexpo(p2)
            if kope!='':
                p2=opemat(p2,kope=kope)
            self.e2.ksym=p2

        if kshow:
            self.s()

    def reducePow(self,op='LR',kope='',kshow=True):
        p1=self.left
        if 'L' in op:
            p1=reducePow(p1)
            if kope!='':
                p1=opemat(p1,kope=kope)
            self.e1.ksym=p1
            
        p2=self.right
        if 'R' in op:
            p2=reducePow(p2)
            if kope!='':
                p2=opemat(p2,kope=kope)
            self.e2.ksym=p2

        if kshow:
            self.s()        
    def simplifyrpow(self,kope='',kop='RL',kshow=True):
        

        if 'L' in kop :
            self.e1.simplifyrpow(kshow=False)
        if 'R' in kop :
            self.e2.simplifyrpow(kshow=False)
        if kope!='':
            kres1=self.e1.ksym
            kres2=self.e2.ksym
            
            kres1=opemat(kres1,kope=kope)
            kres2=opemat(fcc,kope=kope)
            self.e1.ksym=kres1
            self.e2.ksym=kres2
        self.s()
    def simplify_cero(self):
        kres=self.left-self.right
        kres=opemat(kres)
        self.e1.ksym=kres
        self.e2.ksym=0
        self.s()
    

    def primefactor(self, *args,kshow=True):
         
        op='RL' 
        p1=self.L
        p2=self.R
        
        if 'L' in args:
            op='L'
        if 'R' in args:
            op='R'
        if 'L' in op:  
            p1=primefactor(p1)
            self.e1.ksym=p1 
        if 'R' in op: 
            p2=primefactor(p2)
            self.e2.ksym=p2
         
        if kshow:
            self.s()
        
    
    def factor(self,*args, kop='RL',kshow=True):
        
        if 'L' in kop:
            self.e1.factor(*args,kshow=False)

        if 'R' in kop:
            self.e2.factor(*args,kshow=False)
        
        
             
        if kshow:
            self.s()
    def apart(self,kop='RL',kshow=True):
        
        if 'L' in kop:
            self.e1.apart(kshow=False)

        if 'R' in kop:
            self.e2.apart(kshow=False)
        
        
             
        if kshow:
            self.s()
    def inverse(self,kop='RL',kshow=True):
        
        if 'L' in kop:
            self.e1.inverse(kshow=False)

        if 'R' in kop:
            self.e2.inverse(kshow=False)
        
        
             
        if kshow:
            self.s()
            
            
    def factors(self,kfac, kop='RL',kshow=True):
        
        if 'L' in kop:
            p1=self.L
            p1=factors(p1,kfac)
            self.e1.ksym=p1

        if 'R' in kop:
            p2=self.R
            p2=factors(p2,kfac)
            self.e2.ksym=p2
        
        
             
        if kshow:
            self.s()        
        
    def dfactor(self,var,var1,op='',side='LR'):
        if 'L' in side:
            p1=dfactor(self.e1.ksym,var=var,var1=var1,op=op)
            self.e1.ksym=p1
        if 'R' in side:
            p2=dfactor(self.e2.ksym,var=var,var1=var1,op=op)
            self.e2.ksym=p2
        if kshow:
            self.s()
        
    def tsimplify(self, kop='RL',kshow=True,kope=''):
        if 'L' in kop:
            self.e1.ksym=expand_log(self.e1.ksym)
        if 'L' in kop:
            self.e2.ksym=expand_log(self.e2.ksym)
        if kshow:
            self.s()
        
    def trinom2binom(self,sexpr,kop='RL',kshow=True):
        if 'L' in kop:
            self.e1.trinom2binom(sexpr,kshow=False)
        if 'R' in kop:
            self.e2.trinom2binom(sexpr,kshow=False)
        if kshow:
            self.s()
    
    def tfactor(self, kop='RL',kshow=True,kope=''):
        if 'L' in kop:
            self.e1.tfactor(kope=kope,kshow=False)
        if 'R' in kop:
            self.e2.tfactor(kope=kope,kshow=False)
        if kshow:
            self.s()
    
    def texpand(self, kop='RL',kshow=True,kope=''):
        if 'L' in kop:
            self.e1.texpand(kope=kope,kshow=False)
        if 'R' in kop:
            self.e2.texpand(kope=kope,kshow=False)
        if kshow:
            self.s()
        

        
    # def Add(self,*args,kshow=True):
        # vecpara=['simplify','factor','expand'] 
        # kval1,kval2,kname,op=reparte_param(vecpara,*args)
         	    
        # p1=self.left
        # p2=self.right
        # if 'L' in op:
                # p1=fes(self.left+kval1,*args)
        # if 'R' in op:
                # p2=fes(self.right+kval2,*args)         
        # if kname!='':
            # QQ=MyEqEq(p1,p2)
            # return QQ
        # else:
            # self.e1.ksym=premath(p1,*args) 
            # self.e2.ksym=premath(p2,*args)
            # if kshow:
                # self.s()
        
        
        
        
    # def Substrac(self,*args,kshow=True):
         
        # vecpara=['simplify','factor','expand']
        # kval1,kval2,kname,op=reparte_param(vecpara,*args)
         	    
        # p1=self.left
        # p2=self.right
        # if 'L' in op:
                # p1=fes(self.left-kval1,*args)
        # if 'R' in op:
                # p2=fes(self.right-kval2,*args)         
        # if kname!='':
            # QQ=MyEqEq(p1,p2)
            # return QQ
        # else:
            # self.e1.ksym=premath(p1,*args) 
            # self.e2.ksym=premath(p2,*args)
            # if kshow:
                # self.s()
        
    # def Mul(self,*args,kshow=True):
         
        # vecpara=['simplify','factor','expand']
        # kval1,kval2,kname,op=reparte_param(vecpara,*args)
         	    
        # p1=self.left
        # p2=self.right
        # if 'L' in op:
                # p1=fes(self.left*kval1,*args)
                # p1=premath(p1,*args)
        # if 'R' in op:
                # p2=fes(self.right*kval2,*args)
                # p2=premath(p2,*args)                
        # if kname!='':
            # QQ=MyEqEq(p1,p2)
            # return QQ
        # else:
            # self.e1.ksym=p1
            # self.e2.ksym=p2
            # if kshow:
                # self.s()
    
    # def Div(self,*args,kshow=True,unevaluate=False):
        # vecpara=['simplify','factor','expand','nosimplify']
        # kval1,kval2,kname,op=reparte_param(vecpara,*args) 
         
            
                
        # p1=self.left
        # p2=self.right
         
         	    
                 
        # if 'L' in op:
             
                
            # p1=simplify(Div(p1,kval1))
            # p1=premath(p1,*args)
            # if unevaluate:
                # p1=Divunevaluate(p1)
 
        # if 'R' in op:
            # p2=simplify(Div(p2,kval2))
            # p2=premath(p2,*args)
            # if unevaluate:
                # p2=Divunevaluate(p2)
                 
         
        # if kname!='':
            # QQ=MyEqEq(p1,p2)
            # return QQ
        # if self.type=='IQ' and signo(kval)==-1:
            # self.e1.ksym=p2
            # self.e2.ksym=p1
        # else:  
            # self.e1.ksym=p1
            # self.e2.ksym=p2
         
        # if kshow:
            # self.s()

    # def Pow(self, *args,simplify=True):
        # kval=args[0]
        # ops,kshow=getops(*args)
        # p1,p2=self.getdata()
         
        # if 'L' in ops:
            # p1=ppow(p1,kval)
            # p1=premath(p1,*args)
            # if simplify:
                # p1=mulexpo(p1)
                # p1=simplifyexpo(p1)
                
        # if 'R' in ops:
            # p2=ppow(p2,kval)
            # p2=premath(p2,*args)
            # if simplify:
                # p2=mulexpo(p2)
                # p2=simplifyexpo(p2)
        # premath      
        # self.e1.ksym=p1
        # self.e2.ksym=p2
        # self.s()
        
    def Add(self,*args):         
            P1,P2=getdatain(*args)
     
            P1=self.e1.ksym+P1
            P2=self.e2.ksym+P2
            if 'simplify' in args:
                P1=simplify(P1)
                P2=simplify(P2)
            self.e1.ksym=P1
            self.e2.ksym=P2
            if 'noshow'not in args:
                self.s()
    def Substrac(self,*args):         
            P1,P2=getdatain(*args)
     
            P1=self.e1.ksym-P1
            P2=self.e2.ksym-P2
            if 'simplify' in args:
                P1=simplify(P1)
                P2=simplify(P2)
            self.e1.ksym=P1
            self.e2.ksym=P2
            if 'noshow'not in args:
                self.s()
    def Subs(self,*args):         
            P1,P2=getdatain(*args)
     
            P1=self.e1.ksym-P1
            P2=self.e2.ksym-P2
            if 'simplify' in args:
                P1=simplify(P1)
                P2=simplify(P2)
            self.e1.ksym=P1
            self.e2.ksym=P2
            if 'noshow'not in args:
                self.s()
    def Mul(self,*args):         
            P1,P2=getdatain(*args)
     
            P1=self.e1.ksym*P1
            P2=self.e2.ksym*P2
            if 'simplify' in args:
                P1=simplify(P1)
                P2=simplify(P2)
            self.e1.ksym=P1
            self.e2.ksym=P2
            if 'noshow'not in args:
                self.s()            

                
    def Div(self,*args):
            P1,P2=getdatain(*args)
     
            P1=self.e1.ksym/P1
            P2=self.e2.ksym/P2
            if 'simplify' in args:
                P1=simplify(P1)
                P2=simplify(P2)
            self.e1.ksym=P1
            self.e2.ksym=P2
            if 'noshow'not in args:
                self.s()
                
    def Pow(self,*args):
            P1,P2=getdatain(*args)
     
            P1=self.e1.ksym**P1
            P2=self.e2.ksym**P2
            if 'simplify' in args:
                P1=simplify(P1)
                P2=simplify(P2)
            self.e1.ksym=P1
            self.e2.ksym=P2
            if 'noshow'not in args:
                self.s()        

    def Rpow(self,*args):
            rr=args[0]
            P1=self.L
            if not 'R' in args:
               
                self.e1.ksym=rsimplify(rpow(P1,rr))
            P2=self.R
            if not 'L' in args:
                self.e2.ksym=rsimplify(rpow(P2,rr))
            
            
            if 'noshow' not in args:
                self.s()
        
    def lexpand(self,kop='RL',kshow=True):
        if 'L' in kop:
            self.e1.ksym=lexpand(self.L)
        if 'R' in kop:
            self.e2.ksym=lexpand(self.R)
        if kshow:
            self.s()
        
     
    
    def lmul2lpow(self, kop='RL',kshow=True):
        p1=self.e1.ksym
        p2=self.e2.ksym
        if 'L' in kop:
            p1=lmul2lpow(p1)
        if 'R' in kop:
            p2=lmul2lpow(p2)
        self.e1.ksym=p1
        self.e2.ksym=p2       
        if kshow:
            self.s()
        
    def lfactor(self, kop='RL',kshow=True):
        if 'L' in kop:
            self.e1.lfactor(kshow=False)
        if 'R' in kop:
            self.e2.lfactor(kshow=False)
        self.s()
    
    def exp(self, kop='RL',kshow=True):
        if 'L' in kop:
            self.e1.ksym=exp(self.L)
        if 'R' in kop:
            self.e2.ksym=exp(self.R)
        self.s()
    def lexponent(self,op='LR'):
        if 'L' in op:
            p1=self.e1.ksym
            expr=p1
            expr=lexponent(expr)
            self.e1.ksym=expr 
        if 'R' in op:
            p2=self.e2.ksym
            expr=p2
            expr=lexponent(expr)
            self.e2.ksym=expr    
 
        self.s()
    def log(self, kop='RL',kshow=True):
        if 'L' in kop:
            kres=self.e1.ksym
            self.e1.ksym=log(kres)
            
        if 'R' in kop:
            kres=self.e2.ksym
            self.e2.ksym=log(kres)
        self.s()
    
    def factorSec(self, ksym,kop='RL',kshow=True):
        if 'L' in kop:
            
            self.e1.factorSec(ksym,kshow=False)
        if 'R' in kop:
            self.e2.factorSec(ksym,kshow=False)
        self.s()
        
    def linfactor(self, ksym,kop='RL',kshow=True):
        if 'L' in kop:
            p1=self.left
            p1=linfactor(p1,ksym)
            self.e1.ksym=p1
        if 'L' in kop:
            p2=self.right
            p2=linfactor(p2,ksym)
            self.e2.ksym=p2        
        if kshow:
            self.s()
        
    def listR(self,op=''):
        kres=fpoly(self.right,'list')
        if  op!='':
            return kres[op]        
        else:
            return kres
    def listL(self,op=''):
        kres=fpoly(self.left,'list')
        if  op!='':
            return kres[op]        
        else:
            return kres        
    
    def list(self,kop='R'):
        kres=fpoly(self.right,'list')
        if kop=='L':
            kres=fpoly(self.left,'list')
        return kres
        
        
    def clearUp(self,vv,kshow=True):
        p1=0
        p2=0
        for i in self.e1.list():
            mm=fpoly(i,'free')
            if vv in mm:
                p1+=i
            else:
                p2-=i
        for i in self.e2.list():
            mm=fpoly(i,'free')
            if vv in mm:
                p1-=i
            else:
                p2+=i
        self.e1=MyEq(p1,'e1',kshow=False)
        self.e2=MyEq(p2,'e2',kshow=False)
        if kshow:
            self.s()
        
    def clearAlone(self,vv):
        self.clearUp(vv=vv,kshow=False)
        mm=self.e1.list()
        done=True
        if len(mm)==2:
            if type(self.e1.ksym)==Mul:
                if vv in self.e1.list():
                    kfac=1
                    for i in self.e1.list():
                        if i!=vv:
                            self.Div(i,kshow=False)
                            self.s()
                            done=False
        if done:
            self.s()
        
    
    def sin2cos(self, angu, korden=2, kope='', kop='RL'):
        if 'L' in kop:
            self.e1.set(kpow(sin(angu), 3), (1 - kpow(cos(angu), 2)) * cos(angu), kshow=False)
            kres = self.e1.ksym
            kres = sin2cos(kres, angu=angu, korden=korden, kope=kope)
            self.e1.update(kres)
        if 'R' in kop:
            self.e2.set(kpow(sin(angu), 3), (1 - kpow(cos(angu), 2)) * cos(angu), kshow=False)
            kres = self.e2.ksym
            kres = sin2cos(kres, angu=angu, korden=korden, kope=kope)
            self.e2.update(kres)
        self.s()

    def cos2sin(self, angu, korden=2, kope='', kop='RL'):
        if 'L' in kop:
            self.e1.set(kpow(cos(angu), 3), (1 - kpow(sin(angu), 2)) * sin(angu), kshow=False)
            kres = self.e1.ksym
            kres = cos2sin(kres, angu=angu, korden=korden, kope=kope)
            self.e1.update(kres)
        if 'R' in kop:
            self.e2.set(kpow(cos(angu), 3), (1 - kpow(sin(angu), 2)) * sin(angu), kshow=False)
            kres = self.e2.ksym
            kres = cos2sin(kres, angu=angu, korden=korden, kope=kope)
            self.e2.update(kres)

            self.s()
                
    #########################################
    #   MyEqEq      evalue 
    ########################################## 

    def evalue(self,val1,val2,kshow=True):
        QQ=self.xcopy()
        QQ.e1.ksym=QQ.left.subs(val1,val2)
        QQ.e2.ksym=QQ.right.subs(val1,val2)
        
        if kshow:
            QQ.s()
            
        else:
            QQ.s(kshow=False)
            return QQ.right
    def eval(self,**kwargs):
        
        QQ=self.xcopy()
        QQ.set(**kwargs)

    
    #########################################
    #   MyEqEq      solve  
    ##########################################
    def solve_compare_coef(self,*args,var=''):
        lvar=[]
        if var=='':
            var=self.var
        if len(args)==1 and type(args[0])==list:
            lvar=args[0]
        if len(args)>1:
            for i in args:
                lvar.append(i)
        p1=self.L
        cp1=coef_list(p1,var)
        p2=self.R
        cp2=coef_list(p2,var)
        vece=[]
        for i,j in zip(cp1,cp2):
            vece.append(i-j)
        solu=solve(vece,lvar)
        display(Math(latex(solu)))
        return ganswer(solu,'value')
    
    
    
    
        
        
        
        
    def solve_coef_list(self,*args,var2=x):
        r'''
        solve variable from  two polinomies with tha same coefficient ,grade orden
        
        example:
        **********
        a*x*x+b*x+c= 3*x*x+2*x+7+y
         
        return:  MyEq
        a=3
        b=2
        c=7+y

        '''
        vecv=[]
        vecq=[]
        for i in args:
            vecv.append(i)
        m1=coef_list( self.left,var2)
        m2=coef_list( self.right,var2)
        qq1=len(m1)
        qq2=len(m2)

        if qq1<qq2:
            qf=qq2-qq1
            vf=mzero(qf)
            m1=vf+m1
        if qq2<qq1:
            qf=qq1-qq2
            vf=mzero(qf)
            m2=vf+m2 
        

        for i,j in zip(m1,m2):
            if (i-j)!=0:
                vecq.append(i-j)
        if len(vecq)>len( vecv):
            vecq=vecq[0:len(vecv)]
        vecq=vecq+vecv
        kres=solver(*vecq)
        eres=[]
        for i,j in zip(vecv,kres):
            ee=MyEq(j,str(i))
            eres.append(ee)
        return eres
    
    
    def alone(self,*args,kshow=True):
        if len(args)==0:
            helplib('alone')
            return
    
        var=args[0]
         
        try:
            expr=self.L-self.R
            ee=MyEq(expr,'ee',var=var,kshow=False)
                
            kk=ee.solve(var,kshow=False)
            self.e1.ksym=var
            self.e2.ksym=kk.ksym
            if kshow:
                self.s()
        except:
            if kshow:
                self.s()
            
        
    
    def alone2(self,*args,kshow=True,**kwargs):
        
        QQ=self.xcopy()
        sres1=str(QQ.e1.ksym)
        sres2=str(QQ.e2.ksym)
        equ=QQ.e1.ksym-QQ.e2.ksym
        ee=MyEq(equ,'ee',kshow=False)
        ksolu=ee.solve(*args,kshow=False)

        self.e1.ksym=args[0]
        self.e2.ksym=ksolu.ksym
        if kshow:
            self.s()
        
        
    def kisolve(self):
        p1=self.L
        p2=self.R
        ps=self.psymb
        var=self.var
        P1=eval(str(p1)+ps+str(p2))
        return solve_univariate_inequality(P1,var)
        
         
    def isolve(self,kname='',kshow=True):
        expr=self.kisolve()
        if kname=='':
            if kshow:
                display(Math(latex(expr)))
            return  expr
        else:
            Qr=MyEq(expr,kname,var=x,ktype='lq')
            return Qr
    def solve(self,*args,kshow=True,**kwargs):
        if len(args)==1:
            e1=MyEq(self.L-self.R,'e1',kshow=false)
            return e1
        done=False
        cvar=args[0]
        svar=str(cvar)
        allpara=[]
        if len(args)>1:
            allpara=args[1::]
            
        kres=self.L-self.R
        if len(kwargs)>1:        
            kres=real_subs(kres,**kwargs)
         
        ee=MyEq(kres,'ee',var=self.var,kshow=False)
        if 'all' in allpara:
            return ee.solve(cvar,'all',kshow=False)
            
        if 'noimg ' in allpara or 'nonimaginary' in allpara:
            return ee.solve(cvar,'noimg',kshow=False)
        
        ksolu=ee.solve(cvar,kshow=kshow)
         
        ee2= MyEq(ksolu,svar,kshow=False)
        if 'update' in allpara:
            self.set(ee2)
        return ee2            
    
    # def solve(self,*args,kshow=True,**kwargs):
    
        # done=False
        # cvar=args[0]
        # svar=str(cvar)
        # allpara=[]
        # if len(args)>1:
            # allpara=args[1::]
            
        # kres=self.L-self.R
        # if len(kwargs)>1:        
            # kres=real_subs(kres,**kwargs)
         
        # ee=MyEq(kres,'ee',var=self.var,kshow=False)
        # if 'all' in allpara:
            # return ee.solve(cvar,'all',kshow=False)
            
        # if 'noimg ' in allpara or 'nonimaginary' in allpara:
            # return ee.solve(cvar,'noimg',kshow=False)
        
        # ksolu=ee.solve(cvar,kshow=kshow)
         
        # ee2= MyEq(ksolu,svar,kshow=False)
        # if 'update' in allpara:
            # self.set(ee2)
        # return ee2
        
    
    
    # def solve(self,*args,kshow=True,**kwargs):
        # if len(args)==0:
            # helplib('solve')
            # return
        # ovar=args[0]
        # args=args[1:len(args)]
        
        # expr=self.L-self.R
        # if ovar=='' or type(ovar)==tuple:
            # if ovar=='':
                # solu=solve(expr)
            # else:
                # solu=solve(expr,ovar)
            
            # if type(solu)==dict:
                # key,value=unpack(solu)
            # elif type(solu)==list:
                # key,value=unpack(solu[0])    
            # vecsol=[]
            # for i,j in zip(key,value):
                # ee=MyEq(j,kname=str(i))
                # vecsol.append(ee)   
            # if len(vecsol)==1:
                # kres= vecsol[0]
            # else:
                # kres=vecsol
              
            
                     
             
                
        # else:    
            # if len(args)>0 and 'all' in args:
                # kres=solve(expr,ovar)
                # if 'nonimaginary' in args:
                     
                    # kres2=[x for x in kres if 'I' not in str(x) ]
                    # return kres2
                # else:    
                    # return kres
            # else:    
                # solu=ksolve(expr,ovar,*args)
                # solu=real_subs(solu,**kwargs)
                # dvalue=False
                # if len(args)>0:
                    # if 'positive' in args:
                        # solu=solu*signo(solu)
                    # if 'value' in args:
                        # dvalue=True
                    # if 'float' in args:
                        # try:
                            # solu=float(solu)
                        # except:
                            # pass
                # if dvalue:
                    # return solu
                # else:
                    # if type(solu)==dict:
                        # key,value=unpack(solu)
                        # if len(value)==1:
                            # solu=value[0]
                          
                    # ee=MyEq(solu,kname=ovar.name,kshow=kshow)
                    # return ee
                    
            
               
        
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
        x=self.vmain
        ee=MyEq(x-self.e2.ksym,kshow=False)
        svar2=ee.solve_if_and(str(svar),eqv=0,kope=kope,korden=korden,kshow=False,**kwargs)
         
        self.e2.set(svar,svar2,kshow=False)
 
        self.s()     
     
    def autosolve(self,kvar,**kwargs):
        QQ=self.xcopy()
        p1=QQ.L
        p1=real_subs(p1,**kwargs)
        p2=QQ.R
        p2=real_subs(p2,**kwargs)
        QQ.e1.ksym=p1
        QQ.e2.ksym=p2
        solu=QQ.solve(kvar,kshow=False)
        ee=MyEq(solu,kname=str(kvar))
        p1=self.L
        p1=subsubs(p1,str(kvar),solu)
        self.e1.ksym=p1
        p2=self.R
        p2=subsubs(p2,str(kvar),solu)
        self.e2.ksym=p2
        
        self.s()
        return ee
    
    def autofind(self,*args,**kwargs):
        QQ=self.xcopy()
        sres1=str(QQ.e1.ksym)
        sres2=str(QQ.e2.ksym)
        for key, value in kwargs.items():
                sres1=sres1.replace(key,str(value))
                sres2=sres2.replace(key,str(value))
        expr1=parse_expr(sres1)-parse_expr(sres2)
        vecres=solve(expr1,args[0])
         
        kres=vecres[0]
        if 'float' in args:
            try:
                kres=float(kres)
            except:
                pass
         
        var=args[0]
        p1=self.L
        p1=supersubs(p1,var,kres)
        self.e1.ksym=p1
        p2=self.R
        p2=supersubs(p2,var,kres)
        self.e2.ksym=p2
        self.s()  
    def solve_if(self,ksym,**kwargs):
        
    
        QQ=self.xcopy()
        p1=QQ.left
        p2=QQ.right
        if len(kwargs)>0:             
            p1=subskwargs(p1,**kwargs) 
            p2=subskwargs(p2,**kwargs)
        kname=str(ksym)
        ee=MyEq(p1-p2,'ee',kshow=False)
        nee=ee.solve(ksym,kshow=False)
        kres=MyEq(nee,kname)
        return kres
        
       
    def solve_if_and_up(self,*args,**kwargs):
        QQ=self.xcopy()
        sres1=str(QQ.e1.ksym)
        sres2=str(QQ.e2.ksym)
        for key, value in kwargs.items():
                sres1=sres1.replace(key,str(value))
                sres2=sres2.replace(key,str(value))
                 
        QQ.e1.ksym=parse_expr(sres1)
        QQ.e2.ksym=parse_expr(sres2)
        kname=''
        qq=len(args) 
        if qq==1:
            vsym=args[0]
        elif qq==2:
            vsym=args[0]
            kname=args[0]
        else:
            QQ.s()
            return
        if kname=='':
            QQ2=QQ.xcopy()
            QQ2.alone(vsym)
            nval=QQ2.solve(vsym)
            self.set(vsym,nval)
            
        else:
            QQ2=QQ.xcopy()
            ee=MyEq(QQ.solve(vsym),kname)
            self.set(vsym,ee.ksym)

            return ee    
    
    def solveset(self,var,**kwargs): # solve kmain in Q and set as vmain
        kres=self.R-self.L 
        kres=subskwargs(kres,**kwargs)
        ee=MyEq(kres,'ee',kshow=False)
        ksolu=ee.solve(var,kshow=False)
        ee2=MyEq(ksolu.ksym,kname=str(var))
        self.upgrade(ee2)
        return ee2
 
        
    def solvemain(self): # solve vmain in Q and set final Eq
        kres=self.solve(self.vmain)
        self.e1.ksym=self.vmain
        self.e2.ksym=kres
         
        self.s()

    def msimplify(self):
        return self.solvemain()
        
    def toMyEq(self, kname):
        kname = MyEq(self.e1.ksym - self.e2.ksym, kname=kname)
        return kname




    def lock(self,op='RL'):
        if 'R' in op:
            ksym=self.R
            self.e2.ksym=UnevaluatedExpr(ksym)
        if 'L' in op:
            ksym=self.L
            self.e1.ksym=UnevaluatedExpr(ksym)
            
    #########################################
    #   MyEqEq      uprade  
    ##########################################

    def replace(self,*args,kshow=True,**kwargs):
         
        if len(kwargs)>0:
               
            p1=self.L
            p1=subskwargs1(p1,**kwargs)
            self.e1.ksym=p1
            p2=self.R 
            p2=subskwargs1(p2,**kwargs)
            self.e2.ksym=p2
            self.eQ=Eq(p1,p2)
             
        if len(args)>0:
            for i in arg:
                p1=self.L
                p2=self.R
                
                if type(i)==MyEq:
                    sname=i.name
                    value=i.ksym
                else:
                    sname=str(i)
                    value=i
                p1=p1.subs(sname,value)
                p2=p2.subs(sname,value)
            self.e1.ksym=p1
            self.e2.ksym=p2
            self.eQ=Eq(p1,p2)
        if kshow:
            self.s()



    def upgrade(self,*args,kshow=True):
        
        self.e1.upgrade(*args,kshow=False)
        self.e2.upgrade(*args,kshow=False)
        self.updateEq()
        self.s()

    def solveandset(self,*args,kshow=True,name='',**kwargs):
        expr1=self.L
        expr2=self.R
        var1=self.var1
        var2=self.var2
        
        QQ=self.xcopy()
        kvar=args[0]
        QQ.alone(kvar,kshow=False)
        p2=QQ.R
        p2=subskwargs(p2,**kwargs)
        kfloat=False
        kupdate=False
        for i in args:
            if i=='float':
                kfloat=True
            if i=='update':
                kupdate=True
        if kfloat:
            try:
                p2=float(p2)
            except:
                pass
        kvar=MyEq(p2,kname=str(kvar),kshow=False)
        kvar.s()
        if kupdate:
            self.upgrade(kvar)
        return kvar
        
        
    def solve_set_if(self,var,korden=0,kope='',kshow=True,**kwargs):
        self.upgrade_if(var=var,korden=korden,kope=kope,kshow=kshow,**kwargs)
          
    def upgrade_if(self,var,korden=0,kope='',kshow=False,**kwargs):
        Q2=self.xcopy()
        Q2.set(kshow=False,**kwargs)
        Q2.simplify(kshow=False)
        if  monodata(Q1.left,'isexp') or monodata(Q2.right,'isexp'):
            Q2.log(kshow=False)
            Q2.lexpand(kshow=False)
        ee2=MyEq(Q2.right-Q2.left,'ee2',kshow=False,kope=kope)
         
        try:
            ee=ee2.ssolve(str(var),kshow=False,kope=kope)
            ee.opemat(kope=kope,kshow=False)
            self.eQ=self.eQ.subs(var,ee)
        except:
            try:
                ee=ee2.solve(var,str(var),korden=korden,kope=kope,kshow=kshow)
            except:    
                ee=ee2.solve(var,str(var),kope=kope)
                return
        if kshow:
            ee.s()
        self.upgrade(ee,kshow=False)
        self.s()
      
    def clear_exp_QQ(self,kshow=True):
    
        p1,p2=clear_exp_QQ(self.left,self.right)
        self.e1.ksym=p1
        self.e2.ksym=p2
        if kshow:
            self.s()
        
    
        
         
        
    def clear(self,ksym):
        self.alone(ksym=ksym) 
        
    def nofunc(self): #  convert all x(t) in x  dx in 1 dx2 in 1 in the Eq
        vv=self.var0
        vf=self.ft
        vdf=self.df
        vd2f=self.d2f
        vold=[vd2f,vdf,vf]
        vnew=[1,1,vv]
        
        p1=str(self.left)
        p2=str(self.right)
         
        for i,j in zip(vold,vnew):
            p1=p1.replace(str(i),str(j))
            p2=p2.replace(str(i),str(j))
            

        self.e1.ksym=parse_expr(p1)
        self.e2.ksym=parse_expr(p2)
        self.s()
        return self.var0
         
    
    def prima_func(self): #  convert all x(t) in x  dx in 1 dx2 in 1 in the Eq
        vv=self.f
        vf=self.ft
        vdf=self.df
        vd2f=self.d2f
        vold=[vd2f,vdf,vf]
        vnew=[self.p2f,self.pf,vv]
        
 
        p1=str(self.left)
        p2=str(self.right)
         
        for i,j in zip(vold,vnew):
            p1=p1.replace(str(i),str(j))
            p2=p2.replace(str(i),str(j))
            

        self.e1.ksym=parse_expr(p1)
        self.e2.ksym=parse_expr(p2)
        self.s()    
    #########################################
    #   MyEqEq     dsolve 
    ##########################################
    def dsolve(self,*args,op=''):
        expr=self.eQ
        kcond=self.ics
        var2=self.var2
        vmain=self.vmain
        
        # if kcond!='':
            # Q2=MyEqEq(dsolve(expr,ics=kcond),var2=t,vmain=x,kshow=False)
        # else:
            # Q2=MyEqEq(dsolve(expr),var2=t,vmain=x,kshow=False)
            
        if kcond!='':
            kres= dsolve(expr,ics=kcond) 
        else:
            kres=dsolve(expr)

        p1=kres.lhs
        p2=kres.rhs
        kname=str(self.vmain)
        if len(args)==1 and type(args[0])==str:
            kname=args[0]
            nval=symbols(kname)
            nval=MyEq(p2,kname,var2=self.var2)
            return nval
        else:
            nval=symbols(kname)
            nval=MyEq(p2,kname,var2=self.var2)
            return nval
        

    def nodiff(self):
        df=symbol2diff(self.vmain,self.var2)
        d2f=symbol2diff2(self.vmain,self.var2)
        p1=self.left
        p1=p1.subs(d2f,1)
        p1=p1.subs(df,1)
        p2=self.right
        p2=p2.subs(d2f,1)
        p2=p2.subs(df,1)
        self.e1.ksym=p1
        self.e2.ksym=p2
        self.s()
    def derivarespect(self,var1,D1,var2='',D2='',kshow=True):
        p1=derivarespect(self.L,var1,D1)
        if var2=='':
            var2=var1
            D2=D1
        p2=derivarespect(self.R,var2,D2)
        self.e1.ksym=p1
        self.e2.ksym=p2
        if kshow:
            self.s()    
    def chainrule(self,p1,p2,p3,kshow=True):
        PL=self.L
        PL=PL.subs(p1/p2,p3)
        PR=self.R
        PR=PR.subs(p1/p2,p3)
        self.e1.ksym=PL
        self.e2.ksym=PR
        if kshow:
            self.s()
            
    def diffQ(self,*args,niceview=True,update=True):
        ''' 
        input: diff((),())
        return  original eQ
        
        input: diff((x,t),())
        return diff(Q.L, in t and  x=func) =  Q.R )
        
        input: diff((),(t,y))
        return Q.L = diff(Q.R, in t and  y=func)
        
        input: diff((t,x),(t,y))
        return diff(Q.L, in t and  x=func) = diff(Q.R, in t and  y=func)
        
        input: diff((t,x,y),(t,x,y))
        return diff(Q.L, in t and  x,y =func in t) = diff(Q.R, in t and  x,y=func in t)
        
        args:
            update=True, self.update, 
            update=False  return other MQ
            
            niceview=True
            return example Diferential(x(t),t) = dx/dt but  value is original
        '''    
        exp1=self.e1.ksym
        exp2=self.e2.ksym
         
        done1=False
        done2=False
        
        if args[0]!=():
            datax=args[0]

            vart1=datax[0]
            varx1=datax[1]
            vary1=''
               
            if len(datax)==3:
                vary1=datax[2]
            exp11=diffuntion(exp1,vart1,varx1,vary1)
            dif1=diffuntion(exp1,vart1,varx1,vary1)
            exp1=exp11
            done1=True
        if args[1]!=():
             
            datay=args[1]

            vart2=datay[0]
            varx2=datay[1]
            vary2=''
              
            if len(datay)==3:
                vary2=datay[2]
                
            exp22=diffuntion(exp2,vart2,varx2,vary2)
            dif2=diffuntion(exp2,vart2,varx2,vary2)
            exp2=exp22
            done2=True
        if update:
            self.e1.ksym=exp1
            self.e2.ksym=exp2
            rfin=self
        else:
            rfin=MQ(exp1,exp2,kshow=False)
        if niceview:
            nexp1=exp1
            if done1:
                nexp1=viewnicediff(exp1,vart1,varx1,vary1)
            nexp2=exp2
            if done2:
                nexp2=viewnicediff(exp2,vart2,varx2,vary2)    
            if update==False:
                MQ(nexp1,nexp2)
                return MQ(dif1,dif2,kshow=False)
            else:
                MQ(nexp1,nexp2)

        else:    
            rfin.s() 
    
    
    def normalizediff(self,sf,var):
         
        
        p1=self.L
        P1=normalizediff(p1,sf=sf,var=var)
         
        p2=self.R
        P2=normalizediff(p2,sf=sf,var=var)
        sf=symbols(sf)
        self.e1.ksym=P1
        self.e2.ksym=P2
        self.s()
          

        
    def fulldiff(self):
         
        e1= MyFullDiff(self.e1)
        self.e1.ksym=unisymbols(e1.ksym)
        self.e1.diffvar=get_diffvar(e1)
        e2= MyFullDiff(self.e2)
        self.e2.ksym=unisymbols(e2.ksym)
        self.e2.diffvar=get_diffvar(e2)
        self.s()
        
    
    
    def diff2(self,*args): 
        var=self.var
        p1=self.L
        P1=functiondiffk(p1,var,*args)
        p2=self.R
        P2=functiondiffk(p2,var,*args)
        self.e1.ksym=P1
        self.e2.ksym=P2
        self.s() 
        
    def diff(self,*args):
        var=self.var
        p1=self.L
        ff1=prime2sympy(p1,var,*args)
        dff1=diff(ff1,var)
        p2=self.R
        ff2=prime2sympy(p2,var,*args)
        dff2=diff(ff2,var)
        self.e1.ksym=sympyEq2prime(dff1,var,*args)
        self.e2.ksym=sympyEq2prime(dff2,var,*args)
        self.s()        
        
    def hightdiff(self,vec1='',vec2='',*args):
        var=self.var
        if vec1!='':
            ksym1=self.L
            if type(vec1)!=list:
                vec1=[vec1]
            kres=hightdiff(ksym1,var,*vec1)
            if 'nview' in args:
                kres=sympyEq2prime(kres,var,*vec1)
            self.e1.ksym=kres
        if vec2!='':
            ksym2=self.R
            if type(vec2)!=list:
                vec2=[vec2]
            kres=hightdiff(ksym2,var,*vec2)
            if 'nview' in args:
                kres=sympyEq2prime(kres,var,*vec2)
            self.e2.ksym=kres 
        self.s()        
    def maximize(self,kope=''):
        var2=self.var2
        kres=diff(self.right,self.var2)
        kres=simplify(kres)
        kres=reduFac(kres)
        try:
            e1=MyEq(kres,kshow=False)
            e1.reduFac(kshow=False)
            kres=e1.ksym
        except:
            pass
        kres=solve(kres,self.var2)
        kres2=[]
        for i in kres:
            sres=str(i)
            if not 'I' in sres:
                kres2.append(i)
        kres=kres2
        if len(kres)==1:
            for i in kres:
                MyEq(i,str(var2))
                self.evalue(var2,i)
        else:
            cc=0
            cc2=0
            qq=len(kres)
            fres=self.evalue(var2,kres[0],kshow=False)
            for i in range(qq):
                ik=kres[i]
                val=self.evalue(var2,ik,kshow=False)
                val2=val 
                try:
                    val=float(val)
                except:
                    pass
                if val>fres:
                    fres=val
                    cc2=cc
                    val3=val
                cc+=1
                
             
            MyEq(kres[cc2],str(var2),kope=kope)
            self.evalue(var2,opemat(kres[cc2]),kope)    
        

    def apply_chainRuler(self):
        kres1=self.left
        try: 
            kres1=simplify_chain_rule(kres1)
            self.e1.ksym=kres1
        except:
            pass
            
        kres2=self.right    
        try: 
            kres2=simplify_chain_rule(kres2)
            self.e2.ksym=kres2
        except:
            pass    
        self.s()
        
    
    def coef_list(self,var=x):
        p1=self.L
        p2=self.R
        P=p1-p2
        vec=coef_list(P,var)
        return vec


        
    #########################################
    #   MyEqEq     Integral
    ##########################################
    def mode_integral(self,kname='',kside='LR'):
        if 'L' in kside:
            p1=self.left
            if type(p1)==Derivative:
                p1= integrate(self.left) 
            else:
                p1=Integral(p1,self.var0)
            self.e1.ksym=p1
             
        if 'R' in kside :
            p2=self.right
            if type(p2)==Derivative:
                p2= integrate(self.left) 
            else:
                p2=Integral(p2,self.var2)
            self.e2.ksym=p2
            
        if kname!='':
            kname2=kname
            kname2=MyEq(self.right,kname,var2=self.var2)
            return kname2
        else:    
            self.s()    
     
    def applyIntegral(self ):

        kname=''
        p1=self.left
        if type(p1)==Derivative:
            p1=Integral(integrate(Q2.left))
        else:
            p1=Integral(p1,self.vmain)
        p2=self.right    
        if type(p2)==Derivative:
            p2=Integral(integrate(Q1.right))
        else:
            p2=Integral(p2,self.var2)

        self.e1.ksym=p1
        self.e2.ksym=p2
        self.s()
        
        
    def separe_dif(self):
        var1=self.vmain
        name1='d'+str(var1)
        
        
        var2=self.var2
        name2='d'+str(var2)
        
        d1=symbols(name1)
        d2=symbols(name2)

        self.e1.ksym=d1/self.right
        self.e1.var2=var1
        self.e2.ksym=d2 
        self.e2.var2=var2
        self.type='eS'
        self.s()
    
    def eQinteger(self,var1,var2):
        d1='d'+str(var1)
        d2='d'+str(var2)
         
        d1=symbols(d1)
        d2=symbols(d2)
        p1=self.L
        try:
            p1=p1.subs(d1,1)
        except:
            pass
        p2=self.R
        try:
            p2=p2.subs(d2,1)
        except:
            pass
        if Is_Number(p1):
            Q1=MyEq(p1*var1,'Q1',var=var1,kshow=False)
        else:
            Q1=MyIntg(p1,'Q1',var=var1,kshow=False)
            
        if Is_Number(p2):
            Q2=MyEq(p2*var2+C1 ,'Q2',var=var2,kshow=False)
        else:
            Q2=MyIntg(p2,'Q2',var=var2,kshow=False)

        return MQ(Q1,Q2)    
    
    # def integrate(self,var1=y,var2=x):
        # P1=self.L
        # P2=self.R
        # v1,v2=varDiff(str(var1),str(var2))
        # P1=P1.subs(v1,1)
        # P2=P2.subs(v2,1)
        # sv1='d'+str(v1)
        # sv2='d'+str(v2)
        # P1=P1.subs(sv1,1)
        # P2=P2.subs(sv2,1)
        # ee1=MyIntg(P1,'ee1',var=var1,kshow=False)
        # ee2=MyIntg(P2,'ee2',var=var2,kshow=False)
        # self.e1.ksym=ee1.ksym
        # self.e1.type='I'
        # self.e1.var=var1
        # self.e2.ksym=ee2.ksym
        # self.e2.var=var2
        # self.e2.type='I'
        # self.s()
        
    
    def factorinte(self,*args):
        kfac=args[0]
        op='RL'
        if 'L'in args:
            op='L'
        if 'R' in args:
            op='R'
        p1=self.L 
        if 'L' in op:
            kres=factorinte(p1,kfac)
            self.e1.ksym=kres
        p2=self.R
        if 'R' in op:
            kres=factorinte(p2,kfac)
            self.e2.ksym=kres
        self.s()
            
            
    def integral(self,var1='',var2=''):
        sd1=str(var1)+"'"
        sd2=str(var2)+"'"
        p1=self.L
        if var1=='':
            iP1=integrate(p1)
        else:
            p1=p1.subs(sd1,1)
            iP1=creteInt(p1,var1)
        p2=self.R
        if var2=='':
            iP2=integrate(p2)
        else:
            p2=p2.subs(sd2,1)
            iP2=creteInt(p2,var2)

        
         
        try:
            ip1=ip1.subs('C1',0)
        except:
            pass
        

        self.e1.ksym=iP1
        self.e2.ksym=iP2
        
         
        self.e1.type='I'
        # self.e1.var=var1
        # self.e2.ksym=ee2.ksym
        # self.e2.var=var2
        self.e2.type='I'
        self.s()   
            
            
            
    
    def IntegralEq(self,v1='',v2=''):
         
        x=self.vmain
        t=self.var2
        dx=symbol2diff(x,t)
        d2x=symbol2diff2(x,t)
        p1=self.left
        p1=p1.subs(dx,1)
        p2=self.right
        p2=p2.subs(dx,1)
        v11=self.vmain
        
        if v1!='':
            v11=v1
            
        v22=self.var2
        if v2!='':
            v22=v2
        v11=symbols(str(v11))
        v22=symbols(str(v22))
        p1=func2symbol(p1,self.vmain,self.var2)
        p2=func2symbol(p2,self.vmain,self.var2)
        if Is_Poly(p1):
            mm=0
            for i in fpoly(p1,'list'):
               kk=MyEq(i,'',var2=v11,ktype='I',kshow=False)
               mm=mm+kk.ksym
            p1=mm
        else:    
            kk=MyEq(p1,'p1',var2=v11,ktype='I',kshow=False)
            p1=kk.ksym
        if Is_Poly(p2):
            mm=0
            for i in fpoly(p2,'list'):
               kk=MyEq(i,'',var2=v11,ktype='I',kshow=False)
               mm=mm+kk.ksym
            p2=mm
        else:    
            kk=MyEq(p2,'p2',var2=v22,ktype='I',kshow=False)
            p2=kk.ksym
            
         
        self.e1.ksym=p1 
        self.e2.ksym=p2 
         
        self.s()
    
    # def Integral(self,*args):

        # var21=symbols(args[0])
        # var22=symbols(args[1])
        # # clean
        # nd1= 'd'+args[0]
        # nd2= 'd'+args[1]
        # e1s=str(self.left)
        # e1s=e1s.replace(nd1,'1')
        # e1s=e1s.replace(nd2,'1')
         
        # self.e1.ksym=parse_expr(e1s)
        
        # e2s=str(self.right)
        # e2s=e2s.replace(nd1,'1')
        # e2s=e2s.replace(nd2,'1')
         
        # self.e2.ksym=parse_expr(e2s)
        # self.s(kshow=False)
        
        # p1=self.left
        # p2=self.right
        # x1=''
        # x2=0
        # if len(args)==4:
            # x1=args[2]
            # x2=args[3]
            # ne1=MyEq(p1,'e1',var2=var21,x1=x1,x2=x2,ktype='I',kshow=False)
            # ne2=MyEq(p2,'e2',var2=var22,x1=x1,x2=x2,ktype='I',kshow=False) 
        # else:    
        
            # ne1=MyEq(p1,'e1',var2=var21,ktype='Ix2x2',kshow=False)
            # ne2=MyEq(p2,'e2',var2=var22,ktype='I',kshow=False)
             
        
        # self.e1=ne1
        # self.e2=ne2
        # self.s()
        
        
    def dointegral(self,*args):
        data1=''
        data2=''
        var1=x
        var2=t
        x1=0
        x2=x
        y1=0
        y2=t
        done1=True
        done2=True
        if  len(args)==2 and type(args[0])==Symbol:
             
            var1=args[0]
            x2=var1 
            var2=args[1]
            y2=var2  
        for i in args:
            if type(i)==Symbol:
                if done1:
                    var1=i
                    done1=False
                else:
                    var2=i
                    
            if type(i)==tuple:
                if done2:
                    var1,x1,x2=i  
                    done2=False
                else:
                    var2,y1,y2=i
        p1=sympify(self.L)
        p2=parse_expr(str(self.e2.ksym))
        self.e1.ksym=Integral(p1,(var1,x1,x2))
        self.e2.ksym=Integral(p2,(var2,y1,y2))
           
        self.s()
        
        

        
        
    def doit(self,op='RL'):
        p1=self.L
        if 'L' in op:
            p1=p1.doit()
        p2=self.R
        if 'R' in op:
            p2=p2.doit()
            try:
                if self.e2.type=='I':
                    newC=self.vecCC[self.CC]
                    self.CC=self.CC+1
                    p2=p2+newC
            except:
                pass
        self.updatee(p1,p2)
        self.s()        
         
    def doitI(self,*args):
        exp1=self.e1.ksym
        var1=self.e1.var
        kres1=exp1.doit()
        exp2=self.e2.ksym
        var2=self.e2.var
        kres2=exp2.doit()
        self.e1.ksym=kres1
        self.e1.type='P'
        self.e2.ksym=kres2+C1
        self.e2.type='P'
              
        self.s()

    #########################################        
    #  algebrate transformation
    #########################################

    def transformada(self,expr0,ssym,kope='LR'):
        if 'L' in kope:
            p1=self.L
            try:
                p1=transformada(p1,expr0=expr0,ssym=ssym,kope=kope)
            except:
                pass
            self.e1.ksym=p1
        if 'R' in kope:
            p2=self.R
            try:
                p2=transformada(p2,expr0=expr0,ssym=ssym,kope=kope)
            except:
                pass
            self.e2.ksym=p2    

        self.s()
        
    def format(self,qn,kshow):
        p1=self.L
        p2=self.R
        try:
            kres1=format(p1,qn)
            self.e1.ksym=kres1    
        except:
            pass
        try:
            kres3=format(p2,qn)
            self.e2.ksym=kres2    
        except:
            pass
        
        if  kshow:
            self.s()
    
    def subsubs(self,*args,force=False,op='LR'):
        if len(args)==0:
            helplib('subsubs')
            return
        expr1=args[0]
        expr2=args[1]
        
        if 'L' in op:
            expr=self.e1.ksym
            expr=subsubs(expr,expr1,expr2,force=force)
            self.e1.ksym=expr
             
        if 'R' in op:
            expr=self.e2.ksym
            expr=subsubs(expr,expr1,expr2,force=force)
            self.e2.ksym=expr
        self.s()    

    def subsnumber(self,val1,val2):
     
         
        p1=self.left
        p1=subsnumber(p1,val1,val2)
        self.e1.ksym=p1
    
        p2=self.right
        p2=subsnumber(p2,val1,val2)
        self.e2.ksym=p2
        self.s()
        
    def subs(self,*args,kshow=True):
        val1,val2=args[0],args[1]
         
        p1=self.left
        try:
            p1=p1.subs(val1,val2)
            self.e1.ksym=p1
        except:
            pass
            
            
        p2=self.right
        try:
            p2=p2.subs(val1,val2)
            self.e2.ksym=p2
        except:
            pass
        if kshow: 
            self.s()
    def sortdegree(self,var='LR',kope=''):
        '''
        input (x**(a*b))   ---->   return(x**a)**b
        input (x**(a*b),b)   ---->   return(x**b)**a
        '''
        kside='LR'
        if type(var)==str:
            kside=var
         
            
        if  type(var)!=str:
            var2=var
        else:
            var2=x
        if 'L' in kside:
            p1=self.left
            p1=sortdegree(p1,var2=var2)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=sortdegree(p2,var2=var2)
            self.e2.ksym=p2
        self.s()   
    
    def powexpand(self,args):
        '''
        input (x**(a*b))   ---->   return(x**a)**b
        input (x**(a*b),b)   ---->   return(x**b)**a
        '''
        op=''
        kside='LR'
        for i in args:
            if type(i)==str:
                kside=i
            if type(i)==Symbol:
                op=i
        if 'L' in kside:
            p1=self.left
            p1=powexpand(p1,op=op)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=powexpand(p2,op=op)
            self.e2.ksym=p2
        self.s()    
    
    def opematexp(self,kside='LR',kope=''):
        '''
        apply opemat only in exponent monomie
        '''
        if 'L' in kside:
            p1=self.left
            p1=opematexp(p1,kope=kope)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=opematexp(p2,kope=kope)
            self.e2.ksym=p2
        self.s()
        
    def mulexpo(self,kside='LR'):
        '''
        input (x**a)**b
        return(x**(a*b))
        '''
        if 'L' in kside:
            p1=self.left
            p1=mulexpo(p1)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=mulexpo(p2)
            self.e2.ksym=p2
        self.s()
        
     
    def packexp(self,kside='LR'):
        '''
        input (x**a)**b
        return(x**(a*b))
        '''
        if 'L' in kside:
            p1=self.left
            p1=packexp(p1)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=packexp(p2)
            self.e2.ksym=p2
        self.s()
    
    
    
    def base2frac(self,kside='LR'):
         
        if 'L' in kside:
            p1=self.left
            p1=base2frac(p1)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=base2frac(p2)
            self.e2.ksym=p2
        self.s()


    
    def expandexp(self,kside='LR'):
        '''
        input (x**a)**b
        return(x**(a*b))
        '''
        if 'L' in kside:
            p1=self.left
            p1=expandexp(p1)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=expandexp(p2)
            self.e2.ksym=p2
        self.s()
        
    def factorexpo(self,*args):
        '''
        input (x**a)**b
        return(x**(a*b))
        '''
         
        if 'right' not in args: 
             
            p1=self.left
            sp1=str(factorexpo(p1))
             
            self.e1.ksym=parse_expr(sp1,evaluate=False)
             
        if 'left' not in args: 
            p2=self.right
            sp2=str(factorexpo(p2))
             
            self.e2.ksym=parse_expr(sp2,evaluate=False)
            
        self.s()
    def disjoinbase(self,kside='LR'):
        P1=self.e1.ksym 
        P2=self.e2.ksym
        if 'L' in kside:
            P1=separebase(P1)
        if 'R' in kside:
            P2=separebase(P2)
        self.e1.ksym=P1
        self.e2.ksym=P2
        self.s() 
        
    def joinbase(self,kside='LR'):
        '''
        input (x**a)**b
        return(x**(a*b))
        '''
        if 'L' in kside:
            p1=self.left
            p1=joinbase(p1)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=joinbase(p2)
            self.e2.ksym=p2
        self.s()
    def disjoinexpo(self,kside='LR'):
        if 'L' in kside:
            p1=self.left
            p1=disjoinexpo(p1)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=disjoinexpo(p2)
            self.e2.ksym=p2
        self.s()
    def joinexpo(self,kside='LR'):
        if 'L' in kside:
            p1=self.left
            p1=joinexpo(p1)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=joinexpo(p2)
            self.e2.ksym=p2
        self.s()    
    def separebase(self,kside='LR'):
        '''
        input (x**a)**b
        return(x**(a*b))
        '''
        if 'L' in kside:
            p1=self.left
            p1=separebase(p1)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=separebase(p2)
            self.e2.ksym=p2
        self.s()
        
    def sum2mulexpo(self,kside='LR'):
        '''
        input (x**a)**b
        return(x**(a*b))
        '''
        if 'L' in kside:
            p1=self.left
            p1=sum2mulexpo(p1)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=sum2mulexpo(p2)
            self.e2.ksym=p2
        self.s()    
    def factorize(self,kval,kside='LR'):
        if 'L' in kside:
            p1=self.left
            p1=factorize(p1,kval)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=factorize(p2,kval)
            self.e2.ksym=p2
        self.s()
        
    def factoriza(self,kval,kside='LR'):
        if 'L' in kside:
            p1=self.left
            p1=factoriza(p1,kval)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=factoriza(p2,kval)
            self.e2.ksym=p2
        self.s()    
        
    def killsqrtpow(self,kside='LR'):
        if 'L' in kside:
            p1=self.left
            p1=killsqrtpow(p1)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=killsqrtpow(p2)
            self.e2.ksym=p2
        self.s() 
    
    
    def killroot(self):
        p1=self.left
        p11=simplify_root_exp(p1)
        p2=self.right
        p22=simplify_root_exp(p2)
        self.e1.ksym=p11
        self.e2.ksym=p22   

        self.s()
    def killRootPow(self):
        p1=self.left
        p11=kill_RootPow(p1)
        if str(p1)!=str(p11):
            self.e1.ksym=p11
        p2=self.right
        p22=kill_RootPow(p2)
        if str(p2)!=str(p22):
            self.e2.ksym=p22
        self.s()    

        
    def par_frac(self,var=''):
        self.e1.par_frac(var=var,kshow=False)
        self.e2.par_frac(var=var,kshow=False)
        self.s()
        
    def apply(self,sexpr):
        p1=self.e1.ksym
        kres1=apply(p1,sexpr)
        self.e1.ksym=kres1
        p2=self.e2.ksym
        kres2=apply(p2,sexpr)
        self.e2.ksym=kres2
        self.s()    
    def simple_linear(self):
        if self.e1.Is_Mono() and self.e2.Is_Mono():
            p1=self.e1.get_nume()
            p2=self.e1.get_deno()
            p3=self.e2.get_nume()
            p4=self.e2.get_deno()
            self.e1.ksym=p1/p3
            self.e2.ksym=p2*p4
            mm=fpoly(self.left,'list')
            if -1 in mm:
                self.Mul(-1)
            else:
                self.s()    
    
        
    def crossmul(self,kshow=True):
        p1=factor(self.L)
        n1=numer(p1)
        d1=denom(p1)
        p2=factor(self.R)
        n2=numer(p2)
        d2=denom(p2)
        self.e1.ksym=n1*d2
        self.e2.ksym=n2*d1
        if kshow:
            self.s()
            
        
       
       
    def crossMul(self,*args):
        P1=self.L
        P2=self.R
        if Is_Div(P1): 
            pa1,pa2=fraction(P1)
        else:
            pa1,pa2=P1,1
        if Is_Div(P2): 
            pb1,pb2=fraction(P2)
        else:
            pb1,pb2=P2,1    
            
         
        Px=unisymbols(pa1*pb2)
        Py=unisymbols(pb1*pa2)
        
        
        
        Px=premath(Px,*args)
        Py=premath(Py,*args)

        
        self.e1.ksym=Px
        self.e2.ksym=Py
        
        self.s()



tau,psi,omega=symbols('tau psi omega')
varT=[tau,psi,omega] 

       
def eQtools(ee,op=0):
    nvar=varT[op]
    return MQ(nvar,ee)
def real_subs6(QQ,**kwargs):
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
        if  type(value)==MyEq:
            vvalue.append(value.ksym)
        else:    
            vvalue.append(value)
        
        kres=kres.subs(parse_expr(key),value)
    for i,j in zip(mvar,vvar):
        kres=kres.subs(i,j)
    return (kres)  

def real_subs(expr,**kwargs):
    ''' 
    QQ= symbols function
    ** kwargs c7=6,g4=z..etc..
    RETURN real substitucion when variable have underscore name like 'c_7' 'g_4'    
    '''
    if len(kwargs)>0:
        key,value=unpack(kwargs)
        kres=expr 
        for i,j in zip(key,value):
            jj=j
            if type(j)==MyEq:
                jj=j.ksym
            kres=kres.subs(i,j)
            if len(i)>1:
                newi=i[0]+'_'+i[1::]
                try:
                    kres=kres.subs(newi,j)
                except:
                    pass
            
        return kres
    else:
        return expr

def symbol_subs(expr,**kwargs):
    ''' 
    QQ= symbols function
    ** kwargs c7=6,g4=z..etc..
    RETURN real substitucion when variable have underscore name like 'c_7' 'g_4'    
    '''
    if len(kwargs)>0:
        key,value=unpack(kwargs)
        kres=expr
        sres=str(kres)
        for i,j in zip(key,value):
            jj=j
            sjj=str(j)
            if type(j)==MyEq:
                jj=j.ksym
                sjj=str(j.ksym)
            kres=kres.subs(i,j)
            sres=sres.replace(str(i),str(j))
             
            if len(i)>1:
                newi=i[0]+'_'+i[1::]
                try:
                    kres=kres.subs(newi,j)
                    sres=sres.replace(str(newi),str(j))
                except:
                    pass
            
        return parse_expr(sres,evaluate=False)
    else:
        return expr        

def varDiff(*args):
    mm=[]
    for i in args:
         
        k=symbols(i+"'")
        mm.append(k)
    return(mm) 
def varDiff2(*args):
    mm=[]
    for i in args:
         
        k=symbols(i+"''")
        mm.append(k)
    return(mm)   
    
def nice_iq_symbol(ksymbol):
        if ksymbol=='>=' or ksymbol=='=>':
            return '≥'
        elif  ksymbol=='<=' or ksymbol=='=<':
            return '≤'
        else:
            return ksymbol


def apply(expr,sexpr):
    sres=str(expr)
    sres=sexpr+'('+sres+')'
    return parse_expr(sres)
    
def getops(*args):
    ops='LR'
    kshow=True
    if len(args)>0:
        if 'L' in args and not 'R' in args:
            ops='L'
        elif 'R' in args and not 'L' in args:
            ops='R'
        else:
            ops='LR'
        if 'noshow' in args:
            kshow=False
    return ops,kshow


def datacoeff(data,var=x):
    if type(data)==list:
        return data
    elif type(data)==MyEq:
        return coef_list(data.ksym,var)
    else:
        return coef_list(data,var) 
        
def solvecomparecoef(*args):
    val1=args[0]
    if type(val1)==MyEqEq:
        var=args[1]
        vsol=args[2]
        p1=coef_list(val1.L,var)
        p2=coef_list(val1.R,var)
        if type(vsol)!=list:
            qq=len(args)
            Vsol=[]
            for i in args[2::]:
                Vsol.append(i)
        else:
            Vsol=vsol
    else:
        var=args[2]
        p1=datacoeff(args[0],var)
        p2=datacoeff(args[1],var)
        vsol=args[3]
        if type(vsol)!=list:
            qq=len(args)
            Vsol=[]
            for i in args[3::]:
                Vsol.append(i)
        else:
            Vsol=vsol
            
    

    kres=  solve_compare_coef(p1,p2,Vsol)
    for i,j in zip(Vsol,kres):
        display(Math((latex(i)+'= '+latex(j))))
        #display(Math(latex(i)+'= '+latex(j)))
    if len(kres)==1:
        kres=kres[0]

    return kres
    
def edit_primitive(eQ,F):
    expr1=eQ.get_primitive()
    Q=MyEqEq(F,expr1)
    return Q
    
def reparte_param(vecop,*args):
    vecexpr=[]
    kname=''
    op='RL'
    for i in args:
        if type(i)==str and i not in vecop:
            kname=i
        if type(i)!=str:
            vecexpr.append(i)
    if 'R' in args:
        op='R'
    if 'L' in args:
        op='L'
    if len(vecexpr)==1:
        kres=vecexpr[0]
    if type(kres)==MyEqEq:
        expr1,expr2=kres.L,kres.R 
    elif type(kres)==MyEq:
        expr1,expr2=kres.ksym,kres.ksym
    else:
        expr1,expr2=kres,kres
        
    return expr1,expr2,kname,op
    
        
        
def fes(expr,*args):
   
    if 'factor' in args:
        expr=factor(expr)
    if 'simplify' in args:
        expr=simplify(expr)
    if 'expand' in args:    
        expr=expand(expr)
    return expr

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


def argsmap (expr,kname='A',deep=2,data=False,side=''):
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
    if side!='':
        mapL=[side+i for i in mappos]
        mappos=mapL
        
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
class MySys(MyEqEq):
    def __init__(self, *args,kshow=True):
        S=[]
        for i in args:
            S.append(i)
            i.insys=self
        self.S=S    
        cc=1
        for i in self.S:
            p1=i.e1.ksym
            p2=i.e2.ksym
            ps=i.symb
            ps2='Q'+str(cc)+' :'
            display(Math(ps2+latex(p1)+' '+ps+' '+latex(p2)))
            cc=cc+1
    def __call__(self,*args, **kwargs):
        pass
    def __repr__(self):
        return ''
    def s2(self):
        for i in self.S:
            i.s()
    
    def s(self):
        cc=1
        for i in self.S:
            p1=i.e1.ksym
            p2=i.e2.ksym
            ps=i.symb
            ps2='Q'+str(cc)+' :'
            display(Math(ps2+latex(p1)+' '+ps+' '+latex(p2)))
            cc=cc+1    
    def factor(self,kshow=True):
        S2=[]
        for i in range(len(self.S)):
            kres=self.S[i]
            kres.factor(kshow=False)
            S2.append(kres)   
        self.S=S2
        if kshow:
            self.s()
    def subs(self,*args,kshow=True):
        S2=[]
        vecs=self.S 
        for i in range(len(self.S)):
            kres=vecs[i]
            kres.subs(*args,kshow=False)
            if kres.L!=kres.R:
                S2.append(kres)   
        self.S=S2
        if kshow:
            self.s()
    def simplesolve(self,*args):
        S2=self.S
        if len(args)>0:
            for i in args:
                S2.append(i)
        return simplesolve(*S2) 

class MyCelule(MyEqEq):

    def __init__(self, *args,kshow=True):
        op=['MQ','EQ']
        vece=[]
        vecn=[]
        kname=''
        self.sp1=args[0]
        self.sp2=args[1]
        if len(args)>2:
            for i in range(2,len(args)):
                kres=args[i]
                vece.append(kres)
                vecn.append(kres.name)
        self.vece=vece
        self.vecn=vecn
        #self.Q=MQ(eval(self.sp1),eval(self.sp2))
        for i in self.vece:
            i.link=self
        self.e1=MyEq(eval(self.sp1),'e1',kshow=False)
        self.e2=MyEq(eval(self.sp2),'e2',kshow=False)
        p1=self.e1.ksym
        p2=self.e2.ksym
        ps='='
        if kshow:
            display(Eq(self.e1.ksym,self.e2.ksym))
    def s(self):
        self.e1=MyEq(eval(self.sp1),'e1',var=var,kshow=False)
        self.e2=MyEq(eval(self.sp2),'e2',var=var,kshow=False)
 
        ps='='
        display(Math(latex(p1)+' '+ps+' '+latex(p2)))
                
def getdatain(*args):
    kdata=[ i for i in args if type(i)!=str]
    expr=kdata[0]     
    if type(expr)==MyEqEq:
        return expr.L,expr.R
    elif type(expr)==MyEq:
        return expr.ksym,expr.ksym
    else:
        return expr,expr

def doinRL(P1,P2,thisL,thisR,*args):
    if 'R' not in args: 
        P1=thisL
        P1=premath(P1,*args)
     
    if 'L' not in args: 
        P2=thisR
        P2=premath(P2,*args)
     
    return P1,P2 

def get_ksym(obj):
    if type(obj)==MyEqEq:
        return(obj.L-obj.R)
    elif type(obj)==MyEq:
        return obj.ksym
    elif type(obj)==str:
        return parse_expr(obj)
    else:
        return obj
def obj2MyEq(obj,var=x):
    ksym=get_ksym(obj)
    ee=MyEq(ksym,'ee',var=var,kshow=False)
    return ee 

def datadiff(var,*args):
    '''
       var=main var t
       args=[x,y]= x(t),y(t)
       Vs=['x','y']
       Vn=[x,y]
       Vf=[x(t),y(t)]
       Vd=[dxt,dyt]
       Vd2=[dx2t,dy2t]
       Pd=[x',y']
       Pd2=[x'',y'']
    '''
    
    Vs,Vn,Vf,Vd,Vd2,Pd,Pd2 =[],[],[],[],[],[],[]
    P1,P2=[],[]
    for i in args:
        Vs.append(str(i))
        Vn.append(i)
    for i in Vs:
        Vf.append(Function(i)(var))
    for i in Vf:
        Vd.append(i.diff(var))
        Vd2.append(i.diff(var,var))
    for i in Vn:
        Pd.append(symboldiff(i))
        Pd2.append(symboldiff2(i))
        
    return Vs,Vn,Vf,Vd,Vd2,Pd,Pd2 

def hightfunc(expr,var,*args):
    '''
    expr=pi*h*r*r*t/3
    var=t
    args=h,r 
    return pi*t*h(t)*r(t)**2/3 
    '''
    Vs,Vn,Vf,Vd,Vd2,Pd,Pd2=datadiff(var,*args)
    for i,j in zip(Vn,Vf):
        expr=subsubs(expr,i,j)
    return expr 

def hightdiff(expr,var,*args):
    '''
    expr=pi*h*r*r*t/3
    var=t
    args=h,r 
    return 2*pi*t*h(t)*r(t)*Derivative(r(t), t)/3 + pi*t*r(t)**2*Derivative(h(t), t)/3 
           + pi*h(t)*r(t)**2/3
    '''

    Vs,Vn,Vf,Vd,Vd2,Pd,Pd2=datadiff(var,*args)
    for i,j in zip(Vn,Vf):
        expr=subsubs(expr,i,j)
    return diff(expr,t)    
    
def sympyEq2prime(expr,var,*args):
    '''
    return expr diff sympy expr whith diff prime symbols
    '''

    Vs,Vn,Vf,Vd,Vd2,Pd,Pd2=datadiff(var,*args)
    for i,j in zip(Vd2,Pd2):
        expr=expr.subs(i,j)
    for i,j in zip(Vd,Pd):
        expr=expr.subs(i,j)
    for i,j in zip(Vf,Vn):
        expr=expr.subs(i,j)   
    return expr

def symboldiff(*args):
    if len(args)==1:
        svar=str(args[0])
        dsvar=svar+"'"
        return symbols(dsvar)
    else:
        mm=[]
        for i in args:
            mm.append(symboldiff(i))
        return mm
    
def symboldiff2(*args):
    if len(args)==1:
        svar=str(args[0])
        dsvar=svar+"''"
        return symbols(dsvar)
    else:
        mm=[]
        for i in args:
            mm.append(symboldiff2(i))
        return mm 

def functiondiffk(expr,var,*args):
    '''
    if Area= b*h/2 in function of time...
    then A(t)=b(t)*h(t)/2...
    input :  functiondiff(expr, var, varF1, varF2,varF3...etc)
    in this case:
     A=functiondiff(A,t,A) = dA
     R=functiondiff(b*h/2,t,b,h) = b*dh/2 + db*h/2 
 
    '''
    if len(args)==0:
        return diff(expr,var)
    else:    
        v1=[]
        for i in args:
            v1.append(i)
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
            svar='d'+alphaname(i)
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