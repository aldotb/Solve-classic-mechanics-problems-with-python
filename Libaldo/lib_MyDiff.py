
from sympy import symbols
from functools import reduce
from IPython.display import  Math  ,display 
import matplotlib.pyplot as plt
from libaldo_math2 import *
#from libaldo_show import *
from  lib_MyEq import * 
from lib_MyEqEq import * 
import copy
s=symbols('s')

class MyDiff(MyEqEq):
    def __init__(self,expr1,expr2,v1,v2,ics=''):
        self.expr1=expr1
        self.expr2=expr2
        self.v1=v1
        self.v2=v2
        self.ics=ics
        self.tics=False
        self.eQics=''
        self.symb='='
         
        if ics!='':
            self.tics=True
            self.eQics=parse_expr(get_icsexp(ics,v1,v2))
        self.eQdiff=EqDiff(expr1,expr2,v1,v2,kshow=False)
        eQdiff=self.eQdiff    
        self.Ln=eQdiff.lhs
        self.Rn=eQdiff.rhs
        self.Lp=sympyEq2prime(self.Ln,v1,v2)
        self.Rp=sympyEq2prime(self.Rn,v1,v2)
        self.e1=MyEq(0,'e1',var=x,kshow=False)
        self.e2=MyEq(0,'e2',var=x,kshow=False)
        self.e1.ksym=self.Ln
        self.e2.ksym=self.Rn
        p1=self.Lp
        p2=self.Rp
        ps=self.symb
        display(Math(latex(p1)+' '+ps+' '+latex(p2)))
    def s(self):
        v1=self.v1
        v2=self.v2
        p1=sympyEq2prime(self.e1.ksym,v1,v2)
         
        p2=sympyEq2prime(self.e2.ksym,v1,v2)
         
        ps=self.symb
        display(Math(latex(p1)+' '+ps+' '+latex(p2)))

    def dsolve(self,*args):
        eQdiff=Eq(self.e1.ksym,self.e2.ksym)
        if self.tics:
            eQics=self.eQics
            ss=dsolve(eQdiff,ics=eQics)
        else:
            ss=dsolve(eQdiff)
        if len(args)>0:
            if type(ss)==list:
                vecres=[]
                for i in ss:
                    expr1=args[0]
                    expr2=i.rhs
                    vecres.append(MyEqEq(expr1,expr2))
                return vecres
            else:    
                expr1=args[0]
                expr2=ss.rhs
                if type(expr1)==str:
                    return MyEq(expr2,kname=expr1,var=self.v1)
                else: 
                    return MyEqEq(expr1,expr2)
        else:    
            return ss

class MyEqDiff(MyEqEq):
    def __init__(self,*args,ics='',vfunc=[],mark=True):
        
        self.type='MD'
        if type(args[0])==MyEqEq:
            expr=args[0]
            self.exp1=expr.L 
            self.exp2=expr.R
            self.ode=Eq(expr.L,expr.R)
            self.var=args[1]
            self.var1=args[2]
            self.func=Function(str(self.var1))(self.var)
            if len(args)==4:
                self.var2=args[3]
                self.func2=Function(str(self.var2))(self.var)
            self.e1=MyEq(self.exp1,'e1',kshow=False,vfunc=vfunc)
            self.e2=MyEq(self.exp2,'e2',kshow=False,vfunc=vfunc)    
            self.mark=mark
             
            self.ics2='' 
            self.equalityDiff=Eq(expr.L,expr.R)  
        else:    
                
            if len(args)>4:
                self.type=2
                self.sexp1,self.sexp2,self.var,self.var1,self.var2=args
             
            else:
                self.sexp1,self.sexp2,self.var,self.var1=args
                self.var2=''

            self.mark=mark
            self.ics=ics        
                
            self.equalityDiff=EqDiff(*args,kshow=False)
            self.exp1=self.equalityDiff.lhs
            self.exp2=self.equalityDiff.rhs 
            self.e1=MyEq(self.exp1,'e1',kshow=False,vfunc=vfunc)
            self.e2=MyEq(self.exp2,'e2',kshow=False,vfunc=vfunc)   
            var1=self.var1
            var2=self.var2
            var=self.var        
            self.func=Function(str(var1))(var)
            self.func1=Function(str(var1))(var)
            self.func2=Function(str(var2))(var)        
            self.ode=self.equalityDiff
            self.ics2=''
            if  ics!='':
                self.ics2=parse_expr(get_icsexp(self.ics,self.var,self.var1))
             
        if self.mark:
            self.dview()
        else:    
            display(Math(latex(self.equalityDiff)))
        
          
    def __call__(self,*args, **kwargs):
        eqDiff=self.equalityDiff
        display(Math(lates(eqDiff)))
        
    def __repr__(self):
        kres = self.equalityDiff
        return kres
        
    def _latex(self, obj):
        return latex(self.equalityDiff)    
    def __str__(self):
         
        return str(self.__repr__())
    
    @property
    def R(self):
        return self.equalityDiff.rhs
        
    @property
    def L(self):
        return self.equalityDiff.lhs 
        
    @property
    def right(self):
        return self.equalityDiff.rhs
        
    @property
    def left(self):
        return self.equalityDiff.lhs     
    

    def s(self):
        self.exp1=self.e1.ksym
        self.exp2=self.e2.ksym
        self.ode=Eq(self.exp1,self.exp2)
        self.equalityDiff=Eq(self.exp1,self.exp2)
        if self.mark:
            self.dview()
        else:    
            display(Math(latex(self.equalityDiff)))
    
    def dviewF(self,kret=False):
        var=self.var
        var1=self.var1
        var2=self.var2
        if var2!='':
            dexpr=easy_diffviewF(self.ode,var,var1,var2)
        else:    
            dexpr=easy_diffviewF(self.ode,var,var1)
        modiexp=  dexpr.replace('*','.')  
        display(Math(dexpr))
        if kret:
            return dexpr

    def dview(self):
        
        display(Math(latex(self.diff2mark())))
        
    def diff2mark(self):
        '''
            expr= Eq() equation equallity
            vx= independ variable x
            vy= dependent var y(x)
            
            diff2mark(diff1=diff2,y,x)
            
        '''
        self.e1.ksym=self.exp1
        self.e2.ksym=self.exp2
        expr=self.equalityDiff
        vx=self.var
        vy=symbols(str(self.var1))
        dy=symboldiff(self.var1)
        dy2=symboldiff2(self.var1)
        Y=Function(str(vy))(vx)

        expr=expr.subs(Y.diff(vx,vx),dy2)
        expr=expr.subs(Y.diff(),dy)
        expr=expr.subs(Y,vy)
        if self.var2!='':
             
            vz=self.var2
            dz=symboldiff(self.var2)
            dz2=symboldiff2(self.var2)
            Z=Function(str(vz))(vx)

            expr=expr.subs(Z.diff(vx,vx),dz2)
            expr=expr.subs(Z.diff(),dz)
            expr=expr.subs(Z,vz)
            
        return expr
        
    def convert2MQ(self):
        kres=[]
        Px=[]
        Py=[]

        var=self.var
        var1=self.var1
        var2=self.var2
        dx=diffsymbol(str(var1))
        dx2=diffsymbol2(str(var1))
        sdX2=str(var1)+"''"
        sdX =str(var1)+"'"

        sres= self.dview(kret=True,kshow=False) 
        if sdX2 in sres:
            sres=sres.replace(sdX2,str(dx2))
            Px.append(dx2)
        if sdX in sres:
            sres=sres.replace(sdX,str(dx)) 
            Px.append(dx)

        Px.append(dx)    
        if var2!='':
            dy2=diffsymbol2(str(var2))
            sdY2=str(var2)+"''"
            if sdY2 in sres:
                sres=sres.replace(sdY2,str(dy2)) 
                Py.append(dy2)
            dy=diffsymbol(str(var2))
            sdY=str(var2)+"'"
            if sdY in sres:
                sres=sres.replace(sdY,str(dy)) 
                Py.append(dy) 

        p1,p2=sres.split('=')        

        QQ=MQ(parse_expr(p1),parse_expr(p2))
        kres.append(QQ)
        if len(Px)>0:
            for i in Px:
                kres.append(i)
        if len(Py)>0:
            for i in Py:
                kres.append(i)
        
        return kres   




    def vecdata(self):
        vec=[]
        if self.type=='MD':
            vec.append([self.exp1,self.var,self.var1])
            vec.append([self.exp2,self.var,self.var1])
             
        else:
            vec.append([self.exp1,self.var,self.var1,self.var2])
            vec.append([self.exp2,self.var,self.var1,self.var2])
        return vec
    
         
            
    ###########################################
    #               Update                    #
    ###########################################

    def __add__(self, other):
        p1=self.L
        p2=self.R
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
        p1=self.L
        p2=self.R
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
        p1=self.L
        p2=self.R
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
        p1=self.L
        p2=self.R
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
        p1=self.L
        p2=self.R
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
        p1=self.L
        p2=self.R
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
        p1=self.L
        p2=self.R
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
        p1=self.L
        p2=self.R
        if type(other)==MyEqEq:
            p1=p1/other.left
            p2=p2/other.right
        elif type(other)==MyEq:
            p1=p1/other.ksym
            p2=p2/other.ksym

        else:
            p1=other/p1 
            p2=other/p2 
         
        return Eq(p1,p2)     
         
    def update_expr(self):
        self.exp1=self.equalityDiff.lhs
        self.exp2=self.equalityDiff.rhs
        self.e1=MyEq(self.exp1,'e1',kshow=False,vfunc=vfunc)
        self.e2=MyEq(self.exp2,'e2',kshow=False,vfunc=vfunc) 
     
    def update_equalityDiff(self):
        self.equalityDiff=Eq(self.exp1,self.exp2)
        
    def update_sexp(self):
        if self.type=='MD':
            vec=[self.exp1,self.var,self.var1]
            self.sexp1=str_diff(*vec)
            vec=[self.exp2,self.var,self.var1]
            self.sexp2=str_diff(*vec)
        else:
            vec=[self.exp1,self.var,self.var1,self.var2]
            self.sexp1=str_diff(*vec)
            vec=[self.exp2,self.var,self.var1,self.var2]
            self.sexp2=str_diff(*vec)
    def Add(self, kval,kname='',kshow=True ,kop='LR'):
         
        p1=self.e1.ksym
        p2=self.e2.ksym
         
        kval=traducediff(kval,self.var,self.var1) 	    
                 
        if 'L' in kop:
                p1=p1+kval
                self.e1.ksym=p1
        if 'R' in kop:
                p2=p2+kval
                self.e2.ksym=p2
        if kname!='':
            QQ=MyEqEq(p1,p2)
            return QQ
        else:
            self.exp1=p1
            self.exp2=p2
            self.update_equalityDiff()
            if kshow:
                self.s()
        
    def Substrac(self, kval,kname='', kop='RL',kshow=True ):
        p1=self.e1.ksym
        p2=self.e2.ksym
         
        kval=traducediff(kval,self.var,self.var1) 	    
                 
        if 'L' in kop:
                p1=p1-kval
                self.e1.ksym=p1
        if 'R' in kop:
                p2=p2-kval
                self.e2.ksym=p2
        if kname!='':
            QQ=MyEqEq(p1,p2)
            return QQ
        else:
            self.exp1=p1
            self.exp2=p2
            self.update_equalityDiff()
            if kshow:
                self.s()
        
    def Mul(self, kval,kname='', kop='RL',kshow=True ):
        p1=self.e1.ksym
        p2=self.e2.ksym
         
        kval=traducediff(kval,self.var,self.var1) 	    
                 
        if 'L' in kop:
                p1=p1*kval
                self.e1.ksym=p1
        if 'R' in kop:
                p2=p2*kval
                self.e2.ksym=p2
        if kname!='':
            QQ=MyEqEq(p1,p2)
            return QQ
        else:
            self.exp1=p1
            self.exp2=p2
            self.update_equalityDiff()
            if kshow:
                self.s()

    def Div(self, kval,kname='', kop='RL',kshow=True ):
        p1=self.e1.ksym
        p2=self.e2.ksym
         
        kval=traducediff(kval,self.var,self.var1) 	    
                 
        if 'L' in kop:
                p1=p1/kval
                self.e1.ksym=p1
        if 'R' in kop:
                p2=p2/kval
                self.e2.ksym=p2
        if kname!='':
            QQ=MyEqEq(p1,p2)
            return QQ
        else:
            self.exp1=p1
            self.exp2=p2
            self.update_equalityDiff()
            if kshow:
                self.s()

    def Pow(self, kval,kname='', kop='RL',kshow=True ):
        p1=self.e1.ksym
        p2=self.e2.ksym
         
        kval=traducediff(kval,self.var,self.var1) 	    
                 
        if 'L' in kop:
                p1=p1**kval
                self.e1.ksym=p1
        if 'R' in kop:
                p2=p2**kval
                self.e2.ksym=p2
        if kname!='':
            QQ=MyEqEq(p1,p2)
            return QQ
        else:
            self.exp1=p1
            self.exp2=p2
            self.update_equalityDiff()
            if kshow:
                self.s()

    def Rpow(self, kval=2,kname='', kop='RL',kshow=True ):
        p1=self.e1.ksym
        p2=self.e2.ksym
         
        kval=traducediff(kval,self.var,self.var1) 	    
                 
        if 'L' in kop:
                p1=p1**cfrac(1,kval)
                self.e1.ksym=p1
        if 'R' in kop:
                p2=p2**cfrac(1,kval)
                self.e2.ksym=p2
        if kname!='':
            QQ=MyEqEq(p1,p2)
            return QQ
        else:
            self.exp1=p1
            self.exp2=p2
            self.update_equalityDiff()
            if kshow:
                self.s()
        
    def showprimitive(self):
        vec1,vec2=self.vecdata()
        Eqstr=str_diff(*vec1)+'='+str_diff(*vec2)
        return Eqstr
    
    def solvediff(self,*args,kshow=True):
        if self.ics!='':
            kres= dsolve(self.equalityDiff,ics=parse_expr(get_icsexp(self.ics,self.var,self.var1)))
        else:
            kres= dsolve(self.equalityDiff)
        if kshow:
            display(Math(latex(kres)))        
        if 'F' in args:
            return kres.rhs
        if 'Eq' in args:
            ee=MyEq(kres.rhs,str(kres.lhs),var=self.var)
            return ee
        else:    
            return kres    
            
    def dsolve(self,*args):
    
        L=4		
        Id=''
        methW=False
        myeq=False
        kname=str(self.var1)
        for i in args:
            if i=='Wolfram':
                methW=True
            if type(i)==str and len(i)>10:
                Id=id
             
              
			
        if methW:
            return self.dsolveWolfram(id=Id)
        else:    
            if self.ics!='':
                kres=dsolve(self.ode,ics=self.ics)
                 
            else:
                kres=dsolve(self.ode)
                  
        if 'F' in args:
             
            ee=MyEq(kres.rhs,kname=str(self.var1),var=self.var)
            return ee
        elif type(kres)==list:
            for i in kres:
                display(Math(latex(i)))
            return ganswer(kres,'value')
        elif 'Eq' in args:
            QQ=MQ(self.var1,kres.rhs)
            return QQ
        else:
            return kres    
    def setL(self,expr):
        self.e1.ksym=expr
        self.s()
    def setR(self,expr):
        self.e2.ksym=expr
        self.s()    
        
    def set(self,swargs):
        P=swargs.split(',')
        p1=[]
        p1=[]
        for i in P:
            P2=i.split('=')
            p1.append(P2[0])
            p1.append(P2[1])
        sexp1=self.sexp1
        sexp2=self.sexp2        
        for i,j in zip(p1,p2):
            sexp1=sexp1.replace(i,j)
            sexp2=sexp2.replace(i,j)
        
        self.exp1=pru(sexp1,self.var1,self.var2)
        self.exp2=pru(sexp2,self.var1,self.var2)
        

        self.update_equalityDiff()
        self.s() 
    
    def replacediff2(self,expr):
        if type(expr)==MyEq:
            expr=expr.ksym
        var1=self.var1
        var=self.var
        F=Function(str(var1))(var)
        p1=self.exp1 
        p2=self.exp2  
        self.exp1=p1.subs( F.diff(var,var),expr)
        self.exp2=p2.subs( F.diff(var,var),expr)
        self.update_equalityDiff()
        self.s()


    def replacediff(self,expr):
        if type(expr)==MyEq:
            expr=expr.ksym
         
        var1=self.var1
        var=self.var
        F=Function(str(var1))(var)
        p1=self.exp1 
        p2=self.exp2  
        self.exp1=p1.subs( F.diff(var),expr)
        self.exp2=p2.subs( F.diff(var),expr)
        self.update_equalityDiff()
         
         
        self.s()
    
    def replacefunc(self,expr):
        if type(expr)==MyEq:
            expr=expr.ksym
         
        var1=self.var1
        var=self.var
        F=Function(str(var1))(var)
        p1=self.exp1 
        p2=self.exp2  
        self.exp1=p1.subs( F,expr)
        self.exp2=p2.subs( F,expr)
        self.update_equalityDiff()
         
         
        self.s()
    
    
    def simplifyexp(self,op='LR' ,kshow=True):
        p1=self.left
        if 'L' in op:
            p1=simplifyexp(p1)
            self.exp1.ksym=p1
        p2=self.right
        if 'R' in op:
            p2=simplifyexp(p2)
            self.exp2.ksym=p2

        if kshow:
            self.s()
    def basefactor(self,op='LR',kshow=True):
        return self.simplifyexp(op=op,kshow=kshow)
        
        
    def simplifyexp(self,op='LR',kshow=True):
        p1=self.left
        if 'L' in op:
            p1=simplifyexp(p1)
            self.exp1.ksym=p1
        p2=self.right
        if 'R' in op:
            p2=simplifyexp(p2)
            self.exp2.ksym=p2

        if kshow:
            self.s() 
    def div2mulexp(self,op='LR',kshow=True):
        p1=self.left
        if 'L' in op:
            p1=div2mulexp(p1)

            self.exp1.ksym=p1
            
        p2=self.right
        if 'R' in op:
            p2=div2mulexp(p2)
            
            self.exp2.ksym=p2

        if kshow:
            self.s()

    def reducePow(self,op='LR',kshow=True):
        p1=self.left
        if 'L' in op:
            p1=reducePow(p1)
            
            self.exp1.ksym=p1
            
        p2=self.right
        if 'R' in op:
            p2=reducePow(p2)
            
            self.exp2.ksym=p2

        if kshow:
            self.s()        
    def simplifyrpow(self,kop='RL',kshow=True):
        

        if 'L' in kop :
            self.exp1.simplifyrpow(kshow=False)
        if 'R' in kop :
            self.exp2.simplifyrpow(kshow=False)
        
        self.s()
    def simplify_cero(self,):
        kres=self.left-self.right
        kres=opemat(kres)
        self.exp1.ksym=kres
        self.exp2.ksym=0
        self.s()
        
    def expandexp(self,kop='LR',kshow=True):
        op=''
        if 'e' in kop:
            op='e'
        p1=self.left
        if 'L' in kop:
            p1=expandexp(p1,op=op)
            self.exp1.ksym=p1
        p2=self.right
        if 'R' in kop:
            p2=expandexp(p2,op=op)
            self.exp2.ksym=p2
         
        if kshow:
            self.s()
            
            
    def simplifybase(self,kop='LR',kshow=True):

        p1=self.left
        if 'L' in kop:
            p1=simplifybase(p1)
            self.exp1.ksym=p1
        p2=self.right
        if 'R' in kop:
            p2=simplifybase(p2)
            self.exp2.ksym=p2
         
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
            self.exp1.ksym=p1
        p2=self.right
        if 'R' in kop:
            p2=powexpand(p2,op=op)
            self.exp2.ksym=p2
         
        if kshow:
            self.s()

    def lexpand(self):
        self.e1.ksym=expand_log(self.e1.ksym,force=True)
        self.e2.ksym=expand_log(self.e2.ksym,force=True)
        self.s()
        
        
    def mulexpo(self,kop='LR',kshow=True):
        

        p1=self.left
        if 'L' in kop:
            p1=mulexpo(p1)
            self.exp1.ksym=p1
        p2=self.right
        if 'R' in kop:
            p2=mulexpo(p2)
            self.exp2.ksym=p2
         
        if kshow:
            self.s()
            
    def factor(self, kop='RL',kshow=True ):
        if 'L' in kop:
            self.exp1.factor()
        if 'R' in kop:
            self.exp2.factor()
        self.s()
    def simplifyexp(self,op='LR',kshow=True):
        p1=self.left
        if 'L' in op:
            p1=simplifyexp(p1)
            self.exp1.ksym=p1
        p2=self.right
        if 'R' in op:
            p2=simplifyexp(p2)
            self.exp2.ksym=p2

        if kshow:
            self.s()
    def basefactor(self,op='LR',kshow=True):
        return self.simplifyexp(op=op,kshow=kshow)
        
        
    def simplifyexp(self,op='LR',kshow=True):
        p1=self.left
        if 'L' in op:
            p1=simplifyexp(p1)
            self.exp1.ksym=p1
        p2=self.right
        if 'R' in op:
            p2=simplifyexp(p2)
            self.exp2.ksym=p2

        if kshow:
            self.s() 
    def div2mulexp(self,op='LR',kshow=True):
        p1=self.left
        if 'L' in op:
            p1=div2mulexp(p1)
            if kope!='':
                p1=opemat(p1)
            self.exp1.ksym=p1
            
        p2=self.right
        if 'R' in op:
            p2=div2mulexp(p2)
            if kope!='':
                p2=opemat(p2)
            self.exp2.ksym=p2

        if kshow:
            self.s()

    def reducePow(self,op='LR',kshow=True):
        p1=self.left
        if 'L' in op:
            p1=reducePow(p1)
            if kope!='':
                p1=opemat(p1)
            self.exp1.ksym=p1
            
        p2=self.right
        if 'R' in op:
            p2=reducePow(p2)
            if kope!='':
                p2=opemat(p2)
            self.exp2.ksym=p2

        if kshow:
            self.s()        
    def simplifyrpow(self,kop='RL',kshow=True):
        

        if 'L' in kop :
            self.exp1.simplifyrpow(kshow=False)
        if 'R' in kop :
            self.exp2.simplifyrpow(kshow=False)
        if kope!='':
            kres1=self.exp1.ksym
            kres2=self.exp2.ksym
            
            kres1=opemat(kres1)
            kres2=opemat(fcc)
            self.exp1.ksym=kres1
            self.exp2.ksym=kres2
        self.s()
    def expand(self, kop='RL',kshow=True):
        if 'L' in kop:
            self.exp1=expand(self.exp1)
            self.e1.ksym=self.exp1
            
        if 'R' in kop:
            self.exp2=expand(self.exp2)
            self.e2.ksym=self.exp2
        self.s()    
    def simplify_cero(self ):
        kres=self.left-self.right
        kres=opemat(kres)
        self.exp1.ksym=kres
        self.exp2.ksym=0
        self.s()

    def factorSec(self, ksym,kop='RL',kshow=True):
        if 'L' in kop:
            
            self.exp1.factorSec(ksym,kshow=False)
        if 'R' in kop:
            self.exp2.factorSec(ksym,kshow=False)
        self.update_equalityDiff()
        self.s(kshow)
    
    def toMyEqeEq(self):
        return MQ(self.exp1,self.exp2)
        
    def LaplaceEq(self):
        ee=LaplaceEq(self.sexp1,self.sexp2,self.var1,self.var,ics=self.ics)
        Ls=symbols('L_s')
        sLs=str(Ls)
        p1=ee.L
        sp1=str(p1)
        sp1=sp1.replace('L(s)',sLs)
        p1=parse_expr(sp1)
        p2=ee.R
        sp2=str(p2)
        sp2=sp2.replace('L(s)',sLs)
        p2=parse_expr(sp2)
        qq=MQ(p1,p2,kshow=False)
        return qq
    
    def transLaplace(self):
        return transLaplace(self)
         
      
 
F,f=symbols('F f')


def get_icsexp(*args):
    kres=args[0]
    var2=args[1]
    vecv=args[2:len(args)]
    print(vecv)
    Ssym=kres.replace('=',':')
    vec=Ssym.split(',')
    mm=[]
    for j in vecv:
        var1=j
        for i in vec:
            smm=i
            p1=i.find('(')
            p2=i.find(')')
            val=i[0:p1+1]

            if "''" in val:
                sP1=i[0:p2-1]

                sP2=str(var1)+'.diff('+str(var1)+'('+str(var2)+'), '+str(var2)+', '+str(var2)+').subs('+str(var2)+' ,'

                smm=smm.replace(sP1,sP2)

                Ssym=Ssym.replace(i,smm)
                #mm.append(smm) 
            if "'" in i:

                sP1=i[0:p2-1]

                sP2=str(var1)+'.diff('+str(var2)+').subs('+str(var2)+' ,'

                smm=smm.replace(sP1,sP2)  

                Ssym=Ssym.replace(i,smm)
                #mm.append(smm)  

def get_icskwargs(*args):
    kres=args[0]
    var2=args[1]
    vecv=args[2:len(args)]
    print(vecv)
    Ssym=kres
    vec=Ssym.split(',')
    mm=[]
    for j in vecv:
        var1=j
        for i in vec:
            smm=i
            p1=i.find('(')
            p2=i.find(')')
            val=i[0:p1+1]

            if "''" in val:
                sP1=i[0:p2-1]

                sP2=str(var1)+'.diff('+str(var1)+'('+str(var2)+'), '+str(var2)+', '+str(var2)+').subs('+str(var2)+' ,'

                smm=smm.replace(sP1,sP2)

                Ssym=Ssym.replace(i,smm)
                #mm.append(smm) 
            if "'" in i:

                sP1=i[0:p2-1]

                sP2=str(var1)+'.diff('+str(var2)+').subs('+str(var2)+' ,'

                smm=smm.replace(sP1,sP2)  

                Ssym=Ssym.replace(i,smm)
                #mm.append(smm)  

    
    Ssym=Ssym.replace(' ,',', ')
    Ssym='['+Ssym+']'
    return Ssym    
    # Laplace 
def LaplaceEq(kfunc,kfunc2,var1,var2,var3=s,ics='',kshow=True):
    

    
    L= Function('L')(var3)
    
    sL=str(L)
    kwtru=False    
    if ics!='':
        mList=ics.split(',')
        oldv=[]
        newv=[]
        kwtru=True
        for i in mList:
            kyt,vv=i.split('=')
            if str(var1)+"'(" in kyt:
                p1=kyt.find('(')
                p2=kyt.find(')')
                newt=kyt[p1+1:p2]
                kyt= 'Subs(Derivative('+str(var1)+'('+str(var2)+'), '+str(var2)+'), '+str(var2)+', '+newt+')'
                          
            oldv.append(kyt)
            newv.append(vv)
    P1=pru(kfunc,var1,var2)
    Lp1=laplace_transform(P1,var2, var3)
    if type(Lp1)==tuple:
        Lpp=Lp1[0]
        if 'Subs' in str(Lpp):
            Lp1=Lp1[1]
        else:    
            Lp1=Lp1[0]
         
            
       
    sL1=str(Lp1)
     
    
    P2=pru(kfunc2,var1,var2)
    Lp2=laplace_transform(P2,var2, var3)
     
    
    if type(Lp2)==tuple:
         
        Lpp=Lp2[0]
        if 'Subs' in str(Lpp):
            Lp2=Lp2[1]
        else:    
            Lp2=Lp2[0]
    sL2=str(Lp2)
      
    sL2=sL2.replace('1','1/s')
    sL2=sL2.replace('u','exp(-s)/s')
    sL=str(L)
    sLP='LaplaceTransform('+str(var1)+'('+str(var2)+'), '+str(var2)+', '+str(var3)+')' 
    sL1=sL1.replace(sLP,sL)
    sL2=sL2.replace(sLP,sL)
    if kwtru:
        for i,j in zip(oldv,newv):
            sL1=sL1.replace(i,j)
            sL2=sL2.replace(i,j)
    Lp1=parse_expr(sL1)
    #Lp1=clean_LC(Lp1,var1,var2)
    
    Lp2=parse_expr(sL2)
    #Lp2=clean_LC(Lp2,var1,var2)
    return MyEqEq(Lp1,Lp2,var1=var1,var2=var2,kshow=kshow)


s=symbols('s')
def inversa_Ls(Q,kname='',omega=False):
     
    if type(Q)==MyEq:
        s = Q.var 
        t = Q.var1
        kres=inverse_laplace_transform(Q.ksym, s, t)
           
    else:
        kres=inverse_laplace_transform(Q, var, var1)
    
    if omega==False:
        sres = str(kres)
        kfind = 'Heaviside('+str(t)+')' 
        sres = sres.replace(kfind,'1')
        
        kres = parse_expr(sres)
        try:
            kres = apart(expand(simplify(kres)))
        except:
            pass
    if kname=='':
        return kres
    else:    
        nname=str(kname)   
        return MyEq(kres,nname,ktype='F')
        
		 
			
def pru(kfunc,var1,var2,var3=''):
    sfunc=str(kfunc)
    F=Function(str(var1))(var2)
    svar=str(var1)
    sdvar=svar+"'"
    sd2var=svar+"''"
    sF=str(var1)+'('+str(var2)+')'
    sdiff=str(diff(F,var2))
    s2diff=str(diff(F,var2,var2))
    sfunc=sfunc.replace(sd2var,s2diff)
    sfunc=sfunc.replace(sdvar,sdiff)
    return parse_expr(sfunc)


def creaF(kfunc,kfunc2,var1,var2):
    p1=pru(kfunc,var1,var2)
    p2=pru(kfunc2,var1,var2)
     
    return  p1,p2
 
def clean_LC(kres,var1,var2):
    sres=str(kres)
    clean1='Subs(Derivative('+str(var1)+'('+str(var2)+'), '+str(var2)+'), '+str(var2)+', 0)'
    sres=sres.replace(clean1,'0')
    return parse_expr(sres)
        
def solve_LC(obj,var=s):
    kname='L(s)'
    if type(obj)==MyEqEq:
        kres=obj.L-obj.R
    elif type(obj)==MyEq:
        kres=obj.ksym
    else:
        kres=obj
        
    ee=MyEq(kres,'ee',kshow=False)
    ss=ee.solve(kname,kshow=False)
    ss.expand(kshow=False)
    ss.simplify(kshow=False)
    ss.factor(kshow=False)
    try:
        ss.ksym=apart(ss.ksym)
    except:
        pass
    #ss=apart(ss)
    Ls=symbols('Ls')
     
    return MQ(Ls,ss.ksym)  
    
def func_diff(var1,var2):
    try:
        var1=Function(var1.name)
    except:
        vname=str(var1)
        var1=Function(vname)
        
    return var1(var2).diff(var2)
    
def fdiff(var1,var2):
    try:
        var1=Function(var1.name)
    except:
        vname=str(var1)
        var1=Function(vname)
        
    return var1(var2).diff(var2)


    
def func_diff2(var1,var2):
    var1=Function(var1.name) 
    return var1(var2).diff(var2,var2)
def fdiff2(var1,var2):
    var1=Function(var1.name) 
    return var1(var2).diff(var2,var2)   
    
    
    
def func_var2diff(var1,var2):
    var1=Function(var1.name)(var2)
    return var1
def newfunc(var1,var2):
    var1=Function(var1.name)(var2)
    return var1
    
    
def eQdiff(var1,var2,exp1,exp2):
    df='d'+str(var1)
    df= func_diff(var1,var2)
    df2='d'+str(var1)+'2'
    df2= func_diff2(var1,var2)
    f=func_var2diff(var1,var2)
     
    v2='d'+str(var1)+'2'
    exp1=exp1.replace(v2,str(df2))
    exp2=exp2.replace(v2,str(df2))
    v1='d'+str(var1)
    exp1=exp1.replace(v1,str(df))
    exp2=exp2.replace(v1,str(df))
    kres1=parse_expr(exp1)
    kres2=parse_expr(exp2)
    return kres1,kres2

def eQics(var1,var2,exp1):
    df='d'+str(var1)
    df= func_diff(var1,var2)
    df2='d'+str(var1)+'2'
    df2= func_diff2(var1,var2)
    f=func_var2diff(var1,var2)
     
    v2='d'+str(var1)+'2'
    exp1=exp1.replace(v2,str(df2))
    exp2=exp2.replace(v2,str(df2))
    v1='d'+str(var1)
    exp1=exp1.replace(v1,str(df))
    exp2=exp2.replace(v1,str(df))
    kres1=parse_expr(exp1)
    kres2=parse_expr(exp2)
    return kres1,kres2
    
def get_ics(*args):
    fd1=func_var2diff(var1,var2)
     
    if svar[0]=='d':
        if svar[1]==str(var1):
            if svar[2]=='(':
                p2=svar.find(')')
                val1=svar[3:p2]
                val2=svar[p2+1::]
                val3='diff('+str(var1)+'('+str(var2)+'),'+str(var2)+').subs('+str(var2)+','+val1+')'+val2
            else:
                p2=svar.find(')')
                val1=svar[4:p2]
                val2=svar[p2+1::]
                val3='diff('+str(var1)+'('+str(var2)+'),'+str(var2)+','+str(var2)+').subs('+str(var2)+','+val1+')'+val2
             
            return val3
    else:
        return svar    


    ###  Wolfram Alpha 

def pyDiff2wolframDiff(expr,v1,v2): 
    r""" 
    translate sympy differential  equation sintaxis
    to wolgram alpha mathemayicas sintaxis in order
    can use in mathemayica API
    why..??  because sympy sometimes can not solve differential ecuation
    example:
        Derivative(y(x), (x, 2)) + Derivative(y(x), x) + y + x return
                   y''+y'+y+x
                   
        Input(expresi diferential equatin, var depen, var indep.
        
    """
    sres=str(expr)
    
    oexp='Derivative('+str(v1)+'('+str(v2)+'), ('+str(v2)+', 2))'
    nexp=str(v1)+"''"
    sres=sres.replace(oexp,nexp)
    
    oexp='Derivative('+str(v1)+'('+str(v2)+'), '+str(v2)+')'
    nexp=str(v1)+"'"
    sres=sres.replace(oexp,nexp)
    
    oexp=str(v1)+'('+str(v2)+')'
    nexp=str(v1)
    sres=sres.replace(oexp,nexp)
    
    oexp='**'
    nexp='^'
    
    sres=sres.replace(oexp,nexp)
    return sres

   
def get_icsexp(sexpp,var1,var2):
    xx=symbols(str(var1))
    ff=Function(str(var2))(var1)
    kres=sexpp
    
    vecv=[var1]
     
    Ssym=kres.replace('=',':')
    vec=Ssym.split(',')
    mm=[]
    for j in vecv:
        
        
        for i in vec:
            smm=i
            p1=i.find('(')
            p2=i.find(')')
            val=i[0:p1+1]

            if "''" in val:
                sP1=i[0:p2-1]
                 
                sP2='Derivative('+str(var2)+'('+str(var1)+'), '+str(var1)+' ,'+str(var1)+').subs('+str(var1)+', '

                smm=smm.replace(sP1,sP2)

                Ssym=Ssym.replace(i,smm)
                #mm.append(smm) 
            elif "'" in i:

                sP1=i[0:p2-1]

                sP2='Derivative('+str(var2)+'('+str(var1)+'), '+str(var1)+').subs('+str(var1)+', '

                smm=smm.replace(sP1,sP2)  

                Ssym=Ssym.replace(i,smm)
            else:
                sP1=i[0:p2-1]

                sP2=str(var2)+'('+str(var1)+').subs('+str(var1)+','

                smm=smm.replace(sP1,sP2)  

                Ssym=Ssym.replace(i,smm)
                

    Ssym=Ssym.replace(' ,',', ')
    Ssym='{'+Ssym+'}'
    return Ssym
             
def traduce_Lc_Diff(expr,var1,var2):
    oldd2='d'+str(var1)+'2'
    newd2=str(var1)+"''"
    
    expr=expr.replace(newd2,oldd2)
    
    oldd1='d'+str(var1) 
    newd1=str(var1)+"'"
    
    expr=expr.replace(newd1,oldd1)
    return expr
    
def traduce_diff_Lc(expr,var1,var2):
    newd2='d'+str(var1)+'2'
    oldd2=str(var1)+"''"
    
    expr=expr.replace(newd2,oldd2)
    
    newd1='d'+str(var1) 
    oldd1=str(var1)+"'"
    
    expr=expr.replace(newd1,oldd1)
    return expr  


 

def str_diff(*args):
    expr=str(args[0])
    t=args[1]
    vecvar=args[2:len(args)]
    for i in vecvar:
        f=Function(str(i))(t)
        sd2f=str(diff(f,t,t))
        sd1f=str(diff(f,t))
        snd2=str(i)+"''"
        snd1=str(i)+"'"
        expr=expr.replace(sd2f,snd2)
        expr=expr.replace(sd1f,snd1)
    return expr    
        

def EqDiff(*args,kshow=True):
    newa1=[args[0]]
    for i in range (2,len(args)):
        newa1.append(args[i])
    newa2=args[1:len(args)]
    p1=traducediff(*newa1)
    p2=traducediff(*newa2)
    kres=Eq(p1,p2)
    if kshow:
        display(Math(latex(kres)))
    return kres
   
    def toMyEqeEq(self):
        return MQ(self.exp1,self.exp2)
   
def traduce_single_ics(sval,v1,v2):
    partes=sval.split('=')
    p2=partes[1]
    ktype=0
    p1=partes[0]
    if "''" in p1:
        ktype=2
    elif"'" in p1:
        ktype=1
    else:
        ktype=0
    P1=p1.find("(")
    P2=len(p1)-1
    vnum=p1[P1+1:P2]
    
    X=str(v2)
    T=str(v1)
    if ktype==2:
        sres= 'Subs(Derivative('+X+'('+T+'), '+T+','+T+'), '+T+', '+vnum+')'
    elif ktype==1:
        sres= 'Subs(Derivative('+X+'('+T+'), '+T+'), '+T+', '+vnum+')'
    else:
        sres=p1
        
    return sres,p2
##  LAPLACE

def tLaplace(Obj):
    L=symbols('L')
    mm=laplace_transform(Obj.exp1-Obj.exp2,Obj.var,s)
    smm=str(mm)
    X=str(Obj.var1)
    T=str(Obj.var)
    sLap='LaplaceTransform('+X+'('+T+'), '+T+', s)'
    smm=smm.replace(sLap,'L')
    if Obj.ics!='':
        sics=Obj.ics
        partes=sics.split(',')
        for i in partes:
            p1,p2=traduce_single_ics(i,Obj.var,Obj.var1)
            smm=smm.replace(p1,p2)
    smm=smm[1::]
    pp=smm.find(',')
    smm=smm[0:pp] 
    kres=parse_expr(smm)

    return kres
def laplace(*args):
    var=t
    expr=args[0]
    kname=latex2sympy('\mathcal{L}')
    if len(args)==1:
        if type(expr)==MyEqDiff:
            kres=tLaplace(expr)
            if 'L' in str(kres):
                L=symbols('L')
                L=csolve(kres,L)
                return L
            else:
                return kres
                 
        elif type(expr)==MyEqEq:
                f1=expr.L
                f2=expr.R
                var=expr.var
                P=laplacetransform(f1-f2,var)

        elif type(expr)==MyEq:
            f1=expr.ksym            
            var=expr.var
            P=laplacetransform(f1,var)

        else:
            if not 't' in str(expr) and 'x' in str(expr):
                var=x
            f1=expr
            P=laplacetransform(f1,var)
             
    else:
        expr=args[0]
        var=args[1]
        P=laplacetransform(expr,var)
     
    return P
def Laplace(*args):
    kres=laplace(*args)
    P=MQ(latex2sympy('\\mathcal{L}'),kres)
    return  kres


    
def transLaplace(Obj,kshow=True):
    L=symbols('L')
    mm=laplace_transform(Obj.exp1-Obj.exp2,Obj.var,s)
    smm=str(mm)
    X=str(Obj.var1)
    T=str(Obj.var)
    sLap='LaplaceTransform('+X+'('+T+'), '+T+', s)'
    smm=smm.replace(sLap,'L')
    if Obj.ics!='':
        sics=Obj.ics
        partes=sics.split(',')
        for i in partes:
            p1,p2=traduce_single_ics(i,Obj.var,Obj.var1)
            smm=smm.replace(p1,p2)
    '''
    if smm[0]=='(':
        smm=smm[1::]
    pp=smm.find(',')
    smm=smm[0:pp] 
    '''
    kres=parse_expr(smm)
    eL=MyEq(kres,'eL',var=L,ktype='LA')
      
    return eL    
def diff2laplace(Obj,kshow=True):
    L=symbols('L')
    mm=laplace_transform(Obj.exp1-Obj.exp2,Obj.var,s)
    smm=str(mm)
    X=str(Obj.var1)
    T=str(Obj.var)
    sLap='LaplaceTransform('+X+'('+T+'), '+T+', s)'
    smm=smm.replace(sLap,'L')
    if Obj.ics!='':
        sics=Obj.ics
        partes=sics.split(',')
        for i in partes:
            p1,p2=traduce_single_ics(i,Obj.var,Obj.var1)
            smm=smm.replace(p1,p2)
      
    if smm[0]=='(':  
        smm=smm[1::]
        pp=smm.find(',')
        smm=smm[0:pp] 
     
    kres=parse_expr(smm)
    eL=MyEq(kres,'eL',var=L,ktype='LA',kshow=False)
    kres=eL.solve(L,kshow=False)
    kres.var=s
    P=MQ(latex2sympy('\\mathcal{L}'),kres)
    return kres
      
          
def laplacetransform(func,var=t,var2=s):
    kres=laplace_transform(func,var,var2)
    try:
        return kres[0]
    except:
        return kres

def func2laplace(*args,kshow=True):
    var=t
    expr=args[0]
    kname=latex2sympy('\mathcal{L}')
    if len(args)==1:
        if type(expr)==MyEqDiff:
            f1=expr.L
            f2=expr.R
            var=expr.var
            ee=MyEq(f1-f2,'ee',var=var,kshow=False)
            return func2laplace(ee) 
           
             
        elif type(expr)==MyEqEq:
            f1=expr.L
            f2=expr.R
            var=expr.var
            P=laplacetransform(f1-f2,var)
             
        elif type(expr)==MyEq:
            f1=expr.ksym            
            var=expr.var
            P=laplacetransform(f1,var)
                 
        else:
            if not 't' in str(expr) and 'x' in str(expr):
                var=x
            f1=expr
            P=laplacetransform(f1,var)
             
    else:
        expr=args[0]
        var=args[1]
        P=laplacetransform(expr,var)
    if kshow:
        MQ(kname,P)
    return P   

def translaplace(expr,*args):
    Lp=func2laplace(expr)
    kname=''
    vecold=['multt']
    vecnew=[]
    for i in args:
        if i in vecold:
            vecnew.append(i)
        else:
            kname=i 
    if kname=='':
        return 

def ilaplace(expr,var=s,var1=t,omega=False):

    kres=inverse_laplace_transform(expr, var, var1)
    
    if omega==False:
        sres = str(kres)
        kfind = 'Heaviside('+str(t)+')' 
        sres = sres.replace(kfind,'1')
        
        kres = parse_expr(sres)
        try:
            kres = apart(expand(simplify(kres)))
        except:
            pass
    return kres 
    
    
def Ilaplace(expr,var=s,var1=t,omega=False):
    if type(expr)==MyEq:
        expr=expr.ksym
    expr= ilaplace(expr,var=var,var1=var1,omega=omega)  
    ee=MyEq(expr,'F('+str(var1)+')',var=var1)
    return ee

def dsolvesys(Q1,Q2):
    ode=[Q1.ode, Q2.ode]
    func= [Q1.func1,Q2.func1] 
    ics3=''
    if Q1.ics2!='':
        ics3=Q1.ics2
        if Q2.ics2!='':
            ics3={**Q1.ics2,**Q2.ics2}
    if Q2.ics2!='':
        ics3=Q2.ics2
        if Q1.ics2!='':
            ics3={**Q2.ics2,**Q1.ics2} 
     
    if ics3!='':
        kres= dsolve(ode, func,ics=ics3)
    else:
        kres= dsolve(ode, func)
    ff1=kres[0].lhs
    kres1=kres[0].rhs
    var=Q1.var
    ff2=kres[1].lhs
    kres2=kres[1].rhs
    kname1=str(ff1)
    kname2=str(ff2)
    ee1=MyEq(kres1,kname=kname1,var=var)
    ee2=MyEq(kres2,kname=kname2,var=var)
    return ee1,ee2
    

#####################
# view from regular Derivative(x(t)) to x'

def easy_diffviewF(expr,var,var1,var2=''):
    sexpr=preeasyF(expr,var=var,var1=var1)
    if var2!='':
        sexpr=preeasyF(sexpr,var=var,var1=var2)
    if 'Eq(' in sexpr:
        sexpr=sexpr[3:-1]
        sexpr=sexpr.replace(","," = ")
      
    return sexpr 

def easy_diffview(expr,var,var1,var2=''):
    sexpr=preeasy(expr,var=var,var1=var1)
    if var2!='':
        sexpr=preeasy(sexpr,var=var,var1=var2)
    if 'Eq(' in sexpr:
        sexpr=sexpr[3:-1]
        sexpr=sexpr.replace(","," = ")
     
    return sexpr

def preeasyF(expr,var,var1):
     
    X=str(var1)
    T=str(var)
    F=X+'('+T+')'
    odX ='Derivative('+F+', '+T+')'
     
    odX2='Derivative('+F+', ('+T+', 2))'
    ndX=X+"'"
    ndX2=X+"''"
    
    sexpr=str(expr)
    sexpr=sexpr.replace(odX2,ndX2)
    sexpr=sexpr.replace(odX,ndX)
    
    return sexpr
    return sexpr

def preeasy(expr,var,var1):
     
    X=str(var1)
    T=str(var)
    F=X+'('+T+')'
    odX ='Derivative('+F+', '+T+')'
     
    odX2='Derivative('+F+', ('+T+', 2))'
           
    ndX=X+"'"
    ndX2=X+"''"
    
    sexpr=str(expr)
    sexpr=sexpr.replace(odX2,ndX2)
    sexpr=sexpr.replace(odX,ndX)
    sexpr=sexpr.replace(F,X)

    return sexpr
    return sexpr

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





def subsdiff(expr,vy,vx):
    dy=symboldiff(y)
    dy2=symboldiff2(y)
    Y=Function(str(vy))(vx)

    expr=expr.subs(Y.diff(vx,vx),dy2)
    expr=expr.subs(Y.diff(),dy)
    expr=expr.subs(Y,vy)
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

def diff_withrespect(*args):
    qq=len(args)
    expr=args[0]
    sres=str(expr)
    var=args[1]
    if len(args)>2:
        vecf=args[2:qq]
        svar=[str(i) for i in vecf]
        sfun=[str(i)+'('+str(var)+')' for i in vecf]
        for i,j in zip(svar,sfun):
            sres=sres.replace(i,j)
    kres=parse_expr(sres)
    return kres.diff(var)
    
        
def Qdiff_withrespect(self,*args):
    p1=self.L
    args1=[p1]
    for i in args:
        args1.append(i)
    kres1=diff_withrespect(*args1)
    p2=self.R
    args2=[p2]
    for i in args:
        args2.append(i)
     
    kres2=diff_withrespect(*args2)
    QQ=MQ( kres1,kres2,kshow=False)
   
    argst=[QQ]
    for i in args:
        argst.append(i)
    kres= MyEqDiff(*argst)
    return kres

def strdiff(v1,v2):
    v1=Function(str(v1))(v2)
    return str(v1.diff(v2))
def strdiff2(v1,v2):
    v1=Function(str(v1))(v2)
    return str(v1.diff(v2,v2)) 



#  remplazar Icss in LaPlace

def vecics(sexpr,kss):
    return sexpr.split(kss)
def typediff(sexpr):
    if "''" in sexpr:
        return 2
    elif"'" in sexpr:
        return 1
    else:
        return 0
def subsdiffexpr(expr,vt,vx,ics1,ics2):  #,ics1='x(0)',ics2=0
    sexp=str(expr)
    ktype=typediff(ics1)
    if ktype==0:
        sexp=sexp.replace(ics1,ics2)
    elif ktype==1:
        xx=str(vx)
        tt=str(vt)
        olde='Subs(Derivative('+xx+'('+tt+'), '+tt+'), '+tt+', '+'0'+')'
        sexp=sexp.replace(olde,ics2)
    else:
        pass
    return  parse_expr(sexp) 

def subsics(expr,vt,vx,icss):
    '''
    expr= 4*s**2*LaplaceTransform(y(t), t, s) - 4*s*y(0) + LaplaceTransform(y(t), t, s) - 
        4*Subs(Derivative(y(t), t), t, 0) + 2/s
    ics=x(0)=0,x'(0)=0 
    subsics(Ls,t,y,Q.ics)   return   4*s**2*LaplaceTransform(y(t), t, s) + LaplaceTransform(y(t), t, s) + 2/s

        
    '''
    sexpr=str(expr)
    mm=vecics(icss,',')
    for i in mm:
        ics1,ics2=vecics(i,'=')
        expr=subsdiffexpr(expr,vt,vx,ics1,ics2)
    return expr  


def solvelaplace(obj,v1=t,v2=x):
    '''
    obj=4*s**2*LaplaceTransform(y(t), t, s) + LaplaceTransform(y(t), t, s) - 2 + 2/s
    Lt=solvelaplace(Ls,t,y)  return Lt= '(2*s - 2)/(4*s**3 + s)'
    '''
    
    
    tt=str(v1)
    xx=str(v2)
    Lsexp='LaplaceTransform('+xx+'('+tt+'), '+tt+', s)'
    V=parse_expr(Lsexp)
    if type(obj)==MyEqEq:
        kres=obj.L-obj.R
    elif type(obj)==MyEq:
        kres=obj.ksym
    else:
        kres=obj 
    ee=MyEq(obj,'ee',var=s,kshow=False)
    Lt=ee.solve(V,kshow=False)
    ee2=MyEq(Lt,'L_t',var=s)
    return ee2
    
#  laplace of Mul de tLaplace
def Is_tPow(expr): # True if expr=t ,t**2,t**3.... else false
    if expr==t:
        return True
    elif Is_Pow(expr):
        bb=getbase(expr)
        if bb==t:
            return True
        else:
            return False
    else:
        return False    
 
def tpow(expr): # get monomie t**n in a Mul Polynimie retru t, expo of t and thes rest of multwithout t
    mm=expr.args
     
    for i in mm:
        if  Is_tPow(i) or  i=='t':
            kres=simplify(expr/i)
            if i==t:
                return t,1,kres
            else:
                ee=getexpo(i)
                return t,ee,kres 
                
def translaplacemult(expr): # Laplacen trans of mult t**n
    if type(expr)==MyEq:
        expr=expr.ksym
    tt,ee,expr2=tpow(expr)
    Lp=laplace(expr2)
    ee=MyEq(Lp,'ee',var=s,kshow=False)
    eed=ee.diff()
    return eed

def translaplacedivt(expr): # Laplacen trans of mult t**n
    if type(expr)==MyEq:
        kres=expr.ksym
    expr2=kres*t
    kres=laplace(expr2)
 
    kres2= simplify(integrate(kres,s))

    return kres2    
    
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


def eQlaplace(*args):
    var=t
    kname=''
    expr=args[0]
    for i in args[1::]:
        if type(i)==str:
            kname=i
        else:
            var=i
    func=expr        
    if type(expr)==MyEq:        
        func=expr.ksym
    kres=translaplacemult(func)
    if kname=='':
        return kres
    else:
        return MyEq(kres,kname=kname,var=s)
        
def eQinverselaplace(Q,kname='',omega=False):
    s,t=symbols('s t') 
    if type(Q)==MyEq:
        kres=inverse_laplace_transform(Q.ksym, s, t)       
    else:
        kres=inverse_laplace_transform(Q, s, t)
    
    if omega==False:
        sres = str(kres)
        kfind = 'Heaviside('+str(t)+')' 
        sres = sres.replace(kfind,'1')
        
        kres = parse_expr(sres)
        try:
            kres = apart(expand(simplify(kres)))
        except:
            pass
    if kname=='':
        return kres
    else:    
        return MyEq(kres,kname=kname,var=t)        