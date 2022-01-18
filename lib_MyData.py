from sympy import *
from libaldo_show import *
from IPython.display import display, Math  

import numpy as np
from libaldo_math2 import *
from lib_MyEq  import *
from libaldo_show import *
   
class MyData:
    def __init__(self,vmain=[],vsolve=[],vsecu=[],datasys='',kope='',**kwargs):
        if datasys!='':
            vmain=[]
            vsolve=[]
            for i in datasys:
                vmain.append(parse_expr(i.name))
                vsolve.append(opemat(i.ksym,kope=kope))
            self.vmain=vmain
            self.vsolve=vsolve
                
        elif len(kwargs)>0:
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
        

    def make_vsecu(self):
        mm=self.kList()
        vmain2=self.vmain
        vsolve2=self.vsolve
        vsecu2=[]
        for i in mm:
            if i not in vmain2 and i not in vsolve2:
                vsecu2.append(i)
                
        self.vsecu=vsecu2
    def get_vsecu(self):
        self.make_vsecu()
        return self.vsecu
    def add_primaria(self,kval):
        pp=self.vmain
        if type(kval)=='list':
            pp=pp+kval 
        else:    
            pp.append(kval)
        self.vmain=pp
        return pp
        
        
    def add(self,vksym,vkval,krevalue=True):
        self.add_data(vksym=vksym,vkval=vkval,krevalue=krevalue)
        
    def store(self,vksym,vkval,krevalue=True,kshow=True):
        if type(vksym)!=list:
            ksym=[vksym]
            kval=[vkval]
        else:
            ksym= vksym 
            kval= vkval
            
        for i,j in zip(ksym,kval):
            if Is_Number(j):
                self.add_data(i,j,kshow=False)
            else:    
                self.addS(i,j,kshow=False)
        if kshow:    
            self.s()    
        
    def upTriang(self,T3,kang,kshow=True):
        if kang not in self.vsecu:
            vsym=[sin(kang),cos(kang),tan(kang)]
            vval=[T3.sin(),T3.cos(),T3.tan()]
            for i,j  in zip(vsym,vval):
                self.vmain.append(i)
                self.vsolve.append(j)
            self.vsecu.append(kang)
        if kshow:
            self.s()
        
    def add_data(self,vksym,vkval,krevalue=True,kshow=True):
        
        if type(vksym)!=list:
            ksym=[vksym]
            kval=[vkval]
        else:
            ksym= vksym 
            kval= vkval 
            
        Ms=self.dataS
        
        Mv=self.dataV
        qqq=self.qq
        for i,j in zip(ksym,kval):   
         
            Ms.append(i)
            Mv.append(j)
            qqq+=1
            
        self.dataS=Ms

        self.dataV=Mv
        

        self.qq=qqq
        self.sort_val()
       
         
         
        
        if krevalue:
            self.revalue(kshow=kshow)
        else:
            for i,j in zip(self.dataS,self.dataV):
                sE([i,' = ',j])

    def addS(self,vksym,vkval,kshow=True):
        Ms=self.dataS
        Mv=self.dataV
        Mq=self.qq
        Mvmain=self.vmain
        Mvsecu=self.vsecu
        
        Ms.append(vksym)
        Mv.append(vkval)
        Mq+=1
         
          
        if Is_dependient(vkval,Mvmain):
                if vksym not in Mvsecu:
                    Mvsecu.append(vksym) 
        self.dataS=Ms
        self.dataV=Mv
        self.qq=Mq
        self.vsecu=Mvsecu
            
        self.fix_prim_val()
        if kshow:
            self.s()

    def fix_secu(self):
        Ms=self.dataS
        Mv=self.dataV
        Mq=self.qq
        Mvmain=self.vmain
        Mvsecu=self.vsecu
        for i,j in zip(Ms,Mv):
            if i not in Mvmain:
                if not Is_Number(j):
                    if Is_dependient(j,Mvmain):
                        if i not in Mvsecu:
                            Mvsecu.append(i) 
        self.vsecu=Mvsecu

    def fix_prim_val(self):
        self.fix_secu()
        Ms=self.dataS
        Mv=self.dataV
        Mq=self.qq
        Mvmain=self.vmain
        Mvsecu=self.vsecu
        for i in range(Mq):
            if Ms[i] not in Mvmain:
                if Ms[i] not in Mvsecu:
                    for k in Mvsecu:
                        kres1=self.get_value(k)
                        kres2=Mv[i]
                        kres3=kres2.subs(k,kres1)
                        Mv[i]=kres3
                        
        self.fix_secu()                
                        
    def show_data(self):
        for i,j in zip(self.vmain,self.vsolve):
            sE([i,' = ',j])  
            
    def s(self):
        self.show_data()
    def revalue(self,kshow=True):
        self.sort_val()
        qqq = self.qq 
        if qqq>1:
            for i in range(qqq-1,0,-1):
                ks=self.dataS[i]
                kv=self.dataV[i]
                for j in range(0,i):
                    try:
                        oldv=self.dataV[j]
                        newv=oldv.subs(ks,kv)
                        newv=strSubs(newv,ks,kv)
                        self.dataV[j]=newv
                    except:
                        done=True
        if kshow:
            self.s()
    def get_value(self,kval):
            kres=kval
            for i,j in zip(self.dataS,self.dataV):
                if i==kres:
                    return j    
            return(kres)
    def clean_idems(self):
        Ms=self.dataS
        Mv=self.dataV
        qqq=0
        Ms2,Mv2=[],[] 
        for i,j in zip(Ms,Mv):
            if i!=j:
                Ms2.append(i)
                Mv2.append(j)
                qqq+=1
            
        Ms= Ms2 
        Mv= Mv2 
        self.dataS=Ms
        self.dataV=Mv
        self.qq=qqq
        
    def sort_val(self):
        Ms=self.dataS
        Mv=self.dataV
        
        Ms1,Ms2,Mv1,Mv2=[],[],[],[]
        self.make_vsecu()
        contenedor=self.vmain+self.vsecu
        for i,j in zip(Ms,Mv):
            if Is_Number(j) or j in contenedor or Is_dependient(j,contenedor):
                Ms1.append(i)
                Mv1.append(j)
            else:
                Ms2.append(i)
                Mv2.append(j)
        Ms= Ms2+Ms1
        Mv= Mv2+Mv1
        self.dataS=Ms
        self.dataV=Mv
    def kList(self):
        qqq=self.qq
        Ms=self.dataS 
        Mv=self.dataV 
        mm=[]
        for i in range(qqq):
            vv1=fpoly(Ms[i],'free')
            for j in vv1:
                if j not in mm:
                    mm.append(j)
            vv2=fpoly(Mv[i],'free')
            for j in vv2:
                if j not in mm:
                    mm.append(j)        
        return mm
    def recheck_cate(self):
        mV=self.dataV 
        mS=self.dataS 
        mQ=self.qq 
        mM=self.vmain
        for i,j in zip(mS,mV):
            if Is_dependient(j,mM):
                if i not in mM:
                    mM.append(i)
            
        self.vmain=mM
    def simplify(self):
        kres=[]
        Mv=self.dataV
        for i in Mv:
            kres.append(simplify(i))
        self.dataV=kres
        self.sort_val()
        self.show_data()
        
    def expand(self):
        kres=[]
        Mv=self.dataV
        for i in Mv:
            kres.append(expand(i))
        self.dataV=kres
        self.sort_val()
        self.show_data()    
        
    def factor(self):
        kres=[]
        Mv=self.dataV
        for i in Mv:
            kres.append(factor(i))
        self.dataV=kres
        self.sort_val()
        self.show_data()
        
    def texpand(self):
        kres=[]
        Mv=self.dataV
        for i in Mv:
            kres.append(expand_trig(i))
        self.dataV=kres
        self.sort_val()
        self.show_data()  
        
    def tsimplify(self):
        kres=[]
        Mv=self.dataV
        for i in Mv:
            kres.append(trigsimp(i))
        self.dataV=kres
        self.sort_val()
        self.show_data()
        
    def opemat(self,kope=''):
        kres=[]
        Mv=self.dataV
        for i in Mv:
            kres.append(opemat(i,kope=kope))
        self.dataV=kres
        self.sort_val()
        self.show_data()    
def categoria(ksym,kval): # used by Is_dependient
    kres=ksym
    for i in kval:
        kres=kres.subs(i,1)
    try:
        kres=opemat(kres,'v')
        return kres
    except:
        return kres
        
def Is_dependient(ksym,kval):  # Return True if value contain data from sample 
    mm=categoria(ksym,kval)
    return Is_Number(mm)        
    
def MsetValue(vEq=[],op1='',op2=''):
    mm=[]
    for i in vEq:
        eqq=i
        eqq.setValue(op1,op2)
        mm.append(eqq)    