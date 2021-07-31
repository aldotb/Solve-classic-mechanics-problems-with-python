

import sympy
import lib_algebraEq as ae 
from libaldo_math import * 
from libaldo_show import *
import copy
import pandas as pd
    
    
class polyclass:
    def __init__(self,p1='',p2=''):
        if p2=='':
            self.q1=p1.args[0]
            self.q2=p1.args[1]
            self.Q=ae.Equation(self.q1,self.q2)
            self.q3=p1.args[0]-p1.args[1]
            
             
        else:    
            self.q1=p1
            self.q2=p2
            self.Q=ae.Equation(p1,p2)
            self.q3=p1-p2
        self.H=self.Q
        self.V=[]
        
    def showq(self):
        return(self.Q)
    
    
    def re_eQ(self,p1,p2):
        kres=ae.Equation(p1,p2)
        self.update(kres)
        return(kres)
    
    def update(self,newQ):
        self.H=self.Q
        self.Q=newQ
        self.q1=newQ.args[0]
        self.q2=newQ.args[1]
        self.q3=newQ.args[0]-newQ.args[1]
    
    def updateQ(self):
        self.Q=ae.Equation(self.q1,self.q2)
        self.q3=self.q1-self.q2 
        return(self.Q)
        
    def opemat(self,kope='',op1='',op2=''):
        if op1=='1':
            kres=ae.Equation(opemat(self.q1,kope),self.q2)
        elif op1=='2':
            kres=ae.Equation(self.q1,opemat(self.q2,kope))
        else:
            kres=self.Q
            kres=opemat(kres,kope)
        if op1==1 or op2==1:
            self.update(kres)
        return(kres)
    
        
    def ope4( self,op1='',op2=0,kd='0',kope=''):
        kres=self.Q
        if op1 in 'SMDP':
            if op1=='S':
                kres=kres+op2
            if op1=='M':
                kres=kres*op2     
            if op1=='D':
                kres=kres/op2
            
            if op1=='L':
                kres=kres/op2
            
            if op1=='P':
                kres=pow(kres,op2)
            try:
                kres=opemat(kres,kd)
            except :
                kres=opemat(kres,kope)
                 
            if kd==1:
                self.update(kres)        
                 
        return(kres)
        
    def ope1( self,op1='',op2=0,kd='0',kope=''):
        kres=self.q1
        if op1 in 'SMDP':
            if op1=='S':
                kres=kres+op2
            if op1=='M':
                kres=kres*op2     
            if op1=='D':
                kres=kres/op2
            if op1=='P':
                kres=pow(kres,op2)
            try:
                kres=opemat(kres,kd)
            except :
                kres=opemat(kres,kope)
                 
            if kd==1:
                newEq=ae.Equation(kres,self.q2)
                self.update(newEq) 
                return(self.Q)                
                #return(kres)
        return(ae.Equation(kres,self.q2))
                
    def ope2( self,op1='',op2=0,kd='0',kope=''):
        kres=self.q2
        if op1 in 'SMDP':
            if op1=='S':
                kres=kres+op2
            if op1=='M':
                kres=kres*op2     
            if op1=='D':
                kres=kres/op2
            if op1=='P':
                kres=pow(kres,op2)
            try:
                kres=opemat(kres,kd)
            except :
                kres=opemat(kres,kope)
                 
            if kd==1:
                newEq=ae.Equation(self.q1,kres)
                self.update(newEq) 
                return(self.Q)                
                #return(kres)
        return(ae.Equation(self.q1,kres))
        
    def swap(self):
        p2=self.q1
        p1=self.q2
        kres=ae.Equation(p1,p2)
        self.Q=kres
        return(kres)
    
    def solve(self,op1,kd=''):
        kres=self.q1-self.q2
        if kd!='':
            kres=csolve(kres,op1,kd)
        else:
            kres=csolve(kres,op1)
        return(kres)
    def get_exp(self,op1,kope=''):
        
        kres=self.Q
        kres=kres.args[int(op1)-1]
        if type(kres)==Pow:
            kres=kres.args[1]
        kres=opemat(kres,kope)    
        return(kres)
    
    def simpPow(self,op1):
        kres1=self.q1
        kres2=symbols(self.q2)
        if op1==1:
            kres1=poly(kres1,'simpow')
        if op1==2:
            kres=type(kres2)  
        #kres=ae.Equation(kres1,kres2) 
        return(kres)    
    
    def subs(self,op1,op2='',kd=''):         
        kres=self.Q
        if op2=='' or op2=='up':
            kres1=(self.q1).subs(op1)
            kres2=(self.q2).subs(op1)
            kres=ae.Equation(kres1,kres2)
            if op2=='up':
                kd=1
             
        else:
            kres=kres.subs(op1,op2)

        if kd!='':
                self.update(kres)
        return(kres)        
        
    def subsubs(self,op1,op2,kd=''):
        kres=self.Q
        for k1,k2 in zip(op1,op2):
            kres=kres.subs(k1,k2)
        if kd!='':
                self.update(kres)
        return(kres)
    
    def ksubs(self,op1,op2,kd=''):
        kres=self.Q
        try:
            kres=kres.subs(op1,op2) 
        except:
            kres=self.Q
        if kd!='':
                self.update(kres)
        return(kres)
        
    def evalue_if(self,op1,op2,kope=''):
        kres1=self.q1
        kres1=kres1.subs(op1,op2)
        kres1=opemat(kres1,kope)
        
        kres2=self.q2
        kres2=kres2.subs(op1,op2)
        kres2=opemat(kres2,kope)
        
        kres=(kres1==kres2)
        print(kres)
        return(ae.Equation(kres1,kres2))
        
    def simpFac(self,kd=''):
        m1,m2=self.q1,self.q2
        kfac= gcd(m1,m2)
        kres1=simplify(m1/kfac)
        kres2=simplify(m2/kfac)
        kres=ae.Equation(kres1,kres2)
        if kd!='':
                self.update(kres)
        return(kres)
    
    def simp_log_LE(self,kd=''):
        kres=ae.Equation(elel(self.q1),elel(self.q2))
        if kd!='':
                self.update(kres)
        return(kres)
        
    def simp_log_S2M(self,kd=''):
        kres=ae.Equation(opelog2(self.q1,'f'),opelog2(self.q2,'f'))
        if kd!='':
                self.update(kres)
        return(kres)
        
    def dothis(self,ksym1,ksym2,op1=0):
        ksym1=ksym2
        print(ksym1)
        
        
    def polymath(self,op1='',op2=0):
        if op1=='':
            monomath()
        else:    
            done1=True
            done2=True
            kres=self.Q
            kres1=self.q1
            kres2=self.q2
            if '1' not in op1:
                
                kres1=monomath(kres1,op1)
            if '2' not in op1:
                kres2=monomath(kres2,op1)
                        
            kres=ae.Equation(kres1,kres2)
            if op2==1:
                    self.update(kres)
            return(kres)
        
    def make_grp(self,ksym,op2='',side=0):
        kres1=make_grp(self.q1,ksym)
        kres2=make_grp(self.q2,ksym)
        if side==1:
            kres2= self.q2 
        elif side==2:
            kres1= self.q1
        kres=ae.Equation(kres1,kres2)
        if op2==1:
            self.update(kres)
        return(kres)
        
    def make_grp_sec (self,kvec,op2='',side=0):
        kresq=self
        for i in kvec:
            kres=kresq.make_grp(i,op2=op2,side=side)
            kresq=polyclass(kres)
            #print(kres)
            #self.Q=ae.Equation(kres.q1,kres.q2)
        if op2==1:
             
            self.Q=kresq.Q
            self.q1=kresq.q1
            self.q2=kresq.q2
        return(kres)
    
    def multitask(self,vec,kd=''):
        kres1=self.q1
        kres2=self.q2
        op1=vec[0]
        if vec[0]==1:
            kres1=multitask(kres1,vec[1::])
        elif vec[0]==2:    
            kres2=multitask(kres2,vec[1::])
        else:
            kres1=multitask(kres1,vec)
            kres2=multitask(kres2,vec)
        
        kres=ae.Equation(kres1,kres2)
        if kd!='':
                self.update(kres)
        return(kres)    
            
        
    def free(self,op=''):
        kres1=self.q1
        kres2=self.q2
        try:
            mm=fpoly(kres1,'free')
            if op==1:
                return(mm)
        except:
            mm=[]
            if op==1:
                return(mm)
        try:    
            mm2=fpoly(kres2,'free')
            if op==2:
                return(mm2)
            
            for i in mm2:
                if i not in mm:
                    mm.append(i)
        except:
            if op=='N':
                return(len(mm))
            else:    
                return(mm)
        if op=='N':
                return(len(mm))
        else:    
            return(mm)
#  ********  Unific variable ****

    def unifique_exp(self):
        self.q1=parse_expr(str(self.q1))
        self.q2=parse_expr(str(self.q2))
        self.updateQ()
        
class Mpolyclass:
    def __init__(self,eQ):
        self.Q=eQ
        self.mQ= [[i.q1,i.q2] for i in self.Q]
        self.q1=[i.q1 for i in self.Q]
        self.q2=[i.q2 for i in self.Q]
        self.q3=[i.q3 for i in self.Q]
        self.w1=[poly_weight(i) for i in self.q1]
        self.w2=[poly_weight(i) for i in self.q2]
        self.qq=len(self.mQ)
        mm=[]
        mwt=[]
        for i in range (len(self.mQ)):
            kres=self.q1[i]-self.q2[i]
            mm.append(kres)
            mwt.append(poly_weight(kres))
            
        self.wt=mwt
            
        
        self.mtemp=0
         
    def getmono(self,x,y):
        kres=self.mQ[x][y]
        return(kres)

    def getEq(self,op1):
        kres=self.Q[op1]
        return kres.Q
        
    def getPoly(self,op1):
        kres=self.Q[op1]
        return kres   
        
    def redf_Q(self):
        self.Q=[polyclass(i[0],i[1]) for i in self.mQ]
        self.qq=len(self.mQ)
        self.ee=[opemat(i.q1-i.q2,'s') for i in self.Q]
        
    
    def redf_mQ(self):
        self.mQ= [[i.q1,i.q2] for i in self.Q]
        self.ee=[opemat(i.q1-i.q2,'s') for i in self.Q]
        self.qq=len(self.mQ)
        
    def add_Q(self,ksym):
        kres=self.Q 
        kres.append(ksym)
        self.Q=kres
        self.redf_mQ()
    
    def add_mQ(self,op1,op2=0):
        kres=self.mQ 
        kres.append([op1,op2])
        self.mQ=kres
        self.redf_Q()    
        
    def swap(self,op):
        kres1=self.q1[op]
        self.q1[op]=self.q2[op]
        self.q2[op]=kres1
        self.up_class()
         
    
    def up_Df(self,kdf):
        k1=list(kdf["E1"].values)
        k2=list(kdf["E2"].values)
        mm=[]
        for i in range(self.qq):
            self.mQ[i]=(k1[i],k2[i])
        self.redf_Q()
        
    def up_class(self): # con respecto a los valores de q1 y q2
        mq1=[]
        mq2=[]
        mQ=[]
        for i in range(self.qq):
            if not(self.q1[i]==[] or self.q2[i]==[]):
                mq1.append(self.q1[i])
                mq2.append(self.q2[i])
                mQ.append(ae.Equation(self.q1[i],self.q2[i]))
        
        self.q1=mq1
        self.q2=mq2
        self.Q=mQ
        self.qq=len(mq1)
        qq1=self.q1
        qq2=self.q2
        qq=self.qq
        self.w1=[poly_weight(i) for i in qq1]
        self.w2=[poly_weight(i) for i in qq2]
           
        mm=[]
         
        mwt=[]
        for i in range (qq):
            kres=self.q1[i]-self.q2[i]
            mm.append(kres)
            mwt.append(poly_weight(kres))
        self.q3=mm    
        self.wt=mwt
        
        
    def up_from_df(self,df):
        self.q1=df['E1'].tolist()
        self.q2=df['E2'].tolist()
        self.up_class()
        
         
    def up_follow(self,vec):
        kq1=[]
        kq2=[]
        for i in vec:
            kq1.append(self.q1[i])
            kq2.append(self.q2[i])
        self.q1=kq1
        self.q2=kq2
        self.up_class()
        return(self.make_IA2())
        
    def del_eQ(self,op):
        mm=[]
        for i in range(self.qq):
            if i!=op:
                mm.append(self.mQ[i])
        self.mQ=mm
        self.redf_Q()    
                
#  *** Show Mpolyclass *******
   
    def showQ(self):
        for i in range(self.qq):
            sE([self.q1[i],'=',self.q2[i] ])
    
    def s(self): # shoet name of showQ
        return(self.showQ()) 
            
    def showE(self,kope=''):
        for i in self.Qp:
            kres=opemat(i,kope)
            sE([kres])
            
  
    def make_IA(self):
        vecpu=self.iMat('W')
        vecpt=self.iMat('W',True)
        vecst=self.iMat('free',True)
        meq1,meq2,mpu1,mpu2,mpt,msym=[],[],[],[],[],[]
        for i in range (self.qq):
            meq1.append(self.mQ[i][0])
            meq2.append(self.mQ[i][1])
            mpu1.append(vecpu[i][0])
            mpu2.append(vecpu[i][1])
            mpt.append(vecpt[i])
            msym.append(vecst[i])
        
        data={'E1':meq1,
              'W1':mpu1,
              'E2':meq2,
              'W2':mpu2,
              'Wt':mpt,
              'Ksym':msym
              }    
        df=pd.DataFrame(data)
        return(df)
        
    def make_IA2(self):
        msym= []
        vecst=self.iMat('free')
        for i in range (self.qq):
            msym.append(vecst[i])    
            
        data={'E1':self.q1,
              'W1':self.w1,
              'E2':self.q2,
              'W2':self.w2,
              'Wt':self.wt,
              'Ksym':msym
              }    
        df=pd.DataFrame(data)
        return(df)
                
    def make_IAS(self):
        data={'S1':self.iMat('free',1),
              'W1':self.iMat('W',1),
              'N1':self.iMat('n',1),
              'S2':self.iMat('free',2),
              'W2':self.iMat('W',2),
              'N2':self.iMat('n',2),
              'St':self.iMat('free'),
              'Wt':self.iMat('W'),
              'Nt':self.iMat('n')
               }
        df=pd.DataFrame(data)
        return(df)       
                  
    def  get_data_pd(self):
        df=self.make_IA2()
        kq1=df['E1'].tolist()
        return(kq1)
        
        
    #********* Information*******
    
    def get_freelist(self,op1='',op2=''):
        if op1=='':
            mm=[]
            for i in range(self.qq):
                kres=self.Q[i]
                for j in kres.free():
                    if j not in mm:
                        mm.append(j)
            return mm             
        else:                    
            kres=self.getPoly(op1)
            return(kres.free(op=op2))
     
    def mat_EqLin(self,op=''):
        cc=0
        mm=[]
        for i in self.mQ:
            kres1=i[0]
            kres2=i[1]
            if op==1:
                mm.append([cc,int(type(kres1)==Symbol),int(type(kres2)==Symbol)])
            else:
                mm.append([int(type(kres1)==Symbol),int(type(kres2)==Symbol)])
            cc+=1
        return(mm)
    
	
	
    def re_LinForm(self):
        mm=	self.mat_EqLin(1)
        for i in mm:
            kstr=str(i[1])+str(i[2])
            if kstr=='01':
                nEq=i[0]
                self.swap(nEq)
    
    def is_solu(self,op):
        kres1=self.q1[op]
        kres2=self.q2[op]
        if type(kres1)==Symbol and (type(kres2)==Integer or type(kres2)==float or type(kres2)==int or type(kres2)==Rational):
            return True
        return False
    
    def full_solve(self):
        qq=self.qq
        done=True
        for i in range(qq):
            if not self.is_solu(i):
                done=False
        return done
        
#********* algoritmo *******
    def opemat(self,kope='',vec=[]):
        qq=self.qq
        for i in range(qq):
            kres1=self.q1[i]
            kres2=self.q2[i]
            if i in vec or vec==[]: 
                    kres1=opemat(kres1,kope)
                    kres2=opemat(kres2,kope)
           
            self.q1[i]=kres1
            self.q2[i]=kres2    
        self.up_class()
        
    def ope4( self,op1='',op2='',vec=[]):
        for i in range(self.qq):
            if i in vec or vec==[]:
                kres1=self.q1[i]
                kres2=self.q2[i]
            
                if op1 in 'SMDP':
                    if op1=='S':
                        kres1=kres1+op2
                        kres2=kres2+op2
                    if op1=='M':
                        kres1=kres1*op2
                        kres2=kres2*op2
                    if op1=='D':
                        kres1=kres1/op2
                        kres2=kres2/op2
                    
                    if op1=='P':
                        kres1=pow(kres1,op2)
                        kres2=pow(kres2,op2)
                    
                         
                    self.q1[i]=kres1      
                    self.q2[i]=kres2     
                    self.up_class() 
           
        self.s()
            
    def multitask(self,keq=[],vec=[],kshow=False):
        qq=self.qq
        for i in range(qq):
            if i in keq or keq==[]:
                op=i 
                kres1=self.q1[op]
                kres2=self.q2[op]
                op1=vec[0]
                if vec[0]==1:
                    kres1=multitask(kres1,vec[1::])
                elif vec[0]==2:    
                    kres2=multitask(kres2,vec[1::])
                else:


                    kres1=multitask(kres1,vec)
                    kres2=multitask(kres2,vec)
                self.q1[op]=kres1
                self.q2[op]=kres2
                self.up_class()
        if kshow:
            return(self.showQ())
            
    def multitask2(self,op='',vec=[]):
        if type(op)==int and vec!=[]:
            kP=polyclass(self.q1[op],self.q2[op])
            kres=kP.multitask(vec=vec)
            self.q1[op]=kres.args[0]
            self.q2[op]=kres.args[1]
            self.up_class()
                   
                 
    def null_off(self):
        mm=self.ee
        kres=[]
        for i in range(self.qq):
            if opemat(mm[i],'s')!=0:
                kres.apend(self.Q[i])
        self.Q=kres
        self.redf_mQ()
        
                
            
    def subsubs(self,ksym,kval):
        qq=self.qq
        for i in range(qq):
            kres=self.q1[i]
            self.q1[i]=kres.subs(ksym,kval)
            kres=self.q2[i]
            self.q2[i]=kres.subs(ksym,kval)
        self.up_class()
        return(self.s())
        
    def info_mat(self,op='',op2=''):
        mm=[]

        if op2=='':
            if op=='W':        
                for i in range(self.qq):
                    mm.append(poly_weight(self.q1[i]-self.q2[i]))
                return(mm)    
            elif op=='free':
                for i in range(self.qq):
                    mm.append(fpoly(self.q1[i]-self.q2[i],'free'))
                return(mm)
            elif op=='n':
                for i in range(self.qq):
                    mm.append(fpoly(self.q1[i]-self.q2[i],'n'))
                return(mm) 
            else:
                return 0
        else:
            mm=[]
            if op2==1:
                kqq=self.q1
            else :
                kqq=self.q2
            if op=='W':
                mm=[poly_weight(x) for x in kqq]
            elif op=='free':
                mm=[fpoly(x,'free') for x in kqq]
            else:
                mm=[fpoly(x,'n') for x in kqq]
                
            return(mm)    
           
      
    def iMat(self,op='',op2=''):
        kres=self.info_mat(op=op,op2=op2)
         
        return(kres)

    def change_eq(self,ksym,op,kd=''):
        kres=csolve(self.q3[op],ksym)
        self.q1[op]=ksym
        self.q2[op]=kres
        self.up_class() 
        df=self.make_IA2()
         
        return(df)
            
 
    
    def sortEq(self,vsort):
        df=self.make_IA2()
        df2=df.sort_values(by=vsort)
        self.up_from_df(df2)
        return(df2)
        
    def sort_w(self):
        self.sortEq(['W1','W2'])
    
    def solve2Lin(self,vecs,vecp):
        ksym1,ksym2=vecs
        p1,p2=vecp
        kval1,kval2=self.q2[p1],self.q2[p1]
        eq1=kval1.subs(ksym2,kval2)
        kres1=csolve(ksym1-eq1,ksym1)
        eq2=kval2.subs(ksym1,kres1)
        kres2=csolve(ksym2-eq2,ksym2)
        self.q2[p1],self.q2[p1]=kres1,kres2
        self.up_class() 
        df=self.make_IA2()
        return(df)
    
    # algoritmos de autosolve
    
    def solve_from(self,nEq,ksym,kd=''):
        
        kEq=self.q1[nEq]-self.q2[nEq]
        try:
            if kd=='':
                return(csolve(kEq,ksym))
            else:
                return(csolve(kEq,ksym,kd))
                
      
        except:
            return(self.Q[nEq].Q)
        
    def csolves(self,op='',kreturn=False,kope='',kdisp=False):
        kEq=[]
        ksym=[]
        ksum=0
        for i in range(self.qq):
            Eqq=Eq(self.q1[i],self.q2[i])
            kEq.append(Eqq)
            ksum+=self.q1[i]+self.q2[i]
        if op!='':
            k2=op
        else :       
            ksym=fpoly(ksum,'free')
            k2=tuple(ksym)
        k1= tuple(kEq) 
         
            
        kres=[]
        try:
            kres=solve(k1,k2)
            if kreturn:
                mm=[]
                for i in list(kres.items()):
                    mm.append(opemat(i[1],kope))
                if kdisp:
                    cc=0
                    for i in mm:
                        sE([str(op[cc]),'=',i])
                        cc+=1
                        
                return(mm)
            else:    
                return(kres)
        except:    
            return(k1,k2)
    
    def fix_bando(self,op=''):
        m1=[]
        m2=[]
        qq=self.qq
        qq1=self.q1
        qq2=self.q2
        for i in range(qq):
            e1=qq1[i]
            e2=qq2[i]
            if len(fpoly(e2,'free'))==1:
                m2.append(e1)
                m1.append(e2)
            else:
                m1.append(e1)
                m2.append(e2)
             
        self.q1=m1
        self.q2=m2
        self.up_class()
        if op!='':
            return(self.s())
    
    def fix_alone(self,op=''):
        for i in range(self.qq):
            self.q1[i],self.q2[i]=self.redu_p1(i)
        self.up_class()
        if op!='':
            return(self.s())
    
    def fix_free(self,op=''):
        for i in range(self.qq):
            if len(fpoly(self.q1[i],'free'))==0:
                self.swap(i)
        self.up_class()
        if op!='':
            return(self.s())
        
    def fix_first(self,op=0):
        self.fix_free()
        self.fix_bando()
        self.fix_alone()
        self.apli_simp_solu(op=op)
        return(self.s())
        
        
    
    def redu_p1(self,op=0):
    
        ksym1=self.q1[op]
        ksym2=self.q2[op]
        n1=fpoly(ksym1,'n')
        q1=len(fpoly(ksym1,'free'))
        if q1==1:
            ksym=fpoly(ksym1,'free')[0]
            kres=csolve(ksym2-ksym1,ksym)
            rs1=ksym
            rs2=kres
            return(rs1,rs2)    
        else:
            return(ksym1,ksym2)
            
    def redu_p2(self,op=0):
    
        ksym1=self.q1[op]
        ksym2=self.q2[op]
        n1=fpoly(ksym1,'n')
        q1=len(fpoly(ksym1,'free'))
        if q1==1:
            ksym=fpoly(ksym1,'free')[0]
            kres=csolve(ksym2-ksym1,ksym)
            rs1=ksym
            rs2=kres
            return(rs1,rs2)    
        else:
            return(ksym1,ksym2) 
            
    def simplify_idem(self,op=''): # busca terminos iguales para eliminarlos
        for i in range(self.qq):
            kres1=self.q1[i]
            kres2=self.q2[i]
            try:
                kvec1=fpoly(kres1,'list')
                kvec2=fpoly(kres2,'list')
                for j in kvec1:
                    if j in kvec2:
                        kres1=kres1-j
                        kres2=kres2-j
                        
                self.q1[i]=kres1
                self.q2[i]=kres2
            except:
                done=1
        self.up_class()
        if op!='':
            return(self.showQ())
            
    def redu_makein(self,ksym,knum,kall=True):
        kval=self.q2[knum]
        for i in range(self.qq):
            if i!=knum:
                Eq1=self.q2[i]
                Eq2=Eq1.subs(ksym,kval)
                self.q2[i]=Eq2
        self.up_class()
    
    def first_sol(self):
        for i in range(self.qq):
            Eq1=self.q3[i]
            if  len(fpoly(Eq1,'free'))==1:
                ksym=fpoly(Eq1,'free')[0]
                kres=csolve(Eq1,ksym)
                self.q1[i]=ksym
                self.q2[i]=kres
                self.makein(ksym,i,kall=False)
        self.up_class()
        
    def apli_simp_solu(self,op=0,kshow=True):
        ksym=self.q1[op]
        kval=self.q2[op]
        for i in range(self.qq):
            if i!=op:
                kres1=self.q1[i]
                kres2=self.q2[i]
                try:
                    kres1=kres1.subs(ksym,kval)
                except:
                        done=False    
                try:
                    kres2=kres2.subs(ksym,kval)
                except:
                        done=False
                self.q1[i]=kres1
                self.q2[i]=kres2
        self.up_class()
        if kshow:
            self.s()
    
    def appli_solu1(self):
        qq=self.qq
        for i in range(qq):
            if self.is_solu(i):
                self.apli_simp_solu(i)
                self.up_class()
                return(True)
                
                
            
    def change_rep(self):
        qq=self.qq
        for i in range(qq):
            
            ksym=self.q1[i]
            kval=self.q2[i]
            for j in range(i+1,qq):
                if self.q1[i]==self.q1[j]:
                    self.q1[j]=self.q2[i]
                    self.up_class()
                    return(True)
                        
                    
    def reduce3(self):
        qq=self.qq
        for i in range(qq-1):
            if self.is_solu(i):
                ksym=self.q1[i]
                kval=self.q2[i]
                for j in range(i+1,qq):
                    kres1=self.q1[j]
                    kres2=self.q2[j]
                    try:
                        kres2=kres2.subs(ksym,kval)
                    except:
                        done=False
                    self.q2[j]=kres2
                    if not self.is_solu(j):
                        kres1=kres1.subs(ksym,kval)
                        self.q1[j]=kres1
                
                self.up_class()    
                     
    def pro_N_2right(self):
        qq=self.qq
        for i in range(qq):
            kres1=self.q1[i]    
            if  (type(kres1)==Integer or type(kres1)==float or type(kres1)==int or type(kres1)==Rational):
                kres2=self.q2[i]
                ksymL=fpoly(kres2,'free')
                ksym=ksymL[0]
                kval=csolve(kres2-kres1,ksym)
                self.q1[i]=ksym
                self.q2[i]=kval
                self.up_class()
            
        
    def subsolve(self,ksym,kitem,kremp=False): # solve ksym in eqNum kitem and remplaze
        kres1=self.q1[kitem]
        kres2=self.q2[kitem]
        kres=kres1-kres2
        nval=csolve(kres,ksym)
        self.q1[kitem]=ksym
        self.q2[kitem]=nval
        self.up_class()
        if kremp:
            self.apli_simp_solu(op=kitem)
        else:
            return(self.s())
            
    def subsolve2(self,ksym,kitem,kremp=True,kshow=True): # solve ksym in eqNum kitem and remplaze
        kres1=self.q1[kitem]
        kres2=self.q2[kitem]
        kres=kres1-kres2
        nval=csolve(kres,ksym)
        self.q1[kitem]=ksym
        self.q2[kitem]=nval
        self.up_class()
        if kremp:
            self.apli_simp_solu(op=kitem,kshow=False)
        if kshow:
            self.showQ()
    
    def try_solve(self,vksym,vkitem,kremp=True,tkshow=True):
        for ksym,kitem in zip(vksym,vkitem):
            self.subsolve2(ksym,kitem,kremp=kremp,kshow=False)
        if tkshow:
            self.showQ()
            
            
        
#****  end Class  ******
    def matSolve(self,ksym,kor,kremp=True,kshow=False):
        qq=len(ksym)
        for i in range(qq):
            self.subsolve2(ksym[i],kor[i],kremp=kremp,kshow=kshow)
        return(self.s())    
              
def is_eQ_lin(psym):
    done=False
    kres1=psym.q1
    kres2=psym.q2
    if type(kres1)==Symbol and  len(fpoly(kres2,'free'))==1:
        done=True
    return(done)    
        
               
def show_mat_eq(kEq):
    for i in kEq:
        sE([i.Q])            

def join_symbols(ksym):
    kres=ksym
    vsym=fpoly(ksym,'free')
    vsyms=fpoly(ksym,'frees')
    for i,j in zip(vsym,vsyms):
        kres=kres.subs(parse_expr(j),i)
    return(kres)
        
def makemono(ksym,op):
    kres=ksym
    if op[0]=='P':
        if op[1]=='(':
            newsym=parse_expr(op[2:-1])
            kres=kpow(kres,newsym)
            kres=join_symbols(kres)
            #kres=opemat(kres,'es')
        else:    
            kexp=int(op[1::])
            kres=kpow(kres,kexp)
            #kres=opemat(kres,'es')
        return(kres)
        
    elif op[0]=='R':
        if op[1]=='(':
            newsym=parse_expr(op[2:-1])
            kres=rpow(kres,newsym)
            kres=join_symbols(kres)
            #kres=opemat(kres,'es')
        else:    
            kexp=int(op[1::])
            kres=rpow(kres,kexp)
            #kres=opemat(kres,'es')
        return(kres)    
    
    elif op[0]=='S':
        if op[1]=='(':
            newsym=parse_expr(op[2:-1])
            kres=kres+newsym 
            kres=join_symbols(kres)
            #kres=opemat(kres,'es')
        else:    
            kexp=int(op[1::])
            kres=kres+kexp
            #kres=opemat(kres,'es')
        return(kres)
        
    elif op[0]=='M':
        if op[1]=='(':
            newsym=parse_expr(op[2:-1])
            kres=kres*newsym 
            kres=join_symbols(kres)
            #kres=opemat(kres,'es')
        else:    
            kexp=int(op[1::])
            kres=kres*kexp
            #kres=opemat(kres,'es')
        return(kres)
    
    elif op[0]=='F':
        if op[1]=='(':
            newsym=parse_expr(op[2:-1])
            kres=kres/newsym 
            kres=join_symbols(kres)
            #kres=opemat(kres,'es')
        else:    
            kexp=int(op[1::])
            kres=kres/kexp
            #kres=opemat(kres,'es')
        return(kres)    
    
    elif op[0]=='N' or op[0]=='D':
         
        newop=op[1::]
        knum=numer(kres)
        kdeno=denom(kres)
        if op[0]=='N':
            knum=makemono(knum,newop)
            
                    
        elif op[0]=='D':
            kdeno=makemono(kdeno,newop)
            
        kres=knum/kdeno
        kres=join_symbols(kres)
        #kres=opemat(kres,'es')         
        return(kres)
    
    elif call_opemat(op):
        kres=opemat(kres,op)
        return(kres)
    
    else:
        return(ksym)
        
def multitask(ksym,vec):
    kres=ksym
    for i in vec:
        kres=makemono(kres,i)
    return(kres)        
    
#****  Funcion sympy equation args composition     ******

def call_opemat(ksym):
    vec='sfetx'
    kres=True
    for i in ksym:
        if not i in vec:
            kres=False
    return(kres)

def short_tree_exp(ksym): # 'STE(ksym)' short called function 
    kres=srepr(ksym)
     
    nomL=['Symbol','Integer','Add','Mul','Pow','Rational','Float','.','1','2','3','4','5','6','7','8','9','00','000','0000','00000','000000','0000000','00000000','000000000']
    nomC=['S','I','A','M','P','R','F','','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0']
    for i,j in zip(nomL,nomC):
        kres=kres.replace(i,j)
    kres=kres.replace(' ','')
    nomL=['00','000','0000','00000','000000','0000000','00000000','000000000']
    nomC=['0','0','0','0','0','0','0','0']
    for i,j in zip(nomL,nomC):
        kres=kres.replace(i,j)
    kres=kres.replace('00','0')
    kres=kres.replace(',precision=0','')
    kres=kres.replace(',positive=True','')
    kres=kres.replace('(0)','')
    kres=kres.replace(",positive=True", "")
     
    kres=kres.replace("'","")
    return(kres)
 
def short_type_name(ksym): # 'TT(ksym)' short called function 
    Nsym=type(ksym)
    if Nsym == Symbol :
        return('S')
    elif Nsym == Integer :
        return('I')
    elif Nsym == Add :
        return('A')
    elif Nsym == Mul :
        return('M')
    elif Nsym == Pow :
        return('P')
    elif Nsym == Rational :
        return('R')
    elif Nsym == Float :
        return('F')
    else:
        return('N')    

def poly_weight(ksym):
    kres=STE(ksym)
    ktot=0
    
    ktot+=0*kres.count('I)')
    ktot+=1*kres.count('S(')
    ktot+=4*kres.count('A(')
    ktot+=3*kres.count('M(')
    ktot+=8*kres.count('P(')
    ktot+=10*kres.count('R(')
    ktot+=3*kres.count('F(')
    return ktot
   
def TT(ksym):
    return(short_type_name(ksym))
def STE(ksym):    
    return(short_tree_exp(ksym))
    
#  Algoritmos de autosolucion




    
    
     