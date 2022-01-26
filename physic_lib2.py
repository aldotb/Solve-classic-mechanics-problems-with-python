
from  libaldo_math2 import *
from  lib_MyEq import *
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import Image, display
 
 

from sympy import symbols
m,m1,m2,m3,m4,M,g,x,x1,x2,y,y1,y2,X,Y,a,a1,a2,a3,v,v1,v2,M1,M2,M3= symbols('m m1 m2 m3 m4 M g x x1 x2 y y1 y2 X Y a a1 a2 a3 v v1 v2 M1 M2 M3')
w,w1,w2,aw,aw1,aw2,F,F1,F2,Rx,Ry,r,r1,r2,R,ax,ay,V= symbols('w w1 w2 aw aw1 aw2 F F1 F2 R_x R_y r r1 r2 R ax ay V')
mu,mu1,mu2,fr,fr1,fr2,N1,N2,Nm, L,L1,L2,h,h1,h2,b,H= symbols('mu  mu1 mu2  fr fr1 fr2 N1 N2 Nm  L L1 L2 h h1 h2 b H')
alpha,tetha,ac,at,alpha1,alpha2=symbols('alpha theta ac at alpha1 alpha2')
t,vx,vy,vxy,Po,Ti,I_n,In,a_w,Fc,a_c,a_t,  Io=symbols('t vx vy vxy Po Ti I_n In a_w Fc a_c a_t    I_o')
T,T1,T2,t1,t2,t3=symbols('T T1 T2 t1 t2 t3')
K,X1,X2,X0,d,W,P=symbols('K X1 X2 X0 d W P')
rho,z,z1,z2,A,p=symbols('rho z  z1 z2 A p') 
xo,yo,zo,beta,xi,xf=symbols('xo yo zo beta x_i x_f')
yp,xp,pp=symbols("y' x' p'") 
# diff variables

dm,ds,dx,dy,dz,dt,dr,dh,dL,da,dA,dv,dV=symbols('dm ds dx dy dz dt dr dh dL da dA dv dV')   
 
e1=MyEq(0,kname='e1',kshow=False,ktype='C',init=False)
e2=MyEq(0,kname='e2',kshow=False,ktype='C',init=False)  
e3=MyEq(0,kname='e3',kshow=False,ktype='C',init=False)  
e4=MyEq(0,kname='e4',kshow=False,ktype='C',init=False)  
e5=MyEq(0,kname='e5',kshow=False,ktype='C',init=False)   
e6=MyEq(0,kname='e6',kshow=False,ktype='C',init=False)   
e7=MyEq(0,kname='e7',kshow=False,ktype='C',init=False) 
e8=MyEq(0,kname='e8',kshow=False,ktype='C',init=False)
e9=MyEq(0,kname='e9',kshow=False,ktype='C',init=False)
class mparticle:
    def __init__(self,x1=x1,x2=x2,y1=y1,y2=y2,v1=v1,v2=v2,m=m,a=a,g=g,v=v,w=w,ac=ac,s='r',t=t,r=r,r1=r1,r2=r2,vx=vx,vy=vy,vxy=vxy,Ti=Ti,I_n=In,In=In,aw=aw,w1=w1,w2=w2,a_c=a_c,a_t=a_t,ax=ax,ay=ay,typeI='',mu=mu, Nm=Nm,Itype='p',xI=0):
        
        self.m=m   # mass      
        self.g=g   # gravity   
        self.x1=get_real_value(x1) # x1      start Point >(x1,y1),v1D:\Libaldo\SolveLib\Fin_F\physic_lib.py
        self.y1=get_real_value(y1) # y1 
        self.v1=get_real_value(v1) # v1
        self.x2=get_real_value(x2) # x2      finish  Point >(x2,y12),v2
        self.y2=get_real_value(y2) # y2
        self.v2=get_real_value(v2) # v2
        self.v=get_real_value(v)   # Initial Velocity ,  important in parabolic 
        self.a=sex2rad_i(a,s) # in parabolic is the shooting angle
        self.ac=get_real_value(ac) # ac=X acceleration, set ac=0, in parabolic case
        self.r=get_real_value(r)   # Radio  in circular movements       
        self.r1=get_real_value(r1)
        self.r2=get_real_value(r2)
        self.F=[]  #matriz de fuerzas
        self.cG=[] # Centro de Gravedad
        self.cI=[] # Inerce points
        self.Fc=[] # Fuerzas circulares
        self.kvalue=([],[])
        self.kvaluen=([],[])
        self.t=get_real_value(t)
        self.M=[]
        self.CG=[]
        self.vx=get_real_value(vx)
        self.vy=get_real_value(vy)
        self.vxy=vxy
        self.w=get_real_value(w)
        self.Po=Po
        self.Ti=Ti
        self.I_n=I_n
        self.In=In        
        self.aw=get_real_value(aw)
        self.a_w=aw
        self.a_c=a_c
        self.a_t=a_t
        self.xI=xI 
        self.Ima=[]
        self.w1=get_real_value(w1)
        self.w2=get_real_value(w2)
        self.ax=get_real_value(ax)
        self.ay=get_real_value(ay)
        self.typeI=typeI
        self.mu=mu
        self.Nm=get_real_value(Nm)
        
         


    def add_geoforza(self,kmg,kR,a1,a2=0):  # ********   Manage Variables    **********
        # add_geoforza(Fval,polar ang,angulo whit radio Polar)
        xx,yy= polar2xy(kR,a1,op='')
        self.add_forza(kmg,a1+a2,xx,yy)
        #print(kmg,a2+a3,xx,yy)
    
    def add_uni_forza(self,U=1,x1=0,y1=0,x2=0 ,y2=0):  #  (unitare Value,x1,y1,x2,y2)
        # add_geoforza(Fval,polar ang,angulo whit radio Polar)
        Fx=U*(x2-x1)
        self.add_forza(Fx,0)
        Fy=U*(y2-y1)
        self.add_forza(Fy,pi/2)
        #print(kmg,a2+a3,xx,yy)
    
    def add_forzaP(self,kmg,kang,kR,aRx):
        # add_add_pol_forza(Fvalue,angle whit x,Radio,angle Radio x)
        xx,yy= polar2xy(kR,aRx)
        self.add_forza(kmg,kang,xx,yy)
        #print(kmg,a2+a3,xx,yy)
    
    def addF(self,*args):
        mm=[]
        mmpar=[]
        mmF=self.F
        if type(args[0])!=tuple:
            for i in args:
                mmpar.append(i)
            mm.append(mmpar)
        else:
            for i in args:
                mm.append(i)
        for i in mm:
            if len(i)==2:
                kval,kang=i
                x=0
                y=0
                mmF.append([kval,kang,x,y])
            else:
                kval,kang,x,y=i
                mmF.append([kval,kang,x,y])
                

         
         
                
        self.F=mmF
    def add_forza(self, kval,kang,x=0,y=0,s='r'): 
        # fl=self.F
        # mm=[ x[0]  for x in fl]
        # mm=unisymbols(mm)
        # if unisymbols(kval) not in mm or krep==True:
             
        # if s=='s':
            # kang=sex2rad(kang)
        mm=self.F
         
        mm.append([kval,kang,x,y])        
        self.F=mm
            
            
    def add_weight(self): # is tha same that add_forza(m*g,-pi/2)
        mm=self.m 
        gg=self.g
        ww=mm*gg
        self.add_forza(ww,-pi/2)
        
    def add_Inerce(self,mm='',rr='',L=0,ktype='p'):
        if mm!='':
            self.m=mm
        if rr!='':
            self.r=rr
        if ktype!='p':
            kres=self.setIner(ktype)
            
        kres= self.setIner(ktype)+mm*L*L   
        vecCi=self.cI
        vecCi.append(kres)
        self.cI=vecCi
    
    def addCG(self,*args):
        mm=[]
        mmpar=[]
        mmcG=self.cG
        if type(args[0])!=tuple:
            for i in args:
                mmpar.append(i)
            mm.append(mmpar)
        else:
            for i in args:
                mm.append(i)
        for i in mm:
            kval,x,y=i
            mmcG.append([kval,x,y]) 
        self.cG=mmcG   
        
        
        
        
    def revalue(self,kval):
        mm=unisymbols(self.kvalue)
        mmsym=mm[0]
        mmval=mm[1]
        qq=len(mmsym)
        if qq==0:
            return(kval)
        else:
            kres=unisymbols(kval)
            for i in range(qq):
                nnsym=unisymbols(mmsym[i])
                nnval=unisymbols(mmval[i])
                val=kres/10
                val=val*10
                if type(val)!=float:
                    kres=ksubs(kres,nnsym,nnval)
                    #kres=fpoly(kres,'forze_subs',nnsym,nnval)
            return(kres)
            
    def store_val(self,nomv,valv,kdisp=False,kpositive=False,kope='',ksymbolic=True,kreturn=False,krevalue=False):
        if not self.ksym_in_store(nomv):
            mm=self.kvalue
            qq=len(mm[0])
            if qq>0:
                if not(type(valv)==float or type(valv)==int):
                            valv=self.evalue_val(valv)
                
                
            mat1=[nomv]
            mat2=[valv]
            
            
            
            mm1=mat1+mm[0] 
            mm2=mat2+mm[1] 
            self.kvalue=(mm1,mm2)
            kres=valv
            if krevalue:
                kres=self.revalue(kres)
             
            kres=opemat(kres,kope)
            if kpositive:
                if kres.compare(0)==-1:
                    kres=-1*kres
            if kdisp:
                sR=str(nomv)+' ='        
                display(Math(sR+latex(kres)))
            if kreturn:
                return(kres)
    
    def set(self,nomv,valv): 
        mm=self.F
        qq=len(mm)
        for i in range(qq):
            for ss,vv in zip(nomv,valv):
                if str(ss)=='alpha' or str(ss)=='beta':
                    mm[i][0]=mm[i][0].subs(sin(ss),sin(vv))
                    mm[i][0]=mm[i][0].subs(cos(ss),cos(vv))
                    mm[i][0]=mm[i][0].subs(tan(ss),tan(vv))
                    mm[i][1]=mm[i][1].subs(ss,vv)
                else:
                    mm[i][0]=mm[i][0].subs(ss,vv)
            
        self.F=mm    
         
        
    def setValue(self,nomv,valv,kdisp=False,kpositive=False,kope='',kunisym=True,kreturn=False,krevalue=True,vprimaria=True):
        if not self.ksym_in_store(nomv):
            mm=self.kvalue
            qq=len(mm[0])
            if qq>0:
                if not(type(valv)==float or type(valv)==int):
                            valv=self.evalue_val(valv)
            if vprimaria:    
                if str(nomv)=='m':
                        self.m=valv
                else:
                    self.m.subs(nomv,valv)
                    
                if str(nomv)=='g':
                        self.g=valv
                else:
                    self.g.subs(nomv,valv)        
                if str(nomv)=='x1':
                        self.x1=valv
                if str(nomv)=='x2':
                        self.x2=valv
                if str(nomv)=='y1':
                        self.y1=valv
                if str(nomv)=='y2':
                        self.y2=valv        
                if str(nomv)=='v1':
                        self.v1=valv
                if str(nomv)=='v2':
                        self.v2=valv
                else:
                    self.v2.subs(nomv,valv)        
                if str(nomv)=='t':
                        self.t=valv
                if str(nomv)=='a':
                        self.a=valv        
                if str(nomv)=='ac':
                        self.ac=valv
                else:
                    kk=self.ac
                    self.ac=kk.subs(unisymbols(nomv),valv)        
                if str(nomv)=='v':
                        self.v=valv
                if str(nomv)=='w':
                        self.w=valv
                if str(nomv)=='r':
                        self.r=valv        
                if str(nomv)=='w1':
                        self.w1=valv
                if str(nomv)=='w2':
                        self.w2=valv        
                    
                    
                    
                    
            mat1=[nomv]
            mat2=[valv]
            
            
            
            mm1=mat1+mm[0] 
            mm2=mat2+mm[1] 
            self.kvalue=(mm1,mm2)
            kres=valv
            if kunisym:
                kres=unisymbols(kres)
            if krevalue:
                kres=self.revalue(kres)
             
            kres=opemat(kres,kope)
            if kpositive:
                if kres.compare(0)==-1:
                    kres=-1*kres
            if kdisp:
                sR=str(nomv)+' ='        
                display(Math(sR+latex(kres)))
            if kreturn:
                return(kres)
        
    
    def update_kvalue(self):
        qq=len(self.kvalue[0])
        for i in range(qq):
            ks=self.kvalue[0][i]
            kv=self.kvalue[1][i]
            for j in range(qq):
                
                oldVal=self.kvalue[1][j]
                self.kvalue[1][j]=oldVal.subs(ks,kv)
             
    
    def evalue_val(self,kv, kope=''):
        if type(kv)==MyEq:
            eq1=kv.ksym
        else:    
            eq1=kv

        mm=self.kvalue
        mm1=mm[0] 
        mm2=mm[1]
        qq=len(mm1)
        for i in range (qq):
            ksym=mm1[i]
            kval=mm2[i]
            eq1=eq1.subs(ksym,kval)
        kres=eq1
        kres=opemat(kres,kope=kope)
        return(kres)
        
    def define_val(self,ksym,kval,kend=False):
    
        mat1=[ksym]
        
        mat2=[kval]
        kkend=kend
        self.add_value(mat1,mat2,kend=kkend)
        return(kval)
    
    def superstore_val(self,knomv,kvalv,kope='',kdisp=True,keval=True):
            qq2=len(knomv)
            for i in range(qq2):
                nomv=knomv[i]
                valv=kvalv[i]
                self.store_val(nomv,valv)

            if kdisp:
                self.disp_solu()
                
    def redefine(self):
        mm=self.kvalue
        mmsym=mm[0]
        mmval=mm[1]
        qq=len(mmsym)
        for i in range(qq):
            kksym=mmsym[i]
            mmval=self.value(mmsym)
            mm[1][i]=mmval          
        self.kvalue=mm        
                      
    def revalue(self,kval):
        mm=self.kvalue
        mmsym=mm[0]
        mmval=mm[1]
        qq=len(mmsym)
        if qq==0:
            return(kval)
        else:
            kres=kval
            for i in range(qq):
                nnsym=mmsym[i]
                nnval=mmval[i]
                 
                if len(fpoly(kres,'free'))>0:
                    try:
                        kres=fpoly(kres,'forze_subs',nnsym,nnval)
                    except:
                        try:
                            kres=kres.subs(nnsym,nnval)
                        except:
                            kres=kres
            return(kres)    
                        
    def value(self,ksymb,keval=True,kope=''):
        mm=self.kvalue
        mmsym=mm[0]
        mmval=mm[1]
        qq=len(mmsym)
        if qq>0:       
            for i in range(qq):
                nnsym=mmsym[i]
                nnval=mmval[i]                
                if ksymb==nnsym:
                    kres=nnval
                    break
                else:
                    kres=ksymb            
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)    
        return(kres)

    def disp_solu(self,kunic=False):
            mm=self.kvalue
            mmsym=mm[0]
            mmval=mm[1]
            qq=len(mmsym)
            if qq==0:
                return() 
            else:
                sR= '= ' 
                for i in range(qq):
                    nnsym=mmsym[i]
                    nnval=mmval[i]
                    if kunic:
                        if len(nnval.free_symbols)==0:
                            display(Math(latex(nnsym)+sR+latex(nnval)))
                    else:
                        display(Math(latex(nnsym)+sR+latex(nnval)))
                
    def change_store_value(self, ksim,kval):
        mm=self.kvalue
        mmsym=mm[0]
        mmval=mm[1]
        qq=len(mmsym)
        if qq==0:
            return(kval)
        else:

            for i in range(qq):
                nnsym=mmsym[i]
                nnval=mmval[i]
                if nnsym==ksim:
                    mmval[i]=kval
            mm=[mmsym,mmval]
            self.kvalue=mm
            
    def drop_val(self,ksym):
        mm=self.kvalue
        mmsym=mm[0]
        mmval=mm[1]
        qq=len(mmsym)
        if qq==0:
            return()
        else:
            tsym=[]
            tval=[]
            
            for i in range(qq):
            
                nsym=mmsym[i]
                nval=mmval[i]
                if nsym!=ksym:
                    tsym.append(nsym)
                    tval.append(nval)
            
            self.kvalue=(tsym,tval)   

    def store_solve(self,ksym,Eq1):
        kres=supersolve(ksym,Eq1,kdisp=False,kreturn=True)
        qq=len(kres)
        for i in range(qq):
            self.store_val(ksym[i],kres[i])

    def copy_symbal(self,obj2,ksymb):
        self.store_val(ksymb,obj2.value(ksymb))
        
    def solvestore(self,ksym,Eq1,kdisp=True):
        if not self.ksym_in_store(ksym):
            kres=csolve(Eq1,ksym)
            self.store_val(ksym,kres,kdisp=kdisp)
        
    def ksym_in_store(self,ksym):
        allsymb=self.kvalue[0]
        if ksym in allsymb:
            return(True)
        else:
            return(False)
            
    def setval2(self,ksym,newval):
        mm=self.F
        qq=len(mm)
        for i in range(qq):
            olds=mm[i][0]
            news=olds.subs(ksym,newval)
            mm[i][0]=news
        self.F=mm
    
    
         


    def show_store(self,kunic=False): # ********   show()    ********** 
        mm=self.kvalue
        mmsym=mm[0]
        mmval=mm[1]
        qq=len(mmsym)
        if qq==0:
            return() 
        else:
            sR= '= ' 
            for i in range(qq):
                nnsym=mmsym[i]
                nnval=mmval[i]
                if kunic:
                    if len(nnval.free_symbols)==0:
                        display(Math(latex(nnsym)+sR+latex(nnval)))
                else:
                    display(Math(latex(nnsym)+sR+latex(nnval)))
  
 
    
    def x_res(self,kope='',keval=True,kunisym=True): # ********   Static Result sum F 
        kres=self.resultante(ktype='x',kope=kope)
        if kunisym:
            kres=unisymbols(kres)
        if keval:
            kres=self.revalue(kres)
           
        kres=opemat(kres,kope=kope)
        return(kres)
        
        
    def x_res_relative(self,kang,kope='s'):
         
        fx= self.x_res()
        fy= self.y_res()
        kres=fx *cos(kang)+fy  *sin(kang)
        kres=opemat(kres,kope=kope)
        return kres
        
    def xy_res(self,kope='',keval=True,kunisym=True):
        kres=self.resultante(ktype='xy',kope=kope)
        if kunisym:
            kres=unisymbols(kres)
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)
        return(kres) 

    def y_res(self,kope='',keval=True,kunisym=True):
        kres=self.resultante(ktype='y',kope=kope)
        if kunisym:
            kres=unisymbols(kres)
        if keval:
            kres=self.revalue(kres)
           
        kres=opemat(kres,kope=kope)
        return(kres)
        
    def y_res_relative(self,kang=0,kope='s'):
        fx= self.x_res()
        fy= self.y_res()
        kres=fy *cos(kang)-fx *sin(kang)
        kres=opemat(kres,kope=kope)
        return kres                
    def tan_res(self):
        
        xres=self.x_res()
        yres=self.y_res()
        p1=Point(0, 0)
        p2=Point(self.x_res(), self.y_res())
        L = Line(p1, p2)
        return(L.slope)
        
    def ang_res(self,kope=''):
        tanr=self.tan_res()
        kres=atan(tanr)
        kres=opemat(kres,kope=kope)
        return(kres)    
        
    def vec_displace(self):
        xx=self.x2-self.x1
        yy=self.y2-self.y1
        return([xx,yy])
        
    def displacement(self,op='',kope=''): 

        if op=='x':
            xx=self.x2-self.x1
            kres=xx
        
        
        elif op=='y':
            yy=self.y2-self.y1
            kres=yy
        else:
            xx=self.x2-self.x1
            yy=self.y2-self.y1
            kres=rpow(xx*xx+yy*yy)
            
        kres=opemat(kres,kope=kope)
        return(kres)
         
    
    def y_direRes_from_x1(self,xx2='',kope=''):
        xx1=self.x1
        yy1=self.y1
        if xx2=='':
            xx2=self.x2
        
        Ll=xx2-xx1
        ttg=self.tan_res()
        kres=yy1+Ll*ttg
        kres=opemat(kres,kope=kope)
        return(kres)    
        y    
        
    def resultante(self,keval=True,ktype='x',kope=''): # usado x x res y res
            mm=self.F 
            qq=len(mm)
            kresx=0
            kresy=0
            for i in range(qq):
                vec=1*mm[i]
                kk=1*vec[0]
                ka=1*vec[1]
                kval=kk*cos(ka)
                kresx+=kval
                kval=kk*sin(ka)
                kresy+=kval
            if   ktype=='x':
                    kres=kresx
            elif ktype=='y':
                    kres=kresy
            else:  
                    kres=rpow(kpow(kresx,2)+kpow(kresy,2))
            if keval:
                kres=self.revalue(kres)    
            kres=opemat(kres,kope=kope)
            return(kres)
    def mov_Eq(self,ksym,kope=''):
        aa=self.a 
        vv=self.v
        gg=self.g
        xx=self.x1
        
        kres= ksym*tan(aa)-gg*ksym*ksym/(2*vv*vv*cos(aa)*cos(aa))-xx
        
        kres=unisymbols(kres)    
        kres=opemat(kres,kope=kope)
        return(kres)
        
    def eQ_trayec(self,ksym,kname='',kope=''):
        kres= self.mov_Eq(ksym=ksym,kope=kope)
        if kname!='':
            ee=MyEq(kres,kname)
            return ee
        else:
            return(kres)
        
    def Fx_direc(self,alpha,kope=''):
        fx=self.x_res()
        fy=self.y_res()
        kres=fy*sin(alpha)+fx*cos(alpha)
        kres=opemat(kres,kope=kope)
        return kres
    def Fy_direc(self,alpha,kope=''):
        fx=self.x_res()
        fy=self.y_res()
        kres=fy*cos(alpha)-fx*sin(alpha)
        kres=opemat(kres,kope=kope)
        return kres
    def Fxy_direc(self,alpha,kope=''):
        fx=self.Fx_direc(alpha=alpha,kope=kope)
        fy=self.Fy_direc(alpha=alpha,kope=kope)
        kres=get_hipo(fx,fy)
        kres=opemat(kres,kope=kope)
        return kres
    def tan_dezplacement(self,kop='tan'):
        xx1,xx2,yy1,yy2=self.x1,self.x2,self.y1,self.y2
        ktan=(yy2-yy1)/(xx2-xx1)
        if kop=='tan':
            return ktan
        else:
            return atan(ktan)        
 
    def geotorque(self,angle,R,kope='',keval=True):
        xx,yy= polar2xy(kR,a1,op='')
        return self.torque(x=xx,y=yy,kope=kope,keval=keval)
        
    def torque(self,x=0,y=0,kope='',keval=True):
        mm=unisymbols(self.F )
        qq=len(mm)
        kt=0
 
        for i in range(qq):
            vec=mm[i]
            mk=vec[0]
            xx=vec[2]
            yy=vec[3]
            aa=vec[1]
            dx=xx-x
            dy=yy-y
            fx=-1*mk*cos(aa)*dy
            kt+=fx
            fy= mk*sin(aa)*dx
            kt+=fy
            
        kres=unisymbols(kt)
        if keval:
            kres=self.revalue(unisymbols(kres))      
        kres=opemat(kres,kope=kope)
        return(kres) 
        
 

    def kine_ac(self,kname='',kope='',keval=True): # return (vf-vi)/t 
        vv1=unisymbols(self.v1)
        vv2=unisymbols(self.v2)
        tt=unisymbols(self.t)
        kres=  (vv2-vv1)/tt
        if keval:
            kres=self.revalue(unisymbols(kres))      
        kres=opemat(kres,kope=kope)
        if kname!='':
            ee=MyEq(kres,kname)
            return ee
        else:
            return(kres)
    
    def kine_vf_t(self,kname='',kope='',keval=True): # return ac*t +vi
        vv1=unisymbols(self.v1)
        aac=unisymbols(self.ac)
        tt=unisymbols(self.t)
        kres=  aac*tt+vv1
        if keval:
            kres=self.revalue(unisymbols(kres))      
        kres=opemat(kres,kope=kope)
        if kname!='':
            ee=MyEq(kres,kname)
            return ee
        else:
            return(kres)
    def kine_vf_a_x(self,kname='',kope='',keval=True,ktype='P'): # return ac*t +vi
        vv1=unisymbols(self.v1)
        aac=unisymbols(self.ac)
        x11= unisymbols(self.x2-self.x1)
        kres=  rpow(2*aac*x11-vv1*vv1)
        if keval:
            kres=self.revalue(unisymbols(kres))      
        kres=opemat(kres,kope=kope)
        if kname!='':
            ee=MyEq(kres,kname,ktype=ktype)
            return ee
        else:
            return(kres) 
    def kine_vf_x(self,kname='',kope='',keval=True): # return ac*t +vi
        xx1=unisymbols(self.x1)
        xx2=unisymbols(self.x2)
        aac=unisymbols(self.ac)
        vv1=unisymbols(self.v1)
        kres=  vv1*vv1+2*aac*(xx2-xx1)
        if keval:
            kres=self.revalue(unisymbols(kres))      
        kres=opemat(kres,kope=kope)
        if kname!='':
            ee=MyEq(kres,kname)
            return ee
        else:
            return(kres)
        
    def kine_x2(self,kname='',kope='',keval=True): # return ac*t +vi
        xx1=unisymbols(self.x1)
        tt=unisymbols(self.t)
        aac=unisymbols(self.ac)
        vv1=unisymbols(self.v1)
        kres=  xx1+vv1*tt+aac*tt*tt/2
        if keval:
            kres=self.revalue(unisymbols(kres))      
        kres=opemat(kres,kope=kope)
        if kname!='':
            ee=MyEq(kres,kname)
            return ee
        else:
            return(kres)    
     
    def kine_t(self,kname='',kope='',keval=True): # return ac*t +vi
        vv1=unisymbols(self.v1)
        aac=unisymbols(self.ac)
        vv2=unisymbols(self.v2)
        kres=  (vv2-vv1)/aac
        if keval:
            kres=self.revalue(unisymbols(kres))      
        kres=opemat(kres,kope=kope)
        if kname!='':
            ee=MyEq(kres,kname)
            return ee
        else:
            return(kres)  
        
        
 

    def simple_ac(self,ktype='x',kope='',keval=False):
        if ktype=='y':
            fr=self.y_res()
        elif ktype=='xy':
            fr=self.xy_res()    
        else:
            fr=self.x_res()
        mm=self.m
        kres=fr/mm
        kres=opemat(kres,kope=kope)
        if keval:
            mmval=self.kvalue
            kres=s_subs(kres,mmval)
        mtype=['','x','y','xy']
        if ktype not in mtype:
            ee=MyEq(unisymbols(kres),ktype)
            return ee
        else:    
            return(unisymbols(kres))
        
    def accele_due_Vels(self,t='',kope='',keval=True):
        if t=='':
            t=self.t 
            
        kres=(self.v2-self.v1)/t  
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)
        return(kres)
        
        
    def work(self,op='',kope=''):
        if op=='x':
            xx=self.x2-self.x1
            kres=xx*self.xy_res()
        elif op=='y':
            yy=self.y2-self.y1
            kres=yy*self.xy_res()
        else:
            xxyy=self.displacement()
            kres=xxyy*self.xy_res()
            
        kres=opemat(kres,kope=kope)
        return(kres)
        
        
    def work_X(self):
         
        Fr=self.resultante(ktype='x')
        d1=self.x1
        d2=self.x2
        dd=d2-d1
        kres=Fr*dd
        return(kres)
        
    def work_Y(self):
         
        Fr=self.resultante(ktype='y')
        d1=self.y1
        d2=self.y2
        dd=d2-d1
        kres=Fr*dd
        return(kres)
        
        
    def work_due_forza(self, op='',kope=''):
        if op=='x':
            Fr=self.resultante(ktype='x')
            d1=self.x1
            d2=self.x2
            dd=d2-d1
            kres=Fr*dd
        elif op=='y':
            Fr=self.resultante(ktype='y')
            d1=self.y1
            d2=self.y2
            dd=d2-d1
            kres=Fr*dd
        else:
            wx=self.work_due_forza('x')
            wy=self.work_due_forza('y')
            kres=rpow(kpow(wx,2)+kpow(wy,2))
        kres=opemat(kres,kope)
        return(kres) 
     
    def add_central_force(self,ksym):
        self.Ti=ksym

    def add_cg(self, km,x=0,y=0):
        
        mm=self.cG
        mm.append([km,x,y])
        self.cG=mm
        
    def get_cg(self,ktype='',kope='',kshow=False):
        return self.cgravedad(ktype=ktype,kope=kope,kshow=kshow)
        
    def cgravedad(self,ktype='',kname='',kope='',kshow=False):
         
        if ktype=='help':
            listOfImageNames = ['Phelp/cg1.png','Phelp/cg2.png',]
            for imageName in listOfImageNames:
                display(Image(filename=imageName))
            return
        mm=self.cG
         
        Tm=0
        Tx=0
        Ty=0

        for i in mm:
            Tm+=i[0]
            Tx+=i[0]*i[1]
            Ty+=i[0]*i[2]
        Cx=opemat(Tx/Tm,kope=kope)
        Cy=opemat(Ty/Tm,kope=kope)
        
        kres=(Cx,Cy)
        if ktype=='x':
            kres=Cx
        if ktype=='y':
            kres=Cy
        if kname!='':
            ee=MyEq(kres,kname=kname)
            return ee
            
        return kres
            
           
 
    def energia_rot(self,*args,In=''):
        ww1=self.w1
        ww2=self.w2
        ii=self.getInerTotal()
        if In!='':
            ii=In
        er1=ii*ww1*ww1/2
        er2=ii*ww2*ww2/2
        kname=''
        kop=0
        for i in args:
            if type(i)==str:
                kname=i
            if type(i)==int:
                kop=i

        if kop==1:
            kres= er1
        elif kop==2:
            kres= er2
        else:
            kres=er2-er1
        if kname=='':
            return kres
        else:
            return MyEq(kres,kname=kname)
        
    def energia(self,ktype=' ',keval=True,kope=''):
        mm=self.m
        gg=self.g
        vv1=self.v1
        h1=self.y1
        k1=(mm*vv1**2)/2
        p1=mm*gg*h1
        E1=k1+p1
        vv2=self.v2
        h2=self.y2
        k2=(mm*vv2**2)/2
        p2=mm*gg*h2
        E2=k2+p2
        if ktype=='p1':
            kres=p1
        elif ktype=='p2':
            kres=p2    
        elif ktype=='k1':
            kres=k1
        elif ktype=='k2':
            kres=k2
        elif ktype=='1':
            kres=E1
        elif ktype=='2':
            kres=E2
        elif ktype=='P':
            kres=p1-p2
        elif ktype=='K':
            kres=k1-k2     
        else:
            kres=E2-E1
            
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)
        return(kres)
        
        
    def eRotation(self,ktype='',keval=True,kope=''):
        '''
         ktype 'p1' 'p2' 'ro1' 'ro2' '1' '2' 'P' 'R' 
        ''' 
        return self.energiaRot(ktype=ktype,keval=keval,kope=kope)
        
        
    def energiaRot(self,ktype=' ',keval=True,kope=''):
        '''
         ktype 'p1' 'p2' 'ro1' 'ro2' '1' '2' 'P' 'R':
         
        '''    
        mm=self.m
        gg=self.g
        Ii=self.In
        if self.get_InerTotal()!=0:
            Ii=self.get_InerTotal()
         
        ww1=self.w1
        h1=self.y1
         
        p1=mm*gg*h1
        ro1=Ii*ww1*ww1/2
        E1=ro1+p1
        
        ww2=self.w2
        h2=self.y2
         
        p2=mm*gg*h2
        ro2=Ii*ww2*ww2/2
        E2=ro2+p2
        if ktype=='p1':
            kres=p1
        elif ktype=='p2':
            kres=p2    
        elif ktype=='ro1':
            kres=ro1
        elif ktype=='ro2':
            kres=ro2
        elif ktype=='1':
            kres=E1
        elif ktype=='2':
            kres=E2
        elif ktype=='P':
            kres=p1-p2
        elif ktype=='R':
            kres=ro1-ro2     
        else:
            kres=E2-E1
            
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)
        return(kres)
        
    def Etotal(self,op=''):
        kres=0
        if '1' in op or '2' in op:
            if '1' in op:
                if 'p' in op:
                    kres+=self.energia('p1')
                if 'k' in op:
                    kres+=self.energia('k1')
                if 'r' in op:
                    kres+=self.energiaRot('ro1')
            else: 
                 
                if 'p' in op:
                    kres+=self.energia('p2')
                if 'k' in op:
                    kres+=self.energia('k2')
                if 'r' in op:
                    kres+=self.energiaRot('ro2')        
        if 'K' in op:
            kres+=self.energia('K')
        if 'P' in op:
            kres+=self.energia('P')
        if 'R' in op:
            kres+=self.energiaRot('R')
        if op=='':
            kres=self.energia('p2')+self.energia('k2')+self.energiaRot('ro2')
            kres=kres-self.energia('p1')-self.energia('k1')-self.energiaRot('ro1')
        return kres    
        
    def work_x(self, kope=''):
        Fr=self.resultante('x')
        d1=self.x1
        d2=self.x2
        dd=d2-d1
        kres=Fr*dd
        kres=opemat(kres,kope)
        return(kres)       
        
 
    def x_vel(self,t=t,kope='',keval=True):
        acc=self.ac
        vv=self.v
        aa=self.a
        vv1=vv*cos(aa)
        tt=t
        kres=tt*acc+vv1
        kres=opemat(kres,kope=kope)
        if keval:
            kres=self.revalue(kres)
        
        return(kres)
        
    def y_vel(self,t=t,kope='',keval=True):
        
        vv=self.v
        aa=self.a
        vv1=vv*sin(aa)
        gg=self.g
        kres=-t*gg+vv1
        kres=opemat(kres,kope=kope)
        if keval:
            kres=self.revalue(kres)
        return(kres)
        
    def y_velfinal(self,kope='',keval=True):
        yy=self.y1-self.y2
        gg=self.g
        aa=self.a
        vv=self.v
        zz=gg*yy/(vv*vv)
        ss=sin(aa)
        ss2=ss*ss
        kres=-vv*rpow(ss2+2*zz,2)
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)    
        return(kres)

    def y_v2(self,kope='',keval=True):
        kres=self.y_vel(0)+2*self.g*(self.y2-self.y1)
        sres=str(kres)
        if sres[0]=='-':
            kres=-1*rpow(-1*kres)
        else:
            kres=rpow(kres)
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)
        return(kres)

    def vel_final(self,kope='',keval=True):
        
        gg=self.g
        aa=self.a
        vv=self.v
        zz=gg*yy/(vv*vv)
        ss=sin(aa)
        ss2=ss*ss
        kres=vv*rpow(1+2*zz,2)
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)    
        return(kres)
        
    def tan_velfinal(self,kope='',keval=True):
        
        gg=self.g
        aa=self.a
        vv=self.v
        zz=gg*yy/(vv*vv)
        ss=sin(aa)
        ss2=ss*ss
        kres=-rpow(ss2+22*zz,2)/cos(aa) 
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)    
        return(kres)    
        
    def xy_vel(self,tt='',kope='',keval=True):
        if tt=='':
            tt=t
        vvx=self.x_vel(t=tt)
        vvy=self.y_vel(t=tt)
        kres=pow(pow(vvx,2) +  pow(vvy,2),S('1/2'))
        if keval:
            kres=self.revalue(kres)
        return(kres)

    def x_pos(self, t=t,kope='',keval=True):

        tt=1*t
        vvx=self.x_vel(t=0)
        aac=self.ac
        xx1=self.x1
        kres=xx1+vvx*tt+aac*tt*tt/2
        kres=opemat(kres,kope=kope)
        if keval:
            kres=self.revalue(kres)
        return(kres)
    
    def y_pos(self, t=t,kope='',keval=True):
        tt=1*t
        vvy=self.y_vel(t=0)
        gg=self.g
        yy1=self.y1
        kres=yy1+vvy*tt-gg*tt*tt/2
        kres=opemat(kres,kope=kope)
        if keval:
            kres=self.revalue(kres)
         
        return(kres)
    
    def y_max(self,t=t,keval=True,kope='',krelative=True):
        gg=self.g
        aa=self.a
        vv=self.v
        kres=powsimp(vv**2)*(powsimp(sin(aa))**2)/(2*gg)
        kres=opemat(kres,kope=kope)
        yy1=self.y1
        if krelative:
            kres=kres+yy1
        if keval:
            kres=self.revalue(kres)
        return(kres)
  
    def x_max(self,kope='',keval=True):
        yy=self.y1-self.y2
        gg=self.g
        aa=self.a
        vv=self.v
        zz=gg*yy/(vv*vv)
        ss=sin(aa)
        ss2=ss*ss
        rres=rpow(ss2+2*zz,2)
        
        kres=(vv**2)*(ss+   (rpow(ss2+2*zz,2))    )*cos(aa)/gg
        
        
        
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)    
        return(kres)    

    def t_fly(self,kope='',keval=True):
        vv=self.v
        gg=self.g 
        aa=self.a
        yy1=self.y1
        yy2=self.y2
        if yy1==0 and yy2==0:
            kres=2*vv*sin(aa)/gg
        else:
            kk=2*gg*yy1 - 2*gg*yy2 + vv**2*sin(aa)**2
            kkk=opemat(rpow(kk,2),'r')
            kres= (vv*sin(aa)+kkk)/gg
        

        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)    
        return(kres)

    def tanVel(self,t=t,keval=True,kope='',ktan=True): 
        return self.tan_angle_in_t(t=t,keval=keval,kope=kope,ktan=ktan)
    
    def tan_to_x_max(self,keval=True,kope=''):
        yy=self.y1-self.y2
        gg=self.g
        vv=self.v
        kres=vv/rpow(vv*vv+2*yy*gg)
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)    
        return(kres)
         
        
    def tan_angle_in_t(self,t=t,keval=True,kope='',ktan=True):
        tt=t
        vyy=self.y_vel(t=tt)
        vxx=self.x_vel()
        if ktan:
            kres=(vyy/vxx)
        else:
            kres=atan(vyy/vxx)
        kres=opemat(kres,kope=kope)
        if keval:
            kres=self.revalue(kres)
        return(kres)
        
    def vel_comp(self,op='',kope=''): # find mag vec velocity 
        vvx=kpow(self.vx,2)
        vvy=kpow(self.vy,2)
        vvxy=kpow(self.vxy,2)
        if op=='x':
            kres=rpow(vvxy-vvy)
        elif op=='y':
            kres=rpow(vvxy-vvx)
        else :
            kres=rpow(vxy+vvy)
        kres=opemat(kres,kope)
        return(kres)
    
    def simple_x_pos(self,v1='',ac='',t='',kope='',keval=True):
        if v1=='':
            v1=self.v1
        if ac=='':
            ac=self.ac
        if t=='':
            t=self.t
        kres=v1*t+ac*t*t/2 
        
        kres=opemat(kres,kope=kope)
        if keval:
            kres=self.revalue(kres)
        return(kres)
    
    def x_time(self,t=t):
        return csolve(self.x_pos(t)-self.x2,t)
        
    def y_time(self,t=t):
        return csolveR(self.y_pos(t)-self.y2,t)    
        

    def fradial_ang(self,kope=''):
        mm=self.m
        aa=self.aw
        kres=mm*aa 
        kres=opemat(kres,kope=kope)
        return kres
    
    def fradial_w(self,kope=''):
        mm=self.m
        ww=self.w
        rr=self.r
        kres=mm*ww*ww*rr 
        kres=opemat(kres,kope=kope)
        rr=self.r
        return kres    
        
    def fcentripete(self,kang,v='',w='', r='',m='',g='',kop='Fct',Fin=1):
        
        
        if r=='':
            r=self.r
        if m=='':
            m=self.m
        if g=='':
            ww=self.g
        if v=='' and w=='':
            v=self.v
        if v=='' and w!='':
            v=w*r
        if kop=='Fcn':
            kres=Fin*m*v*v/r
        if kop=='Fcx':
            kres=kres=self.fcen_res(kang)
        if kop=='Fct':
            kres=self.fcen_res(kang)+Fin*m*v*v/r
        return kres     
            
    
    def fcentripeta_res(self,kangle='',keval=True,kope=''):
        F1=self.fcentripeta(kangle=kangle,keval=keval,kope=kope)
        F2=self.Fc
        return F2+F1
        
    def set_forza_central(self,fc=0):
        self.Fc=fc
    def fradial_res(self,kang=0,kope=''):
        fext=self.fcentripeta( kang)
        finer=self.Fc
        mm=self.m
        rr=self.r
        vv=self.v
        ftang=mm*vv*vv/rr
        kres=finer+fext-ftang
        kres=opemat(kres,kope=kope)
        return kres
    
        
    def add_forzac(self, kval,kang,x=0,y=0,s='r'): 
            
            if s=='s':
                kang=sex2rad(kang)
            mm=self.Fc
            qq=len(mm)  
            if qq>0:
                veri=[]
                for i in range(qq):
                    veri.append(mm[i][0])
                if kval not in veri:    
                    mm.append([kval,kang,x,y])
            else: 
                mm.append([kval,kang,x,y])        
            self.Fc=mm
      
    
        
    def fcentripeta(self,kname='',kangle='',keval=True,kope='',kdire=1):
        '''
            fcentripeta()     , return Radial Force - Forza due Velocity
            fcentripeta(alpha), return Sum all Radial Force in alpha pos
            fcentripeta('v')  , return  fcentrip in velocity function m*v*v/R 
            fcentripeta('w')  , return  fcentrip in ang velocity  function m*v*v*R 
        '''    
        
        if kangle=='':
            
            kres=kdire*self.simple_Fcentripeta()
             
        
        else: 
            fx=self.x_res()*cos(kangle)
            fy=self.y_res()*sin(kangle)
            
            kres= fx-fy    
        kres=unisymbols(kres)
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)
        if kname!='':
            ee=MyEq(kres,kname=kname)
            return ee
        else:
            return(kres)
        
            

    def fcirc_res_due_vtan(self,ang):
        F1=self.fcentripeta(ang)
        F2=self.Ti
        F3=F2+F1
        ve=self.v
        mm=self.m
        rr=self.r
        kres=F3- mm*ve*ve/rr
        return kres
    def fcirc_res_due_w(self,ang):
        F1=self.fcentripeta(ang)
        F2=self.Ti
        F3=F2+F1
        ve=self.w
        mm=self.m
        rr=self.r
        kres=F3- mm*ve*ve*rr
        return kres    
   
            
    def Fcentrip_v(self,kope='',keval=True):
        mm=self.m 
        rr=self.r
        vv=self.v
        kres=unisymbols(mm*vv*vv/rr)
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)
        return kres
    
    def Fcentrip_w(self,kope='',keval=True):
        mm=self.m 
        rr=self.r
        ww=self.w
        kres=unisymbols(mm*ww*ww*rr)
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)
        return kres
        
        
     
    def c_res1(self,aac=0):
        xr=self.x_res()
        yr=self.y_res()
        
        ffc=-yr*sin(aac)+xr*cos(aac)
        ffc=self.revalue(ffc)
        fft=-yr*cos(aac)+xr*sin(aac)
        fft=self.revalue(fft)
        return sE(['ffc=',ffc,' fft=',fft])
    def c_res(self,aac=''):
        return self.fcen_res(acc=acc)
        
    def fcen_res(self,aac=''):  #return centripete force in pos aac angle due forces added
        if aac=='':
            beta=symbols('beta')
            aac=beta
        xr=self.x_res()
        yr=self.y_res()
        ffc=-yr*sin(aac)+xr*cos(aac)
        ffc=self.revalue(ffc)
        return  ffc+self.Ti
        
    def Fsum_radial(self,aac=0):  #return centripete force in pos aac angle due forces added
        if aac=='':
            beta=symbols('beta')
            aac=beta
        xr=self.x_res()
        yr=self.y_res()
        ffc=-yr*sin(aac)+xr*cos(aac)
        ffc=self.revalue(ffc)
        return  ffc+self.Ti     
        
    def ftan_res(self,aac=0,kabs=False):  #return tangential force in pos aac angle due forces added
        xr=self.x_res()
        yr=self.y_res()
        
        fft=-yr*cos(aac)+xr*sin(aac)
        fft=self.revalue(fft)
        if kabs:
             
            return  abs(fft)
             
        else:
             
            return  fft    
            
    def simple_vtan(self):
        ww=self.w
        rr=self.r
        kres=ww*rr
        kres=self.revalue(kres)
        
        return(kres)
        
        
    def fcen_v(self):
        return self.fcen_due_vel_tang()
        
    def fcen_due_vel_tang(self):
        return self.simple_Fcentripeta('v')
    
    def fcen_w(self):
        return self.simple_Fcentripeta('w')
    
    def fcen_due_vel_ang(self):
        return self.simple_Fcentripeta('w')
        
    def simple_Fcentripeta(self,ktype='w'):
        mm=self.m
        ww=self.w
        rr=self.r
        vv=self.v
        kres=mm*ww*ww*rr
        if ktype=='v':
            kres=mm*vv*vv/rr
            

        kres=self.revalue(kres)
        
        return(kres)
    
    
    def Fc_dueV(self,kope='',keval=True): # retun m*v*v/r
        mm=self.m 
        vv=self.v 
        rr=self.r 
        kres=mm*vv*vv/rr
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)
        return kres
    
    def xdueW(self,kope='',keval=True):
        mm=self.m 
        ww=self.w
        rr=self.r 
        kres=mm*ww*ww*rr
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)
        return kres
    
   
         
      
        
    def equiTriangle(self,ksym,kangle,ksym2,ksym3,kstore=False):
        kres1=ksym*cos(kangle)
        kres1=self.revalue(kres1)
        kres2=ksym*sin(kangle)
        kres2=self.revalue(kres2)
        ktext=' >> not stored '
        if kstore:
            self.store_val(ksym2,kres1)
            self.store_val(ksym3,kres2)
            ktext='  >> was stored '
        sE([ksym2,'=',kres1,ktext])
        sE([ksym3,'=',kres2,ktext])
            
#    help
    def help(self,op1=''):
        if op1=='':
            sE([' type .help("kine") to get help with kinematic help'])
            sE([' type .help("static") to get help with static help'])
        
        if op1=='kine':
            display(Image(filename="Phelp/help_k1.png"))
            display(Image(filename="Phelp/help_k2.png"))
            
         
 
    def simplePow(self,t=t,kope='',keval=True,ktype=''):
        tt=t
        ww=self.energia(ktype)
        kres=ww/tt
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope)
        return kres
    def get_InerTotal(self,kname='',var2=''):
        '''
            Return Sum   Inercia of all particles in the system
            before we need add each inerce point whit:
            
            add_Inerce(self,mm='',rr='',L=0,ktype='p'): 
            add_Inerce(mass,radioof mass , L= steiner dist,type Inerce..'p','r'..):
        ''' 
        mm=[]
        Ivec=self.cI
        for i in Ivec:
            mm.append(i)
        mms=sum(mm)
        if kname=='':
            return mms
        else:
            if var2!='':
                return MyEq(mms,kname,var2=var2,ktype='F')
            else:
                return MyEq(mms,kname)
        
        
    def getInerTotal(self,kname='',var2=''):
        mm=[]
        Ivec=self.cI
        for i in Ivec:
            mm.append(i)
        mms=sum(mm)
        if kname=='':
            return mms
        else:
            if var2!='':
                return MyEq(mms,kname,var2=var2,ktype='F')
            else:
                return MyEq(mms,kname)
         
     
        
    
        
    def setIner(self,ktype='',xy=0):
        if ktype=='help':
            listOfImageNames = ['Phelp/inertia.png']
            for imageName in listOfImageNames:
                display(Image(filename=imageName))
        else:
            if type(ktype)!=str:
                self.I_n=ktype
            else:
                mm=self.m 
                rr=kpow(self.r,2)
                K=mm*rr
                vecs=['p','r','v','g','c','s','q']
                vecv=[K,K/12,K/3,K,K/2,2*K/5,K/3]
                self.I_n= vecv[vecs.index(ktype)]+mm*xy*xy 
            return self.I_n
         
        
            
    def ang_moment(self,op='',kope='',keval=True): 
        Ia=self.I_n
        
        if op=='w':
            kres=kpow(self.w,2)*Ia/2
        else:
            kres=Ia*self.aw
            
            
         
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope)
        return kres
     
    def Mo_In_due_geoI(self,kope='',keval=True):
        Inn=unisymbols(self.I_n)
        aw=unisymbols(self.a_w)
         
        kres=Inn*aw
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope)
        return kres
    
    def Mo_In_due_velang(self,kope='',keval=True):
        mm=unisymbols(self.m)
        aw=unisymbols(self.a_w)
        rr=unisymbols(self.r)
         
        kres=mm*aw*rr*rr
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope)
        return kres
    
    def Mo_In_due_torque(self,kope='',keval=True):
        tt=self.torque()
        ww=self.a_w
        kres=tt/ww
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope)
        return kres
        
    
    def eCir(self,vv): # Internal use
        Ii=self.I_n
        rr=self.r
        kres=frs(Ii*vv*vv/(rr*rr),2)
        return unisymbols(kres)
        
    def energiaC(self,kv='',keval=True,kope=''):

        if kv=='':
            kres=self.eCir(self.v2)-self.eCir(self.v1)
        elif kv=='1':
            kres=self.eCir(self.v1)
        elif kv=='2':
            kres=self.eCir(self.v2)
        else :
            kres=self.eCir(kv)
    
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope)
        return kres
        
    def energiaCw(self,kv='',keval=True,kope=''):

        if kv=='':
            kres=self.eCir(self.w2*self.r)-self.eCir(self.w1*self.r)
        elif kv=='1':
            kres=self.eCir(self.w1*self.r)
        elif kv=='2':
            kres=self.eCir(self.w2*self.r)
        else :
            kres=self.eCir(kv)
    
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope)
        return kres    
    
    def Ener_wrotation(self,kop='',keval=True,kope=''):
        Ii=self.I_n
        ec1=(Ii/2)*kpow(self.w1,2)
        ec2=(Ii/2)*kpow(self.w2,2)
        if kop=='1':
            kres=ec1
        elif kop=='2':
            kres=ec2
        else:
            kres=ec2-ec1
        
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope)
        return kres
            
        
    def torqueC(self,keval=True,kope=''):
        Ii=self.I_n
        aw=self.a_w
        rr=self.r
        if self.ac!=ac:
            aw=self.ac/rr
        kres=Ii*aw
            
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope)
        return kres
      
    def help(self,*args):
        if len(args)==0:
            mm=[0,1,2,3]
        else:
            mm=args
        istOfImageNames = ['Phelp/HelpPh0.png','Phelp/HelpPh1.png','Phelp/HelpPh2.png','Phelp/HelpPh3.png']
        for kop in mm:
            display(Image(filename=istOfImageNames[kop]))

         
       
def ToInQ(P,kname='',kope=''):
    if kname == '':
        return MyEq(P.torque()-P.ang_moment(), 'T_oI_m', kope=kope)
    else:
        return MyEq(P.torque()-P.ang_moment(),kname, kope=kope)   

        
def InToQ(P,kname='',kope=''):
    return ToInQ(P=P,kname=kname,kope=kope)

def WoEn(P,kname='',kope=''):
    if kname == '':
        return MyEq(P.energia()-P.x_res()*(P.x2-P.x1) , 'W_oE_e', kope=kope)
    else:
        return MyEq(P.energia()-P.x_res()*(P.x2-P.x1) ,kname, kope=kope)        

def EnWo(P,kname='',kope=''):
    return WoEn(P=P,kname=kname,kope=kope)    


def Eq_aceleration(P,var2=t,kname='',kop='x'):
    if type(P)==mparticle:
        if kop=='y':
            kres=P.y_res()/P.m
        else:
            kres=P.x_res()/P.m
    else:
        kres=P
    
    ee1=MyEq(kres,sD('x',2),ktype='d2',var2=t)
    ee2=MyDiff(kres,kvar=var2,kname=kname,ktype='square')
    return ee1,ee2
    
def accx2diffxQ(ksym,kname='',kope='',Fx=0,m=0):

    if type(ksym)==mparticle:
        Fx=ksym.x_res()
        m=ksym.m
        kres=Fx/m
    else:
        kres=ksym
            
     
    kres=opemat(kres,kope)
    ee1=MyEq(kres,sD('x',2),ktype='d2')
    ee2=MyDiff(kres,t,kname=kname,ktype='square')
    return [ee1,ee2]
def accy2diffyQ(*args,kope=''):
     
    qq=len(args)
    if qq==1:
        P=args[0]
        Fy=P.y_res()
        mm=P.m
        kres=Fy/m
        kres=opemat(kres,kope)
        ee=MyEq(kres,sD('y',2),type='d2')
        return ee        

def sD(kstr,kop=1): #put dot or double dot over symbol name
    if kop==2:
        sd='d\ddot{'+kstr+'_t}'
    else:
        sd='d\dot{'+kstr+'_t}'
    return sd



    
def sym2vec(mm):
    kres=''
    mm2=mm.split()
    for i in mm2:
        kres+=str(i)+'vec'+' ' 
    return( symbols(kres))


class mspring:
    def __init__(self,K=K,X1=X1,X2=X2,X0=X0):
                 
        self.K=K # posicion 1
        self.X0=X0
        self.X1=X1 # altura 1     
        self.X2=X2 # altura 1
        
        
    def xL(self):
        x1=self.X1
        x2=self.X2
        return x2-x1
        
    def Fx(self):
        kk=self.K 
        return kk*self.xL() # Function used by physic_lab
        
    

def polar2xy(R,alpha,op=''):
    xx=R*cos(alpha)
    yy=R*sin(alpha)
    if op=='x':
        return xx
    elif op=='y':
        return yy
    else:
        return (xx,yy)
        
#  Plot 

def eqParametricPlot(e1,e2,ksym,x1,x2,xx):
        xxx=np.linspace(x1,x2,xx)
        X=e1.evalueArray(ksym,xxx)
        Y=e2.evalueArray(ksym,xxx)
        plt.plot(X,Y,scalex=True,)
        if max(Y)>0 and min(Y)<0:
             
          x = [int(min(X)),int(max(X))]
          y = [0,0]
          xres=(x1+x2)/2
          plt.fill_between(x, y, color='0.1')
          
       
def get_Iner_func(ktype='',mm=0,rr=0): # generar dormula de Inercia
    
    kres=0
    if ktype=='':
        return kres
    
    elif ktype=='0': # puntual
        kres=mm*rr*rr
        return kres
        
    elif ktype=='1': # varilla de radio rr que gira  en su centro
        kres=mm*rr*rr/12
        return kres    
        
    elif ktype=='2': # varilla de radio rr que gira en su extremo
        kres=mm*rr*rr/3
        return kres
        
    else:
        return 0
        
def reset_ee(*args):
     
    for i in args:
        i.init=False    










        


    