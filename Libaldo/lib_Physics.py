from sympy import symbols

from IPython.display import Image, display,Math,Markdown
from sympy import *
import numpy as np
import inspect
from functools import wraps 

 
from lib_Mathbasic import *
from lib_Mathematica import *
from lib_Algorith import *
from lib_Exponencial import *
 
  
from lib_MyEq import *
from lib_MyEqEq import *
from lib_MyFunctions import *
from lib_MyIntegral import *
from Libaldomath import show_assignment
'''
| Letra | Descripción     | Copiar |
| :---- | :-------------- | :----- |
| α     | alfa minúscula  | **α**  |
| β     | beta            | **β**  |
| θ     | theta           | **θ**  |
| φ     | phi             | **φ**  |
| ω     | omega           | **ω**  |
| Δ     | delta mayúscula | **Δ**  |
| λ     | lambda          | **λ**  |
''' 

t=symbols('t',positive=True)
m,m1,m2,m3,m4,M,g,x,x0,x1,x2,y,y0,y1,y2,X,Y,a,a1,a2,a3,v,v1,v2,M1,M2,M3, V ,V1 ,V2= symbols('m m1 m2 m3 m4 M g x xo x1 x2 y y0 y1 y2 X Y a a1 a2 a3 v v1 v2 M1 M2 M3 V V1 V2 ')
w,w1,w2,aw, aw1,aw2,F,F1,F2,Rx,Ry,r,r1,r2,R, ax, ax1, ax2,ay, ay1,ay2= symbols('w w1 w2 aw aw1 aw2 F F1 F2 Rx Ry r r1 r2 R ax ax1 ax2 ay ay1 ay2')
 

mu,mu1,mu2,fr,fr1,fr2,f1, f2, f3,N1,N2,N3,Nm, L,L1,L2,h,h1,h2,b,H= symbols('mu  mu1 mu2  f_r fr1 fr2 f1 f2 f3 N1 N2 N3 Nm  L L1 L2 h h1 h2 b H')
Io,R1,R2,R3,Io1,Io2=symbols('Io R1 R2 R3 Io_1 Io_2')

alpha,tetha,ac,at,alpha1,alpha2=symbols('alpha theta a_c a_t alpha1 alpha2')
vx,vy,vxy,Po,Ti,I_n,In,a_w,Fc,a_t=symbols('vx vy vxy Po T_i I_n In a_w F_c  a_t ')
T,T1,T2,t1,t2,t3=symbols('T T1 T2 t1 t2 t3')
K,k,X1,X2,X0,d,W,P=symbols('K k X1 X2 X0 d W P')
rho,z,z1,z2,A,p=symbols('rho z  z1 z2 A p') 
xo,yo,zo,beta,xi,xf=symbols('x_o y_o z_o beta x_i x_f')
yp,xp,pp=symbols("y' x' p'") 
vv,Ma=symbols('vv  M_a') # velocidad al cuadrado

# diff variables
dt=symbols('d_t') 
dm,ds,dx,dy,dz,dt,dr,dh,dL,da,dA,dv,dV,dM,=symbols('dm ds dx dy dz dt dr dh dL da dA dv dV dM')   
aw=symbols('ậ') 
from sympy.physics.mechanics import *

anam=str(alpha)
anam1=str(alpha)+'d'
anam2=str(alpha)+'d2'
 
ax,ay,q,E=symbols('a_x a_y q E') 
class mparticle:
    '''
    P=mparticle()
    P.cmass() : return (x,y) position masss gravity
    
    '''
    def __init__(self,*args,x1=x1,x2=x2,y1=y1,y2=y2,v1=v1,v2=v2,m=m,a=0,ang=0,g=g,v=v,w=w,ac=0,s='r',t=t,r=r,r1=r1,r2=r2,vx='',vy='',vxy=vxy,Ti=Ti,In=In,aw=aw,w1=w1,w2=w2 ,at=at,ax=ax,ay='',typeI='',mu=mu, Nm=Nm,Itype='p',xI=0,dire=0,direx=0,direy=0,vv='',deltat=0.01,posx=0,posy=0,dt=dt,Ma=Ma,type='P',q=q,K=K,Fr=0,tipo=''):
        self.name='P'
        self.tipo=tipo
        if len(args)==1:
            self.name=args[0]
        self.m=m   # mass      
        self.g=g   # gravity   
        self.x1= x1 # x1      start Point >(x1,y1),v1D:\Libaldo\SolveLib\Fin_F\physic_lib.py
        self.y1=y1 # y1 
        self.v1=v1 # v1
        self.x2=x2 # x2      finish  Point >(x2,y12),v2
        self.y2=y2 # y2
        self.v2=v2 # v2
        self.v=v   # Initial Velocity ,  important in parabolic 
        if ang!=0:
            self.ang=ang
            self.a=ang
        elif a!=0:
            self.ang=a 
            self.a=a 
        else:
            self.a=0
            self.ang=0
           
        self.ac=ac # ac=X acceleration, set ac=0, in parabolic case
        self.r=r   # Radio  in circular movements       
        self.r1=r1
        self.r2=r2
        self.F=[]  #matriz de fuerzas
        self.cG=[] # Centro de Gravedad
        self.cI=[] # Inerce points
        self.Fc=[] # Fuerzas circulares
        self.kvalue=([],[])
        self.kvaluen=([],[])
        self.t=t
        self.M=[]
        self.CG=[]
        self.vx=vx
         
        self.vy=vy
             
        self.vxy=vxy
        self.w=w
        self.Po=Po
        self.Ti=Ti
        self.I_n=Io
        self.In=In
        self.Io=Io        
        self.aw=aw
        self.aw=aw
        self.ac=ac
        self.at=at
        self.xI=xI 
        self.Ima=[]
        self.w1=w1
        self.w2=w2
        self.ax=ax
        self.ay=ay
        if ay=='':
            self.ay=self.g
        self.typeI=typeI
        self.mu=mu
        self.Nm=Nm
        self.dire=dire
        self.direx=direx
        self.direy=direy
        self.vv=vv
        self.deltat=deltat
        self.posx=posx        
        self.posy=posy
        self.dt=dt
        self.Ma=Ma
        self.type=type
        self.Q=q
        self.K=K
        self.Fr=Fr 
 
        
 
        
    '''
    def add_vector2(self,x=0,y=0,vx=0,vy=0,ang=0):
        
        x=0,y=0,vx=0,vy=0
        if vx==0 and vy==0:
            self.add_forza(mag,ang,x=x,y=y)
        else:
            self.add_forza(vx,0,x=x,y=y)
            self.add_forza(vy,pi/2,x=x,y=y)
    '''        
    def add_vector(self,*args,ang=0):
        '''
           add_vector(x1,y1,x2,y2)
           add_vector(x1,y1,magnitud,ang=alpha) 
        '''
        if len(args)==2:
            vx=args[0]
            vy=args[1]
            x=0
            y=0
            self.add_forza(vx,0,x=x,y=y)
            self.add_forza(vy,pi/2,x=x,y=y)
            
             
        elif len(args)==4:
            x,y,vx,vy=args
            self.add_forza(vx,0,x=x,y=y)
            self.add_forza(vy,pi/2,x=x,y=y)
        
        else:  
            self.add_forza(args[2],ang,x=args[0],y=args[1])
            
 
    def add_Normal(self,n=1):
        if n==1:
            self.add_forza(N1,pi/2)
            display(Math('Normal :'+latex(N1)+', '+latex(pi/2)))
        elif n==2:
            self.add_forza(N2,pi/2)
            display(Math('Normal :'+latex(N2)+', '+latex(pi/2)))
        else:
            self.add_forza(N3,pi/2)
            display(Math('Normal :'+latex(N3)+', '+latex(pi/2)))
            
    def add_Weight(self):
        self.add_forza(self.m*self.g,-pi/2)
        display(Math('Weight :'+latex(self.m*self.g)+', '+latex(pi/2)))
             
             
    def add_forza(self, *args):

        mm = self.F

        # caso 1: F,alpha  o  F,alpha,x,y
        if not isinstance(args[0], (tuple,list)):

            if len(args) == 2:
                F,alpha = args
                x,y = 0,0

            elif len(args) == 4:
                F,alpha,x,y = args

            else:
                raise ValueError("Use (F,alpha) o (F,alpha,x,y)")

            mm.append([F,alpha,x,y])

        # caso 2: tuplas
        else:

            for t in args:

                if len(t)==2:
                    F,alpha = t
                    x,y = 0,0

                elif len(t)==4:
                    F,alpha,x,y = t

                else:
                    raise ValueError("Tupla debe ser (F,alpha) o (F,alpha,x,y)")

                mm.append([F,alpha,x,y])

        self.F = mm
    
    #####  Resultantes  #####        
    def xres(self,*args,**kwargs):
        name, alpha = pick2data(*args)

        kres = 0
        for val, ang, *_ in self.F:
            kres += val*cos(ang - alpha)
            kres = real_subs(kres, **kwargs)

        kres = tsimplify(kres)
        return out2data(kres,name)


    def yres(self,*args,**kwargs):
        name, alpha = pick2data(*args)

        kres = 0
        for val, ang, *_ in self.F:
            kres += val*sin(ang - alpha)
            kres = real_subs(kres, **kwargs)

        kres = tsimplify(kres)
        return out2data(kres,name)    

    

        


    def fieldforce(self,x,y,*args,**kwargs):
        qq=self.Q
        x1=x
        y1=y
        x2=self.x1
        y2=self.y1        
        dd2=(x2-x1)**2+(y2-y1)**2
        kres=K*qq/dd2
        kres2=kres*signo(kres)
        if signo(qq) ==-1:
            fy=cfrac(kres2*(y2-y1),sqrt(dd2))
            fx=cfrac(kres2*(x2-x1),sqrt(dd2))
        else:
            fy=-cfrac(kres2*(y2-y1),sqrt(dd2))
            fx=-cfrac(kres2*(x2-x1),sqrt(dd2))
        fx=real_subs(fx,**kwargs) 
        fy=real_subs(fy,**kwargs)
        if 'modulo' in args:
            return rsimplify(sqrt(fx*fx+fy*fy))
        elif 'angulo' in args:
            return atan(fy/fx)
        elif 'matrix' in args:
            return Matrix([fx,fy]).T
        else:    
            return fx,fy
        
    def xfieldforce(self,x,y,**kwargs):
        fxx,fyy=self.fieldforce(x,y,**kwargs)
        return fxx
    def yfieldforce(self,x,y,**kwargs):
        fxx,fyy=self.fieldforce(x,y,**kwargs)
        return fyy        
        
        

    
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
        
        
        
        
 
 
     
        
     
 
    def xy_accel(self,*args):
        mm=self.m
        kname=''
        kdire=0
        for i in args:
            if type(i)==str:
                kname=i 
            else:
                kdire=i
        kres=self.fx_relative(kdire)
        if kname!='':
            return MyEq(kres/mm,kname=kname)
        else:
            return unisymbols(kres/mm) 
 
        
 
        
    def modres(self,kname=''):
        
        kres=self.get_resultante()
        try:
            kres=killsqrtpow(kres)
        except:
            pass
        if kname!='':
            return MyEq(kres,kname)
        else:
            return kres
            
 
        
        
    def xy_res(self,**kwargs):
        kres=sqrt((self.xres())**2+(self.yres())**2)
        kres=real_subs(kres,**kwargs)
        return kres
        # kres=self.resultante(ktype='xy',kope=kope)
        # if kunisym:
            # kres=unisymbols(kres)
        # if keval:
            # kres=self.revalue(kres)
        #  
        # return(kres)
        
  
            
  
 
        
  
     
    def Fxy(self,kname='Fxy',kshow=True):
        fx=self.Fx(kshow=False)
        fy=self.Fy(kshow=False)
        fxy=rsimplify(sqrt(fx.ksym**2+fy.ksym**2))
        ee=MyEq(fxy,kname=kname,kshow=kshow)
        return ee
        
    def tan_res(self):
        
        xres=self.xres()
        yres=self.yres()
        p1=Point(0, 0)
        p2=Point(self.xres(), self.yres())
        L = Line(p1, p2)
        return(L.slope)
        
    def ang_res(self,kope=''):
        tanr=self.tan_res()
        kres=atan(tanr)
        
         
        return(kres)    
    def rendervector(self,alpha=pi/3):
        datav=pdatagraf(self,alpha=alpha)
        name=self.name
        rendervector(name,*datav,tipo=self.tipo) 
        
    def vec_displace(self):  #lib_Trajectory.py
        xx=self.x2-self.x1
        yy=self.y2-self.y1
        return([xx,yy])
        
    def displacement(self,op='',kope=''):   #lib_Trajectory.py

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
            
         
        return(kres)
         
    
    def y_direRes_from_x1(self,xx2='',kope=''):  #lib_Trajectory.py
        xx1=self.x1
        yy1=self.y1
        if xx2=='':
            xx2=self.x2
        
        Ll=xx2-xx1
        ttg=self.tan_res()
        kres=yy1+Ll*ttg
         
        return(kres)    
              
    def get_resultante(self):  
        xx=self.xres()
        yy=self.yres()
        kres=get_hipo(xx,yy)
        return kres
    
    @show_assignment
    def cgravity(self,*args,c):
        return self.gcenter(*args)
    @show_assignment    
    def cmass(self,*args,**kwargs):
        return self.gcenter(*args,**kwargs)  
    @show_assignment    
    def gcenter(self,*args,**kwargs):
        tipo=''
         
        mm=self.cG

        Tm=0
        Tx=0
        Ty=0

        for i in mm:
            Tm+=i[0]
            Tx+=i[0]*i[1]
            Ty+=i[0]*i[2]

        Cx= cf(Tx,Tm)
        Cy= cf(Ty,Tm)
         
        kres=(Cx,Cy)
        if 'x' in args:
            return Cx
        elif 'y' in args:
            return Cy
        else:
            return (Cx,Cy)  
        
    def get_cg(self,show=True):
        
        return self.cgravedad(ktype=ktype,kope=kope,kshow=kshow)    
        
    def resultante(self,keval=True,ktype='x',kope='',**kwargs): # usado x x res y res
            if len(kwargs)>0:
                mm=subsvec(self.F,**kwargs)
            else:    
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
             
            return(kres)
    def mov_Eq(self,ksym,kope=''):
        aa=self.a 
        vv=self.v
        gg=self.g
        xx=self.x1
        
        kres= ksym*tan(aa)-gg*ksym*ksym/(2*vv*vv*cos(aa)*cos(aa))-xx
        
        kres=unisymbols(kres)    
         
        return(kres)
        
    def eQ_trayec(self,ksym,kname='',kope=''): #lib_Trajectory.py
        kres= self.mov_Eq(ksym=ksym,kope=kope)
        if kname!='':
            ee=MyEq(kres,kname)
            return ee
        else:
            return(kres)
        
    def Fx_direc(self,alpha,kope=''): #lib_Trajectory.py
        fx=self.xres()
        fy=self.yres()
        kres=fy*sin(alpha)+fx*cos(alpha)
         
        return kres
    def Fy_direc(self,alpha,kope=''): #lib_Trajectory.py
        fx=self.xres()
        fy=self.yres()
        kres=fy*cos(alpha)-fx*sin(alpha)
         
        return kres
    def Fxy_direc(self,alpha,kope=''): #lib_Trajectory.py
        fx=self.Fx_direc(alpha=alpha,kope=kope)
        fy=self.Fy_direc(alpha=alpha,kope=kope)
        kres=get_hipo(fx,fy)
         
        return kres
    def tan_dezplacement(self,kop='tan'): #lib_Trajectory.py
        xx1,xx2,yy1,yy2=self.x1,self.x2,self.y1,self.y2
        ktan=(yy2-yy1)/(xx2-xx1)
        if kop=='tan':
            return ktan
        else:
            return atan(ktan)        
 
    def geotorque(self,angle,R,kope='',keval=True):
        xx,yy= polar2xy(kR,a1,op='')
        return self.torque(x=xx,y=yy,kope=kope,keval=keval)
    

    def To(self,*args,**kwargs):
        name,xx,yy=pick3data(*args) 
 
        
        
        ktorque=factor(expand(-1*self.torque(xx,yy)))
        ktorque=real_subs(ktorque,**kwargs)
        ktorque=simplify(ktorque)
        
        return out2data(ktorque,name)
 
             
                    
    def Mo(self,kname='Mo',**kwargs):
        kres=self.I_n
        kres=kres*self.aw
        
         
        if len(kwargs)>0:
            for key, value in kwargs.items():
                kres=real_subs(kres,**kwargs)
        ee= MyEq(kres,kname=kname)
        return ee 

    def MoAng(self,kname='Ma', **kwargs):
        kres=self.I_n
        kres=kres*self.w
        
         
        if len(kwargs)>0:
            for key, value in kwargs.items():
                kres=real_subs(kres,**kwargs)
        ee= MyEq(kres,kname=kname)
        return ee  

    
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
         
        return(kres) 
        
 
     
   
 
    def acdynamic(self,*args):
        name,angle=pick2data(*args)
        kres=cf(self.xres(angle),self.m)
        return out2data(kres,name)
        
    faccel=acdynamic
    dynamicac=acdynamic
    
    
    def ac_due_Forza(kdire='',kname=''):  #lib_Trajectory.py
        if kdire=='':
            kdire=self.dire
        
        frr=self.Fx_relative(dire)
        kres=frr/self.m 
        
        if kname!='':
            ee=MyEq(kres,kname)
            return ee
            
        return kres
    def Xt(self):
        x1=self.x1
         
        v1=self.v1
        alpha=self.a
        ac=self.ac
        kres=x1+v1*t*cos(alpha)+ac*t*t/2
        ee=MyEq(kres,'Xt',var=t)
        return ee
    
    
    
    def relative_aceleration(self,*args,**kwargs): #lib_Trajectory.py
        kname=''
        angle=0
        for i in args:
            for i in args:
                if type(i)==str:
                    kname=i 
                else:
                    angle=i 
        F=self.x_res_relative(angle)
        mm=self.m 
        return smartreturn(F/mm,kname)
        
        
              
        
    def simple_ac(self,kdire='',ktype='x',kope='',keval=False): #lib_Trajectory.py
        if ktype=='y':
            fr=self.yres()
        elif ktype=='xy':
            fr=self.xy_res()    
        else:
            fr=self.xres()
        mm=self.m
        kres=fr/mm
         
        if keval:
            mmval=self.kvalue
            kres=s_subs(kres,mmval)
        mtype=['','x','y','xy']
        if ktype not in mtype:
            ee=MyEq(unisymbols(kres),ktype)
            return ee
        else:    
            return(unisymbols(kres))
        
    def accele_due_Vels(self,t='',kope='',keval=True): #lib_Trajectory.py
        if t=='':
            t=self.t 
            
        kres=(self.v2-self.v1)/t  
        if keval:
            kres=self.revalue(kres)
         
        return(kres)
        
        
    def work(self,**kwargs):
         
        return(self.xwork()+self.ywork())
        
   
 
    def xwork(self,name='',**kwargs):
        x1=self.x1
        x2=self.x2
        L=x2-x1
        F=self.xres()
        kres=F*L
        kres=real_subs(kres,**kwargs)
        if name!='':
            return MyEq(kres,name)
        else:    
            return kres
    def ywork(self,name='',**kwargs):
        y1=self.y1
        y2=self.y2
        L=y2-y1
        F=self.yres()
        kres=F*L
        kres=real_subs(kres,**kwargs)
        if name!='':
            return MyEq(kres,name)
        else:    
            return kres 

        
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
         
        return(kres) 
     
    def add_central_force(self,ksym):
        self.Ti=ksym
        
    def add_rforce(self,ksym):
        self.Ti=ksym
        
    def add_cg(self, km,x=0,y=0):
        
        mm=self.cG
        mm.append([km,x,y])
        self.cG=mm
        
    def get_cg(self,ktype='',kope='',kshow=False):
        return self.cgravedad(ktype=ktype,kope=kope,kshow=kshow)
  
            
           
   
    def eQenergy(self, name='', op='', **kwargs):
        expr = self.energy(op, **kwargs)
        if name:
            return MyEq(expr, name)
        return expr
    def energy(self,op='',**kwargs):
        '''
        K=kinetic,P=potential,R=rotation, 
        1: initial, 2: final
        op='','KP','KP1','R1'  etc
        '''
        kres=0
        mm=self.m
        gg=self.g
        v1=self.v1
        v2=self.v2
        h1=self.y1
        h2=self.y2
        w1=self.w1
        w2=self.w2
        rr=self.r
        Io=self.Io
        k1=mm*v1*v1/2
        k2=mm*v2*v2/2
        p1=mm*gg*h1
        p2=mm*gg*h2
        er1=Io*w1*w1/2
        er2=Io*w2*w2/2
        kres=0
        K=0
        P=0
        R=0
        if '1' in op:
            K=k1
            P=p1
            R=er1
        elif '2' in op:
            K=k2
            P=p2
            R=er2
        else:
            K=k2-k1
            P=p2-p1
            R=er2-er1
        kres=0    
        if 'K' in op:
            kres=kres+K
        if 'P' in op:
            kres=kres+P
        if 'R' in op:
            kres=kres+R
        if op=='':
            kres=P+K
        if kwargs:
            kres=subskwargs(kres,**kwargs)
        return kres    
        
         
        
        
    def eQtraslation(self,*args,**kwargs):
        name,angle=pick2data(*args)
        F=self.xres(angle)
        aa=self.ac
        mm=self.m
        kres= F-aa*mm  
  
        kres=real_subs(kres,**kwargs)
        if name=='':
            return kres
        else:
            return MyEq(kres,name)    
        
    def eQrotation(self,*args,In=''):
        name,x1,x2=pick3data(*args)
        To=self.To(x1,x2)
        if In=='':
            In=self.In
        aw=self.aw    
        kres=To-aw*In
        if name!='':
            return MyEq(kres,name)
        else:
            return kres         
    
        
    def eQmovradial(self,ang=0,negative=False):
        fx=self.Fx_relative(ang,kshow=False)
        mm=self.m
        rr=self.r
        ww=self.w
        if negative:
            return MQ(mm*rr*ww*ww,fx)
        else:    
            return MQ(-1*mm*rr*ww*ww,fx)   
    def eQmovtangencial(self,ang=0,negative=False):
        fy=self.Fy_relative(ang,kshow=False)
        mm=self.m
        rr=self.r
        aw=self.aw
        if negative:
            return MQ(1*mm*rr*aw,fy)
        else:    
            return MQ(-1*mm*rr*aw,fy)
        
    def eRotation(self,ktype='',keval=True,kope=''):
        '''
         ktype 'p1' 'p2' 'ro1' 'ro2' '1' '2' 'P' 'R' 
        ''' 
        return self.energiaRot(ktype=ktype,keval=keval,kope=kope)
        
        
    def energiaRot(self,ktype='R',keval=True,kope=''):
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
         
        return(kres)
    
   
    
    
   
    def simple_Ek(self):
        mm=self.m 
        vv=self.v 
        kres=mm*vv*vv/2
        return kres
        
    def simple_Er(self):
        ii=self.getInerTotal()
        ww=self.w
        kres=ii*ww*ww/2
        return kres
        
 
                
        
    ####################################
    ###  Kunematic functions ###########
    ####################################


    def kinematichelp(self):

        display(Markdown("""
    # Projectile Kinematics

    ### Position
    - `xpos(t)` : x position
    - `ypos(t)` : y position

    ### Velocity
    - `xvel(t)` : velocity x
    - `yvel(t)` : velocity y
    - `speed(t)` : speed magnitude

    ### Trajectory
    - `eQtrajectory()` : equation y(x)

    ### Extremes
    - `xmax()` : horizontal range
    - `ymax()` : maximum height
    - `thmax()` : time to max height
    - `tfly()` : flight time

    ### Acceleration
    - `atangencial(t)` : tangential acceleration
    - `anormal(t)` : normal acceleration

    ### Angles
    - `angle(t)` : trajectory angle
    - `angleshot()` : shooting angle
    """))
    def helpkinematic(self):
        print('xpos(t),xmax(),ypos(t),hmax(),xvel(t),yvel(t),vfinal(t),vangle(t),tfly(), xaceleration(t)')
    
    def xpos(self,t=t,**kwargs):  #lib_Trajectory.py
        x1=self.x1
        v1=self.v1
        alpha=self.ang
        ax=self.ax

        x2= x1 + v1*cos(alpha)*t + ax*t**2/2
        x2=real_subs(x2,**kwargs)
        return x2 
        
    def xmax(self,**kwargs):  #lib_Trajectory.py
        v1=self.v1
        alpha=self.ang
        ax=self.ax
        tf=self.tfly()
        xm = v1*cos(alpha)*tf + (ax*tf**2)/2
        x1=self.x1
        x2=real_subs(xm,**kwargs)
        return x1+x2
        
    def ypos(self,t=t,**kwargs):  #lib_Trajectory.py
        v1=self.v1
        alpha=self.ang
        ay=self.ay
        y1=self.y1
        
        y2=y1 + v1*sin(alpha)*t + ay*t**2/2

        y2=real_subs(y2,**kwargs)
        return y2
        
    def ymax(self,**kwargs):
        return self.hmax(**kwargs)
        
    def hmax(self,**kwargs): #lib_Trajectory.py
        y1=self.y1
        v1=self.v1
        a =self.a 
        g =self.g 
        alpha=self.ang
        hm = v1**2*sin(alpha)**2/(-2*g)
        hm=real_subs(hm,**kwargs)
        return y1+hm  
        
    def xvel(self,t=t,**kwargs):  #lib_Trajectory.py
        ax=self.ax
        alpha=self.ang
        v1=self.v1
        
        v2=v1*cos(alpha) + ax*t
        v2=real_subs(v2,**kwargs)
        return v2
    
    def yvel(self,t=t,**kwargs):  #lib_Trajectory.py
        ay=self.ay
        alpha=self.ang
        v1=self.v1
        
        v2=v1*sin(alpha) + ay*t
        v2=real_subs(v2,**kwargs)
        return v2    
 
    def vfinal(self,tt=t,**kwargs):  #lib_Trajectory.py
        vx=self.xvel(tt=tt,**kwargs)   
        vy=self.yvel(tt=tt,**kwargs)
        vf=sqrt(vx*vx+vy*vy)
        vf=real_subs(vf,**kwargs)
        return vf
        
    def angle(self,t=t,**kwargs):  #lib_Trajectory.py
        vx=self.xvel(t)
        vy=self.yvel(t)
        vx=real_subs(vx,**kwargs)
        vy=real_subs(vy,**kwargs)
        
        return atan(cf(vy,vx)) 

    def tfly(self,**kwargs):  #lib_Trajectory.py
        v1=self.v1
        g=self.g
        alpha=self.ang 
        y1=self.y1
        
        tf = (-v1*sin(alpha) - sqrt(v1**2*sin(alpha)**2 - 2*g*y1))/g
        tf=real_subs(tf,**kwargs)
        return tf
        
    def thmax(self,**kwargs):
        v1 = self.v1
        alpha = self.ang
        g = self.g
        
        th = -v1*sin(alpha)/g
        th = real_subs(th,**kwargs)
        return th
        
    def xaceleration(self,tt=t,**kwargs):  #lib_Trajectory.py
        v1=self.v1
        v2=self.v2
        kres= cf(v2-v1,tt)
        kres=real_subs(kres,**kwargs)
        return kres
    ###  end Kinematic functions ###########
     
    def atangencial(self,t=t,**kwargs):
        ax = self.ax
        ay = self.ay
        
        vx = self.xvel(t)
        vy = self.yvel(t)
        
        v = sqrt(vx**2 + vy**2)
        
        at = (ax*vx + ay*vy)/v
        
        at = real_subs(at,**kwargs)
        return at
        
    def anormal(self,t=t,**kwargs):
        ax = self.ax
        ay = self.ay
        
        at = self.atangencial(t)
        
        a = sqrt(ax**2 + ay**2)
        
        an = sqrt(a**2 - at**2)
        
        an = real_subs(an,**kwargs)
        return an    
    def eQtrajectory(self,name='',**kwargs):

        x1 = self.x1
        y1 = self.y1
        v1 = self.v1
        a  = self.ang
        g  = self.g

        y = y1 + (x-x1)*tan(a) + g*(x-x1)**2/(2*v1**2*cos(a)**2)

        y = real_subs(y,**kwargs)
        if name=='':
            return y
        else:
            return MyEq(y,name)
    def angleshot(self,**kwargs):

        x1 = self.x1
        y1 = self.y1
        v1 = self.v1
        g  = self.g

        dx = x - x1
        dy = y - y1

        A = g*dx**2/(2*v1**2)
        B = dx
        C = A - dy

        T1 = (-B + sqrt(B**2 - 4*A*C))/(2*A)
        T2 = (-B - sqrt(B**2 - 4*A*C))/(2*A)

        ang1 = atan(T1)
        ang2 = atan(T2)

        ang1 = real_subs(ang1,**kwargs)
        ang2 = real_subs(ang2,**kwargs)

        return ang1,ang2            
           
    def parabolicEq(self,**kwargs):
        g,v1,a=self.g,self.v1,self.a
        kres=x*tan(a)-g*x*x/(2*v1*v1*(cos(a)**2))
        kres=real_subs(kres,**kwargs)
        ee=MyEq(kres,'y',var=x)
        return ee
    def anglevel(self,t,**kwargs):
        th = atan(self.yvel(t)/self.xvel(t))
        return real_subs(th,**kwargs)
        
    def telemetric(self,t):

        return {
            'x': self.xpos(t),
            'y': self.ypos(t),
            'vx': self.xvel(t),
            'vy': self.yvel(t),
            'v': self.speed(t)
        }    
    def speed(self,t,**kwargs):
        v = sqrt(self.xvel(t)**2 + self.yvel(t)**2)
        return real_subs(v,**kwargs) 
    def circul_desplace(self, t=t,kope='',keval=True):  #lib_Trajectory.py
        ww1=self.w1
        rr=self.r
        tt=1*t
        aaw=self.aw
        kres=ww1*rr*tt+rr*aaw*tt*tt/2
         
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
            
        return(kres)
         
        
    def tan_angle_in_t(self,t=t,keval=True,kope='',ktan=True):
        tt=t
        vyy=self.y_vel(t=tt)
        vxx=self.x_vel()
        if ktan:
            kres=(vyy/vxx)
        else:
            kres=atan(vyy/vxx)
         
        if keval:
            kres=self.revalue(kres)
        return(kres)
        
    def vel_comp(self,op='',kope=''):  #lib_Trajectory.py
        vvx=kpow(self.vx,2)
        vvy=kpow(self.vy,2)
        vvxy=kpow(self.vxy,2)
        if op=='x':
            kres=rpow(vvxy-vvy)
        elif op=='y':
            kres=rpow(vvxy-vvx)
        else :
            kres=rpow(vxy+vvy)
         
        return(kres)
    
    def simple_x_pos(self,v1='',ac='',t='',kope='',keval=True):  #lib_Trajectory.py
        if v1=='':
            v1=self.v1
        if ac=='':
            ac=self.ac
        if t=='':
            t=self.t
        kres=v1*t+ac*t*t/2 
        
         
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
         
        return kres
    
    def fradial_w(self,kope=''):
        mm=self.m
        ww=self.w
        rr=self.r
        kres=mm*ww*ww*rr 
         
        rr=self.r
        return kres    
        
    
            
    
    def fcentripeta_res(self,angle=0,Fr=''):
        '''
        Fr= vector that have direccion whit object and the rotation center , example T
        '''
        
        if Fr=='':
            Fr=self.Fr
        xres=self.xres_relative(angle)
        return Fr-xres
        
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
    def fcentripeta_res(self,angle=0,Fr=''):
        '''
        Fr= vector that have direccion whit object and the rotation center , example T
        '''
        
        if Fr=='':
            Fr=self.Fr
        xres=self.xres_relative(angle)
        return Fr-xres  
    def eQcentripeta(self,*args,angle=0):
        '''
        return Fc=m*(w*r)**2/2
        included 'v' in args: return Fc=m*v*v/2
        '''
        p1=self.fcentripeta_res(angle=angle,Fr=self.Fr)
        if '-1' in args:
            p1=-1*p1
        p2=self.m*(self.w)**2*self.r
        if 'v' in args:
            p2=self.m*self.v**2/self.r
        return MQ(p1,p2)
        
    def fcentripeta(self,*args):
        '''
            fcentripeta()     , return Radial Force - Forza due Velocity
            fcentripeta(alpha), return Sum all Radial Force in alpha pos
            fcentripeta('v')  , return  fcentrip in velocity function m*v*v/R 
            fcentripeta('w')  , return  fcentrip in ang velocity  function m*v*v*R 
        '''    
        angle=''
        dire=1
        name=''
        op=['simplify','factor','expand']
        for data in args:
            if type(data)==str:
                if not data in op:
                    name=data
            else:
                if data==-1:
                    dire=-1
                else:
                    angle=data
     
        if angle=='':
            kres=dire*self.simple_Fcentripeta()
        else: 
            fx=self.xres()*cos(angle)
            fy=self.yres()*sin(angle)
            kres= fx-fy    
        kres=unisymbols(kres)
        kres=transympy(kres,*args)
        if name!='':
            ee=MyEq(kres,name)
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
         
        return kres
    
    def Fcentrip_w(self,kope='',keval=True):
        mm=self.m 
        rr=self.r
        ww=self.w
        kres=unisymbols(mm*ww*ww*rr)
        if keval:
            kres=self.revalue(kres)
         
        return kres
        
    def fcEq(self,a='',kname='eFc'): # return fcentripeta accel centripeta Equation   
        mm=self.m
        ve=self.v
        if self.vv=='':
            vv=ve*ve 
        else:
            vv=self.vv
        
        e1=MyEq(self.fcen_res(a)-mm*vv/2,kname)
        return e1
        
        
    def c_res1(self,aac=0):
        xr=self.xres()
        yr=self.yres()
        
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
        xr=self.xres()
        yr=self.yres()
        ffc=-yr*sin(aac)+xr*cos(aac)
        ffc=self.revalue(ffc)
        return  ffc+self.Ti
        
    def Fsum_radial(self,aac=0):  #return centripete force in pos aac angle due forces added
        if aac=='':
            beta=symbols('beta')
            aac=beta
        xr=self.xres()
        yr=self.yres()
        ffc=-yr*sin(aac)+xr*cos(aac)
        ffc=self.revalue(ffc)
        return  ffc+self.Ti     
        
    def ftan_res(self,aac=0,kabs=False):  #return tangential force in pos aac angle due forces added
        xr=self.xres()
        yr=self.yres()
        
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
        
    def setFcentripeta(self,F):
        self.Fr=F
        if type(F)==MyEq:
            self.Fr=F()
        display(self.Fr)
        
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
         
        return kres
    
    def xdueW(self,kope='',keval=True):
        mm=self.m 
        ww=self.w
        rr=self.r 
        kres=mm*ww*ww*rr
        if keval:
            kres=self.revalue(kres)
         
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
         
     
        
    
        
    def setIner(self,kval):
        self.I_n=kval
        return  self.I_n
         
         
        
            
    def ang_moment(self,op='',kope='',keval=True): 
        Ia=self.I_n
        
        if op=='w':
            kres=kpow(self.w,2)*Ia/2
        else:
            kres=Ia*self.aw
            
            
         
        if keval:
            kres=self.revalue(kres)
         
        return kres
     
    def Mo_In_due_geoI(self,kope='',keval=True):
        Inn=unisymbols(self.I_n)
        aw=unisymbols(self.a_w)
         
        kres=Inn*aw
        if keval:
            kres=self.revalue(kres)
         
        return kres
    
    def Mo_In_due_velang(self,kope='',keval=True):
        mm=unisymbols(self.m)
        aw=unisymbols(self.a_w)
        rr=unisymbols(self.r)
         
        kres=mm*aw*rr*rr
        if keval:
            kres=self.revalue(kres)
         
        return kres
    
    def Mo_In_due_torque(self,kope='',keval=True):
        tt=self.torque()
        ww=self.a_w
        kres=tt/ww
        if keval:
            kres=self.revalue(kres)
         
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
         
        return kres
      
    def help(self,*args):
        if len(args)==0:
            mm=[0,1,2,3]
        else:
            mm=args
        istOfImageNames = ['Phelp/HelpPh0.png','Phelp/HelpPh1.png','Phelp/HelpPh2.png','Phelp/HelpPh3.png']
        for kop in mm:
            display(Image(filename=istOfImageNames[kop]))
    def eQkinematic(self,**kwargs):
        v1=self.v1
        v2=self.v2
        ac=self.ac
        x1=self.x1
        x2=self.x2
        t=self.t
        p1=x1+v1*t+ac*t*t/2
        p1=real_subs(p1,**kwargs)
        p2=p2
        p3=real_subs(p2,**kwargs)
        QQ=MQ(p1,p2)
         
        
             
    def eQdynamicX(self,alpha=0,negative=False):
        F=self.Fx_relative(alpha,kshow=False)
        mm=self.m
        ac=self.ac
        if negative:
            return MQ(F,-1*mm*ac)
        else:
            return MQ(F,mm*ac)

    def eQdynamicY(self,alpha=0,negative=False):
        F=self.Fy_relative(alpha,kshow=False)
        mm=self.m
        ac=self.ac
        if negative:
            return MQ(F,-1*mm*ac)
        else:
            return MQ(F,mm*ac)
    def eQmechanic(self):
        kres=[]
        are=self.dire
        if are==0:
            fx=self.xres()
            fy=self.yres()
        else:
            fx=self.fx_relative(are)
            fy=self.fy_relative(are)
        veriV=fx+fy
        all_var=fpoly(veriV,'list')
        Nm=0
        Nmm=''
        
        if N1 in all_var:
            Nm=N1
            Nmm=str(N1)
        if N2 in all_var:
            Nm=N2
            Nmm=str(N2)
        N0=MyEq(csolve(fy,Nm),kname=Nmm)
        kres.append(N0)
        
        if self.mu!=0:
            fr=symbols('fr')
            fr=MyEq(Nm*self.mu,kname='fr')
            fx.upgrade(fr,kshow=False)
            kres.append(fr)
        if self.ac!=0:
            ac0=symbols(self.ac.name)
             
            ac0=MyEq(fx/self.m,kname=self.ac.name)
            kres.append(ac0)
        else:
            if T in all_var:
                T0=symbols('T')
                T0=fx.solve(T,'T')
                kres.append(T0)
            elif F in all_var:
                F0=symbols('F')
                F0=fx.solve(F,'F')
                kres.append(F0)
        return kres            
    
    def eQ_dVx(self,kname='',**kwargs): #retuen dV eQ dynamic due Forze in x
        mm=self.m
        #Fx=self.xres(**kwargs)
        if self.dire==0:
            fx=self.xres(**kwargs)
        else:
            fx=self.x_res_relative(self.dire)
            fx=real_subs(fx,**kwargs)
        tt=self.dt
        vv=self.v
        kres=fx*tt/mm+vv
        kres=real_subs(kres,**kwargs)
        if kname=='':
            return kres
        else:    
            ee=MyEq(kres,kname)
            return ee

    def eQ_dx_ac(self,kname=''):
         
        tt=self.dt
        vv=self.v
        acc=self.ac
        kres=vv*tt+acc*tt*tt/2
        if kname=='':
            kname='dx'
        ee=MyEq(kres,kname)
        return ee    

    def eQ_dx(self,kname='',**kwargs):    #retuen dx eQ dynamic due Forze in x
        mm=self.m 
        tt=self.dt
        vv=self.v
        if self.dire==0:
            fx=self.xres(**kwargs)
        else:
            fx=self.x_res_relative(self.dire)
            fx=real_subs(fx,**kwargs)
            
        acc=fx/mm
        kres=vv*tt+acc*tt*tt/2
        kres=real_subs(kres,**kwargs)
        if kname=='':
            return kres
        else:    
            ee=MyEq(kres,kname)
            return ee
            

       
       


def ToInQ(P,kname='',kope=''):
    if kname == '':
        return MyEq(P.torque()-P.ang_moment(), 'T_oI_m', kope=kope)
    else:
        return MyEq(P.torque()-P.ang_moment(),kname, kope=kope)   

        
def InToQ(P,kname='',kope=''):
    return ToInQ(P=P,kname=kname,kope=kope)

def WoEn(P,kname='',kope=''):
    if kname == '':
        return MyEq(P.energia()-P.xres()*(P.x2-P.x1) , 'W_oE_e', kope=kope)
    else:
        return MyEq(P.energia()-P.xres()*(P.x2-P.x1) ,kname, kope=kope)        

def EnWo(P,kname='',kope=''):
    return WoEn(P=P,kname=kname,kope=kope)    


def Eq_aceleration(P,var2=t,kname='',kop='x'):
    if type(P)==mparticle:
        if kop=='y':
            kres=P.yres()/P.m
        else:
            kres=P.xres()/P.m
    else:
        kres=P
    
    ee1=MyEq(kres,sD('x',2),ktype='d2',var2=t)
    ee2=MyDiff(kres,kvar=var2,kname=kname,ktype='square')
    return ee1,ee2


    
def accx2diffxQ(ksym,kname='',kope='',Fx=0,m=0):

    if type(ksym)==mparticle:
        Fx=ksym.xres()
        m=ksym.m
        kres=Fx/m
    else:
        kres=ksym
            
     
     
    ee1=MyEq(kres,sD('x',2),ktype='d2')
    ee2=MyDiff(kres,t,kname=kname,ktype='square')
    return [ee1,ee2]
def accy2diffyQ(*args,kope=''):
     
    qq=len(args)
    if qq==1:
        P=args[0]
        Fy=P.yres()
        mm=P.m
        kres=Fy/m
         
        ee=MyEq(kres,sD('y',2),type='d2')
        return ee        

 

    
def sym2vec(mm):
    kres=''
    mm2=mm.split()
    for i in mm2:
        kres+=str(i)+'vec'+' ' 
    return( symbols(kres))

'''
class mspring:
    def __init__(self,K=k,X1=x1,X2=x2,X0=x0):
                 
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
        
'''   

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

def DimEq(val='',kname=''):
    return setDeQ(val=val,kname=kname)
def setDeQ(val='',kname=''):
    if val=="Area" :
        kres= L*L       
    elif val=="Vol":
        kres= L*L*L 
    elif val=="veloc":
        kres= L/T
    elif val=="accel":
        kres= L/T/T 
    elif val=="Forz":
        kres= M*L/T/T
    elif val=="Work":
        kres= M*L*L/T/T 
    elif val=="Energy":
        kres= M*L*L/T/T 
    elif val=="Pot":
        kres= M*L*L/T/T/T 
    elif val=="Pres":
        kres= M/L/T/T
    elif val=="Mom":
        kres= M*L*L/T/T
    elif val=="Freq":
        kres=  1/T
    elif val=="Perd":
        kres=  T
    elif val=="Dens":
        kres=  M/L/L/L 
    elif val=="AngV":
        kres=  1/T
    else:
        print('options : "Area","Vol","veloc","accel","Forz","Work","Energy","Pot","Pres","Mom","Freq","Perd","Dens","AngV"')
        return
    if kname=='':
        return kres
    else:
        return MyEq(kres,kname)
    
def subsvec(vec,**kwargs): # chage Obj.F data
    newk=vec
    newv=[]
    for i in newk:
        newp=[]
        for j in i:
            j=real_subs(j,**kwargs)
            newp.append(j)
        newv.append(newp)
    return newv   



def basicope(kval,*args):
    if 'simplify' in args:
        kval=simplify(kval)
    if 'factor' in args:
        kval=factor(kval) 
    if 'expand' in args:
        kval=factor(kval)
    return kval        



def smartreturn(kres,name='',var=t):
    if name=='':
        return float2int(kres)
    else:
        ee=MyEq(float2int(kres),kname=name,var=var)
        return ee
        
        

# ELECTRIC CHARGUES, FIELD, POTENTIAL, COULUMB....
# electric Field due several charges

def fx_electricfield(x,y,*args,**kwargs):
    kres=0
    for data in args:
        kres=kres+data.xfieldforce(x,y)
    kres=real_subs(kres,**kwargs)
    return kres
    
def fy_electricfield(x,y,*args,**kwargs):
    kres=0
    for data in args:
        kres=kres+data.yfieldforce(x,y)
    kres=real_subs(kres,**kwargs)
    return kres    

def electrifield(x,y,*args,**kwargs):
    '''
    electrifield(x,y,*args,**kwargs)
        x=xpos
        y=ypos
        *args= Q1,Q2,Q3.....
        *args include:
            'modulo': return sqrt(fx*fx+fy*fy)
            'matriz': return Matrix([fx,fy]).T = [fx,fy]
            'angulo': return atan(fy/fx)
            'fx': return fx
            'fy': return fy
            none of the above: return (fx,fy) 
            
    also see:  electrifield, electricPot, electricWork        
    '''    
    ops=['modulo','matrix','angulo','fx','fy']
    vec=[]
    for data in args:
        if not data in ops:
            vec.append(data)
    fx=fx_electricfield(x,y,*vec,**kwargs)
    fy=fy_electricfield(x,y,*vec,**kwargs)
    fx=real_subs(fx,**kwargs) 
    fy=real_subs(fy,**kwargs)
    fx=rsimplify(fx)
    fy=rsimplify(fy)
    if 'modulo' in args:
        return rsimplify(sqrt(fx*fx+fy*fy))
    elif 'angulo' in args:
        return atan(fy/fx)
    elif 'matrix' in args:
        return Matrix([fx,fy]).T
    elif 'fx' in args:
        return fx
    elif 'fy' in args:
        return fy
    else:    
        return fx,fy 


def electricPot(x,y,*args,K=9*10**9):
    '''
    return electric Potencial in th point(x,y) due *args mparticle cahargues
    example:
        Q1=mparticle(q=pow10(6,-8),x1=cfrac(-3,10),y1=0)
        Q2=mparticle(q=pow10(8,-6),x1=cfrac(5,10),y1=0)
        Q1(0,0)..........P(0.3,0)...............Q2(0.8,0) finf P in(0.3,0)
        electricPot(0,0,Q1,Q2) 
        also see:  electrifield, electricPot, electricWork
    '''
    x1=x
    y1=y
    kP=0
    for data in args:
        x2=data.x1
        y2=data.y1
        Qq=data.Q
        Rr=rsimplify(sqrt((x2-x1)**2+(y2-y1)**2))
        kP=kP+Qq*K/Rr
    return kP
    
    
def electricWork(Qq,*args,x2=0,y2=0,K=9*10**9):
    '''
    Qq=mparticle(q=carga, x1=xini,y1=yini)
    *args= mparticle chragues who create Potential
    x2,y2 = posfinal
    Qq, makw Work from (x1,y1) to (x2,y2)
    also see:  electrifield, electricPot, electricWork
    '''
    
    x1=Qq.x1
    y1=Qq.y1
    Va=0
    Vb=0
    for data in args:
        Va+=electricPot(x1,y1,data)
        Vb+=electricPot(x2,y2,data)      
    qres=Qq.Q*(Vb-Va)
    return qres

def pdatagraf(self,alpha=pi/6):

    alpha2=alpha
    vecf = self.F
    dataf = []

    for data in vecf:

        nom = str(data[0])
        dire = data[1]

        if signo(data[0]) == -1:
            nom = str(-1*data[0])
            dire = simplify(dire + pi)

        dire = dire.subs('alpha', alpha2)

        dataf.append((nom, dire))

    return dataf 
 

def rendervector(nombre,*args,subs={},tipo=''):

    cx = -0.2
    cy = -0.2

    L = 1.6

    colors=["red","blue","green","orange","purple","brown","black"]

    fig,ax=plt.subplots(figsize=(5,5))


    # objeto central
    if tipo=='block':

        rect=plt.Rectangle(
            (cx-0.3,cy-0.3),
            0.6,
            0.6,
            fill=False,
            linewidth=2
        )

        ax.add_patch(rect)

    if tipo=='ball':

        circ=plt.Circle(
            (cx,cy),
            0.3,
            fill=False,
            linewidth=2
        )

        ax.add_patch(circ)


    ax.text(cx,cy,nombre,ha='center',va='center')


    # vectores
    for i,data in enumerate(args):

        if len(data)==2:
            name,ang=data
            ang_draw=ang
        else:
            name,ang,ang_draw=data

        try:
            ang=float(sp.N(ang_draw.subs(subs)))
        except:
            ang=float(ang_draw)

        vx=L*np.cos(ang)
        vy=L*np.sin(ang)

        color=colors[i%len(colors)]

        ax.arrow(
            cx,cy,
            vx,vy,
            width=0.035,
            head_width=0.18,
            length_includes_head=True,
            color=color
        )

        ax.text(
            cx+vx*1.2,
            cy+vy*1.2,
            name,
            fontsize=14,
            ha='center',
            color=color
        )


    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_xlim(-2,2)
    ax.set_ylim(-2,2)

    plt.show()    