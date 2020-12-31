from sympy import *
from lib_subclass import *

# Symbols Variable used in Library
# --------------------------------
x1,y1,x2,y2,v1,v2,m,mu,a,g,r,w,v,Ti,ac,alpha=symbols('x1 y1 x2 y2 v1 v2 m mu a g r w v Ti ac alpha')

# Symbols Variable used in Problems aplications
# --------------------------------

# mass symbols
# -------------
m1,m2,m3 =symbols('m1 m2 m3 ')

# Distance symbols
# -------------
L,L1,L2,h1,h2,h3,d1,d2,d3=symbols('L L1 L2 h1 h2 h3 d1 d2 d3 ')


# Angle symbols
# -------------
beta,theta=symbols(' beta theta')

# Ratio symbols
# -------------
r1,r2=symbols('r1 r2')


# Forces symbols
# -------------
F1,F2,Fx,Fy, R1,R2,Rx,Ry,N1,N2,T1,T2,T3,fr, fr1,fr2=symbols('F1 F2 Fx Fy  R1 R2 Rx Ry N1 N2 T1 T2 T3 fr  fr1 fr2')

# Velocity Aceleration symbols
# -------------
V1,V2,a1,a2 =symbols('V1 V2 a1 a2')

# Energy Work
# -------------
E1,E2,W1,W2=symbols('E1 E2 W1 W2') 

# Other symbols
# -------------
mu1,mu2=symbols('mu1,mu2')

class kparticle:
    def __init__(self,x1=x1,y1=y1,x2=x2,y2=y2,v1=v1,v2=v2,m=m,mu=mu,a=alpha,s='r',g=g,r=r,w=w,v=v,Ti=Ti,ac=ac):
                 
        self.x1=x1 # posicion 1
        self.y1=y1 # altura 1
        self.x2=x2 # posicion 2
        self.y2=y2 # altura 2
        self.v1=v1 # velocidad en 1
        self.v2=v2 # velocidad en 2
        self.m=m   # masa 
        self.mu=mu # coeff rozamiento
        self.a=sex2rad_i(a,s) # angulo de direccion
        self.F=[]  #matriz de furzas
        self.g=g   # gravedad
        self.r=r   # radio de Giro
        self.w=sex2rad_i(w,s)  # velocidad angular
        self.v=v   # velocidad tangencial
        self.cG=[] # Centro de Gravedad
        self.Fc=[] # Fuerzas circulares
        self.Ti=Ti
        self.ac=ac
        self.kvalue=([],[])
    

    def add_value(self,mat1,mat2,kend=False):
        mm=self.kvalue
        if kend:
            mm1=mm[0]+mat1 
            mm2=mm[1]+mat2
        else:
            mm1=mat1+mm[0] 
            mm2=mat2+mm[1] 
        self.kvalue=(mm1,mm2)
        
    def add_forza(self, kval,kang,x=0,y=0,s='r'):
        if s=='s':
            kang=sex2rad(kang)
        mm=self.F
        mm.append([kval,kang,x,y])
        self.F=mm
    
    def add_cg(self, km,x=0,y=0):
        
        mm=self.cG
        mm.append([km,x,y])
        self.cG=mm
        
    def add_Fcentripeta(self, km,kang,inceter=True):
        
        mm=self.Fc
        if incenter: 
            kpun=0
        else:
            kpun=1
        mm.append([km,x,y])
        self.cG=mm    
        
    def resultante(self,s='s',kabs=False,ktype='x'):
        mm=self.F 
        qq=len(mm)
        kresx=0
        kresy=0
        for i in range(qq):
            vec=mm[i]
            kk=vec[0]
            ka=vec[1]
            kval=kk*cos(ka)
            kresx+=kval
            kval=kk*sin(ka)
            kresy+=kval
        if ktype=='x':
            kres=kresx
        elif ktype=='y':
            kres=kresy
        else:
            kres=sqrt(kresx**2+kresy**2)
        if s=='e':
           kres=kres.evalf()
        if s=='t':
           kres=trigsimp(kres) 
        return(kres)
        
    def torque(self,x=0,y=0,s='s'):
        mm=self.F 
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
            
        kres=kt    
         
        if s=='e':
           kres=kres.evalf()
        return(kres) 
        
    def changeMatF(self,k1,k2):
        mm=self.F 
        qq=len(mm)
        for i in range(qq):
            vec=mm[i]
            kk=vec[0]
            k2=kk.subs(k1,k2)  
            mm[i][0]=k2
        self.F=mm

    def simple_ac(self,s='s',ktype='x',keval=False):
        if ktype=='y':
            fr=self.resultante(ktype='y')
        else:
            fr=self.resultante()
        mm=self.m
        kres=fr/mm
        if keval:
            mmval=self.kvalue
            kres=s_subs(kres,mmval)
        
        if s=='e':
           kres=kres.evalf()
        return(kres)

    
    def cgravedad(self,s='s',kabs=False,ktype='x'):
        mm=self.cG 
        qq=len(mm)
        kresm=0
        kresx=0
        kresy=0
        for i in range(qq):
            vec=mm[i]
            km=vec[0]
            kx=vec[1]
            ky=vec[2]
            kresm+=km
            kresx+=km*kx
            kresy+=km*ky
        kresx=kresx/kresm
        kresy=kresy/kresm
        
        if ktype=='y':
            kres=kresy
        
        else:
            kres=kresx
        if s=='e':
           kres=kres.evalf()
        if s=='t':
           kres=trigsimp(kres) 
        return(kres)
    
    def vel_ang(self,ktype='v',s='s'):
        if ktype=='w':
            kres=self.w
        else:
            rr=self.r
            vv=self.v
            kres=vv/rr
        if s=='e':
           kres=kres.evalf()
        return(kres)   
    
    def vel_tan(self,ktype='w',s='s'):
        if ktype=='v':
            kres=self.v
        else:
            rr=self.r
            ww=self.w
            kres=ww/rr
        if s=='e':
           kres=kres.evalf()
        return(kres)
     
    def simple_acentripeta(self,ktype='w',s='s'):
        rr=self.r
        if ktype=='v':
            vv=self.v
            kres=vv*vv/rr
        else:
            ww=self.w
            kres=ww*ww*rr
        if s=='e':
           kres=kres.evalf()
        return(kres) 

    def simple_fcentripeta(self,s='s',ktype='v'):
        mm=self.m
        rr=self.r
        if ktype=='w':
            ww=self.w
            kres=mm*rr*ww**2
        else:
            vv=self.v
            kres=mm*(vv**2)/rr
        
        if s=='e':
           kres=kres.evalf()
        return(kres)


    def w_tangencial(self,aw=0, s='r'):
        mm=self.m
        gg=self.g
        aa=aw
        if s=='s':
            aa=sex2rad(aa)
        kres=mm*gg*cos(aa)
        return(kres)
    
    def w_incenter(self,aw=0, s='r'):
        mm=self.m
        gg=self.g
        aa=aw
        if s=='s':
            aa=sex2rad(aa)
        kres=mm*gg*sin(aa)
        return(kres)
    
    def res_incenter(self,aw=0,s='r',ktype=True):
        tt=self.Ti
        ss=s
        aa=aw
        fw=0
        if(ktype):
            fw=self.w_incenter(aw=aa, s=ss)
        kres=tt+fw
        return(kres)
    
    def forza_ac(self,s='s',ktype='f'):
        mm=self.m
        if ktype=='a':
            fx=self.resultante()
            kres=fx/mm
        
        else:
            aa=self.ac
            kres=mm*aa
        
        if s=='e':
           kres=kres.evalf()
        return(kres)
        
    def energia(self,s='s',ktype=' ',keval=False):
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
        else:
            kres=E2-E1
            
        if keval:
            mmval=self.kvalue
            kres=s_subs(kres,mmval)    
        if s=='e':
           kres=kres.evalf()
        return(kres)

    def work_fd(self, s='s'):
        Fr=self.resultante()
        d1=self.x1
        d2=self.x2
        dd=d2-d1
        kres=Fr*dd
        if s=='e':
           kres=kres.evalf()
        return(kres)
                
    def res_x(self,keval=False):
        kres=self.resultante()
        if keval:
            mmval=self.kvalue
            kres=s_subs(kres,mmval)
            
        return(kres)
    def res_y(self,keval=False):
        kres=self.resultante(ktype='y')
        if keval:
            mmval=self.kvalue
            kres=s_subs(kres,mmval)
        return(kres)
        
        
    def din_ac(self):
        mm=self.m
        ff=self.rex()
        kres=ff/mm
        return(kres)
        
    def cchange(self,kval,ktype=' '):
        if ktype=='mu':
            self.mu=kval
            sR='mu=' 
            display(Math(sR+latex(self.mu)))
        if ktype=='ac':
            self.ac=kval
            sR='aceleration=' 
            display(Math(sR+latex(self.ac)))
    def get_value(self,kkval):
        self.kvalue=kkval
    
    def Eq_eval(self,Eq1,ksimplify=False):
        mmval=self.kvalue
        kres=s_subs(Eq1,mmval)
        if ksimplify:
            kres=simplify(kres)
        return(kres)
        
    def value(self,Eq1,ksimplify=False):
        mmval=self.kvalue
        kres=s_subs(Eq1,mmval)
        if ksimplify:
            kres=simplify(kres)
        return(kres)

    def symval(self,kv):
        mm=self.kvalue
        mm1=mm[0] 
        mm2=mm[1]
        qq=len(mm1)
        for i in range(qq):
            if kv==mm1[i]:
                return(mm2[i])
        return(kv)

    def define_val(self,ksym,kval,kend=False):  
        mat1=[ksym]
        mat2=[kval]
        kkend=kend
        self.add_value(mat1,mat2,kend=kkend)
        return(kval)
    
 





        
'''
   
    def add_torque(self, kval,kang,kdis,s='r'):
        if s=='s':
            kang=sex2rad(kang)
        mm=self.Tr
        mm.append([kval,kang,kdis])
        self.Tr=mm
        
    def torque(self,s='s',kabs=False):
        mm=self.Tr 
        qq=len(mm)
        kres=0
        for i in range(qq):
            vec=mm[i]
            kv=vec[0]
            ka=vec[1]
            kd=vec[2]
            kval=kd*kv*sin(ka)
            kres+=kval

        if kabs==True:
            kres=-1*kres
        if s=='e':
           kres=kres.evalf()
        return(kres)
    
    
    def resultanteX(self,s='s',kabs=False,peso=False):
        mm=self.F 
        qq=len(mm)
        kres=0
        for i in range(qq):
            vec=mm[i]
            kk=vec[0]
            ka=vec[1]
            kval=kk*cos(ka)
            kres+=kval
            
        kmu=self.mu
        nn=self.n
        gg=self.g
        if peso:
            mm=self.m
            gg=self.g
            aa=self.a
            pp=mm*gg*cos(aa)
            kres+=pp
        if nn is not 0:
            mm=self.m
            gg=self.g
            aa=self.a
            uu=self.mu
            pp=nn*uu*mm*gg*cos(aa)
            kres+=pp
        if kabs==True:
            kres=-1*kres
        if s=='e':
           kres=kres.evalf()
        return(kres)    
    
    def resultanteY(self,s='s',kabs=False,peso=False):
        mm=self.F 
        qq=len(mm)
        kres=0
        for i in range(qq):
            vec=mm[i]
            kk=vec[0]
            ka=vec[1]
            kval=kk*sin(ka)
            kres+=kval
            
        kmu=self.mu
        nn=self.n
        gg=self.g
        if peso:
            mm=self.m
            gg=self.g
            aa=self.a
            pp=mm*gg*sin(aa)
            kres+=pp
        if nn is not 0:
            mm=self.m
            gg=self.g
            aa=self.a
            uu=self.mu
            pp=nn*uu*mm*gg*sin(aa)
            kres+=pp
        if kabs==True:
            kres=-1*kres
        if s=='e':
           kres=kres.evalf()
        return(kres)        

    def eKinetic(self,s='s',p=1,kabs=False):
        if p==1:
            vv=self.v1
        else:
            vv=self.v2
        mm=self.m
        kres=(mm*vv**2)/2
        if kabs==True:
            kres=-1*kres
        if s=='e':
           kres=kres.evalf()
        return(kres)
    
    def ePotencial(self,s='s',p=1,kabs=False):
        if p==1:
            hh=self.y1
        else:
            hh=self.y2
        gg=self.g
        mm=self.m
        kres=mm*gg*hh
        if kabs==True:
            kres=-1*kres
        if s=='e':
           kres=kres.evalf()
        return(kres)    
       
    def energia(self,s='s',p=1,kabs=False):
        kres=0
        if p==1:
            kres=self.eKinetic(p=1)+self.ePotencial(p=1)
        else:
            kres=self.eKinetic(p=2)+self.ePotencial(p=2)
        if kabs==True:
            kres=-1*kres
        if s=='e':
           kres=kres.evalf()
        return(kres)    
    
    def EnergiaTotal(self,s='s',kabs=False):
        kres=self.energia(p=2)-self.energia(p=1)
        if kabs==True:
            kres=-1*kres
        if s=='f':
           kres=simplify(kres)     
        if s=='e':
           kres=kres.evalf()
        return(kres) 
            
            
    def trabajo(self,s='s',kabs=False):
        Fx=self.resultanteX()
        L=self.x2-self.x1
        kres=Fx*L
        if kabs==True:
            kres=-1*kres
        if s=='e':
           kres=kres.evalf()
        return(kres)
    
    def Normal(self,s='s',kabs=False):
        nn=self.n
        if nn ==0:
            return(0) 
        else:
            mm=self.m
            gg=self.g
            aa=self.a
            kres=mm*gg*sin(aa)
        if kabs==True:
            kres=-1*kres
        if s=='e':
           kres=kres.evalf()
        return(kres)
    
    def Frozamiento(self,s='s',kabs=False):
        nn=self.n
        if nn ==0:
            return(0)
        else:
            uu=self.mu
            Nm=self.Normal()
            kres=uu*Nm*nn
        if kabs==True:
            kres=-1*kres
        if s=='e':
           kres=kres.evalf()
        return(kres)
    def Fuerza(self,s='s',kabs=False):
        dE=self.EnergiaTotal()
        L=self.x2-self.x1
        kres=dE/L
        if kabs==True:
            kres=-1*kres
        if s=='e':
           kres=kres.evalf()
        return(kres)
    
    def acentripeta(self,ktype='v',s='s',):
        rr=self.r
        if ktype=='w':
            ww=self.w
            kres=(ww**2)*rr
        elif ktype=='x':
            mm=self.m
            fx=self.resultante(ktype='x')
            kres=fx/mm
        elif ktype=='y':
            mm=self.m
            fy=self.resultante(ktype='y')
            kres=fy/mm    
        else:
            vv=self.v
            kres=(vv**2)/rr
        if s=='e':
           kres=kres.evalf()
        return(kres)    
            
    def angulot(self,kt,s='s'):
        vvw=self.vw
        kres=vvw*kt
        if s=='e':
           kres=kres.evalf()
        return(kres)
    
    def peso(self,s='s',ktype='xy'):
        mm=self.m
        gg=self.g
        aa=self.a
        if ktype=='y':
            kres=mm*gg*sin(aa) 
        elif ktype=='x':
            kres=mm*gg*cos(aa)
        else:
             kres=mm*gg
        if s=='e':
           kres=kres.evalf()
        return(kres)
    
    def fcentripeta(self,s='s',ktype='w'):
        ftype=ktype
        aac=self.acentripeta(ktype=ftype)
        mm=self.m
        kres=mm*aac
        if ktype=='v':
            rr=self.r
            vv=self.vt
            kres=mm*(vv**2)/rr
        if s=='e':
           kres=kres.evalf()
        return(kres)
        
    def fpeso(self,kang,s='r'):
        pp=self.peso()
        if s=='s':
            kang=sex2rad(kang)
        kres=pp*cos(kang)
        return(kres)
    
    def resultante(self,s='s',kabs=False,ktype='x'):
        mm=self.F 
        qq=len(mm)
        kresx=0
        kresy=0
        for i in range(qq):
            vec=mm[i]
            kk=vec[0]
            ka=vec[1]
            kval=kk*cos(ka)
            kresx+=kval
            kval=kk*sin(ka)
            kresy+=kval
        if ktype=='x':
            kres=kresx
        elif ktype=='y':
            kres=kresy
        else:
            kres=sqrt(kresx**2+kresy**2)
        if s=='e':
           kres=kres.evalf()
        if s=='t':
           kres=trigsimp(kres) 
        return(kres)
    
    def simple_ac(self,s='s',ktype='x'):
        if ktype=='y':
            fr=self.resultante(ktype='y')
        else:
            fr=self.resultante()
        mm=self.m
        kres=fr/mm
        if s=='e':
           kres=kres.evalf()
        return(kres)  
    
    def solve_normal(self,N1,s='s'):
        eq1=self.resultante(ktype='y')
        kres=csolve(N1-eq1,N1)
        if s=='e':
           kres=kres.evalf()
        return(kres)
        
    def creplace_M(self,k1,k2):
        mm=self.F 
        qq=len(mm)
        for i in range(qq):
            vec=mm[i]
            kk=vec[0]
            if kk==k1:
                mm[i][0]=k2
        self.F=mm
 
'''