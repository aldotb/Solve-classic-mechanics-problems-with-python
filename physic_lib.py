from sympy import *
import traceback
import PySimpleGUI as sg 
from libaldo_math import *
from polyclass  import *
import libaldo_show as ae
from lib_algebraEq import *
from IPython.display import display, Math

m,g,x1,x2,y1,y2,v1,v2,g, t, a,ac,v,r,vx,vy,vxy,w,Po=symbols('m g x1 x2 y1 y2 v1 v2 g t a  ac v r vx vy vxy w Po')
 

class mparticle:
    def __init__(self,x1=x1,x2=x2,y1=y1,y2=y2,v1=v1,v2=v2,m=m,a=a,g=g,v=v,ac=ac,s='r',t=t,r=r,vx=vx,vy=vy,vxy=vxy,w=w):
                 
        self.x1=x1 # posicion 1
        self.y1=y1 # altura 1
        self.x2=x2 # posicion 2
        self.y2=y2 # altura 2
        self.v1=v1 # velocidad en 1
        self.v2=v2 # velocidad en 2
        self.m=m   # masa 
         
        self.a=sex2rad_i(a,s) # angulo de direccion
        self.F=[]  #matriz de fuerzas
        self.g=g   # gravedad
        
        self.v=v   # velocidad tangencial o velocidad inicial cinemat
        self.cG=[] # Centro de Gravedad
        self.Fc=[] # Fuerzas circulares
        self.r=r 
        self.ac=ac
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
# ***************************************** 
# ********   Manage Variables    **********
# *****************************************
    def add_forza(self, kval,kang,x=0,y=0,s='r'): 
        
        if s=='s':
            kang=sex2rad(kang)
        mm=self.F
        qq=len(mm)
        if qq>0:
            veri=[]
            for i in range(qq):
                veri.append(mm[i][0])
                 
            mm.append([kval,kang,x,y])
        else: 
            mm.append([kval,kang,x,y])        
        self.F=mm
        
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
    
           
    def setValue(self,nomv,valv,kdisp=False,kpositive=False,kope='',kunisym=True,kreturn=False,krevalue=True):
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
        
    def evalue_val(self,kv, kope=''):
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
    
    
         
# ********   show()    ********** 

    def show_store(self,kunic=False):
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
  
  
# ***************************************** 
# ********   Static Result sum F    **********
# *****************************************    
    
    def x_res(self,kope='',keval=True,kunisym=True):
        kres=self.resultante(kope=kope)
        
        if kunisym:
            kres=unisymbols(kres)
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope=kope)
        return(kres)
     
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
    
    def tan_res(self):
        
        xres=self.x_res()
        yres=self.y_res()
        p1=Point(0, 0)
        p2=Point(self.x_res(), self.y_res())
        L = Line(p1, p2)
        return(L.slope)
        
    def ang_res(self):
        tanr=self.tan_res()
        kres=atan(tanr)
        return(kres)    
        
    def vec_displace(self):
        xx=self.x2-self.x1
        yy=self.y2-self.y1
        return([xx,yy])
        
    def displacement(self):
        xx=self.x2-self.x1
        yy=self.y2-self.y1
        kres=rpow(xx*xx+yy*yy)   
        return(kres)
        
    def resultante(self,keval=True,ktype='x',kope=''): # usado x x res y res
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
            
# ********   Torque    **********

    def torque(self,x=0,y=0,kope='',keval=True):
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
        if keval:
            kres=self.revalue(kres)        
        kres=opemat(kres,kope=kope)
        return(kres) 
        
        
# **************************************** 
# ********    SIMPLE FUNC      ***********
# ****************************************

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
        
        return(unisymbols(kres))
        
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
     

# ********   C Gravedad    ********** 

    def add_cg(self, km,x=0,y=0):
        
        mm=self.cG
        mm.append([km,x,y])
        self.cG=mm
        
    def cgravedad(self,ktype='',kope='',kshow=False):
        mm=self.cG 
        Tm=0
        Tx=0
        Ty=0

        for i in mm:
            Tm+=i[0]
            Tx+=i[0]*i[1]
            Ty+=i[0]*i[2]
        kres= Ty/Tm
        kresy=opemat(kres,kope)
        kres= Tx/Tm
        kresx=opemat(kres,kope)
        
        if ktype=='y':
            return kresy
        elif ktype=='x':
            return kresx
        else:
            if kshow:
                show_res(kresx,'Cx')
                show_res(kresy,'Cy')
            return [kresx, kresy]

# ********    Energy work      ***********

    def energia(self,ktype=' ',keval=False,kope=''):
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
            mmval=self.kvalue
            kres=s_subs(kres,mmval)    
        kres=opemat(kres,kope)
        return(kres)
                
    def work_x(self, kope=''):
        Fr=self.resultante('x')
        d1=self.x1
        d2=self.x2
        dd=d2-d1
        kres=Fr*dd
        kres=opemat(kres,kope)
        return(kres)       
        
# ********    kinematic      ***********

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
        tt=t
        vvx=self.x_vel(t=0)
        aac=self.ac
        xx1=self.x1
        kres=xx1+vvx*tt+aac*tt*tt/2
        kres=opemat(kres,kope=kope)
        if keval:
            kres=self.revalue(kres)
        return(kres)
    
    def y_pos(self, t=t,kope='',keval=True):
        tt=t
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
        gg=self.g
        aa=self.a
        vv=self.v
        kres=(vv**2)*(sin(2*aa))/gg
        kres=opemat(kres,'x')
        kres=opemat(kres,kope=kope)
        if keval:
            kres=self.revalue(kres)
        return(kres)    

    def t_fly(self,kope='',keval=True):
        xx=self.x_max()
        vx=self.x_vel()
        kres=xx/vx

        kres=opemat(kres,kope=kope)
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
        
        
# ********   Dynamic circular      ***********

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
            
    def fcentripeta(self,kangle='',ktype='', keval=True,kope=''):
        if kangle=='':
            kres=self.Ti
            if keval:
                kres=self.revalue(kres)
            
        else:
            kres=self.Ti
            kres=kres.subs(alpha,kangle)
        if ktype=='w':
            ww=self.w
            rr=self.r
            kwr=ww*rr
            kres=kres.subs(v,kwr)
            
# **************************************** 
# *******  CIRCULAR DYNAMIC    *********** 
# ****************************************              
    
     
    def c_res(self,aac=0):
        xr=self.x_res()
        yr=self.y_res()
        
        ffc=-yr*sin(aac)+xr*cos(aac)
        ffc=self.revalue(ffc)
        fft=-yr*cos(aac)+xr*sin(aac)
        fft=self.revalue(fft)
        return sE(['ffc=',ffc,' fft=',fft])
    
    def fcen_res(self,aac=0):  #return centripete force in pos aac angle due forces added
        xr=self.x_res()
        yr=self.y_res()
        ffc=-yr*sin(aac)+xr*cos(aac)
        ffc=self.revalue(ffc)
        return  ffc 
        
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
            

# **************************************** 
# *******          Power      *********** 
# ****************************************             
 
    def simplePow(self,t=t,kope='',keval=True,ktype=''):
        tt=t
        ww=self.energia(ktype)
        kres=ww/tt
        if keval:
            kres=self.revalue(kres)
        kres=opemat(kres,kope)
        return kres
        
             