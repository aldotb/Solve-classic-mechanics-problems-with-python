from sympy import symbols

import numpy as np
import matplotlib.pyplot as plt
from IPython.display import Image, display
from lib_Mathematica import *
from lib_Mathbasic import *
from lib_Algorith import *

from lib_Exponencial import *
from lib_Tools import *

from lib_MyEq import *
from lib_MyEqEq import *
from lib_Physics import *
from MyMatrix import *
from sympy.vector import CoordSys3D

class MyVector(mparticle):
    def __init__(self,*args,origin='',x='',y='',z='',F=1):
        N = CoordSys3D('N')
        self.F=F
        if len(args)==3:
             
            self.x1=args[0]
            self.y1=args[1]
            self.z1=args[2]
            self.V=float2int(self.x1)*N.i+float2int(self.y1)*N.j+float2int(self.z1)*N.k  
        elif len(args)==2:
            p1=x_pol(args[0],args[1])
            p2=y_pol(args[0],args[1])
            self.x1=p1
            self.y1=p2
            self.z1=0
            self.V=float2int(self.x1)*N.i+float2int(self.y1)*N.j+float2int(self.z1)*N.k  
             
        else:
            self.V=args[0]
            self.x1=self.x
            self.y1=self.y
            self.z1=self.z
        self.xx=self.x1
        self.yy=self.y1
        self.zz=self.z1        
        self.O=(0,0,0)
        if origin!='':
            self.O=origin
        if F!=1:
            R=self.V.magnitude()
            self.x1=F*self.x1/R
            self.y1=F*self.y1/R
            self.z1=F*self.z1/R
            self.V=float2int(self.x1)*N.i+float2int(self.y1)*N.j+float2int(self.z1)*N.k
    def __call__(self,*args, **kwargs):
        if len(kwargs)>0:
            kres=self.V
            kres=real_subs(kres,**kwargs)
            return kres
        return self.V
         
            
    
    def __repr__(self):
        kk=str(self.V)
        return kk

        
    def _latex(self, obj):
        return latex(self.V)    
    def __str__(self):
        return str(self.__repr__())
    def __add__(self, other):
        V=self.V
        L1=veccomp(V)
        L2=veccomp(other)
        Ls=AddList(L1,L2)
        return List2Vec(Ls)
         

    def __radd__(self, other):
        V=self.V
        L1=veccomp(V)
        L2=veccomp(other)
        Ls=AddList(L2,L1)
        return List2Vec(Ls)

    def __sub__(self, other):
        V=self.V
        L1=veccomp(V)
        L2=veccomp(other)
        Ls=SubList(L1,L2)
        return List2Vec(Ls)

    def __rsub__(self, other):
        V=self.V
        L1=veccomp(V)
        L2=veccomp(other)
        Ls=SubList(L2,L1)
        return List2Vec(Ls)
    
    def __mul__(self, other):
        p1=self.V          
        p2=other
        return p1*p2
    
    def __rmul__(self, other):
        p1=self.V          
        p2=other
        return p2*p1
        
        
    def __truediv__(self, other):
        p1=self.V          
        p2=other
        return p1/p2

    def __rtruediv__(self, other):
        p1=self.V          
        p2=other
        return p2/p1
    @property
    def x(self):
        vect=self.V
        var,val=unpack(vect.components)
        var2=[str(i) for i in var]
        try:
            value_index = var2.index('N.'+'i')
            return unpack(vect.components)[1][value_index] 
        except:
            return 0
    @property
    def y(self):
        vect=self.V
        var,val=unpack(vect.components)
        var2=[str(i) for i in var]
        try:
            value_index = var2.index('N.'+'j')
            return unpack(vect.components)[1][value_index] 
        except:
            return 0
    @property
    def z(self):
        vect=self.V
        var,val=unpack(vect.components)
        var2=[str(i) for i in var]
        try:
            value_index = var2.index('N.'+'k')
            return unpack(vect.components)[1][value_index] 
        except:
            return 0
    @property
    def vx(self):
        N = CoordSys3D('N')
        return self.x1*N.i
    @property
    def vy(self):
        N = CoordSys3D('N')
        return self.y1*N.j
    @property
    def vz(self):
        N = CoordSys3D('N')
        return self.z1*N.k
        
    @property
    def module(self):
        Vec=self.V
        return Vec.magnitude()
    
    @property
    def unitary(self):
        vec=self.V
        var,val=unpack(vec.components)
        kres=0
        for i in val:
            kres=kres+i**2
        kres= rsimplify(simplify(expand(sqrt(kres))))
        return vec/kres
    def one(self):
        vec=self.V
        var,val=unpack(vec.components)
        kres=0
        for i in val:
            kres=kres+i**2
        kres= rsimplify(simplify(expand(sqrt(kres))))
        return vec/kres
    def Add(self,vec):
        V=self.V
        L1=veccomp(V)
        L2=veccomp(other)
        Ls=AddList(L1,L2)
        kres=List2Vec(Ls)
        self.V=kres
        self.x1,self.y1,self.z1=veccomp(kres)
        return kres
    def veccomp(self):
        kres=self.V
        return veccomp(kres)
    
    @property
    def axes(self):
        vecc=self.veccomp()
        xyz=(x,y,z)
        kres=[]
        for i,j in zip(xyz,vecc):
            if j!=0:
                kres.append(i)
            else:
                kres.append(0)
        return kres
    @property
    def saxes(self):
        kres=self.axes 
        kres2=''
        for i in kres:
            if i!=0:
                kres2=kres2+str(i)
        return kres2        
    def torque(self):
        vec=self.veccomp()
        pos=self.O
        kres=0
        for i,j in zip(vec,pos):
            kres=kres+i*j
        return kres
    def axes_torque(self):
        vec=self.veccomp()
        vec2=self.O
        kres=[]
        vec3=[x,y,z]
        for i,j,k in zip(vec,vec2,vec3) :
            if i+j==0:
                kres.append(k)
            else:
                kres.append(0)
        return kres 
    def moment(self):
        N = CoordSys3D('N')
        Uv=[N.i,N.j,N.k]
        H=MyMat(Uv,self.MatrixO(),self.Matrix(),kshow=False)
        return H.D
    def R(self):
        vecr=MyVector(*self.O)
        return vecr
    def Rvec(self):
        N = CoordSys3D('N')
        O=self.O
        return [O[0]*N.i,O[1]*N.j,O[2]*N.k]
    def Rmodule(self):
        N = CoordSys3D('N')
        O=self.O
        vec=MyVector(O[0],O[1],O[2])
        return vec.module
    def Matrix(self):
        return vec2matx(self.V,op='')
    def MatrixO(self):
        return vec2matx(self.O,op='')
    def s(self):
        sexpr=0
        p1=self.x1
        sexpr=sexpr+p1*i
        p2=self.y1
        sexpr=sexpr+p2*j
        p3=self.z1
        sexpr=sexpr+p3*k
        display(Math(latex(sexpr)))

def veccomp(expr):
    if type(expr)==MyVector:
        return expr.x,expr.y,expr.z
    else:
        vec=MyVector(expr)
        return vec.x,vec.y,vec.z
def obj2MyVec(obj):
    if type(obj)==MyVector:
        return obj
    else:
        return MyVector(obj)
def MyVec2obj(obj):
    if type(obj)==MyVector:
        return obj.V
    else:
        return obj
    
def AddList(L1,L2):
    vec=[]
    for i,j in zip(L1,L2):
        vec.append(i+j)
    return vec
def SubList(L1,L2):
    vec=[]
    for i,j in zip(L1,L2):
        vec.append(i-j)
    return vec
def List2MyVec(L):
    return MyVector(*L)
def List2Vec(L):
    vec=List2MyVec(L)
    return vec.V 

def vec2matx(*args,op=''):
     
    if len(args)==1:
        vec=args[0]
        if type(vec)==MyVector:
            x1,y1,z1=vec.x,vec.y,vec.z
        elif type(vec)==tuple:
            vec2=MyVector(*vec)
            x1,y1,z1=vec2.x,vec2.y,vec2.z
        else:
            vec2=MyVector(vec)
            x1,y1,z1=vec2.x,vec2.y,vec2.z
         
        if op=='T':
            kres= MyMat(x1,y1,z1,3,1,kshow=False)
        else:
            kres= MyMat(x1,y1,z1,1,3,kshow=False)
        
        return kres.M
        
        
def find_unitary(vec):
    disp('vec.one()')
    x,y,z=symbols('x y z')
    disp('if \; vec= \;xi\;+yj\;+zk \; then:')
    R=sqrt(x*x+y*y+z*z)
    disp('R=',R)
     
    disp('VecUnitary \; is :', x/sqrt(x*x+y*y+z*z),'+',y/sqrt(x*x+y*y+z*z),'+',z/sqrt(x*x+y*y+z*z))
    x,y,z=veccomp(vec)
    disp('finding  \;   components \; i,j,k \;of \;vec')
    disp('x=',x,',y=',y,',z=',z)
    R=sqrt(x*x+y*y+z*z)
    disp('R=',R)
    disp('componets \;unitary \;are \;: x=',x/R,',y=',y/R,',z=',z/R)
    disp(vec.one())
    disp('type \; Vec.one()')


def find_resultant(*args):
    ss='\;'
    disp(' find V1+V2+V3 ')
    x,y,z=symbols('x y z')
    disp('for each vector include in args: we will separate by type axes i, j, k' )
    disp('Rx, Ry, Rz = 0,0,0')
    Rx,Ry,Rz=0,0,0
    disp('for i in args:')
    for i in args:
        disp(ss*6,'x ,y ,z =',i.veccomp())
        x,y,z=i.veccomp()
        Rx=Rx+x
        Ry=Ry+y
        Rz=Rz+z
        disp(ss*6,'Rx ,Ry ,Rz =',Rx,',',ss,Ry,',',ss,Rz)
        disp('next..') 
    disp(  'kres = MyVector(Rx,Ry,Rz)')  
    kres=MyVector(Rx,Ry,Rz)
    disp(kres) 


from IPython.display import Image
from Libaldomath import *
from lib_MyVec import *
from MyMatrix import *
sys.path.insert(0, 'imgHelp/')
def find_moment(vec):
    ss='\;'
    disp('whith this sample Vector')
    if vec.F==1:
        if vec.O=='': 
            disp(ss*6,'V= MyVector(',vec.xx,',',vec.yy,',',vec.zz,')')
        else:
            disp(ss*6,'V= MyVector(',vec.xx,',',vec.yy,',',vec.zz,'origin=',vec.O,')')
    else:
        if vec.O=='': 
            disp(ss*6,'V= MyVector(',vec.xx,',',vec.yy,',',vec.zz,'F=',vec.F,')')
        else:
            disp(ss*6,'V= MyVector(',vec.xx,',',vec.yy,',',vec.zz,'origin=',vec.O,'F=',vec.F,')')
    Image(filename = "Libaldo/torquevector.png")
    display(Image(filename = "Libaldo/torquevector.png"))
    Uvec=vec.one()
    direction=MyVector(Uvec)
    x1,y1,z1=direction.veccomp()
    disp('direction =',vec.xx,ss,',',vec.yy,ss,',',vec.zz)
    disp("N = CoordSys3D('N')")
    N = CoordSys3D('N')
    
    disp("Uv = [N.i,N.j,N.k]")
    Uv=[N.i,N.j,N.k]
    disp(' creatibg matrix to find determinant')
    disp('H = MyMat(Uv,vec.MatrixO() ,vec.Matrix())')
    
    H=MyMat(Uv,vec.MatrixO(),vec.Matrix(),kshow=False)
    disp( ' result = H.D')
    disp(H.D)
    disp('call: vec.moment()')    