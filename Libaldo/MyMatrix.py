from sympy import *
import copy
from libaldo_math2 import *
from libaldo_algorith import *
from mathbasic import *
from mathexponencial import *
from lib_tools import *
 

a00,a01,a02,a03,a04,a05=symbols('a00 a01 a02 a03 a04 a05')
a10,a11,a12,a13,a14,a15=symbols('a10 a11 a12 a13 a14 a15')
a20,a21,a22,a23,a24,a25=symbols('a20 a21 a22 a23 a24 a25')
a30,a31,a32,a33,a34,a35=symbols('a30 a31 a32 a33 a34 a35')
a40,a41,a42,a43,a44,a45=symbols('a40 a41 a42 a43 a44 a45')
a50,a51,a52,a53,a54,a55=symbols('a50 a51 a52 a53 a54 a55')

b00,b01,b02,b03,b04,b05=symbols('b00 b01 b02 b03 b04 b05')
b10,b11,b12,b13,b14,b15=symbols('b10 b11 b12 b13 b14 b15')
b20,b21,b22,b23,b24,b25=symbols('b20 b21 b22 b23 b24 b25')
b30,b31,b32,b33,b34,b35=symbols('b30 b31 b32 b33 b34 b35')
b40,b41,b42,b43,b44,b45=symbols('b40 b41 b42 b43 b44 b45')
b50,b51,b52,b53,b54,b55=symbols('b50 b51 b52 b53 b54 b55')

c00,c01,c02,c03,c04,c05=symbols('c00 c01 c02 c03 c04 c05')
c10,c11,c12,c13,c14,c15=symbols('c10 c11 c12 c13 c14 c15')
c20,c21,c22,c23,c24,c25=symbols('c20 c21 c22 c23 c24 c25')
c30,c31,c32,c33,c34,c35=symbols('c30 c31 c32 c33 c34 c35')
c40,c41,c42,c43,c44,c45=symbols('c40 c41 c42 c43 c44 c45')
c50,c51,c52,c53,c54,c55=symbols('c50 c51 c52 c53 c54 c55')

MatMasterA=Matrix([(a00,a01,a02,a03,a04,a05),(a10,a11,a12,a13,a14,a15),(a20,a21,a22,a23,a24,a25),(a30,a31,a32,a33,a34,a35),(a40,a41,a42,a43,a44,a45),(a50,a51,a52,a53,a54,a55)]) 
MatMasterB=Matrix([(b00,b01,b02,b03,b04,b05),(b10,b11,b12,b13,b14,b15),(b20,b21,b22,b23,b24,b25),(b30,b31,b32,b33,b34,b35),(b40,b41,b42,b43,b44,b45),(b50,b51,b52,b53,b54,b55)]) 
MatMasterC=Matrix([(c00,c01,c02,c03,c04,c05),(c10,c11,c12,c13,c14,c15),(c20,c21,c22,c23,c24,c25),(c30,c31,c32,c33,c34,c35),(c40,c41,c42,c43,c44,c45),(c50,c51,c52,c53,c54,c55)]) 

e1,e2,e3,e4,e5,e6=symbols('e1 e2 e3 e4 e5 e6')
vecee=[e1,e2,e3,e4,e5,e6]
def getdataMat(expr):
    if type(expr)==MyMat:
        return expr.M
    else:
        return expr
    
def prematrix(*args):
    
    cc=2
    M=[]
    for i in range(args[0]):
        MM=[]
        for j in range(args[1]):
            MM.append(args[cc])
            cc=cc+1
        M.append(MM)
    kres=M
    print(M)
                    
                
   
class MyMat:
    def __init__(self, *args,kshow=True):
        try: 
            self.type='M'
            qq=len(args)
            kcol=args[qq-1]
            krow=args[qq-2]
             
            mm=[]
            cc=0
            for i in range(krow):
                smat=[]
                for j in range(kcol):
                    smat.append(args[cc])
                    cc=cc+1
                  
                mm.append(smat)
             
            M=Matrix(mm)
        except:
            expr=args[0] 
            if len(args)==1:
                expr=args[0]
                if type(expr)==MyMat:
                    M=expr.M
                elif type(expr)==Matrix:
                    M=expr
            elif type(expr)==str:
                op=expr
                
                if op=='zero' or op=='zeros' or op=='ceros' or op=='cero':
                    p1=args[1]

                    if len(args)==3:

                        p2=args[2]
                        M=zeros(p1,p2)

                    else:
                        M=zeros(p1)
                if op=='eye' or op=='eyes':
                    
                    p1=args[1]
                    M=eye(p1)
                if op=='one' or op=='ones':

                    p1=args[1]

                    if len(args)==3:

                        p2=args[2]
                        M=ones(p1,p2)

                    else:
                        M=ones(p1)
                    
            else:    
                mm=[]
                for i in args:
                    mm.append(i)
                M=Matrix(mm)
            
        self.M=M
        if kshow:     
            display(self.M)
    def __call__(self,*args,**kwargs):
        if len(args)==0 and len(kwargs)==0:
            return  self.M
        if len(args)>0:
            p1=args[0]
            p2=args[1]
            M=self.M
            return M[p1,p2]

    ## Basisc ope
    def __add__(self, other):        
        M=self.M
        p2=getdataMat(other)
        return M + p2
        
    def __radd__(self, other):
        M=self.M
        p2=getdataMat(other)
        return p2 + M

    def __sub__(self, other):
        M=self.M
        p2=getdataMat(other)
        return M - p2
        

    def __rsub__(self, other):
        M=self.M
        p2=getdataMat(other)
        return p2 - M
        
    def __mul__(self, other):
        """ Returns the vector addition of self and other """
        M=self.M
        p2=getdataMat(other)
        return M * p2
         
    def __rmul__(self, other):
        """ Returns the vector addition of self and other """
        M=self.M
        p2=getdataMat(other)
        return p2 * M
         
    def s(self):
        M=self.M
        display(M)
    def row(self,*args):
        M=self.M
        p1=args[0]
        if len(args)==1:
            return M[:,p1]
        else:
            p2=args[1]
            if p2>self.qrow:
                p2=self.qrow
            return M[:,p1:p2+1]
            
    def col(self,*args):
        M=self.M
        p1=args[0]
        if len(args)==1:
            return M[p1,:]
        else:
            p2=args[1]
            if p2>self.qcol:
                p2=self.qcol
            return M[p1:p2+1,:]
    
    def submatrix(self,p1,p2):
        return self.subMat(p1,p2)

        
    def subMat(self,p1,p2):
        M=self.M
        if type(p1)!=list:
            if type(p2)!=list:
                kres=M[p1,p2]
            else:
                p21,p22=p2
                kres=M[p1,p21:p22+1]
        else:
            p11,p12=p1
            if type(p2)!=list:
                kres=M[p11:p12+1,p2]
            else:
                p21,p22=p2
                kres=M[p11:p12,P21:p22+1]
        return kres
        
    def complement(self,p1,p2):
        M=self.M
        return M.minor_submatrix(p1,p2)
        
    def subdet(self,p1,p2):
        kres =self.submatrix(p1,p2) 
        return kres.det()
       
    # propeties
    @property
    def qcol(self):
        M=self.M
        qc,qr=M.shape
        return qr
    @property
    def qrow(self):
        M=self.M
        qc,qr=M.shape
        return qc
    
    @property
    def shape(self):
        M=self.M
        return M.shape
    
    def row(self,yy):
        M=self.M
        return M[yy,:]
    def col(self,xx):
        M=self.M
        return M[:,xx]    
        
    ## Transformation
    def set(self,*args,**kwargs):
        if len(args)==3:
            M=self.M
             
            M[args[0],args[1]]=simplify(args[2])
            
            self.s()
        if len(kwargs)>0:
            qcol=self.qcol
            qrow=self.qrow
            M2=zeros(qrow,qcol)
             
            M=self.M
            for i in range(qrow):
                for j in range(qcol):
                    if 'noeval' in args:
                        M2[i,j]=symbol_subs(M[i,j],**kwargs)
                         
                    else:    
                        M2[i,j]=simplify(unisymbols(real_subs(M[i,j],**kwargs)))
            self.M=M2
            self.s()

        
    # Transpose
    @property
    def T(self):
        kres=self.M
        return kres.T
    
    # Inverse
    @property
    def I(self):
        M=self.M
        try:
            MI= M.inv()
        except:
            MI=0
        return MI
    # Determinant
    @property
    def D(self):
        M=self.M
        return M.det()
    
    def Abs(self):
        qcol=self.qcol
        qrow=self.qrow
        M2=zeros(qrow,qcol)
         
        M=self.M
        M2d=M2.M
        for i in range(qrow):
            for j in range(qcol):
                M2d[i,j]=abs(M[i,j])
            M2.M=M2d    
        return M2
    def pow(self,expr):
        M=self.M
        return M**expr
        
    def Pow(self,expr):
        M=self.M
        self.M = M**expr
        self.s()
        
    def simplesolve(self,*args):
        M=self.M
        Mdata= M[:]
        if len(args)>0:
            for i in args:
                Mdata.append(i)
        return simplesolve(*Mdata)
    def xcopy(self):
        return copy.deepcopy(self)

    def add(self,expr):
        M=self.M
        H,W=M.shape
        M2=mconstant(expr,H,W)
        M3=M+M2
        return M3
        
     
    def substrac(self,expr):
        M=self.M
        H,W=M.shape
        M2=mconstant(expr,H,W)
        M3=M-M2
        return M3
        
    def Add(self,expr):
        M3=self.add(expr)
        self.M=M3 
        self.s()
        
    def Substrac(self,expr):
        M3=self.substrac(expr)
        self.M=M3 
        self.s() 
        
    def mul(self,expr):
        M=self.M
        M2=M*expr
        return M2
    
    def Mul(self,expr):
        M=self.M
        M2=M*expr
        self.M=M2
        self.s()
        
    def div(self,expr):
        M=self.M
        M2=M/expr
        return M2
    
    def Div(self,expr):
        M=self.M
        M2=M/expr
        self.M=M2
        self.s()
    def pow(self,expr):
        M=self.M
        M2=M**expr
        return M2  
    
    def Pow(self,expr):
        M=self.M
        M2=M**expr
        self.M=M2
        self.s()

    def sum(self):
        M=self.M
        return sum(M)
    
    def vecdata(self):
        M=self.M
        return M[:]
        
    def getfactor(self):
        vecM=self.vecdata()
        kres=gcd(vecM)
        return kres
    # sympy
    def simplify(self):
        M=self.M
        return simplify(M)
    def expand(self):
        M=self.M
        return expand(M)
    def factor(self):
        M=self.M
        return factor(M)

    
    def tsimplify(self):
        M=self.M
        return trigsimp(M)
         
    def evalif(self,sexp1,sexp2):
        '''
        self.evalif('3','i=j') change self.M[i,j]=3 if i=j
        '''
        M=self.M
        M2=self.M
        H,W=M.shape
        for i in range(H):
            for j in range(W):
                if eval(sexp2):
                    M2[i,j]=eval(sexp1)
        return M2
        
    def setif(self,sexp1,sexp2):
        M=self.M
        M2=self.M
        H,W=M.shape
        for i in range(H):
            for j in range(W):
                if eval(sexp2):
                    M2[i,j]=eval(sexp1)
        self.M=M2
        self.s()
    def Adj(self):
        expr=self.M
        H,W=expr.shape
        M2=copy.deepcopy(expr)
        for i in range(H):
            for j in range(W):
                if i!=j:
                    M2[i,j]=expr[j,i]
        return M2

    def traz(self):
        M=self.M
        H,W=M.shape
        kres=0
        for i in range(H):
            kres=kres+M[i,i]
        return kres 
    def Hstack(self,expr):
        if type(expr)==MyMat:
            expr=expr.M
        M=self.M
        M2=Matrix.hstack(M,expr)
        self.M=M2
        self.s()
        
    def Vstack(self,expr):
        if type(expr)==MyMat:
            expr=expr.M
        M=self.M
        M2=Matrix.vstack(M,expr)
        self.M=M2
        self.s()    
            
def deter(expr):
    if type(expr)==MyMat:
        return expr.D 
    else:
        return expr.det()

def mconstant(expr,H,W):
    M=zeros(H,W)
    for i in range(H):
        for j in range(W):
            M[i,j]=expr
    return M

def getlist(mmat):
        if type(mmat)==MyMat:
            mmat=mmat.M
        yy,xx=mmat.shape
         
        kres=[]
        for i in range(yy):
            for j in range(xx):
                kres.append(mmat[i,j])
        display(kres)        
        return kres
    
def MatSolveSys(xvar,xval):
    if type(xvar)==MyMat:
        xvar=xvar.M        
    if type(xval)==MyMat:
        xval=xval.M
    kres=xvar.LUsolve(xval)
    kres2=kres.T
    val=kres2[:]
     
    ksolu=getlist(kres2)

    return ksolu

def prematrix(*args):
    
    cc=2
    M=[]
    for i in range(args[0]):
        MM=[]
        for j in range(args[1]):
            MM.append(args[cc])
            cc=cc+1
        M.append(MM)
    kres=M
    sres=str(M)
    sres=sres.replace('[','(')
    sres=sres.replace(']',')')
    sres="=MyMat"+sres
    return sres


def trans_simetric(mmat,var='a'):
    yy,xx=mmat.shape
    for i in range(1,yy):
        for j in range(0,i):
            mmat[j,i]=mmat[i,j]
    m2=createMatrix(yy,xx,var,kshow=False)
    for i in range(yy):
        mmat[i,i]=m2(i,i)
        
    return mmat
def trans_asimetric(mmat,var='a'):
    yy,xx=mmat.shape
    for i in range(1,yy):
        for j in range(0,i):
            mmat[j,i]=-1*mmat[i,j]
    for i in range(yy):
        mmat[i,i]=0
    m2=createMatrix(yy,xx,var,kshow=False)
    for i in range(yy):
        mmat[i,i]=0    
    return mmat


def createMatrix(*args,kshow=True):
    nrow=args[0]
    ncol=args[1]
    vecop=['b','c','simetric','asimetric']
    MM=copy.deepcopy(MatMasterA)
    var='a'
    if 'b' in args:
        MM=copy.deepcopy(MatMasterB)
        var='b'
    if 'c' in args:
        MM=copy.deepcopy(MatMasterC)
        var='c'

    if 'simetric' in args  or 'symmetric' in args or 'symetric'  in args: 
        MM=trans_simetric(MM,var)
    if 'asimetric' in args  or 'asymmetric' in args or 'asymetric'  in args: 
        MM=trans_asimetric(MM,var)    
    newm=MM[0:nrow,0:ncol]
    return MyMat(newm,kshow=kshow)




 
from vectors import Vector
x,x1,x2,x3,y,y1,y2,y3,z,z1,z2,z3,i,j,k=symbols('x x1 x2 x3 y y1  y2 y3 z z1 z2 z3 i j k')
Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz=symbols('Ax Ay Az Bx By Bz Cx Cy Cz')

class MyVec(MyMat):
    def __init__(self, *args,kshow=True):
  
        mm=[]
        self.type='V'
        if len(args)==1:
            self.P = Matrix(args[0])
            self.M = Matrix(args[0])    
            self.O = Matrix([0,0,0])
            
        elif len(args)==2 and type(args[0])!=list:
             
            self.O = Matrix(args[0])
            self.P = Matrix(args[1])
            mm=args[1]
            mm2=[]
            for i,j in zip(mm,self.O):
                mm2.append(i-j)
                
                
            self.M = Matrix(mm2)
             
             
        else:    
            for i in args:
                mm.append([i])
            self.P =  Matrix([mm[0][0],mm[1][0],mm[2][0]])
            self.M = Matrix(mm)    
            self.O = Matrix([0,0,0])
             
        x1,y1,z1=self.M
        self.Vector=Vector( x1,y1,z1)   
        if kshow:     
            display(self.Vector)
    def __call__(self,*args,**kwargs):
        if len(args)==0 and len(kwargs)==0:
            return  self.M
    def __abs__(self):
        return self.module

    def __add__(self, u):
        if isinstance(u, MyVec):
            return  self.M + u.M
        else:
            return self.M + u
        
    def __radd__(self, u):
        if isinstance(u, MyVec):
            return  u.M + self.M  
        else:
            return u+ self.M  
    

    def __mul__(self, u):
        if isinstance(u, MyVec):
            return self.M*u.M
        else:
            return  self.M*u
    def __rmul__(self, u):
        if isinstance(u, MyVec):
            return u.M*self.M
        else:
            return u*self.M
    
    def __or__(self, u):
        return acos(self * u / (self.module * u.module))

    def __pow__(self, u):
        return self.cross(u)

     
    def __repr__(self):
        return str(Matrix((self.vector.coeff(i),self.vector.coeff(j),self.vector.coeff(z))))
        '''
                      return 'Vector(' + str(self.vector.coeff(i)) + ', ' \
            + str(self.vector.coeff(j)) + ', ' \
            + str(self.vector.coeff(k)) + ')'
        '''
    def __str__(self):
        return str(self.M)

    def __sub__(self, u):
        if isinstance(u, MyVec):
            return  self.M - u.M
        else:
            return self.M - u
    def __rsub__(self, u):
        if isinstance(u,MyVec):
            return  u.M - self.M  
        else:
            return u - self.M    

    def __xor__(self, u):
        return asin((self ** u).module / (self.module * u.module))   
    
    def s(self):
        M=self.M
        O=self.O
        self.P=M+O
        display(M)
    def angle(self,V1):
        v1=self.Vector
        v2=V1.Vector
        return v1|v2
    
    @property
    def module(self):
        return self.Vector.module
    
    
    def dot(self,v2):
        M=self.M
        if isinstance(v2, MyVec):
            return  M.dot(v2.M)
        elif type(v2)==Matrix:
            return M.dot(v2)
        else:
            p2=Matrix(v2)
            return M.dot(p2)
        
        
    
    def cross(self,v2):
         
        x1=self.x
        y1=self.y
        z1=self.z
        if isinstance(v2, MyVec):
            x2=v2.x
            y2=v2.y
            z2=v2.z
        else:
            x2=v2[0]
            y2=v2[1]
            z2=v2[2]

        kres=kres=([-x2*z1 + y1*z2] ,[-x1*z2 + x2*z1], [x1*y2 - x2*y1])

        kres= Matrix(kres)
        return kres.T
        
        return self.Vector.cross(MyVec2Vec(v2))
    @property
    def x(self):
        return self.M[0]
    @property
    def y(self):
        return self.M[1]
    @property
    def z(self):
        return self.M[2]
    @property
    def T(self):
        return self.M.T
    @property
    def List(self):
        return list(self.M.T)
    @property
    def anglex(self):
        return acos(self.x/self.module)
    @property
    def xangle(self):
        return acos(self.x/self.module)
        
    @property
    def angley(self):
        return acos(self.y/self.module)
    
    @property
    def yangle(self):
        return acos(self.y/self.module)
    
    @property
    def anglez(self):
        return acos(self.z/self.module)

    @property
    def zangle(self):
        return acos(self.z/self.module)
        
    def area_paralell(self,V2):
        kres=self.cross(V2)
        return kres.module
    
    def Point3D(self):
        return Point3D(self.P)

def createVector3D(var='A'):
    s1=var+'x'
    s2=var+'y'
    s3=var+'z'
    s1,s2,s3=symbols(s1+' '+s2+' '+s3)
    return MyVec(s1,s2,s3)

    