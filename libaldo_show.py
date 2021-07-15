# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 07:35:49 2021

@author: LibaldoHouse
"""

from sympy import *
from IPython.display import display, Math


def multishow(vecs):
    for eeq in vecs:
        show_eq(eeq)
        
def show_eq(Eq1):
    kres=Eq1
    display(Math(latex(kres)))
    
def show_text(kd=''):
    sR=kd
    display(Math(sR))

def show_res(Eq1,kd='',kr=''):
    kres=Eq1
    sR=kd+' ='
    display(Math(sR+latex(kres)))
    if kr!='':
        return(kres)
        
def show_if(Eq1):
    kres=Eq1
    sR='if  \; \; \; \;' 
    display(Math(sR+latex(kres)))
def show_ifthen(Eq1):
    kres=Eq1
    sR1='if  \; \; \; \;'
    sR2='\; \; \; \; then'   
    display(Math(sR1+latex(kres)+sR2))
       
def fraseconvert(lst):
    return (lst[0].split())
def show_expP(vecexp):
    kexp=''
    for i in vecexp:
        if type(i)==str:
            kexp=kexp+ i + '\;'
        else:
            kexp=kexp+latex(i)+'\;'
    display(Math(kexp))
    
def show_exp(vecexp):
    kexp=''
    for i in vecexp:
        if type(i)==str:
            vecfrase=fraseconvert([i])
            for j in vecfrase:
                kexp=kexp+j + '\;'
        else:
            kexp=kexp+latex(i)+'\;'
    display(Math(kexp))     
def sE(vecexp): # short funcion of show_exp
    show_exp(vecexp)
    
def show_sval(ksym):
    vsym=str(ksym)+'='
    show_exp([vsym,ksym])
    
def define_eQ(kstr,ksym):
    kexp=[]
    for vs,vv in zip(kstr,ksym):
        kexp.append(vs+'=')
        kexp.append(vv)
        kexp.append(', ')
        kexp.append('\; \; \;')
        kexp.append(' ')
    show_exp(kexp)
    return(ksym)    
        
        
        
def show_eques(vecvec):
    kexp=''
    for i in vecvec:
        kexp=str(i)+' \;'+'='+' \;'+latex(i)+', \;'
    display(Math(kexp))
    

def kdisplay(skres):
        display(Math(latex(skres)))
        

def supersotore(kclass,kvar,kval):
    mm=kclass.kvalue
    

def csubs(Eq1,ks,kv,kd=' ',s='s'):
    kres=Eq1.subs(ks,kv)
    sR=kd+' ='
    if s=='e':
        kres=kres.evalf()
    if s=='t':
        kres=trigsimp(kres)
    display(Math(sR+latex(kres)))    
    return(kres)



def super_subs(Eq1,tup1,tup2):
    qq=len(tup2)
    for i in range(qq):
        ss=tup1[i]
        sv=tup2[i]
        Eq2=Eq1.subs(ss,sv)
        Eq1=Eq2
    return(Eq1)

def s_subs(Eq1,keval):
    tup1=keval[0]
    tup2=keval[1]
    kres=super_subs(Eq1,tup1,tup2)

    return(kres)     
