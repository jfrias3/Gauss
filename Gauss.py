# -*- coding: utf-8 -*-
from functools import *
import matplotlib.pyplot as plt
from random import randint
from math import sqrt



def even(n):
    f = lambda x : x/2 == x//2
    
    return f(n)
    
def coll(n):#image under "collatz map"
    if (even(n)):
        n = n//2
        
    else:
        n = 3*n+1
    return n

def Primes(n):#Only positive
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)
    return factors

def qres(n):#quadratic residues
    qres = []
    i = 1
    while (i<n):
        quad = lambda x : (x**2)%n == i%n
        if (len(list(filter(quad,range(1,n+1)))) >0):
            qres.append(i)
        i+=1
    return qres

def L(n,p):#Legendre Symbol, where p is prime!!! I will include 2 even though techn I can't yet. This will  be nicer for the generalizations of this function.
    if (n/p == n//p):
        return 0
    else:
        if(p==2):
            if (n%8 == 1 or n%8 == -1):
                return 1
            if (n%8 == 3 or n%8 == -3):
                return -1
        if (n in qres(p)):
            return 1
        return -1

def J(m,n):# Jacobi Symbol, generalization of L. n is odd.
    M = list(map(lambda x : L(m,x),Primes(n)))
    result = reduce(lambda x,y : x*y, M)
    return result

def K(m,n):#Kronecker Symbol, a further generalization of the Legendre Symbol. All integers work.
    if (n==0):
        if (m == 1 or m == -1):
            return 1
        return 0
    if (n == -1):
        if (m < 0):
            return -1
        return 1
    if (m < 0):
        return -J(m,n)
    return J(m,n)
    
def pop(n):# returns the set (a list, but for my purposes it doesn't matter) of Gaussian integers with a given norm n. At least I think so. This might be new work; I hope someone has done this before but idk..
    i = GaussInt(0,1)
    result = []
    l = len(qres(n))
    for  y in qres(n):
        x = list(map(lambda x : K(x,n)*K(x*x+y,n),range(1,l+1)))
        a = reduce(lambda x,y: x+y, x)
        z = GaussInt(a, int(sqrt(n-a**2)))
        result.append(z)
        result.append(i*z)
        result.append((i**2)*z)
        result.append((i**3)*z)
        result.append(~z)
        result.append(i*(~z))
        result.append((i**2)*(~z))
        result.append((i**3)*(~z))
    q = []
    counter = 0
    for g in result:
        if counter == 0:
            counter += 1
            q.append(g)
        if g not in q:
            counter = 0
            q.append(g)
    result = q
    return result




class GaussInt(object):# Gaussian Integer class
    
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.FN = x*x+y*y
        self.list = [x,y]
        self.title = str(x) + ' + ' + str(y)+'i'
    
    def __str__(self):
        if (self.y == 1):
            return (str(self.x)+' + i')
        return self.title
    
    def __radd__(self,other):
        return self + other
    
    def __add__(self,other):
        return GaussInt(self.x+other.x,self.y + other.y)
    
    def __lt__(self,other):
        return self.FN < other.FN
    def __eq__(self,other):
        return self.title == other.title
    def __mul__(self,other):
        if (type(other) == int):
            return GaussInt(self.x * other, self.y*other)
        return GaussInt(self.x*other.x - self.y*other.y, self.x*other.y+self.y*other.x)
   
    __rmul__ = __mul__
    
    def __invert__(self):
        return GaussInt(self.x, -self.y)
    
    
    def __pow__(self, n):
        x = self
        f = lambda x : x > 1
        while(f(n)):
            x = x*self
            n -= 1
        return x
    
    def cplot(self,iter):# still wonky, but soon I'll figure out the runtime issus
        val = self.FN
        XLIST = []
        YLIST = []
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        while (iter>0):
            listy = pop(val)
            val = coll(val)
            XLIST = list(map(lambda f: f.x, listy))
            YLIST = list(map(lambda f: f.y, listy))
            iter -= 1
            
            plt.plot(XLIST,YLIST, 'ro')
            #for xy in zip(XLIST, YLIST):                                       
            #   ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='data') 
            
    def randomize(n):
        x = GaussInt(randint(0,n),randint(0,n))
        return x