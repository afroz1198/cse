#!/usr/bin/env python
# coding: utf-8


MOD=1000000007

def binaryExponentiationRecursive(x,n):
    if n==0:
        return 1;
    elif n%2==1: # n is odd
        return x*binaryExponentiationRecursive(x*x,(n-1)/2);
    else: # n is even
        return binaryExponentiationRecursive(x*x,n/2);
        



def binaryExponentiationIterative(x,n):    
    result=1
    while(n>0):
        if(n % 2 ==1):
            result=result * x
            n = n-1
            continue
        x=x*x
        n=n/2
    return result

def modularExponentiationRecursive(x,n,C):
    if n==0:
        return 1;
    elif n%2==1:
        return x*modularExponentiationRecursive((x*x)%C,(n-1)/2,C)
    else:
        return modularExponentiationRecursive((x*x)%C,n/2,C)



def modularExponentiationIterative(x,n,C):
    result = 1
    while n>0:
        if n%2 == 1:
            result = (result*x)%C
            n = n-1
        else:
            x = (x*x)%C
            n = n/2
    return result



def fibRecursive(n):
    if n==0:
        return 0
    elif n == 1:
        return 1
    return fibRecursive(n-1)+fibRecursive(n-2)



def fibIterative(n):
    """
    returns f_n%MOD
    
    """
    #bottom-up/iterative approach
    n1,n2 = 0,1
    for i in range(1,n):
        n3 = n1+n2
        n1 = n2
        n2 = n3
    return n3%MOD



def fibFastDoubling(n):
    """
    returns (f_n%MOD,f_n+1%MOD)
    
    """
    if n==0:
        return (0,1)
    # go k to 0 without doing any calculations,then 0 to k with calculations -> stack magic
    # also called recursive/top-down approach
    (f_n,f_nplus1) = fibFastDoubling(n>>1) # n>>1 == n//2 ,cause we have doubled n -> 2n ,check proof
    f_2n = (f_n*(2*f_nplus1 - f_n))%MOD
    f_2nplus1 = (f_nplus1**2 + f_n**2)%MOD
    
    if n&1:
        return (f_2nplus1,f_2n+f_2nplus1)
    else:        
        return (f_2n,f_2nplus1)

def matmul(M, N):
    """
    2x2 Matrix multiplication
    
    Parmeters
    ---------
    M: Matrix of size (2,2)
    N: Matrix of size (2,2)
    
    """
    # List to store matrix multiplication result
    R = [[0,0],[0,0]]
 
    for i in range(0, 2): 
        for j in range(0, 2):
            for k in range(0, 2): 
                R[i][j] += (M[i][k] * N[k][j])%MOD
     
    return R

def matBinExpo(mat,n):
    """
    Binary Exponentiation on 2x2 matrix 
    
    Parmeters
    ---------
    mat: Matrix of size 2x2
    n: exponent
    """
    
    res = mat
    while n>0:
        if n%2 == 1:
            res = matmul(res,mat)
            n = n-1
            continue
        mat = matmul(mat,mat)
        n = n/2
    return res

def fibUsingmatBinExpo(n):
    """
    returns f_n%MOD
    
    """
    transform_matrix = matBinExpo([[1,1],[1,0]],n-2)
    return transform_matrix[0][0]



def gcdBrute(a,b):
    m = min(a,b)
    gcd = 0
    for i in reversed(range(m)):
        if a%i ==0 and b%i ==0:
            gcd = i
            return gcd




def gcdEuclid(a,b):
    if b==0:
        return a
    else:
        return gcdEuclid(b,a%b)



import math
def gcdExtendedEuclid(A,B):
    """
    Implements extended Euclidean Algo (Check Page-4 of Ref[1] / Ref[2])
    A*x+B*y = GCD(A,B)
    returns (gcd,x,y)
    
    gcd: GCD(A,B)
    x: Modular Multiplicative inerse of (A%B) when A,B are coprime (i.e GCD(A,B)=1)
    y: Modular Multiplicative inerse of (B%A) when A,B are coprime (i.e GCD(A,B)=1)
    
    Parameters
    ----------
    A : value1
    B : value2
    
    """
    if B == 0:
        gcd=A
        x=1
        y=0
    else:
        gcd,x,y = gcdExtendedEuclid(B,A%B)
        temp = x
        x = y
        y = temp - math.floor(A/B)*y
    return gcd,x,y



def gcdExtendedEuclidImproved(A,B):
    """
    Implements extended Euclidean Algo (Check Page-4 of Ref[1] / Ref[2])
    A*x+B*y = GCD(A,B)
    returns (gcd,x,y)
    
    gcd: GCD(A,B)
    x: Modular Multiplicative inerse of (A%B) when A,B are coprime (i.e GCD(A,B)=1)
    y: Modular Multiplicative inerse of (B%A) when A,B are coprime (i.e GCD(A,B)=1)
    
    Parameters
    ----------
    A : value1
    B : value2
    
    """
    if B == 0:
        gcd=A
        x=1
        y=0
    else:
        Q = A//B
        gcd,x,y = gcdExtendedEuclidImproved(B,A-Q*B)
        temp = x
        x = y
        y = temp -Q*y
    return gcd,x,y




def gcdExtendedEuclidImprovedIterative(A,B):
    """
    Implements extended Euclidean Algo (Check Page-4 of Ref[1] / Ref[2])
    A*x+B*y = GCD(A,B)
    returns (gcd,x,y)
    
    a1: GCD(A,B)
    x: Modular Multiplicative inerse of (A%B) when A,B are coprime (i.e GCD(A,B)=1)
    y: Modular Multiplicative inerse of (B%A) when A,B are coprime (i.e GCD(A,B)=1)
    
    Parameters
    ----------
    A : value1
    B : value2
    
    """
    x=1;y=0;x1=0;y1=1;
    while B>0:
        q = A//B
        x,x1 = x1,x-q*x1
        y,y1 = y1,y-q*y1
        A,B = B,A-q*B
    return (A,x,y)




def modInverseBrute(A,C):
    """
    returns Modular Multiplicative Inverse (A*B)%C=1 where B=A**-1 and B should be in [1,C-1]
    uses BruteForce Approach
    
    Parameters
    ----------
    A : value
    C : Mod value
    """
    A = A%C
    for B in range(1,C):
        if (A*B)%C == 1:
            return B
        


def modInverseUsingEEA(A,C):
    """
    returns Modular Multiplicative Inverse (A*B)%C=1 where B=A**-1 and B should be in [1,C-1]
    uses extended Euclidean algo 
    
    Parameters
    ----------
    A : value
    C : Mod value
    
    """
    _,x,_ = gcdExtendedEuclidImprovedIterative(A,C)
    return (x%C+C)%C # because x could be negative (Check Page-3 of Ref[1]/ Ref[2]-> Extended Euclidean algorithm)
    
    

def modInverseFermat(A,C):
    return modularExponentiationIterative(A,C-2,C)



def modInversesTill(m):
    inv = [1]*m
    inv[1] = 1
    for i in range(2,m):
         inv[i] = m-(m//i)*inv[m%i]%m
    return inv


# ##  Refereces
# 1. Notes (Modular Arithmetic)
# 2. https://www.hackerearth.com/practice/math/number-theory/basic-number-theory-1/tutorial/
# 3. https://en.wikipedia.org/wiki/Fermat%27s_little_theorem
# 4. https://www.khanacademy.org/computing/computer-science/cryptography/modarithmetic/a/the-euclidean-algorithm





