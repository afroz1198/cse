#!/usr/bin/env python
# coding: utf-8



def checkPrimeBrute(n):
    count =0
    for i in range(1,n+1):
        if n%i == 0:
            count += 1
    if count==2:
        return True
    else:
        return False



def checkPrimeImproved(n):
    count =0 
    d=1
    while d*d<=n:
        if n%d == 0:
            if d*d == n: # sqrt(25) = 5*5,so count only one factor,avoid redundancy
                count += 1 
            else:
                count += 2 # Lemma 1 ,counting  i,n/i
        d += 1
    if count==2:
        return True
    else:
        return False



def sieveMask(n):
    """
    Sieve of Eratosthenes
    returns Prime number boolean Mask of length n
    
    Parameters
    ----------
    n: integer value
    
    """
    isPrime = [True]*(n+1)
    isPrime[0],isPrime[1] = False,False
    d = 2
    while d*d<=n:
        if isPrime[d]:
            # Mark all the multiples of d as composite numbers (# sieve of eratosthenes)
            j = d*d
            while j<=n:
                isPrime[j] = False
                j += d  
        d += 1
    return isPrime


def factorize(n):
    """
    returns factors/prime factors of n
    
    """
    factors = []
    d=2
    while d*d <= n: 
        while n%d == 0:
            factors.append(d)
            n //= d # similar to lcm reduction
        d += 1
    if n != 1:
        factors.append(n)
    return factors
        


def factorizeWheel(n):
    factors = []
    for d in [2,3,5]:
        while n%d == 0:
            factors.append(d)
            n /= d
    increments = [4, 2, 4, 2, 4, 6, 2, 6] #increments to get next prime number
    i=0 # access increments values by this index 
    d=7 #start from 7,as 5 is completed
    while d*d<=n:
        while n%d == 0:
            factors.append(d)
            n /= d
        i += 1
        if i == 8:
            i = 0
        d += increments[i]
    if n>1:
        factors.append(n)
    return factors


def buildMinPrimeStore(n):
    """
    returns factors/prime factors of n using sieve of eratosthenes
    
    """
    minPrime = [0]*(n+1)
    minPrime[0],minPrime[1] = 1,1
    d=2
    while d*d <= n:
        if minPrime[d] == 0:
            # Mark all the multiples of i with minPrime as "i" (# sieve of eratosthenes)
            j=d*d
            while j<=n:
                if minPrime[j] == 0:                    
                    minPrime[j]=d
                j += d
        d += 1
    # all prime numbers are divisible by 1,itself.1 is not prime,so adding "itself"
    for i in range(2,n+1):
        if minPrime[i] == 0:
            minPrime[i] = i
    return minPrime



def factorizeSieve(n): 
    minPrime = buildMinPrimeStore(n)
    factors=[]
    while n>1:
        factors.append(minPrime[n])
        n //= minPrime[n]
    return factors



from collections import Counter
from functools import reduce
def NoOfDivisors(n):
    """
    Based on above fact
    returns no_of divisors including 1,itself
    """
    factors = factorizeWheel(n)
    counted = Counter(factors)
    for item in counted:
        counted[item] += 1
    return reduce(lambda x,y:x*y,counted.values(),1)


from collections import Counter
def sumGeometricSeries(a,r,n):
    """
    returns sum of Geometric series 
    a + ar+ ar^2 + ...+ ar^n=a(r^{n+1}-1)/{r-1} 
    """
    if r==1:
        raise ArithmeticError("r should be > 1")
    numerator = a*(r**(n+1) -1)
    denomerator = r-1
    return numerator//denomerator
def sumOfDivisors(n):
    """
    based above Facts
    returns sum of all divisors including 1,itself
    """
    factors = factorizeWheel(n)
    counted = Counter(factors)
    result = 1
    for item in counted.items():
        result *= sumGeometricSeries(1,item[0],item[1])
    return result



def phiWheel(n):
    result = n
    for d in [2,3,5]:
        if n%d == 0:
            while n%d == 0:
                n /= d
            result -= result//d
    increments = [4, 2, 4, 2, 4, 6, 2, 6]
    i=0 # access increments values by this index 
    d=7 #start from 7,as 5 is completed
    while d*d<=n:
        if n%d == 0:
            while n%d == 0:
                n /= d
            result -= result//d
        i += 1
        if i == 8:
            i = 0
        d += increments[i]
    if n>1:
        result -= result//n
    return result


def phi(n):
    result=n
    d=2
    while d*d <= n: 
        if n%d == 0:
            while n%d == 0:
                n //=d # similar to lcm reduction
            result -= result//d
        d += 1
    if n > 1:
        result -= result//n
    return result


# ## References
# 1. https://www.hackerearth.com/practice/math/number-theory/basic-number-theory-2/tutorial/
# 2.https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
# 3. https://cp-algorithms.com/algebra/sieve-of-eratosthenes.html
# 4. https://cp-algorithms.com/algebra/prime-sieve-linear.html
# 5. https://cp-algorithms.com/algebra/factorization.html
# 




