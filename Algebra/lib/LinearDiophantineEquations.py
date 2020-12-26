#!/usr/bin/env python
# coding: utf-8

# ## Linear DioPhantine equations 



from BasicMath1 import gcdExtendedEuclidImprovedIterative


def isDegenerateCase(a,b,c):
    return True if a==0 and b==0 else False


def findAnySolution(a,b,c):
    if isDegenerateCase(a,b,c):
        return False
    (g,x_g,y_g) = gcdExtendedEuclidImprovedIterative(a,b)
    if c%g:
        return False
    x_0 = x_g*(c//g)
    y_0 = y_g*(c//g)
    if a<0:
        x_0 = -x_0
    if b<0:
        y_0 = -y_0
    return (g,x_0,y_0)


def shiftSolution(x_0,y_0,a,b,k):
    return x_0+k*b, y_0-k*a

def findAllSolutions(a,b,c,minx,maxx,miny,maxy):
    if a==0 or b==0:
        return 1;
    if not findAnySolution(a,b,c):
        return 0;
    g,x,y = findAnySolution(a,b,c)
    
    # from eq(4),eq(5),simplifies shiftSolution equations
    a //= g;
    b //= g;
    
    sign_a = 1 if a>0 else -1
    sign_b = 1 if b>0 else -1
    
    x,y = shiftSolution(x,y,a,b,(minx-x)//b)
    if x < minx:
        x,y = shiftSolution(x,y,a,b,sign_b)
    if x > maxx:
        return 0;
    lx1 = x
    
    x,y = shiftSolution(x,y,a,b,(maxx-x)//b)
    if x > maxx:
        shiftSolution(x,y,a,b,-sign_b)
    rx1 = x
    
    x,y = shiftSolution(x,y,a,b,-(miny-y)//a)
    if y < miny:
        x,y = shiftSolution(x,y,a,b,-sign_a)
    if y > maxy:        
        return 0;
    lx2 = x
    
    x,y = shiftSolution(x,y,a,b,-(maxy-y)//a)
    if y > maxy:
        x,y  = shiftSolution(x,y,a,b,sign_a)
    rx2 = x
    
    if lx2 > rx2:
        temp = lx2;
        lx2 = rx2;
        rx2 = temp
    # Intersection
    lx = max(lx1,lx2)
    rx = min(rx1,rx2)
    
    if lx > rx:
        return 0;
    return (rx-lx)//abs(b)+1;   
    


