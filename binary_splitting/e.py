'''
binary splitting test

e^x = sum x^n/n! from n=0 to n=inf

binary splitting method

sum p_n/q_n from n=a to n=b (a<b) (all p and q are integers, qn != 0)
summation result P(a,b)/Q(a,b)

'''

import math

def standard_summation():
    MAX_INDEX = 50
    summation = 0.0
    n = 0
    while n <= MAX_INDEX:
        newsummation = summation + 1.0/math.factorial(n)
        if newsummation == summation:
            print('stopped at term index',n)
            break
        summation = newsummation
        n += 1
    print(summation)

# binary splitting routine for e
# [a,b] inclusive index range
# call computePQ(0,n) to compute e
def computePQ(a,b):
    assert a <= b
    if a == b: # compute directly
        return (1,math.factorial(a))
    elif a+1 == b: # compute directly
        # 1/a! + 1/((a+1)a!)
        # = (a+2)/((a+1)!)
        n,d = a+2,math.factorial(a+1)
        g = math.gcd(n,d)
        return (n//g,d//g)
    else:
        # split interval and combine results
        mid = (a+b)//2
        p1,q1 = computePQ(a,mid) # p1/q1
        p2,q2 = computePQ(mid+1,b) # p2/q2
        # sum = (p1q2+p2q1)/(q1q2)
        n,d = p1*q2+p2*q1,q1*q2
        g = math.gcd(n,d)
        return (n//g,d//g)

P,Q = computePQ(0,17) # 18 terms to get full precision with IEEE754 double
ps,qs = str(P),str(Q)
print(len(ps),'digits, P =',ps)
print(len(qs),'digits, Q =',qs)
print('gcd(P,Q) =',math.gcd(P,Q))
print('e ~= P/Q =',P/Q)
print('P/Q == math.e ?',P/Q==math.e)
print('e =',math.e)
print('|P/Q-math.e| =',abs(P/Q-math.e))
