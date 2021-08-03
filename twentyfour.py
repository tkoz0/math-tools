import pyperclip
import sys
from fractions import Fraction

target = Fraction(sys.argv[1])
inputs = list(map(Fraction,sys.argv[2:]))
assert len(inputs) > 0

# binary operators that can always be used
# division uses / because Fraction is being used for more flexibility
ops = \
{
    '+': lambda x,y : x+y,
    '-': lambda x,y : x-y,
    '*': lambda x,y : x*y,
    '/': lambda x,y : x/y
}

# returns a dictionary mapping possible value to a parse tree for computing it
# the parse tree is represented as a Fraction or recursive list:
# [operation, left, right], left/right are recursive or Fraction
# operation is a string key from ops
def search(nums):
    assert len(nums) > 0
    if len(nums) == 1:
        return {nums[0]:nums[0]}
    # use the bits to determine set inclusion for the partitioning
    S = dict()
    for subsetnum in range(1,2**len(nums)-1):
        A = [] # 0
        B = [] # 1
        for i in range(len(nums)):
            if subsetnum % 2 == 0:
                A.append(nums[i])
            else:
                B.append(nums[i])
            subsetnum = subsetnum // 2
        A_vals = search(A)
        B_vals = search(B)
        for a in A_vals:
            for b in B_vals:
                for op in ops:
                    # try/except because zero division is possible
                    try:
                        r = ops[op](a,b)
                    except:
                        continue
                    if r in S: continue # already found a way to compute it
                    S[r] = [op,A_vals[a],B_vals[b]]
    return S

def t2s(t,init=True): # tree to string
    if type(t) != list:
        return str(t)
    r = t2s(t[1],False)+t[0]+t2s(t[2],False)
    return r if init else '('+r+')'

result = search(inputs)

if target in result:
    print('found solution')
    print(target,'=',t2s(result[target]))
    pyperclip.copy(t2s(result[target]))
else:
    print('no solution')
    pyperclip.copy('none')
