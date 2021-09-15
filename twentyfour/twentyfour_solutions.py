import json
import sys
import twentyfour
from fractions import Fraction

def solve(nums):
    nums = list(Fraction(num) for num in nums)
    values = twentyfour.search(nums)
    return values[Fraction(24)] if Fraction(24) in values else []

solutions = dict()

def i2c(i):
    if i < 10: return str(i)
    return {10:'T',11:'J',12:'Q',13:'K'}[i]

for a in range(1,14):
    for b in range(a,14):
        for c in range(b,14):
            for d in range(c,14):
                sys.stderr.write(str([a,b,c,d])+'\n')
                key = ''.join(map(i2c,[a,b,c,d]))
                sols = solve([a,b,c,d])
                solutions[key] = list(map(twentyfour.tree2str,sols))
                sys.stderr.write('key %s, %d solutions\n'%(key,len(sols)))

sys.stdout.write(json.dumps(solutions,separators=(',',':')))
