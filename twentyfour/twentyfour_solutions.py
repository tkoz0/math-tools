import json
import sys
import twentyfour

def solve(nums):
    values = twentyfour.search(nums)
    return values[24] if 24 in values else 'none'

solutions = dict()

for a in range(1,10):
    for b in range(a,10):
        for c in range(b,10):
            for d in range(c,10):
                solution = twentyfour.t2s(solve([a,b,c,d]))
                solutions[''.join(map(str,[a,b,c,d]))] = solution

sys.stdout.write(json.dumps(solutions,separators=(',',':')))
