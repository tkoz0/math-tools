import twentyfour

def solution(nums):
    values = twentyfour.search(nums)
    return values[24] if 24 in values else 'none'

for a in range(1,10):
    for b in range(a,10):
        for c in range(b,10):
            for d in range(c,10):
                print(a,b,c,d,twentyfour.t2s(solution([a,b,c,d])))
