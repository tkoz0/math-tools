# usage: row_reduce <r> <c> <g | gj>
# g for gauss elim (to upper triangle) and gj for gauss jordan
import sys
from fractions import Fraction
R = int(sys.argv[1])
C = int(sys.argv[2])
type = sys.argv[3]
assert type == 'g' or type == 'gj'
assert R > 0 and C > 0
M = []
for _ in range(R):
    line = input().split()
    assert len(line) == C
    M.append(list(map(Fraction,line)))

def arrow(text):
    return f'\\xrightarrow{{{text}}}'

def pmatrix(M):
    M = [[str(c) for c in r] for r in M]
    col_len = [0]*C
    for r in M:
        for i,c in enumerate(r):
            col_len[i] = max(col_len[i],len(c))
    pad = lambda s, c : (col_len[c]-len(s))*' ' + s
    return '\\begin{pmatrix}\n' + \
           '\n'.join(' & '.join(pad(r[c],c) for c in range(C)) + ' \\\\'
                     for r in M) + \
           '\n\\end{pmatrix} \\\\'

def elim_g(): # gauss elim to upper triangle
    global M
    print('\\\\\\text{Gaussian Elimination}\\\\\\\\')
    r = 0
    for c in range(C):
        if r == R: break # done
        if M[r][c] == 0: # swap pivot row
            pr = -1
            for rr in range(r+1,R):
                if M[rr][c] != 0:
                    pr = rr
                    break
            if pr == -1: continue # move to next column
            print(arrow('R_{%d}\\leftrightarrow R_{%d}'%(r+1,pr+1)))
            M[r], M[pr] = M[pr], M[r]
            print(pmatrix(M))
        assert M[r][c] != 0 # pivot number
        mult = 1/M[r][c]
        if M[r][c] != 1:
            print(arrow('(%s)R_{%d}\\to R_{%d}'%(str(mult),r+1,r+1)))
            M[r] = [M[r][cc]*mult for cc in range(C)]
            assert M[r][c] == 1
            print(pmatrix(M))
        for rr in range(r+1,R):
            if M[rr][c] == 0: continue # already zero
            mult = -M[rr][c] # add mult times row r to row rr
            print(arrow('R_{%d}+(%s)R_{%d}\\to R_{%d}'
                        %(rr+1,str(mult),r+1,rr+1)))
            M[rr] = [M[rr][cc]+mult*M[r][cc] for cc in range(C)]
            assert M[rr][c] == 0
            print(pmatrix(M))
        r += 1

def elim_gj(): # gauss jordan elimination
    global M
    print('\\\\\\text{Gauss-Jordan Elimination}\\\\\\\\')
    for r in range(R-1,-1,-1): # eliminate upward using pivot rows
        c = 0
        while c < C and M[r][c] == 0: c += 1 # find pivot 1
        if c == C: continue # cannot find 1 on this row
        for rr in range(r-1,-1,-1): # eliminate above rows
            mult = -M[rr][c]
            if mult == 0: continue # no need if it is already zero
            print(arrow('R_{%d}+(%s)R_{%d}\\to R_{%d}'
                        %(rr+1,str(mult),r+1,rr+1)))
            M[rr] = [M[rr][cc]+mult*M[r][cc] for cc in range(C)]
            print(pmatrix(M))

print('\\documentclass{article}')
print('\\usepackage{fullpage}')
print('\\usepackage{mathtools}')
print('\\usepackage{parskip}')
print('\\begin{document}')
print('$')

print(arrow('\\text{initial}'))
print(pmatrix(M))

elim_g()
if type == 'gj': elim_gj()

print('$')
print('\\end{document}')
