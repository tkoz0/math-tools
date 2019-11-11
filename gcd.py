# usage: gcd <a> <b> [extended]

import sys
a, b = int(sys.argv[1]), int(sys.argv[2])
if len(sys.argv) > 3 and sys.argv[3] == 'extended':
    extended = True
else: extended = False

print('\\documentclass{article}')
print('\\usepackage{fullpage}')
print('\\usepackage{mathtools}')
print('\\usepackage{parskip}')
print('\\begin{document}')

print(f'$a = {a}$ \\\\')
print(f'$b = {b}$ \\\\')
print('\\\\')

# adjust numbers so both are nonnegative and a >= b
if a < 0: a = -a
if b < 0: b = -b
if a < b: a,b = b,a

x, y = a, b

if extended:
    print(f'Let $x = {x}$ and $y = {y}$ \\\\')
    print('\\\\')

def linear_str(x,y): # string for linear combination of x and y
    s = ''
    if x == -1: s += '-'
    elif x != 1: s += str(x)
    s += 'x '
    if y < 0: s += '- '
    else: s += '+ '
    y = abs(y)
    if y == 1: s += 'y'
    else: s += str(y) + 'y'
    return s

if a == 0 and b == 0:
    print('Both numbers are zero, thus $\\gcd(a,b)$ does not exist.')

else:

    if extended: print('\\begin{tabular}{l|l|l}')
    else: print('\\begin{tabular}{l|l}')

    print('$\\gcd(a,b)$ & \\text{division} ',end='')
    if extended: print('& \\text{bezout} ',end='')
    print('\\\\')
    print('\\hline')

    # linear combination for bezout coefficients
    ax, ay = 1, 0
    bx, by = 0, 1

    s1 = linear_str(ax,ay)
    s2 = linear_str(bx,by)

    if extended:
        print(f'& & ${a} = {s1}$ \\\\')
        print(f'& & ${b} = {s2}$ \\\\')
        print('& & \\\\')

    print(f'$= \\gcd({a},{b})$ & ',end='')

    while b != 0: # iterate euclid gcd

        q,r = divmod(a,b) # a = q * b + r
        # r = a - q*b = (ax*x + ay*y) - q*(bx*x + by*y) = (ax - q*bx)*x + (ay - q*by)*y
        rx, ry = ax - q*bx, ay - q*by

        # division algorithm coefficients
        print(f'${a} = {q} \\times {b} + {r}$ ',end='')

        if extended:
            #s1 = f'({ax} - {q}({bx}))x' # (ax - q*bx)*x
            #s2 = f'({ay} - {q}({by}))y' # (ay - q*by)*y
            print(f'& ${r} = {a} - {q} * {b}$ \\\\')
            print(f'& & $= ({s1}) - {q} \\times ({s2})$ \\\\')
            s1 = s2
            s2 = linear_str(rx,ry)
            print(f'& & $= ({ax}-{q}({bx}))x + ({ay}-{q}({by}))y$ \\\\')
            print(f'& & $= {s2}$',end='')

        print('\\\\')

        print(f'$= \\gcd({b},{r})$ & ',end='') # begin next line

        # setup values for next iteration
        a,b = b,r
        ax, ay = bx, by
        bx, by = rx, ry

    # print final result
    if extended: print('& ',end='')
    print('\\\\')
    print(f'$= \\fbox{{{a}}}$ & ',end='')
    if extended:
        print(f'& \\fbox{{${a} = {s1}$}} \\\\')
        print(f'& & $= {x}({ax}) + {y}({bx})$ ',end='')
        pass
    print('\\\\')

    print('\\end{tabular}')

print('\\end{document}')
