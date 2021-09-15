import json
import pyperclip
import sys

sys.stderr.write('input 4 digits and a solution is printed\n')

with open('twentyfour_solutions_small.json') as infile:
    solutions = json.loads(infile.read())

chars = '0123456789TJQK'
comp = dict()
for i,c in enumerate(chars):
    comp[c] = i

while True:
    try:
        line = input()
    except:
        break
    try:
        digits = ''.join(sorted(filter(lambda x : x in chars, line)))
        digits = ''.join(sorted(digits,key=lambda k:comp[k]))
        if digits in solutions:
            print(solutions[digits])
            pyperclip.copy(solutions[digits])
        else:
            print('invalid input')
    except:
        print('input error')
