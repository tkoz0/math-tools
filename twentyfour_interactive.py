import json
import pyperclip
import sys

sys.stderr.write('input 4 digits and a solution is printed\n')

with open('twentyfour_solutions.json') as infile:
    solutions = json.loads(infile.read())

while True:
    try:
        line = input()
    except:
        break
    digits = ''.join(sorted(filter(lambda x : x.isdigit(), line)))
    if digits in solutions:
        print(solutions[digits])
        pyperclip.copy(solutions[digits])
