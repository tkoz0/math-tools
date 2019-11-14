from fractions import Fraction
from functools import reduce
from typing import List
from typing import Union
import copy

class DimensionError(Exception): pass
class NotInvertibleError(Exception): pass

# RationalMatrix: matrix whose entries are rational numbers
# Invariant: self.M is a list of at least 1 row. Each row has equal
# positive length and consists of entries that are of type Fraction

class RationalMatrix:

    # check representation invariant
    def _check_rep(self):
        assert type(self.M) == list and len(self.M) > 0
        rlen = len(self.M[0])
        for r in self.M:
            assert len(r) == rlen
            for c in r:
                assert type(c) == Fraction

    # 1st argument is RationalMatrix: creates a copy
    # 1st argument is int and no 2nd argument: creates a zeroed square matrix
    # 2 arguments are int: creates a zeroed r by c matrix
    def __init__(self, r:Union[int,'RationalMatrix'], c:int = None):
        if type(r) == type(self):
            self.M = copy.deepcopy(r.M)
        if c == None: c = r
        if type(r) != int or type(c) != int: raise TypeError()
        if r < 1 or c < 1: raise ValueError()
        self.M = [[Fraction(0)]*c for _ in range(r)]

    # check matrix dimensions
    def rows(self) -> int: return len(self.M)
    def cols(self) -> int: return len(self.M[0])
    def isSquare(self) -> bool: return self.rows() == self.cols()

    # helper functions
    def _delta(r,c): return 1 if r == c else 0 # kronecker delta
    def _parity(self,n): return 1 if n % 2 == 0 else -1 # (-1)^n
    # checks that each matrix entry satisfies a function of the position
    def _checkM(self,f): # f: r,c,v -> bool
        for r in range(self.rows()): # return false on first unsatisfied entry
            for c in range(self.cols()):
                if not f(r,c,self.M[r][c]):
                    return False
        return True # all satisfy
    # matrix split for determinant computation with cofactor expansion
    def _det_split(self,r,c):
        M = [r[:c]+r[c+1:] for r in self.M[0:r]+self.M[r+1:]]
        D = RationalMatrix(1)
        D.M = M
        return D

    # check matrix types, they use lambdas mapping r,c,v to boolean
    # the value v at each position r,c must satisfy some condition
    def isIdentity(self) -> bool: # entries equal to kronecker delta
        return self.isSquare() and self._checkM(lambda r,c,v: v == _delta(r,c))
    def isDiagonal(self) -> bool: # r != c implies v == 0
        return self.isSquare() and self._checkM(lambda r,c,v: r == c or v == 0)
    def isUpperTriangle(self) -> bool: # r > c implies v == 0
        return self.isSquare() and self._checkM(lambda r,c,v: r <= c or v == 0)
    def isLowerTriangle(self) -> bool: # r < c implies v == 0
        return self.isSquare() and self._checkM(lambda r,c,v: r >= c or v == 0)
    def isStrictUpperTriangle(self) -> bool: # r >= c implies v == 0
        return self.isSquare() and self._checkM(lambda r,c,v: r < c or v == 0)
    def isStrictLowerTriangle(self) -> bool: # r <= c implies v == 0
        return self.isSquare() and self._checkM(lambda r,c,v: r > c or v == 0)

    # modify a single entry
    def set(self, r:int, c:int, v:Union[int,Fraction]):
        if type(v) == Fraction: self.M[r][c] = v
        elif type(v) == int: self.M[r][c] = Fraction(v)
        else: raise TypeError()

    # access matrix, getRow and getMatrix expose representation
    def get(self,r:int,c:int) -> Fraction: return self.M[r][c]
    def getRow(self,r:int) -> List[Fraction]: return self.M[r]
    def getCol(self,c:int) -> List[Fraction]:
        return [self.M[r][c] for r in range(self.rows())]
    def getMatrix(self): return self.M

    # returns a new matrix that is equal to the transpose
    def transpose(self) -> 'RationalMatrix': # returns a new matrix
        T = RationalMatrix(1)
        T.M = [ [M[c][r] for c in range(self.rows())]
                for r in range(self.cols())]
        return T

    # augments the matrix (such as for solving a linear system) by adding a
    # number to the end of each row
    def augment(self,b): # appends entries to each row
        if len(b) != self.rows(): raise ValueError()
        for n in b: # check types
            if type(n) != Fraction and type(n) != int:
                raise ValueError()
        for r in range(self.rows()): # append to each row
            if type(b[r]) == Fraction: self.M[r].append(b[r])
            else: self.M[r].append(Fraction(b[r]))

    # compute determinant using cofactor expansion (slow)
    def detCofactor(self) -> Fraction:
        if not self.isSquare(): raise DimensionError()
        n = self.rows() # for convenience, since matrix must be square
        if n == 1: return self.M[0][0] # base case
        # find row/col with most zeroes to optimize computation
        # count zeroes in each row and column
        row0 = [sum(1 if n == 0 else 0 for n in r) for r in self.M]
        col0 = [sum(1 if self.M[r][c] == 0 else 0 for r in range(n))
                for c in range(n)]
        use_row = True # expand over row or column
        best_index = 0
        most_zeroes = 0
        for r in range(1,n): # find row with more zeroes
            if row0[r] > most_zeroes: best_index, most_zeroes = r, row0[r]
        for c in range(n): # find column with more zeroes
            if col0[c] > most_zeroes:
                use_row = False
                best_index, most_zeroes = c, col0[c]
        if use_row: # optimize computation by using row/col with
            return sum(self._parity(best_index+c) * self.M[best_index][c]
                       * self._det_split(best_index,c).detCofactor()
                       for c in range(n))
        else: return sum(self._parity(r+best_index) * self.M[r][best_index]
                         * self._det_split(r,best_index).detCofactor()
                         for r in range(n))

    # compute determinant using row reduction (faster)
    def detRowReduce(self) -> Fraction:
        if not self.isSquare(): raise DimensionError()
        M = self.clone() # row reduce a copy
        n = M.rows()
        multiplier = Fraction(1) # multiply at end, due to swapping rows
        for r in range(n): # row reduce
            if M.M[r][r] == 0: # swap in a pivot row
                pr = -1 # pivot row
                for rr in range(r+1,n): # find pivot row
                    if M.M[rr][r] != 0:
                        pr = rr
                        break
                if pr == -1: return M.M[r][r] # cannot find pivot row
                # (return existing zero instead of making new fraction)
                else: # swap and negate multiplier
                    M.swapRows(r,pr)
                    multiplier *= -1
            # zero out column for upper triangle form
            for rr in range(r+1,n): # subtract with leading M.M[rr][r]
                M.addRowMult(-M.M[rr][r]/M.M[r][r],r,rr)
        # multiplier times product of diagonal
        return multiplier * reduce(lambda x,y: x*y,
                                   [M.M[r][r] for r in range(n)])

    # return a new matrix that is the inverse, exception if not invertible
    def inverse(self) -> 'RationalMatrix':
        if not self.isSquare(): raise DimensionError()
        n = self.rows()
        # convert (M|N) = (original|I) --> (M|N) = (I|inverse)
        M = self.clone()
        N = RationalMatrix.identity(n)
        for r in range(n): # row reduce to M being an upper triangle matrix
            if M.M[r][r] == 0: # swap in a pivot row
                pr = -1 # pivot row
                for rr in range(r+1,n):
                    if M.M[rr][r] != 0:
                        pr = rr
                        break
                if pr == -1: raise NotInvertibleError() # cannot find pivot row
                else: # swap in pivot row
                    M.swapRows(r,pr)
                    N.swapRows(r,pr)
            row_mult = 1/M.M[r][r]
            M.multRow(row_mult,r) # have leading element become 1
            N.multRow(row_mult,r)
            for rr in range(r+1,n): # zero out column for upper triangle form
                row_mult = -M.M[rr][r]
                M.addRowMult(row_mult,r,rr)
                N.addRowMult(row_mult,r,rr)
        for r in range(n-1,-1,-1): # reduce upward to get M = I
            for rr in range(r-1,-1,-1):
                row_mult = -M.M[rr][r]
                ###### TODO MAY BE ABLE TO ELIMINATE NEXT LINE #######
                M.addRowMult(row_mult,r,rr)
                N.addRowMult(row_mult,r,rr)
        return N

    # row manipulations, these modify the matrix
    def swapRows(self,r1:int,r2:int):
        self.M[r1], self.M[r2] = self.M[r2], self.M[r1]
    def multRow(self,a:Union[int,Fraction],r:int):
        if type(a) != int and type(a) != Fraction: raise TypeError()
        self.M[r] = [a*c for c in self.M[r]]
    def addRowMult(self,a:Union[int,Fraction],r1:int,r2:int): # add a*r1 to r2
        if type(a) != int and type(a) != Fraction: raise TypeError()
        for c in range(self.cols()):
            self.M[r2][c] += a * self.M[r1][c]

    # string representation: each row as [n1 n2 ... nc], columns are aligned
    def __str__(self):
        M = [list(map(str,r)) for r in self.M]
        col_len = [0]*self.cols() # compute column widths needed for alignment
        for r in M:
            for i,c in enumerate(r):
                col_len[i] = max(col_len[i],len(c))
        pad = lambda s, c : (col_len[c]-len(s))*' ' + s
        return '\n'.join('['+' '.join(pad(r[c],c)
                                      for c in range(self.cols()))+']'
                         for r in M)

    # return a copy of this matrix
    def clone(self):
        M = RationalMatrix(1)
        M.M = copy.deepcopy(self.M)
        return M

    # some types of matrices
    def identity(n:int): # identity matrix
        M = RationalMatrix(n)
        for i in range(n): M.set(i,i,1)
        return M
    def diagonal(e:list): # diagonal matrix with entries
        if type(e) != list: raise TypeError()
        M = RationalMatrix(len(e))
        for i,n in enumerate(e):
            M.set(i,i,n)
        return M

    def scale(self,s:Union[int,Fraction]): # scalar multiple
        if type(s) != int and type(s) != Fraction: raise TypeError()
        self.M = [[s * self.M[r][c] for c in range(self.cols())]
                  for r in range(self.rows())]

    # operator overloading
    def __add__(self,M):
        if self.rows() != M.rows() or self.cols() != M.cols():
            raise DimensionError()
        N = RationalMatrix(1)
        N.M = [[self.M[r][c] + M.M[r][c] for c in range(self.cols())]
               for r in range(self.rows())]
        return N
    def __sub__(self,M):
        if self.rows() != M.rows() or self.cols() != M.cols():
            raise DimensionError()
        N = RationalMatrix(1)
        N.M = [[self.M[r][c] - M.M[r][c] for c in range(self.cols())]
               for r in range(self.rows())]
        return N
    def __mul__(self,M):
        if self.cols() != M.rows(): raise DimensionError()
        N = RationalMatrix(1)
        N.M = [[sum(self.M[r][i] * M.M[i][c] for i in range(self.cols()))
                for c in range(M.cols())]
               for r in range(self.rows())]
        return N
    def __pow__(self,p:int): # TODO possibly use diagonalization
        if type(p) != int: raise TypeError()
        if not self.isSquare(): raise DimensionError()
        n = self.rows()
        if self.isDiagonal(): # fast case
            return RationalMatrix.diagonal([self.M[r][r]**p for r in range(n)])
        if p == 0: return RationalMatrix.identity(self.rows())
        elif p >= 1:
            M = self.clone()
            for _ in range(1,p): M *= self
            return M
        elif p <= -1:
            INV = self.inverse()
            M = INV.clone()
            for _ in range(1,-p): M *= INV
            return M
    def __eq__(self,M):
        return self.rows() == M.rows() and \
            reduce(lambda x,y: x and y, [self.M[r] == M.M[r]
                   for r in range(self.rows())])

def _fix(A):
    for r in range(A.rows()):
        for c in range(A.cols()):
            A.M[r][c] = Fraction(A.M[r][c])

if __name__ == '__main__':
    A = RationalMatrix(1)
    A.set(0,0,5)
    print(A)
    print('det =',A.detRowReduce(), '(should be 5)')
    print('det =',A.detCofactor(), '(should be 5)')
    A.M = [[1,2],
           [3,4]]
    _fix(A)
    print(A)
    print('det =',A.detRowReduce(), '(should be -2)')
    print('det =',A.detCofactor(), '(should be -2)')
    A.M = [[1,2,3],
           [2,-2,2],
           [-1,0,5]]
    _fix(A)
    print(A)
    print('det =',A.detRowReduce(), '(should be -40)')
    print('det =',A.detCofactor(), '(should be -40)')
    A.M = [[1,-1,1,-1],
           [5,0,6,7],
           [3,-1,-1,-1],
           [1,0,1,1]]
    _fix(A)
    print(A)
    print('det =',A.detRowReduce(), '(should be 6)')
    print('det =',A.detCofactor(), '(should be 6)')
    print('begin inverse')
    D = A.inverse()
    print(D)
    print(D*A)
    print('end inverse')
    A.M = [[Fraction(-2,3),Fraction(-1,4)],
         [Fraction(-1,5),Fraction(1,2)]]
    print(A)
    print('det =',A.detRowReduce(), '(should be -23/60)')
    print('det =',A.detCofactor(), '(should be -23/60)')
    A.M = [[1,-1,1,-1],
           [5,0,6,7],
           [3,-1,-1,-1],
           [-1,-2,-6,-9]] # last row = r1 - r2 + r3
    _fix(A)
    print(A)
    print('det =',A.detRowReduce(), '(should be 0)')
    print('det =',A.detCofactor(), '(should be 0)')
    try:
        Z = A.inverse()
        print(Z)
    except NotInvertibleError:
        print('cannot invert')
    print('identity')
    print(RationalMatrix.identity(5))
    print('diagonal')
    print(RationalMatrix.diagonal([1,-2,3,-4,5]))
    print('diagonal power')
    print(RationalMatrix.diagonal([1,2,3])**-3)
    print('power')
    A.M = [[1,2],
           [3,4]]
    _fix(A)
    print(A**0)
    print(A**1)
    print(A**2)
    print(A**3)
    B = RationalMatrix(2)
    B.M = [[7,10],
           [15,22]]
    _fix(B)
    print(A**2 == B)
    print('row operations')
    I = RationalMatrix.identity(2)
    print(I)
    I.swapRows(0,1)
    print(I)
    I.multRow(2,1)
    print(I)
    I.multRow(-3,0)
    print(I)
    I.addRowMult(-1,1,0)
    print(I)
    I.addRowMult(2,0,1)
    print(I)
    print('inverse')
    print(A)
    print(A.inverse())
    print(A*A.inverse())
    print(A == A.inverse())
    print('inverse powers')
    print(A)
    print(A**0)
    print(A**-1)
    print(A**-2)
    print('identity again')
    I = RationalMatrix.identity(3)
    print(I.inverse())
    print(I*I == I)
    print(I*I.inverse())
