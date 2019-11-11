from fractions import Fraction
from functools import reduce
from typing import List

class DimensionError(Exception): pass
class NotInvertible(Exception): pass

class Matrix:
    # Matrix: a representation of a matrix whose entries are rational numbers
    # Invariant: self.M is a list of at least 1 row. Each row has the same
    # positive length and consists of entries that are of type Fraction
    def _check_rep(self):
        assert type(self.M) == list and len(self.M) > 0
        rlen = len(self.M[0])
        for r in self.M:
            assert len(r) == rlen
            for c in r:
                assert type(c) == Fraction
    # Initializes the matrix with zeroes. If 1 argument is given, the matrix is
    # a square matrix.
    def __init__(self,r,c=None):
        if type(r) == type(self): self.M = r.M
        if c == None: c = r
        if type(r) != int or type(c) != int: raise TypeError()
        if r < 1 or c < 1: raise ValueError()
        self.M = [[Fraction(0)]*c for _ in range(r)]
    def rows(self) -> int: return len(self.M)
    def cols(self) -> int: return len(self.M[0])
    def is_square(self) -> bool: return self.rows() == self.cols()
    def is_identity(self) -> bool:
        if not self.is_square(): return False
        for r in range(self.rows()):
            for c in range(self.cols()):
                if self.M[r][c] != 1 if r == c else 0: return False
        return True
    def is_diagonal(self) -> bool:
        if not self.is_square(): return False
        for r in range(self.rows()):
            for c in range(self.cols()):
                if r != c and self.M[r][c] != 0: return False
        return True
    def set(self,r,c,n):
        if type(n) == Fraction: self.M[r][c] = n
        elif type(n) == int: self.M[r][c] = Fraction(n)
        else: raise TypeError()
    def get(self,r,c) -> Fraction: return self.M[r][c]
    def get_row(self,r) -> List[Fraction]: return self.M[r]
    def transpose(self) -> 'Matrix': # returns a new matrix
        T = RationalMatrix(self.cols(),self.rows())
        for r in range(self.rows()):
            for c in range(self.cols()):
                T.M[c][r] = self.M[r][c]
        return T
    def augment(self,b): # appends entries to each row
        if type(b) != list: raise TypeError()
        if len(b) != self.rows(): raise ValueError()
        for n in b:
            if type(n) != Fraction and type(n) != int: raise ValueError()
        for r in range(self.rows()):
            if type(b[r]) == Fraction: self.M[r].append(b[r])
            else: self.M[r].append(Fraction(b[r]))
    def _det_split(self,r,c): # split sub matrix for determinant computation
        N = self.M[0:r] + self.M[r+1:]
        for i in range(len(N)): N[i] = N[i][:c] + N[i][c+1:]
        R = Matrix(1)
        R.M = N
        return R
    def _parity(self,n): return 1 if n % 2 == 0 else -1
    # compute determinant using cofactor expansion (slow)
    def det(self) -> Fraction:
        if not self.is_square(): raise DimensionError()
        n = self.rows() # for convenience, since matrix must be square
        if n == 1: return self.M[0][0] # base case
        # find row/col with most zeroes to optimize computation
        row0 = [sum(1 if n == 0 else 0 for n in r) for r in self.M]
        col0 = [sum(1 if self.M[r][c] == 0 else 0 for r in range(n))
                for c in range(n)]
        use_row = True
        best_index = 0
        most_zeroes = 0
        for r in range(1,n):
            if row0[r] > most_zeroes: best_index, most_zeroes = r, row0[r]
        for c in range(n):
            if col0[c] > most_zeroes:
                use_row = False
                best_index, most_zeroes = c, col0[c]
        if use_row:
            return sum(self._parity(best_index+c) * self.M[best_index][c]
                       * self._det_split(best_index,c).det()
                       for c in range(n))
        else: return sum(self._parity(r+best_index) * self.M[r][best_index]
                         * self._det_split(r,best_index).det()
                         for r in range(n))
    def det_fast(self) -> Fraction:
        if not self.is_square(): raise DimensionError()
        M = self.clone()
        multiplier, n = Fraction(1), M.rows()
        for r in range(n): # row reduce
            if M.M[r][r] == 0: # swap in a pivot row
                pr = -1 # pivot row
                for rr in range(r+1,n):
                    if M.M[rr][r] != 0:
                        pr = rr
                        break
                if pr == -1: return Fraction(0) # cannot find pivot row
                else: # swap and negate multiplier
                    M.swap_rows(r,pr)
                    multiplier *= -1
            for rr in range(r+1,n): # zero out column for upper triangle form
                M.add_row_mult(-M.M[rr][r]/M.M[r][r],r,rr)
        return multiplier * reduce(lambda x,y: x*y,
                                   [M.M[r][r] for r in range(n)])
    def inverse(self) -> 'Matrix':
        if not self.is_square(): raise DimensionError()
        n = self.rows()
        M = self.clone()
        N = Matrix.identity(n) # perform same row operations to get inverse
        for r in range(n): # row reduce to M being an upper triangle matrix
            if M.M[r][r] == 0: # swap in a pivot row
                pr = -1 # pivot row
                for rr in range(r+1,n):
                    if M.M[rr][r] != 0:
                        pr = rr
                        break
                if pr == -1: raise NotInvertible() # cannot find pivot row
                else: # swap and negate multiplier
                    M.swap_rows(r,pr)
                    N.swap_rows(r,pr)
            row_mult = 1/M.M[r][r]
            M.mult_row(row_mult,r) # have leading element become 1
            N.mult_row(row_mult,r)
            for rr in range(r+1,n): # zero out column for upper triangle form
                row_mult = -M.M[rr][r]
                M.add_row_mult(row_mult,r,rr)
                N.add_row_mult(row_mult,r,rr)
        for r in range(n-1,-1,-1): # reduce upward to get M = I
            for rr in range(r-1,-1,-1):
                row_mult = -M.M[rr][r]
                M.add_row_mult(row_mult,r,rr)
                N.add_row_mult(row_mult,r,rr)
        return N
    def swap_rows(self,r1,r2):
        self.M[r1], self.M[r2] = self.M[r2], self.M[r1]
    def mult_row(self,a,r):
        if type(a) != int and type(a) != Fraction: raise TypeError()
        for c in range(self.cols()):
            self.M[r][c] *= a
    def add_row_mult(self,a,r1,r2):
        if type(a) != int and type(a) != Fraction: raise TypeError()
        for c in range(self.cols()):
            self.M[r2][c] += a * self.M[r1][c]
    def __str__(self):
        M = [list(map(str,r)) for r in self.M]
        col_len = [0]*self.cols()
        for r in M:
            for i,c in enumerate(r):
                col_len[i] = max(col_len[i],len(c))
        pad = lambda s, c : (col_len[c]-len(s))*' ' + s
        return '\n'.join('['+' '.join(pad(r[c],c)
                                      for c in range(self.cols()))+']'
                         for r in M)
    def clone(self):
        M = Matrix(1)
        M.M = [r[:] for r in self.M]
        return M
    def identity(n): # identity matrix
        M = Matrix(n)
        for i in range(n): M.set(i,i,1)
        return M
    def diagonal(e): # diagonal matrix with entries
        if type(e) != list: raise TypeError()
        M = Matrix(len(e))
        for i,n in enumerate(e):
            M.set(i,i,n)
        return M
    def scale(self,s): # scalar multiple
        if type(s) != int and type(s) != Fraction: raise TypeError()
        self.M = [[s * self.M[r][c] for c in range(self.cols())]
                  for r in range(self.rows())]
    def __add__(self,M):
        if self.rows() != M.rows() or self.cols() != M.cols():
            raise DimensionError()
        N = Matrix(1)
        N.M = [[self.M[r][c] + M.M[r][c] for c in range(self.cols())]
               for r in range(self.rows())]
    def __sub__(self,M):
        N = M.clone()
        N.scale(-1)
        return self.__add__(N)
    def __mul__(self,M):
        if self.cols() != M.rows(): raise DimensionError()
        N = Matrix(1)
        N.M = [[sum(self.M[r][i] * M.M[i][c] for i in range(self.cols()))
                for c in range(M.cols())]
               for r in range(self.rows())]
        return N
    def __pow__(self,p):
        if type(p) != int: raise TypeError()
        if not self.is_square(): raise DimensionError()
        n = self.rows()
        if self.is_diagonal():
            return Matrix.diagonal([self.M[r][r]**p for r in range(n)])
        if p == 0: return Matrix.identity(self.rows())
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
    A = Matrix(1)
    A.set(0,0,5)
    print(A)
    print('det =',A.det(), '(should be 5)')
    print('det =',A.det_fast(), '(should be 5)')
    A.M = [[1,2],
           [3,4]]
    _fix(A)
    print(A)
    print('det =',A.det(), '(should be -2)')
    print('det =',A.det_fast(), '(should be -2)')
    A.M = [[1,2,3],
           [2,-2,2],
           [-1,0,5]]
    _fix(A)
    print(A)
    print('det =',A.det(), '(should be -40)')
    print('det =',A.det_fast(), '(should be -40)')
    A.M = [[1,-1,1,-1],
           [5,0,6,7],
           [3,-1,-1,-1],
           [1,0,1,1]]
    _fix(A)
    print(A)
    print('det =',A.det(), '(should be 6)')
    print('det =',A.det_fast(), '(should be 6)')
    print('begin inverse')
    D = A.inverse()
    print(D)
    print(D*A)
    print('end inverse')
    A.M = [[Fraction(-2,3),Fraction(-1,4)],
         [Fraction(-1,5),Fraction(1,2)]]
    print(A)
    print('det =',A.det(), '(should be -23/60)')
    print('det =',A.det_fast(), '(should be -23/60)')
    A.M = [[1,-1,1,-1],
           [5,0,6,7],
           [3,-1,-1,-1],
           [-1,-2,-6,-9]] # last row = r1 - r2 + r3
    _fix(A)
    print(A)
    print('det =',A.det(), '(should be 0)')
    print('det =',A.det_fast(), '(should be 0)')
    try:
        Z = A.inverse()
        print(Z)
    except NotInvertible:
        print('cannot invert')
    print('identity')
    print(Matrix.identity(5))
    print('diagonal')
    print(Matrix.diagonal([1,-2,3,-4,5]))
    print('diagonal power')
    print(Matrix.diagonal([1,2,3])**-3)
    print('power')
    A.M = [[1,2],
           [3,4]]
    _fix(A)
    print(A**0)
    print(A**1)
    print(A**2)
    print(A**3)
    B = Matrix(2)
    B.M = [[7,10],
           [15,22]]
    _fix(B)
    print(A**2 == B)
    print('row operations')
    I = Matrix.identity(2)
    print(I)
    I.swap_rows(0,1)
    print(I)
    I.mult_row(2,1)
    print(I)
    I.mult_row(-3,0)
    print(I)
    I.add_row_mult(-1,1,0)
    print(I)
    I.add_row_mult(2,0,1)
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
    I = Matrix.identity(3)
    print(I.inverse())
    print(I*I == I)
    print(I*I.inverse())
