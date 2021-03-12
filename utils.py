"""
Perform Gaussian Elimation on the input matrix.

Elements of this matrix need to support +, *, /,
and identity checking with 0.

If the matrix is non-square, reduction is performed on as many
rows/columns as possible.

@param A
    Input matrix, of any dimension

@return
    A', input matrix transformed to reduced row echelon form

"""
def gaussian_elimination(A):
    num_rows = len(A)
    num_cols = len(A[0])
    
    r = 0
    c = 0

    # forward scan, creating matrix in row echelon form
    while r < num_rows and c < num_cols:
        # if this column has a nonzero entry below row r,
        # row-operate so that everything below row 'r'
        # is 0 for column 'c'
        
        # find the nonzero entry
        pivot_row = r
        pivot = A[pivot_row][c]
        while pivot == 0 and pivot_row < num_rows - 1:
            pivot_row += 1
            pivot = A[pivot_row][c]
        if pivot == 0:
            # column c is entirely of zeros; move to next col
            c += 1
        else:
            # swap pivot_row and r
            tmp = A[r]
            A[r] = A[pivot_row]
            A[pivot_row] = tmp

            # reduce
            # for each i>r, subtract a multiplier of row r
            #   onto row i, making A[i][c] == 0
            for i in range(r+1, num_rows):
                if A[i][c] == 0:
                    continue

                multiplier = A[i][c] / A[r][c]
                for j in range(c, num_cols):
                    A[i][j] -= multiplier * A[r][j]

            # make row r's leading nonzero entry == 1
            divisor = A[r][c]
            for j in range(c, num_cols):
                A[r][j] /= divisor

            c += 1
            r += 1


    # backward scan, reducing to unique reduced row echelon form
    # for each column, find the last row which is non-zero, and
    # then zero-out all the earlier rows
    # we can start at (r, c) from before
    # every column from last to first, the last non-zero row must be strictly decreasing
    r -= 1
    c -= 1
    while r >= 0 and c >= 0:
        i = r
        last_non_zero = A[i][c]
        while i > 0 and last_non_zero == 0:
            i -= 1
            last_non_zero = A[i][c]
        if last_non_zero == 0:
            c -= 1  # go to next column
        else:
            for j in range(i):
                # zero out A[j][c] by subtracting row i for all previous rows
                # row i only has non-zero entries after column c
                multiplier = A[j][c] / A[i][c]
                for k in range(c, num_cols):
                    A[j][k] -= multiplier * A[i][k]
                
            r -= 1
            c -= 1


    return A


"""
Solves the exact linear system Ax = b using Gaussian Elimination.
A should be a n-by-n square matrix and b a 1-by-n vector.

Throws ValueError if there's no solution or if there are multiple solutions.
"""
def solve_lin_sys(A, b):
    num_rows = len(A)
    if len(b) != num_rows:
        raise ValueError(f"Unmatched dimension in linear system: {num_rows} and {len(b)}")

    if len(A[0]) != num_rows:
        raise ValueError("Expected A to be a square matrix.")

    A_aug = [A[i] + [b[i]] for i in range(num_rows)]

    A_reduced = gaussian_elimination(A_aug)

    solution = [None] * num_rows

    for i in range(num_rows):
        # if there's an unique solution, A_reduced[i][i] should be 1
        if A_reduced[i][i] == 0:
            if b[i] == 0:
                raise ValueError("A is not full rank, no unique solution.")
            else:
                raise ValueError("No solution exist for this system.")
        if A_reduced[i][i] != 1:
            raise ValueError("Unknown error happened, corrupted result from "
                             + "Gaussian Elimination.")

        solution[i] = A_reduced[i][-1]

    return solution


'''
Divides f by g.
Solves the equation
    f(x) = g(x)q(x) + r(x)
for some r with degree strictly less than deg(g).

@param f
    divident polynomial, represented as a list of coefficients (lower degree first)
@param g
    divisor polynomial, same format as f

@return
    (q, r), where q is the quotient polynomial
    and r is the reminder, same formats as f
'''
def poly_div(f, g):
    n = len(f)
    k = len(g)

    reminder = f
    quotient = []

    # each loop iter correspond to one long division step
    for x_pow in range(n-k, -1, -1):
        quotient_term = reminder[-1] / g[-1]
        quotient.append(quotient_term)
        # multiply g by divisor_term * x^{x_pow}, then
        # subtract from reminder to get new reminder
        for i in range(k):
            reminder[x_pow + i] -= quotient_term * g[i]

        # last term of reminder should be 0 now; trim
        reminder = reminder[:-1]

    # long division gets back quotient with leading coefficient first
    quotient = quotient[::-1]

    # trim zeros in reminder
    while len(reminder) > 0 and reminder[-1] == 0:
        reminder = reminder[:-1]

    return (quotient, reminder)


''' Evaluate f(x). '''
def poly_eval(f, x):
    return sum([f[i] * x**i for i in range(len(f))])


''' Given a list of n points x, y,
finds the unique polynomial f of degree n-1 
such that f(x_i) = y_i '''
def poly_interpolate(X, Y):
    n = len(X)
    if len(Y) != n:
        raise ValueError("Expected X, Y to be of the same length.")

    # f(x_i) = y_i gives us a set of n linear equations
    A = [[X[i]**j for j in range(n)] for i in range(n)]
    b = Y
    f = solve_lin_sys(A, b)
    return f
