"""
Perform Gaussian Elimation on the input matrix.

Elements of this matrix need to support +, *, /,
and identity checking with 0.

If the matrix is non-square, reduction is performed on as many
rows/columns as possible.

@param A
    Input matrix, of any dimension

@return
    (A', row_map), where A' is the reduced row echelon form
    of A, and row_map is a dictionary giving us the row permutations.
    In particular, row_map[i] = j means that A'[i] was originally
    in row j of A. This is useful for solving systems of linear equations.

"""
def gaussian_elimination(A):
    num_rows = len(A)
    num_cols = len(A[0])
    
    r = 0
    c = 0

    final_row_map_to_orig = dict([(i, i) for i in range(num_rows)])

    # forward scan, creating matrix in row echelon form
    while r < num_rows and c < num_cols:
        # if this column has a nonzero entry below row r,
        # row-operate so that everything below row 'r'
        # is 0 for column 'c'
        
        # find the nonzero entry
        pivot_row = r
        pivot = A[pivot_row][c]
        while pivot == 0 and pivot_row < num_rows:
            pivot = A[pivot_row][c]
            pivot_row += 1
        if pivot == 0:
            # column c is entirely of zeros; move to next col
            c += 1
        else:
            # swap pivot_row and r
            tmp = A[r]
            A[r] = A[pivot_row]
            A[pivot_row] = tmp

            tmp = final_row_map_to_orig[pivot_row]
            final_row_map_to_orig[pivot_row] = final_row_map_to_orig[r]
            final_row_map_to_orig[r] = tmp

            # reduce
            # for each i>r, subtract a multiplier of row r
            #   onto row i, making A[i][c] == 0
            pivot_elem = A[r][c]
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
        while i >= 0 and last_non_zero == 0:
            last_non_zero = A[i][c]
            i -= 1
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


    return A, final_row_map_to_orig