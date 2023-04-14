import numpy
import numpy as np

"""
1.converting to row echelon form
2.converting to reduced row echelon form
    a.making all the leading non-zero numbers equal to 1
    b.making the columns containing the leading coefficients to be equal to zero except the leading coefficient itself
"""


class LinearEqnSolver:
    def __init__(self, equations, r, c):
        self.equations = equations
        self.r = r
        self.c = c

    def swap_rows(self, r1: int, r2: int):
        # assigning the sub-matrix of the 2 rows to the reverse order or it
        self.equations[[r1, r2], :] = self.equations[[r2, r1], :]
        print(f"R{r1} <-> R{r2}")
        self.display()

    def sort(self):
        """
        sorting the matrix so that the rows with more leading 0s go down
        :return: None
        """
        for h in range(self.r - 1, -1, -1):
            most_leading_zeroes_index = 0
            most_leading_zeroes = 0

            for i in range(0, h + 1):
                num_leading_zeroes = 0
                # let's count the number of leading zeroes
                for j in range(0, self.c - 1):
                    # if we find an element that isn't a zero, we stop incrementing the number of leading zeroes
                    if self.equations[i][j] != 0:
                        break
                    else:
                        num_leading_zeroes += 1

                if num_leading_zeroes > most_leading_zeroes:
                    most_leading_zeroes = num_leading_zeroes
                    most_leading_zeroes_index = i

            if most_leading_zeroes > 0:
                if h != most_leading_zeroes_index:  # swapping with same row is redundant
                    self.swap_rows(h, most_leading_zeroes_index)
            else:  # if the most leading zeroes in any row is not more than 0, then we don't need to sort
                break

    def check_contradiction(self):
        """
        checking whether there are rows with all zeroes except the last element-->which is shows that the system is contradicgtory
        "No Solutions" if there is a contradiction
        :return: True if there is a contradiction, else False
        """
        for i in range(0, self.r):
            j = 0
            while j < self.c - 1:
                if self.equations[i][j] != 0:
                    break
                j += 1
            # do we have all zero elements in this row, except the last element?
            if j == self.c - 1 and self.equations[i][j] != 0:
                return True

        return False

    def check_sig(self):
        """
         edge case: not enough info i.e., number of significant coefficients in each row > number of rows
         which means you can't find the values of all the variables
         "Infinitely many solutions"
        :return: True if (more vars > number of rows with all coeff=0) else False
        """
        non_sig = 0  # variable to keep track of the number of rows with all coeff=0
        for i in range(0, self.r):
            j = 0
            while j < self.c - 1:
                if self.equations[i][j] != 0:
                    break
                j += 1
            # do we have all coeff=0 in this row, then j would end on (c-1)
            if j == self.c - 1:
                non_sig += 1

        if (self.r - non_sig) < (self.c - 1):
            print("sig<vars Infinite solutions")
            return True

    def display(self):
        print(np.array_str(self.equations))  # , precision=3, suppress_small=True))

    def add(self, place: int, r1: int, k: float, r2: int):
        for j in range(0, self.c):
            self.equations[place][j] = self.equations[r1][j] + k * self.equations[r2][j]

    def mult(self, r1: int, k: float):
        for j in range(0, self.c):
            self.equations[r1][j] = self.equations[r1][j] * k

    def solve(self):
        self.display()

        # 1
        for i in range(0, self.r - 1):
            self.sort()
            for j in range(0, self.c - 1):
                # is variable != 0? then make everything below this leading coefficient equal to zero
                if self.equations[i][j] != 0:
                    for I in range(i + 1, self.r):
                        if self.equations[I][j] != 0:  # it has found a number !=0 below it
                            # rowI = rowI + [rowi * (1 / element below it) * -leading coefficient(of i))]
                            print(f"R{I} = R{I} + (-{self.equations[I][j]} /{self.equations[i][j]})*R{i}")
                            self.add(I, I, -1 * self.equations[I][j] / self.equations[i][j], i)
                            self.display()
                    break  # break from the inner loop

        self.sort()
        self.check_contradiction()
        self.check_sig()
        print("1. converted to row echelon form\n")

        # 2.a
        for i in range(0, self.r):
            for j in range(0, self.c - 1):
                if self.equations[i][j] != 0:  # checking for non-zero leading element
                    if self.equations[i][j] != 1:  # if leading element equal to 1, we don't do anything to that row
                        print(f"R{i}=(1/{self.equations[i][j]})*R{i}")
                        self.mult(i, 1.0 / self.equations[i][j])  # we want the leading element to equal 1
                        self.display()
                    break

        print("\n 2.a all non-zero leading numbers = 1")
        self.display()

        # 2.b
        for i in range(1, self.r):
            for j in range(0, self.c - 1):
                if self.equations[i][j] == 1:  # if the current element is a leading 1
                    for I in range(i - 1, -1, -1):  # we go on making the elements above equal to zero
                        print(f"R{I} = R{I} + (-{self.equations[I][j]} /{self.equations[i][j]})*R{i}")
                        self.add(I, I, self.equations[I][j] * -1, i)
                        self.display()

                    break


"""
p+q+2r+s+2t=15
p+4q-r+2s+t=8
3p-2q+r+3s+4t=22
-p-q+3r+s+2t=5
2p+3q+4r+4s-t=30
"""
eqns = np.array([
    [1, 1, 2, 1, 2, 15],
    [1, 4, -1, 2, 1, 8],
    [3, -2, 1, 3, 4, 22],
    [-1, -1, 3, 1, 2, 5],
    [2, 3, 4, 4, -1, 30]
])  # only pass in float values
a = LinearEqnSolver(eqns, 5, 6)
a.solve()
