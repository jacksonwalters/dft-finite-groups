import os
from sage.matrix.constructor import Matrix

def load_csv_as_matrix(filename):
    """
    load a .csv file which is a matrix of integers
    """
    import csv
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile)
        matrix = [list(map(int, row)) for row in reader]
    return Matrix(matrix)

def write_matrix_to_csv(P,filename):
    """
    dump the eigenvector matrix into a .csv file
    """
    with open(filename, "w") as f:
        for row in P:
            f.write(",".join(map(str, row)) + "\n")

def save_array(array, path, filename):
    """
    save an array as a comma separated value file by converting the elements to strings
    """
    if not os.path.exists(path):
        os.makedirs(path)  # Create the directory if it doesn't exist
    full_path = os.path.join(path, filename)
    with open(full_path, 'w') as f:
        for i, element in enumerate(array):
            f.write(str(element))
            if i < len(array) - 1:
                f.write(",")
        f.write("\n")

def load_array(path, filename):
    """
    Load an array from a comma-separated value file and convert elements to strings.
    """
    full_path = os.path.join(path, filename)
    with open(full_path, 'r') as f:
        line = f.readline().strip()
        return line.split(",")
    
#given a list of eigenvalues in string representation, convert to elements of algebraic closure
import re
def str_repn_to_alg_closure(F_bar, eigenvalues_F_bar_list):
    eigenvalues_F_bar = []
    for eig in eigenvalues_F_bar_list:
        monomials = eig.split('+')
        pattern = r'(\d*)\*?z(\d+)(\^(\d+))?|\b(\d+)\b'
        eig_F_bar = 0
        for monomial in monomials:
            match = re.match(pattern, monomial.strip())
            if match:
                if match.group(5):  # Case for constants (no 'z' part)
                    coefficient = int(match.group(5))
                    variable_number = 1
                    exponent = 0
                else:  # Case for terms with 'z' and potentially an exponent
                    coefficient = int(match.group(1)) if match.group(1) else 1  # Default to 1 if no coefficient
                    variable_number = int(match.group(2))  # The integer after 'z'
                    exponent = int(match.group(4)) if match.group(4) else 1  # Default to 1 if no exponent
                monomial_F_bar = F_bar(coefficient)*F_bar.gen(variable_number)**exponent
                eig_F_bar += monomial_F_bar
            else:
                print(f"Monomial: {monomial.strip()} - No match found")
        eigenvalues_F_bar.append(eig_F_bar)
    return eigenvalues_F_bar