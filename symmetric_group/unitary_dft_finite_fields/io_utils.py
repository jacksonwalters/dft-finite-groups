import os

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