import numpy as np

def process_files_fc3(file, neighbors):
    all_fc3_matrices, all_lattice_disp, all_iatom_indices = [], [], []
    with open(file, 'r') as f1: 
        for _ in range(3):  
            line = f1.readline()
        while len(all_fc3_matrices) < 4 * neighbors:
            fc3_matrix, lattice_disp, iatom_indices = [], [], []
            # to store iatom indices
            for _ in range(3):
                line = f1.readline()
                if not line.strip():
                    break
                # print("all_iatoms length:", len(iatom_indices), "values:", line) 
                iatom_indices.extend(list(map(float, line.strip().split()[0])))
            if len(iatom_indices) == 3:
                all_iatom_indices.append(iatom_indices)                
            # to store lattice displacement vectors 
            for _ in range(3):  
                line = f1.readline()
                if not line.strip():
                    break
                # print("all_lattice_disp length:", len(lattice_disp), "values:", list(map(float, line.strip().split()))) 
                lattice_disp.append(list(map(float, line.strip().split())))
            if len(lattice_disp) == 3:
                all_lattice_disp.append(lattice_disp)
            # to store fc3            
            for _ in range(9):  
                line = f1.readline()
                if not line.strip():
                    break
                # print("all_fc3_matrices length:", len(fc3_matrix), "values:", list(map(float, line.strip().split()))) 
                fc3_matrix.append(list(map(float, line.strip().split())))
            if len(fc3_matrix) == 9:
                all_fc3_matrices.append(fc3_matrix)  

            if len(all_fc3_matrices)%neighbors==0:
                # print("within if")
                for _ in range(1):  
                    print(f1.readline())
            if len(all_fc3_matrices) >= 4*neighbors:
                break

    return all_fc3_matrices, all_lattice_disp, all_iatom_indices

