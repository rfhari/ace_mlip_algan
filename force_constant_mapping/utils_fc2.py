def process_files_fc2(file, neighbors):
    matrices, neighbor_index, lattice_disp = [], [], []
    with open(file, 'r') as f1: 
        for _ in range(3):  
            line = f1.readline()
        line1 = f1.readline()
        line2 = f1.readline()
        numbers = list(map(float, line2.strip().split()))
        formatted_numbers = [round(x, 15) for x in numbers]
        # print(line1.strip().split()[0], line1)
        neighbor_index.append(line1.strip().split()[0])
        lattice_disp.append(formatted_numbers)
        # print(formatted_numbers)
        while len(matrices) <= (4 * neighbors):
            matrix = []
            for _ in range(3):  
                line = f1.readline()
                if not line.strip():
                    break
                # print("matrix length:", len(matrices), "values:", list(map(float, line.strip().split()))) 
                matrix.append(list(map(float, line.strip().split())))
            if len(matrix) == 3:
                matrices.append(matrix)
            if len(matrices) >=500:
                break
            if len(matrices)%neighbors==0:
                f1.readline() 
                line1 = f1.readline()
                line2 = f1.readline()
                numbers = list(map(float, line2.strip().split()))
                formatted_numbers = [round(x, 15) for x in numbers]
                print(line1.strip().split()[0], line1)
                neighbor_index.append(line1.strip().split()[0])
                lattice_disp.append(formatted_numbers)
                # print(formatted_numbers)
            else:               
                line1 = f1.readline()
                line2 = f1.readline()
                numbers = list(map(float, line2.strip().split()))
                formatted_numbers = [round(x, 15) for x in numbers]
                print(line1.strip().split()[0], line1)
                neighbor_index.append(line1.strip().split()[0])
                lattice_disp.append(formatted_numbers)
                print(formatted_numbers)


    return matrices, lattice_disp, neighbor_index

