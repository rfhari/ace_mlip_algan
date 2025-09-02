import numpy as np
from utils_fc2 import * 

temp = "300K"
fraction_al = 1.0
fraction_ga = 1 - fraction_al

aln_input_filename = "input_fcs/" + temp + "/infile_6.6_AlN.forceconstant"
gan_input_filename = "input_fcs/" + temp + "/infile_GaN.forceconstant"
vca_output_filename = f'vca_fcs/{temp}/vca_{fraction_al}Al_infile_6.6_AlN.forceconstant'

fc2_gan, lattice_disp_gan, neighbor_index_gan = process_files_fc2(gan_input_filename, neighbors=125)
fc2_aln, lattice_disp_aln, neighbor_index_aln = process_files_fc2(aln_input_filename, neighbors=125)

n_neighbors = 125
all_found, ordered_f2_gan = [], []

# -------------------- finding match pairs --------------------------

for iatom1 in range(1, 5):
    print("iatom1:", iatom1)
    for ii in range(len(lattice_disp_aln)):
        all_found_for_ii = []
        lower_bound_i, upper_bound_i = n_neighbors * (iatom1-1), n_neighbors * iatom1
        if lower_bound_i <= ii < upper_bound_i:
            if (neighbor_index_aln[ii] != neighbor_index_gan[ii]) or (lattice_disp_aln[ii] != lattice_disp_gan[ii]):
                print("Not an exact match:", ii)
                found_index = None
                for j, gan_disp in enumerate(lattice_disp_gan):
                    if (lower_bound_i <= j < upper_bound_i) and (gan_disp == lattice_disp_aln[ii]) and (neighbor_index_aln[ii] == neighbor_index_gan[j]):
                        found_index = j
                        all_found_for_ii.append(found_index)
                if len(all_found_for_ii)>1:
                    print("wrong matches:", ii, all_found_for_ii)
                else:
                    print("correct match, all_found_for_ii:", ii, all_found_for_ii)
                all_found.extend(all_found_for_ii)
                print("aln:", lattice_disp_aln[ii], neighbor_index_aln[ii])
                if found_index is not None:
                    ordered_f2_gan.append(fc2_gan[found_index])
                    print("gan:", lattice_disp_gan[found_index], neighbor_index_gan[found_index])
                    print(f"Match found after search: AlN index {ii} matches GaN index {found_index}")
                else:
                    print("No matching GaN lattice_disp found for this AlN entry.")
            elif (neighbor_index_aln[ii] == neighbor_index_gan[ii]) and (lattice_disp_aln[ii] == lattice_disp_gan[ii]):
                if ii not in all_found:
                    all_found.extend([ii])
                    ordered_f2_gan.append(fc2_gan[ii])
                    print(f"Exact Match found without search: AlN index {ii} matches GaN index {ii}")
                    print("gan:", lattice_disp_gan[ii], neighbor_index_gan[ii])
                    print("aln:", lattice_disp_aln[ii], neighbor_index_aln[ii])
                else:
                    print("Conflict in matching")
                    break
    print("len of ordered fc2 gan:", len(ordered_f2_gan))


# -------------------- calculating interpolated values -------------------------- 
vca_fc2 = []

for ii in range(len(fc2_aln)):
    vca_fc2.append(np.asarray(fc2_aln[ii]) * fraction_al + np.asarray(ordered_f2_gan[ii]) * fraction_ga)

# -------------------- writing interpolated values --------------------------
ptr = 0
reminder_ind = 0
line_to_read = [0, 1, 2]

with open(gan_input_filename, 'r') as infile:
    gan_lines = infile.readlines()

with open(aln_input_filename, 'r') as infile, open(vca_output_filename, 'w') as outfile:
    for lineno, line in enumerate(infile, start=1):
        if lineno in [629, 1255, 1881]:
            reminder_ind+=1
            if reminder_ind==1:
                line_to_read = [1, 2, 3]
            elif reminder_ind==2:
                line_to_read = [2, 3, 4]
            elif reminder_ind==3:
                line_to_read = [3, 4, 0]
        print("line_to_read:", line_to_read)
        if lineno <= 2507:
            if lineno >= 6 and ((lineno - 6) % 5) in line_to_read and lineno not in [629, 630, 1255, 1266, 1881, 1882]: # hard coded 
                print("line number:", lineno, "line:", line)
                if reminder_ind==0:
                    index = (lineno - 6) % 5
                elif (reminder_ind==3 and ((lineno - 6) % 5)==0):
                    index = 2 
                else:
                    index = (lineno - 6) % 5 - reminder_ind
                print("line number:", lineno, "ptr:", ptr, "index:", index, "reminder_ind:", reminder_ind)
                formatted = " ".join(f"{item:8.15f}" for item in vca_fc2[ptr][index])
                modified_line = "   " + formatted + "\n"
                outfile.write(modified_line)
                if index==2:
                    print("changing current ptr:", ptr)
                    ptr = ptr + 1
            else:
                outfile.write(line)
        
        # interpolating born charge based second-order-force-constants (hard coded line no. currently)
        elif 2508 <= lineno <= 2510:
            dieletric_value = list(map(str, line.strip().split()[3:])) 
            aln_value = list(map(float, line.strip().split()[:3]))
            gan_value = list(map(float, gan_lines[lineno-1].strip().split()[:3]))
            average_value = np.asarray(aln_value) * fraction_al + np.asarray(gan_value) * fraction_ga
            print('dielectric:', dieletric_value, 'aln/gan values:', aln_value, gan_value)
            value_to_write = list(average_value) + dieletric_value
            outfile.write("       " + " ".join(map(str, value_to_write)) + "\n")
            
        elif lineno == 2511:
            text_v = list(map(str, line.strip().split()[1:]))
            v1 = float(line.strip().split()[0])
            v2 = float(gan_lines[lineno-1].strip().split()[0])
            print("coupling values:", v1, v2)
            to_write = [str(v1 * fraction_al + v2 * fraction_ga)] + text_v
            outfile.write("  " + " ".join(map(str, to_write)) + "\n")
        elif 2513 <= lineno <= 2539:
            v1 = float(line.strip().split()[0])
            v2 = float(gan_lines[lineno-1].strip().split()[0])
            print("irreps of Born Charges:", v1, v2)
            value = v1 * fraction_al + v2 * fraction_ga
            to_write = f"{value:.15f}"
            outfile.write("       " + to_write + "\n")
        elif 2576 <= lineno <= 2587:
            try:
                text_value = list(map(str, line.strip().split()[3:])) 
            except:
                text_value = []
            aln_value = list(map(float, line.strip().split()[:3]))
            gan_value = list(map(float, gan_lines[lineno-1].strip().split()[:3]))
            average_value = np.asarray(aln_value) * fraction_al + np.asarray(gan_value) * fraction_ga
            average_value = [f"{element:.15f}" for element in average_value]
            print('born charge:', text_value, 'aln/gan values:', aln_value, gan_value)
            value_to_write = average_value + text_value
            outfile.write("       " + " ".join(map(str, value_to_write)) + "\n")
        else:
            outfile.write(line)
        