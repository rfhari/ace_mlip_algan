import numpy as np
from utils_fc3 import *

temp = "300K"
fraction_al = 0.01
fraction_ga = 1 - fraction_al

input_aln_filename = "input_fcs/" + temp + "/infile_6.6_AlN.forceconstant_thirdorder"
intput_gan_filename = "input_fcs/" + temp + "/infile_GaN.forceconstant_thirdorder"
vca_output_filename = f'vca_fcs/{temp}/vca_{fraction_al}_infile_6.6_AlN.forceconstant_thirdorder'

gan_fc3_matrices, gan_lattice_disp, gan_iatom_indices = process_files_fc3(intput_gan_filename, neighbors=361)
aln_fc3_matrices, aln_lattice_disp, aln_iatom_indices = process_files_fc3(input_aln_filename, neighbors=361)

# get matched pairs 
all_matches = []
n_neighbors = 361
for ii in range (len(aln_iatom_indices)):
    match_indices = []
    index_start = (ii // n_neighbors) * n_neighbors + 0
    index_end = (ii // n_neighbors) * n_neighbors + n_neighbors
    for jj in range (len(gan_iatom_indices)):
        if index_start <= jj < index_end and (aln_iatom_indices[ii] == gan_iatom_indices [jj]) and (aln_lattice_disp[ii] == gan_lattice_disp[jj]):
            match_indices.append(jj)
        if len(match_indices) > 1:
            print("wrong matches:", ii, match_indices)
    print("ii:", ii, "matches:", match_indices, "search range:", index_start, index_end)
    all_matches.extend(match_indices)

# --------- verify if there are any differences ----------------- 
# for ii in range(len(all_matches)):
#     l1 = np.asarray(gan_iatom_indices[all_matches[ii]])
#     l2 = np.asarray(aln_iatom_indices[ii])
#     d1 = np.asarray(gan_lattice_disp[all_matches[ii]])
#     d2 = np.asarray(aln_lattice_disp[ii])
#     difference1 = l1 - l2
#     difference2 = d1 - d2
#     if not np.allclose(difference1, [0, 0, 0]) or not np.allclose(d1, d2):
#         print("not matched index:", ii, "aln index:", l2, "gan index:", l1, "aln ld:", d2, "gan ld:", d1)
#     else:
#         print("matched index:", ii, "difference:", difference1, difference2)

# --------- vca FC3 calculations -----------------
vca_fc3 = []

for ii in range(len(aln_fc3_matrices)):
    vca_fc3.append(np.asarray(aln_fc3_matrices[ii]) * fraction_al + np.asarray(gan_fc3_matrices[all_matches[ii]]) * fraction_ga)

# ---------- write to fc3 format -----------------

ptr = 0
extra_line_positions = {5419, 10835, 16251}  
index = 0
offset = 0  

with open(input_aln_filename, "r") as infile, open(vca_output_filename, "w") as outfile:
    for orig_line_no, line in enumerate(infile, start=1):
        if orig_line_no in extra_line_positions:
            offset += 1
        normalized_line_no = orig_line_no - offset
        if normalized_line_no >= 10 and ((normalized_line_no - 10) % 15) < 9 and orig_line_no not in extra_line_positions:
            # print("lineno:", orig_line_no, "ptr:", ptr, "index:", index)
            formatted = " ".join(f"{item:8.15f}" for item in vca_fc3[ptr][index])
            modified_line = "   " + formatted + "\n"
            outfile.write(modified_line)
            if index<8:
                index+=1
            elif index==8:
                index=0
                ptr+=1
        else:
            # print("lineno:", orig_line_no, "line:", line)
            outfile.write(line)