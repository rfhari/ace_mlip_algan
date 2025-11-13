from pyace import PyACECalculator
import ase
from ase.io import read, write
import glob
import numpy as np
from phonopy.interface.calculator import read_crystal_structure
from phono3py import Phono3py
from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5
from tqdm import tqdm
from ase import Atoms
from phonopy.file_IO import parse_BORN

# temp = str(300)
# unitcell = read_crystal_structure("new_minimized_aln_"+temp+"K.unitcell", interface_mode='lammps')
unitcell = read_crystal_structure('minimized_aln.unitcell', interface_mode='lammps')
print(unitcell[0])


ph3 = Phono3py(
    unitcell[0],
    supercell_matrix=[3,3,2],
    phonon_supercell_matrix=[6,6,4],
)

ph3.generate_displacements(cutoff_pair_distance=2.5)
ph3.save("phono3py_disp.yaml")
print("Generated displacements")

calc = PyACECalculator('hari_AlGaN_v0.yaml')
forces=[]
# calc FC2
for sc in tqdm(ph3.phonon_supercells_with_displacements):
    atoms = Atoms(sc.symbols, cell=sc.cell, positions=sc.positions, pbc=True)
    atoms.set_calculator(calc)
    f = atoms.get_forces()
    forces.append(f)

ph3.phonon_forces = np.array(forces)    

nac_params = parse_BORN(ph3.primitive, filename="BORN")
ph3.nac_params = nac_params


ph3.produce_fc2()
write_fc2_to_hdf5(
    ph3.fc2,
    p2s_map=ph3.phonon_primitive.p2s_map,
    physical_unit="eV/angstrom^2",
)

# # convert fc3

print("Computing forces for FC3")

forces = []
nat = len(ph3.supercells_with_displacements[0])
for sc in tqdm(ph3.supercells_with_displacements):
    if sc is not None:
        atoms = Atoms(sc.symbols, cell=sc.cell, positions=sc.positions, pbc=True)
        atoms.set_calculator(calc)
        f = atoms.get_forces()
    forces.append(f)

# compute 3rd order force constants
ph3.forces = np.array(forces)

nac_params = parse_BORN(ph3.primitive, filename="BORN")
ph3.nac_params = nac_params


print("Computing FC3")
ph3.produce_fc3()
write_fc3_to_hdf5(
    ph3.fc3,
    p2s_map=ph3.primitive.p2s_map,
)
