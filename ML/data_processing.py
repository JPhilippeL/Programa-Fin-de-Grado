from rdkit import Chem
import os
import numpy as np
from molgraph.chemistry import MolecularGraphEncoder, AtomFeaturizer


atom_encoder = AtomFeaturizer()
encoder = MolecularGraphEncoder(atom_encoder)

def load_sdf_directory_to_molgraphs(sdf_dir):
    molgraphs, filenames = [], []
    for fname in sorted(os.listdir(sdf_dir)):
        if fname.endswith(".sdf"):
            mol = Chem.MolFromMolFile(os.path.join(sdf_dir, fname), sanitize=True)
            if mol:
                molgraphs.append(encoder(mol))
                filenames.append(os.path.splitext(fname)[0])
    return molgraphs, filenames



def load_target_file(target_file, molecule_names):
    target_map = {}
    with open(target_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                target_map[parts[0]] = float(parts[1])

    targets = []
    for name in molecule_names:
        targets.append(target_map.get(name, 0.0))  # valor por defecto si falta
    return np.array(targets)
