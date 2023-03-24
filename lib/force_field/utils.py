from rdkit import Chem


def get_submol_rad_n(mol, radius, atom):
    env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom.GetIdx(), useHs=True)
    amap = {}
    sub_mol = Chem.PathToSubmol(mol, env, atomMap=amap)
    return sub_mol, amap, env


def get_type_from_cache(molecule, atom, chem_envs_cache, hash_cut):
    sub_mol, sub_mol_amap, sub_env = get_submol_rad_n(molecule, hash_cut, atom)  # 3 is good enough for hash
    atom_env_hash = Chem.MolToSmiles(sub_mol, rootedAtAtom=sub_mol_amap[atom.GetIdx()], canonical=False)
    if chem_envs_cache.get(atom.GetSymbol()) is None:
        chem_envs_cache[atom.GetSymbol()] = {}
    atom_type = chem_envs_cache[atom.GetSymbol()].get(atom_env_hash)
    return atom_type, atom_env_hash
