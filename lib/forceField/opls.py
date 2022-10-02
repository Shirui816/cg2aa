import os
import re
import warnings

from rdkit import Chem


def _warning(message, category=None, filename=None, lineno=None, file=None, line=None):
    print("WARNING: ", message)


warnings.showwarning = _warning

this_dir, this_filename = os.path.split(__file__)
opls_tpl_path = os.path.join(this_dir, "assets", "opls", "STaGE_opls_tomoltemplate_opls.txt")
opls_bonded_itp = os.path.join(this_dir, "assets", "opls", "bonded.itp")  # from gromacs


# lopls_tpl_path = 'lopls_tomoltemplate.txt'


def generate_feature_defn(fpath):
    info = {}
    feature_defn = ''
    with open(fpath, 'r') as infile:
        feat_index = 0
        for line in [line.strip() for line in infile if line.strip()]:
            if line[0] != '*':
                el, atomname, typename, patt, lttype, chg, desc = [el.strip() for el in line.split("|")]
                info[lttype] = (atomname, typename, chg, desc)
                if el == 'nom':
                    continue
                feature_defn += \
                    """
                    DefineFeature {0} {1}
                    Family {2}{3}
                    EndFeature""".format(lttype, patt, feat_index, atomname)
                # add new charge dictionary entry
                feat_index += 1
    return info, feature_defn


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


def ff(molecule, chem_envs_cache=None, chem_envs=None, large=500, hash_cut=3, **kwargs):
    info, feature_defn = generate_feature_defn(opls_tpl_path)
    factory = Chem.ChemicalFeatures.BuildFeatureFactoryFromString(feature_defn)
    radius = kwargs.get('radius') or 7
    if molecule.GetNumAtoms() > large:
        warnings.warn(f"Molecule too large, chemical env cut off for FF is set to {radius}, "
                      "cut off method may cause inaccuracy.")
        if radius < 7:
            warnings.warn(f"Chemical cutoff {radius} < 7, "
                          "please make sure fragment larger than the adj-Aromatic ring!")
        for _i, atom in enumerate(molecule.GetAtoms()):
            # print(_i, mol.GetNumAtoms())
            # check cache first
            atom_type, atom_env_hash = get_type_from_cache(molecule, atom, chem_envs_cache, hash_cut)
            if atom_type is not None:
                atom.SetProp('AtomType', atom_type)
                continue
            sub_mol, sub_mol_amap, sub_env = get_submol_rad_n(molecule, radius, atom)  # no need to sanitize
            for _atom in sub_mol.GetAtoms():
                if _atom.GetIsAromatic():
                    if not _atom.IsInRing():
                        _atom.SetIsAromatic(False)
            res = Chem.SanitizeMol(sub_mol, catchErrors=True)
            sub_smi = Chem.MolToSmiles(sub_mol, rootedAtAtom=sub_mol_amap[atom.GetIdx()], canonical=False)
            if res is not Chem.rdmolops.SanitizeFlags.SANITIZE_NONE:
                warnings.warn(
                    f"Molecule fragment {sub_smi} is not SANITIZED: {res}"
                    "this may cause errors, re-check env cutoff and molecular size."
                )
            features = factory.GetFeaturesForMol(sub_mol)
            [sub_mol.GetAtomWithIdx(f.GetAtomIds()[0]).SetProp('AtomType', f.GetType()) for f in features]

            # sub_mol_amap = {atom_idx: sub_mol_idx}
            try:
                sm_atom_id = sub_mol_amap.get(atom.GetIdx())
                if sm_atom_id is not None:
                    target_type = sub_mol.GetAtomWithIdx(sm_atom_id).GetProp('AtomType')
                    atom.SetProp('AtomType', target_type)
                    chem_envs_cache[atom.GetSymbol()][atom_env_hash] = target_type
            except KeyError:
                pass
    else:
        features = factory.GetFeaturesForMol(molecule)
        [molecule.GetAtomWithIdx(f.GetAtomIds()[0]).SetProp('AtomType', f.GetType()) for f in features]

    defaults = kwargs.get('defaults')
    fc = 0
    # check after auto-type
    for atom in molecule.GetAtoms():
        try:
            # print("Atom {0} has {1}".format(at.GetIdx(), at.GetProp('AtomType')))
            _ = atom.GetProp("AtomType")
            atom.SetIntProp("ff_failed", 0)
        except KeyError:
            fc += 1
            sub_mol, sub_mol_amap, sub_env = get_submol_rad_n(molecule, hash_cut, atom)  # default chem_env is 3
            chem_env = Chem.MolToSmiles(sub_mol, rootedAtAtom=sub_mol_amap[atom.GetIdx()], canonical=False)
            atom.SetIntProp("ff_failed", 1)
            _s = None
            if defaults is not None:
                if not defaults.get(atom.GetSymbol()) is None:
                    atom_type = defaults.get(atom.GetSymbol()).get(chem_env)
                    if atom_type is not None:
                        atom.SetProp("AtomType", atom_type)
                        _s = info[atom_type][0] + ", " + info[atom_type][3]
                        atom.SetIntProp("ff_failed", 0)  # not failed if there is default type
                        chem_envs_cache[atom.GetSymbol()][chem_env] = atom_type
            atom.SetProp("ChemEnv", chem_env)
            if chem_envs.get(atom.GetSymbol()) is None:
                chem_envs[atom.GetSymbol()] = {}
                chem_envs[atom.GetSymbol()][chem_env] = _s

    failed = False

    if failed:
        return None

    sum_of_charge = 0
    for atom in molecule.GetAtoms():
        names = info[atom.GetProp("AtomType")]
        et = names[0]
        if atom.GetIntProp("ff_failed") == 1:
            et = atom.GetProp("ChemEnv")
        atom.SetProp("ElementType", et)
        atom.SetProp("OplsType", names[1])
        atom.SetProp("Charge", names[2])
        sum_of_charge += float(names[2])
    if abs(sum_of_charge) > 0.01:
        warnings.warn("*** Total Charge is %.6f (%d atoms)!" % (sum_of_charge, molecule.GetNumAtoms()))

    bonded_itp = open(opls_bonded_itp, 'r').read()
    # TODO: make a reader and hash map, but for small systems re is fine
    bonds = []
    angles = []
    dihedrals = []
    for bond in molecule.GetBonds():  # bonds
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        atom_i = molecule.GetAtomWithIdx(i)
        atom_j = molecule.GetAtomWithIdx(j)

        bonds.append('%s-%s %d %d' % (atom_i.GetProp("ElementType"),
                                      atom_j.GetProp("ElementType"),
                                      i, j))

        for nbr_i in atom_i.GetNeighbors():  # a bond (i, j) and neighbor of i and j define a dihedral
            ni = nbr_i.GetIdx()
            if ni == j:
                continue
            for nbr_j in atom_j.GetNeighbors():
                nj = nbr_j.GetIdx()
                if nj == i or ni == nj:  # avoid 3-ring
                    continue
                dihedral = [nbr_i.GetProp("ElementType"),
                            atom_i.GetProp("ElementType"),
                            atom_j.GetProp("ElementType"),
                            nbr_j.GetProp("ElementType")]
                if re.search(r'\s*'.join(dihedral) + r'\s*[1-9]+', bonded_itp) or re.search(
                        r'\s*'.join(dihedral[::-1]) + r'\s*[1-9]+', bonded_itp):
                    dihedrals.append('%s %d %d %d %d' % ('-'.join(dihedral), ni, i, j, nj))

    for atom in molecule.GetAtoms():  # any atom can be a center of an angle
        for nbr_i in atom.GetNeighbors():
            for nbr_j in atom.GetNeighbors():
                if nbr_i.GetIdx() <= nbr_j.GetIdx():
                    continue
                angle = [nbr_i.GetProp("ElementType"), atom.GetProp("ElementType"), nbr_j.GetProp("ElementType")]
                if re.search(r'\s*'.join(angle) + r'\s*[1-9]+', bonded_itp) or re.search(
                        r'\s*'.join(angle[::-1]) + r'\s*[1-9]+', bonded_itp):
                    angles.append('%s %d %d %d' % ('-'.join(angle), nbr_i.GetIdx(), atom.GetIdx(), nbr_j.GetIdx()))

    return molecule, bonds, angles, dihedrals
