from rdkit import Chem
from rdkit.Chem import AllChem


class Position(object):
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


class FakeConf():
    def __init__(self, num_atoms):
        self.x = {}

    def set_pos(self, i, pos):
        self.x[i] = pos

    def GetAtomPosition(self, idx):
        return self.x.get(idx)


def generate_pos_fragment(molecule):
    conf = FakeConf(molecule.GetNumAtoms())
    atom_map = {}
    for atom in molecule.GetAtoms():
        monomer_id = atom.GetIntProp('molecule_id')
        if not atom_map.get(monomer_id):
            atom_map[monomer_id] = []
        atom_map[monomer_id].append(atom.GetIdx())
    for monomer_id in atom_map:  # Get current monomer
        m = Chem.RWMol()
        bonds = set()
        local_map1 = {}
        local_map2 = {}
        count = 0
        for atom_id in atom_map[monomer_id]:
            local_map1[atom_id] = count
            local_map2[count] = atom_id  # map monomer to polymer
            count += 1
        for atom_id in atom_map[monomer_id]:
            atom = molecule.GetAtomWithIdx(atom_id)
            m.AddAtom(atom)
            for bond in atom.GetBonds():
                a, b, t = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType()
                if molecule.GetAtomWithIdx(a).GetIntProp("molecule_id") != monomer_id:
                    continue
                if molecule.GetAtomWithIdx(b).GetIntProp("molecule_id") != monomer_id:
                    continue
                if a < b:
                    bonds.add((a, b, t))
                else:
                    bonds.add((b, a, t))
        for bond in bonds:
            m.AddBond(local_map1[bond[0]], local_map1[bond[1]], bond[2])
        for atom in m.GetAtoms():  # avoiding breaking aromatic rings
            atom.SetIsAromatic(0)
        Chem.SanitizeMol(m)
        _mh = AllChem.AddHs(m)
        AllChem.EmbedMolecule(_mh, useRandomCoords=True)
        _conf = _mh.GetConformer(0)
        for atom in _mh.GetAtoms():
            if local_map2.get(atom.GetIdx()) is None:
                continue
            p = _conf.GetAtomPosition(atom.GetIdx())
            conf.set_pos(local_map2.get(atom.GetIdx()), Position(p.x, p.y, p.z))
    return conf
