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


def generate_pos_fragment(molecule, meta):
    conf = FakeConf(molecule.GetNumAtoms())
    atom_map = {}
    for atom in molecule.GetAtoms():
        monomer_id = atom.GetIntProp('molecule_id')
        if not atom_map.get(monomer_id):
            atom_map[monomer_id] = []
        atom_map[monomer_id].append(atom.GetIdx())
    adj_dict = dict(meta.adjacency())
    for m_id in meta.nodes:
        # generate position of A by its neighbor monomers
        # to obtain better monomer-monomer connections
        n_ids = list(adj_dict[m_id].keys())
        fragment = Chem.RWMol()
        fragment_ids = [m_id] + n_ids
        bonds = set()
        local_map1 = {}
        local_map2 = {}
        count = 0
        broke = set()
        for monomer_id in fragment_ids:  # Get current monomer
            for atom_id in atom_map[monomer_id]:
                local_map1[atom_id] = count
                local_map2[count] = atom_id  # map monomer to polymer
                count += 1
        for i in range(count):
            atom_id = local_map2[i]
            atom = molecule.GetAtomWithIdx(atom_id)
            fragment.AddAtom(atom)
            for bond in atom.GetBonds():
                a, b, t = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType()
                if molecule.GetAtomWithIdx(a).GetIntProp("molecule_id") not in fragment_ids:
                    broke.add(atom_id)
                    continue
                if molecule.GetAtomWithIdx(b).GetIntProp("molecule_id") not in fragment_ids:
                    broke.add(atom_id)
                    continue
                if a < b:
                    bonds.add((a, b, t))
                else:
                    bonds.add((b, a, t))
        for bond in bonds:
            fragment.AddBond(local_map1[bond[0]], local_map1[bond[1]], bond[2])
        for atom in fragment.GetAtoms():  # avoiding breaking aromatic rings
            if atom.GetIsAromatic():
                if not atom.IsInRing():
                    atom.SetIsAromatic(0)
            if local_map2.get(atom.GetIdx()) in broke:
                atom.SetIsAromatic(0)
                for btom in atom.GetNeighbors():
                    btom.SetIsAromatic(0)
                    _bond = fragment.GetBondBetweenAtoms(atom.GetIdx(), btom.GetIdx())
                    _bond.SetBondType(Chem.rdchem.BondType.SINGLE)
        Chem.SanitizeMol(fragment, Chem.SanitizeFlags.SANITIZE_ADJUSTHS)
        res = Chem.SanitizeMol(fragment, catchErrors=True)
        if not res is Chem.rdmolops.SanitizeFlags.SANITIZE_NONE:
            fragment = Chem.MolFromSmiles(Chem.MolToSmiles(fragment))
            # I don't know why yet
            # but not important, for the truncated monomers are not used
            # raise ValueError(f"{res}, {Chem.MolToSmiles(fragment)}")
        _mh = AllChem.AddHs(fragment)
        AllChem.EmbedMolecule(_mh, useRandomCoords=True)
        _conf = _mh.GetConformer(0)
        for atom in _mh.GetAtoms():
            if local_map2.get(atom.GetIdx()) is None:
                continue
            m_a_id = local_map2.get(atom.GetIdx())
            p = _conf.GetAtomPosition(atom.GetIdx())
            if molecule.GetAtomWithIdx(m_a_id).GetIntProp("molecule_id") != m_id:
                continue
            conf.set_pos(m_a_id, Position(p.x, p.y, p.z))
    return conf
