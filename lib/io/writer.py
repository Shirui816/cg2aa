import warnings


def _warning(message, category=None, filename=None, lineno=None, file=None, line=None):
    print("WARNING: ", message)


warnings.showwarning = _warning

from rdkit import Chem
from rdkit.Chem import AllChem

template = '''<?xml version ="1.0" encoding ="UTF-8" ?>
<{program}_xml version="{version}">
<configuration time_step="0" dimensions="3" natoms="{n_atoms:d}" >
<box lx="{lx:.8f}" ly="{ly:.8f}" lz="{lz:.8f}" xy="{xy:8f}" xz="{xz:8f}" yz="{yz:8f}"/>
<position num="{n_atoms:d}">
{positions}</position>
<type num="{n_atoms:d}">
{types}</type>
<opls_type num="{n_atoms:d}">
{opls_type}</opls_type>
<monomer_id num="{n_atoms:d}">
{monomer_id}</monomer_id>
<charge num="{n_atoms:d}">
{charge}</charge>
<mass num="{n_atoms:d}">
{mass}</mass>
<bond num="{n_bonds:d}">
{bond}
</bond>
<angle num="{n_angles:d}">
{angle}
</angle>
<dihedral num="{n_dihedrals:d}">
{dihedral}
</dihedral>
</configuration>
</{program}_xml>'''

LARGE = 500


class position(object):
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


class fake_conf():
    def __init__(self, num_atoms):
        self.x = {}

    def set_pos(self, i, pos):
        self.x[i] = pos

    def GetAtomPosition(self, idx):
        return self.x.get(idx)


def generate_pos_fragment(molecule):
    conf = fake_conf(molecule.GetNumAtoms())
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
            conf.set_pos(local_map2.get(atom.GetIdx()), position(p.x, p.y, p.z))
    return conf


def write_xml(molecule, box, bonds, angles, dihedrals, postfix, program='galamost', version='1.3'):
    if molecule.GetNumAtoms() > LARGE:
        warnings.warn(f"*** Num of atoms {molecule.GetNumAtoms()} is greater than {LARGE}, using fragment method.")
        conf = generate_pos_fragment(molecule)
        if len(conf.x) != molecule.GetNumAtoms():
            warnings.warn(f"*** configuration generation error! {len(conf.x)} != {molecule.GetNumAtoms()}")
    else:
        conf_id = AllChem.EmbedMolecule(molecule, useRandomCoords=True)
        if conf_id == -1:
            conf = generate_pos_fragment(molecule)
            if len(conf.x) != molecule.GetNumAtoms():
                warnings.warn(f"*** configuration generation error! {len(conf.x)} != {molecule.GetNumAtoms()}")
        else:
            conf = molecule.GetConformer(conf_id)
    n_atoms = molecule.GetNumAtoms()
    n_bonds = molecule.GetNumBonds()
    mass = types = opls_type = positions = charge = monomer_id = ''
    for atom in molecule.GetAtoms():
        mass += '%.6f\n' % atom.GetMass()
        types += '%s\n' % atom.GetProp("ElementType")
        opls_type += '%s\n' % atom.GetProp("OplsType")
        pos = conf.GetAtomPosition(atom.GetIdx())
        positions += '%.6f %.6f %.6f\n' % (pos.x, pos.y, pos.z)
        charge += '%.6f\n' % float(atom.GetProp("Charge"))
        monomer_id += '%d\n' % atom.GetIntProp("molecule_id")
    n_angles = len(angles)
    n_dihedrals = len(dihedrals)
    angle = '\n'.join(angles)
    dihedral = '\n'.join(dihedrals)
    bond = '\n'.join(bonds)
    lx, ly, lz, xy, xz, yz = box
    o = open('out_%s.xml' % postfix, 'w')
    o.write(
        template.format(
            n_atoms=n_atoms, n_bonds=n_bonds, mass=mass, types=types, opls_type=opls_type, positions=positions,
            bond=bond, charge=charge, angle=angle, dihedral=dihedral, n_angles=n_angles, n_dihedrals=n_dihedrals,
            monomer_id=monomer_id, program=program, version=version, lx=lx, ly=ly, lz=lz, xy=xy, xz=xz, yz=yz
        ))
    o.close()
    return
