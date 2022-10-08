import warnings

from rdkit.Chem import AllChem

from lib.io.utils import generate_pos_fragment


def _warning(message, category=None, filename=None, lineno=None, file=None, line=None):
    print("WARNING: ", message)


warnings.showwarning = _warning

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
