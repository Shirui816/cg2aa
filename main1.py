__doc__ = r"""
An example for epoxy.
"""

import os.path
from concurrent.futures import ProcessPoolExecutor as Executor  # cpu-bound jobs
# from concurrent.futures import ThreadPoolExecutor as Executor  # io-bound
from multiprocessing import Manager
from sys import argv

from rdkit import Chem

from lib.cg_topology import read_cg_topology
from lib.force_field.opls import ff
from lib.force_field.printMissingType import print_missing_type
from lib.io.xml_parser import XmlParser
from lib.io.xml_writer import write_xml
from lib.reactor import Reactor
from lib.utils import set_molecule_id_for_h

# ctl = ...
# get molecules from control.in

# Read this part from file
molecules = {
    'M0': {'smiles': 'N#COc1ccc(cc1)C(C)(C)c2ccc(OC#N)cc2', 'pdb': None},
    'M1': {'smiles': 'O1CC1COc2ccc(cc2)C(C)(C)c3ccc(cc3)OCC(O)COc4ccc(cc4)C(C)(C)c5ccc(cc5)OCC6CO6', 'pdb': None},
    # E.g. N=NC=O (not real molecule) for testing missing types
}

reaction_template = {
    'tri': {
        'cg_reactant_list': [('M0', 'M0', 'M0')],
        'smarts': "[C:1]#[N:2].[C:3]#[N:4].[C:5]#[N:6]>>[c:1]1[n:2][c:3][n:4][c:5][n:6]1",
        'prod_idx': [0]
    },
    'di': {
        'cg_reactant_list': [('M0', 'M1')],
        'smarts': "[#6:8][O:1][C:2]#[N:3].[#6:4][C:5]1[C:6][O:7]1>>[#6:4][C:5]1[C:6][N:3]([#6:8])[C:2](=[O:1])[O:7]1",
        'prod_idx': [0]
    }
}

defaults = {
    'O': {
        "O(c(c(c)[H])c(c)[H])C#N": "@atom:408",
    }
}

# End


# this part is for ANY cases
for key in molecules:
    if not molecules[key]['pdb'] is None:
        mol = Chem.RemoveAllHs(Chem.MolFromPDBFile(molecules[key]['pdb']))
        molecules[key]['smiles'] = Chem.MolToSmiles(mol)


# for key in molecules:
#     mol = Chem.MolFromSmiles(molecules[key]['smiles'])
#     mol = Chem.AddHs(mol)
#     AllChem.EmbedMolecule(mol)
#     pdb = Chem.MolToPDBFile(mol, '%s.pdb' % key)


def processing(i, mol, cg_mol, box, mt, ch, meta, defaults, radius=7):
    f_name = f"out_{i:06d}.xml"
    if os.path.isfile(f_name):  # skip existing files
        return 2
    ret = ff(mol, chem_envs=mt, chem_envs_cache=ch, large=500, defaults=defaults, radius=radius)
    if ret is not None:
        plm, bonds, angles, dihedrals, impropers = ret
        for m in meta.nodes:  # set monomer id
            molecule = meta.nodes[m]
            for idx in molecule['atom_idx'].values():
                atom = plm.GetAtomWithIdx(idx)
                atom.SetIntProp('molecule_id', int(m))
        plm = set_molecule_id_for_h(plm)
        write_xml(plm, cg_mol, box, bonds, angles, dihedrals, impropers, '%06d' % i)
        return 1
    return 0
    # mol = Chem.AddHs(mol)
    # AllChem.EmbedMolecule(mol)
    # pdb = Chem.MolToPDBFile(mol, '%06d.pdb' % i)


# Serial ver.
# def main(mols, box, meta, default_types=None):
#     cache = dict()
#     missing_types = dict()
#     res = {}
#     for i, molecule in enumerate(mols):
#         args = (i, molecule, box, missing_types, cache, meta[i], default_types)
#         res[i] = processing(*args)
#     print_missing_type(missing_types)
#     for i in res:
#         if res[i].result() == 0:
#             print(f"Molecule {i} is not generated, please re-check force field parameters.")


def main(mols, box, cg_mols, meta, default_types=None, draw=False):
    all_elements = set()
    for _mol in mols:
        for atom in _mol.GetAtoms():
            all_elements.add(atom.GetSymbol())
    m = Manager()  # worry not, it's atomic
    cache = m.dict()
    missing_types = m.dict()
    for ele in all_elements:
        cache[ele] = m.dict()
        missing_types[ele] = m.dict()
    futures = {}
    with Executor() as e:
        for i, molecule in enumerate(mols):
            args = (i, molecule, cg_mols[i], box, missing_types, cache, meta[i], default_types)
            futures[i] = e.submit(processing, *args)
    missing_types = dict(missing_types)
    for ele in all_elements:
        missing_types[ele] = dict(missing_types[ele])
    print_missing_type(missing_types, draw=draw)
    for i in futures:
        if futures[i].exception() is not None:  # catch errors first
            print(f"Errors in {i}th molecule: {futures[i].result()}")
        elif futures[i].result() == 0:
            print(f"Molecule {i} is not generated, please re-check force field parameters.")


# end

if __name__ == "__main__":
    xml = XmlParser(argv[1])
    box = (xml.box.lx, xml.box.ly, xml.box.lz, xml.box.xy, xml.box.xz, xml.box.yz)
    box = tuple(map(float, box))

    cg_sys, cg_mols = read_cg_topology(xml, molecules)
    reactor = Reactor(molecules, reaction_template)

    # find tri / di reactions for PFR
    # should be read from xml.data['reaction'] directly
    reactions = []
    for r in xml.data['di']:
        reactions.append(('di', r[0], r[1]))
    for r in xml.data['tri']:
        reactions.append(('tri', r[0], r[1], r[2]))
    # end

    aa_mols, meta = reactor.process(cg_mols, reactions)
    [Chem.SanitizeMol(_) for _ in aa_mols]
    aa_mols_h = [Chem.AddHs(m) for m in aa_mols]
    print(f"{len(aa_mols_h)} molecules!")
    main(aa_mols_h, box, cg_mols, meta, default_types=defaults, draw=True)
