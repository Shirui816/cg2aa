__doc__ = r"""
An example for pfr chains.
"""

import os.path
from concurrent.futures import ProcessPoolExecutor as Executor  # cpu-bound jobs
from concurrent.futures import ThreadPoolExecutor as Executor  # io-bound
from multiprocessing import Manager
from sys import argv

from rdkit import Chem

from lib.cgTopology import read_cg_topology
from lib.forceField.opls import ff
from lib.forceField.printMissingType import print_missing_type
from lib.io.xmlParser import XmlParser
from lib.io.xmlWriter import write_xml
from lib.reactor import Reactor
from lib.utils import set_molecule_id_for_h

# ctl = ...
# get molecules from control.in

# Read this part from file
molecules = {
    'A': {'smiles': 'Oc1ccccc1', 'pdb': None},
    'B': {'smiles': 'N=NCC=O', 'pdb': None},  # not real molecule, for testing missing types
}

reaction_template = {
    'tri': {
        'cg_reactant_list': [('A', 'B', 'A')],
        'smarts': "Oc1cccc[c:1]1.[C:2]=[O:3].Oc1cccc[c:4]1>>Oc1cccc[c:1]1[C:2][c:4]2ccccc2O.[O:3]",
        'prod_idx': [0]
    },
    'di': {
        'cg_reactant_list': [('A', 'B')],
        'smarts': "Oc1cccc[c:1]1.[C:2]=[O:3]>>Oc1cccc[c:1]1[C:2][O:3]",
        'prod_idx': [0]
    }
}

defaults = {
    'N': {
        # comment out some types to test missing type
        # "N(=NC(C)([H])[H])[H]": "@atom:nnn",
        "N(=N[H])C(C(c)(c)[H])([H])[H]": "@atom:nnn",
        "N(=N[H])C(C(O)(c)[H])([H])[H]": "@atom:nnn",
    },
    'H': {
        # '[H]N=NC': "@atom:nnn",
    },
    'C': {
        "C(N=N[H])(C(c(c)c)(c(c)c)[H])([H])[H]": "@atom:nnn",
        "C(N=N[H])(C(O[H])(c(c)c)[H])([H])[H]": "@atom:nnn",
        "C(N=N[H])(C(c(c)c)(O[H])[H])([H])[H]": "@atom:nnn",
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


def processing(i, mol, box, mt, ch, meta, defaults, radius=7):
    f_name = f"out_{i:06d}.xml"
    if os.path.isfile(f_name):  # skip existing files
        return 2
    ret = ff(mol, chem_envs=mt, chem_envs_cache=ch, large=500, defaults=defaults, radius=radius)
    if ret is not None:
        plm, bonds, angles, dihedrals = ret
        for m in meta.nodes:  # set monomer id
            molecule = meta.nodes[m]
            for idx in molecule['atom_idx'].values():
                atom = plm.GetAtomWithIdx(idx)
                atom.SetIntProp('molecule_id', int(m))
        plm = set_molecule_id_for_h(plm)
        write_xml(plm, box, bonds, angles, dihedrals, '%06d' % i)
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


def main(mols, box, meta, default_types=None, draw=False):
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
            args = (i, molecule, box, missing_types, cache, meta[i], default_types)
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
    for monomer in cg_sys.nodes:
        _type = cg_sys.nodes[monomer]['type']
        if _type in ['B']:
            if len(cg_sys.adj[monomer]) == 2:
                a, b = list(cg_sys.adj[monomer])
                reactions.append(('tri', a, b, monomer))
            elif len(cg_sys.adj[monomer]) == 1:
                a = list(cg_sys.adj[monomer])[0]
                reactions.append(('di', a, monomer))
    # end

    reactor.process(cg_mols, reactions)
    aa_mols = reactor.aa_molecules
    meta = reactor.meta
    [Chem.SanitizeMol(_) for _ in aa_mols]
    aa_mols_h = [Chem.AddHs(m) for m in aa_mols]
    print(f"{len(aa_mols_h)} molecules!")
    main(aa_mols_h, box, meta, default_types=defaults, draw=True)
