__doc__ = r"""
An example for pfr chains.
"""

from sys import argv

from rdkit import Chem

from lib.cgTopology import read_cg_topology
from lib.forceField.opls import ff
from lib.io.writer import write_xml
from lib.io.xmlParser import XmlParser
from lib.reactor import Reactor
from lib.utils import set_molecule_id_for_h

xml = XmlParser(argv[1])
box = (xml.box.lx, xml.box.ly, xml.box.lz, xml.box.xy, xml.box.xz, xml.box.yz)
box = tuple(map(float, box))

# ctl = ...
# get molecules from control.in


molecules = {
    'A': {'smiles': 'Oc1ccccc1', 'pdb': None},
    'B': {'smiles': 'C=O', 'pdb': None}
}

for key in molecules:
    if not molecules[key]['pdb'] is None:
        mol = Chem.RemoveAllHs(Chem.MolFromPDBFile(molecules[key]['pdb']))
        molecules[key]['smiles'] = Chem.MolToSmiles(mol)

# for key in molecules:
#     mol = Chem.MolFromSmiles(molecules[key]['smiles'])
#     mol = Chem.AddHs(mol)
#     AllChem.EmbedMolecule(mol)
#     pdb = Chem.MolToPDBFile(mol, '%s.pdb' % key)


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

cg_sys, cg_mols = read_cg_topology(xml, molecules)
reactor = Reactor(molecules, reaction_template)

# find tri / di reactions for PFR
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

reactor.process(cg_mols, reactions)
aa_mols = reactor.aa_molecules
[Chem.SanitizeMol(_) for _ in aa_mols]
aa_mols_h = [Chem.AddHs(m) for m in aa_mols]
aa_mols_h = [set_molecule_id_for_h(mh) for mh in aa_mols_h]

defaults = None


def processing(i):
    molecule = aa_mols_h[i]
    plm_h_ff, bonds, angles, dihedrals = ff(molecule, defaults=defaults, radius=None)
    write_xml(plm_h_ff, box, bonds, angles, dihedrals, '%06d' % i)
    # mol = Chem.AddHs(molecule)
    # AllChem.EmbedMolecule(mol)
    # pdb = Chem.MolToPDBFile(mol, '%06d.pdb' % i)


if __name__ == "__main__":
    for i in range(len(aa_mols_h)):
        processing(i)
