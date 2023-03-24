import networkx as nx


def read_cg_topology(cg_system, monomers):
    cg_sys = nx.Graph()
    for bond in cg_system.data['bond']:
        bond_type, i, j = str(bond[0]), int(bond[1]), int(bond[2])
        type_i, type_j = cg_system.data['type'][i], cg_system.data['type'][j]
        cg_sys.add_node(i, type=type_i, reaction_type=monomers[type_i].get('reaction_type'),
                        smiles=monomers[type_i]['smiles'], x=cg_system.data['position'][i])
        cg_sys.add_node(j, type=type_j, reaction_type=monomers[type_j].get('reaction_type'),
                        smiles=monomers[type_j]['smiles'], x=cg_system.data['position'][j])
        cg_sys.add_edge(i, j, bond_type=bond_type)
    cg_molecules = [cg_sys.subgraph(c).copy() for c in nx.connected_components(cg_sys)]
    return cg_sys, cg_molecules
