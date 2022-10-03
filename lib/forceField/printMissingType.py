from rdkit import Chem

from lib.draw_molecule_ascii import print_mol_ascii


def print_missing_type(chem_envs, draw=False):
    flag = True
    for ele in chem_envs:
        for chem_env in chem_envs[ele]:
            atom_type = chem_envs[ele][chem_env]
            if atom_type is not None:
                flag = False
                print(
                    f"Chemical env for *{ele} -IN- *{chem_env} is NOT FOUND in database\n"
                    f"Manually set as {atom_type}\n"
                )
    for ele in chem_envs:
        for chem_env in chem_envs[ele]:
            atom_type = chem_envs[ele][chem_env]
            if atom_type is None:
                flag = False
                if draw:
                    try:
                        s = print_mol_ascii(Chem.MolFromSmiles(chem_env, sanitize=False))
                        print(s)
                    except:
                        print("Drawing failed!")
                print(
                    f"Chemical env for *{ele} -IN- *{chem_env} is "
                    "NOT FOUND in database and there is no default type."
                )
                print(f"Please set default type with map key ``{chem_env}''")

    if flag:
        print("All types found in database.")
