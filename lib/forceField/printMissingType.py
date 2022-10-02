def print_missing_type(chem_envs):
    if bool(chem_envs):
        for ele in chem_envs:
            for chem_env in chem_envs[ele]:
                atom_type = chem_envs[ele][chem_env]
                if atom_type is not None:
                    print(
                        f"Chemical env for *{ele} -IN- *{chem_env} is NOT FOUND in database\n"
                        f"Manually set as {atom_type}\n"
                    )
        for ele in chem_envs:
            for chem_env in chem_envs[ele]:
                atom_type = chem_envs[ele][chem_env]
                if atom_type is None:
                    failed = True
                    print(
                        f"Chemical env for *{ele} -IN- *{chem_env} is "
                        "NOT FOUND in database and there is no default type."
                    )
                    print(f"Please set default type with map key ``{chem_env}''")
    else:
        print("All types found in database.")
