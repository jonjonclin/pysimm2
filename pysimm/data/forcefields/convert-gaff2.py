import json

element_names_by_mass = {1: 'H', 4: 'He', 7: 'Li', 9: 'Be', 11: 'B', 12: 'C',
                         14: 'N', 16: 'O', 19: 'F', 20: 'Ne', 23: 'Na',
                         24: 'Mg', 27: 'Al', 28: 'Si', 31: 'P', 32: 'S',
                         35: 'Cl', 39: 'K', 40: 'Ca', 80: 'Br', 127: 'I'}

r_to_sigma = 2./2.**(1/6)  # conversion factor from r to sigma in LAMMPS

def mass_to_element(mass):
    """Convert atomic mass to element symbol using a predefined mapping."""
    rounded_mass = round(mass)
    return element_names_by_mass.get(rounded_mass, 'X')

def parse_gaff_dat(filename):
    data = {
        "title": "",
        "particle_types": {},
        "bond_types": {},
        "angle_types": {},
        "dihedral_types": {},
        "improper_types": {}
    }
    with open(filename, 'r') as f:
        lines = [line.rstrip('\n') for line in f]

    idx = 0
    # 1. Title
    while idx < len(lines) and not lines[idx].strip():
        idx += 1
    data["title"] = lines[idx].strip()
    idx += 1

    # 2. Atom types
    while idx < len(lines):
        line = lines[idx].strip()
        if not line:
            idx += 1
            break
        # Each line: SYMBOL, MASS, POLARIZABILITY
        p_type = line[:2].strip()
        mass = float(line[2:14].strip())  # F10.4
        # pol = float(line[14:24].strip())  # F10.4
        desc = line[24:].strip()  # rest of the line
        elem = mass_to_element(mass)
        data["particle_types"][p_type] = {"name": p_type, "tag": p_type, "elem": elem,"mass": mass, "desc": desc}
        idx += 1

    # 3. Hydrophilic atom types, line can have up to 20 atom types, skipped
    line = lines[idx].strip()
    # print(f"Hydrophilic line: {line}")
    # atoms = [line[i:i+2].strip() for i in range(0, len(line), 4) if line[i:i+2].strip()]
    # data["hydrophilic"].extend(atoms)
    idx += 1

    # 4. Bonds
    while idx < len(lines):
        line = lines[idx].strip()
        if not line:
            idx += 1
            break
        # Extracting the BOND PARAMETERS
        a1 = line[0:2].strip()  # A2
        a2 = line[3:5].strip()  # 1X, A2 (index 2 is the skipped character)
        if a1 > a2:
            a1, a2 = a2, a1
        k_r = float(line[5:15].strip())  # F10.2
        req = float(line[15:25].strip())  # F10.2
        key = f"{a1}-{a2}"
        name = f"{a1},{a2}"
        rname = f"{a2},{a1}"
        data["bond_types"][key] = {"name": name, "tag": name, "rname": rname, "k": k_r, "r0": req}
        idx += 1

    # 5. Angles
    while idx < len(lines):
        line = lines[idx].strip()
        if not line:
            idx += 1
            break
        # Extracting the ANGLE PARAMETERS
        a1 = line[0:2].strip()
        a2 = line[3:5].strip()
        a3 = line[6:8].strip()
        if a1 > a3:
            a1, a3 = a3, a1
        k_t = float(line[9:19].strip())  # F10.2
        teq = float(line[19:29].strip()) # The second float follows immediately
        key = f"{a1}-{a2}-{a3}"
        name = f"{a1},{a2},{a3}"
        rname = f"{a3},{a2},{a1}"
        data["angle_types"][key] = {"name": name, "tag": name, "rname": rname, "k": k_t, "theta0": teq}
        idx += 1

    # 6. Dihedrals
    while idx < len(lines):
        line = lines[idx].strip()
        if not line:
            idx += 1
            break
        # Extracting the DIHEDRAL PARAMETERS
        a1 = line[0:2].strip()
        a2 = line[3:5].strip()
        a3 = line[6:8].strip()
        a4 = line[9:11].strip()

        if a2 > a3:
            a1, a2, a3, a4 = a4, a3, a2, a1

        idivf = int(line[11:15].strip())  # I4

        pk = float(line[15:30].strip())  # F15.2])
        phase = float(line[30:45].strip())  # F15.2
        pn = float(line[45:60].strip())  # F15.2
        key = f"{a1}-{a2}-{a3}-{a4}"
        name = f"{a1},{a2},{a3},{a4}"
        rname = f"{a4},{a3},{a2},{a1}"

        k_d = pk / idivf  # force constant divided by idivf
        # Dihedrals may have multiple terms for the same key
        # m: # of terms, d: phases, n: periodicities, k: force constants
        if key not in data["dihedral_types"]:
            data["dihedral_types"][key] = {"name": name, "tag": name, "rname": rname, "m": 0, "k":[], "d":[], "n":[]}

        data["dihedral_types"][key]["m"] += 1
        data["dihedral_types"][key]["k"].append(k_d)
        data["dihedral_types"][key]["d"].append(phase)
        data["dihedral_types"][key]["n"].append(pn) 
        idx += 1

    # 7. Impropers
    while idx < len(lines):
        line = lines[idx].strip()
        if not line:
            idx += 1
            break
        # Extracting the IMPROPER PARAMETERS
        a1 = line[0:2].strip()
        a2 = line[3:5].strip()
        a3 = line[6:8].strip()
        a4 = line[9:11].strip()

        sorted_atoms = sorted([a1, a2, a4])
        a1, a2, a4 = sorted_atoms[0], sorted_atoms[1], sorted_atoms[2]
        key = f"{a1}-{a2}-{a3}-{a4}"
        # idivf = int(line[12:16].strip()) # not used in impropers

        pk = float(line[15:30].strip())
        phase = float(line[30:45].strip())
        pn = float(line[45:60].strip())
        name = f"{a1},{a2},{a3},{a4}"
        rname = f"{a4},{a3},{a2},{a1}"
        # Impropers do not have multiple terms for the same key
        data["improper_types"][key] = {"name": name, "tag": name, "rname": rname, "k":pk, "d":phase, "n":pn}
        idx += 1

    # 8. H-bonds; no longer used in GAFF2, but included for completeness
    while idx < len(lines):
        line = lines[idx].strip()
        if not line:
            idx += 1
            break
        idx += 1

    # 9. Equivalences; no longer used in GAFF2, but included for completeness
    while idx < len(lines):
        line = lines[idx].strip()
        if not line:
            idx += 1
            break
        idx += 1

    # 10. Nonbonded parameters
    # Find LABEL and KINDNB
    line = lines[idx].strip()
    # label = line[0:4]  
    # kindnb = line[10:12] 
    idx += 1
    # KINDNB .EQ. 'RE'
    while idx < len(lines):
        line = lines[idx].strip()
        if not line:
            idx += 1
            continue
        if line.startswith('END'):
            break
        
        p_type = line[0:2].strip()
        sigma = r_to_sigma * float(line[10:20].strip()) 
        epsilon = float(line[20:28].strip())
        # data['nonbonded'][ltynb] = {"r": r, "edep": edep}
        data["particle_types"][p_type].update({"sigma": sigma, "epsilon": epsilon})
        idx += 1

    return data

def main():
    files = ['gaff211.dat', 'gaff221.dat']
    for fname in files:
        data = parse_gaff_dat(fname)
        forcefield_name = fname.replace('.dat', '')
        data_out = {
            "ff_name": forcefield_name,
            "ff_class": 1,
            "pair_style": 'lj',
            "bond_style": 'harmonic',
            "angle_style": 'harmonic',
            "dihedral_style": 'fourier',
            "improper_style": 'cvff',
            "mix_rule": 'arithmetic',
        }
        for k, v in data.items():
            data_out[k] = list(v.values()) if isinstance(v, dict) else v
        outname = fname.replace('.dat', '.json')
        with open(outname, 'w') as f:
            json.dump(data_out, f, indent=2)
        print(f"Saved {outname}")

if __name__ == "__main__":
    main()