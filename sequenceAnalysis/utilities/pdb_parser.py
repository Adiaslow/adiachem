import pandas as pd


def parse_pdb_to_dataframe(file_path):
    data = []
    with open(file_path, 'r') as pdb_file:
        pdb_records = pdb_file.read().split("\n")
        for record in pdb_records:
            if record.startswith("ATOM"):
                atom_info = record.split()
                atom_serial = int(atom_info[1])
                atom_name = atom_info[2]
                residue_name = atom_info[3]
                chain_id = atom_info[4]
                if len(chain_id) != 1:
                    residue_number = int(chain_id[1:])
                    chain_id = [0]
                    x_coord = float(atom_info[5])
                    y_coord = float(atom_info[6])
                    z_coord = float(atom_info[7])
                    occupancy = float(atom_info[8])
                    temp_factor = float(atom_info[9])
                    atom_type = atom_info[-1]
                else:
                    residue_number = int(atom_info[5])
                    chain_id = atom_info[4]
                    x_coord = float(atom_info[6])
                    y_coord = float(atom_info[7])
                    z_coord = float(atom_info[8])
                    occupancy = float(atom_info[9])
                    temp_factor = float(atom_info[10])
                    atom_type = atom_info[-1]
                data.append([atom_serial, atom_name, residue_name,
                             chain_id, residue_number, x_coord, y_coord, z_coord, occupancy, temp_factor, atom_type])
    columns = ['AtomSerial', 'AtomName', 'ResidueName',
               'ChainID', 'ResidueNumber', 'X', 'Y', 'Z','Occupancy', 'TempFactor', 'AtomType']
    df = pd.DataFrame(data, columns=columns)
    return df

def dataframe_to_pdb(df, output_file):
    with open(output_file, 'w') as pdb_file:
        for _, row in df.iterrows():
            atom_line = (
                f"ATOM  {row['AtomSerial']:>5} "
                f"{row['AtomName']:>4}{row['ResidueName']:>3} "
                f"{row['ChainID']}{row['ResidueNumber']:>4}    "
                f"{row['X']:>8.3f}{row['Y']:>8.3f}{row['Z']:>8.3f}"
                f"{row['Occupancy']:>6.2f}{row['SASA']:>6.2f}"
                f"          {row['AtomType']:>2}"
            )
            if 'SASA' in row:
                atom_line += f"{row['SASA']:>6.2f}"
            pdb_file.write(atom_line + "\n")
    pdb_file.close()
    print(f"PDB file '{output_file}' successfully created.")