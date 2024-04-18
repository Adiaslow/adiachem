"""
Module for working with protein structures.
"""

import numpy as np

class Atom:
    """Class for representing an atom in a protein structure.

    Attributes:
        serial (int): Serial number of the atom.
        name (str): Name of the atom.
        alt_loc (str): Alternate location indicator.
        res_name (str): Name of the residue to which the atom belongs.
        chain_id (str): Chain identifier.
        res_seq (int): Residue sequence number.
        icode (str): Insertion code.
        coord (numpy.ndarray): Coordinates of the atom.
        occupancy (float): Occupancy of the atom.
        temp_factor (float): Temperature factor of the atom.
        element (str): Element symbol.
        charge (str): Charge of the atom.

    Methods:
        transform: Transform the coordinates of the atom using a given matrix.
    """
    def __init__(self, serial, name, alt_loc, res_name, chain_id, res_seq, icode, coord, occupancy, temp_factor, element, charge):
        """Initialize the Atom object.

        Args:
            serial (int): Serial number of the atom.
            name (str): Name of the atom.
            alt_loc (str): Alternate location indicator.
            res_name (str): Name of the residue to which the atom belongs.
            chain_id (str): Chain identifier.
            res_seq (int): Residue sequence number.
            icode (str): Insertion code.
            coord (numpy.ndarray): Coordinates of the atom.
            occupancy (float): Occupancy of the atom.
            temp_factor (float): Temperature factor of the atom.
            element (str): Element symbol.
            charge (str): Charge of the atom.

        Returns:
            None

        Raises:
            None
        """
        self.serial = serial
        self.name = name
        self.alt_loc = alt_loc
        self.res_name = res_name
        self.chain_id = chain_id
        self.res_seq = res_seq
        self.icode = icode
        self.coord = coord
        self.occupancy = occupancy
        self.temp_factor = temp_factor
        self.element = element
        self.charge = charge

    def transform(self, matrix):
        """Transform the coordinates of the atom using a given matrix.

        Args:
            matrix (numpy.ndarray): Transformation matrix.

        Returns:
            None

        Raises:
            None
        """
        self.coord = np.dot(matrix, np.append(self.coord, 1))[:3]

class Residue:
    """Class for representing a residue in a protein structure.

    Attributes:
        res_name (str): Name of the residue.
        chain_id (str): Chain identifier.
        res_seq (int): Residue sequence number.
        icode (str): Insertion code.
        atoms (list): List of atoms in the residue.

    Methods:
        add_atom: Add an atom to the residue.
        transform: Transform the coordinates of the residue using a given matrix.
    """
    def __init__(self, res_name, chain_id, res_seq, icode):
        """Initialize the Residue object.

        Args:
            res_name (str): Name of the residue.
            chain_id (str): Chain identifier.
            res_seq (int): Residue sequence number.
            icode (str): Insertion code.

        Returns:
            None

        Raises:
            None
        """
        self.res_name = res_name
        self.chain_id = chain_id
        self.res_seq = res_seq
        self.icode = icode
        self.atoms = []

    def add_atom(self, atom):
        """Add an atom to the residue.

        Args:
            atom (Atom): Atom object to add to the residue.

        Returns:
            None

        Raises:
            None
        """
        self.atoms.append(atom)

    def transform(self, matrix):
        """Transform the coordinates of the residue using a given matrix.

        Args:
            matrix (numpy.ndarray): Transformation matrix.

        Returns:
            None

        Raises:
            None
        """
        for atom in self.atoms:
            atom.transform(matrix)

class Chain:
    """Class for representing a chain in a protein structure.

    Attributes:
        chain_id (str): Chain identifier.
        residues (list): List of residues in the chain.

    Methods:
        add_residue: Add a residue to the chain.
        transform: Transform the coordinates of the chain using a given matrix.
    """
    def __init__(self, chain_id):
        """Initialize the Chain object.

        Args:
            chain_id (str): Chain identifier.

        Returns:
            None

        Raises:
            None
        """
        self.chain_id = chain_id
        self.residues = []

    def add_residue(self, residue):
        """Add a residue to the chain.

        Args:
            residue (Residue): Residue object to add to the chain.

        Returns:
            None

        Raises:
            None
        """
        self.residues.append(residue)

    def transform(self, matrix):
        """Transform the coordinates of the chain using a given matrix.

        Args:
            matrix (numpy.ndarray): Transformation matrix.

        Returns:
            None

        Raises:
            None
        """
        for residue in self.residues:
            residue.transform(matrix)

class Protein:
    """Class for representing a protein structure.

    Attributes:
        chains (list): List of chains in the protein.

    Methods:
        load_pdb: Load a PDB file into the protein object.
        _parse_atom_record: Parse an ATOM record in a PDB file.
        _find_or_create_chain: Find or create a chain in the protein.
        _find_or_create_residue: Find or create a residue in a chain.
        transform_all_chains: Transform the coordinates of all chains in the protein using a given matrix.
        transform_chain: Transform the coordinates of a chain in the protein using a given matrix.
        transform_residue: Transform the coordinates of a residue in the protein using a given matrix.
        transform_atom: Transform the coordinates of an atom in the protein using a given matrix.
        save_pdb: Save the protein structure to a PDB file.
        _format_atom_record: Format an ATOM record for a given atom.
    """
    def __init__(self, pdb_file):
        """Initialize the Protein object.

        Args:
            pdb_file (str): Path to the PDB file to load.

        Returns:
            None

        Raises:
            None
        """
        self.chains = []
        self.load_pdb(pdb_file)

    def load_pdb(self, pdb_file):
        """Load a PDB file into the protein object.

        Args:
            pdb_file (str): Path to the PDB file to load.

        Returns:
            None

        Raises:
            None
        """
        with open(pdb_file, 'r') as file:
            for line in file:
                if line.startswith('ATOM'):
                    self._parse_atom_record(line)

    def _parse_atom_record(self, line):
        """Parse an ATOM record in a PDB file.

        Args:
            line (str): ATOM record line from a PDB file.

        Returns:
            None

        Raises:
            None
        """
        serial = int(line[6:11].strip())
        name = line[12:16].strip()
        alt_loc = line[16].strip()
        res_name = line[17:20].strip()
        chain_id = line[21].strip()
        res_seq = int(line[22:26].strip())
        icode = line[26].strip()
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        coord = np.array([x, y, z])
        occupancy = float(line[54:60].strip())
        temp_factor = float(line[60:66].strip())
        element = line[76:78].strip()
        charge = line[78:80].strip()

        atom = Atom(serial, name, alt_loc, res_name, chain_id, res_seq, icode, coord, occupancy, temp_factor, element, charge)

        chain = self._find_or_create_chain(chain_id)
        residue = self._find_or_create_residue(res_name, chain_id, res_seq, icode)
        residue.add_atom(atom)

    def _find_or_create_chain(self, chain_id):
        """Find or create a chain in the protein.

        Args:
            chain_id (str): Chain identifier.

        Returns:
            Chain: Chain object.

        Raises:
            None
        """
        for chain in self.chains:
            if chain.chain_id == chain_id:
                return chain
        chain = Chain(chain_id)
        self.chains.append(chain)
        return chain

    def _find_or_create_residue(self, res_name, chain_id, res_seq, icode):
        """Find or create a residue in a chain.

        Args:
            res_name (str): Residue name.
            chain_id (str): Chain identifier.
            res_seq (int): Residue sequence number.
            icode (str): Insertion code.

        Returns:
            Residue: Residue object.

        Raises:
            None
        """
        chain = self._find_or_create_chain(chain_id)
        for residue in chain.residues:
            if residue.res_seq == res_seq and residue.icode == icode:
                return residue
        residue = Residue(res_name, chain_id, res_seq, icode)
        chain.add_residue(residue)
        return residue

    def transform_all_chains(self, matrix):
        """Transform the coordinates of all chains in the protein using a given matrix.

        Args:
            matrix (numpy.ndarray): Transformation matrix.

        Returns:
            None

        Raises:
            None
        """
        for chain in self.chains:
            chain.transform(matrix)

    def transform_chain(self, chain_id, matrix):
        """Transform the coordinates of a chain in the protein using a given matrix.

        Args:
            chain_id (str): Chain identifier.
            matrix (numpy.ndarray): Transformation matrix.

        Returns:
            None

        Raises:
            None
        """
        chain = self._find_or_create_chain(chain_id)
        chain.transform(matrix)

    def transform_residue(self, chain_id, res_seq, icode, matrix):
        """Transform the coordinates of a residue in the protein using a given matrix.

        Args:
            chain_id (str): Chain identifier.
            res_seq (int): Residue sequence number.
            icode (str): Insertion code.
            matrix (numpy.ndarray): Transformation matrix.

        Returns:
            None

        Raises:
            None
        """"
        chain = self._find_or_create_chain(chain_id)
        for residue in chain.residues:
            if residue.res_seq == res_seq and residue.icode == icode:
                residue.transform(matrix)
                break

    def transform_atom(self, chain_id, res_seq, icode, atom_name, matrix):
        """Transform the coordinates of an atom in the protein using a given matrix.

        Args:
            chain_id (str): Chain identifier.
            res_seq (int): Residue sequence number.
            icode (str): Insertion code.
            atom_name (str): Atom name.
            matrix (numpy.ndarray): Transformation matrix.

        Returns:
            None

        Raises:
            None
        """
        chain = self._find_or_create_chain(chain_id)
        for residue in chain.residues:
            if residue.res_seq == res_seq and residue.icode == icode:
                for atom in residue.atoms:
                    if atom.name == atom_name:
                        atom.transform(matrix)
                        break
                break

    def save_pdb(self, output_file):
        """Save the protein structure to a PDB file.

        Args:
            output_file (str): Output file path.

        Returns:
            None

        Raises:
            None
        """
        with open(output_file, 'w') as file:
            for chain in self.chains:
                for residue in chain.residues:
                    for atom in residue.atoms:
                        file.write(self._format_atom_record(atom))
                        file.write('\n')

    def _format_atom_record(self, atom):
        """Format an ATOM record line for a given atom.

        Args:
            atom (Atom): Atom object.

        Returns:
            str: Formatted ATOM record line.

        Raises:
            None
        """
        record = 'ATOM  {:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n'.format(
            atom.serial, atom.name, atom.alt_loc, atom.res_name, atom.chain_id, atom.res_seq, atom.icode,
            atom.coord[0], atom.coord[1], atom.coord[2], atom.occupancy, atom.temp_factor, atom.element, atom.charge
        )
        return record
