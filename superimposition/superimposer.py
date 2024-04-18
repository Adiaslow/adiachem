"""
Module for superimposing protein structures.
"""

class ProteinSuperimposer:
    """ Class for superimposing protein structures.

        Parameters:
            reference_structure (Structure): Reference protein structure.
            test_structure (Structure): Test protein structure.

        Attributes:
            reference_structure (Structure): Reference protein structure.
            test_structure (Structure): Test protein structure.
            sequence_aligner (SequenceAligner): Sequence aligner object.

        Methods:
            superimpose(): Superimpose the test structure onto the reference structure.
            _extract_sequences(): Extract sequences from the reference and test models.
            _extract_aligned_coordinates(): Extract coordinates of aligned residues from the reference and test models.
            _calculate_transformation(): Calculate the optimal rotation and translation for superimposition.
            _apply_transformation(): Apply the transformation to the test structure.
            _convert_to_one_letter(): Convert sequence to one-letter code.
            _traceback_path(): Traceback the path of the optimal alignment.
            visualize_alignment(): Visualize the alignment of two sequences.
            visualize_alignment_matrix(): Visualize the alignment matrix.
    """
    def __init__(self, reference_structure, test_structure, ref_superimpose_chains, test_superimpose_chains, ref_rmsd_chains, test_rmsd_chains):
        """ Initialize the ProteinSuperimposer object.

            Args:
                reference_structure (Structure): Reference protein structure.
                test_structure (Structure): Test protein structure.

            Returns:
                None.

            Raises:
                None.
        """
        self.reference_structure = reference_structure
        self.test_structure = test_structure
        self.ref_superimpose_chains = ref_superimpose_chains
        self.test_superimpose_chains = test_superimpose_chains
        self.ref_rmsd_chains = ref_rmsd_chains
        self.test_rmsd_chains = test_rmsd_chains
        self.sequence_aligner = SequenceAligner()

    def superimpose(self):
        """ Superimpose the test structure onto the reference structure.

            Args:
                None.

            Returns:
                superimposed_structure (Structure): Superimposed test structure.

            Raises:
                None.
        """
        best_rmsd = np.inf
        best_rotation_angle = None

        for i in range(4):
            for ref_model, test_model in zip(self.reference_structure, self.test_structure):
                ref_seq, test_seq = self._extract_sequences(ref_model, test_model, self.ref_superimpose_chains, self.test_superimpose_chains)
                ref_aligned_seq, test_aligned_seq = self.sequence_aligner.align(ref_seq, test_seq)

                ref_coords, test_coords = self._extract_aligned_coordinates(ref_model, test_model, ref_aligned_seq, test_aligned_seq, self.ref_superimpose_chains, self.test_superimpose_chains)
                ref_main_axis = self._find_axis_of_symmetry_by_convex_hull(ref_coords, num_sections=20, plot=False)
                test_main_axis = self._find_axis_of_symmetry_by_convex_hull(test_coords, num_sections=20, plot=False)

                # self._align_with_axis(ref_model, ref_main_axis)
                self._align_with_axis(test_model, ref_main_axis, test_main_axis)

                ref_coords_aligned, test_coords_aligned = self._extract_aligned_coordinates(ref_model, test_model, ref_aligned_seq, test_aligned_seq, self.ref_superimpose_chains, self.test_superimpose_chains)
                rotation_matrix, translation_vector = self._calculate_transformation(ref_coords_aligned, test_coords_aligned)

                identity_matrix = np.eye(3)

                self._apply_transformation(test_model, identity_matrix, translation_vector)

                angle = (np.pi / 2) * i
                axis_of_rotation = test_main_axis
                self._apply_rotation(test_model, axis_of_rotation, angle)

                # Fine-tune the rotation about the main axis
                for angle_deg in np.arange(-10, 11, 0.1):
                    angle_rad = np.deg2rad(angle_deg)
                    self._apply_rotation(test_model, axis_of_rotation, angle_deg)
                    _, test_coords_rotated = self._extract_aligned_coordinates(ref_model, test_model, ref_aligned_seq, test_aligned_seq, self.ref_superimpose_chains, self.test_superimpose_chains)
                    rotation_rmsd = self._rmsd(ref_coords_aligned, test_coords_rotated)
                    if rotation_rmsd < best_rmsd:
                        best_rmsd = rotation_rmsd
                        best_rotation_angle = angle_rad

                self._apply_rotation(test_model, axis_of_rotation, best_rotation_angle)

                ref_coords_final, test_coords_final = self._extract_aligned_coordinates(ref_model, test_model, ref_aligned_seq, test_aligned_seq, self.ref_rmsd_chains, self.test_rmsd_chains)
                final_rmsd = self._rmsd(ref_coords_final, test_coords_final)
                print(f"Final RMSD (rotation {i * 90} degrees): {final_rmsd}")

        print(f"Best RMSD: {best_rmsd}")
        # self._plot_superimposition(ref_model, test_model)
        return self.test_structure

    def _extract_chain_coordinates(self, model, chains):
        """ Get the coordinates of the specified chains in the model.

            Args:
                model (Structure): Protein model.
                chains (list): List of chain identifiers.

            Returns:
                coords (np.array): Coordinates of the specified chains.

            Raises:
                None.
        """
        coords = []
        for chain_id in chains:
            chain = model[chain_id]
            for residue in chain:
                for atom in residue:
                    coords.append(atom.get_coord())
        return np.array(coords)

    def _get_all_atom_coordinates(self, model):
        """ Get the coordinates of all atoms in the model.

            Args:
                model (Model): Protein model.

            Returns:
                all_coords (np.array): Coordinates of all atoms in the model.

            Raises:
                None.
        """
        all_coords = []
        for chain in model:
            for residue in chain:
                for atom in residue.get_atoms():
                    all_coords.append(atom.get_coord())
        return np.array(all_coords)

    def _extract_sequences(self, ref_model, test_model, ref_chains, test_chains):
        """ Extract sequences of the specified chains from the reference and test models.

            Args:
                ref_model (Structure): Reference protein model.
                test_model (Structure): Test protein model.
                ref_chains (list): List of chain identifiers for the reference model.
                test_chains (list): List of chain identifiers for the test model.

            Returns:
                ref_seq (str): Reference sequence.
                test_seq (str): Test sequence.

            Raises:
                None.
        """
        ref_seq = ""
        test_seq = ""

        for ref_chain_id, test_chain_id in zip(ref_chains, test_chains):
            ref_chain = ref_model[str(ref_chain_id)]  # Convert to string
            test_chain = test_model[str(test_chain_id)]  # Convert to string

            for ref_res, test_res in zip(ref_chain, test_chain):
                ref_seq += ref_res.get_resname()
                test_seq += test_res.get_resname()

        if len(ref_seq) > len(test_seq):
            test_seq += "-" * (len(ref_seq) - len(test_seq))
        elif len(ref_seq) < len(test_seq):
            ref_seq += "-" * (len(test_seq) - len(ref_seq))
        return ref_seq, test_seq

    def _extract_aligned_coordinates(self, ref_model, test_model, ref_aligned_seq, test_aligned_seq, ref_chains, test_chains):
        """ Extract aligned coordinates from the reference and test models.

            Args:
                ref_model (Structure): Reference protein model.
                test_model (Structure): Test protein model.
                ref_aligned_seq (str): Reference aligned sequence.
                test_aligned_seq (str): Test aligned sequence.
                ref_chains (list): List of chain identifiers for the reference model.
                test_chains (list): List of chain identifiers for the test model.

            Returns:
                ref_coords (np.array): Reference coordinates.
                test_coords (np.array): Test coordinates.

            Raises:
                None.
        """
        ref_coords = []
        test_coords = []
        ref_residues = list(ref_model.get_residues())
        test_residues = list(test_model.get_residues())
        ref_index = 0
        test_index = 0
        for ref_chain_id, test_chain_id in zip(ref_chains, test_chains):
            ref_chain = ref_model[str(ref_chain_id)]
            test_chain = test_model[str(test_chain_id)]
            for ref_res, test_res, ref_aa, test_aa in zip(ref_chain, test_chain, ref_aligned_seq, test_aligned_seq):
                if ref_aa != "-" and test_aa != "-":
                    if ref_index < len(ref_residues) and test_index < len(test_residues):
                        ref_res = ref_residues[ref_index]
                        test_res = test_residues[test_index]
                        if "CA" in ref_res and "CA" in test_res:
                            ref_coords.append(ref_res["CA"].get_coord())
                            test_coords.append(test_res["CA"].get_coord())
                    else:
                        break
                if ref_aa != "-":
                    ref_index += 1
                if test_aa != "-":
                    test_index += 1

        return np.array(ref_coords), np.array(test_coords)

    def _calculate_transformation(self, ref_coords, test_coords):
        """ Calculate the optimal rotation and translation for superimposition.

            Args:
                ref_coords (np.array): Reference coordinates.
                test_coords (np.array): Test coordinates.

            Returns:
                rotation_matrix (np.array): Optimal rotation matrix.
                translation_vector (np.array): Translation vector.

            Raises:
                None.
        """
        ref_centroid = self._calculate_centroid(ref_coords)
        test_centroid = self._calculate_centroid(test_coords)
        centered_ref_coords = ref_coords - ref_centroid
        centered_test_coords = test_coords - test_centroid
        covariance_matrix = self._calculate_covariance_matrix(centered_ref_coords, centered_test_coords)
        rotation_matrix = self._calculate_optimal_rotation(covariance_matrix)
        translation_vector = self._calculate_translation_vector(ref_centroid, test_centroid, rotation_matrix)
        return rotation_matrix, translation_vector

    def _find_axis_of_symmetry(self, model, axis=0, plot=False):
        """ Find the main axis of symmetry using PCA.

            Args:
                model (Structure): Protein model.

            Returns:
                main_axis (np.array): Main axis of symmetry.

            Raises:
                None.
        """
        # Extract coordinates of all atoms in the model
        coords = []
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coords.append(atom.get_coord())

        coords = np.array(coords)

        # Perform PCA
        pca = PCA(n_components=3)
        pca.fit(coords)
        main_axis = pca.components_[axis]
        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            arrow_scale = 20
            # Get the first principal axis (main axis of symmetry)

            # plt.scatter(coords[:, 0], coords[:, 1], coords[:, 2])
            plt.plot(coords[:, 1], coords[:, 2])
            plt.quiver(0, 0, 0, pca.components_[0][0], pca.components_[0][1], pca.components_[0][2], color='red', zorder=2)
            plt.quiver(0, 0, 0, pca.components_[1][0], pca.components_[1][1], pca.components_[1][2], color='orange', zorder=2)
            plt.quiver(0, 0, 0, pca.components_[2][0], pca.components_[2][1], pca.components_[2][2], color='olive', zorder=2)
            plt.show()

        return main_axis

    def _find_axis_of_symmetry_by_convex_hull(self, coords, num_sections=30, plot=False):
        """ Find the axis that runs through the center of the pore.

            Args:
                coords (np.array): Aligned coordinates.
                num_sections (int): Number of sections to divide the protein into.
                plot (bool): Whether to plot the results.

            Returns:
                pore_axis (np.array): Axis that runs through the center of the pore.

            Raises:
                None.
        """
        # Perform PCA on the atom coordinates to identify the longest axis
        pca = PCA(n_components=3)
        pca.fit(coords)
        longest_axis = pca.components_[1]

        # Project the coordinates onto the longest axis
        projected_coords = np.dot(coords, longest_axis)

        # Divide the protein into sections along the longest axis
        min_proj, max_proj = np.min(projected_coords), np.max(projected_coords)
        section_boundaries = np.linspace(min_proj, max_proj, num_sections + 1)

        # Calculate the centroid of each section
        section_centroids = []
        for i in range(num_sections):
            section_mask = (projected_coords >= section_boundaries[i]) & (projected_coords < section_boundaries[i + 1])
            section_coords = coords[section_mask]
            section_centroid = np.mean(section_coords, axis=0)
            section_centroids.append(section_centroid)

        section_centroids = np.array(section_centroids)

        # Perform PCA on the section centroids
        pca_sections = PCA(n_components=3)
        pca_sections.fit(section_centroids)
        pore_axis = pca_sections.components_[1]

        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot(coords[:, 0], coords[:, 1], coords[:, 2])
            # ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c='b', marker='o', alpha=0.1)
            ax.scatter(section_centroids[:, 0], section_centroids[:, 1], section_centroids[:, 2], c='r', marker='o')
            ax.quiver(np.mean(section_centroids[:, 0]), np.mean(section_centroids[:, 1]), np.mean(section_centroids[:, 2]),
                    pore_axis[0], pore_axis[1], pore_axis[2], color='r', length=20)
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            plt.show()

        return pore_axis

    def _align_with_axis(self, model, ref_axis, test_axis, plot=False):
        """ Align the model with the main axis of symmetry.

            Args:
                model (Structure): Protein model.
                ref_axis (np.array): Reference axis of symmetry.
                test_axis (np.array): Test axis of symmetry.
                plot (bool): Whether to plot the results.

            Returns:
                None.

            Raises:
                None.
        """
        angle_between_axes = np.arccos(np.dot(ref_axis, test_axis))
        angle_between_axes = np.degrees(angle_between_axes)
        self._apply_rotation(model, ref_axis, angle_between_axes)

    def _apply_transformation(self, model, rotation_matrix, translation_vector):
        """ Apply the rotation and translation to the model.

            Args:
                model (Structure): Protein model to be transformed.
                rotation_matrix (np.array): Rotation matrix.
                translation_vector (np.array): Translation vector.

            Returns:
                None.

            Raises:
                None.
        """
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coord = atom.get_coord()
                    transformed_coord = np.dot(coord, rotation_matrix) + translation_vector
                    atom.set_coord(transformed_coord)

    def _apply_rotation(self, model, rotation_axis, angle):
        """ Apply the rotation to the model.

            Args:
                model (Structure): Protein model to be transformed.
                rotation_axis (np.array): Axis of rotation.
                angle (float): Angle of rotation in radians.

            Returns:
                None.

            Raises:
                None.
        """
        # Normalize the rotation axis
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)

        # Create the rotation matrix using the axis and angle
        cos_theta = np.cos(angle)
        sin_theta = np.sin(angle)
        ux, uy, uz = rotation_axis
        rotation_matrix = np.array([
            [cos_theta + ux**2 * (1 - cos_theta), ux * uy * (1 - cos_theta) - uz * sin_theta, ux * uz * (1 - cos_theta) + uy * sin_theta],
            [uy * ux * (1 - cos_theta) + uz * sin_theta, cos_theta + uy**2 * (1 - cos_theta), uy * uz * (1 - cos_theta) - ux * sin_theta],
            [uz * ux * (1 - cos_theta) - uy * sin_theta, uz * uy * (1 - cos_theta) + ux * sin_theta, cos_theta + uz**2 * (1 - cos_theta)]
        ])

        # Apply the rotation to each atom in the model
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coord = atom.get_coord()
                    transformed_coord = np.dot(rotation_matrix, coord)
                    atom.set_coord(transformed_coord)

    def _apply_translation(self, model, translation_vector):
        """ Apply the translation to the model.

            Args:
                model (Structure): Protein model to be transformed.
                translation_vector (np.array): Translation vector.

            Returns:
                None.

            Raises:
                None.
        """
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coord = atom.get_coord()
                    transformed_coord = coord - translation_vector
                    atom.set_coord(transformed_coord)

    def _plot_superimposition(self, ref_model, superimposed_model):
        """ Plot the superimposition of the reference and superimposed models.

            Args:
                ref_model (Model): Reference model.
                superimposed_model (Model): Superimposed model.

            Returns:
                None.

            Raises:
                None.
        """
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Plot the reference model
        ref_coords = self._get_all_atom_coordinates(ref_model)
        ax.plot(ref_coords[:, 0], ref_coords[:, 1], ref_coords[:, 2], label='Reference', alpha=1)

        # Plot the superimposed model
        superimposed_coords = self._get_all_atom_coordinates(superimposed_model)
        ax.plot(superimposed_coords[:, 0], superimposed_coords[:, 1], superimposed_coords[:, 2], label='Superimposed', alpha=1)

        # Set labels and title
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Protein Superimposition')
        ax.legend()

        # Adjust the view angle
        ax.view_init(elev=20, azim=90)

        plt.tight_layout()
        plt.show()
