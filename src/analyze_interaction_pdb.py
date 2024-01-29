from Bio.PDB import PDBParser, DSSP
import numpy as np

class AnalyzeInteractionPdb:

    def __init__(self, file_name):
        self.file_name = file_name

    # 読み込んだファイルから以下の相互作用を解析する
    # ジスルフィド結合、疎水性相互作用、イオン結合、水素結合
    def analyze_interactions(self):

        # ファイルの読み込み
        pdb_parser = PDBParser()
        pdb_path = "pdb/" + self.file_name + ".pdb"
        structure = pdb_parser.get_structure("PDB_ID", pdb_path)

        # ファイルの情報読み取り
        pdb_molecule = structure.header['compound']['1']['molecule']
        pdb_ec = structure.header['source']['1']['organism_scientific']
        pdb_organism_scientific = structure.header['source']['1']['organism_scientific']
        pdb_organism_taxid = structure.header['source']['1']['organism_taxid']

        print(f"{pdb_molecule} ({pdb_ec}) in {pdb_organism_scientific} ({pdb_organism_taxid})")

        try:
            # ジスルフィド結合の数
            pdb_disulfile_bonds = self.find_disulfide_bonds(structure)
            print(f"disulfide : {len(pdb_disulfile_bonds)}")
        except Exception as e:
            # This block will handle any other exceptions
            print(f"An error occurred: {e}")
        try:
            # 疎水性相互作用の数
            pdb_internal_hydrophobic_interactions = self.find_internal_hydrophobic_interactions(structure)
            print(f"hydrophobic : {len(pdb_internal_hydrophobic_interactions)}")
        except Exception as e:
            # This block will handle any other exceptions
            print(f"An error occurred: {e}")
        try:
            # イオン結合
            pdb_ionic_interactions = self.find_ionic_interactions(structure)
            print(f"ionic : {len(pdb_ionic_interactions)}")
        except Exception as e:
            # This block will handle any other exceptions
            print(f"An error occurred: {e}")
        try:
            # 水素結合
            pdb_hydrogen_bonds = self.find_hydrogen_bonds(structure)
            print(f"hydrogen : {len(pdb_hydrogen_bonds)}")
        except Exception as e:
            # This block will handle any other exceptions
            print(f"An error occurred: {e}")

    
    # ジスルフィド結合
    def find_disulfide_bonds(self, structure):
        """異なる鎖にあるシステイン残基間の距離を計算し、条件に合うものを見つける"""
        cys_atoms = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() == "CYS":
                        # システインの硫黄原子を探す
                        for atom in residue.get_atoms():
                            if atom.get_name() == "SG":
                                cys_atoms.append((chain, atom))
                                break
                            
        distances = []
        for i, (chain1, atom1) in enumerate(cys_atoms):
            for (chain2, atom2) in cys_atoms[i+1:]:
                # 異なる鎖のシステインを考慮
                distance = self.calc_distance(atom1, atom2)
                if 2.05 <= distance <= 2.1:
                    distances.append((chain1, atom1.get_id(), chain2, atom2.get_id(), distance))
        return distances
    
    # 疎水性相互作用
    def find_internal_hydrophobic_interactions(self, structure, threshold=4.0, accessible_threshold=25):
        hydrophobic_residues = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'GLY'}
        model = next(structure.get_models())  # 最初のモデルを使用
        dssp = DSSP(model, "pdb/2ejv.pdb", dssp='mkdssp')
        print(dssp)

        internal_atoms = []
        for chain in model:
            print("chain")
            for residue in chain:
                if residue.get_resname() in hydrophobic_residues:
                    # DSSPが認識できる形式のresidue IDを取得
                    dssp_key = (chain.id, residue.id)
                    if dssp_key in dssp and dssp[dssp_key][3] < accessible_threshold:  # ASAをチェック
                        for atom in residue.get_atoms():
                            if atom.element != 'H':  # 水素以外の原子を考慮
                                internal_atoms.append(atom)
                                print(atom)

        interactions = []
        for i, atom1 in enumerate(internal_atoms):
            for atom2 in internal_atoms[i+1:]:
                distance = self.calc_distance(atom1, atom2)
                if distance <= threshold:
                    interactions.append((atom1, atom2))
                    print(f"{atom1}, {atom2}")

        return interactions
    
    # イオン結合
    def find_ionic_interactions(self, structure, threshold=4.0):
        """イオン結合を持ちうるアミノ酸残基の組み合わせを探す"""
        positive_residues = {'ARG', 'LYS', 'HIS'}
        negative_residues = {'ASP', 'GLU'}

        positive_atoms = []
        negative_atoms = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    res_name = residue.get_resname()
                    if res_name in positive_residues or res_name in negative_residues:
                        for atom in residue.get_atoms():
                            if atom.element in ['N', 'O']:  # 窒素と酸素原子を対象とする
                                if res_name in positive_residues:
                                    positive_atoms.append(atom)
                                else:
                                    negative_atoms.append(atom)

        interactions = []
        for atom_pos in positive_atoms:
            for atom_neg in negative_atoms:
                distance = self.calc_distance(atom_pos, atom_neg)
                if distance <= threshold:
                    interactions.append((atom_pos, atom_neg))

        return interactions
    
    # 水素結合
    def find_hydrogen_bonds(self, structure, distance_threshold=3.5):
        hydrogen_bonds = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if atom.element == 'N':
                            for residue2 in chain:
                                if residue != residue2:
                                    for atom2 in residue2:
                                        if atom2.element == 'O':
                                            distance = self.calc_distance(atom, atom2)
                                            if distance <= distance_threshold:
                                                hydrogen_bonds.append((atom, atom2))
        return hydrogen_bonds
    
    # 原子間距離計算
    def calc_distance(self, atom1, atom2):
        """二つの原子間の距離を計算する"""
        diff_vector = atom1.coord - atom2.coord
        return np.sqrt(np.sum(diff_vector * diff_vector))

ins1 = AnalyzeInteractionPdb("2dfv")
ins2 = AnalyzeInteractionPdb("3gfb")

ins2.analyze_interactions()