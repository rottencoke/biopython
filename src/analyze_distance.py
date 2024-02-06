from Bio.PDB import PDBParser
import numpy as np

class AnalyzeDistances:

    def __init__(self, pdb_file):
        self.pdb_file = "pdb/0205/" + pdb_file + ".pdb"
    
    # PDBファイルからたんぱくの立体構造を解析して、各残基の不斉炭素間の距離を行列にする
    def calculate_distance_matrix(self):
        parser = PDBParser()
        structure = parser.get_structure("protein", self.pdb_file)

        # ファイルの情報読み取り
        pdb_molecule = structure.header['compound']['1']['molecule']
        # pdb_ec = structure.header['compound']['1']['ec']
        pdb_organism_scientific = structure.header['source']['1']['organism_scientific']
        pdb_organism_taxid = structure.header['source']['1']['organism_taxid']
        # pdb_strain = structure.header['source']['1']['strain']

        print(f"解析PDBファイル : {pdb_molecule} in {pdb_organism_scientific} ({pdb_organism_taxid}txid)")

        model = next(structure.get_models())  # 最初のモデルを使用

        # residues = [residue for residue in model.get_residues() if residue.has_id('CA')]
        residues = [residue for residue in model.get_residues() if residue.id[0] == ' ']
        num_residues = self.__count_seqres_amino_acids()
        distance_matrix = np.zeros((num_residues, num_residues))
        print(f"作成した行列の長さ : {len(distance_matrix)}, 配列の長さ : {num_residues}")

        for i, residue1 in enumerate(residues):
            for j, residue2 in enumerate(residues):
                distance = residue1['CA'] - residue2['CA']
                distance_matrix[i, j] = distance

        return distance_matrix

    def convert_distance_matrix(self, distance_matrix, threshold=5.0):
        # 5 Å 以下の距離を持つ要素を1に、それ以外を0に変換
        converted_matrix = np.where(distance_matrix <= threshold, 1, "_")
        return converted_matrix
    
    def __count_seqres_amino_acids(self):
        # アミノ酸の総数を保持する変数
        total_amino_acids = 0
        # PDBファイルを開く
        with open(self.pdb_file, 'r') as file:
            for line in file:
                # SEQRESレコードを見つける
                if line.startswith("SEQRES"):
                    # SEQRESレコードからアミノ酸の数を抽出（スペースで分割し、11番目の要素以降がアミノ酸）
                    amino_acids = line.split()[4:]  # アミノ酸は4番目の要素から始まります
                    total_amino_acids += len(amino_acids)
        return total_amino_acids