import numpy as np

# self.obj_blast_result
# self.matrix_distance_raw (PDBファイルで解析したタンパク質の各アミノ残基間の距離を行列にしたもの)
# self.matrix_distance (読み込んだ行列をquery_sequenceをもとに変更したもの)
# self.hydrophobic_amino_acids (疎水性相互作用を行うアミノ酸の種類)
# self.query_from
class AnalyzeInteractions:
    def __init__(self, obj_blast_result, matrix_distance_raw):
        self.obj_blast_result = obj_blast_result
        self.matrix_distance_raw = matrix_distance_raw
        self.__initialize_matrix_distance()

        # 疎水性相互作用を作るアミノ酸
        self.hydrophobic_amino_acids = {'W', 'F', 'M', 'L', 'I', 'V'}

    # blast_resultを読み込んで、距離の行列を変更する
    def __initialize_matrix_distance(self):

        # query_sequenceの読み込み
        arr_query_sequence = self.obj_blast_result['arr_query_sequence']

        # query_fromの読み込み
        self.query_from = int(self.obj_blast_result['query_from']) - 1

        print(f"query from : {self.query_from}")

        # 操作する行列
        matrix_distance_edited = self.matrix_distance_raw

        # 削除した位置を記録する配列
        arr_pos_insert = []

        # query_sequenceの長さで繰り返し、欠失があればその位置に99で埋まった行と列を追加
        for index, query_ac in enumerate(arr_query_sequence):

            if query_ac == '-':
                
                # 行列上の位置を計算
                pos = self.query_from + index + len(arr_pos_insert)

                value_to_add = 99

                arr_add = np.full(len(matrix_distance_edited), value_to_add)
                matrix_distance_edited = np.insert(matrix_distance_edited, pos + 1, arr_add, axis=0)

                arr_add = np.full(len(matrix_distance_edited), value_to_add)
                matrix_distance_edited = np.insert(matrix_distance_edited, pos + 1, arr_add, axis=1)

                # 追加した位置を保存
                arr_pos_insert.append(pos) 
        
        print(f"query側の欠失の位置 : {arr_pos_insert}")

        # 距離の行列を定義
        self.matrix_distance = matrix_distance_edited

        print(f"行列の長さ{len(self.matrix_distance)}")

    # 指定されたアミノ酸の組み合わせの相互作用をカウントする
    def count_interactions(self, amino_acids, threshold):

        arr_hit_sequence = self.obj_blast_result['arr_hit_sequence']
        print(f"hit_sequenceの長さ{len(arr_hit_sequence)}")
        count = 0
        for i in range(len(arr_hit_sequence)):
            for j in range(i+1, len(arr_hit_sequence)):
                # 距離行列で定義された閾値以下の距離にある組み合わせを数える
                if self.matrix_distance[i + self.query_from][j + self.query_from] <= threshold:
                    if (arr_hit_sequence[i], arr_hit_sequence[j]) in amino_acids or (arr_hit_sequence[j], arr_hit_sequence[i]) in amino_acids:
                        count += 1
        return count
    
    # 疎水性相互作用を持つアミノ酸ペアの数をカウントする（改変版）
    def count_interactions_hydrophobic(self):
        arr_hit_sequence = self.obj_blast_result['arr_hit_sequence']
        count = 0
        for i in range(len(arr_hit_sequence)):
            for j in range(i+1, len(arr_hit_sequence)):
                # アミノ酸ごとの閾値を設定
                threshold = self.__get_threshold(arr_hit_sequence[i], arr_hit_sequence[j])

                # 距離行列で定義された閾値以下の距離にある疎水性アミノ酸を数える
                if self.matrix_distance[i + self.query_from][j + self.query_from] <= threshold:
                    if arr_hit_sequence[i] in self.hydrophobic_amino_acids and arr_hit_sequence[j] in self.hydrophobic_amino_acids:
                        count += 1
        return count
    
    # 疎水性相互作用を持つアミノ酸残基の総数をカウントする（改変版）
    def count_hydrophobic_residues_unique(self):
        arr_hit_sequence = self.obj_blast_result['arr_hit_sequence']
        hydrophobic_residues_counted = [False] * len(arr_hit_sequence)

        for i in range(len(arr_hit_sequence)):
            for j in range(i+1, len(arr_hit_sequence)):
                threshold = self.__get_threshold(arr_hit_sequence[i], arr_hit_sequence[j])
                if self.matrix_distance[i + self.query_from][j + self.query_from] <= threshold:
                    if arr_hit_sequence[i] in self.hydrophobic_amino_acids and arr_hit_sequence[j] in self.hydrophobic_amino_acids:
                        hydrophobic_residues_counted[i] = True
                        hydrophobic_residues_counted[j] = True

        count_unique_hydrophobic_residues = sum(hydrophobic_residues_counted)
        return count_unique_hydrophobic_residues
    
    def identify_hydrophobic_cores(self):
        arr_hit_sequence = self.obj_blast_result['arr_hit_sequence']
        matrix_distance = self.matrix_distance
        hydrophobic_cores = []

        for i, aa in enumerate(arr_hit_sequence):
            if aa in self.hydrophobic_amino_acids and aa != 'M':  # 開始コドンのMは除外
                new_core = True
                for core in hydrophobic_cores:
                    if any(self.__get_threshold(aa, arr_hit_sequence[j]) >= matrix_distance[i][j] for j in core):
                        core.add(i)
                        new_core = False
                        break
                if new_core:
                    hydrophobic_cores.append({i})

        # アミノ酸の数が3つ未満のコアを削除
        hydrophobic_cores = [core for core in hydrophobic_cores if len(core) >= 3]

        # コアの詳細を含む辞書を作成
        hydrophobic_cores_details = [
            {'positions': core, 'amino_acids': [arr_hit_sequence[i] for i in core]} for core in hydrophobic_cores
        ]

        number_of_cores = len(hydrophobic_cores)
        cores_counts = self.__count_amino_acids_in_cores(hydrophobic_cores_details)
        return number_of_cores, cores_counts

    # 閾値を決定するヘルパーメソッド（改変版）
    def __get_threshold(self, amino_acid1, amino_acid2):
        set_fw = {'F', 'W'}
        set_vil = {'V', 'I', 'L'}

        # (F, W)のどちらかと(V, I, L)のどちらかの組み合わせ、または、Mを含む場合の閾値設定
        if (amino_acid1 in set_fw and amino_acid2 in set_vil) or (amino_acid2 in set_fw and amino_acid1 in set_vil) or 'M' in {amino_acid1, amino_acid2}:
            return 7
        elif amino_acid1 in set_fw or amino_acid2 in set_fw:
            return 8
        elif amino_acid1 in set_vil or amino_acid2 in set_vil:
            return 7
        else:
            # ここにはデフォルトの閾値を設定する
            return 0
        
    def __count_amino_acids_in_cores(self, hydrophobic_cores):
        aromatic_amino_acids = {'F', 'W', 'Y', 'H'}
        aliphatic_amino_acids = {'L', 'I', 'V', 'M'}

        cores_counts = []
        for core in hydrophobic_cores:
            aromatic_count = sum(aa in aromatic_amino_acids for aa in core['amino_acids'])
            aliphatic_count = sum(aa in aliphatic_amino_acids for aa in core['amino_acids'])
            core_count = {'aromatic_count': aromatic_count, 'aliphatic_count': aliphatic_count}
            cores_counts.append(core_count)

        return cores_counts