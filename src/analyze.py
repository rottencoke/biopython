from analyze_distance import AnalyzeDistances
from analyze_interaction import AnalyzeInteractions
from analyze_hydrophobicity import AnalyzeHydrophobicity
from datetime import datetime
import csv
import numpy as np

# PDBファイル名とblast_resultファイル名
file_pdb = "1bkg"
file_blast_result = "aspartate aminotransferase_20240130_170132"

# PDBファイルを指定して、AnalyzeDistancesのインスタンスを作成
ins_distance = AnalyzeDistances(file_pdb)

# 距離行列の計算
matrix_distance = ins_distance.calculate_distance_matrix()

# csvに保存するデータをリストで作成
matrix_interactions = []

# アミノ酸の組み合わせと距離の閾値
amino_pair_disulfide = {('C', 'C')}
amino_pair_ionic = {('R', 'D'), ('R', 'E'), ('K', 'D'), ('K', 'E'), ('H', 'D'), ('H', 'E')}
threshold = 5.0

# blast_resultのパス
path_blast_result = f"blast_result/{file_blast_result}.csv"

# blast_resultのcsvの対象のヘッダーの列
row_query_from = row_gaps = row_align_len = row_bit_score = 0
row_hit_sequence = row_query_sequence = row_genus_and_species = row_strain = row_topt = row_protein = ""

# blast_resultの結果の数
num_data = 0

# blast_result読み込み
with open(path_blast_result, newline='') as csvfile:
    reader = csv.reader(csvfile)
    headers = next(reader)  # ヘッダ行を読み込む

    # ヘッダ行から目的の列のインデックスを見つける
    row_genus_and_species = headers.index("genus_and_species")
    row_strain = headers.index("strain")
    row_topt = headers.index("Topt_ave[℃]")
    row_protein = headers.index("protein")
    row_bit_score = headers.index("bit_score")
    row_gaps = headers.index("gaps")
    row_align_len = headers.index("align_len")
    row_query_from = headers.index("query_from")
    row_hit_sequence = headers.index("hit_sequence")
    row_query_sequence = headers.index("query_sequence")

    # ヘッダー以降の行数（つまりデータ数）を数える
    for row in reader:
        num_data += 1

print(f"データ数 : {num_data}")

# fastaファイルの数だけ繰り返す
for index in range(num_data):

    # blast resultの情報を取得
    query_from = gaps = align_len = bit_score = 0
    hit_sequence = query_sequence = genus_and_species = strain = topt = protein = ""
    
    with open(path_blast_result, newline='') as csvfile:
        reader = csv.reader(csvfile)
        headers = next(reader)  # ヘッダ行を読み込む

        for j, row in enumerate(reader):
            if j == index:

                # 目的の行が見つかった場合、該当の列を取得
                genus_and_species = row[row_genus_and_species]
                strain = row[row_strain]
                topt = row[row_topt]
                protein = row[row_protein]
                bit_score = row[row_bit_score]
                gaps = row[row_gaps]
                align_len = row[row_align_len]
                query_from = row[row_query_from]
                hit_sequence = row[row_hit_sequence]
                query_sequence = row[row_query_sequence]

                break
    
    # リスト (配列) に
    arr_hit_sequence = list(hit_sequence)
    arr_query_sequence = list(query_sequence)

    # analyze_interactionに渡すblast_resultのオブジェクト
    obj_blast_result = {
        'arr_hit_sequence': arr_hit_sequence,
        'arr_query_sequence': arr_query_sequence,
        'query_from': query_from
    }

    print(f"{index} : {genus_and_species}")

    # 行列を元にfastaファイルで読み込んだ配列を解析するインスタンスを作成
    ins_interaction = AnalyzeInteractions(obj_blast_result, matrix_distance)

    # 各相互作用の解析
    num_hydrophobic = ins_interaction.count_interactions_hydrophobic()
    num_hydrophobic_pos = ins_interaction.count_hydrophobic_residues_unique()
    num_hydrophobic_core, hydrophobic_cores = ins_interaction.identify_hydrophobic_cores()

    result_hydrophobic = []
    for item in hydrophobic_cores:
        result_hydrophobic.append(item['aromatic_count'])
        result_hydrophobic.append(item['aliphatic_count'])

    print(f"疎水性相互作用 : {num_hydrophobic}個")
    print(f"コア数 : {num_hydrophobic_core}個")
    print(hydrophobic_cores)

    # 疎水度計算
    ins_hydrophobicity = AnalyzeHydrophobicity(hit_sequence)
    value_hydrophobicity = ins_hydrophobicity.calculate_hydrophobicity()
    
    # データを配列にまとめる
    arr_data = np.array([genus_and_species, strain, topt, protein, bit_score, align_len, gaps, num_hydrophobic, num_hydrophobic_pos, value_hydrophobicity])

    arr_result_hydrophobic = np.array(result_hydrophobic)

    arr_data_combined = np.concatenate((arr_data, arr_result_hydrophobic))

    # 配列を行列に保存
    matrix_interactions.append(arr_data_combined)


# 現在の日時を取得
now = datetime.now()

# 日時を 'mmdd_hhmmss' 形式でフォーマット
date_formatted = now.strftime("%m%d_%H%M%S")

# 作成ファイル名
name_file = f"result/analysis_{date_formatted}_{file_pdb}.csv"

# ヘッダー（各データの名前）
headers = ['genus_and_species', 'strain', 'Topt_ave[℃]', 'protein', 'bit_score', 'align_len', 'gaps', 'hydrophobic', 'pos_hydrophobic', 'hydrophobicity']

# CSVファイルに書き込む
with open(name_file, 'w', newline='', encoding='utf-8') as file:
    writer = csv.writer(file)

    # ヘッダーを書き込む
    writer.writerow(headers)

    # データを書き込む
    writer.writerows(matrix_interactions)

print(f"結果を保存 : {name_file}")