class AnalyzeHydrophobicity:

    def __init__(self, amino_acid_sequence):

        self.amino_acid_sequence = amino_acid_sequence

        self.values = {
            'A': 0.267, 'C': 1.806, 'D': -2.303, 'E': -2.442,
            'F': 0.427, 'G': 0.160, 'H': -2.189, 'I': 0.971,
            'K': -2.996, 'L': 0.623, 'M': 0.136, 'N': -1.988,
            'P': -0.451, 'Q': -1.814, 'R': -2.749, 'S': -0.119,
            'T': -0.083, 'V': 0.721, 'W': -0.875, 'Y': 0.386
        }
    
    def calculate_hydrophobicity(self):

        # 合計値を計算
        total = sum(self.values[aa] for aa in self.amino_acid_sequence if aa in self.values)

        return total
    