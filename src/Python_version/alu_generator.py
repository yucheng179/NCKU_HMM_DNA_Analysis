# alu_generator.py
import numpy as np

# 範例 Alu 序列（部分）：來自人類 Alu consensus
ALU_SEQUENCE = "AGCTCGGAGGCTGAGGCAGGAGAATCGCTTGAACCCGGGAGGCGGAGGTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAA"

# 映射表：將 ACGT 轉換為數字（0,1,2,3）
BASE_TO_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

def generate_alu_sequence():
    """
    將 ALU_SEQUENCE 轉換成對應的觀測值序列（0,1,2,3）
    """
    return [BASE_TO_INDEX[base] for base in ALU_SEQUENCE if base in BASE_TO_INDEX]

def get_alu_emission_distribution():
    """
    回傳 Alu 狀態的 emission 機率分布（根據 ALU_SEQUENCE 統計）
    """
    counts = np.zeros(4)
    for base in ALU_SEQUENCE:
        if base in BASE_TO_INDEX:
            counts[BASE_TO_INDEX[base]] += 1
    return counts / np.sum(counts)
