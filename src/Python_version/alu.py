# alu.py
import random

class AluGenerator:
    def __init__(self):
        # Alu 共通序列元件（這段可再擴充，或從多個 Alu 中挑出 conserved 的部分）
        self.alu_sequence = list("AGCTCCGGTGCGTGACGCGGAGTCGAGGACGGCAGGCGGAGGAGTGGGCGCGG")

    def generate(self):
        # 回傳 Alu 序列轉成 observation index（A=0, C=1, G=2, T=3）
        base_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        return [base_to_index.get(base, 0) for base in self.alu_sequence]

    def length(self):
        return len(self.alu_sequence)


if __name__ == '__main__':
    generator = AluGenerator()
    print(generator.generate())
