def count_acgt_from_file(filepath):
    with open(filepath, 'r') as file:
        sequence = file.read().strip().upper()

    total = len(sequence)
    count_A = sequence.count('A')
    count_C = sequence.count('C')
    count_G = sequence.count('G')
    count_T = sequence.count('T')

    # 計算比例
    ratio_A = count_A / total
    ratio_C = count_C / total
    ratio_G = count_G / total
    ratio_T = count_T / total

    print(f"總長度：{total}")
    print(f"A: {count_A} ({ratio_A:.3%})")
    print(f"C: {count_C} ({ratio_C:.3%})")
    print(f"G: {count_G} ({ratio_G:.3%})")
    print(f"T: {count_T} ({ratio_T:.3%})")

    return {
        'A': ratio_A,
        'C': ratio_C,
        'G': ratio_G,
        'T': ratio_T
    }

# 使用方式
if __name__ == "__main__":
    ratios = count_acgt_from_file("train_data.txt")
