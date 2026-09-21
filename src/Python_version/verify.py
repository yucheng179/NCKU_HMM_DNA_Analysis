import numpy as np
from model import A, B, pi

def load_verify_data(filepath):
    """
    從 verify_data.txt 載入序列（只接受 A, C, G, T 字母）
    """
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    with open(filepath, "r") as f:
        content = f.read().strip().upper()
        return [mapping[char] for char in content if char in mapping]

def forward_probability(observations, A, B, pi):
    """
    使用 Forward 演算法計算觀察序列的生成機率
    """
    N = len(pi)
    T = len(observations)

    alpha = np.zeros((T, N))
    alpha[0] = pi * B[:, observations[0]]
    
    for t in range(1, T):
        for j in range(N):
            alpha[t, j] = np.sum(alpha[t - 1] * A[:, j]) * B[j, observations[t]]
    
    return np.sum(alpha[-1])

# 主程式
if __name__ == "__main__":
    observations = load_verify_data("verify_data.txt")
    if len(observations) == 0:
        print("Error: verify_data.txt is empty or contains invalid characters.")
        exit()

    print(f"Loaded {len(observations)} observations from verify_data.txt")

    hmm_prob = forward_probability(observations, A, B, pi)
    natural_prob = (1 / 4) ** len(observations)

    print(f"\nHMM Model Generated Probability: {hmm_prob:.10e}")
    print(f"Natural Random Generation Probability: {natural_prob:.10e}")
    print(f"Log2 Probability Ratio (HMM vs Random): {np.log2(hmm_prob / natural_prob):.4f}")
