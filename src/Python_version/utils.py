# utils.py
import numpy as np

def log_sum_exp(arr):
    """ 計算 log(∑exp(x_i)) 來避免數值溢出 """
    max_val = np.max(arr)
    return max_val + np.log(np.sum(np.exp(arr - max_val)))

def normalize(matrix):
    """ 正規化矩陣，使每列總和為 1 """
    return matrix / np.sum(matrix, axis=1, keepdims=True)

def initialize_hmm(num_states, num_observations):
    """ 隨機初始化 HMM 參數 """
    np.random.seed(42)  # 固定種子，確保結果可重現
    pi = np.random.rand(num_states)
    A = np.random.rand(num_states, num_states)
    B = np.random.rand(num_states, num_observations)

    # 正規化
    pi /= np.sum(pi)
    A = normalize(A)
    B = normalize(B)

    return pi, A, B
