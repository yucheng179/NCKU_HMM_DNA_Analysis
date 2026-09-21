# forward_backward.py
import numpy as np

def forward_algorithm(hmm, observations):
    num_obs = len(observations)
    alpha = np.zeros((num_obs, hmm.num_states))

    # 初始化 alpha
    alpha[0] = hmm.pi * hmm.B[:, observations[0]]
    alpha[0] /= np.sum(alpha[0]) + 1e-10  # 避免溢出
    
    print(f"Alpha[0]: {alpha[0]}")

    # 遞推 alpha
    for t in range(1, num_obs):
        for j in range(hmm.num_states):
            alpha[t, j] = np.sum(alpha[t-1] * hmm.A[:, j]) * hmm.B[j, observations[t]]
        
        alpha[t] /= np.sum(alpha[t]) + 1e-10  # 正規化
        #if t % 100 == 0:  # 每 100 步列印一次
            #print(f"Alpha[{t}]: {alpha[t]}")

    return alpha

def backward_algorithm(hmm, observations):
    num_obs = len(observations)
    beta = np.zeros((num_obs, hmm.num_states))

    # 初始化 beta
    beta[-1] = 1
    print(f"Beta[-1]: {beta[-1]}")

    # 遞推 beta
    for t in range(num_obs - 2, -1, -1):
        for i in range(hmm.num_states):
            beta[t, i] = np.sum(hmm.A[i] * hmm.B[:, observations[t+1]] * beta[t+1])

        beta[t] /= np.sum(beta[t]) + 1e-10  # 正規化
        #if t % 100 == 0:
            #print(f"Beta[{t}]: {beta[t]}")

    return beta


def compute_gamma_xi(hmm, observations, alpha, beta):
    gamma = (alpha * beta) / (np.sum(alpha * beta, axis=1, keepdims=True) + 1e-10)  # 避免除以 0

    xi = np.zeros((len(observations) - 1, hmm.num_states, hmm.num_states))
    for t in range(len(observations) - 1):
        denominator = np.sum(alpha[t] * hmm.A * hmm.B[:, observations[t + 1]] * beta[t + 1]) + 1e-10
        xi[t] = (alpha[t][:, None] * hmm.A * hmm.B[:, observations[t + 1]] * beta[t + 1]) / denominator

    # 檢查 gamma 是否過於小
    if np.all(gamma < 1e-6):
        print("Warning: gamma is too small, model might not be learning.")
    
    return gamma, xi

