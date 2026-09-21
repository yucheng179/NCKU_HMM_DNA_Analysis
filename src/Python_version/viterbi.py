# viterbi.py
import numpy as np

def viterbi_algorithm(hmm, observations):
    """ 找到最可能的隱狀態序列 """
    T = len(observations)
    N = hmm.num_states
    delta = np.zeros((T, N))
    psi = np.zeros((T, N), dtype=int)

    # 初始化
    delta[0] = hmm.pi * hmm.B[:, observations[0]]

    # 遞推
    for t in range(1, T):
        for j in range(N):
            delta[t, j] = np.max(delta[t-1] * hmm.A[:, j]) * hmm.B[j, observations[t]]
            psi[t, j] = np.argmax(delta[t-1] * hmm.A[:, j])

    # 回溯
    states = np.zeros(T, dtype=int)
    states[-1] = np.argmax(delta[-1])
    for t in range(T-2, -1, -1):
        states[t] = psi[t+1, states[t+1]]

    return states
