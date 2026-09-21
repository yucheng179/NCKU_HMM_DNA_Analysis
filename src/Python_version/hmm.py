import numpy as np

import numpy as np

class HMM:
    def __init__(self, num_states, num_observations):
        self.num_states = num_states
        self.num_observations = num_observations

        # 初始化狀態機率 π (避免過度偏向某個狀態)
        self.pi = np.array([0.5, 0.5])  # 設定為均等機率
        self.pi = np.clip(self.pi, 0.05, 0.95)  # 限制範圍
        self.pi /= np.sum(self.pi)  # 確保正規化

        # 轉移機率矩陣 A (讓 HMM 不會頻繁切換狀態)
        self.A = np.array([[0.95, 0.05],
                           [0.05, 0.95]])  # 這樣 HMM 比較傾向維持同一個狀態
        self.A /= np.sum(self.A, axis=1, keepdims=True)  # 確保正規化

        # 發射機率矩陣 B (讓每個狀態的碱基比例合理)
        self.B = np.array([[0.05, 0.45, 0.45, 0.05],  # high CG 狀態均等發射
                           [0.45, 0.05, 0.05, 0.45]])  # low CG 狀態 A 和 T 比較多
        self.B /= np.sum(self.B, axis=1, keepdims=True)  # 確保正規化

        # Debug
        print(f"Initialized pi:\n{self.pi}")
        print(f"Initialized A:\n{self.A}")
        print(f"Initialized B:\n{self.B}")



    def update_parameters(self, gamma, xi, observations, epsilon=1e-10):
        """ M-Step: 更新 HMM 參數，並檢查數值異常 """
        
        # 確保 observations 是 NumPy 陣列
        observations = np.asarray(observations)  

        # Debug: Print gamma 及 xi 的範圍
        print(f"Gamma range: min={np.min(gamma)}, max={np.max(gamma)}")
        print(f"Xi range: min={np.min(xi)}, max={np.max(xi)}")

        # 更新初始機率 π
        pi_denominator = np.sum(gamma[0]) + epsilon
        self.pi = gamma[0] / pi_denominator
        print(f"Updated pi: {self.pi}")

        # 更新轉移機率 A
        denominator_A = np.sum(gamma[:-1], axis=0, keepdims=True).T + epsilon
        self.A = np.sum(xi, axis=0) / denominator_A
        
        # 平滑調整，讓 HMM 更傾向維持同一狀態
        self.A = self.A * 0.9 + np.array([[0.95, 0.05], 
                                            [0.05, 0.95]]) * 0.1  
        
        # 限制範圍，避免機率太極端
        self.A = np.clip(self.A, 0.05, 0.95)
        
        # 再次正規化，確保每列加總為 1
        self.A /= np.sum(self.A, axis=1, keepdims=True)
        
        print(f"Updated A:\n{self.A}")
        # 更新發射機率 B
        epsilon = 1e-6
        for i in range(self.num_states):
            for k in range(self.num_observations):
                mask = (observations == k).astype(float)
                denominator_B = np.sum(gamma[:, i]) + epsilon
                self.B[i, k] = (np.sum(gamma[:, i] * mask) + epsilon) / (denominator_B + epsilon)


        self.B /= np.sum(self.B, axis=1, keepdims=True)  # 確保正規化
        print(f"Updated B:\n{self.B}")
