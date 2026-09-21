# data_loader.py
def load_data(filename):
    """
    讀取 DNA 觀測序列並轉換為數字索引
    A -> 0, C -> 1, G -> 2, T -> 3
    """
    base_mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    with open(filename, "r") as file:
        sequence = file.read().strip()

    # 將 DNA 序列轉換為數字索引
    observations = [base_mapping[base] for base in sequence if base in base_mapping]
    
    return observations
