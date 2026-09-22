# 隱馬可夫模型於 DNA 序列之分析

> Hidden Markov Models for Analysis of DNA Sequence Data

國立成功大學資訊工程學系畢業專題。此研究以隱馬可夫模型（Hidden Markov Model, HMM）描述人類 DNA 序列的統計結構，從可觀察的鹼基序列 `A / C / G / T` 推測不可直接觀察的區域狀態，並比較 HMM 與 0、1、2 階馬可夫模型對測試序列的解釋能力。

- **指導教授：** 賀保羅教授
- **專題成員：** 曾暐翔、謝承佑
- **資料來源：** NCBI RefSeq - *Homo sapiens* chromosome 6 標準參考序列
- **主要語言：** C、Python、MATLAB
- **開發／測試環境：** Visual Studio Code、Microsoft Windows

## 研究動機

DNA 並非由四種鹼基均勻且彼此獨立地隨機產生。不同基因組區域可能呈現不同的 GC 含量、轉移模式與重複序列結構，因此單一的全域機率分布不一定足以描述整條序列。

本專題的核心假設是：觀察到的鹼基由一組隨序列位置改變、但無法直接看到的狀態所生成。HMM 可同時表示「區域之間如何切換」與「每個區域傾向產生哪些鹼基」，適合用來探索 DNA 中的區段性統計差異。

## 研究目標

1. 建立可處理 `A / C / G / T` 離散觀測值的 HMM。
2. 以 Viterbi 演算法推測最可能的隱藏狀態路徑。
3. 以 Forward-Backward 與 Baum-Welch（EM）演算法估計模型參數。
4. 針對長 DNA 序列處理機率下溢與訓練穩定性問題。
5. 將 HMM 與隨機模型、0 階、1 階及 2 階馬可夫模型進行比較。
6. 探索 GC 結構與 Alu 重複序列對模型辨識能力的影響。

## 模型設計

### 觀測值與隱藏狀態

| 類型 | 設計 |
|---|---|
| 觀測值 | 四種 DNA 鹼基：A、C、G、T |
| High | GC 含量較高的區域 |
| Mid | GC 含量介於高、低之間的區域 |
| Low | GC 含量較低的區域 |
| Alu | 由外部規則辨識的 Alu-like 重複序列區段 |

模型以初始狀態機率 `π`、狀態轉移矩陣 `A` 與發射機率矩陣 `B` 表示。High、Mid、Low 為可由資料學習的隱藏狀態；Alu 則是結合 AluY、AluS、AluJ 共識序列與近似字串比對所加入的外部決定性狀態。這是一種混合式設計，而非完全由無監督 HMM 自行辨識 Alu。

### HMM 的三個基本問題

本專題實作對應經典 HMM 架構的三項核心問題：

| 問題 | 本專題採用的方法 | 用途 |
|---|---|---|
| Evaluation | Forward / Forward-Backward | 計算觀測序列在模型下的 likelihood |
| Decoding | Viterbi | 尋找最可能的隱藏狀態路徑 |
| Learning | Baum-Welch / EM | 反覆更新 `π`、`A`、`B`，提升訓練序列的 likelihood |

### 數值穩定性

長序列中大量小於 1 的機率連乘會快速趨近 0。為避免浮點數下溢，主要 C 實作使用 log-domain、log-sum-exp 與極小值 `EPSILON` 進行計算；這也呼應經典 HMM 實作中以 scaling 或 logarithm 維持數值範圍的做法。

## 系統流程

```mermaid
flowchart LR
    A[讀取訓練與驗證 DNA] --> B[Viterbi 初步解碼]
    B --> C[標記 Alu-like 區段]
    C --> D[頻率統計與參數初始化]
    D --> E[Baum-Welch / EM 訓練]
    E --> F[Forward likelihood 驗證]
    F --> G[與隨機及 0/1/2 階模型比較]
    G --> H[狀態比例、轉移熵與圖表輸出]
```

1. **Read Data：** 載入訓練與驗證 DNA 序列。
2. **Estimation：** 使用 Viterbi 進行初始狀態分配，並以頻率統計初始化參數。
3. **Alu annotation：** 以 Alu 共識序列的近似比對標記外部決定性狀態。
4. **Training：** 使用 Baum-Welch 進行多輪 EM 訓練。
5. **Verification：** 以 Forward likelihood 評估驗證序列。
6. **Comparison：** 計算相對均勻隨機模型 `(1/4)^n` 的 improvement rate，並與低階馬可夫模型比較。

## 評估方式

若驗證序列長度為 `n`，均勻隨機模型的生成機率為：

```text
P_random(x) = (1/4)^n
```

模型提升倍數定義為：

```text
Improvement rate = P_model(x) / P_random(x)
```

實作以 log probability 計算，避免直接處理極小機率。除 HMM 外，repository 亦保留 0、1、2 階馬可夫模型，作為由簡至繁的統計基準。

## 實驗結果與觀察

- 在具有明顯 GC-rich 與 GC-poor 區域差異的測試序列上，HMM 對序列結構的解釋能力較佳。
- 當資料含有真實或接近共識序列的 Alu 片段時，加入 Alu 資訊可提高模型對該類序列的生成機率。
- 多組測試中，HMM 表現多半介於 1 階與 2 階馬可夫模型之間，部分資料高於 2 階模型。
- 結果會受到序列長度、資料切分、初始參數、隱藏狀態設計與 Alu 標記方式影響，因此目前成果應視為方法可行性驗證，而非完成生物註解工具的臨床或基因體學基準。

完整圖表與說明請參閱 [`project-report.pdf`](./project-report.pdf) 與 [`project-summary.pdf`](./project-summary.pdf)。

## 實作過程中的問題

| 問題 | 處理方式 |
|---|---|
| 狀態結構不明顯 | 調整 High／Mid／Low 狀態設計，並加入外部 Alu 標註 |
| EM 落入局部最佳解 | 嘗試不同初始參數；以 Viterbi 與頻率統計改善初始化 |
| 模型過擬合 | 調整隱藏狀態數量並進行多次初始化 |
| 長序列機率下溢 | 使用 log-domain 與 log-sum-exp |
| 稀疏事件造成零機率 | 加入 smoothing、最小轉移機率與 `EPSILON` |
| 訓練資料有限 | 控制模型複雜度，並將結果定位為探索性研究 |

Baum-Welch 只能保證 likelihood 不下降並收斂至臨界點，不能保證找到全域最佳解；狀態數量與初始值的選擇因此是本研究的重要限制。

## Repository 結構

```text
.
├── project-report.pdf             # 專題成果海報
├── project-summary.pdf            # 兩頁專題簡介
└── src/
    ├── HMM_Project/               # 主要 C 語言 HMM 實作
    │   ├── Version_1 ... Version_5/  # 開發歷程快照
    │   ├── hmm.c / hmm.h          # Viterbi、Forward-Backward、Baum-Welch
    │   ├── io.c / io.h            # DNA 序列讀取
    │   ├── matrix_utils.*         # 機率矩陣工具
    │   └── sequence_generation.*  # HMM 序列生成
    ├── Python_version/            # Python 原型與小型測試資料
    ├── web/                       # Flask 操作介面原型
    ├── zero_order_Markov.c        # 0 階基準
    ├── first_order_Markov.c       # 1 階基準
    ├── second_order_Markov.c      # 2 階基準
    ├── draw_improvement.m         # 實驗圖表繪製
    └── run_all.sh                 # 批次執行模型
```

## 執行方式

### 1. 準備資料

大型染色體資料未收錄於 GitHub。請將只包含 `A/C/G/T` 的訓練與驗證序列放在：

```text
src/data/train_data.txt
src/data/vertified_data.txt
```

> `vertified_data.txt` 是目前原始碼沿用的檔名。若要更名為 `verified_data.txt`，需同步修改程式中的路徑。

### 2. 編譯主要 C 版本

需要 GCC、Make 與 math library：

```bash
cd src/HMM_Project
make
./hmm_project
```

Windows 使用 MinGW 時可改用 `mingw32-make`，輸出程式可能為 `hmm_project.exe`。

### 3. 編譯基準模型

```bash
cd src
make
./zero_order
./first_order
./second_order
```

### 4. Python 原型

```bash
cd src/Python_version
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
pip install numpy
python train_model.py
python verify.py
```

### 5. Web 介面原型

`src/web/app.py` 會呼叫 `src/web/models/` 中已編譯的模型執行檔。需先自行編譯並放置執行檔，再執行：

```bash
pip install flask
cd src/web
python app.py
```

## 限制與後續方向

- 目前狀態標籤主要代表統計上的 GC 差異，不等同於經生物實驗驗證的功能性基因區段。
- Alu 狀態仰賴既定共識序列與字串距離，並非完全由 HMM 自動學得。
- 尚未建立跨染色體、跨物種或標準註解資料集上的完整 precision／recall 評估。
- 可加入多條獨立訓練序列、交叉驗證、更多初始化與模型選擇準則，降低單次資料切分偏差。
- 可與 CpG island、RepeatMasker 或其他基因註解結果比較，進一步驗證狀態的生物意義。
- 可將資料路徑、狀態數與 EM 超參數改為命令列設定，提升可重現性。

## 參考資料

- Lawrence R. Rabiner, “A Tutorial on Hidden Markov Models and Selected Applications in Speech Recognition,” *Proceedings of the IEEE*, vol. 77, no. 2, 1989.
- NCBI Reference Sequence (RefSeq), *Homo sapiens* chromosome 6 reference sequence.
- 本專題成果文件：[`project-report.pdf`](./project-report.pdf)、[`project-summary.pdf`](./project-summary.pdf)。

---

# Hidden Markov Models for DNA Sequence Analysis

This NCKU Computer Science graduation project investigates whether a discrete Hidden Markov Model can capture regional statistical structure in human DNA. The model observes `A/C/G/T` bases and uses High-, Mid-, and Low-GC hidden states together with an externally annotated Alu state.

The implementation addresses the three classical HMM problems: Forward evaluation, Viterbi decoding, and Baum-Welch/EM parameter learning. The primary C version performs calculations in the log domain to reduce numerical underflow on long sequences. Its likelihood is compared with a uniform random baseline and zero-, first-, and second-order Markov models.

Experiments use human chromosome 6 reference data from NCBI RefSeq. The observed results suggest that HMM performance improves when a sequence contains clear GC-rich/GC-poor regions or recognizable Alu repeats. Across the reported trials, the HMM generally falls between first- and second-order Markov baselines and exceeds the second-order model on some samples. These findings demonstrate feasibility rather than a production-grade genome annotation benchmark.

The repository contains the main C implementation, Python prototypes, development snapshots, baseline Markov models, a MATLAB plotting script, and a Flask interface prototype. Large genomic datasets, virtual environments, compiled executables, and build artifacts are intentionally excluded. See the Chinese sections above for setup instructions, methodology, limitations, and the complete repository map.

This is a team project advised by Professor Paul Ho and completed by Wei-Hsiang Tseng and Cheng-Yu Hsieh.
