#!/bin/bash

# Step 0: 清空改善率表
> improvement_table.txt

# Step 1: 編譯 extract.c
echo "🔧 Compiling extract.c..."
gcc -Wall -O2 -o ./data/extract ./data/extract.c
if [ $? -ne 0 ]; then
    echo "❌ Failed to compile extract.c"
    exit 1
fi

# Step 2: 編譯 Markov 模型
echo "🔧 Compiling Markov models..."
make
if [ $? -ne 0 ]; then
    echo "❌ make (Markov models) failed"
    exit 1
fi

# Step 3: 編譯 HMM project
echo "🔧 Compiling HMM_Project/hmm_project..."
cd HMM_Project
make
if [ $? -ne 0 ]; then
    echo "❌ make (HMM_Project) failed"
    cd ..
    exit 1
fi
cd ..

#for i in {0..2}  # 外層：10 個不同起始位置
#do
    start_pos=$((80000 * 0))
    echo "==============================="
    echo "🧪 Iteration $i (start_pos = $start_pos)"
    echo "==============================="

    # 清空本輪的改善率表
    > improvement_table.txt

    # 對每個資料長度執行實驗
    lengths=(5000 10000 20000 30000 40000 50000 60000 70000 80000)
    for len in "${lengths[@]}"
    do
        echo "🔍 Running for length: $len"

        # 擷取訓練資料
        ./data/extract "$len" "$start_pos"
        if [ $? -ne 0 ]; then
            echo "❌ extract failed for len=$len, start=$start_pos"
            continue
        fi

        # 執行三個 Markov 模型
        make run_zero
        make run_first
        make run_second

        cd HMM_Project
        # 執行 HMM 模型
        ./hmm_project
        cd ..
        echo "✅ Finished length: $len"
        echo
    done

    # 儲存本次迴圈的改善率表
    mv improvement_table.txt "improvement_table_${i}.txt"

done

echo "🎯 All 10 iterations completed. Check improvement_table_*.txt"
