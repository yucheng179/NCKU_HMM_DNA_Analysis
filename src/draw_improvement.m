% draw_improvement.m

% 搜尋 improvement_table_*.txt 檔案
files = dir('improvement_table_*.txt');

% 定義長度與模型名稱
lengths = [5000, 10000, 20000, 40000, 80000];
models = {'Zero-order', 'First-order', 'Second-order', 'Baum-Welch HMM'};

for idx = 1:length(files)
    filename = files(idx).name;
    data = readmatrix(filename);

    % 檢查格式
    if size(data,1) ~= length(lengths) || size(data,2) ~= length(models)
        warning("⚠️ Skipping %s: unexpected shape (%dx%d)", filename, size(data,1), size(data,2));
        continue;
    end

    % 繪圖
    figure;
    hold on;
    for i = 1:size(data,2)
        plot(lengths, data(:,i), '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
            'DisplayName', models{i});
    end

    xlabel('Sequence Length (n)');
    ylabel('Improvement Rate (log_{10} scale)');
    title(sprintf('Improvement Rate vs Length — Trial %d', idx-1));
    legend('Location', 'northeast');
    grid on;
    xticks(lengths);
    xticklabels(string(lengths));
    set(gca, 'YScale', 'log');

    % 儲存圖（覆蓋命名）
    outname = sprintf('improvement_plot_%d.png', idx-1);  % 因為 idx-1 是 trial 編號
    saveas(gcf, outname);
    close;  % 關掉當前圖，避免太多視窗
end
