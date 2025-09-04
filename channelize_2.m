
clear; clc;

%% 1. 系统参数设置
K = 16;          % 信道数（必须=16，匹配流程图）
M = 8;           % 多相滤波器抽头数（原型长度=K×M）
fs = 1e9;        % 基带采样率（110GHz信号需先下变频）
N = K * 200;     % 输入信号长度（K的倍数，保证抽取完整）
L = N / K;       % 每条支路的有效长度（L = N/K = 200）

%% 2. 设计原型低通滤波器
h_prototype = fir1(K*M - 1, 1/K);  % 长度=16×8-1=127的FIR滤波器

%% 3. 多相分解（修正索引越界）
h_poly = zeros(K, M);  % K行×M列的多相滤波器矩阵
for m = 0:K-1          % m: 0~15（信道索引，0-based）
    for n = 0:M-1      % n: 0~7（抽头索引，0-based）
        idx = m + n*K + 1;  % 转换为1-based索引
        if idx <= length(h_prototype)
            h_poly(m+1, n+1) = h_prototype(idx);  % 存储到1-based矩阵
        else
            h_poly(m+1, n+1) = 0;  % 超出范围补0（理论上不会发生）
        end
    end
end
%% 多相分解结果可视化
figure('Name','多相分解结果可视化','Position',[100 100 1000 800]);

% 显示多相滤波器矩阵的热图
subplot(2,1,1);
heatmap(h_poly,'Title','多相滤波器矩阵h_poly','XLabel','抽头索引(n)','YLabel','信道索引(m)');
colormap('jet');
colorbar;
title('多相滤波器矩阵热图显示');

% 绘制每个信道的滤波器抽头
subplot(2,1,2);
colors = lines(K);  % 生成K种不同颜色
for m = 1:K
    stem(0:M-1, h_poly(m,:), 'MarkerFaceColor', colors(m,:), ...
         'Color', colors(m,:), 'DisplayName', sprintf('信道 %d', m));
    hold on;
end
grid on;
xlabel('抽头索引 n');
ylabel('滤波器系数值');
title('各信道的多相滤波器抽头系数');
legend('Location','bestoutside');
hold off;

% 显示数值矩阵（在命令行）
disp('多相滤波器矩阵h_poly的数值：');
disp(h_poly);

% 原型滤波器与多相分解的关系可视化
figure('Name','原型滤波器与多相分解关系','Position',[200 200 900 500]);
subplot(2,1,1);
stem(0:length(h_prototype)-1, h_prototype, 'b', 'MarkerFaceColor', 'b');
grid on;
xlabel('原型滤波器系数索引');
ylabel('系数值');
title('原型低通滤波器h_prototype');

subplot(2,1,2);
for m = 1:K
    % 提取原型滤波器中属于第m个多相分量的系数
    indices = m:K:length(h_prototype);
    stem(indices, h_prototype(indices), 'MarkerFaceColor', colors(m,:), ...
         'Color', colors(m,:), 'DisplayName', sprintf('信道 %d', m));
    hold on;
end
grid on;
xlabel('原型滤波器系数索引');
ylabel('系数值');
title('各信道在原型滤波器中对应的系数位置');
legend('Location','bestoutside');
hold off;

%% 4. 生成输入信号（16个单音，覆盖各信道中心频率）
n = 0:N-1;                % 全局时间索引（0-based，长度N）
X = zeros(1, N);          % 初始化输入信号（1×N）
for m = 0:K-1
    f_m = (m + 0.5) * fs / K;  % 信道m的中心频率（基带全带宽覆盖）
    X = X + exp(1j*2*pi*f_m/fs * n);  % 叠加单音信号
end

%% 5. 多相支路抽取（0-based→1-based索引转换）
x_branch = cell(K, 1);    % 存储16条支路的信号（每条长度L）
for m = 0:K-1
    % 抽取公式：x_m(k) = X(m + K×k)，转换为1-based索引：m+1 + K×k
    x_m = X(m + 1 + K*(0:L-1));  % 长度=L=200
    x_branch{m+1} = x_m;         % 存储为1-based细胞数组
end

%% 5. 多相支路抽取结果可视化
figure('Name','多相支路抽取时域信号','Position',[300 300 1000 800]);
sgtitle('多相支路抽取后的时域信号');

% 绘制前4条支路的时域信号（实部）
subplot(4,1,1);
for m = 1:4
    plot(0:L-1, real(x_branch{m}), 'Color', colors(m,:), ...
         'LineWidth', 1.2, 'DisplayName', sprintf('信道 %d', m));
    hold on;
end
grid on;
xlabel('样本索引');
ylabel('信号幅度（实部）');
title('1-4号信道信号');
legend('Location','best');
hold off;

% 绘制5-8条支路的时域信号（实部）
subplot(4,1,2);
for m = 5:8
    plot(0:L-1, real(x_branch{m}), 'Color', colors(m,:), ...
         'LineWidth', 1.2, 'DisplayName', sprintf('信道 %d', m));
    hold on;
end
grid on;
xlabel('样本索引');
ylabel('信号幅度（实部）');
title('5-8号信道信号');
legend('Location','best');
hold off;

% 绘制9-12条支路的时域信号（实部）
subplot(4,1,3);
for m = 9:12
    plot(0:L-1, real(x_branch{m}), 'Color', colors(m,:), ...
         'LineWidth', 1.2, 'DisplayName', sprintf('信道 %d', m));
    hold on;
end
grid on;
xlabel('样本索引');
ylabel('信号幅度（实部）');
title('9-12号信道信号');
legend('Location','best');
hold off;

% 绘制13-16条支路的时域信号（实部）
subplot(4,1,4);
for m = 13:16
    plot(0:L-1, real(x_branch{m}), 'Color', colors(m,:), ...
         'LineWidth', 1.2, 'DisplayName', sprintf('信道 %d', m));
    hold on;
end
grid on;
xlabel('样本索引');
ylabel('信号幅度（实部）');
title('13-16号信道信号');
legend('Location','best');
hold off;

% 各支路信号频谱分析
figure('Name','多相支路信号频谱','Position',[400 400 1200 800]);
sgtitle('多相支路抽取后的信号频谱');

% 计算每条支路的频谱并显示
for m = 1:K
    x = x_branch{m};
    X_fft = fftshift(fft(x));
    f_axis = (-L/2 : L/2 - 1) * (fs/(K*L));  % 频率轴（Hz）
    
    subplot(4,4,m);
    plot(f_axis/1e6, abs(X_fft), 'Color', colors(m,:), 'LineWidth', 1.2);
    grid on;
    xlabel('频率 (MHz)');
    ylabel('幅度');
    title(sprintf('信道 %d', m));
    xlim([-fs/(2*K) fs/(2*K)]/1e6);  % 限制显示该信道的频率范围
    ylim([0 max(abs(X_fft))*1.1]);
end

% 单条支路的详细分析（以第1条为例）
figure('Name','支路信号详细分析','Position',[500 500 900 600]);
m = 1;  % 选择第1条支路进行详细分析

subplot(2,1,1);
plot(0:L-1, real(x_branch{m}), 'b', 'LineWidth', 1.2);
hold on;
plot(0:L-1, imag(x_branch{m}), 'r--', 'LineWidth', 1.2);
grid on;
xlabel('样本索引');
ylabel('信号幅度');
title(sprintf('第%d条支路信号（实部：蓝色，虚部：红色虚线）', m));
legend('实部','虚部','Location','best');
hold off;

subplot(2,1,2);
x = x_branch{m};
X_fft = fftshift(fft(x));
f_axis = (-L/2 : L/2 - 1) * (fs/(K*L));
plot(f_axis/1e6, 20*log10(abs(X_fft)), 'g', 'LineWidth', 1.2);
grid on;
xlabel('频率 (MHz)');
ylabel('幅度 (dB)');
title(sprintf('第%d条支路信号频谱', m));
xlim([-fs/(2*K) fs/(2*K)]/1e6);
hold off;
%% 6. 混频（复指数搬移）
x_mixed = cell(K, 1);
for m = 0:K-1
    x_m = x_branch{m+1};
    k = 0:L-1;                  % 支路内时间索引（0-based，长度L）
    mix_factor = exp(-1j*pi/2 * m);    % 混频因子
    x_mixed{m+1} = x_m .* mix_factor;
end

%% 6. 混频结果可视化
colors = lines(K);  % 复用之前的颜色方案

% 1. 混频前后时域信号对比（实部）
figure('Name','混频前后时域信号对比','Position',[100 100 1000 800]);
sgtitle('混频前后的时域信号（实部对比）');

% 前8条信道对比
for i = 1:8
    subplot(4,2,i);
    % 原始信号
    plot(0:L-1, real(x_branch{i}), 'b--', 'LineWidth', 1, ...
         'DisplayName', '混频前');
    hold on;
    % 混频后信号
    plot(0:L-1, real(x_mixed{i}), 'Color', colors(i,:), 'LineWidth', 1.2, ...
         'DisplayName', '混频后');
    grid on;
    xlabel('样本索引');
    ylabel('信号幅度（实部）');
    title(sprintf('信道 %d', i));
    legend('Location','best');
    hold off;
end

% 后8条信道对比
figure('Name','混频前后时域信号对比（续）','Position',[200 200 1000 800]);
sgtitle('混频前后的时域信号（实部对比）');

for i = 9:16
    subplot(4,2,i-8);
    % 原始信号
    plot(0:L-1, real(x_branch{i}), 'b--', 'LineWidth', 1, ...
         'DisplayName', '混频前');
    hold on;
    % 混频后信号
    plot(0:L-1, real(x_mixed{i}), 'Color', colors(i,:), 'LineWidth', 1.2, ...
         'DisplayName', '混频后');
    grid on;
    xlabel('样本索引');
    ylabel('信号幅度（实部）');
    title(sprintf('信道 %d', i));
    legend('Location','best');
    hold off;
end

% 2. 混频前后频谱对比
figure('Name','混频前后频谱对比','Position',[300 300 1200 800]);
sgtitle('混频前后的信号频谱对比');

for m = 1:K
    % 计算混频前后的频谱
    X_fft_original = fftshift(fft(x_branch{m}));
    X_fft_mixed = fftshift(fft(x_mixed{m}));
    f_axis = (-L/2 : L/2 - 1) * (fs/(K*L));  % 频率轴（Hz）
    
    subplot(4,4,m);
    % 原始频谱
    plot(f_axis/1e6, abs(X_fft_original), 'b--', 'LineWidth', 1, ...
         'DisplayName', '混频前');
    hold on;
    % 混频频谱
    plot(f_axis/1e6, abs(X_fft_mixed), 'Color', colors(m,:), 'LineWidth', 1.2, ...
         'DisplayName', '混频后');
    grid on;
    xlabel('频率 (MHz)');
    ylabel('幅度');
    title(sprintf('信道 %d', m));
    xlim([-fs/(2*K) fs/(2*K)]/1e6);  % 限制频率范围
    legend('Location','best');
    hold off;
end

% 3. 相位变化分析（以第1、5、9、13条信道为例）
figure('Name','混频相位变化分析','Position',[400 400 1000 800]);
sgtitle('混频前后的相位变化分析');

selected_channels = [1,5,9,13];  % 选择代表性信道
for i = 1:4
    m = selected_channels(i);
    subplot(4,1,i);
    
    % 计算相位
    phase_original = angle(x_branch{m});
    phase_mixed = angle(x_mixed{m});
    phase_diff = phase_mixed - phase_original;  % 相位差
    
    plot(0:L-1, phase_original, 'b--', 'LineWidth', 1, ...
         'DisplayName', '混频前相位');
    hold on;
    plot(0:L-1, phase_mixed, 'r-', 'LineWidth', 1.2, ...
         'DisplayName', '混频后相位');
    plot(0:L-1, phase_diff, 'g-.', 'LineWidth', 1, ...
         'DisplayName', '相位差');
    grid on;
    xlabel('样本索引');
    ylabel('相位（弧度）');
    title(sprintf('信道 %d 的相位变化', m));
    ylim([-pi, pi]);  % 相位范围限制在[-π, π]
    legend('Location','best');
    hold off;
end

% 4. 混频因子可视化
figure('Name','混频因子特性','Position',[500 500 800 500]);
m = 0:K-1;
mix_factors = exp(-1j*pi/2 * m);

subplot(2,1,1);
stem(m, real(mix_factors), 'b', 'MarkerFaceColor', 'b');
hold on;
stem(m, imag(mix_factors), 'r', 'MarkerFaceColor', 'r');
grid on;
xlabel('信道索引 m');
ylabel('幅度');
title('混频因子的实部（蓝色）与虚部（红色）');
legend('实部','虚部','Location','best');
hold off;

subplot(2,1,2);
stem(m, angle(mix_factors), 'g', 'MarkerFaceColor', 'g');
grid on;
xlabel('信道索引 m');
ylabel('相位（弧度）');
title('混频因子的相位特性');
ylim([-pi, pi]);
hold off;
%% 7. 带通滤波（多相滤波器卷积）
x_filtered = cell(K, 1);
for m = 0:K-1
    x_m = x_mixed{m+1};
    h_m = h_poly(m+1, :)';      % 提取第m+1行，转置为列向量（M×1）
    x_filt = conv(x_m, h_m, 'same');  % 同长度卷积（输出长度L）
    x_filtered{m+1} = x_filt;
end


%% 7. 带通滤波结果可视化
colors = lines(K);  % 复用颜色方案
fs_branch = fs/K;   % 支路采样率

% 1. 多相滤波器的频率响应
figure('Name','多相滤波器频率响应','Position',[100 100 1000 800]);
sgtitle('各信道多相滤波器的频率响应');

for m = 1:K
    h_m = h_poly(m, :)';  % 获取第m个滤波器
    [H, f] = freqz(h_m, 1, 1024, fs_branch);  % 计算频率响应
    
    subplot(4,4,m);
    plot(f/1e6, 20*log10(abs(H)), 'Color', colors(m,:), 'LineWidth', 1.2);
    grid on;
    xlabel('频率 (MHz)');
    ylabel('幅度 (dB)');
    title(sprintf('信道 %d 滤波器', m));
    xlim([-fs_branch/2 fs_branch/2]/1e6);
    ylim([-60 5]);
    hold off;
end

% 2. 滤波前后时域信号对比（实部）
figure('Name','滤波前后时域对比','Position',[200 200 1000 800]);
sgtitle('滤波前后的时域信号（实部对比）');

% 前8条信道
for i = 1:8
    subplot(4,2,i);
    % 滤波前信号
    plot(0:L-1, real(x_mixed{i}), 'b--', 'LineWidth', 1, ...
         'DisplayName', '滤波前');
    hold on;
    % 滤波后信号
    plot(0:L-1, real(x_filtered{i}), 'Color', colors(i,:), 'LineWidth', 1.2, ...
         'DisplayName', '滤波后');
    grid on;
    xlabel('样本索引');
    ylabel('信号幅度（实部）');
    title(sprintf('信道 %d', i));
    legend('Location','best');
    hold off;
end

% 后8条信道
figure('Name','滤波前后时域对比（续）','Position',[300 300 1000 800]);
sgtitle('滤波前后的时域信号（实部对比）');

for i = 9:16
    subplot(4,2,i-8);
    % 滤波前信号
    plot(0:L-1, real(x_mixed{i}), 'b--', 'LineWidth', 1, ...
         'DisplayName', '滤波前');
    hold on;
    % 滤波后信号
    plot(0:L-1, real(x_filtered{i}), 'Color', colors(i,:), 'LineWidth', 1.2, ...
         'DisplayName', '滤波后');
    grid on;
    xlabel('样本索引');
    ylabel('信号幅度（实部）');
    title(sprintf('信道 %d', i));
    legend('Location','best');
    hold off;
end

% 3. 滤波前后频谱对比
figure('Name','滤波前后频谱对比','Position',[400 400 1200 800]);
sgtitle('滤波前后的信号频谱对比');

for m = 1:K
    % 计算滤波前后的频谱
    X_fft_mixed = fftshift(fft(x_mixed{m}));
    X_fft_filtered = fftshift(fft(x_filtered{m}));
    f_axis = (-L/2 : L/2 - 1) * (fs_branch/L);  % 支路频率轴（Hz）
    
    subplot(4,4,m);
    % 滤波前频谱
    plot(f_axis/1e6, 20*log10(abs(X_fft_mixed)+eps), 'b--', 'LineWidth', 1, ...
         'DisplayName', '滤波前');
    hold on;
    % 滤波后频谱
    plot(f_axis/1e6, 20*log10(abs(X_fft_filtered)+eps), 'Color', colors(m,:), 'LineWidth', 1.2, ...
         'DisplayName', '滤波后');
    % 叠加滤波器频率响应作为参考
    h_m = h_poly(m, :)';
    [H, f] = freqz(h_m, 1, L, fs_branch);
    plot(f/1e6, 20*log10(abs(H))+30, 'k-.', 'LineWidth', 1, ...  % +30dB偏移以便观察
         'DisplayName', '滤波器响应');
    grid on;
    xlabel('频率 (MHz)');
    ylabel('幅度 (dB)');
    title(sprintf('信道 %d', m));
    xlim([-fs_branch/2 fs_branch/2]/1e6);
    ylim([-60 60]);
    legend('Location','best');
    hold off;
end

% 4. 信号能量变化分析
figure('Name','滤波前后能量对比','Position',[500 500 800 500]);
energy_mixed = zeros(1, K);
energy_filtered = zeros(1, K);

for m = 1:K
    energy_mixed(m) = sum(abs(x_mixed{m}).^2);
    energy_filtered(m) = sum(abs(x_filtered{m}).^2);
end

bar([energy_mixed; energy_filtered]','FaceColor',{'b','r'});
grid on;
xlabel('信道索引');
ylabel('信号能量');
title('滤波前后各信道的信号能量对比');
legend('滤波前','滤波后','Location','best');
xticks(1:K);

% 5. 代表性信道的详细分析（以信道1、5、9、13为例）
figure('Name','滤波效果详细分析','Position',[600 600 1000 800]);
sgtitle('滤波效果详细分析（实部、虚部与相位）');

selected_channels = [1,5,9,13];
for i = 1:4
    m = selected_channels(i);
    subplot(4,1,i);
    
    % 实部对比
    plot(0:L-1, real(x_mixed{m}), 'b--', 'LineWidth', 1, ...
         'DisplayName', '滤波前（实部）');
    hold on;
    plot(0:L-1, real(x_filtered{m}), 'r-', 'LineWidth', 1.2, ...
         'DisplayName', '滤波后（实部）');
    
    % 虚部对比（偏移显示以便区分）
    plot(0:L-1, imag(x_mixed{m})-2, 'g--', 'LineWidth', 1, ...
         'DisplayName', '滤波前（虚部）');
    plot(0:L-1, imag(x_filtered{m})-2, 'm-', 'LineWidth', 1.2, ...
         'DisplayName', '滤波后（虚部）');
    
    grid on;
    xlabel('样本索引');
    ylabel('信号幅度（实部/虚部-2）');
    title(sprintf('信道 %d 的滤波效果', m));
    legend('Location','best');
    ylim([-3 3]);
    hold off;
end

%% 8. 加权（符号+复指数）
x_weighted = cell(K, 1);
for m = 0:K-1
    x_m = x_filtered{m+1};
    k = 0:L-1;
    sign_factor = (-1).^k;          % (-1)^k 符号因子
    exp_factor = exp(1j*pi/2 * m * k);  % 复指数因子
    weight_factor = sign_factor .* exp_factor;
    x_weighted{m+1} = x_m .* weight_factor;
end

%% 8. 加权结果可视化
colors = lines(K);  % 复用颜色方案
fs_branch = fs/K;   % 支路采样率

% 1. 加权前后时域信号对比（实部）
figure('Name','加权前后时域信号对比','Position',[100 100 1000 800]);
sgtitle('加权前后的时域信号（实部对比）');

% 前8条信道对比
for i = 1:8
    subplot(4,2,i);
    % 加权前信号（滤波后）
    plot(0:L-1, real(x_filtered{i}), 'b--', 'LineWidth', 1, ...
         'DisplayName', '加权前（滤波后）');
    hold on;
    % 加权后信号
    plot(0:L-1, real(x_weighted{i}), 'Color', colors(i,:), 'LineWidth', 1.2, ...
         'DisplayName', '加权后');
    grid on;
    xlabel('样本索引');
    ylabel('信号幅度（实部）');
    title(sprintf('信道 %d', i));
    legend('Location','best');
    hold off;
end

% 后8条信道对比
figure('Name','加权前后时域信号对比（续）','Position',[200 200 1000 800]);
sgtitle('加权前后的时域信号（实部对比）');

for i = 9:16
    subplot(4,2,i-8);
    % 加权前信号（滤波后）
    plot(0:L-1, real(x_filtered{i}), 'b--', 'LineWidth', 1, ...
         'DisplayName', '加权前（滤波后）');
    hold on;
    % 加权后信号
    plot(0:L-1, real(x_weighted{i}), 'Color', colors(i,:), 'LineWidth', 1.2, ...
         'DisplayName', '加权后');
    grid on;
    xlabel('样本索引');
    ylabel('信号幅度（实部）');
    title(sprintf('信道 %d', i));
    legend('Location','best');
    hold off;
end

% 2. 加权前后频谱对比
figure('Name','加权前后频谱对比','Position',[300 300 1200 800]);
sgtitle('加权前后的信号频谱对比');

for m = 1:K
    % 计算加权前后的频谱
    X_fft_filtered = fftshift(fft(x_filtered{m}));
    X_fft_weighted = fftshift(fft(x_weighted{m}));
    f_axis = (-L/2 : L/2 - 1) * (fs_branch/L);  % 支路频率轴（Hz）
    
    subplot(4,4,m);
    % 加权前频谱
    plot(f_axis/1e6, 20*log10(abs(X_fft_filtered)+eps), 'b--', 'LineWidth', 1, ...
         'DisplayName', '加权前');
    hold on;
    % 加权后频谱
    plot(f_axis/1e6, 20*log10(abs(X_fft_weighted)+eps), 'Color', colors(m,:), 'LineWidth', 1.2, ...
         'DisplayName', '加权后');
    grid on;
    xlabel('频率 (MHz)');
    ylabel('幅度 (dB)');
    title(sprintf('信道 %d', m));
    xlim([-fs_branch/2 fs_branch/2]/1e6);
    ylim([-60 40]);
    legend('Location','best');
    hold off;
end

% 3. 加权因子特性分析
figure('Name','加权因子特性分析','Position',[400 400 1000 800]);
sgtitle('加权因子（符号因子×复指数因子）特性分析');

k = 0:L-1;  % 支路内时间索引
selected_m = [0,3,7,11];  % 选择4个代表性信道（0-based）

for i = 1:4
    m = selected_m(i);
    % 计算加权因子
    sign_factor = (-1).^k;          % 符号因子
    exp_factor = exp(1j*pi/2 * m * k);  % 复指数因子
    weight_factor = sign_factor .* exp_factor;  % 总加权因子
    
    subplot(4,1,i);
    % 实部
    plot(k, real(weight_factor), 'b-', 'LineWidth', 1.2, ...
         'DisplayName', '实部');
    hold on;
    % 虚部
    plot(k, imag(weight_factor), 'r--', 'LineWidth', 1.2, ...
         'DisplayName', '虚部');
    % 幅度（绝对值）
    plot(k, abs(weight_factor), 'k:', 'LineWidth', 1.5, ...
         'DisplayName', '幅度');
    grid on;
    xlabel('样本索引 k');
    ylabel('值');
    title(sprintf('加权因子（信道 m=%d）', m));
    legend('Location','best');
    ylim([-1.2 1.2]);
    hold off;
end

% 4. 相位变化详细分析（加权前后）
figure('Name','加权前后相位变化','Position',[500 500 1000 800]);
sgtitle('加权前后的信号相位变化');

selected_channels = [1,5,9,13];  % 选择代表性信道（1-based）
for i = 1:4
    m = selected_channels(i);
    subplot(4,1,i);
    
    % 加权前相位（滤波后）
    phase_filtered = angle(x_filtered{m});
    % 加权后相位
    phase_weighted = angle(x_weighted{m});
    % 相位差（理论应为加权因子的相位）
    phase_diff = phase_weighted - phase_filtered;
    % 理论相位（加权因子的相位）
    k = 0:L-1;
    theoretical_phase = angle( (-1).^k .* exp(1j*pi/2*(m-1)*k) );  % m-1是0-based
    
    plot(k, phase_filtered, 'b--', 'LineWidth', 1, ...
         'DisplayName', '加权前相位');
    hold on;
    plot(k, phase_weighted, 'r-', 'LineWidth', 1.2, ...
         'DisplayName', '加权后相位');
    plot(k, phase_diff, 'g-.', 'LineWidth', 1, ...
         'DisplayName', '实际相位差');
    plot(k, theoretical_phase, 'k:', 'LineWidth', 1.5, ...
         'DisplayName', '理论相位差');
    grid on;
    xlabel('样本索引 k');
    ylabel('相位（弧度）');
    title(sprintf('信道 %d 的相位变化', m));
    ylim([-pi pi]);
    legend('Location','best');
    hold off;
end

%% 9. 分组执行4点IDFT（核心修正：行索引范围）
idft_result = zeros(K, L);  % 16行×200列，存储16信道结果
for group = 0:3              % 4组：0~3（每组4信道）
    start_m = group * 4;     % 每组起始信道（0,4,8,12，0-based）
    end_m = start_m + 3;     % 每组结束信道（3,7,11,15，0-based）
    
    % 构建组矩阵：4行×L列（每组4信道）
    group_matrix = zeros(4, L);
    for m = start_m:end_m
        row_idx = m - start_m + 1;  % 组内行索引（1~4，1-based）
        group_matrix(row_idx, :) = x_weighted{m+1};  % 提取第m+1条支路
    end
    
    % 对每组的每一列执行4点IDFT
    for col = 1:L
        col_data = group_matrix(:, col);  % 4×1列数据（1-based）
        idft_col = ifft(col_data, 4);     % 4点IDFT，输出4×1
        % 修正行索引：覆盖start_m+1 ~ end_m+1（共4行，1-based）
        idft_result(start_m+1:end_m+1, col) = idft_col;
    end
end

%% 10. 结果分析（绘制第14信道频谱）
channel_idx = 0;            % 选择信道8（0-based，对应第9行）
channel_signal = idft_result(channel_idx+1, :);  % 1-based索引
N_fft = 1024;               % FFT长度
fft_spec = abs(fft(channel_signal, N_fft)) / N_fft;  % 归一化频谱
freq_axis = (0:N_fft/2-1) * (fs/N_fft) / 1e6;        % 频率轴（MHz）

figure;
plot(freq_axis, fft_spec(1:N_fft/2));
xlabel('频率 (MHz)');
ylabel('归一化幅度');
title(['信道', num2str(channel_idx), '的频谱（采样率', num2str(fs/1e9), 'GHz）']);
grid on;