%% 基于FPGA的数字信道化接收机仿真（2.4.4节）
% 功能：生成正弦信号+线性调频信号，通过多相滤波+IFFT实现16路信道化，输出时域I分量及频域频谱（没有加移频和线性补偿）
% 参考参数：
% - ADC采样频率：1.25GHz
% - 信号1：正弦波，载频1005MHz，功率2dBm
% - 信号2：线性调频（LFM），中心频率830MHz，调频带宽60MHz，功率2dBm
% - 信道化：16路均匀划分，子信道带宽39.0625MHz（1.25e9/32）
% - 原型滤波器：192阶等波纹低通，通带19.0625MHz，阻带39.0625MHz

clear; clc; close all;

%% 1. 基本参数设置
fs = 1.25e9;               % ADC采样频率（1.25GHz）
Tp = 10e-6;                % 脉冲宽度（10μs，确保信号完整采样）
K = 32;                    % 多相滤波器组总路数（后16路与前16路共轭，取前16路有效）
M = 16;                    % 抽取倍数（对应1:16串并转换）
F = K / M;                 % 重叠因子（F=2，满足无混叠条件）
N_prototype = 191;         % 原型低通滤波器阶数
R = 50;                    % 负载电阻（50Ω，用于功率-幅度转换）
P_dBm = 2;                 % 信号功率（2dBm）

% 计算时间向量（确保覆盖2个脉冲周期，避免截断）
t = 0 : 1/fs : 2*Tp - 1/fs;
N = length(t);             % 总采样点数


%% 2. 信号生成（正弦脉冲 + LFM脉冲）
% 2.1 功率转电压幅度（dBm → 电压有效值）
P_mW = 10^(P_dBm / 10);    % dBm转mW
P_W = P_mW / 1000;         % mW转W
A = sqrt(2 * P_W * R);     % 正弦信号电压有效值（考虑50Ω负载）

% 2.2 信号1：正弦脉冲（载频1005MHz）
f1 = 1005e6;               % 信号1载频
signal1 = A * sin(2 * pi * f1 * t) .* rectpuls(t - Tp/2, Tp);  % 带脉冲包络

% 2.3 信号2：线性调频（LFM）脉冲（中心频率830MHz，带宽60MHz）
f2 = 830e6;                % LFM中心频率
B = 60e6;                  % 调频带宽
k = B / Tp;                % 调频斜率
% LFM复数表达式：实部为信号，虚部后续信道化生成
lfm = A * cos(2 * pi * f2 * t + pi * k * t.^2) .* rectpuls(t - Tp/2, Tp);
signal2 = lfm;             % 取实部作为LFM信号


%% 3. 无噪声信号（已移除加性高斯白噪声）
x = signal1 + signal2;  % 最终ADC采样输入信号（无噪声）


%% 4. 设计原型低通滤波器（等波纹法，192阶）
% 滤波器参数（归一化频率：f/(fs/2)）
f_pass = 19.0625e6;        % 通带截止频率（文档2.4.2节）
f_stop = 39.0625e6;        % 阻带起始频率（文档2.4.2节）
Wp = f_pass / (fs / 2);    % 归一化通带
Ws = f_stop / (fs / 2);    % 归一化阻带
Rp = 0.1;                  % 通带波纹（0.1dB）
Rs = 60;                   % 阻带衰减（60dB）

% 等波纹法设计FIR滤波器（firpm函数）
h_prototype = firpm(N_prototype, [0, Wp, Ws, 1], [1, 1, 0, 0], [1, Rs/(Rp)]); %firpm函数的阶数定义为系数个数 - 1
% 验证滤波器幅频响应（可选）
figure; freqz(h_prototype, 1, 1024, fs/1e6); title('原型低通滤波器幅频响应');


%% 5. 多相滤波器组分解（K=32路，2倍插0）
% 5.1 多相分量分解：将原型滤波器分为K路，每路长度ceil(N_prototype/K)
h_poly = reshape(h_prototype, K, [])';  % K=32列，每列是一个多相分量（6阶）
% 5.2 2倍插0（文档2.4.1节：F=2时多相分量2倍插0）
h_poly_up = upsample(h_poly, 2);  % 每个多相分量插0后变为12阶


%% 6. 数字信道化处理（串并转换 → 多相滤波 → IFFT）
% 6.1 M=16抽取对应文档3.3.1节）
x_parallel = x(1:M:end);    % 16倍抽取，得到并行数据（长度N/M）
L = length(x_parallel);     % 并行数据长度

% 6.2 多相滤波：32路并行滤波（每路对应一个多相分量）
x_filtered = zeros(K, L);   % 存储32路滤波后数据
for k = 1:K
    % 第k路多相滤波（卷积实现）
    x_filtered(k, :) = conv(x_parallel, h_poly_up(:, k), 'same');
end

% 6.3 32点IFFT运算（文档2.4.1节：高效信道化核心）
x_ifft = fftshift(ifft(x_filtered, K, 1), 1);  % 沿第1维（信道维度）IFFT，fftshift对齐频率

% 6.4 提取有效信道（前16路，后16路共轭对称，文档2.4.1节）
x_channelized = x_ifft(1:K/2, :);  % 16路信道化输出（复数IQ信号）
I_components = real(x_channelized); % 提取I分量（实部）


%% 7. 结果可视化（时域I分量 + 频域频谱）
% 7.1 时域I分量图（16路信道，2行8列布局）
figure('Position', [100, 100, 1200, 600]);
for ch = 1:16
    subplot(2, 8, ch);
    plot(t(1:M:end)*1e6, I_components(ch, :));  % 时间轴单位：μs
    xlabel('时间（μs）');
    ylabel('I分量幅度');
    title(['信道', num2str(ch)]);
    grid on;
    xlim([0, 2*Tp*1e6]);  % 显示2个脉冲周期
end
sgtitle('16路信道化输出I分量时域波形');

% 7.2 频域频谱图（16路信道，2行8列布局）
% 计算每个信道的频率轴（对应1250~625MHz，文档2.4.2节信道划分）
f_channel = (fs/2 : -fs/K : fs/2 - fs/K*(K/2-1)) / 1e6;  % 16路信道中心频率（MHz）
figure('Position', [100, 100, 1200, 600]);
for ch = 1:16
    subplot(2, 8, ch);
    % 计算当前信道I分量的频谱（FFT）
    fft_I = fftshift(fft(I_components(ch, :)));
    f = linspace(-fs/(2*M), fs/(2*M), L) / 1e6;  % 子信道频率轴（MHz）
    plot(f + f_channel(ch), 20*log10(abs(fft_I)/max(abs(fft_I))));  % 归一化频谱
    xlabel('频率（MHz）');
    ylabel('幅度（dB）');
    title(['信道', num2str(ch), '（中心频率：', num2str(f_channel(ch)), 'MHz）']);
    grid on;
    ylim([-60, 0]);  % 显示0~-60dB范围
end
sgtitle('16路信道化输出I分量频域频谱');


%% 8. 仿真结果说明
fprintf('仿真完成！关键参数：\n');
fprintf('1. 采样频率：%.2f GHz\n', fs/1e9);
fprintf('2. 子信道数量：%d 路\n', K/2);
fprintf('3. 子信道带宽：%.2f MHz\n', fs/K/1e6);
fprintf('4. 原型滤波器阶数：%d 阶\n', N_prototype);