experimentID = [21, 22, 23, 24, 25, 26, 27, 28, 29, 30];
accel_arr = ["ACC\_X", "ACC\_Y", "ACC\_Z"];

fs = 50;


%File Readers
user10 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/ATD2020/RawData/acc_exp21_user10.txt', 21);
user11_1 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/ATD2020/RawData/acc_exp22_user11.txt', 22);
user11_2 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/ATD2020/RawData/acc_exp23_user11.txt', 23);
user12_1 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/ATD2020/RawData/acc_exp24_user12.txt', 24);
user12_2 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/ATD2020/RawData/acc_exp25_user12.txt', 25);
user13_1 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/ATD2020/RawData/acc_exp26_user13.txt', 26);
user13_2 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/ATD2020/RawData/acc_exp27_user13.txt', 27);
user14_1 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/ATD2020/RawData/acc_exp28_user14.txt', 28);
user14_2 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/ATD2020/RawData/acc_exp29_user14.txt', 29);
user15 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/ATD2020/RawData/acc_exp30_user15.txt', 30);


%Activity Index Reader
activityfile = fileread('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/ATD2020/RawData/activity_labels.txt');
activitycells = textscan(activityfile, '%f %s');

%Labels Reader
labelsfile = fileread('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/ATD2020/RawData/labels.txt');
labelsdata = cell2mat(textscan(labelsfile, '%f%f%f%f%f'));
lengthofarr = length(labelsdata);
i = 1;

while i <= lengthofarr
    if ismember(labelsdata(i, 1),experimentID) == 0
        labelsdata(i,:) = [];
        lengthofarr = lengthofarr - 1;
    else
        i = i + 1;
    end
end    

%raw_signal(user10, activitycells, labelsdata, fs, 1);
%raw_signal(user11_1, activitycells, labelsdata, fs, 2);
%raw_signal(user11_2, activitycells, labelsdata, fs, 3);
%raw_signal(user12_1, activitycells, labelsdata, fs, 4);
%raw_signal(user12_2, activitycells, labelsdata, fs, 5);
%raw_signal(user13_1, activitycells, labelsdata, fs, 6);
%raw_signal(user13_2, activitycells, labelsdata, fs, 7);
%raw_signal(user14_1, activitycells, labelsdata, fs, 8);
%raw_signal(user14_2, activitycells, labelsdata, fs, 9);
%raw_signal(user15, activitycells, labelsdata, fs, 10);

%DFT_Hamming(user10, fs, labelsdata, activitycells, 11, accel_arr);
%DFT_Hamming(user11_1, fs, labelsdata, activitycells, 12, accel_arr);
%DFT_Hamming(user11_2, fs, labelsdata, activitycells, 13, accel_arr);
%DFT_Hamming(user12_1, fs, labelsdata, activitycells, 14, accel_arr);
%DFT_Hamming(user12_2, fs, labelsdata, activitycells, 15, accel_arr);
%DFT_Hamming(user13_1, fs, labelsdata, activitycells, 16, accel_arr);
%DFT_Hamming(user13_2, fs, labelsdata, activitycells, 17, accel_arr);
%DFT_Hamming(user14_1, fs, labelsdata, activitycells, 18, accel_arr);
%DFT_Hamming(user14_2, fs, labelsdata, activitycells, 19, accel_arr);
%DFT_Hamming(user15, fs, labelsdata, activitycells, 20, accel_arr);

%DFT_rectangular(user10, fs, labelsdata, activitycells, 11, accel_arr);
%DFT_rectangular(user11_1, fs, labelsdata, activitycells, 12, accel_arr);
%DFT_rectangular(user11_2, fs, labelsdata, activitycells, 13, accel_arr);
%DFT_rectangular(user12_1, fs, labelsdata, activitycells, 14, accel_arr);
%DFT_rectangular(user12_2, fs, labelsdata, activitycells, 15, accel_arr);
%DFT_rectangular(user13_1, fs, labelsdata, activitycells, 16, accel_arr);
%DFT_rectangular(user13_2, fs, labelsdata, activitycells, 70, accel_arr);
%DFT_rectangular(user14_1, fs, labelsdata, activitycells, 18, accel_arr);
%DFT_rectangular(user14_2, fs, labelsdata, activitycells, 19, accel_arr);
%DFT_rectangular(user15, fs, labelsdata, activitycells, 20, accel_arr);

w_total = zeros(10,1);
wup_total = zeros(10,1);
wd_total = zeros(10,1);
w_steps = [];
wup_steps = [];
wd_steps = [];

%WALKING STEPS
[w_total(1,1), w_steps] = store_steps(user10(:,3), labelsdata, experimentID(1), 1, w_steps, w_total(1,1), fs);
[w_total(2,1), w_steps] = store_steps(user11_1(:,3), labelsdata, experimentID(2), 1, w_steps, w_total(2,1), fs);
[w_total(3,1), w_steps] = store_steps(user11_2(:,3), labelsdata, experimentID(3), 1, w_steps, w_total(3,1), fs);
[w_total(4,1), w_steps] = store_steps(user12_1(:,3), labelsdata, experimentID(4), 1, w_steps, w_total(4,1), fs);
[w_total(5,1), w_steps] = store_steps(user12_2(:,3), labelsdata, experimentID(5), 1, w_steps, w_total(5,1), fs);
[w_total(6,1), w_steps] = store_steps(user13_1(:,3), labelsdata, experimentID(6), 1, w_steps, w_total(6,1), fs);
[w_total(7,1), w_steps] = store_steps(user13_2(:,3), labelsdata, experimentID(7), 1, w_steps, w_total(7,1), fs);
[w_total(8,1), w_steps] = store_steps(user14_1(:,3), labelsdata, experimentID(8), 1, w_steps, w_total(8,1), fs);
[w_total(9,1), w_steps] = store_steps(user14_2(:,3), labelsdata, experimentID(9), 1, w_steps, w_total(9,1), fs);
[w_total(10,1), w_steps] = store_steps(user15(:,3), labelsdata, experimentID(10), 1, w_steps, w_total(10,1), fs);

%WALKING_UP STEPS
[wup_total(1,1), wup_steps] = store_steps(user10(:,3), labelsdata, experimentID(1), 2, wup_steps, wup_total(1,1), fs);
[wup_total(2,1), wup_steps] = store_steps(user11_1(:,3), labelsdata, experimentID(2), 2, wup_steps, wup_total(2,1), fs);
[wup_total(3,1), wup_steps] = store_steps(user11_2(:,3), labelsdata, experimentID(3), 2, wup_steps, wup_total(3,1), fs);
[wup_total(4,1), wup_steps] = store_steps(user12_1(:,3), labelsdata, experimentID(4), 2, wup_steps, wup_total(4,1), fs);
[wup_total(5,1), wup_steps] = store_steps(user12_2(:,3), labelsdata, experimentID(5), 2, wup_steps, wup_total(5,1), fs);
[wup_total(6,1), wup_steps] = store_steps(user13_1(:,3), labelsdata, experimentID(6), 2, wup_steps, wup_total(6,1), fs);
[wup_total(7,1), wup_steps] = store_steps(user13_2(:,3), labelsdata, experimentID(7), 2, wup_steps, wup_total(7,1), fs);
[wup_total(8,1), wup_steps] = store_steps(user14_1(:,3), labelsdata, experimentID(8), 2, wup_steps, wup_total(8,1), fs);
[wup_total(9,1), wup_steps] = store_steps(user14_2(:,3), labelsdata, experimentID(9), 2, wup_steps, wup_total(9,1), fs);
[wup_total(10,1), wup_steps] = store_steps(user15(:,3), labelsdata, experimentID(10), 2, wup_steps, wup_total(10,1), fs);

%WALKING_DOWN STEPS
[wd_total(1,1), wd_steps] = store_steps(user10(:,3), labelsdata, experimentID(1), 3, wd_steps, wd_total(1,1), fs);
[wd_total(2,1), wd_steps] = store_steps(user11_1(:,3), labelsdata, experimentID(2), 3, wd_steps, wd_total(2,1), fs);
[wd_total(3,1), wd_steps] = store_steps(user11_2(:,3), labelsdata, experimentID(3), 3, wd_steps, wd_total(3,1), fs);
[wd_total(4,1), wd_steps] = store_steps(user12_1(:,3), labelsdata, experimentID(4), 3, wd_steps, wd_total(4,1), fs);
[wd_total(5,1), wd_steps] = store_steps(user12_2(:,3), labelsdata, experimentID(5), 3, wd_steps, wd_total(5,1), fs);
[wd_total(6,1), wd_steps] = store_steps(user13_1(:,3), labelsdata, experimentID(6), 3, wd_steps, wd_total(6,1), fs);
[wd_total(7,1), wd_steps] = store_steps(user13_2(:,3), labelsdata, experimentID(7), 3, wd_steps, wd_total(7,1), fs);
[wd_total(8,1), wd_steps] = store_steps(user14_1(:,3), labelsdata, experimentID(8), 3, wd_steps, wd_total(8,1), fs);
[wd_total(9,1), wd_steps] = store_steps(user14_2(:,3), labelsdata, experimentID(9), 3, wd_steps, wd_total(9,1), fs);
[wd_total(10,1), wd_steps] = store_steps(user15(:,3), labelsdata, experimentID(10), 3, wd_steps, wd_total(10,1), fs);


valtotal_w = 0;
valtotal_wup = 0;
valtotal_wd = 0;

for i = 1:length(w_total)
    valtotal_w = valtotal_w + w_total(i);
    valtotal_wup = valtotal_wup + wup_total(i);
    valtotal_wd = valtotal_wd + wd_total(i);
end

fprintf("[WALKING] Total por sinal:\n")
disp(w_total)
fprintf("[WALKING_UP] Total por sinal:\n")
disp(wup_total)
fprintf("[WALKING_DOWN] Total por sinal:\n")
disp(wd_total)
fprintf("\n[WALKING] Valor total de passos: %d\n", valtotal_w)
fprintf("[WALKING_UP] Valor total de passos: %d\n", valtotal_wup)
fprintf("[WALKING_DOWN] Valor total de passos: %d\n\n", valtotal_wd)



mean_wsteps = mean(w_total);
mean_wupsteps = mean(wup_total);
mean_wdsteps = mean(wd_total);

fprintf("[WALKING] Media de passos: %d\n", mean_wsteps)
fprintf("[WALKING_UP] Media de passos: %d\n", mean_wupsteps)
fprintf("[WALKING_DOWN] Media de passos: %d\n\n", mean_wdsteps)


dp_wsteps = w_steps-mean_wsteps;
dp_wsteps = (dp_wsteps.^2)/length(dp_wsteps);
dp_wupsteps = wup_steps-mean_wupsteps;
dp_wupsteps = (dp_wupsteps.^2)/length(dp_wupsteps);
dp_wdsteps = wd_steps-mean_wdsteps;;
dp_wdsteps = (dp_wdsteps.^2)/length(dp_wdsteps);

totaldp_wsteps = 0;
totaldp_wupsteps = 0;
totaldp_wdsteps = 0;


for i=1:length(dp_wsteps)
    totaldp_wsteps = totaldp_wsteps + dp_wsteps(i);
end

for i=1:length(dp_wupsteps)
    totaldp_wupsteps = totaldp_wupsteps + dp_wupsteps(i);
end

for i=1:length(dp_wdsteps)
    totaldp_wdsteps = totaldp_wdsteps + dp_wdsteps(i);
end

fprintf("[WALKING] Desvio padrao: %d\n", totaldp_wsteps)
fprintf("[WALKING_UP] Desvio padrao: %d\n", totaldp_wupsteps)
fprintf("[WALKING_DOWN] Desvio padrao: %d\n\n", totaldp_wdsteps)


%Caracteristicas espectrais
hold on
spec_features(user15, activitycells, labelsdata, 1, [1 0 0]) %WALKING U10
spec_features(user15, activitycells, labelsdata, 2, [0 0 1]) %WALKING UP U10
spec_features(user15, activitycells, labelsdata, 3, [1 0 1]) %WALKING DOWN U10
spec_features(user15, activitycells, labelsdata, 4, [0 1 0]) %SITTING U10
spec_features(user15, activitycells, labelsdata, 5, [1 0.5 1]) %STANDING U10
spec_features(user15, activitycells, labelsdata, 6, [0.3 0.6 0.9]) %LAYING U10
spec_features(user15, activitycells, labelsdata, 7, [1 1 0]) %STAND TO SIT U10
spec_features(user15, activitycells, labelsdata, 8, [0 1 1]) %SIT TO STAND U10
spec_features(user15, activitycells, labelsdata, 9, [0.5 1 0.5]) %SIT TO LIE U10
spec_features(user15, activitycells, labelsdata, 10, [0.3 0.3 0.3]) %LIE TO SIT U10
spec_features(user15, activitycells, labelsdata, 11, [0.6 0.6 0.6]) %STAND TO LIE U10
spec_features(user15, activitycells, labelsdata, 12, [0.9 0.9 0.3]) %LIE TO STAND U10
hold off

%STFT
Zaxis = {user10(:,3)};
STFT(Zaxis);
Zaxis = {user11_1(:,3)};
STFT(Zaxis);
Zaxis = {user11_2(:,3)};
STFT(Zaxis);
Zaxis = {user12_1(:,3)};
STFT(Zaxis);
Zaxis = {user12_2(:,3)};
STFT(Zaxis);
Zaxis = {user13_1(:,3)};
STFT(Zaxis);
Zaxis = {user13_2(:,3)};
STFT(Zaxis);
Zaxis = {user14_1(:,3)};
STFT(Zaxis);
Zaxis = {user14_2(:,3)};
STFT(Zaxis);
Zaxis = {user15(:,3)};
STFT(Zaxis);

function raw_signal(data, activity, labels, fs, fig_num)
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);
    id = data(1, 4);
    timevec = [0: length(data)-1]./fs;
    timevec = timevec./60; %conversao para minutos
    figure(fig_num);
    tiledlayout(6,1);
    ax1 = nexttile;
    plot(ax1, timevec, x);
    activityplot(x, ax1, timevec, activity, labels, id);
    title(ax1, 'Accelerometer\_X');
    ylabel(ax1, 'acc\_x');
    xlabel(ax1, 'Time(min)');
    ax2 = nexttile;
    plot(ax2, timevec, y);
    activityplot(y, ax2, timevec, activity, labels, id);
    title(ax2, 'Accelerometer\_Y');
    ylabel(ax2, 'acc\_y');
    xlabel(ax2, 'Time(min)');
    ax3 = nexttile;
    plot(ax3, timevec, z);
    activityplot(z, ax3, timevec, activity, labels, id);
    title(ax3, 'Accelerometer\_Z');
    ylabel(ax3, 'acc\_z');
    xlabel(ax3, 'Time(min)');
end

function activityplot(data, ax, time, activity, labels, id)
    colorArray = ["1172BD", "D95319", "EDB120", "7E2F8E", "77AC30", "4DBEEE", "A2142F", "48C9B0", "27AE60", "A569BD", "5D6D7E", "D35400"];
    hold(ax, 'on')
    index = 1;
    while index <= length(labels)
        if labels(index, 1) == id
            length4 = labels(index,4);
            length5 = labels(index,5);
            t = linspace(time(labels(index,4)), time(labels(index,5)), length5-length4+1);
            hex = hex2dec(colorArray(labels(index, 3)));
            hex = dec2hex(hex);
            hex2 = strcat('#', hex);
            plot(ax, t, data(length4:length5), 'Color', hex2)
            if labels(index, 3) > 6
                text(((time(length5)-time(length4))/2)+time(length4), max(data), activity{2}{labels(index, 3)})
            else
                if labels(index, 3) == 1 || labels(index,3) == 3
                    text(((time(length5)-time(length4))/2)+time(length4), max(data), activity{2}{labels(index, 3)})
                else
                    text(((time(length5)-time(length4))/2)+time(length4), min(data), activity{2}{labels(index, 3)})
                end
            end
        end
        index = index + 1;
    end
    hold(ax, 'off')
end

function DFT_Hamming(data, fs, labels, activity, fig_num, arraccel)
    colorArray = ["1172BD", "D95319", "EDB120"];
    plotn = 211;
    counter = 1;
    color = 1;
    for accel=1:1:3
        index = 1;
        while index < length(labels)
            if counter == 2
                fig_num = fig_num + 1;
                plotn = 211;
                counter = 1;
            end
            while data(1, 4) ~= labels(index, 1) && index < length(labels)
                index = index + 1;
            end
            length1 = labels(index, 4);
            length2 = labels(index, 5);
            N = length2 - length1 + 1;
            if(mod(N,2)==0)
                f = -fs/2:fs/N:fs/2-fs/N;
            else
                f = -fs/2+fs/(2*N):fs/N:fs/2-fs/(2*N);
            end
            T = fftshift(fft(data([length1:length2],accel)));
            T(abs(T)<0.001)=0;
            m_T = abs(T);
            t_Hamming = data([length1:length2], accel).*hamming(N);
            T_HammW = fftshift(fft(t_Hamming));
            m_HammW = abs(T_HammW);
            %time = linspace(0,(N-1)/fs,N);
            figure(fig_num);
            hex = hex2dec(colorArray(color));
            hex = dec2hex(hex);
            hex2 = strcat('#', hex);
            subplot(plotn);
            plot(f, m_T, 'Color', hex2)
            hex = hex2dec(colorArray(color+1));
            hex = dec2hex(hex);
            hex2 = strcat('#', hex);
            hold on
            plot(f, 1.85*m_HammW, 'Color', hex2)
            str = sprintf('[DFT] da atividade %s com Hamming Window (laranja) no acelerometro %s\n', activity{2}{labels(index,3)}, arraccel(accel));
            title(str)
            hold off
            plotn = plotn + 1;
            hex = hex2dec(colorArray(color+2));
            hex = dec2hex(hex);
            hex2 = strcat('#', hex);
            subplot(plotn);
            plot([1:N], 1.85*t_Hamming, 'Color', hex2)
            str = sprintf('Hamming Window da atividade %s no acelerometro %s\n', activity{2}{labels(index,3)}, arraccel(accel));
            title(str)
            counter = counter + 1;
            index = index + 1;
            plotn = plotn + 1;
            color = 1;
        end
    end
end

function DFT_rectangular(data, fs, labels, activity, fig_num, arraccel)
    colorArray = ["1172BD", "D95319", "EDB120"];
    plotn = 211;
    counter = 1;
    color = 1;
    for accel=1:1:3
        index = 1;
        while index < length(labels)
            if counter == 2
                fig_num = fig_num + 1;
                plotn = 211;
                counter = 1;
            end
            while data(1, 4) ~= labels(index, 1) && index < length(labels)
                index = index + 1;
            end
            length1 = labels(index, 4);
            length2 = labels(index, 5);
            N = length2 - length1 + 1;
            if(mod(N,2)==0)
                f = -fs/2:fs/N:fs/2-fs/N;
            else
                f = -fs/2+fs/(2*N):fs/N:fs/2-fs/(2*N);
            end
            T = fftshift(fft(data([length1:length2],accel)));
            T(abs(T)<0.001)=0;
            m_T = abs(T);
            t_rect = data([length1:length2], accel).*rectwin(N);
            T_rectW = fftshift(fft(t_rect));
            m_rectW = abs(T_rectW);
            %time = linspace(0,(N-1)/fs,N);
            figure(fig_num);
            hex = hex2dec(colorArray(color));
            hex = dec2hex(hex);
            hex2 = strcat('#', hex);
            subplot(plotn);
            plot(f, m_T, 'Color', hex2)
            hex = hex2dec(colorArray(color+1));
            hex = dec2hex(hex);
            hex2 = strcat('#', hex);
            hold on
            plot(f, m_rectW, 'Color', hex2)
            str = sprintf('[DFT] da atividade %s com Rectangular Window (laranja) no acelerometro %s\n', activity{2}{labels(index,3)}, arraccel(accel));
            title(str)
            hold off
            plotn = plotn + 1;
            hex = hex2dec(colorArray(color+2));
            hex = dec2hex(hex);
            hex2 = strcat('#', hex);
            subplot(plotn);
            plot([1:N], t_rect, 'Color', hex2)
            str = sprintf('Rectangular Window da atividade %s no acelerometro %s\n', activity{2}{labels(index,3)}, arraccel(accel));
            title(str)
            counter = counter + 1;
            index = index + 1;
            plotn = plotn + 1;
            color = 1;
        end
    end
end

function [total, min] = step_count(data, arr_min, arr_total, fs)
    bm_window = blackman(length(data));
    T_bm = fftshift(fft(detrend(data).*bm_window));
    m_bm = abs(T_bm);
    [~, peaks_f] = findpeaks(m_bm);
    Ts = 1/length(peaks_f(1));
    min = [arr_min; 60/Ts];
    total = arr_total + (length(data)/fs)/Ts;
end

function [arr_total, arr_min] = store_steps(data, labelsdata, id, activity_id, arr_min, arr_total, fs)
    for index=1:length(labelsdata)
        if labelsdata(index, 3) == activity_id && labelsdata(index, 1) == id
            len_max = labelsdata(index, 5);
            len_min = labelsdata(index, 4);
            [arr_total, arr_min] = step_count(data([len_min:len_max]), arr_min, arr_total, fs);
        end
    end
end

function plot = spec_features(data, activity, labels, id, color)
    m_x = 0;
    m_y = 0;
    m_z = 0;
    for i=1:3
        index = 1;
        while index <= length(labels)
            if labels(index, 3) == id && labels(index, 1) == data(1,4)
                min_len = labels(index, 4);
                max_len = labels(index, 5);
                bm_win = blackman(length(data(min_len:max_len, i)));
                m = abs(fftshift(fft(detrend(data(min_len:max_len,i)).*bm_win)));
                [peaks_m, ~] = findpeaks(m);
                if i == 1
                    m_x = peaks_m;
                elseif i == 2
                    m_y = peaks_m;
                else
                    m_z = peaks_m;
                end
            end
            index = index + 1;
        end
    end
    %figure(fig_num)
    lengths = [length(m_x), length(m_y), length(m_z)];
    max = min(lengths);
    S = repmat([70,50,20],max,1);
    s = S(:);
    plot = scatter3(m_x(1:max), m_y(1:max), m_z(1:max), s(1:max), color, 'filled', 'MarkerEdgeColor', 'k');
    str = ["Scatter of activities"];
    title(str)
    grid on
    view(-30,10)
end

function data = reader(input, id)
    f = fileread(input);
    data = cell2mat(textscan(f, '%f%f%f%f'));
    data(:,4) = [];
    data(1, 4) = id;
end