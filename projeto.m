experimentID = [21, 22, 23, 24, 25, 26, 27, 28, 29, 30];

%File Readers
user10 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/RawData/acc_exp21_user10.txt', 21);
user11_1 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/RawData/acc_exp22_user11.txt', 22);
user11_2 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/RawData/acc_exp23_user11.txt', 23);
user12_1 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/RawData/acc_exp24_user12.txt', 24);
user12_2 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/RawData/acc_exp25_user12.txt', 25);
user13_1 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/RawData/acc_exp26_user13.txt', 26);
user13_2 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/RawData/acc_exp27_user13.txt', 27);
user14_1 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/RawData/acc_exp28_user14.txt', 28);
user14_2 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/RawData/acc_exp29_user14.txt', 29);
user15 = reader('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/RawData/acc_exp30_user15.txt', 30);

%Activity Index Reader
activityfile = fileread('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/RawData/activity_labels.txt');
activitycells = textscan(activityfile, '%f %s');
for i = 1:length(activitycells{2})
    if activitycells{1}(i) == 6
        activitycells{2}{i};
    end
end
%activitydata = cell2mat()

%Labels Reader
labelsfile = fileread('/home/loirito/Documents/Documents/1920/2Sem/ATD/Projeto/RawData/labels.txt');
labelsdata = cell2mat(textscan(labelsfile, '%f%f%f%f%f'));
lengthofarr = length(labelsdata);
i = 1;
while i <= lengthofarr
    if ismember(labelsdata(i, 1),experimentID) == 0
        labelsdata(i,:) = [];
        lengthofarr = lengthofarr - 1;
    else i = i + 1;
    end
end

        


raw_signal(user10, activitycells, labelsdata);
%raw_signal(user11_1, activitycells, labelsdata);
%raw_signal(user11_2);
%raw_signal(user12_1);
%raw_signal(user12_2);
%raw_signal(user13_1);
%raw_signal(user13_2);
%raw_signal(user14_1);
%raw_signal(user14_2);
%raw_signal(user15);

function raw_signal(data, activity, labels)
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);
    id = data(1, 4);
    t = linspace(0, 7, length(x));
    tiledlayout(3,1);
    ax1 = nexttile;
    plot(ax1, t, x);
    activityplot(x, ax1, t, activity, labels, id);
    title(ax1, 'Accelerometer\_X');
    ylabel(ax1, 'acc\_x');
    xlabel(ax1, 'Time(min)');
    ax2 = nexttile;
    plot(ax2, t, y);
    activityplot(y, ax2, t, activity, labels, id);
    title(ax2, 'Accelerometer\_Y');
    ylabel(ax2, 'acc\_y');
    xlabel(ax2, 'Time(min)');
    ax3 = nexttile;
    plot(ax3, t, z);
    activityplot(z, ax3, t, activity, labels, id);
    title(ax3, 'Accelerometer\_Z');
    ylabel(ax3, 'acc\_z');
    xlabel(ax3, 'Time(min)');
end

function activityplot(data, ax, time, activity, labels, id)
    colorArray = ["1172BD", "D95319", "EDB120", "7E2F8E", "77AC30", "4DBEEE", "A2142F", "48C9B0", "27AE60", "A569BD", "5D6D7E", "D35400"];
    hold(ax, 'on')
    index = 1;
    while index < length(labels)
        if labels(index, 1) == id
            length4 = labels(index,4);
            length5 = labels(index,5);
            t = linspace(time(labels(index,4)), time(labels(index,5)), length5-length4+1);
            hex = hex2dec(colorArray(labels(index, 3)));
            hex = dec2hex(hex);
            hex2 = strcat('#', hex);
            plot(ax, t, data(length4:length5), 'Color', hex2)
            if labels(index, 3) > 6
                text(((time(length5)-time(length4))/2)+time(length4), 1.5, activity{2}{labels(index, 3)})
            else
                text(((time(length5)-time(length4))/2)+time(length4), -0.5, activity{2}{labels(index, 3)})
            end
        end
        index = index + 1;
    end
    hold(ax, 'off')
end

function data = reader(input, id)
    f = fileread(input);
    data = cell2mat(textscan(f, '%f%f%f%f'));
    data(:,4) = [];
    data(1, 4) = id;
end

