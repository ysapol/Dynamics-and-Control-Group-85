% filename = '';
Data = importdata("MeasuredSignals-1.mat",delimiterIn,headerlinesIn);

for k = [3, 5]
    disp(Data.data(:, k))
    disp(' ')
end