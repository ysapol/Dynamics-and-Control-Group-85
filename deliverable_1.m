f = fopen("MeasuredSignals-1.mat");
data = textscan(f,'%s');
fclose(f);
for i 1:3
variable{i} = (data{i}(2:2:end));
end
disp(variable{i})