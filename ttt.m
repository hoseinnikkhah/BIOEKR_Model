% Read CSV file
data = readmatrix('1.csv');

wave_number = data(2:1977, 1);
alginate = data(2:1977, 2);
alg_ac = data(2:1977, 4);
AC = data(2:1977, 6);
a1 = data(2:1977,8);
a2 = data(2:1977,9);
a3 = data(2:1977,10);
a4 = data(2:1977,11);
a5 = data(2:1977,12);
a6 = data(2:1977,13);
a7 = data(2:1977,14);
a8 = data(2:1977,15);

figure(1);
hold on;

plot(wave_number,alginate);
plot(wave_number,alg_ac);
plot(wave_number,AC);
plot(wave_number,a1);
plot(wave_number,a2);
plot(wave_number,a3);
plot(wave_number,a4);
plot(wave_number,a5);
plot(wave_number,a6);
plot(wave_number,a7);
plot(wave_number,a8);



%plot(time, Bl, '-o', 'LineWidth', 1, 'MarkerSize', 4, 'Color', 'b','DisplayName', 'test 3');