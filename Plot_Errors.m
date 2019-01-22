
sigma_tau = dlmread('sigma.dat');

NN_error = dlmread('NN_Errors.dat');

Pereb_errors = dlmread('Perebor_Errors.dat');


plot(sigma_tau, NN_error, 'b', sigma_tau, Pereb_errors, 'r'); grid on;

xlabel('СКО измерения времени задержки (мкс)');
ylabel('СКО местоопределения (км)');