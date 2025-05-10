% Sampling parameters
Fs = 1000;              % Sampling frequency (in Hz)
T = 1/Fs;               % Sampling period
L = 1000;               % Length of signal (number of samples)
t = (0:L-1) * T;        % Time vector

% Define the signals
f1 = sin(220 * pi * t);  % 110 Hz signal
f2 = sin(440 * pi * t);  % 220 Hz signal
f3 = sin(880 * pi * t);  % 440 Hz signal
f = f1 + f2 + f3;

% Perform FFT for each signal
F1 = fft(f1);
F2 = fft(f2);
F3 = fft(f3);
F = fft(f);

% Compute the two-sided spectrum, then compute the single-sided spectrum
P2_1 = abs(F1/L);          % Two-sided spectrum for f1
P1_1 = P2_1(1:L/2+1);      % Single-sided spectrum for f1
P1_1(2:end-1) = 2*P1_1(2:end-1);

P2_2 = abs(F2/L);          % Two-sided spectrum for f2
P1_2 = P2_2(1:L/2+1);      % Single-sided spectrum for f2
P1_2(2:end-1) = 2*P1_2(2:end-1);

P2_3 = abs(F3/L);          % Two-sided spectrum for f3
P1_3 = P2_3(1:L/2+1);      % Single-sided spectrum for f3
P1_3(2:end-1) = 2*P1_3(2:end-1);

P2_4 = abs(F/L);          % Two-sided spectrum for f3
P1_4 = P2_4(1:L/2+1);      % Single-sided spectrum for f3
P1_4(2:end-1) = 2*P1_4(2:end-1);

% Define the frequency axis
f_axis = Fs * (0:(L/2)) / L;

% Plot the FFT results
figure;
hold on;

plot(f_axis, P1_1, 'r', 'DisplayName', '110 Hz A2');
plot(f_axis, P1_2, 'g', 'DisplayName', '220 Hz A3');
plot(f_axis, P1_3, 'b', 'DisplayName', '440 Hz A4');
%plot(f_axis, P1_4, 'DisplayName', 'A2, A3, A4')

title('Single-Sided Amplitude Spectrum of f_1(t), f_2(t), f_3(t)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 500]);          % Limit x-axis to 1000 Hz to include all signals
ylim([0 1.2]);
legend('show');

hold off;
