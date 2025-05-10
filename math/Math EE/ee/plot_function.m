% Define the x range
t = linspace(0, 0.02, 1000);

% Define the functions
f1 = sin(220 * pi * t);
f2 = sin(440 * pi * t);
f3 = sin(880 * pi * t);
f=f1+f2+f3;
% Plot the functions
figure;
hold on;
%plot(t, f1, 'r', 'LineWidth', 1.5); % Plot f_1(x) in red
%plot(t, f2, 'g', 'LineWidth', 1.5); % Plot f_2(x) in green
%plot(t, f3, 'b', 'LineWidth', 1.5); % Plot f_3(x) in blue
plot(t,f,'LineWidth',1.5);
hold off;

% Set axis limits
ylim([-3.2 3.2]);

% Add labels and title
xlabel('x');
ylabel('f(x)');
title('Plots of f(t)');
legend('f(t) = sin(220 \pi t) + sin(440 \pi t) + sin(880 \pi t)', 'Location', 'best');

% Show grid
grid on;
