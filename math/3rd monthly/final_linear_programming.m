import matlab.io.*

problem = optimproblem('ObjectiveSense', 'minimize' );
%minimize when fat, and supper for normal;
%maximize when protein and breakfast, lunch for normal

%different in different meals
carbohydrate_min = 43;
carbohydrate_max = 83;
protein_min = 15;
protein_max = 65;
fat_min = 21.4;
fat_max = 30.2;

%carbs, protein, fat in 1 gram of 1st, 2nd, 3rd dish
carbohydrates = [0.04, 0.07, 0.04];
protein = [0.03, 0.08, 0.02];
fat = [0, 0.05, 0.04];

x = optimvar('x', 'LowerBound', 0, 'UpperBound', 400);
y = optimvar('y', 'LowerBound', 0, 'UpperBound', 400);
z = optimvar('z', 'LowerBound', 0, 'UpperBound', 400);

%obese objective function: fat(1) * x + fat(2) * y + fat(3) * z 
%fitness objective function: protein(1) * x+ protein(2) * x+protein(3) * z 
%normal objective function: x + y + z

objective = fat(1) * x + fat(2) * y + fat(3) * z;
problem.Objective = objective;

problem.Constraints.carbohydrate_min = carbohydrates(1) * x + carbohydrates(2) * y + carbohydrates(3) * z >= carbohydrate_min;
problem.Constraints.carbohydrate_max = carbohydrates(1) * x + carbohydrates(2) * y + carbohydrates(3) * z <= carbohydrate_max;
problem.Constraints.protein_min = protein(1) * x + protein(2) * y + protein(3) * z >= protein_min;
problem.Constraints.protein_max = protein(1) * x + protein(2) * y + protein(3) * z <= protein_max;
problem.Constraints.fat_min = fat(1) * x + fat(2) * y + fat(3) * z >= fat_min;
problem.Constraints.fat_max = fat(1) * x + fat(2) * y + fat(3) * z <= fat_max;


%constraints of mass of different dishes
problem.Constraints.x_min = x >= 100;
problem.Constraints.y_min = y >= 100;
problem.Constraints.z_min = z >= 100;
problem.Constraints.x_max = x <= 400;
problem.Constraints.y_max = y <= 400;
problem.Constraints.z_max = z <= 400;

[solution, fval] = solve(problem);

disp("Optimal Solution:");
disp("x = " + solution.x);
disp("y = " + solution.y);
disp("z = " + solution.z);
disp("Objective Value: " + fval);

%data plot
num_samples = 100;
x_range = linspace(100, 400, num_samples);
y_range = linspace(100, 400, num_samples);
z_range = linspace(100, 400, num_samples);
[X, Y, Z] = meshgrid(x_range, y_range, z_range);

valid_points = carbohydrates(1) * X + carbohydrates(2) * Y + carbohydrates(3) * Z >= carbohydrate_min & ...
               carbohydrates(1) * X + carbohydrates(2) * Y + carbohydrates(3) * Z <= carbohydrate_max & ...
               protein(1) * X + protein(2) * Y + protein(3) * Z >= protein_min & ...
               protein(1) * X + protein(2) * Y + protein(3) * Z <= protein_max & ...
               fat(1) * X + fat(2) * Y + fat(3) * Z >= fat_min & ...
               fat(1) * X + fat(2) * Y + fat(3) * Z <= fat_max;

K = convhulln([X(valid_points(:)), Y(valid_points(:)), Z(valid_points(:))]);

trisurf(K, X(valid_points), Y(valid_points), Z(valid_points), 'FaceColor', 'blue', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on;

scatter3(solution.x, solution.y, solution.z, 100, 'r', 'filled');

center_x = mean(X(valid_points));
center_y = mean(Y(valid_points));
center_z = mean(Z(valid_points));

scatter3(center_x, center_y, center_z, 100, 'g', 'filled');

text(center_x, center_y, center_z, sprintf('(%0.2f, %0.2f, %0.2f)', center_x, center_y, center_z), 'FontSize', 10, 'Color', 'k');

text(solution.x, solution.y, solution.z, sprintf('(%0.2f, %0.2f, %0.2f)', solution.x, solution.y, solution.z), 'FontSize', 10, 'Color', 'k');

xlabel('x');
ylabel('y');
zlabel('z');
title('Feasible Region, Optimal Solution, and Geometric Center');
legend('Feasible Region', 'Optimal Solution', 'Geometric Center');
grid on;
