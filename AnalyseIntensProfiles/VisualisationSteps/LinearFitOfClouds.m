Blue = MTmax_ForBiggerCytMax1197;
x_Blue = AvInt1197;

Red = Result11766;
x_Red = AverIntens11766;

p_Red = polyfit(x_Red, Red, 1);
p_Blue = polyfit(x_Blue, Blue, 1);

plot(x_Blue, Blue', 'bs', 'MarkerSize', 2);
hold on
plot(x_Red, Red, 'rs', 'MarkerSize', 2.5);
hold on
plot(x_Red, p_Red(1) * x_Red + p_Red(2), 'r-');
plot(x_Blue, p_Blue(1) * x_Blue + p_Blue(2), 'b-');

hold off