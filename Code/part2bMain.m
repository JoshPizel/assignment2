
loop = 8;

current = zeros(1,loop);
for i =1:1:loop
    current(i) = part2b(i);
end

figure(10)
plot(1:1:loop,current);
title('Mesh Density')
xlabel('Points Per Unit Spacing')
ylabel('Current')
ylim([-0.5e-6 0])