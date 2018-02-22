
loop = 10;

current = zeros(1,loop);
for i =2:1:loop
    current(i) = part2c(1/i);
end

figure(11)
plot(1:1:loop,current);
title('Bottleneck Narrowing')
xlabel('Spacing of L*1/x')
ylabel('Current')
xlim([2 loop])