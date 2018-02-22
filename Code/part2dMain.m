
loop = 100;

current = zeros(1,loop);
for i =1:1:loop
    sig = i*1e-2;
    current(i) = part2d(sig);
end

figure(12)
plot(1e-2:1e-2:1,current);
title('Sigma Variations')
xlabel('Sigma')
ylabel('Current')