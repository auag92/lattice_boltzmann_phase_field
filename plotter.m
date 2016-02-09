d = load('exp.dat');

figure
subplot(2,2,1)       % add first plot in 2 x 2 grid
axis('tight')
plot(d(:,2),d(:,3),':g*')
title('DeltaH vs T')

subplot(2,2,2)       % add second plot in 2 x 2 grid
plot(d(:,1),d(:,3),'--ro')       % plot using + markers
title('DeltaH vs. pH')

subplot(2,2,3)       
plot(d(:,1),d(:,2),'--bo')       
title('T vs. pH')
