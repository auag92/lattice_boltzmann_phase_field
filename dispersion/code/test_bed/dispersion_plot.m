My      = 50;
len     = 19;
delta   = zeros(len,2);
vel     = zeros(len,1);
l = 598;
v = zeros(l,3);
for k = 1:len
    fname = sprintf('%d.dat',My)
    A = load(fname);
%     for i = 1:100
%         v(i) = (A(i+200,3) - A(i,3))/(4000*0.55*(A(i,3)+A(i+200,3)));
%     end
    
    for i = 1:l
        i
        v(i,2) = (A(i+1,3) - A(i,3))/(1000*0.02*.5*(A(i+1,3)+A(i,3)));
        v(i,3) = (A(i+1,1) - A(i,1))/(1000*0.02);
        v(i,1) = i;
    end
%     plot(v(:,1),v(:,2),':b*')
%     % mean of nonzero elements
%     n = sum(v~=0);
%     n(n==0) = NaN;
%     v_sort  = sort(v, 'descend');
%     delta(k,2) = sum(v) ./ n;
    delta(k,1) = pi/My;
    i_init  = 450;
    i_final = l;
    for i = i_init : i_final
        delta(k,2) = delta(k,2) + v(i,2);
        vel(k) = vel(k) + v(i,3);
    end
    vel(k)     = vel(k)/(i_final - i_init);
    delta(k,2) = delta(k,2)/(i_final - i_init);
    My = My + 5;
end
figure
plot(delta(:,1),delta(:,2),'--bo');
axis([0 0.07 -inf inf])
xlabel('Frequency');
ylabel('Amplification_Factor');
fname = sprintf('delta.dat');
dlmwrite(fname, delta);