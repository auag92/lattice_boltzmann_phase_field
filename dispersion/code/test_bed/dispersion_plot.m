My      = 50;
len     = 15;
v       = zeros(150,1);
delta   = zeros(len,2);
v_sort  = zeros(150,1);
for k = 1:len
    fname = sprintf('%d.dat',My)
    A = load(fname);
    for i = 1:200
        v(i) = (A(i+100,3) - A(i,3))/(2000*0.5*(A(i+100,3)+A(i,3)));
    end
    % mean of nonzero elements
    n = sum(v~=0);
    n(n==0) = NaN;
    v_sort  = sort(v, 'descend')
    delta(k,2) = sum(v) ./ n;
    delta(k,1) = pi/My;
    My = My + 5;
end
figure
plot(delta(:,1),delta(:,2),'--bo')
fname = sprintf('delta.dat');
dlmwrite(fname, delta);