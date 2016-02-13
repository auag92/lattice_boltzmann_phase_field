My = 145;
fname   = sprintf('%d.dat',My);
l = 598;
A       = load(fname);
v       = zeros(l,2);
for i = 1:l
    i
    v(i,2) = (A(i+1,3) - A(i,3))/(1000*0.02);
    v(i,1) = i;
end
% mean of nonzero elements
n = sum(v~=0);
n(n==0) = NaN;

i_init  = 1600;
i_final = l;
vel = 0;
for i = i_init : i_final
    vel = vel + v(i,2)
end
vel = vel/(i_final-i_init)
plot(v(:,1),v(:,2),':b*')