fname   = sprintf('2.dat')
A       = load(fname);
v       = zeros(201,2);
for i = 1:201
    i
    v(i,2) = (A(i+100,1) - A(i,1))/(2000);
    v(i,1) = i;
end
% mean of nonzero elements
n = sum(v~=0);
n(n==0) = NaN;
sum(v) ./ n
plot(v(:,1),v(:,2),':b*')