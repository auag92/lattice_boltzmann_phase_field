f       = dir('phi_*.dat');
len     = length(f);
left    = zeros(len,1);
right   = zeros(len,1);
delta   = zeros(len,1);
data    = zeros(len,5);
My = 60;
for k = 1:len
    A = load(f(k).name);
    c = 0;
    indx = 1;
    while(c == 0)
        if(A(indx,3) < 0.5)            
            c = 1;
            data(k,4) = A(indx,3);
        else
            indx = indx+My;
        end
    end
    sl = A(indx,3) - A(indx-My,3);
    data(k,1) = (indx-1)/My + (0.5 - A(indx,3))/sl;
    indx = My;
    c = 0;
    while(c == 0)
        if(A(indx,3) < 0.5)
            c = 1;
            data(k,5) = A(indx,3);
        else
            indx = indx+My;
        end
    end
    sl          = A(indx,3) - A(indx-My,3);
    data(k,2)   = (indx-My)/My + (0.5 - A(indx,3))/sl;
    data(k,3)   = (data(k,1)-data(k,2))/2.0;    
end
fname = sprintf('%d.dat',My);
dlmwrite(fname, data);
