f       = dir('phi_*.dat');
len     = length(f);
left    = zeros(len,1);
right   = zeros(len,1);
delta   = zeros(len,1);
data    = zeros(len,5);
My = 75;
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
    left(k,1) = (indx-1)/My;
    data(k,1) = (indx-1)/My;
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
    right(k,1)  = (indx-1)/My;
    data(k,2)   = (indx-1)/My;
    delta(k,1)  = (left(k)-right(k))/2.0;
    data(k,3)   = (left(k)-right(k))/2.0;    
end

