mesh = input(">>>");
% Read txt into cell A
fid = fopen('constants.h','r+');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);

A{1,5}
replace_string = sprintf("#define Mx       %d", mesh);
A{1,5} =  replace_string;
A{1,5}

fid = fopen('constants.h','w');
% Write cell A into txt
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s',A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end
fclose(fid);
