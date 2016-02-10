mesh = 50;
for i = 26:36
  % Read txt into cell A
  fid = fopen('constants.h','r+');
  j = 1;
  tline = fgetl(fid);
  A{j} = tline;
  while ischar(tline)
      j = j+1;
      tline = fgetl(fid);
      A{j} = tline;
  end
  fclose(fid);

  replace_string = sprintf("#define My       %d", mesh);
  A{1,3} =  replace_string;
  replace_string = sprintf("#define ftag %d", i);
  A{1,36} =  replace_string;
  mesh = mesh + 5;

  fid = fopen('constants.h','w');
  % Write cell A into txt
  for j = 1:numel(A)
      if A{j+1} == -1
          fprintf(fid,'%s',A{j});
          break
      else
          fprintf(fid,'%s\n', A{j});
      end
  end
  fclose(fid);
  fname = sprintf('mkdir datafiles%d',i);
  system('echo Hello World!%d',i);
  system(fname);
  system('gcc 2d_binary_lbm_mpi.c -lm');
  system('./a.out');
end
