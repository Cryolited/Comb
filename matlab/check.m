fid = fopen('/home/anatoly/hub/Comb/Anal', 'rb');
if fid == -1                     % проверка корректности открытия
    error('File is not opened');
end 

Anal = fread(fid,  'double');    % чтение 
Anal2 = Anal(1:2:end) + 1j*Anal(2:2:end);
Anal3 = buffer(Anal2, 18).';
%Anal3 = buffer(Anal2, 128);

fclose(fid);

figure
subplot(2,1,1)
mesh(abs(Anal3))
subplot(2,1,2)
mesh(abs(fb_analysis_data))

sum(sum(Anal3 - fb_analysis_data))

%%
fid = fopen('/home/anatoly/hub/Comb/Signal', 'rb');

if fid == -1                     % проверка корректности открытия
    error('File is not opened');
end 

Output = fread(fid,  'double');    % чтение 
Output2 = Output(1:2:end) + 1j*Output(2:2:end);
rez = Output2 - fb_synth_data.';
err = sum(abs(rez))/sum(abs(Output2));

figure
subplot(2,1,1)
plot(real(Output2))
subplot(2,1,2)
plot(real(fb_synth_data))

fclose(fid);   
%%

tic
conv([1:1000],[1:1000]);
toc
tic
filter([1:1000],1,[1:1000]);
toc