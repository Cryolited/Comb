
clearvars
clc
Fs = 512e6;
NFFT = 128;
WIN_OVERLAP_RATIO = 8;

%addpath('src');

%% signal params
t_us = 0.2 * 1e-6;
period_us = 10 * 1e-6;
first_period_part = NFFT / Fs;
second_period_part = NFFT * WIN_OVERLAP_RATIO / Fs;
%first_period_part = 5 * 1e-6;
%second_period_part = 5 * 1e-6;
%% generate harmonic signal
DATA_AMPL           = 30000;
THRESHOLD           = DATA_AMPL/2;
time_sec            = t_us;
f0                  = 0; % freq_shift              
t                   = (0:1:Fs*time_sec -1)/Fs; 
% % simple sin signal
% sig                 = exp(1i*2*pi*1e6*t);
% sig_round           = round(DATA_AMPL*sig);
% %signal              = [zeros(1,first_period_part*Fs) sig_round zeros(1,second_period_part*Fs)];
% complex lfm signal
LFM_dev_hz = 50*1e6;
f0                  = -LFM_dev_hz/2;                
b                   = (LFM_dev_hz/time_sec);   
sig                 = exp(1i*2*pi*(f0*t+b/2*t.^2));
sig_round           = round(DATA_AMPL*sig);
signal              = [zeros(1,first_period_part*Fs) sig_round zeros(1,second_period_part*Fs)];
% signal              = [zeros(1,period_us*Fs) sig_round];
% signal              = [ sig_round zeros(1,period_us*Fs)];
%plot(-Fs/2:Fs/length(signal):(Fs/2 - Fs/length(signal)),abs(fftshift(fft(signal))));
%spectrogram(signal,256,250,[],Fs,'yaxis')


fid = fopen('../I', 'rb'); %/home/anatoly/hub/Comb/I
if fid == -1                     % проверка корректности открытия
    error('File is not opened');
end 

signalI = fread(fid,'int16');    % чтение 
fclose(fid); 

fid2 = fopen('../Q', 'rb');
signalQ = fread(fid2,  'int16');    % чтение 
%Output2 = Output(1:2:end) + 1j*Output(2:2:end);
%rez = Output2 - fb_synth_data.';
fclose(fid2); 



%signal = signalI + 1j*signalQ;
f_signal = (fft(signal));

signal = f_signal; % comment this after
ZZ = buffer(signal,NFFT);
mesh(abs(ZZ))
%% analysis filter Bank


WIN_H_RADIX       = 18;
FB_OVERLAP_RATIO  = 2;
BAND              = 500e6;
FW                = 4e6;

% 1.create npr coeff
c=npr_coeff(NFFT,2*WIN_OVERLAP_RATIO);
coeff = c(:);
f_coeff = abs(fftshift(fft(coeff)));

coeff = f_coeff; % comment this after 

max_coeff_val = max(abs(coeff));
coeff_radix = fix(log2(2^(WIN_H_RADIX-1)/max_coeff_val));
h_fb_win_fxp = round(coeff*2^coeff_radix);

% 2.get signal throwgh the filterBank
fb_analysis_win_max_gain_bit = ceil(max(log2(sum(abs(buffer(h_fb_win_fxp,NFFT)),2))));
round_fir = fb_analysis_win_max_gain_bit;
round_fft = coeff_radix-fb_analysis_win_max_gain_bit;
fb_analysis_data = non_maximally_decimated_fb(signal, ...
NFFT, FB_OVERLAP_RATIO, ...
h_fb_win_fxp, round_fir, round_fft);

%rmpath('src');



%% synthesis filter bank
% 
% max_gain_bit = ceil(max(log2(sum(abs(buffer(h_fb_win_fxp,NFFT)),2))));
% round_fft = 1;
% round_fir = max_gain_bit-round_fft;
% fb_synth_data = npr_synthesis(h_fb_win_fxp, fb_analysis_data, round_fir, round_fft);
% 
% 
% 
% 
% 
% plot(real(signal));
% hold on;
% plot(real(fb_synth_data));
% legend('input','output');
% 
% sconv = abs(xcorr(signal,signal));
% xconv = abs(xcorr(signal,fb_synth_data));
% 
% sconv_max = find(sconv == max(sconv));
% xconv_max = find(xconv == max(xconv));
% 
% delay = xconv_max - sconv_max;
% disp([ 'delay = ' num2str(delay/Fs*1e6) ' us']);
% disp([ 'add zeros = ' num2str(first_period_part + second_period_part)  ' us']);
% 
% 
% 
% 
% 
