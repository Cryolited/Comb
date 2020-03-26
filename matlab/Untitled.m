clear        %%%%%%%%%%%% Исходные данные
close all;

A= sqrt(0.5*erfc(-3:0.2:400)); % ачх фильтра
N=length(A);

n=0:(N/2-1);
A(N-n)=conj(A(2+n));
A(1+N/2)=0;


b(1:800) = 1;   % сигнал
b(801:N) = 0;


%% Первый способ (сейчас)

a = real(fft(A)); % Имп. хар-ка фильтра
a = fftshift(a);
%plot((a));

subplot(2,1,1)
c = filter(a,1,b); % фильтрация(свертка)
plot(abs(c));

subplot(2,1,2)
cy = fft(c);    % спектр вых. сигнала
cy = fftshift(cy);
plot(abs(cy));

%% Второй способ

B = (fft(b)); % спектр сигнала
%plot(real(B));

C =(fftshift(B)) .* fftshift(A); % фильтрация (перемножение)
%C = fftshift(C);
subplot(2,1,2)
plot(abs(C));

subplot(2,1,1)
CY = fft(C);    % вых. сигнал
%CY = fftshift(CY);
plot(abs(CY));

%% Сравнение рез-тов
CYF = abs(((CY)));
CYF = CYF/max(CYF);
c = abs(c/max(c));
figure
plot(c)
hold on
plot(CYF)

sum(CYF - c) % разность от Модулей
