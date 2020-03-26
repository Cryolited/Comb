clear        %%%%%%%%%%%% �������� ������
close all;

A= sqrt(0.5*erfc(-3:0.2:400)); % ��� �������
N=length(A);

n=0:(N/2-1);
A(N-n)=conj(A(2+n));
A(1+N/2)=0;


b(1:800) = 1;   % ������
b(801:N) = 0;


%% ������ ������ (������)

a = real(fft(A)); % ���. ���-�� �������
a = fftshift(a);
%plot((a));

subplot(2,1,1)
c = filter(a,1,b); % ����������(�������)
plot(abs(c));

subplot(2,1,2)
cy = fft(c);    % ������ ���. �������
cy = fftshift(cy);
plot(abs(cy));

%% ������ ������

B = (fft(b)); % ������ �������
%plot(real(B));

C =(fftshift(B)) .* fftshift(A); % ���������� (������������)
%C = fftshift(C);
subplot(2,1,2)
plot(abs(C));

subplot(2,1,1)
CY = fft(C);    % ���. ������
%CY = fftshift(CY);
plot(abs(CY));

%% ��������� ���-���
CYF = abs(((CY)));
CYF = CYF/max(CYF);
c = abs(c/max(c));
figure
plot(c)
hold on
plot(CYF)

sum(CYF - c) % �������� �� �������
