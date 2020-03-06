addpath('src');
profile on;

for q= 1:1000
delay_test;
end

profile off;
profile viewer;

clear
 