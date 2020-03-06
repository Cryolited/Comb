%NPR_SYNTHESIS Near perfect reconstruction synthesis
%  Y = NPR_SYNTHESIS(C,X) synthesizes a wideband signal Y from a
%  number of subbands stored in X. Each subband is a row in X.
%  C is a two dimensional array, containing the filter coefficients.
%  The number of rows in X must be twice the number of rows in C. 
%
% (c) 2007 Wessel Lubberhuizen

function [y_result] = npr_synthesis(coeff, x, round_fir, round_fft, disp_max_val)  
    % number of channels
    N=size(x,1);

    % number of slices
    M=size(x,2);

    coeff_buff = buffer(coeff,N);

    % split into even and odd channels
    y1 = x(:,1:2:end);
    y2 = x(:,2:2:end);

    % apply ifft using fft with inverse input and output data
    y1_iq_swap = imag(y1)+1i.*real(y1);
    y2_iq_swap = imag(y2)+1i.*real(y2);
    y1_fft_inv = fft(y1_iq_swap,[],1);
    y2_fft_inv = fft(y2_iq_swap,[],1);
    y1_fft     = imag(y1_fft_inv)+1i.*real(y1_fft_inv);
    y2_fft     = imag(y2_fft_inv)+1i.*real(y2_fft_inv);
    y1_fft_fxp = fpga_round(y1_fft, round_fft);
    y2_fft_fxp = fpga_round(y2_fft, round_fft);

    % cyclic shift
    y2_sh = [y2_fft_fxp(N/2+1:end,:); y2_fft_fxp(1:N/2,:)];

    % apply channel filters
    y1_filter = zeros(N,M/2);
    y2_filter = zeros(N,M/2);
    for i=1:N
        y1_filter(i,:) = filter(coeff_buff(i,:),1,y1_fft_fxp(i,:));
        y2_filter(i,:) = filter(coeff_buff(i,:),1,y2_sh(i,:));
    end

    y1_result = y1_filter(:).';
    y2_result = y2_filter(:).';
    y1_fir_fxp = fpga_round(y1_result, round_fir);
    y2_fir_fxp = fpga_round(y2_result, round_fir);

    % combine filter results
    y_result = y1_fir_fxp + [zeros(1,N/2) y2_fir_fxp(1:end-N/2)];

    % debug section
    if nargin == 4   % if the number of inputs equals 4
        disp_max_val = 0;
    end
    max_even_fft = max([max(abs(real(y1_fft))) max(abs(imag(y1_fft)))]);
    max_odd_fft = max([max(abs(real(y2_fft))) max(abs(imag(y2_fft)))]);
    max_fft = max([max_even_fft max_odd_fft]);
    if disp_max_val == 1
        disp(['max_fft_value = ' num2str(max_fft)]);
    end

end