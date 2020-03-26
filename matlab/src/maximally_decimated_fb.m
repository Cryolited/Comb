function [ Y ] = maximally_decimated_fb( pulse_sig, NFFT, overlapped_ratio, phase_fft, h_fir, round_fir, round_fft )
	N = NFFT;
	b = h_fir;
	B = buffer(b,N);
	B = flipud(B);
	X = buffer(pulse_sig,N);
	for k = 1:N
        %F(k,:) = filter(B(k,:),1,X(k,:));
        F(k,:) = B(k,:) .* X(k,[1:8]);
		F(k,:) = fpga_round(F(k,:), round_fir);
	end
	for k = 1:1:length(F(1,:))
		F(:,k) = circshift(F(:,k),phase_fft*NFFT/overlapped_ratio);
	end
	Y = F;%fft(F);
	Y = fpga_round(Y, round_fft);
end