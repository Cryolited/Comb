function [ Y_sum ] = non_maximally_decimated_fb( pulse_sig, NFFT, overlapped_ratio, h_fir, round_fir, round_fft )
	pulse_sig_round = pulse_sig(1:NFFT*floor(length(pulse_sig)/NFFT));
	Y_sum = zeros(NFFT,overlapped_ratio*ceil(length(pulse_sig_round)/NFFT));
	for k = 1:1:overlapped_ratio
		pulse_sig_phase_n = pulse_sig_round(1+NFFT/overlapped_ratio*(k-1):end);
		Y_phase_n = maximally_decimated_fb( pulse_sig_phase_n, NFFT, overlapped_ratio, k-1, h_fir, round_fir, round_fft );
		for n = 1:1:length(Y_phase_n(1,:))
			Y_sum(:,k+(n-1)*overlapped_ratio) = Y_phase_n(:,n);
		end
	end
end