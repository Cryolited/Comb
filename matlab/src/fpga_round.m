function [ dout ] = fpga_round( din, shval )
	dout = floor((din)./2.^shval + 0.5 + 1i.*0.5);
end