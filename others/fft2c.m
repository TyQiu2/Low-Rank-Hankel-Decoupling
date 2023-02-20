function res = fft2c(x)

% res = fft2c(x)
% 
% orthonormal forward 2D FFT
%
% (c) Michael Lustig 2005

% res = 1/sqrt(length(x(:)))*fftshift(fft2(ifftshift(x)));
 res=fftshift(fft2(x))/sqrt(length(x(:)));
%  res=(fft2(x))/sqrt(length(x(:)));

