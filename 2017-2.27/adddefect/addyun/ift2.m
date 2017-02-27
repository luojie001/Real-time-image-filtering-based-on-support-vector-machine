function [ g ] = ift2( G,delta_f )
%function g = ift2(G,deltaa_f)
%   Calculate the inverse fourier transform of a defined function or
%   anything
N = size(G,1);
g = ifftshift(ifft2(ifftshift(G)))*(N*delta_f)^2;
end

