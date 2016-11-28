function [ p ] = av( v, sv )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

p = 1./(abs(v)+sv)./(abs(v)+sv);

%p = v.*exp(-2.*v./sv);

end

