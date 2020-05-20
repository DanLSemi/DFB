% File name: fmuller.m
function out = fmuller(x)
% Test functions for Muller¡¯s method
%out = x.^6 - 2;
n1=3.2515;n2=3.456;
r=(n2-n1)/(n1+n2);
t=sqrt(1-r^2);
kappaLg=1;
m=kappaLg/(2*r);
T11 = -(exp(j*x./m)+r^2)/t^2;
T22 = -(r^2+exp(-j*x./m))/t^2;
xi = log((T11+T22)/2+sqrt(((T11+T22)/2).^2-1));
coshmxi=(exp(m*xi)+exp(-m*xi))/2;
sinhmxi=(exp(m*xi)-exp(-m*xi))/2;
tanhmxi=sinhmxi./coshmxi;
tanhxi=(exp(xi)+exp(-xi))./(exp(xi)-exp(-xi));
out = (tanhmxi./tanhxi).*(T22-T11)./(T22+T11)-1;