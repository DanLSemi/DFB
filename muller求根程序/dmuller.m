% File name: dmuller.m
% Driver to test Muller¡¯s method
%clear all
format short
max = 1000;
epsilon = 1e-6;
% starting points
x0 = 0; % original version from the book is wrong about x1, x2, x3
x1 = 1.5;
x2 = 2.5;
% call to Muller¡¯s method
out = muller(@fmuller, x0, x1, x2 , epsilon, max)