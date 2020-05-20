function T = Tmatrix_interface(phi)
% 实际r应该考虑色散，
% 就是失谐量deltaLg表示波长不满足布拉格条件，而是一个范围，那么n应该也是随波长变化的
T11 = exp(j*phi);
T12 = 0;
T21 = 0;
T22 = exp(-j*phi);
T = [T11 T12;T21 T22];