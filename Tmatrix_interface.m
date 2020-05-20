function T = Tmatrix_interface(r)
% 实际r应该考虑色散，
% 就是失谐量deltaLg表示波长不满足布拉格条件，而是一个范围，那么n应该也是随波长变化的
t = sqrt(1-r.^2);
T11 = 1/t;
T12 = r/t;
T21 = r/t;
T22 = 1/t;
T = [T11 T12;T21 T22];