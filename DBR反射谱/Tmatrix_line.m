function T = Tmatrix_interface(phi)
% ʵ��rӦ�ÿ���ɫɢ��
% ����ʧг��deltaLg��ʾ���������㲼��������������һ����Χ����ônӦ��Ҳ���沨���仯��
T11 = exp(j*phi);
T12 = 0;
T21 = 0;
T22 = exp(-j*phi);
T = [T11 T12;T21 T22];