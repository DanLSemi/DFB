function T = Tmatrix_interface(r)
% ʵ��rӦ�ÿ���ɫɢ��
% ����ʧг��deltaLg��ʾ���������㲼��������������һ����Χ����ônӦ��Ҳ���沨���仯��
t = sqrt(1-r.^2);
T11 = 1/t;
T12 = r/t;
T21 = r/t;
T22 = 1/t;
T = [T11 T12;T21 T22];