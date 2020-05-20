% �̶���Lg���ı�r�Ըı��դ����m�Ĳ������侵����ϵ��ͼ
% �ϴ�kappaLg,ȡ4

clear;
%��������
deltaLg = (-10:0.01:10)'; %��һ��ʧг�������������仯��Χ��
kappaLg = 4; %��Lg
lambda=1.55; %��1.55umΪ�����񲨳���Ʒ��侵
m = [3 4 5 9 100] ; %ʵ�ʹ�դ����Ӧ��������
r=kappaLg./(2*m);
t=sqrt(1-r.^2);

% ����ϵ�� rg = (T21/T11)./(1-sinh((m-1)*xi)./(T11*sinh(m*xi)))  ʽA7.18
% ����rg��Ҫ�Ĳ�����������դ���ڵĴ������ �е� T11 T21 T22
% ������դ���ڵĴ������ ���Կ��� 12��+2��+21��+1�� ������Tmatrix_interface/line�����㡿
for a=1:1:length(deltaLg)
    for b=1:1:length(m)
        phi1(a,b) = pi/2 + deltaLg(a)./(2*m(b));
        phi2(a,b) = pi/2 + deltaLg(a)./(2*m(b));
        T{a,b} = Tmatrix_interface(r(b))*Tmatrix_line(phi1(a,b))*...
            Tmatrix_interface(-r(b))*Tmatrix_line(phi2(a,b));
        T11(a,b) = T{a,b}(1);
        T21(a,b) = T{a,b}(2);
        T12(a,b) = T{a,b}(3);
        T22(a,b) = T{a,b}(4);
        k(a,b) = (T11(a,b)+T22(a,b))/2; % ��xi���ʽ�ľ������
        xi(a,b) = log( k(a,b) + sqrt(k(a,b).^2-1) );
        mxi(a,b) = xi(a,b)*m(b);
    end
end
rg = (T21./T11)./(1-sinh(mxi-xi)./(T11.*sinh(mxi)));
%% ���λ�ͼ���ֱ�����ߣ� 
% for n=1:1:length(m)
%     figure(1)
%     plot(deltaLg,abs(rg(:,n)));
%     pause;
%     hold on;
% end     
% hold off;
%%
Ld = plot(deltaLg,abs(rg));
set(Ld,'Linewidth',2);
fsize = 20;
xlabel('\delta L_g','FontSize',fsize);
ylabel('r_g ','FontSize',fsize);
title('DBR����ϵ��','FontSize',fsize);
xlim([-10 10]);ylim([0 1.2]);
text(-3.3,0.7,'m = 3(��),4(��),5(��),9(��),100(��)','FontSize',fsize-3);
% text(-2.7,0.5,'3','FontSize',fsize);
% text(-2.7,0.3,'30','FontSize',fsize);
text(-1.07,1.1,'��L_g= 4','FontSize',fsize-3);
set(gca,'FontSize',fsize); 
