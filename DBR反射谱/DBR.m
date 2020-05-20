% �̶����ϣ�r��,���ӹ�դ����m������kappaLg�Ĳ������侵����ϵ��ͼ

clear;
%��������
deltaLg = (-15:0.01:15)'; %��һ��ʧг�������������仯��Χ��
kappaLg = [0.5 1 2 4]; %��Lg
lambda=1.55; %��1.55umΪ�����񲨳���Ʒ��侵
n1=3.17;n2=3.40; %��InP/InGaAsPΪ��
r=(n2-n1)/(n1+n2);
t=sqrt(1-r^2);
m = round(kappaLg/(2*r)); %ʵ�ʹ�դ����Ӧ��������
% L1=1.55/(4*n1);L2=1.55/(4*n2);
% grating_period = L1+L2;  % ��դ����
% Lg=m*grating_period;

% ����ϵ�� rg = (T21/T11)./(1-sinh((m-1)*xi)./(T11*sinh(m*xi)))  ʽA7.18
% ����rg��Ҫ�Ĳ�����������դ���ڵĴ������ �е� T11 T21 T22
% ������դ���ڵĴ������ ���Կ��� 12��+2��+21��+1�� ������Tmatrix_interface/line�����㡿
for a=1:1:length(deltaLg)
    for b=1:1:length(kappaLg)
        phi1(a,b) = pi/2 + deltaLg(a)./(2*m(b));
        phi2(a,b) = pi/2 + deltaLg(a)./(2*m(b));
        T{a,b} = Tmatrix_interface(r)*Tmatrix_line(phi1(a,b))*...
            Tmatrix_interface(-r)*Tmatrix_line(phi2(a,b));
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
xlim([-16 16]);ylim([0 1.2]);
text(-1.07,1.07,'��L_g= 4','FontSize',fsize);
text(-0.5,0.9,'2','FontSize',fsize);
text(-0.5,0.7,'1','FontSize',fsize);
text(-0.5,0.3,'0.5','FontSize',fsize);
set(gca,'FontSize',fsize); 
