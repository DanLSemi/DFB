%%���������Լ��ٵ����ӣ�������F 3.26������ȫ��ͬ��û�ҵ����жԲ�����ֱ��˵����
clear;
%��ز�������
deltaLg = (-15:0.05:15)'; %����F 3.26
kappaLg = [0.5 1 1.5 2 3 4]; %����F 3.26
n1=3.2515;n2=3.456; %��InGaAsP/InPΪ��
r=(n2-n1)/(n1+n2);
t=sqrt(1-r^2);
m=round(kappaLg/(2*r));
% Ggth_alphai = 4;
% beta_c = beta_r - j*alphai/2
% beta_c = beta_r - j*(Ggth-alphai)/2
LA = 1.55/(4*n1)+1.55/(4*n2);  % ��դ����
g_scan=0:0.25:8;

%ֱ���������й�դ������������ ���DFBģʽ����ֵ����
for b=1:1:length(deltaLg)
for n=1:1:length(kappaLg)
    for a=1:1:length(g_scan)
        phi = pi + deltaLg(b)./m(n)+ j*g_scan(a)./(2*m(n));
        T11(b,n,a) = (exp(j*phi)-r^2)/t^2;
        T22(b,n,a) = -(r^2-exp(-j*phi))/t^2;
        T21(b,n,a) = (r/t^2)*(exp(j*phi)-1);
        T12(b,n,a) = (r/t^2)*(exp(-j*phi)-1);
        xi(b,n,a) = log((T11(b,n,a)+T22(b,n,a))/2 - sqrt(((T11(b,n,a)+T22(b,n,a))/2).^2-1));
        Ddelta(b,n,a)=j*(T22(b,n,a)-T11(b,n,a))./(T22(b,n,a)+T11(b,n,a));
        m_eff(b,n,a)=tanh(m(n)*xi(b,n,a))./tanh(xi(b,n,a));
        % k(b,n,a)=(T22(:,n,a)-T11(:,n,a))./(T22(:,n,a)+T11(:,n,a));         
        Tg11(b,n,a) = (1 + j*m_eff(b,n,a).*Ddelta(b,n,a)).*cosh(m(n).*xi(b,n,a));
        Tg21(b,n,a) = (T21(b,n,a)./T11(b,n,a)).*m_eff(b,n,a).*(1 + j*Ddelta(b,n,a)).*cosh(m(n).*xi(b,n,a));
        Tg12(b,n,a) = (T12(b,n,a)./T22(b,n,a)).*m_eff(b,n,a).*(1 - j*Ddelta(b,n,a)).*cosh(m(n).*xi(b,n,a));
        Tg22(b,n,a) = (1 - j*m_eff(b,n,a).*Ddelta(b,n,a)).*cosh(m(n).*xi(b,n,a));
        %S21(b,n,a)=((cosh(m(n).*xi(b,n,a)) - (sinh(m(n).*xi(b,n,a))./tanh(xi(b,n,a))) + T11(b,n,a).*sinh(m(n).*xi(b,n,a))./sinh(xi(b,n,a)))).^-1;
        %TT1{b,n,a} = [1.1765 -0.62353;-0.62353 1.1765]*[Tg11(b,n,a) Tg12(b,n,a);Tg21(b,n,a) Tg22(b,n,a)]*[1 0;0 1];
        S21(b,n,a)=1./Tg11(b,n,a);
    end
end
end

% ���ȹ�դ
% S21���ֵ��Ӧ�����ľ��ǳ���ģʽ����Ӧ��ʽ�� [��ʽ�Ƶ���¼].docx��
% S21=(abs(Tg11)).^-1;
% S21=(coshmxi - (sinhmxi./tanhxi) + T11.*sinhmxi./sinhxi).^-1;

%% ����Ѱ�����Ƿ�׼ȷ
% figure(2);
% for n=1:1:length(kappaLg)
%     for a=1:1:length(g_scan)
%     plot(deltaLg,abs(S21(:,n,a)));
%     ylim([0 500]);
%     text(-14,485,['��L_g=' num2str(kappaLg(n)) ' ' 'b=' num2str(n)],'FontSize',12);
%     text(-14,460,['(��g-��_i)L_g=' num2str(g_scan(a)) ' ' 'a=' num2str(a)],'FontSize',12);
%     pause
%     end
% end
%% (������һ�ڳ���)�۲�S21����¼��ֵ��
% ����Ѱ�����ֻ�Ǽ�¼���Լ��۲�õ�����ֵ����������ɳ����ж��ĸ�ģʽ������ֵ����Ϊʵ�������е�����
% �۲�1-3��ģ�Ľ����¼��aa�����У�Ҳ�����·����ο�
%       b=1: 22 27 30 32
%       b=2: 15 21 25 27
%       b=3: 11 18 21 23
%       b=4: 9  15 19 21
%       b=5: 6  12 15 18
%       b=6: 4  9  13 15

 aa = zeros(4,length(kappaLg));
 [row col]=size(aa);
 aa(:,1)=[22 27 30 32];aa(:,2)=[15 21 25 27];aa(:,3)=[11 18 21 23];
 aa(:,4)=[9  15 19 21];aa(:,5)=[6  12 15 18];aa(:,6)=[4  9  13 15];
for b=1:1:col
    for a=1:1:row
        maxS21=max(abs(S21(:,b,aa(a,b))));
        [locs{a,b}]=find(maxS21==abs(S21(:,b,aa(a,b))));
    end
end
 
%% ��ͼ(δ���б�׼���������ֻ��ʾ��ͼ)
% ͼ1. S21��ʧг���Ĺ�ϵͼ S21 vs ��Lg ��ע������F3.26 ѡkappaLg=1Ϊ����
% ͼ2. ��ֵ��ʧг���Ĺ�ϵͼ (��g-��)Lg vs ��Lg
figure(1);
Lf=plot(deltaLg,abs(S21(:,2,12))); % ��kappaLg=1ʱ��S21ʾ��ͼ����ӦF3.26a�ϰ���ͼ
set(Lf,'Linewidth',2);
fsize=20;
ylim([0 9]);
xlabel('\delta L_g','FontSize',fsize);
ylabel('S21(\kappa L_g = 1)','FontSize',fsize);
set(gca,'FontSize',fsize); 

%%
% ���������ֵ��������꣬aaa[��ֵindex]��locs0[deltaLg index]����;ʧг���Ӹ��������У�m(1->6)
locs0 = ones(2*row,col);aaa=locs0;
for b=1:1:col % û���뵽�ð취��Ԫ����ֵ��ȥ ��д�����ѭ��
    for a=1:1:row
       % if 
        locs0(a,b)=length(deltaLg)-locs{row+1-a,b};
        locs0(a+row,b)=locs{a,b};
    end
end
aaa(row+1:end,:)=aa(:,:);
aaa(1:row,:)=aa(row:-1:1,:);

figure(2);
%��㣨��ͬkappaLg�������£���ͬģʽ����ֵ�仯��
deltaLg0=sort(deltaLg(locs0)); %Ϊ��֮������߶�����
Ld=plot(deltaLg0',g_scan(aaa'),'h');
hold on
%���ߣ���ͬkappaLg�������£�ĳ��ģʽ����ֵ�����ߣ�
for t=1:1:6
    plot(deltaLg0',g_scan(aaa'),'-k') ;
    hold on
end
set(Ld,'Linewidth',2);
xlabel('\delta L_g','FontSize',fsize);
ylabel('(��g_t_h - \alpha ) L_g','FontSize',fsize);
set(gca,'FontSize',fsize); 

