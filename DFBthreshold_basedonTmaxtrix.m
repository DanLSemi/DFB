%%���������Լ��ٵ����ӣ�������F 3.26������ȫ��ͬ��û�ҵ����жԲ�����ֱ��˵����

%��ز�������
deltaLg = (-15:0.01:15)'; %����F 3.26
kappaLg = 0.5:0.5:4.0; %����F 3.26
n1=3.2515;n2=3.456; %ȡ��InGaAsP/InP������
r=(n2-n1)/(n1+n2);
t=sqrt(1-r^2);
m=kappaLg/(2*r);

%����������DFBģʽ����ֵ�����Ĺ���
for n=1:1:length(kappaLg)
  T11(:,n) = -(exp(j*deltaLg./m(n))+r^2)/t^2;
  T22(:,n) = -(r^2+exp(-j*deltaLg./m(n)))/t^2;
  T21(:,n) = (r/t^2)*(exp(j*deltaLg./m(n))+1);
  xi(:,n) = log((T11(:,n)+T22(:,n))/2+sqrt(((T11(:,n)+T22(:,n))/2).^2-1));
  coshmxi(:,n)=(exp(m(n)*xi(:,n))+exp(-m(n)*xi(:,n)))/2;
  sinhmxi(:,n)=(exp(m(n)*xi(:,n))-exp(-m(n)*xi(:,n)))/2;
  tanhmxi(:,n)=sinhmxi(:,n)./coshmxi(:,n);
  tanhxi(:,n)=(exp(xi(:,n))+exp(-xi(:,n)))./(exp(xi(:,n))-exp(-xi(:,n)));
  k(:,n)=(T22(:,n)-T11(:,n))./(T22(:,n)+T11(:,n));  

end

%S21���ֵ��Ӧ�����ľ��ǳ���ģʽ����Ӧ��ʽ�� [��ʽ�Ƶ���¼].docx��
S21=(coshmxi - (sinhmxi./tanhxi).*k).^-1;

%������ֵ����(��g-��)Lg = ln(1/rg*rg')����Ӧ��ʽ�� [��ʽ�Ƶ���¼].docx��
m_eff=tanhmxi./tanhxi;
thres=-2*log((T21./T11).*m_eff.*(1-k)./(1-m_eff.*k));

%% (��g-��)Lg�����ĸ������ȣ���Ӧ��ͬģʽ����ֵ����˶Բ���ȡ�෴����Ѱ��
% Ѱ�庯���������õò��Ǻܺã��������ֶ�ȥ��һЩ��ֵĵ�
pks=zeros(8,length(kappaLg));
locs=zeros(8,length(kappaLg));
MINPEAK=[-7.1 -6 -5.5 -4.5 -4.2 -7.1 -7.1 -7.1]; %�����ǹ۲���ֵ�����������һ�ڳ��� �˴���ע�Ͳ����У�֮���ֶ���ӵ����� 
for n=1:1:length(kappaLg)
    [pks(:,n) locs(:,n)]=findpeaks(-abs(thres(:,n)),'NPEAKS',8,'MINPEAKHEIGHT',MINPEAK(n));
end
% ����ģʽ����һ������
a=zeros(6,8);b=zeros(6,8);
a(1:6,1:3)=pks(2:7,1:3);
a(1:3,4:8)=pks(1:3,4:8);a(4:6,4:8)=pks(6:8,4:8);
b(1:6,1:3)=locs(2:7,1:3);
b(1:3,4:8)=locs(1:3,4:8);b(4:6,4:8)=locs(6:8,4:8);

%% ����Ѱ�����Ƿ�׼ȷ
% figure(2);
% for n=1:1:length(kappaLg)
%     plot((deltaLg),-abs(thres(:,n)));
%     ylim([-30 0])
%     pause
% end
%% ��ͼ(δ���б�׼���������ֻ��ʾ��ͼ)
% ͼ1. S21��ʧг���Ĺ�ϵͼ S21 vs ��Lg ��ע������F3.26 ѡkappaLg=1Ϊ����
% ͼ2. ��ֵ��ʧг���Ĺ�ϵͼ (��g-��)Lg vs ��Lg
figure(1);
Lf=plot(deltaLg,abs(S21(:,2))); % ��kappaLg=1ʱ��S21ʾ��ͼ����ӦF3.26a�ϰ���ͼ
set(Lf,'Linewidth',2);
fsize=20;
xlabel('\delta L_g','FontSize',fsize);
ylabel('S21(\kappa L_g = 1)','FontSize',fsize);
set(gca,'FontSize',fsize); 

figure(2);
%��㣨��ͬkappaLg�������£���ͬģʽ����ֵ�仯��
Ld=plot(deltaLg(b),abs(a),'h');
hold on
%���ߣ���ͬkappaLg�������£�ĳ��ģʽ����ֵ�����ߣ�
for t=1:1:6
plot(deltaLg(b(t,1:8)),abs(a(t,1:8)),'-k') ;
hold on
end
set(Ld,'Linewidth',2);
xlabel('\delta L_g','FontSize',fsize);
ylabel('(��g_t_h - \alpha ) L_g','FontSize',fsize);
set(gca,'FontSize',fsize); 

