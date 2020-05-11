%%本程序是自己举得例子，和书中F 3.26参数不全相同（没找到书中对参数的直接说明）

%相关参数定义
deltaLg = (-15:0.01:15)'; %按照F 3.26
kappaLg = 0.5:0.5:4.0; %按照F 3.26
n1=3.2515;n2=3.456; %取了InGaAsP/InP的例子
r=(n2-n1)/(n1+n2);
t=sqrt(1-r^2);
m=kappaLg/(2*r);

%传输矩阵法求解DFB模式与阈值函数的过程
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

%S21最大值对应波长的就是出光模式（对应公式见 [公式推导记录].docx）
S21=(coshmxi - (sinhmxi./tanhxi).*k).^-1;

%计算阈值函数(Γg-α)Lg = ln(1/rg*rg')（对应公式见 [公式推导记录].docx）
m_eff=tanhmxi./tanhxi;
thres=-2*log((T21./T11).*m_eff.*(1-k)./(1-m_eff.*k));

%% (Γg-α)Lg函数的各个波谷，对应不同模式的阈值。因此对波谷取相反数再寻峰
% 寻峰函数参数设置得不是很好，所以有手动去除一些奇怪的点
pks=zeros(8,length(kappaLg));
locs=zeros(8,length(kappaLg));
MINPEAK=[-7.1 -6 -5.5 -4.5 -4.2 -7.1 -7.1 -7.1]; %这里是观察阈值函数结果（下一节程序 此处被注释不运行）之后，手动添加的限制 
for n=1:1:length(kappaLg)
    [pks(:,n) locs(:,n)]=findpeaks(-abs(thres(:,n)),'NPEAKS',8,'MINPEAKHEIGHT',MINPEAK(n));
end
% 根据模式进行一定调整
a=zeros(6,8);b=zeros(6,8);
a(1:6,1:3)=pks(2:7,1:3);
a(1:3,4:8)=pks(1:3,4:8);a(4:6,4:8)=pks(6:8,4:8);
b(1:6,1:3)=locs(2:7,1:3);
b(1:3,4:8)=locs(1:3,4:8);b(4:6,4:8)=locs(6:8,4:8);

%% 检验寻峰结果是否准确
% figure(2);
% for n=1:1:length(kappaLg)
%     plot((deltaLg),-abs(thres(:,n)));
%     ylim([-30 0])
%     pause
% end
%% 绘图(未进行标准化输出处理，只是示意图)
% 图1. S21与失谐量的关系图 S21 vs δLg 【注：对照F3.26 选kappaLg=1为例】
% 图2. 阈值与失谐量的关系图 (Γg-α)Lg vs δLg
figure(1);
Lf=plot(deltaLg,abs(S21(:,2))); % 画kappaLg=1时的S21示意图，对应F3.26a上半张图
set(Lf,'Linewidth',2);
fsize=20;
xlabel('\delta L_g','FontSize',fsize);
ylabel('S21(\kappa L_g = 1)','FontSize',fsize);
set(gca,'FontSize',fsize); 

figure(2);
%描点（不同kappaLg的条件下，不同模式的阈值变化）
Ld=plot(deltaLg(b),abs(a),'h');
hold on
%连线（不同kappaLg的条件下，某个模式的阈值连成线）
for t=1:1:6
plot(deltaLg(b(t,1:8)),abs(a(t,1:8)),'-k') ;
hold on
end
set(Ld,'Linewidth',2);
xlabel('\delta L_g','FontSize',fsize);
ylabel('(Γg_t_h - \alpha ) L_g','FontSize',fsize);
set(gca,'FontSize',fsize); 

