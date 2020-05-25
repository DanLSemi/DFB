%%本程序是自己举的例子，和书中F 3.26参数不全相同（没找到书中对参数的直接说明）
clear;
%相关参数定义
deltaLg = (-15:0.05:15)'; %按照F 3.26
kappaLg = [0.5 1 1.5 2 3 4]; %按照F 3.26
n1=3.2515;n2=3.456; %以InGaAsP/InP为例
r=(n2-n1)/(n1+n2);
t=sqrt(1-r^2);
m=round(kappaLg/(2*r));
% Ggth_alphai = 4;
% beta_c = beta_r - j*alphai/2
% beta_c = beta_r - j*(Ggth-alphai)/2
LA = 1.55/(4*n1)+1.55/(4*n2);  % 光栅周期
g_scan=0:0.25:8;

%直接利用书中光栅传输矩阵级联结果 求解DFB模式与阈值增益
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

% 均匀光栅
% S21最大值对应波长的就是出光模式（对应公式见 [公式推导记录].docx）
% S21=(abs(Tg11)).^-1;
% S21=(coshmxi - (sinhmxi./tanhxi) + T11.*sinhmxi./sinhxi).^-1;

%% 检验寻峰结果是否准确
% figure(2);
% for n=1:1:length(kappaLg)
%     for a=1:1:length(g_scan)
%     plot(deltaLg,abs(S21(:,n,a)));
%     ylim([0 500]);
%     text(-14,485,['κL_g=' num2str(kappaLg(n)) ' ' 'b=' num2str(n)],'FontSize',12);
%     text(-14,460,['(Γg-α_i)L_g=' num2str(g_scan(a)) ' ' 'a=' num2str(a)],'FontSize',12);
%     pause
%     end
% end
%% (利用上一节程序)观察S21，记录阈值点
% 以下寻峰程序只是记录下自己观察得到的阈值结果，并非由程序判断哪个模式到了阈值，因为实现起来有点困难
% 观察1-3阶模的结果记录在aa矩阵中，也列在下方供参考
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
 
%% 绘图(未进行标准化输出处理，只是示意图)
% 图1. S21与失谐量的关系图 S21 vs δLg 【注：对照F3.26 选kappaLg=1为例】
% 图2. 阈值与失谐量的关系图 (Γg-α)Lg vs δLg
figure(1);
Lf=plot(deltaLg,abs(S21(:,2,12))); % 画kappaLg=1时的S21示意图，对应F3.26a上半张图
set(Lf,'Linewidth',2);
fsize=20;
ylim([0 9]);
xlabel('\delta L_g','FontSize',fsize);
ylabel('S21(\kappa L_g = 1)','FontSize',fsize);
set(gca,'FontSize',fsize); 

%%
% 调整结果阈值点矩阵坐标，aaa[阈值index]和locs0[deltaLg index]按列;失谐量从负到正；行：m(1->6)
locs0 = ones(2*row,col);aaa=locs0;
for b=1:1:col % 没有想到好办法把元胞赋值上去 先写了这个循环
    for a=1:1:row
       % if 
        locs0(a,b)=length(deltaLg)-locs{row+1-a,b};
        locs0(a+row,b)=locs{a,b};
    end
end
aaa(row+1:end,:)=aa(:,:);
aaa(1:row,:)=aa(row:-1:1,:);

figure(2);
%描点（不同kappaLg的条件下，不同模式的阈值变化）
deltaLg0=sort(deltaLg(locs0)); %为了之后的连线而排序
Ld=plot(deltaLg0',g_scan(aaa'),'h');
hold on
%连线（不同kappaLg的条件下，某个模式的阈值连成线）
for t=1:1:6
    plot(deltaLg0',g_scan(aaa'),'-k') ;
    hold on
end
set(Ld,'Linewidth',2);
xlabel('\delta L_g','FontSize',fsize);
ylabel('(Γg_t_h - \alpha ) L_g','FontSize',fsize);
set(gca,'FontSize',fsize); 

