%%本程序是自己举的例子，和书中F 3.26参数不全相同（没找到书中对参数的直接说明）
% 程序中kappaLg的改变是通过增加光栅对数m实现的
% 阈值增益通过 扫描不同增益 对应的S21获得，因此考虑结果精度，本程序运行时间会较长【练习中的精度为0.25】
% 增益扫描法也适用于寻找其他激光器的模式阈值

tic;
clear;
%相关参数定义
deltaLg = (-15:0.05:15)'; %按照F 3.26
lambda0 = 1.55; % 按布拉格长度1.55um设计光栅
kappaLg = [0.5 1 1.5 2 3 4]; %按照F 3.26
n1=3.2515;n2=3.456; %以InGaAsP/InP为例
r=(n2-n1)/(n1+n2);
t=sqrt(1-r^2);
m=round(kappaLg/(2*r));
% grating_period = lambda0/(4*n1)+lambda0/(4*n2);  % 光栅周期

% 增益项扫描项 g_scan 代表(Γg-α)Lg【参考F3.26】
g_scan = 0:0.25:8; 

% 传输面矩阵
T_interface12 = Tmatrix_interface(r);
T_interface21 = Tmatrix_interface(-r);

% 以下循环中的过程变量(phi1/2,Tline1/2)不存储所有结果，节省工作空间。
for a=1:1:length(g_scan)
    for b=1:1:length(m)
        for c=1:1:length(deltaLg)
            % βc = β + j(Γg-α)/2  【公式说明】
            % Φ = βc*grating_period = pi + deltaLg/m +...
            % j*g_scan/(2*m) 【公式说明】
            phi1 = pi/2 + deltaLg(c)/(2*m(b)) + j*g_scan(a)/(4*m(b)); 
            phi2 = pi/2 + deltaLg(c)/(2*m(b)) + j*g_scan(a)/(4*m(b));
            % 两个phi相等，参考书中Φ+=Π+δLg/2m 与 Φ-=0
            
            % 传输线矩阵
            T_line1 = Tmatrix_line(phi1);
            T_line2 = Tmatrix_line(phi2);
            T1{c,b,a} = T_interface12*T_line2*T_interface21*T_line1; % 前半段重复单元T矩阵
            T2{c,b,a} = T_line1*T_interface12*T_line2*T_interface21; % 后半段重复单元T矩阵
            T_shift{c,b,a} = T_interface12*Tmatrix_line(2*phi2)*T_interface21; %λ/4相移区T矩阵
            Tg{c,b,a} = (T1{c,b,a}^floor(m(b)/2))*T_shift{c,b,a}*(T2{c,b,a}^floor(m(b)/2));
            S21(c,b,a) = 1/Tg{c,b,a}(1);
        end
    end    
end

%% 绘制S21示意图
figure(1);
Lf=plot(deltaLg,abs(S21(:,2,10))); % 画kappaLg=1时的S21示意图，对应F3.26b上半张图【g_scan增益值从展示效果方面选的2.25】
set(Lf,'Linewidth',2);
ylim([0 13])
fsize=20;
xlabel('\delta L_g','FontSize',fsize);
ylabel('S21(\kappa L_g = 1)','FontSize',fsize);
set(gca,'FontSize',fsize); 

%% 依次画线(画出每个kappaLg(即m)下 扫描增益的S21结果)，为了观察寻峰 
% figure(2);
% for b=1:1:length(m)
%     for a=1:1:length(g_scan)     
%         plot(deltaLg,abs(S21(:,b,a)));
%         ylim([0 400]);
%         text(-14,385,['κL_g=' num2str(kappaLg(b)) ' ' 'b=' num2str(b)],'FontSize',12);
%         text(-14,360,['(Γg-α_i)L_g=' num2str(g_scan(a)) ' ' 'a=' num2str(a)],'FontSize',12);
%         pause;
%     end
% end
%% (利用上一节程序)观察S21，记录阈值点
% 以下寻峰程序只是记录下自己观察得到的阈值结果，并非由程序判断哪个模式到了阈值，因为实现起来有点困难
% 观察0-4阶模的结果记录在aa矩阵中，也列在下方供参考
%       b=1: 16 20 23 25 27
%       b=2: 12 16 21 23 25
%       b=3: 9  14 19 21 24
%       b=4: 6  12 17 19 22
%       b=5: 4  9  14 16 19
%       b=6: 2  7  11 13 16
aa = zeros(5,6);
aa(:,1)=[16 20 23 25 27];aa(:,2)=[12 16 21 23 25];aa(:,3)=[9  14 19 21 24];
aa(:,4)=[6  12 17 19 22];aa(:,5)=[4  9  14 16 19];aa(:,6)=[2  7  11 13 16];
for b=1:1:6
    for a=1:1:5
        maxS21=max(abs(S21(:,b,aa(a,b))));
        [locs{a,b}]=find(maxS21==abs(S21(:,b,aa(a,b))));
    end
end
%% 
% 正负模式坐标赋值与调整
locs0 = ones(9,6);aaa=locs0;
for b=1:1:6 % 不写这个循环没有想到好办法把元胞赋值上去
    for a=1:1:4
        locs0(a,b)=length(deltaLg)-locs{6-a,b};
        locs0(a+5,b)=locs{a+1,b};
    end
end
locs0(5,:)=locs{1,:};
aaa(5:9,:)=aa(:,:);
aaa(1:4,:)=aa(5:-1:2,:);
figure(2);

%描点（不同kappaLg的条件下，不同模式对应的阈值）
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


toc;