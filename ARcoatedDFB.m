%%本程序是自己举的例子，和书中F 3.26参数不全相同（没找到书中对参数的直接说明）
% 程序中kappaLg的改变是通过增加光栅对数m实现的
% 阈值增益通过 扫描不同增益 对应的S21获得，因此考虑结果精度，本程序运行时间会较长【练习中的精度为0.25】
% 增益扫描法也适用于寻找其他激光器的模式阈值

tic;
clear;
%相关参数定义
deltaLg = (-15:0.05:15)'; %按照F 3.26，想起码画出三节模
lambda0 = 1.55; % 按布拉格长度1.55um设计光栅
kappaLg = [0.5 1 1.5 2 3 4]; %按照F 3.26
n1=3.2515;n2=3.456; %以InGaAsP/InP为例
r=(n1-n2)/(n1+n2);
t=sqrt(1-r^2);
m=round(kappaLg/(-2*r)); % 考虑到实际光栅对数是正整数
grating_period = lambda0/(4*n1)+lambda0/(4*n2);  % 光栅周期
Lg=m*grating_period;

r_AR = -0.1; r_cleavage = sqrt(0.32); r_HR = 0.9;% 仿 F3.28 的R值（直接定的r，后来没改）
% % 设计增透膜材料折射率为2，结构为λ/4,即传输线经历pi/2，最终r约0.1
% n_AR=2;
% r_n1toAR=(n1-n_AR)/(n1+n_AR);r_ARtoAir=(n_AR-1)/(1+n_AR);
% L_AR=lambda0/(4*n_AR); 
% 
% % 设计增反膜材料折射率为10，结构为λ/4,即传输线经历pi/2，最终R约0.87
% n_HR=10;
% r_AirtoHR=(1-n_HR)/(1+n_HR);r_HRton2=(n2-n_HR)/(n2+n_HR);
% L_HR=lambda0/(4*n_HR); 
% 
% %解理面根据书中设计，R=0.32
% r_cleavage = sqrt(0.32);
% T_cleavage = Tmatrix_interface(r_cleavage);

% 增益项扫描项 g_scan 代表(Γg-α)Lg【参考F3.26】
g_scan = 0:0.01:2.5; 

% 传输面矩阵
T_interface12 = Tmatrix_interface(r);
T_interface21 = Tmatrix_interface(-r);
T_AR = Tmatrix_interface(r_AR);
T_HR = Tmatrix_interface(r_HR);
T_cleavage = Tmatrix_interface(r_cleavage);

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
            T{c,b,a} = T_interface12*T_line2*T_interface21*T_line1; % 重复单元T矩阵
            %Tend{c,b,a} = T_interface12*T_line2*T_interface21*Tmatrix_line(phi1/2); % 最后一个单元T矩阵，书上标的只有λ/8
            
%             T_AR = Tmatrix_interface(r_n1toAR)*Tmatrix_line(pi/2 + L_AR*deltaLg(c)/(Lg(b)))*Tmatrix_interface(r_ARtoAir);
%             %T_ARr = T_AR*T_line2*T_interface21*T_line1; % AR材料之后的第一个周期（第一个光栅单元 反射率有变化）
%             T_HR = Tmatrix_interface(r_HRton2)*Tmatrix_line(pi/2 + L_HR*deltaLg(c)/(Lg(b)))*Tmatrix_interface(r_AirtoHR);
%             T_HRr = T_HR*T_line2*T_interface21*T_line1;% HR材料之后的第一个周期（第一个光栅单元 反射率有变化）
            
            %Tg{c,b,a} = T_AR*(T{c,b,a}^(m(b)))*T_interface12*Tmatrix_line(phi1/2)*T_cleavage; % 解理面
            Tg{c,b,a} = T_AR*(T{c,b,a}^(m(b)))*T_interface12*Tmatrix_line(phi1/2)*T_HR; % HR面
%             Tg{c,b,a} = T_HRr*(T{c,b,a}^(m(b)-1))*T_interface12*Tmatrix_line(phi1/2)*T_AR;
            S21(c,b,a) = 1/Tg{c,b,a}(1);
        end
    end    
end

%% 绘制S21示意图
figure(1);
Lf=plot(deltaLg,abs(S21(:,2,79))); % 画kappaLg=1时的S21示意图，对应F3.26a上半张图【增益值从展示效果方面选的0.78】
set(Lf,'Linewidth',2);
fsize=20;
%xlim([-25 25]);
ylim([0 9]);
xlabel('\delta L_g','FontSize',fsize);
ylabel('S21(\kappa L_g = 1)','FontSize',fsize);
set(gca,'FontSize',fsize); 

%% 依次画线(画出每个kappaLg(即m)下 扫描增益的S21结果)，为了观察寻峰 
% figure(2);
% for b=6:1:length(m)
%     for a=34:1:length(g_scan)     
%         plot(deltaLg,abs(S21(:,b,a)));
%         ylim([0 50]);
%         text(-13.5,40,['κL_g=' num2str(kappaLg(b)) ' ' 'b=' num2str(b)],'FontSize',12);
%         text(-13.5,30,['(Γg-α_i)L_g=' num2str(g_scan(a)) ' ' 'a=' num2str(a)],'FontSize',12);
%         pause;
%     end
% end
%% (利用上一节程序)观察S21，记录阈值点
% 以下寻峰程序只是记录下自己观察得到的阈值结果，并非由程序判断哪个模式到了阈值，因为实现起来有点困难
% 观察1-3阶模的结果记录在aa矩阵中，也列在下方供参考
%  【0表示对应模式没有激射】【只选到三阶模】
%       b=1: 133 179 201 207 213 215 217 0
%       b=2: 93 146 180 192 204 208 0 214
%       b=3: 67 114 156 171 189 195 202 204
%       b=4: 54 92 138 154 175 182 194 198 
%       b=5: 40 61 108 123 150 159 175 0
%       b=6: 33 42 86 98 128 135 155 159
 aa = zeros(6,length(kappaLg));
 [row col]=size(aa);
 aa(:,1)=[133 179 201 207 213 215];aa(:,2)=[93 146 180 192 204 208];aa(:,3)=[67 114 156 171 189 195];
 aa(:,4)=[54 92 138 154 175 182];aa(:,5)=[40 61 108 123 150 159];aa(:,6)=[33 42 86 98 128 135];
for b=1:1:col
    for a=1:1:row
        maxS21=max(abs(S21(:,b,aa(a,b))));
        [locs{a,b}]=find(maxS21==abs(S21(:,b,aa(a,b))));
    end
end
 
%% 无需调整结果阈值点矩阵坐标，直接绘图

figure(2);
%描点1（不同kappaLg的条件下，不同模式对应的阈值）
Ld=plot(deltaLg(cell2mat(locs)),g_scan(aa),'h');
hold on

%连线（不同kappaLg的条件下，某个模式的阈值连成线）
for t=1:1:6
    plot(deltaLg(cell2mat(locs')),g_scan(aa'),'-k') ;
    hold on
end
set(Ld,'Linewidth',2);
xlabel('\delta L_g','FontSize',fsize);
ylabel('(Γg_t_h - \alpha ) L_g','FontSize',fsize);
set(gca,'FontSize',fsize); 

toc;