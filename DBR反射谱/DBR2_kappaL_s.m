% 固定κLg，改变r以改变光栅对数m的布拉格反射镜反射系数图
% 较大kappaLg,取4

clear;
%参数定义
deltaLg = (-10:0.01:10)'; %归一化失谐参数（代表波长变化范围）
kappaLg = 4; %κLg
lambda=1.55; %按1.55um为布拉格波长设计反射镜
m = [3 4 5 9 100] ; %实际光栅周期应该是整数
r=kappaLg./(2*m);
t=sqrt(1-r.^2);

% 反射系数 rg = (T21/T11)./(1-sinh((m-1)*xi)./(T11*sinh(m*xi)))  式A7.18
% 计算rg需要的参数：单个光栅周期的传输矩阵 中的 T11 T21 T22
% 单个光栅周期的传输矩阵 可以看作 12面+2线+21面+1线 【调用Tmatrix_interface/line函数算】
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
        k(a,b) = (T11(a,b)+T22(a,b))/2; % 简化xi表达式的局域变量
        xi(a,b) = log( k(a,b) + sqrt(k(a,b).^2-1) );
        mxi(a,b) = xi(a,b)*m(b);
    end
end
rg = (T21./T11)./(1-sinh(mxi-xi)./(T11.*sinh(mxi)));
%% 依次绘图（分辨各条线） 
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
title('DBR反射系数','FontSize',fsize);
xlim([-10 10]);ylim([0 1.2]);
text(-3.3,0.7,'m = 3(蓝),4(绿),5(红),9(青),100(紫)','FontSize',fsize-3);
% text(-2.7,0.5,'3','FontSize',fsize);
% text(-2.7,0.3,'30','FontSize',fsize);
text(-1.07,1.1,'κL_g= 4','FontSize',fsize-3);
set(gca,'FontSize',fsize); 
