%%���������Լ��ٵ����ӣ�������F 3.26������ȫ��ͬ��û�ҵ����жԲ�����ֱ��˵����
% ������kappaLg�ĸı���ͨ�����ӹ�դ����mʵ�ֵ�
% ��ֵ����ͨ�� ɨ�費ͬ���� ��Ӧ��S21��ã���˿��ǽ�����ȣ�����������ʱ���ϳ�����ϰ�еľ���Ϊ0.25��
% ����ɨ�跨Ҳ������Ѱ��������������ģʽ��ֵ

tic;
clear;
%��ز�������
deltaLg = (-15:0.05:15)'; %����F 3.26�������뻭������ģ
lambda0 = 1.55; % �������񳤶�1.55um��ƹ�դ
kappaLg = [0.5 1 1.5 2 3 4]; %����F 3.26
n1=3.2515;n2=3.456; %��InGaAsP/InPΪ��
r=(n1-n2)/(n1+n2);
t=sqrt(1-r^2);
m=round(kappaLg/(-2*r)); % ���ǵ�ʵ�ʹ�դ������������
grating_period = lambda0/(4*n1)+lambda0/(4*n2);  % ��դ����
Lg=m*grating_period;

r_AR = -0.1; r_cleavage = sqrt(0.32); r_HR = 0.9;% �� F3.28 ��Rֵ��ֱ�Ӷ���r������û�ģ�
% % �����͸Ĥ����������Ϊ2���ṹΪ��/4,�������߾���pi/2������rԼ0.1
% n_AR=2;
% r_n1toAR=(n1-n_AR)/(n1+n_AR);r_ARtoAir=(n_AR-1)/(1+n_AR);
% L_AR=lambda0/(4*n_AR); 
% 
% % �������Ĥ����������Ϊ10���ṹΪ��/4,�������߾���pi/2������RԼ0.87
% n_HR=10;
% r_AirtoHR=(1-n_HR)/(1+n_HR);r_HRton2=(n2-n_HR)/(n2+n_HR);
% L_HR=lambda0/(4*n_HR); 
% 
% %���������������ƣ�R=0.32
% r_cleavage = sqrt(0.32);
% T_cleavage = Tmatrix_interface(r_cleavage);

% ������ɨ���� g_scan ����(��g-��)Lg���ο�F3.26��
g_scan = 0:0.01:2.5; 

% ���������
T_interface12 = Tmatrix_interface(r);
T_interface21 = Tmatrix_interface(-r);
T_AR = Tmatrix_interface(r_AR);
T_HR = Tmatrix_interface(r_HR);
T_cleavage = Tmatrix_interface(r_cleavage);

% ����ѭ���еĹ��̱���(phi1/2,Tline1/2)���洢���н������ʡ�����ռ䡣
for a=1:1:length(g_scan)
    for b=1:1:length(m)
        for c=1:1:length(deltaLg)
            % ��c = �� + j(��g-��)/2  ����ʽ˵����
            % �� = ��c*grating_period = pi + deltaLg/m +...
            % j*g_scan/(2*m) ����ʽ˵����
            phi1 = pi/2 + deltaLg(c)/(2*m(b)) + j*g_scan(a)/(4*m(b)); 
            phi2 = pi/2 + deltaLg(c)/(2*m(b)) + j*g_scan(a)/(4*m(b));
            % ����phi��ȣ��ο����Ц�+=��+��Lg/2m �� ��-=0
            
            % �����߾���
            T_line1 = Tmatrix_line(phi1);
            T_line2 = Tmatrix_line(phi2);
            T{c,b,a} = T_interface12*T_line2*T_interface21*T_line1; % �ظ���ԪT����
            %Tend{c,b,a} = T_interface12*T_line2*T_interface21*Tmatrix_line(phi1/2); % ���һ����ԪT�������ϱ��ֻ�Ц�/8
            
%             T_AR = Tmatrix_interface(r_n1toAR)*Tmatrix_line(pi/2 + L_AR*deltaLg(c)/(Lg(b)))*Tmatrix_interface(r_ARtoAir);
%             %T_ARr = T_AR*T_line2*T_interface21*T_line1; % AR����֮��ĵ�һ�����ڣ���һ����դ��Ԫ �������б仯��
%             T_HR = Tmatrix_interface(r_HRton2)*Tmatrix_line(pi/2 + L_HR*deltaLg(c)/(Lg(b)))*Tmatrix_interface(r_AirtoHR);
%             T_HRr = T_HR*T_line2*T_interface21*T_line1;% HR����֮��ĵ�һ�����ڣ���һ����դ��Ԫ �������б仯��
            
            %Tg{c,b,a} = T_AR*(T{c,b,a}^(m(b)))*T_interface12*Tmatrix_line(phi1/2)*T_cleavage; % ������
            Tg{c,b,a} = T_AR*(T{c,b,a}^(m(b)))*T_interface12*Tmatrix_line(phi1/2)*T_HR; % HR��
%             Tg{c,b,a} = T_HRr*(T{c,b,a}^(m(b)-1))*T_interface12*Tmatrix_line(phi1/2)*T_AR;
            S21(c,b,a) = 1/Tg{c,b,a}(1);
        end
    end    
end

%% ����S21ʾ��ͼ
figure(1);
Lf=plot(deltaLg,abs(S21(:,2,79))); % ��kappaLg=1ʱ��S21ʾ��ͼ����ӦF3.26a�ϰ���ͼ������ֵ��չʾЧ������ѡ��0.78��
set(Lf,'Linewidth',2);
fsize=20;
%xlim([-25 25]);
ylim([0 9]);
xlabel('\delta L_g','FontSize',fsize);
ylabel('S21(\kappa L_g = 1)','FontSize',fsize);
set(gca,'FontSize',fsize); 

%% ���λ���(����ÿ��kappaLg(��m)�� ɨ�������S21���)��Ϊ�˹۲�Ѱ�� 
% figure(2);
% for b=6:1:length(m)
%     for a=34:1:length(g_scan)     
%         plot(deltaLg,abs(S21(:,b,a)));
%         ylim([0 50]);
%         text(-13.5,40,['��L_g=' num2str(kappaLg(b)) ' ' 'b=' num2str(b)],'FontSize',12);
%         text(-13.5,30,['(��g-��_i)L_g=' num2str(g_scan(a)) ' ' 'a=' num2str(a)],'FontSize',12);
%         pause;
%     end
% end
%% (������һ�ڳ���)�۲�S21����¼��ֵ��
% ����Ѱ�����ֻ�Ǽ�¼���Լ��۲�õ�����ֵ����������ɳ����ж��ĸ�ģʽ������ֵ����Ϊʵ�������е�����
% �۲�1-3��ģ�Ľ����¼��aa�����У�Ҳ�����·����ο�
%  ��0��ʾ��Ӧģʽû�м��䡿��ֻѡ������ģ��
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
 
%% ������������ֵ��������ֱ꣬�ӻ�ͼ

figure(2);
%���1����ͬkappaLg�������£���ͬģʽ��Ӧ����ֵ��
Ld=plot(deltaLg(cell2mat(locs)),g_scan(aa),'h');
hold on

%���ߣ���ͬkappaLg�������£�ĳ��ģʽ����ֵ�����ߣ�
for t=1:1:6
    plot(deltaLg(cell2mat(locs')),g_scan(aa'),'-k') ;
    hold on
end
set(Ld,'Linewidth',2);
xlabel('\delta L_g','FontSize',fsize);
ylabel('(��g_t_h - \alpha ) L_g','FontSize',fsize);
set(gca,'FontSize',fsize); 

toc;