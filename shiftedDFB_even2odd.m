%%���������Լ��ٵ����ӣ�������F 3.26������ȫ��ͬ��û�ҵ����жԲ�����ֱ��˵����
% ������kappaLg�ĸı���ͨ�����ӹ�դ����mʵ�ֵ�
% ��ֵ����ͨ�� ɨ�費ͬ���� ��Ӧ��S21��ã���˿��ǽ�����ȣ�����������ʱ���ϳ�����ϰ�еľ���Ϊ0.25��
% ����ɨ�跨Ҳ������Ѱ��������������ģʽ��ֵ

tic;
clear;
%��ز�������
deltaLg = (-15:0.05:15)'; %����F 3.26
lambda0 = 1.55; % �������񳤶�1.55um��ƹ�դ
kappaLg = [0.5 1 1.5 2 3 4]; %����F 3.26
n1=3.2515;n2=3.456; %��InGaAsP/InPΪ��
r=(n2-n1)/(n1+n2);
t=sqrt(1-r^2);
m=[9 17 25 33 49 67];
% grating_period = lambda0/(4*n1)+lambda0/(4*n2);  % ��դ����

% ������ɨ���� g_scan ����(��g-��)Lg���ο�F3.26��
g_scan = 0:0.25:8; 

% ���������
T_interface12 = Tmatrix_interface(r);
T_interface21 = Tmatrix_interface(-r);

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
            T1{c,b,a} = T_interface12*T_line2*T_interface21*T_line1; % ǰ����ظ���ԪT����
            T2{c,b,a} = T_line1*T_interface12*T_line2*T_interface21; % �����ظ���ԪT����
            T_shift{c,b,a} = T_interface12*Tmatrix_line(2*phi2)*T_interface21; %��/4������T����
            Tg{c,b,a} = (T1{c,b,a}^floor(m(b)/2))*T_shift{c,b,a}*(T2{c,b,a}^floor(m(b)/2));
            S21(c,b,a) = 1/Tg{c,b,a}(1);
        end
    end    
end

%% ����S21ʾ��ͼ
figure(1);
Lf=plot(deltaLg,abs(S21(:,2,10))); % ��kappaLg=1ʱ��S21ʾ��ͼ����ӦF3.26b�ϰ���ͼ��g_scan����ֵ��չʾЧ������ѡ��2.25��
set(Lf,'Linewidth',2);
ylim([0 13])
fsize=20;
xlabel('\delta L_g','FontSize',fsize);
ylabel('S21(\kappa L_g = 1)','FontSize',fsize);
set(gca,'FontSize',fsize); 

%% ���λ���(����ÿ��kappaLg(��m)�� ɨ�������S21���)��Ϊ�˹۲�Ѱ�� 
% figure(2);
% for b=1:1:length(m)
%     for a=1:1:length(g_scan)     
%         plot(deltaLg,abs(S21(:,b,a)));
%         ylim([0 400]);
%         text(-14,385,['��L_g=' num2str(kappaLg(b)) ' ' 'b=' num2str(b)],'FontSize',12);
%         text(-14,360,['(��g-��_i)L_g=' num2str(g_scan(a)) ' ' 'a=' num2str(a)],'FontSize',12);
%         pause;
%     end
% end
%% (������һ�ڳ���)�۲�S21����¼��ֵ��
% ����Ѱ�����ֻ�Ǽ�¼���Լ��۲�õ�����ֵ����������ɳ����ж��ĸ�ģʽ������ֵ����Ϊʵ�������е�����
% �۲�0-4��ģ�Ľ����¼��aa�����У�Ҳ�����·����ο�
%       b=1: 16 20 23 25 27
%       b=2: 12 16 21 23 25
%       b=3: 9  14 19 21 24
%       b=4: 6  12 17 19 22
%       b=5: 4  9  14 16 19
%       b=6: 2  7  11 13 16
aa = zeros(5,6);
aa(:,1)=[18 22 26 28 30];aa(:,2)=[12 17 22 24 26];aa(:,3)=[9  14 19 21 24];
aa(:,4)=[6  12 17 19 22];aa(:,5)=[4  9  14 16 19];aa(:,6)=[2  7  11 14 16];
for b=1:1:6
    for a=1:1:5
        maxS21=max(abs(S21(:,b,aa(a,b))));
        [locs{a,b}]=find(maxS21==abs(S21(:,b,aa(a,b))));
    end
end
%% 
% ����ģʽ���긳ֵ�����
locs0 = ones(9,6);aaa=locs0;
for b=1:1:6 % ��д���ѭ��û���뵽�ð취��Ԫ����ֵ��ȥ
    for a=1:1:4
        locs0(a,b)=length(deltaLg)-locs{6-a,b};
        locs0(a+5,b)=locs{a+1,b};
    end
end
locs0(5,:)=locs{1,:};
aaa(5:9,:)=aa(:,:);
aaa(1:4,:)=aa(5:-1:2,:);
figure(2);

%��㣨��ͬkappaLg�������£���ͬģʽ��Ӧ����ֵ��
deltaLg0=sort(deltaLg(locs0)); %Ϊ��֮������߶�����
Ld=plot(deltaLg0',g_scan(aaa'),'h');
hold on

%���ߣ���ͬkappaLg�������£�ĳ��ģʽ����ֵ�����ߣ�
for t=1:1:6
    plot(deltaLg0',g_scan(aaa'),'-k') ;
    hold on
end
set(Ld,'Linewidth',2);
fsize = 20;
xlabel('\delta L_g','FontSize',fsize);
ylabel('(��g_t_h - \alpha ) L_g','FontSize',fsize);
set(gca,'FontSize',fsize); 

toc;