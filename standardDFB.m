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
r=(n2-n1)/(n1+n2);
t=sqrt(1-r^2);
m=round(kappaLg/(2*r)); % ���ǵ�ʵ�ʹ�դ������������
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
            T{c,b,a} = T_interface12*T_line2*T_interface21*T_line1; % �ظ���ԪT����
            
            Tg{c,b,a} = T{c,b,a}^m(b);
            S21(c,b,a) = 1/Tg{c,b,a}(1);
        end
    end    
end

%% ����S21ʾ��ͼ
figure(1);
Lf=plot(deltaLg,abs(S21(:,2,12))); % ��kappaLg=1ʱ��S21ʾ��ͼ����ӦF3.26a�ϰ���ͼ������ֵ��չʾЧ������ѡ��2.75��
set(Lf,'Linewidth',2);
fsize=20;
%xlim([-25 25]);
ylim([0 9]);
xlabel('\delta L_g','FontSize',fsize);
ylabel('S21(\kappa L_g = 1)','FontSize',fsize);
set(gca,'FontSize',fsize); 

%% ���λ���(����ÿ��kappaLg(��m)�� ɨ�������S21���)��Ϊ�˹۲�Ѱ�� 
% figure(2);
% for b=1:1:length(m)
%     for a=1:1:length(g_scan)     
%         plot(deltaLg,abs(S21(:,b,a)));
%         ylim([0 400]);
%         text(-13.5,385,['��L_g=' num2str(kappaLg(b)) ' ' 'b=' num2str(b)],'FontSize',12);
%         text(-13.5,360,['(��g-��_i)L_g=' num2str(g_scan(a)) ' ' 'a=' num2str(a)],'FontSize',12);
%         pause;
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
 
%% ���������ֵ��������꣬aaa[��ֵindex]��locs0[deltaLg index]����;ʧг���Ӹ��������У�m(1->6)
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
%���1����ͬkappaLg�������£���ͬģʽ��Ӧ����ֵ��
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

toc;