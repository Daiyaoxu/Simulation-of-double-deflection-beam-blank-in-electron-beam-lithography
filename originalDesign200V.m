clc;
clear all;


% 设计参数
% 定义电压变化函数的参量
a = 0.32; %指数函数系数
b = 0.32; %三角函数系数

% 定义电极几何参数
Vd = 200;       %极板之间电压为7.8V
plategap= 0.020;    %上下极板间长度20mm
platelength=   0.015;%         0.0125;% 极板自身长度11.32mm
platewidth=    0.0005;%         0.0005;% 极板自身宽度2.03mm
aputurewidth=0.0005;
halfway=platelength+plategap/2;

% 定义电子飞行参数

m0=9.10938215e-31;    % 电子静止质量： 9.109 38215(45)×10 -31 kg
c0=299792458;         % 光速： 299 792 458 m / s
VB=30000   ;          % 加速电压：10Kv
e0=1.602e-19 ;        % 电子带电量：-1.602 × 10−19C
E0=e0/(2*m0*c0^2);    
X20=(1+2*E0*VB)^2/(1+E0*VB);
vz=((2*e0*VB)/(m0*X20))^0.5     %飞行速度
Gamma0=1/(1-(vz/c0)^2);
t01=(platelength/vz)*1e9         %飞行穿过极板的时间
t03=(plategap/vz)*1e9 ;         %飞行穿过极板间距的时间
m00 = m0 / sqrt(1 - (vz^2 / c0^2));% 电子修正质量 
% （这里有个问题就是横向加速度没有在程序里进行修正！！！）

                                                                                                                                                                              
% 定义电极最大加速度
aMAX = ((Vd/platewidth)*e0)/m00;

% 定义时间数组
t001=round(t01, 3);
tt=1e3;        %计算精度
t = linspace(0, 2000, 2000*tt);
indexflyingtime = round(t01, 3)*tt;  %飞过极板时间索引
indexflyinggaptime = round(t03, 3)*tt;  %飞过极板间空缺的时间索引



% 定义加速度函数数组
f = zeros(size(t));  % 创建一个与 t 大小相同的全零数组 

% 设置秒表查看计算时间
tic

% 定义时间范围
data1 = zeros(length(t), 1);  % 时间预分配数组1
data2 = zeros(length(t), 1);  % 加速度预分配数组2
Vtotaldata = zeros(99, 1);  % 预分配数组2

% 设置计算函数值
for i = 1:length(t)
    if t(i) <= 0
        f(i) = 0;
    elseif 0 < t(i) && t(i) <= 300
        f(i) = aMAX*(1+0.2) ;
    elseif 300 < t(i) && t(i) <= 1300
        f(i) = aMAX * ((exp(-a * (-300+t(i))) * cos(b * (-300+t(i))))+0.2);

    elseif 1300 < t(i) && t(i) <= 2000
        f(i) = aMAX * ((1 - exp(-a * (-1300+t(i))) * cos(b * (-1300+t(i))))+0.2);        
    end
    data1(i) = t(i);
    data2(i) = f(i);
end

% 真实时间数组
A=[data1 data2];% 左边时间 右边加速度
column1 = A(:, 1);
multipliedColumn1 = column1 * 1e-9;  % 假设要乘以2
A(:, 1) = multipliedColumn1;

% 脉冲弛豫时间

PINGYI1=0.02*tt;
PINGYI2=0*tt;


% 计算束流空间变化
T0=0.1*tt;  % 计算束流密度密度
Matrixsize1=length(A)/T0;

%% 上极板
X10=zeros(Matrixsize1:1);
X101=zeros(Matrixsize1:1);

V10=zeros(Matrixsize1:1);
Telapse=zeros(Matrixsize1:1);


for j = 1+PINGYI1:T0:1950*tt+PINGYI1

    T1 = A(j:(j+indexflyingtime), :);
    Tflyingtime = T1(:, 1); 
    Aflyingtime = T1(:, 2); 
    VflyingtimefunctionUnit=zeros(length(Tflyingtime)-1,1); 

    % 对于极板间飞行速度的单位变化量
    for j1=1:(length(Tflyingtime)-1)

        A1=Aflyingtime(j1);
        B1=Aflyingtime(j1+1);
        UnitV=(A1+B1)*(1e-9/tt)/2;
        VflyingtimefunctionUnit(j1)=UnitV;
    end

    AA00=[0]; % 修正数组
    VflyingtimefunctionUnit0 = [AA00;VflyingtimefunctionUnit]; % 每个单位时间极板内电子横向速度
    Vflyingtimefunction=zeros(length(VflyingtimefunctionUnit0),1); % 极板电子横向速度函数

    % 对于极板间飞行横向速度的函数
    for k1=1:(length(VflyingtimefunctionUnit0))

        subArray1 = VflyingtimefunctionUnit0(1:k1);
        functionV11=sum(subArray1);
        Vflyingtimefunction(k1)=functionV11;
    end

    Vfinal1= Vflyingtimefunction(end);         % 求解器计算


    % 离开极板时的横向位移
    X1Single = trapz(Tflyingtime-PINGYI1, Vflyingtimefunction);  

    X10((j+T0-1-PINGYI1)/T0)=X1Single;
    V10((j+T0-1-PINGYI1)/T0)=Vfinal1;
    Telapse((j+T0-1-PINGYI1)/T0)=(j+T0-1-PINGYI1);
    
    X101((j+T0-1-PINGYI1)/T0)=X1Single+((plategap/2)/vz)*Vfinal1;
end
Telapse0=Telapse/tt;

Z02=X10.*(vz./V10);
smooth_Z02 = smoothdata(Z02)';
standardZ02=smooth_Z02(1000*tt/T0)
final_Z02=smooth_Z02-standardZ02;
X02=final_Z02./(vz./V10)';


Xtotal1=[Telapse0,X10];
Vtotal1=[Telapse0,V10];

%% 下极板

X20=zeros(Matrixsize1:1);
V20=zeros(Matrixsize1:1);
Telapse2=zeros(Matrixsize1:1);


for j2 = 1+PINGYI2:T0:1450*tt+PINGYI2

    T1 = A(j2:(j2+indexflyingtime),:);
    Tflyingtime = T1(:, 1); 
    Aflyingtime = T1(:, 2); 
    VflyingtimefunctionUnit2=zeros(length(Tflyingtime)-1,1); 
    Vupperplate=V10((j2-PINGYI2+T0-1)/T0);
    Xupperplate=X10((j2-PINGYI2+T0-1)/T0);


    % 对于极板间飞行速度的单位变化量
    for j3=1:(length(Tflyingtime)-1)

        A2=Aflyingtime(j3);
        B2=Aflyingtime(j3+1);
        UnitV=(A2+B2)*(1e-9/tt)/2;
        VflyingtimefunctionUnit2(j3)=UnitV;
    end

    AA00=[0]; % 修正数组
    VflyingtimefunctionUnit02 = [AA00;VflyingtimefunctionUnit2]; % 每个单位时间极板内电子横向速度
    Vflyingtimefunction2=zeros(length(VflyingtimefunctionUnit02),1); % 极板电子横向速度函数

    % 对于极板间飞行横向速度的函数
    for k2=1:(length(VflyingtimefunctionUnit02))

        subArray2 = VflyingtimefunctionUnit02(1:k2);
        functionV12=sum(subArray2);
        Vflyingtimefunction2(k2)=functionV12++Vupperplate;
    end

    Vfinal2= Vflyingtimefunction2(end);         % 求解器计算
    % Vfinal2 = sum(VflyingtimefunctionUnit)    % 求和总速度（删掉）
    % Vfinal3 = trapz(Tflyingtime, Aflyingtime) % Matlab方法计算

    % 离开极板时的横向位移
    Tgap=plategap/vz;
    Xgap=Vupperplate*Tgap;
    X2Single = trapz(Tflyingtime-PINGYI2, Vflyingtimefunction2)+Xupperplate+Xgap;  

    % T2 = Q((j+t02):(j+t01+t02), :);
    % col21 = T2(:, 1);  
    % col22 = T2(:, 2);  
    % V2 = trapz(col21, col22);

    X20((j2-PINGYI2+T0-1)/T0)=X2Single;
    V20((j2-PINGYI2+T0-1)/T0)=Vfinal2;
    Telapse2((j2-PINGYI2+T0-1)/T0)=(j2-PINGYI2+T0-1);
end


Telapse0=Telapse2/tt;

Z03=X20.*(vz./V20);
smooth_Z03 = smoothdata(Z03)';
standardZ03=smooth_Z03(5000);
final_Z03=smooth_Z03-standardZ03;
X03=final_Z03./(vz./V20)';

% Angle1=V20/vz;    % 电子离开极板飞行方向
% Vangle=(vz^2+V20.^2).^0.5;
% Angle2=V20./Vangle; 
% L20=X20./Angle2;
% Z03=(L20.^2-X20.^2).^0.5-0.0133986; %w0.01132+
% Z03=Z03-0.012541;
% X33=Z03.*Angle1;

% Z0=X10./Angle1;   % 轴上位移
%  Z01=Z0-(platelength/2);
%  Xmid= Angle1 .*Z01;  
Xtotal1=[Telapse0,X20,'b'];
Vtotal1=[Telapse0,V20,'r'];


%% 后处理

hold on
%plot(Telapse0,X33,'R') % 绿色双极板
%plot(V10,Angle1,'b')
%plot(t,f,'b')
% plot(Telapse0,X101-X101(1000*tt/T0),'b')
% plot(Telapse0,smooth_Z02,'b')
plot(Telapse0,X10-X10(1000*tt/T0),'b')
%plot(Telapse0,X02-X02(1000*tt/T0),'b')
ylabel('a(t)')
%plot(Telapse0,X03,'r')
% ylim([-4e-7, 4e-7])
% smooth_X22 = smooth(X22); % 应用平滑方法，可以尝试不同的平滑算法
% plot(Telapse0,smooth_X22,'r') % 红色单极板


hold off


xlabel('时间')
ylabel('s(t)')

title('单偏转Beam Blank 电子束移动误差')
grid on





%% Dose 计算

% 定义剂量高斯分布的均值和标准差
FWHM=10;                           % Full Width at Half Maximum (FWHM) 50NM
sigma = FWHM/(2*(2*log(2))^0.5);   % X方向标准差
M=8;                               

% 锯齿波

    % 绘图尺寸 束斑间距 曝光10 个点
    spotgap= 10;            % 束斑间距10纳米
    spotnum=48.5 ;            % 束斑数量
    t02 = 20;               % 每个束斑曝光时间50ns
    xMax=spotnum*spotgap;
    alpha=deg2rad(0);
    Q1=sin(alpha);
    Q2=cos(alpha);

    % 计算量
    t0totoal=spotnum*t02;         %总曝光时间
    Numelectron=t0totoal*tt/T0   %总电子数量
    Tjuchi=t0totoal/spotnum;       %周期
    Ajuchi=spotgap;               %幅值
    totalsize=100+spotgap*spotnum+100 ;       % 绘图尺寸200纳米

    % 设置计算剂量的区域范围和步长
    Q_x_min = 0;           % X方向区域最小值
    Q_x_max = totalsize;   % X方向区域最大值
    Q_y_min = 0;           % Y方向区域最小值
    Q_y_max = totalsize;   % Y方向区域最大值
    step_size = 1;       % 步长0.1 nm（计算精度）

    % 找到束流开启的状态
%     targetvoltage = aputurewidth/2;
%     diff_values = abs(X10 - targetvoltage);
%     index_closest = find(diff_values == min(diff_values), 1);
%     indicesofbeamon = find(X10 == )
%     Topen=Telapse0(indicesofbeamon)

    TBlank=346.7; %开启时间

    % 指定取出的的起始位置和结束位置
    start_index = TBlank*tt/T0
    end_index = (TBlank+t0totoal)*tt/T0

    % 使用冒号运算符生成一段数据
    deltaX10=0;    %偏移量归零
    X101=X02-deltaX10;
    XBeamopen = X101(start_index:end_index);
    XBeamopenlength=length(XBeamopen);


    % 定义数组
    tarray = 1:Numelectron;
    tarray1=tarray/(tt/T0);
    xarry =(tarray*xMax)/((t0totoal)*(tt/T0));

    % 
    Dose_profile=zeros(totalsize/step_size+1, totalsize/step_size+1);




for n=1:1:Numelectron % spot 分段
        erro_x =0;     
        A5=n/Numelectron *100;
         disp([num2str(A5), '%']);

        erro_y =((round(XBeamopen(n),9)) /M)*1e9 ;  %(n-1)*50+n2*tt/T0       % Y方向误差

        % 生成剂量的空间网格，并计算每个位置上的剂量
        [Q_x, Q_y] = meshgrid(Q_x_min:step_size:Q_x_max, Q_y_min:step_size:Q_y_max);  

        dose_values1 =  exp(-0.5 * ((Q_x-(totalsize/2)+ xarry(n)-200).^2 / sigma^2 + ((Q_y-(totalsize/2) -(erro_y)).^2) / sigma^2)); %-5*spotgap  -(erro_y)*1e9  -(erro_y)
        %dose_values2 =  exp(-0.5 * ((Q_x-(totalsize/2)+xarry(n)/t02).^2 / sigma^2 + ((Q_y-(totalsize/2) ).^2) / sigma^2)); %-5*spotgap  -(erro_y)*1e9  -(erro_y)

          Dose_profile=Dose_profile+dose_values1;%+dose_values2;


end

% 
% 

%% 后处理

% 绘制剂量分布曲面
figure;
surf(Q_x,Q_y, Dose_profile)

% 假设您希望颜色范围在 [min_val, max_val] 之间
min_val = 0;    % 替换为您期望的最小值
max_val = 120;   % 替换为您期望的最大值
caxis([min_val, max_val]);

% 添加颜色条
%colorbar;


shading flat % 消除栅格线（语句一）
xlabel('位置 Q_x')
ylabel('位置 Q_y')
zlabel('剂量')
title('二维高斯分布的剂量分布')



% 
% 停止计时器，并获取经过的时间
elapsedTime = toc;

% 显示计算时间
disp(['计算时间：', num2str(elapsedTime), ' 秒']);