clc;
clear all;
%%The only difference is that in the latter case, the lower plate pair is present. In this example, the beam voltage YB= 10 kV, the plate length of a single pair L,=11.32 mm, plate separation D=2.03 mm, location of the beam stop as measured from the top of the upper 
% plate pair is 25.9 mm. In both cases, conjugate blanking is assumed and, for simplicity, the same lens magnification of M rens=0.24 is used. The slight adjustment in 
% object plane position between the two cases is ignored. The 
% 100/o-90% rise/fall time is taken to be 1.0 ns, which corresponds to (a$) = ( l.O,l. 1) in both cases. The blankervoltage swing is VD=7.5 V and the blanker aperture diameter is &,=50 pm. The blanker sensitivity at the beam 
% stop is then 5.7 pm/V. In addition, for the doubledeflection blanker, the gap g=29.2 mm and the electron 
% drift length is ZD=41.615 mm. 

% 设计参数
% 定义电压变化函数的参量
a = 0.01644; %指数函数系数
b = 0.02042; %三角函数系数

% 定义电极几何参数
Vd =           7.5;%         200;       %极板之间电压为7.8V
plategap=               0.080;    %上下极板间长度20mm
platelength=   0.048;%         0.0125;% 极板自身长度11.32mm
platewidth=    0.002;%         0.0005;% 极板自身宽度2.03mm
aputurewidth=0.0005;
halfway=platelength+plategap/2;

% 定义电子飞行参数

m0=9.10938215e-31;    % 电子静止质量： 9.109 38215(45)×10 -31 kg
c0=299792458;         % 光速： 299 792 458 m / s
VB=50000   ;          % 加速电压：10Kv
e0=1.602e-19 ;        % 电子带电量：-1.602 × 10−19C
E0=e0/(2*m0*c0^2);    
X20=(1+2*E0*VB)^2/(1+E0*VB);
vz=((2*e0*VB)/(m0*X20))^0.5;     %飞行速度
Gamma0=1/(1-(vz/c0)^2);
t01=(platelength/vz)*1e9;          %飞行穿过极板的时间
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
    % Vfinal2 = sum(VflyingtimefunctionUnit)    % 求和总速度（删掉）
    % Vfinal3 = trapz(Tflyingtime, Aflyingtime) % Matlab方法计算

    % 离开极板时的横向位移
    X1Single = trapz(Tflyingtime-PINGYI1, Vflyingtimefunction);  

    % T2 = Q((j+t02):(j+t01+t02), :);
    % col21 = T2(:, 1);  
    % col22 = T2(:, 2);  
    % V2 = trapz(col21, col22);

    X10((j+T0-1-PINGYI1)/T0)=X1Single;
    V10((j+T0-1-PINGYI1)/T0)=Vfinal1;
    Telapse((j+T0-1-PINGYI1)/T0)=(j+T0-1-PINGYI1);
    
    X101((j+T0-1-PINGYI1)/T0)=X1Single+((plategap/2)/vz)*Vfinal1;
end
Telapse0=Telapse/tt;

Z02=X10.*(vz./V10);
smooth_Z02 = smoothdata(Z02)';
standardZ02=smooth_Z02(1000*tt/T0);
final_Z02=smooth_Z02-standardZ02;
X02=final_Z02./(vz./V10)';

% Angletantheta=V10/vz;    % 电子离开极板飞行方向
% Vxiebian=(vz^2+V10.^2).^0.5;
% Angle2=V10./Vxiebian; 
% L11=X10./Angle2;
% Z02=(L11.^2-X10.^2).^0.5;% 这里有一部分修正因素
% smooth_Z02 = smooth(Z02);
% X22=Z02.*Angletantheta;
% Z0=X10./Angle1;   % 轴上位移
%  Z01=Z0-(platelength/2);
%  Xmid= Angle1 .*Z01;  
Xtotal1=[Telapse0,X10];
Vtotal1=[Telapse0,V10];



%% 后处理

hold on
%plot(Telapse0,X33,'R') % 绿色双极板
%plot(V10,Angle1,'b')
%plot(t,f,'b')
% plot(Telapse0,X101-2.41149e-5,'b')

plot(Telapse0,X02,'b')
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


% 停止计时器，并获取经过的时间
elapsedTime = toc;

% 显示计算时间
disp(['计算时间：', num2str(elapsedTime), ' 秒']);