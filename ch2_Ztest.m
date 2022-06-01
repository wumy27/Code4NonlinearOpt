clear all;close all;
%% 定义函数
% 用matlabFunction()可以将符号函数转化为匿名函数
syms X Y 
f = -2*cos(X) - cos(Y) - cos(X-Y);      % f =@(X,Y) -2*cos(X) - cos(Y) - cos(X-Y);
fx1 = gradient(f,X);                    % fx1 =@(X,Y) 2*sin(X) + sin(X - Y);
fx2 = gradient(f,Y);                    % fx2 =@(X,Y) sin(Y) - sin(X - Y);
fx1x2 = gradient(fx1,Y);                % fx1x2 =@(X,Y) -cos(X - Y);
fx1x1 = gradient(fx1,X);                % fx1x1 =@(X,Y) 2*cos(X) + cos(X - Y);
fx2x2 = gradient(fx2,Y);                % fx2x2 =@(X,Y) cos(Y) + cos(X - Y);

%% 固定x2=0,对f(x1,0)用二分法
[a, b, epsilon] = deal(-3.54, 2.03, 1e-6);
fye0 =@(X) -2*cos(X) - 1 - cos(X);     % 设定 x2 = 0
xm_list = ch2_BisectionMethod(fye0, a, b, epsilon);
ShowGD('二分法',xm_list,zeros(size(xm_list)));

%% 固定x2=0,对f(x1,0)用割线法
[a, b, epsilon] = deal(-3, 2, 1e-6);
fye0 =@(X) -2*cos(X) - 1 - cos(X);        % 设定 x2 = 0
xkp2_list = ch2_SecantMethod(fye0, a, b, epsilon);
ShowGD('割线法',xkp2_list,zeros(size(xkp2_list)));

%%  固定x2=0,对f(x1,0)用牛顿法
[x0, epsilon] = deal(-2, 1e-6);
fye0 =@(X) -2*cos(X) - 1 - cos(X);     % 设定 x2 = 0
xkp1_list = ch2_NewtonMethod(fye0, x0, epsilon);
ShowGD('牛顿法',xkp1_list,zeros(size(xkp1_list)));

%% 固定x2=0,对f(x1,0)用改进牛顿法
[a, b, epsilon] = deal(-3, 2, 1e-6);
fye0 =@(X) -2*cos(X) - 1 - cos(X);     % 设定 x2 = 0
xkp1_list = ch2_ModifiedNewtonMethod(fye0, a, b, epsilon);
ShowGD('改进牛顿法',xkp1_list,zeros(size(xkp1_list)));
