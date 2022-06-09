clear all;close all;
%% 定义函数
syms X Y
f = -2*cos(X) - cos(Y) - cos(X-Y);
X0 = [-5;-2];
lb = [-8;-8];
ub = [2;2];
epsilon = 1e-6;

%% 绘制函数图像
% 用matlabFunction()可以将符号函数转化为匿名函数
figure;SurfCG(matlabFunction(f));hold on;

%% Steepest descent
X_list = ch7_SteepestDescent(f, X0, lb, ub, epsilon);
x_list = X_list(1,:);y_list = X_list(2,:);  
ShowGD('最速下降法',matlabFunction(f),x_list,y_list);    % 迭代次数 = 10/18/19次
% [两种弱线搜索方法对于不同初始点表现不同]
% [-5;-2] 18/19
% [-2;-2] 18/12  完美线搜索的迭代次数均为10

h = gca;
set(h,'ylim',[-3 3],'xlim',[-10 -2]);
hold off;legend;