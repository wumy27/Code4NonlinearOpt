clear all;close all;
%% 定义函数
% 用matlabFunction()可以将符号函数转化为匿名函数
syms X Y
f = -2*cos(X) - cos(Y) - cos(X-Y);    % 目标函数
X0 = [-5;3];                % 初始迭代点
lb = [-8;-8]; ub = [2;2];   % 重新设置上下限区域
epsilon = 1e-6;

%% 坐标轮转法
X_list = ch5_UnivariateSearch(f, X0, lb, ub, epsilon);
x_list = X_list(1,:);y_list = X_list(2,:);
ShowGD('坐标轮转法',x_list,y_list);

%% Hooke-Jeeves算法: X_list的长度比坐标轮转法更少
X_list = ch5_HooJee(f, X0, lb, ub, epsilon);
x_list = X_list(1,:);y_list = X_list(2,:);
ShowGD('Hooke-Jeeve法',x_list,y_list);

%% Nelder-Mead Simplex算法【比较低效的算法】
% 很依赖定义域，如果初始域过大，算法最后可能找到一个离初始点更远的最优解
X_list = ch5_NMSimplex(f, lb, ub, epsilon);   % 初值是通过程序随机生成
x_list = X_list(1,:);y_list = X_list(2,:);    % 最终迭代次数为80次
ShowGD('Nelder-Mead Simplex法',x_list,y_list);

%% DIRECT
[X_list,P_list,m,t] = ch5_MultiDIRECT(f, lb, ub, -4, epsilon);  
% 如果不给定函数的最优值-4作为判别收敛条件，则收敛速度很慢
x_list = X_list(1,:);y_list = X_list(2,:);    
ShowGD('SIMPLEX法',x_list,y_list);
% 绘制被评估过的点位置
Px = P_list(1,:);Py = P_list(2,:);scatter3(Px,Py,-5.*ones(size(Px)),'bh');  
