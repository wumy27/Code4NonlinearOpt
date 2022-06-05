clear all;close all;
%% 定义函数
% 用matlabFunction()可以将符号函数转化为匿名函数
syms X Y
f = -2*cos(X) - cos(Y) - cos(X-Y);
% X0 = [-5;-2];   % no bug
X0 = [-3;-2];   % [fixed]bug[最终到达一个鞍点]--因为没有保证deltak.'*gammak>0
epsilon = 1e-6;

% 因为是无约束优化问题，所以没有lb和ub

%% Quasi Newton Method - versionA
formula = 'BFGS';  % 'DFP/BFGS/SR1'
X_list = ch10_QuasiNewtonMethod(f, X0, epsilon, formula);
x_list = X_list(1,:);y_list = X_list(2,:);  
ShowGD(['拟牛顿法-',formula],x_list,y_list);   % 三种Hk更新公式的迭代次数为 8/8/8 % 比牛顿法的迭代次数+2

%% Quasi Newton Method - versionB
formula = 'BFGS';  % 'DFP/BFGS/SR1'
X_list = ch10_QuasiNewtonMethodB(f, X0, epsilon, formula);
x_list = X_list(1,:);y_list = X_list(2,:);  
ShowGD(['拟牛顿法versionB-',formula],x_list,y_list);   % 三种Bk更新公式的迭代次数为 7//8 % 比牛顿法的迭代次数+2
