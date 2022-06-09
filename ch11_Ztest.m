clear all;close all;
%% 定义函数
syms X Y
f = -2*cos(X) - cos(Y) - cos(X-Y);
% X0 = [-5;-2]; 
X0 = [-3;-2]; 
epsilon = 1e-6;

% 因为是无约束优化问题，所以没有lb和ub

%% 绘制函数图像
% 用matlabFunction()可以将符号函数转化为匿名函数
figure;SurfCG(matlabFunction(f));hold on;

%% Conjugate Gradient Method
X_list = ch11_ConjugateGradient(f, X0, epsilon);
x_list = X_list(1,:);y_list = X_list(2,:);  
ShowGD('共轭梯度法',matlabFunction(f),x_list,y_list);      % 迭代次数 8

%% Truncated Newton Method
X_list = ch11_TruncatedNewton(f, X0, epsilon);
x_list = X_list(1,:);y_list = X_list(2,:);  
ShowGD('裁剪牛顿法',matlabFunction(f),x_list,y_list);      % 迭代次数 7[完美pk]/13[误差容许pk]

h = gca;
set(h,'ylim',[-3 3],'xlim',[-3 3]);
hold off;legend;
