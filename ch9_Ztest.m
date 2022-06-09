clear all;close all;
%% 定义函数
syms X Y
f = -2*cos(X) - cos(Y) - cos(X-Y);
X0 = [-5;-2];   % no bug
% X0 = [-3;-2];   % no bug
epsilon = 1e-6;
% 因为是无约束优化问题，所以没有lb和ub

%% 绘制函数图像
% 用matlabFunction()可以将符号函数转化为匿名函数
figure;SurfCG(matlabFunction(f));hold on;

%% Newton Method
X_list = ch9_NewtonMethod(f, X0, epsilon);
x_list = X_list(1,:);y_list = X_list(2,:);  
ShowGD('牛顿法',matlabFunction(f),x_list,y_list);      % 两种判定Gk正定的方法，迭代次数相同，在w/w线搜索下迭代次数分别为22/5次

%% Trust Region Method
X_list = ch9_TrustRegionMethod(f, X0, epsilon); % 使用完美线搜索效果更好，此时需要额外设定可行域的范围
x_list = X_list(1,:);y_list = X_list(2,:);  
ShowGD('信赖域法',matlabFunction(f),x_list,y_list);    % 选取p/w/w线搜索，迭代次数分别为4/23/5次  

%% Guass Newton Method
% 高斯牛顿法，针对的目标函数为 \sumf^{2}(x)形式
% 此时可以 降低 求解pk的方程组 的 计算复杂度

h = gca;
set(h,'ylim',[-3 3],'xlim',[-10 -2]);
hold off;legend;