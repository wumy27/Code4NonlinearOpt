clear all;close all;
%% 定义线搜索函数
% 以下算法针对线搜索问题进行构建
% 用matlabFunction()可以将符号函数转化为匿名函数
syms X Y
f = -2*cos(X) - cos(Y) - cos(X-Y);
X0 = [-5;-2];

%% Line Search Algorithm
sstar = ch8_ArmijoLineSearch(f, X0)
sstar = ch8_ModifiedArmijoLS(f, X0)
