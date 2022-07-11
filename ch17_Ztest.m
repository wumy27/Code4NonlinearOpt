clear all;
close all;

%% test ShowGDp2.m
method_name = 'test';
x_list = [0.5:0.3:1.4];
y_list = [0.5:0.3:1.4];
ShowGDp2(method_name,x_list,y_list);

%% EQP problem        % 一步到位，无法用ShowGDp2绘图
% 使用拉格朗日函数的一阶最优性进行求解
% P171 example1
G = [6 12;12 28]; h = [-72; -160]; Aeq = [-1 0]; beq = [-5];
[x_star, lambda_star] = ch17_EQPprog(G, h, Aeq, beq)  % [5.0000 3.5714; -0.8571]

% P171 example2
G = [2 3;3 0]; h = [0; 0]; Aeq = [1 5]; beq = [1];
[x_star, lambda_star] = ch17_EQPprog(G, h, Aeq, beq)  % [-0.7500 0.3500; -0.4500]

%% Reduced gradients  % 一步到位，无法用ShowGDp2绘图
% P171 example2
G = [2 3;3 0]; h = [0; 0]; Aeq = [1 5]; beq = [1]; 
x_star = ch17_RG4EQPprog(G, h, Aeq, beq, 1) % [ -0.7500; 0.3500]

% P173 exercise2
G = [0.06,-0.02,0; -0.02,0.05,0; 0,0,0]; h = [0;0;0];
Aeq = [0.2,0.3,0.15;1,1,1]; beq = [0.2;1];
x_star = ch17_RG4EQPprog(G, h, Aeq, beq, 2)   % [0.1549, 0.2817, 0.5634]
% ground truth: x_star = quadprog(G,h,[],[],Aeq,beq)

%% Reduced gradients for EQP  % 一步到位，无法用ShowGDp2绘图
% P171 example2
G = [2 3;3 0]; h = [0; 0]; Aeq = [1 5]; beq = [1]; 
x_star = ch17_RG4EQPprogv2(G, h, Aeq, beq)  % [ -0.7500; 0.3500]

% P173 exercise2
G = [0.06,-0.02,0; -0.02,0.05,0; 0,0,0]; h = [0;0;0];
Aeq = [0.2,0.3,0.15;1,1,1]; beq = [0.2;1];
x_star = ch17_RG4EQPprogv2(G, h, Aeq, beq)  % [0.1549, 0.2817, 0.5634]

%% 一般线性约束问题的既约梯度法
% 二次型problem1
X = sym('x',[1,2]).';
G = [2 3;3 0]; h = [0; 0]; Aeq = [1 5]; beq = [1]; 
f = 0.5.*X.'*G*X+h.'*X;
X_list = ch17_RG4GLECprog(f, [0;1/5], 1e-6, Aeq, beq)    % 一步优化[初始点需要满足约束/[1;0]]

% 非二次型problem2
X = sym('x',[1,3]).';
f = sum((X-[1 2 3].').^8);
Aeq = [1 1 1]; beq = [6];
X_list = ch17_RG4GLECprog(f, [3;-2;5], 1e-6, Aeq, beq)   % 一步优化[初始点需要满足约束]

%% 一般线性约束问题的投影梯度法
% 二次型problem1
G = [2 3;3 0]; h = [0; 0]; Aeq = [1 5]; beq = [1]; 
X = sym('x',[1,2]).';
f = 0.5.*X.'*G*X+h.'*X;  
X_list = ch17_PG4GLECprog(f, [1;0], 1e-6, Aeq, beq)      % 一步优化[初始点需要满足约束]

% 非二次型problem2
X = sym('x',[1,3]).';
f = sum((X-[1 2 3].').^8);
Aeq = [1 1 1]; beq = [6];
X_list = ch17_PG4GLECprog(f, [3;-2;5], 1e-6, Aeq, beq)   % 一步优化[初始点需要满足约束]
