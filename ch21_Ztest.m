clear all;

%% IQP4pdH √
% example in page227
X = sym('x',[1,2]).';
G = diag([2,6]); h = [0;0];
% 常数项c没有影响，可以考虑删除
Age = [1,5;0,1]; bge = [-1;0];   % 仅第一个为有效约束
% Age = [1,5]; bge = [-1];       % 仅第一个为有效约束
% Age = [1,5;0,1;3,2]; bge = [-1;0;0.5];  % 仅第一个为有效约束
x0 = [0;0];
[X_list,x_star,lambda_star]=ch21_IQP4pdH(G,h,Age,bge,x0);   
x_star                  % 真值为[0.1071;0.1786]
lambda_star

%% RG4IC
% example in Page252
X = sym('x',[1,2]).';
f = X(1)^2+3*X(2)^2;
Age = [1 5;5 -1]; bge = [-1;-0.25];  % c = [X(1)+5*X(2)-1; 5*X(1)-X(2)-0.25];
x0 = [1.0;0.0];     % 初始解需要为可行解，且满足一个不等式约束的等号情况 
[X_list,x_star] = ch21_RG4IC(f,Age,bge,x0);  % [0.0;0.2]&[1.0;0.0] -> 真解[0.1071;0.1786]
x_star                                % 只需要进行两步迭代

%% RG4SB
% personal designed problem
X = sym('x',[1,2]).';
f = (X(1)-2)^4+3*(X(2)-1)^4;
lb = [1;0]; ub = [3;2];
x0 = [1;0];
[X_list,x_star] = ch21_RG4SB(f,lb,ub,x0);
x_star                                    % 42步迭代，真解为[2;1]

%% P-SUMT for inequality constraints
X = sym('x',[1,2]).';
f = X(1)^2+2*X(2)^2;
ceq = []; cin = [X(1)+X(2)-2];
x0 = [0.2;0.2];            % 初始解非可行解/不能为[0;0]-->此时函数值为0
[X_list,x_star] = ch21_P_SUMT4IC(f, ceq, cin, x0, 1e-6);   % 11步迭代
x_star                     % 真解为[1.3333;0.6667]

%% AL-SUMT for inequality constraints
X = sym('x',[1,2]).';
f = X(1)^2+2*X(2)^2;
ceq = []; cin = [X(1)+X(2)-2];
x0 = [0.2;0.2];            % 初始解非可行解/不能为[0;0]-->此时函数值为0
[X_list,x_star] = ch21_AL_SUMT4IC(f, ceq, cin, x0, 1e-6);   % 7步迭代
x_star                     % 真解为[1.3333;0.6667]

%% EIQP4pdH for SQP algorithm [bug]
% example in page227
X = sym('x',[1,2]).';
G = diag([2,6]); h = [0;0];
% 常数项c没有影响，可以考虑删除
Aeq = [1 1]; beq = [8/28];
Age = [1,5;0,1]; bge = [-1;0];
x0 = [0;0];
[X_list,x_star,lambda_star]=ch21_EIQP4pdH(G, h, Aeq, beq, Age, bge, x0) ; 
x_star                  % 真值为[0.1071;0.1786]
lambda_star

%% WHP-SQP for inequality constraints  [bug]
% personal designed problem [not converge]
X = sym('x',[1,2]).';
f = (X(1)-2)^2+3*(X(2)-1)^2;
ceq = [X(1)-2*X(2)];  cin = [X(1)+X(2)-2];
x0 = [0.8;0.4];
[X_list,x_star] = ch21_WHP_SQP4IC(f, ceq, cin, x0, 1e-6);
x_star                                    % 真解为[2;1]

% example2 [solved]
X = sym('x',[1,2]).';
f = X(1)^2+2*X(2)^2;
ceq = [X(1)-2*X(2)]; cin = [X(1)+X(2)-2];
x0 = [0.4;0.2];            % 初始解非可行解/不能为[0;0]-->此时函数值为0
[X_list,x_star] = ch21_WHP_SQP4IC(f, ceq, cin, x0, 1e-6);
x_star                     % 收敛到真解为[1.3333;0.6667]

%% AL-SQP for inequality constraints
% cannot be done

