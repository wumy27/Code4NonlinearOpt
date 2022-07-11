clear all;
close all;

%% 一般非线性约束问题的既约梯度法
% code describtion in Section 18.1[not clear]
% competed by BevanWu
% exercise in page188
X = sym('x',[1,2]).';
f = X(1)^3+X(2)^2;
c = [X(2)-X(1)^2-1];
[X_list, lambda_star] = ch18_RG4NLECprog(f, c, [1;2], 1e-6);    
                      % 初始解[1;2]满足约束结果会好看些，不满足约束初始解[0.2;0.2]所得结果要差些
X_list(:,end)                                                      % 收敛到真解为[0 1]

%% 惩罚函数法 P-SUMT
% problem1  [一步优化]
X = sym('x',[1,3]).';
f = sum((X-[1 2 3].').^8);
c = [[1 1 1]*X-6];                           % Aeq = [1 1 1]; beq = [6];
X_list = ch18_P_SUMT(f, c, [0;0;0], 1e-6);   % 初始解不满足约束
X_list(:,end)                                % 收敛到真解为[1 2 3];

% an example in page187
X = sym('x',[1,2]).';
f = X.'*[1 3/2;3/2 0]*X;                         % = X(1)^2+3*X(1)*X(2);
c = [[1 5]*X-1];  
X_list = ch18_P_SUMT(f, c, [0.2;0.2], 1e-6);     % 初始解不满足约束
X_list(:,end)                                    % 收敛到真解为[-0.75 0.35];

% exercise in page188
X = sym('x',[1,2]).';
f = X(1)^3+X(2)^2;
c = [X(2)-X(1)^2-1];  
X_list = ch18_P_SUMT(f, c, [0.2;0.2], 1e-6);     % 初始解不满足约束
X_list(:,end)                                    % 收敛到真解为[0 1];

%% 惩罚函数法 AL-SUMT
% problem1  [一步优化]
X = sym('x',[1,3]).';
f = sum((X-[1 2 3].').^8);
c = [[1 1 1]*X-6];                            % Aeq = [1 1 1]; beq = [6];
X_list = ch18_AL_SUMT(f, c, [0;0;0], 1e-6);   % 初始解不满足约束
X_list(:,end)                                 % 收敛到真解为[1 2 3];

% an example in page187
X = sym('x',[1,2]).';
f = X.'*[1 3/2;3/2 0]*X;                          % = X(1)^2+3*X(1)*X(2);
c = [[1 5]*X-1];  
X_list = ch18_AL_SUMT(f, c, [0.2;0.2], 1e-6);     % 初始解不满足约束
X_list(:,end)                                     % 收敛到真解为[-0.75 0.35];

% exercise in page188
X = sym('x',[1,2]).';
f = X(1)^3+X(2)^2;
c = [X(2)-X(1)^2-1];  
X_list = ch18_AL_SUMT(f, c, [0.2;0.2], 1e-6);     % 初始解不满足约束
X_list(:,end)                                     % 收敛到真解为[0 1];

%% 精确惩罚函数法
% problem1  [一步优化]
X = sym('x',[1,3]).';
f = sum((X-[1 2 3].').^8);
c = [[1 1 1]*X-6];                                 % Aeq = [1 1 1]; beq = [6];
[X_list, x_star] = ch18_EPF(f, c, [0;0;0], 1e-6);  % 初始解不满足约束
x_star                                             % 收敛到真解为[1 2 3];

% an example in page187
X = sym('x',[1,2]).';
f = X.'*[1 3/2;3/2 0]*X;                              % = X(1)^2+3*X(1)*X(2);
c = [[1 5]*X-1];  
[X_list, x_star] = ch18_EPF(f, c, [0.2;0.2], 1e-6);   % 初始解不满足约束
x_star                                                % 收敛到真解为[-0.75 0.35];

% exercise in page188
X = sym('x',[1,2]).';
f = X(1)^3+X(2)^2;
c = [X(2)-X(1)^2-1];
[X_list, x_star] = ch18_EPF(f, c, [0.2;0.2], 1e-6);   % 初始解不满足约束
x_star                                                % 收敛到真解为[0 1];

