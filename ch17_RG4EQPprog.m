function x_star = ch17_RG4EQPprog(G, h, Aeq, beq, l)
    % Reduced gradient for Equality Constraint Quadratic Programming
    % 目标函数f的设计不需考虑常数项，因为其不影响求得的最优解，仅影响最优函数值
    % 输入的G必须为正定矩阵，才能保证输出的解是最小值点
    % l: 为所需要被n-l个变元代替的变元数量[hyper parameter]
   
    [~,n] = size(Aeq);   % 变元(决策变量)个数

    Atilde = Aeq(:,1:n-l);
    Abar = Aeq(:,n-l+1:end);
    v = zeros(n,1); v(n-l+1:end) = Abar\beq;
    M = [eye(n-l);-Abar\Atilde];

    Gplus = M'*G*M;
    hplus = (v'*G*M+h'*M).'; 

    % new Quadratic Programming Problem donot have Constraints
    % means: xtilde = ch17_EQPprog(Gplus, hplus, [], [])
    % 变元个数可变的符号函数的创建: https://www.ilovematlab.cn/thread-296226-1-1.html
    X = sym('x',[1,n-l]).';
    f = 0.5.*X.'*Gplus*X+hplus.'*X;

    % 无约束优化算法: 不能取全零作为初始点，因为其为一个稳定点
    X_list = ch11_TruncatedNewton(f, rand(n-l,1), 1e-6); 
%     X_list = ch9_NewtonMethod(f, rand(n-l,1), 1e-6);
    xtilde = X_list(:,end);
    
    % 输出最优解
    % method-1
%     xbar = Abar\(beq-Atilde*xtilde);
%     x_star = [xtilde;xbar];
    % method-2
    x_star = v+M*xtilde;
    
end