function [x_star,lambda_star] = ch17_EQPprog(G, h, Aeq, beq)
    % Equality Constraint Quadratic Programming
    % form: Aeq*x = beq // not same as A*x+b=0 [beq = -b]
    % One Order Optimity: Solve [G,-Aeq';-Aeq,0]*[x_star;lambda_star]=[-h;-beq];
    % 目标函数f的设计不需考虑常数项，因为其不影响求得的最优解，仅影响最优函数值
    % 输入的G必须为正定矩阵，才能保证输出的解是最小值点

    % method-1
%     sol_star = [G,-Aeq.';-Aeq,0]\[-h;-beq];
%     x_star = sol_star(1:end-1);
%     lambda_star = sol_star(end);
    
    % method-2  % 计算复杂度更小
    lambda_star = (Aeq/G*Aeq.')\(Aeq/G*h+beq);
    x_star = G\(Aeq.'*lambda_star-h); 
    
    if norm(lambda_star) == inf  % 针对超定方程组
        lambda_star = lsqminnorm((Aeq/G*Aeq.'),(Aeq/G*h+beq));
        x_star = lsqminnorm(G,Aeq.'*lambda_star-h); 
    end

end