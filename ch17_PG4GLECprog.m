function X_list = ch17_PG4GLECprog(f, x0, epsilon, Aeq, beq)
    % projected-gradient algorithm for linear equality constraints 
    % min f(x) s.t. Aeq*x=beq
    % todo: (17.5.2)公式计算无效

    % initialization
    X_list = x0;
    n = length(x0);           % 解的维度

    % 设计 f'(x)[vector]
    df = jacobian(f,symvar(f));      % f的一阶导f'[vector]    

    while true
        % 因为只在Z空间进行更新,故不许用到beq的信息,即ybar不需要更新
        xk = X_list(:,end);
        gk = eval(subs(df,symvar(df),xk.')).';      % 计算当前迭代点的一阶信息

        % 投影矩阵1
        P = eye(n)-Aeq.'/(Aeq*Aeq.')*Aeq;               
%         % 投影矩阵2【需要求取f的Hessian近似矩阵/正定】
%         % 【当Bk==eye时，可以完成优化，但如果用拟牛顿法近似Bk，无法收敛】
%         Bk = QuasiNewton();
%         P = Aeq/Bk*Aeq.';
%         P = eye(n)-Aeq.'/P*Aeq;
%         P = Bk\P/Bk;   
        % 计算更新方向 
        pk = -P*gk;

        % 一阶收敛条件[以往版本放在最后，此函数若放最后会不收敛]
        if norm(pk)<epsilon  % 投影梯度
            break;
        end

        sstar = ch8_ModifiedArmijoLS(f, xk, pk);          % 求解更新步长s 
        xkp1 = xk + sstar.*pk; X_list = [X_list, xkp1];   % 更新xkp1
    end
end