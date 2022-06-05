function X_list = ch10_QuasiNewtonMethod(f, x0, epsilon, formula)
    % 该版本拟牛顿法是用Hk对Hessian的逆矩阵进行近似
    X_list = x0;
    % 设计 f'(x)[vector]
    df = jacobian(f,symvar(f));      % f的一阶导f'[vector]
    % 用对称正定矩阵近似 f''(x)^{-1}[matrix]
    Hk = eye(length(symvar(f)));     % 初始化Hk---Hk近似Gk^{-1}
    while true
        xk = X_list(:,end);
        % 计算当前迭代点的一阶信息
        gk = eval(subs(df,symvar(df),xk.')).';   % =f'(xk) 要做一个转置
        % 计算更新方向
        pk = -Hk*gk;          % Hk一定正定  
        % 求解更新步长s 
        sstar = ch8_ModifiedArmijoLS(f, xk, pk);
        % 更新xkp1
        xkp1 = xk + sstar.*pk; X_list = [X_list, xkp1];
        % 计算梯度差与解间隔
        gkp1 = eval(subs(df,symvar(df),xkp1.')).';
        gammak = gkp1 - gk;
        deltak = xkp1 - xk;
        
        %% 通过求解Hkp1*gammak = deltak得到一个正定矩阵Hkp1更新Hk
        if deltak.'*gammak > 0   % 当其为正数时，才能保证更新公式给出的Hk正定
            %[not work] Hk = deltak/gammak;   % 无法得到正定矩阵
            if strcmp(formula,'DFP')       % DFP formula   % 迭代次数 7
                Hk = Hk - (Hk*gammak*(gammak.')*Hk)/((gammak.')*Hk*gammak) + (deltak*deltak.')/((deltak.')*gammak);
            elseif strcmp(formula,'BFGS')  % BFGS formula  % 迭代次数 7
                Hk = Hk - (Hk*gammak*(deltak.')+deltak*gammak.'*Hk)/((deltak.')*gammak) + ( 1 + (gammak.'*Hk*gammak)./(deltak.'*gammak) ).*(deltak*deltak.')./(deltak.'*gammak);
            elseif strcmp(formula,'SR1')   % 一阶对称更新公式  % 迭代次数 8
                vk = deltak-Hk*gammak;
                Hk = Hk + vk*vk.'/(vk.'*gammak);
            else
                disp(['sorry, formula-',formula,' has not been implemented.']);
                return 
                % future code
            end
        else
            Hk = eye(length(symvar(f)));  % 随机给定一个正定矩阵更新Hk
        end
        
        %% 一阶收敛条件
        if norm(eval(subs(df,symvar(df),xkp1.')))<epsilon  % =f'(xkp1)
            break;
        end
    end
end