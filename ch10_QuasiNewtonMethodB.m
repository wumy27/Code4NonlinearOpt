function X_list = ch10_QuasiNewtonMethodB(f, x0, epsilon, formula)
    % 该版本拟牛顿法是用Bk对Hessian矩阵进行近似
    X_list = x0;
    % 设计 f'(x)[vector]
    df = jacobian(f,symvar(f));      % f的一阶导f'[vector]
    % 用对称正定矩阵近似 f''(x)[matrix]
    Bk = eye(length(symvar(f)));     % 初始化Bk---Bk近似Gk
    while true
        xk = X_list(:,end);
        % 计算当前迭代点的一阶信息
        gk = eval(subs(df,symvar(df),xk.')).';   % =f'(xk) 要做一个转置
        % 计算更新方向
        pk = -Bk\gk;          % Bk一定正定  
        % 求解更新步长s 
        sstar = ch8_ModifiedArmijoLS(f, xk, pk);
        % 更新xkp1
        xkp1 = xk + sstar.*pk; X_list = [X_list, xkp1];
        % 计算梯度差与解间隔
        gkp1 = eval(subs(df,symvar(df),xkp1.')).';
        gammak = gkp1 - gk;
        deltak = xkp1 - xk;
        
        %% 通过求解Bkp1*deltak = gammak得到一个正定矩阵Bkp1更新Bk
        if deltak.'*gammak > 0   % 当其为正数时，才能保证更新公式给出的Bk正定
            %[not work] Bk = gammak/deltak;   % 无法得到正定矩阵
            if strcmp(formula,'DFP')       % DFP formula   % 迭代次数 7
                Bk = Bk - (Bk*deltak*(deltak.')*Bk)/((deltak.')*Bk*deltak) + (gammak*gammak.')/((deltak.')*gammak);
            elseif strcmp(formula,'BFGS')  % BFGS formula  % 迭代次数 7
                Bk = Bk - (Bk*deltak*(gammak.')+gammak*deltak.'*Bk)/((deltak.')*gammak) + ( 1 + (deltak.'*Bk*deltak)./(deltak.'*gammak) ).*(gammak*gammak.')./(deltak.'*gammak);
            elseif strcmp(formula,'SR1')   % 一阶对称更新公式  % 迭代次数 8
                wk = gammak-Bk*deltak;
                Bk = Bk + wk*wk.'/(wk.'*deltak);
            else
                disp(['sorry, formula-',formula,' has not been implemented.']);
                return 
                % future code
            end
        else
            Bk = eye(length(symvar(f)));  % 随机给定一个正定矩阵更新Bk
        end
        
        %% 一阶收敛条件
        if norm(eval(subs(df,symvar(df),xkp1.')))<epsilon  % =f'(xkp1)
            break;
        end
    end
end