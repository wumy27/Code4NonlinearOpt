function X_list = ch7_SteepestDescent(f, x0, lb0, ub0, epsilon)
    syms s
    X_list = x0;
    % 设计函数f'(x)
    df = jacobian(f,symvar(f));      % f的一阶导函数f'[vector]
    while true
        xk = X_list(:,end);
        % 计算当前迭代点的负梯度
        pk = -eval(subs(df,symvar(df),xk.')).';   % =-f'(xk) 要做一个转置

        %% 求解更新步长s
        % ch2 一维优化算法  
        % 待最小化函数 
%         phi = matlabFunction(subs(f,symvar(f),[xk+s.*pk].'));   
%         % 计算步长s的容许上下界
%         lbANDub = [(lb0 - xk)./pk,(ub0 - xk)./pk];
%         ub = min(max(lbANDub,[],2));
%         lb = max(min(lbANDub,[],2));
%         s_list = ch2_BisectionMethod(phi, lb, ub, epsilon);    % 完美线搜索
%         sstar = s_list(end);
        
        % ch8 线搜索   
%         sstar = ch8_ArmijoLineSearch(f, xk);    % 弱线搜索
        sstar = ch8_ModifiedArmijoLS(f, xk);    % 弱线搜索[收敛速度与完美线搜索相差很小]
        
        % 更新xkp1
        X_list = [X_list, xk + sstar.*pk];
        
        if norm(eval(subs(df,symvar(df),X_list(:,end).')))<epsilon  % =f'(X_list(:,end))
            break;
        end
    end
end