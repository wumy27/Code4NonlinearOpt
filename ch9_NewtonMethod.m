function X_list = ch9_NewtonMethod(f, x0, epsilon)
    X_list = x0;
    % 设计 f'(x) & f''(x)
    df = jacobian(f,symvar(f));      % f的一阶导f'[vector]
    d2f = hessian(f,symvar(f));      % f的二阶导f''[matrix]
    while true
        xk = X_list(:,end);
        % 计算当前迭代点的一阶和二阶信息
        gk = eval(subs(df,symvar(df),xk.')).';   % =f'(xk) 要做一个转置
        Gk = eval(subs(d2f,symvar(d2f),xk.'));    % 对称矩阵
        
        %% 判断Gk是否正定
        % method-1  22/5[ArmijoLineSearch/ModifiedArmijoLS]
%         if all(eig(Gk) > epsilon)
%             pk = -Gk\gk;     % 求解方程组[可能需要其他手段]
%         else     % 非正定，选择负梯度方向
%             pk = -gk;
%         end
        % method-2  22/5[ArmijoLineSearch/ModifiedArmijoLS]
        try L = chol(Gk).';  % Matrix is symmetric positive definite.
            z = -L\gk;
            pk = L.'\z;
        catch ME            % Matrix is not symmetric positive definite
            pk = -gk;
        end
        
        % 弱线搜索求解更新步长s 
%         sstar = ch8_ArmijoLineSearch(f, xk, pk);   
        sstar = ch8_ModifiedArmijoLS(f, xk, pk);  
        % 更新xkp1
        X_list = [X_list, xk + sstar.*pk];
        
        if norm(eval(subs(df,symvar(df),X_list(:,end).')))<epsilon  % =f'(X_list(:,end))
            break;
        end
    end
end