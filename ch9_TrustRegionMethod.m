function X_list = ch9_TrustRegionMethod(f, x0, epsilon)
    X_list = x0;
    % 设计 f'(x) & f''(x)
    df = jacobian(f,symvar(f));      % f的一阶导f'[vector]
    d2f = hessian(f,symvar(f));      % f的二阶导f''[matrix]
    while true
        xk = X_list(:,end);
        % 计算当前迭代点的一阶和二阶信息
        gk = eval(subs(df,symvar(df),xk.')).';   % =f'(xk) 要做一个转置
        Gk = eval(subs(d2f,symvar(d2f),xk.'));   % 对称矩阵
        
        %% 判断Gk是否正定
        try L = chol(Gk).';  % Matrix is symmetric positive definite.
            z = -L\gk;
            pk = L.'\z;
        catch ME            % Matrix is not symmetric positive definite
%             disp(['Matrix is not symmetric positive definite']);
            mu = 1/norm(-Gk).*rand(1,10);    % 选取10[超参数]个mu值
            pi_list = -mu.*gk+mu.^2.*(Gk*gk)-mu.^3.*(Gk^2*gk);  % 每一列是一个pi
            kappai_list = [];      % kappa是一个数值
            for i = 1:10  % 总共10个pi
                kappai_list = [kappai_list,(pi_list(:,i).'*Gk*pi_list(:,i))./(pi_list(:,i).'*pi_list(:,i))];
            end
            [~,minInd]=min(kappai_list);
            pk = pi_list(:,minInd);
        end
        
        % 求解更新步长s 
        %% 完美线搜索   % 迭代次数 4
        sstar = ch8_PerfectLineSearch(f, xk, pk); 
        
        %% weak line search
%         sstar = ch8_ArmijoLineSearch(f, xk, pk);   % 迭代次数 23
%         sstar = ch8_ModifiedArmijoLS(f, xk, pk);   % 迭代次数 5
        % 更新xkp1
        X_list = [X_list, xk + sstar.*pk];
        
        if norm(eval(subs(df,symvar(df),X_list(:,end).')))<epsilon  % =f'(X_list(:,end))
            break;
        end
    end
end