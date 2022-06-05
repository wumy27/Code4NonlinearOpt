function X_list = ch11_TruncatedNewton(f, x0, epsilon)
    % initialization
    X_list = x0;
    C = 1.5;
    df = jacobian(f,symvar(f));      % f的一阶导f'[vector]
    d2f = hessian(f,symvar(f));      % f的二阶导f''[Matrix]
    k = 0;                           % 记录迭代数
    while true
        k = k + 1;
        xk = X_list(:,end);          % end = k
        gk = eval(subs(df,symvar(df),xk.')).';
        Gk = eval(subs(d2f,symvar(d2f),xk.'));

        vk = min([C*norm(gk),1/k]);  % 阈值会随着k增大而减小
        
        % 求取更新方向
        % method-1 完美更新方向
%         pk = -Gk\gk;   % 迭代次数 7
        % method-2 有一定误差容许的更新方向
        [pk,~] = ch11_CG4LinearEqts( Gk, gk, xk, vk);  % 迭代次数 13

        % 求取更新步长--必须使用完美线搜索
        sstar = ch8_PerfectLineSearch(f,xk,pk);  

        X_list = [X_list, xk+sstar.*pk];
        if norm(gk)<epsilon
            break;
        end
    end
end