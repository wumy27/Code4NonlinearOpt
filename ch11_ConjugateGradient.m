function X_list = ch11_ConjugateGradient(f, x0, epsilon)
    % f'(x)
    df = jacobian(f,symvar(f));      % f的一阶导f'[vector]
    n = length(symvar(f));           % 变元个数
    % initialization
    X_list = x0;
    xk = X_list(:,end);
    gk = eval(subs(df,symvar(df),xk.')).';
    pk = -gk;
    
    while true
        [~,k] = size(X_list);   % k为当前的迭代轮次
        if k>=2   % k = 1不存在x/g/pkp1
            xk = xkp1;
            gk = gkp1;
            pk = pkp1;
        end
        % 必须进行完美线搜索
        sstar = ch8_PerfectLineSearch(f,xk,pk);
        % 更新解与梯度
        xkp1 = xk + sstar.*pk; X_list = [X_list,xkp1];
        gkp1 = eval(subs(df,symvar(df),xkp1.')).';
        % 更新方向
        if mod(k-1,n)~=0    % 判断当前迭代轮次k是否为n的倍数
            % calculate beta & pkp1 from (11.1.2)
            % 两个计算beta的公式后，计算得到的xk有不同
            beta = gkp1.'*gkp1/(gk.'*gk);        % Fletcher-Reeves公式
%             beta = gkp1.'*(gkp1-gk)/(gk.'*gk);   % Polak-Ribere公式
            pkp1 = -gkp1 + beta.*pk;
        else
            pkp1 = -gkp1;
        end
        % 一阶收敛条件
        if norm(gkp1)<epsilon
            break;
        end
end