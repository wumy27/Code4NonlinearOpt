function X_list = ch19_WHP_SQP(f, c, x0, epsilon)

    % initialization
    X_list = x0;
    lambda_k = ones(size(c));

    % function derivate
    Bk = eye(length(symvar(f)));
    df = jacobian(f,symvar(f)).';    % g(x)
    d2f = hessian(f,symvar(f));      % G(x)
    dc = jacobian(c,symvar(c)).';    % A
    d2c = hessian(c,symvar(c));   
        
    % Lagrange function
    lambda = sym('lambda',[1,length(c)]).'; 
    L = f - lambda.'*c;
    dL = jacobian(L, symvar(f)).';   % L'[vector]
    d2L = d2f - lambda.' * d2c;      % G_hat
    
    % 误差度量函数
    kappa = 1e5;      % 惩罚参数
    T = norm( df - lambda.' * dc) + kappa * norm( c ); % g = df

    while true
        xk = X_list(:,end);
        gk = eval(subs(dL, [symvar(f), lambda.'], [xk.',lambda_k.'])); 
 
        % 解eqp子问题
        % qp_f = 1/2.* (p.' * Bk * p) + eval(subs(df,symvar(f),xk.')).' * p;
        % qp_c = c + dc.'*p;
        G = Bk; h = eval(subs(df,symvar(f),xk.')); Aeq = eval(subs(dc,symvar(c),xk.')).'; beq = -eval(subs(c,symvar(c),xk.'));
        [pk,lambda_k] = ch17_EQPprog(G, h, Aeq, beq);

        % 线搜索
        sstar = ch8_ModifiedArmijoLS(subs(L,lambda,lambda_k.'), xk, pk);    % 求解更新步长s 
        xkp1 = xk + sstar.*pk;
        X_list = [X_list, xkp1];

        % 更新Bk
        gkp1 = eval(subs(dL, [symvar(f), lambda.'], [xkp1.',lambda_k.']));
        gammak = gkp1 - gk;   % 计算梯度差
        deltak = xkp1 - xk;   % 计算解间隔
        % 通过求解Bkp1*deltak = gammak得到一个正定矩阵Bkp1更新Bk
        if deltak.'*gammak > 0   % 当其为正数时，才能保证更新公式给出的Bk正定
            % BFGS formula
            Bk = Bk - (Bk*deltak*(gammak.')+gammak*deltak.'*Bk)/((deltak.')*gammak) + ( 1 + (deltak.'*Bk*deltak)./(deltak.'*gammak) ).*(gammak*gammak.')./(deltak.'*gammak);
        else
            % 书本P178 (17.3.11)给出用etak代替gammak的公式
            for theta = 1:-0.1:0
                etak = (1-theta).*gammak+theta.*Bk*deltak;
                if deltak.'*etak>0    % 只能保证略微大于0
                    Bk = Bk - (Bk*deltak*(etak.')+etak*deltak.'*Bk)/((deltak.')*etak) + ( 1 + (deltak.'*Bk*deltak)./(deltak.'*etak) ).*(etak*etak.')./(deltak.'*etak);
                    break;
                end
            end
        end
        
        % 最优性判定
        if (eval(subs(T,[symvar(f), lambda.'],[xkp1.',lambda_k.']))<epsilon)|(abs(xkp1-xk)<epsilon*1e2)  % 用判定公式可以收敛，加上迭代差更快收敛
            break;
        end
    end
end