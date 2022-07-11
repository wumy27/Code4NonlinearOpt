function [X_list,lambda_star] = ch17_RG4GLECprog(f, x0, epsilon, Aeq, beq)
    % reduced-gradient algorithm for general linear equality constrainted problem 
    % min f(x) s.t. Aeq*x=beq

    % initialization
    X_list = x0;
    lambda_star = zeros(size(beq));

    % 设计 f'(x)[vector]
    df = jacobian(f,symvar(f));      % f的一阶导f'[vector]
    % 用对称正定矩阵近似 f''(x)[matrix]
    Bk = eye(length(symvar(f)));     % 初始化Bk---Bk近似Gk
    
    % calculate Y Z
    [Q,R]= qr(Aeq.'); L=R.';   % Q=Q;
    [l,~]=size(L);                   % l: size of L
    Y = Q(:,1:l); Z = Q(:,l+1:end);  % range/null space

    while true
        % 因为只在Z空间进行更新,故不许用到beq的信息,即ybar不需要更新
        xk = X_list(:,end);
        gk = eval(subs(df,symvar(df),xk.')).';      % 计算当前迭代点的一阶信息
        Z
        Bk
        zbar = -(Z.'*Bk*Z)\(Z.'*gk);
        pk = Z*zbar;                                % 计算更新方向
        lambda_star = (Y.'*Aeq.')\(Y.'*gk+Y.'*Bk*pk);
        sstar = ch8_ModifiedArmijoLS(f, xk, pk);    % 求解更新步长s 
        xkp1 = xk + sstar.*pk; X_list = [X_list, xkp1];   % 更新xkp1

        % 更新Bk 
        gkp1 = eval(subs(df,symvar(df),xkp1.')).';
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
        % 一阶收敛条件
        if norm(Z.'*gkp1)<epsilon  % 既约梯度
            break;
        end  
    end
end