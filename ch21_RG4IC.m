function [X_list,xkp1] = ch21_RG4IC(f,Age,bge,x0)
    % reduced-gradients for nonquadratic objective function and linear inequality constraints
    % min f(x) s.t. Age*x+bge>=0
    % x0需要为可行解，且满足一个不等式约束的等号情况
    % xkp1 == x_star

    % initialization
    epsilon = 1e-6;           % 直接设置误差限（数值算法一定有误差）
    X_list = x0;
    lambda = zeros(size(bge));
    Bk = eye(length(symvar(f)));   % estimate d2f
    df = jacobian(f,symvar(f)).';  

    while true 
        xk = X_list(:,end);
        gk = eval(subs(df,symvar(f),xk.'));

        % 选择有效集
        active_ind = find(abs(Age*xk+bge)<epsilon & lambda>=0).';
        
        % 计算列空间与零空间
        Aeq = Age(active_ind,:);   % each row is (nabla c_{i}).'
        [Q,R]= qr(Aeq.');
        [l,~]=size(R.');                 % l: size of R.'
        Y = Q(:,1:l); Z = Q(:,l+1:end);  % range/null space
        
        % 计算更新方向
        zbar = -(Z.'*Bk*Z)\(Z.'*gk);
        pk = Z*zbar;                

        lambda = (Y.'*Aeq.')\(Y.'*gk+Y.'*Bk*pk);
        try
            sstar = ch8_ModifiedArmijoLS(f, xk, pk);    % 求解更新步长s 
        catch    % 可能会出现zero division的情形
            break;
        end
        xkp1 = xk + sstar.*pk;             % 更新xkp1
        X_list = [X_list, xkp1];                    

        % 更新Bk 
        gkp1 = eval(subs(df, symvar(f), xkp1.'));
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
        if ((norm(Z.'*gkp1)<epsilon)|(abs(xkp1-xk)<epsilon))  % 非线性约束问题的既约梯度无法收敛
            break;
        end
    end
end
