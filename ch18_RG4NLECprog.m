function [X_list,lambda_k] = ch18_RG4NLECprog(f, c, x0, epsilon)
    % reduced-gradient algorithm for general nonlinear equality constrainted problem 
    % min f(x) s.t. ci(x)=0, i=1,...,m 
    
    % can refer to :
    % https://github.com/ishank011/grgdescent/blob/master/reduced_gradient.py
    % https://github.com/el-mouatasim/Reduced-Gradient-Algorithm

    % initialization
    X_list = x0;
    lambda_k = ones(size(c));

    % Lagrange function
    lambda = sym('lambda',[1,length(c)]).'; 
    L = f - lambda.'*c;              

    % Lagrange gradient and hessian
    dL = jacobian(L, symvar(f)).';    % L'[vector]
    d2L = hessian(L, symvar(f));      % L''[matrix]
    % 用对称正定矩阵近似 L''[matrix]
    Bk = eye(length(symvar(f)));      % 初始化Bk---Bk近似d2L

    while true
        % calculate xk, gk
        xk = X_list(:,end);
        gk = eval(subs(dL, [symvar(f), lambda.'], [xk.',lambda_k.']));      % 计算当前迭代点的一阶信息
        % calculate Y Z
        Aeq = eval(subs(jacobian(c),symvar(c),xk.'));  % each row is (nabla c_{i}).'
        [Q,R]= qr(Aeq.');
        [l,~]=size(R.');                 % l: size of R.'
        Y = Q(:,1:l); Z = Q(:,l+1:end);  % range/null space

        % 因为只在Z空间进行更新, ybar不需要更新
        zbar = -(Z.'*Bk*Z)\(Z.'*gk);
        pk = Z*zbar;                                % 计算更新方向
        lambda_k = (Y.'*Aeq.')\(Y.'*gk+Y.'*Bk*pk);
        try
            sstar = ch8_ModifiedArmijoLS(subs(L,lambda,lambda_k.'), xk, pk);    % 求解更新步长s 
        catch    % 可能会出现zero division的情形
            break;
        end
        pk_hat = -Aeq\eval(subs(c,symvar(c),(xk + sstar.*pk).'));       % 计算回溯步
        xkp1 = xk + sstar.*pk + pk_hat;             % 更新xkp1
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
        % 一阶收敛条件
        if ((norm(Z.'*gkp1)<epsilon)|(abs(xkp1-xk)<epsilon))  % 非线性约束问题的既约梯度无法收敛
            break;
        end  
    end
end