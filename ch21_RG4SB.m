function [X_list,xkp1] = ch21_RG4SB(f,lb,ub,x0)
    % reduced-gradient algorithm for simple bounds.
    % min f(x) s.t. lb<=x<=ub
    % x0 should be feasible
    % xkp1 == x_star

    % 针对f的形式不同，所得到最终解的质量也不一样，（可能是因为Bk逼近的原因）
    % 减小epsilon能够提高解的质量

    % initialization
    epsilon = 1e-12;           % 直接设置误差限（数值算法一定有误差）
    X_list = x0;
    Bk = eye(length(symvar(f)));   % estimate d2f
    df = jacobian(f,symvar(f)).';  

    while true 
        xk = X_list(:,end);
        gk = eval(subs(df,symvar(f),xk.'));
        
        % 计算零空间
        Z = eye(length(symvar(f)));
%         active_ind = [find(abs(xk-lb)<epsilon & gk>0).',...
%             find(abs(xk-ub)<epsilon & gk<0).'];  % 21.2.1
        active_ind = [find(abs(xk-lb)<epsilon & gk>0).',...
            find(abs(xk-ub)<epsilon & gk<0).'];  % 21.2.1
        inactive_ind = [1:length(ub)]; inactive_ind(active_ind)=[]; 
        Z(:,active_ind)=[];
        
        % 计算更新方向
        zbar = -(Z.'*Bk*Z)\(Z.'*gk);
        pk = Z*zbar;                
        
        % 求解更新步长s 
        try
            sstar = ch8_ModifiedArmijoLS(f, xk, pk);    
        catch    % 可能会出现zero division的情形
            break;
        end

        % 计算容许步长
        sigma_list = [];
        for i = inactive_ind
            if pk(i)>0
                sigma_list = [sigma_list,(ub(i)-xk(i))/pk(i)];
            elseif pk(i)<0
                sigma_list = [sigma_list,(lb(i)-xk(i))/pk(i)];
            else
                continue;
            end
        end
        sstar = min([sstar,sigma_list]);

        % 更新xkp1
        xkp1 = xk + sstar.*pk;             
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