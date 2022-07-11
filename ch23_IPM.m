function X_list = ch23_IPM(f,c,x0,epsilon)
    % min f(x)  s.t. c(x)>=0
    % = min f(x)-r*sum(log(w)) s.t. c(x)-w=0,w>=0,r>0
    
    % symvar for decision variables
    % symvar(f)                   % decision variables
    w = sym('w',[1,length(c)]).';   % slack variables
    lambda = sym('lambda',[1,length(c)]).';   % lagrangian multipliers
    r = sym('r').';   % penalty function parameters(r decrease, P up)
    v = sym('v',[1,length(c)]).';   % parameter vector for M

    % jacobian and hessian for f & c
    df = jacobian(f,symvar(f)).';   
    dc = jacobian(c,symvar(f));   
%     d2f = hessian(f,symvar(f));
%     d2c = hessian(c,symvar(f));    

    % 参数表达式
%     G_tilde = d2f - lambda.'*d2c;    
    g = df;        % df
    A = dc;        % row = jacobian(ci)
    W = diag(w);
    Lmda = diag(lambda);

    % initialization
    X_list = x0
    wk = ones(size(c));
    lambda_k = ones(size(c));
    Bk = eye(length(symvar(f)));            % 估计矩阵G_tilde
    rk = 1e3*eval(subs(f,symvar(f),x0.')); 
    beta = 0.5;                             % decay parameters
    vk = ones(size(c));                     % approx lambda_k    
    e = ones(size(c));                      % all 1 vector
    % error function
    tau = norm(c-w)^2 + norm(g-A.'*lambda)^2 + norm(W*Lmda*e-r*e)^2;
    tau_minus = eval(subs(tau,[symvar(f), w.',lambda.',r],[x0.',wk.',lambda_k.',rk]));

    % 阈值设置
    lambda_min = 0.5;   % 正数

    % lagrangian function para: x,w,lambda
    L = f - r*sum(log(w)) - lambda.'*(c-w);
    dL = jacobian(L, symvar(f)).';   % L'[vector] == d(f - lambda.'*c)/dx    

    % enhanced lagrangian function of barrier nonlinear programming
    M = f - r*sum(log(w)) - (c-w).'*v + norm(c-w)^2/r;
    

    while true
        xk = X_list(:,end);
        gkp0 = eval(subs(dL, [symvar(f), lambda.'], [xk.',lambda_k.']));   % for Bk
        % obtain delta_x,delta_lambda,delta_w
        Ak = eval(subs(A,symvar(f),xk.'));
        gk = eval(subs(g,symvar(f),xk.'));
        Lmdak = diag(lambda_k);
        Wk = diag(wk);
        ck = eval(subs(c,symvar(f),xk.'));
        
%         % 直接使用精确G_tilde [not work!]
%         Bk = eval(subs(hessian(L,symvar(f)),[symvar(f),lambda.'],[xk.',lambda_k.']));

        sol = [Bk,-Ak.';-Ak,-Wk/Lmdak]\[-gk+Ak.'*lambda_k;ck-rk.*Lmdak\e];
        delta_x = sol(1:length(symvar(f)));
        delta_lambda = sol(length(symvar(f))+1:end);
        delta_w = rk.*Lmdak\e - wk - Wk/Lmdak*delta_lambda;
        
        % get shat
        s_list = [];
        for i = 1:length(c)
            if delta_w(i)<0
                s_list = [s_list,-0.9*wk(i)/delta_w(i)];
            end
        end
        if length(s_list) == 0
            shat = 1;    % 设置为1
        else
            shat = min(s_list);
        end
        
        % get lambda_plus
        lambda_plus = max([lambda_min.*ones(size(lambda_k)),lambda_k+delta_lambda],[],2);

        % line search to update xkp1, wkp1
        w_x = [wk;xk];   % w1,w2,x1,x2
        pk = [delta_w;delta_x];
        
        sstar = ch8_ModifiedArmijoLS(subs(M,[v.',r],[vk.',rk]), w_x, pk);
        sstar = min([shat,sstar]);
        % 满足约束的步长
%         sstar = ch23_WLS4C(subs(M,[v.',r],[vk.',rk]),w_x,pk,shat)      % [not work]
%         sstar = ch23_WLS4Cv2(subs(M,[v.',r],[vk.',rk]), w_x, pk, shat)   % [收敛到另外的解]
        xkp1 = xk+sstar.*delta_x
        wkp1 = wk+sstar.*delta_w;
        X_list = [X_list,xkp1];
        wk = wkp1;

        % lambda_k = lambda_k; rk = rk; vk = vk;
        if eval(subs(tau,[symvar(f),w.',lambda.',r],[xkp1.',wkp1.',lambda_plus.',rk]))<tau_minus
            tau_minus = eval(subs(tau,[symvar(f),w.',lambda.',r],[xkp1.',wkp1.',lambda_plus.',rk]));
            lambda_k = lambda_plus; rk = beta*rk; vk = lambda_k;
        end

        % update Bk[ to do ]
        gkp1 = eval(subs(dL, [symvar(f), lambda.'], [xkp1.',lambda_k.']));
        gammak = gkp1 - gkp0;   % 计算梯度差
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

        % 收敛条件
        if (norm(eval(subs(tau,[symvar(f),w.',lambda.',r],[xkp1.',wkp1.',lambda_k.',0])))<epsilon)|(norm(xkp1-xk)<epsilon)
            break;
        end
    end
end


