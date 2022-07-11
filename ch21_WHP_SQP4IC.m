function [X_list,xkp1] = ch21_WHP_SQP4IC(f, ceq, cin, x0, epsilon)
    % f: 目标函数
    % ceq/cin: (不)等式约束函数[vector] 
    % min f(x) s.t. ceq_i(x)=0, i=1,...,l cin_j(x)>=0, j=l+1,...,m
    % x0: 初始点
    % epsilon: 误差限
    % xkp1 == x_star

    % initialization
    X_list = x0;
    r = 1e-3./eval(subs(f,symvar(f),x0.')); % penalty function parameters(r decrease, P up)
    beta = 0.1;                             % decay parameters
    % parameter vector to approximate lambda
    veq = ones(size(ceq)); vin = ones(size(cin));    % must be non-negetive  
    cin_ = min([zeros(size(cin)),cin],[],2);         % for convergence condition

    % function derivate
    Bk = eye(length(symvar(f)));         % estimate G(x)
    df = jacobian(f,symvar(f)).';        % g(x)
%     d2f = hessian(f,symvar(f));          % G(x)
    dceq = jacobian(ceq,symvar(f));    
    dcin = jacobian(cin,symvar(f));
%     d2ceq = hessian(ceq,symvar(f));  
%     d2cin = hessian(ceq,symvar(f));
    
    while true
        xk = X_list(:,end);
        P = f + 1/r*( sum((ceq-r/2.*veq).^2) + sum(-(cin-r/2.*vin).^3) );    % al_penalty function
        dP = jacobian(P,symvar(f)).';                                        % nabla P

        gk = eval(subs(dP, symvar(f), xk.')); 
 
        % 解iqp子问题
        % method-1: use EIQP4pdH [bug]
        % qp_f = 1/2.* (p.' * Bk * p) + eval(subs(df,symvar(f),xk.')).' * p;
        % qp_c = [ceq + dceq.'*p;-ceq - dceq.'*p;cin + dcin.'*p];
        G = Bk; h = eval(subs(df,symvar(f),xk.')); 
        Age = eval(subs(dcin,symvar(f),xk.')); Aeq = eval(subs(dceq,symvar(f),xk.'));
        bge = eval(subs(cin,symvar(f),xk.'));  beq = eval(subs(ceq,symvar(f),xk.'));
        [~,pk,~] = ch21_EIQP4pdH(G, h, Aeq, beq, Age, bge, zeros(size(xk)));

%         % method-2: use P_SUMT4IC
%         qp_f = 1/2.* (symvar(f) * Bk * symvar(f).') + eval(subs(df,symvar(f),xk.')).' * symvar(f).';
%         [~,pk] = ch21_P_SUMT4IC(qp_f, ceq, cin, xk, epsilon);

%         % method-3: matlab quadprog
%         G = Bk; h = eval(subs(df,symvar(f),xk.')); 
%         Age = eval(subs(dcin,symvar(f),xk.')); bge = eval(subs(cin,symvar(f),xk.'));
%         Aeq = eval(subs(dceq,symvar(f),xk.')); beq = eval(subs(ceq,symvar(f),xk.'));
%         pk = quadprog(G,h,-Age,-bge,Aeq,beq);

        % 线搜索求解更新步长s     
        sstar = ch8_ModifiedArmijoLS(P, xk, pk);
        xkp1 = xk + sstar.*pk
        X_list = [X_list, xkp1];

        % 更新Bk
        gkp1 = eval(subs(dP, symvar(f), xkp1.')); 
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
 
        % update veq
        veq = veq - 2.*eval(subs(ceq,symvar(ceq),xkp1.'))./r;   
        % update vin
        ck = eval(subs(cin,symvar(cin),xkp1.'));
        change_ind = find(ck<r/2.*vin);
        vin_ = zeros(size(vin));
        vin_(change_ind) = vin(change_ind) - 2/r.*ck(change_ind);
        vin = vin_;

        r = r*beta;                                  % decay r

        % 最优性判定：将误差度量函数T改为用约束的收敛来判定
        if norm(eval(subs(ceq,symvar(f),xkp1.'))) + norm(eval(subs(cin_,symvar(f),xkp1.')))<epsilon
            break;
        end
    end
end