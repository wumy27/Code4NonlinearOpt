function X_list = ch19_AL_SQP(f, c, x0, epsilon)
    
    lambda = sym('lambda',[1,length(c)]);
    syms r

    % T parameter:x,lambda
    kappa = 1e5;      % 惩罚参数
    df = jacobian(f,symvar(f)).';
    dc = jacobian(c,symvar(c)).';
    T = norm( df - lambda.' * dc) + kappa * norm( c ); 
    
    % initialization
    X_list = x0;
    lambda_k = ones(size(c));               % lagrangian multipliers
    rk = 1e-3./eval(subs(f,symvar(f),x0.')); % penalty function parameters(r decrease, P up)
    Hk = eye(length(symvar(f)));
    beta = 0.1;                             % decay parameters
    uk = zeros(size(c));                    % equal to lambda_k

    Tminus = eval(subs(T,[symvar(f),lambda.'],[x0.',lambda_k.']));
    
    % al_penalty function(parameter:x,lambda,r)  
    M = f + 1/r.*sum((c-r/2.*lambda).^2);  

    % A and G
    A = jacobian(c,symvar(c));    % A[not change] 
    G = df;                       % G == df
    
    while true
        xk = X_list(:,end);
        ck = eval(subs(c,symvar(f),xk.'));
        gk = eval(subs(G,symvar(f),xk.'));
        Ak = eval(subs(A,symvar(f),xk.'));

        % obtain uk, pk
        uk = (rk/2*eye(length(uk))+Ak*Hk*Ak.')\(Ak*Hk*gk-ck+rk/2*lambda_k);
        pk = Hk*(Ak.'*uk-gk);
        
        % line search 
        sstar = ch8_ModifiedArmijoLS(subs(M,[lambda.',r],[lambda_k.',rk]),xk,pk);
        xkp1 = xk + sstar.*pk;
        X_list = [X_list, xkp1];

        if eval(subs(T,[symvar(f),lambda.'],[xk.',uk.']))<Tminus
            rk = beta*rk; lambda_k = uk; 
            Tminus = eval(subs(T,[symvar(f),lambda.'],[xkp1.',lambda_k.']));
        % else 
        %     rk = rk; lambda_k = lambda_k;
        end

        % update Hk
        % 计算梯度差与解间隔
        gkp1 = eval(subs(G,symvar(df),xkp1.'));
        gammak = gkp1 - gk;
        deltak = xkp1 - xk;
        
        % 通过求解Hkp1*gammak = deltak得到一个正定矩阵Hkp1更新Hk
        if deltak.'*gammak > 0   % 当其为正数时，才能保证更新公式给出的Hk正定
        % BFGS formula 
            Hk = Hk - (Hk*gammak*(deltak.')+deltak*gammak.'*Hk)/((deltak.')*gammak) + ...
                ( 1 + (gammak.'*Hk*gammak)./(deltak.'*gammak) ).*(deltak*deltak.')./(deltak.'*gammak);
        else
            % 两种方法的迭代次数分别为42/103
            % method 1:随机给定一个正定矩阵更新Hk
            Hk = eye(length(symvar(f)));  

%             % method 2:[原版本为Bk]改编自书本P178 (17.3.11)给出用etak代替gammak的公式 % 可能不收敛
%             for theta = 1:-0.1:0
%                 etak = (1-theta).*deltak+theta.*Hk*gammak;
%                 if gammak.'*etak>0    % 只能保证略微大于0
%                     Hk = Hk - (Hk*gammak*(etak.')+etak*gammak.'*Hk)/((gammak.')*etak) + ( 1 + (gammak.'*Hk*gammak)./(gammak.'*etak) ).*(etak*etak.')./(gammak.'*etak);
%                     break;
%                 end
%             end

        end

        if (eval(subs(T,[symvar(f),lambda.'],[xkp1.',lambda_k.']))<epsilon)|(abs(xkp1-xk)<epsilon*1e2) % 用判定公式可以收敛，加上迭代差更快收敛[很难收敛]
            break;
        end

    end

end