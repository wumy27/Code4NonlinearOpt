function [X_list,xkp1] = ch21_AL_SUMT4IC(f, ceq, cin, x0, epsilon)
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

    while true
        xk = X_list(:,end);
        P = f + 1/r*( sum((ceq-r/2.*veq).^2) + sum(-(cin-r/2.*vin).^3) );    % al_penalty function

        % solve unconstraint optimization problem
        try
            x_list = ch9_NewtonMethod(P, xk, epsilon); xkp1 = x_list(:,end);  
        catch    % 可能会出现zero division的情形
            break;
        end
        X_list = [X_list,xkp1];

        % update veq
        veq = veq - 2.*eval(subs(ceq,symvar(ceq),xkp1.'))./r;   
        % update vin
        ck = eval(subs(cin,symvar(cin),xkp1.'));
        change_ind = find(ck<r/2.*vin);
        vin_ = zeros(size(vin));
        vin_(change_ind) = vin(change_ind) - 2/r.*ck(change_ind);
        vin = vin_;

        r = r*beta;                                  % decay r

        if norm(eval(subs(ceq,symvar(f),xkp1.'))) + norm(eval(subs(cin_,symvar(f),xkp1.')))<epsilon
            break;
        end
    end
end