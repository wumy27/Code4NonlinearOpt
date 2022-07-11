function [X_list,xkp1] = ch21_P_SUMT4IC(f, ceq, cin, x0, epsilon)
    % f: 目标函数
    % ceq/cin: (不)等式约束函数[vector] 
    % min f(x) s.t. ceq_i(x)=0, i=1,...,l cin_j(x)>=0, j=l+1,...,m
    % x0: 初始点
    % epsilon: 误差限
    % xkp1 == x_star

    X_list = x0;
    r = 1e-3./eval(subs(f,symvar(f),x0.'));   % penalty function parameters(r decrease, P up)
    beta = 0.1;                               % decay parameters
    cin_ = min([zeros(size(cin)),cin],[],2);  % for penalty function & convegence condition

    while true
        xk = X_list(:,end);

        % penalty function
%         P = f + 1/r*(sum(ceq.^2) + sum(cin_.^2));  % // 对用牛顿法求解出现bug
        P = f + 1/r*(sum(ceq.^2) + sum(-cin.^3));    % cin>=0,误差很小;cin<0误差很大

        % solve unconstraint optimization problem
        try
            x_list = ch9_NewtonMethod(P, xk, epsilon); xkp1 = x_list(:,end);  
        catch    % 可能会出现zero division的情形
            break;
        end
        X_list = [X_list,xkp1];
        
        r = r*beta;             % decay r
            
        if norm(eval(subs(ceq,symvar(f),xkp1.'))) + norm(eval(subs(cin_,symvar(f),xkp1.')))<epsilon
            break;
        end
    end
end