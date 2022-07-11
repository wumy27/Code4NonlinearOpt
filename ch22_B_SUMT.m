function X_list = ch22_B_SUMT(f, c, x0, epsilon)
    % f: 目标函数
    % c: 不等式约束函数[vector] c_i(x)>=0, i=1,...,l
    % x0: 初始点[需要为不严格满足约束的可行解]
    % epsilon: 误差限

    X_list = x0
    r = 1e3*eval(subs(f,symvar(f),x0.'));   % penalty function parameters(r up, P up)
    beta = 0.5;                             % decay parameters
    while true
        xk = X_list(:,end);
        % penalty function
        B = f + r*sum(1./c);                              % form 1
%         B = f - r*sum(log(c));                            % form 2[c<0无定义]

        % solve unconstraint optimization problem
        try
            x_list = ch9_NewtonMethod(B, xk, epsilon);    % iter list
            xkp1 = x_list(:,end)                         % result
        catch    % 可能会出现zero division的情形
            break;
        end
        X_list = [X_list,xkp1]; 
        r = r*beta;                                   % decay r
        % lagrange multipliers
        lambda = r./(eval(subs(c,symvar(f),xkp1.')).^2);   % form 1
%         lambda = r./eval(subs(c,symvar(f),xkp1.'));        % form 2
        
        if max(abs(lambda.*eval(subs(c,symvar(c),xkp1.'))))<epsilon
            break;
        end
    end
end