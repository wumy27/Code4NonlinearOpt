function X_list = ch18_AL_SUMT(f, c, x0, epsilon)
    % f: 目标函数
    % c: 等式约束函数[matrix] c_i(x)=0, i=1,...,l
    % v: parameter vector
    % x0: 初始点
    % epsilon: 误差限

    X_list = x0;
    r = 1e-3./eval(subs(f,symvar(f),x0.')); % penalty function parameters(r decrease, P up)
    beta = 0.1;                             % decay parameters
    v = ones(size(c));                      % parameter vector
    while true
        xk = X_list(:,end);
        P = f + 1/r.*sum((c-r/2.*v).^2);    % al_penalty function
        
        % solve unconstraint optimization problem
        x_list = ch9_NewtonMethod(P, xk, epsilon);        % iter list
        xkp1 = x_list(:,end);   % result
        X_list = [X_list,xkp1]; 
        v = v - 2.*eval(subs(c,symvar(c),xkp1.'))./r;      % update v
        r = r*beta;                                  % decay r
        if norm(eval(subs(c,symvar(c),xkp1.')))<epsilon
            break;
        end
    end
end