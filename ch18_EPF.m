function [X_list,x_star,lambda_star] = ch18_EPF(f, c, x0, epsilon)
    % exact penalty function
    % c: vector of ci
    
    % penalty function parameters(r decrease, P up)
    r = 1e-3./eval(subs(f,symvar(f),x0.')); 

    % first form of E  % 非光滑函数，难以求导，算法失败
%     E = f + 1/r * sum(abs(c));
    
    % second form of E % 只能选择E的形式2
    A = jacobian(c,symvar(c));  % matrix of jacobian
    g = gradient(f, symvar(f));
    E = f + c.'/(A*A.')*A*g + 1/r.*c.'*c;

    X_list = ch9_NewtonMethod(E, x0, epsilon);
    x_star = X_list(:,end);
    lambda_star = (A*A.')\A*g;
end