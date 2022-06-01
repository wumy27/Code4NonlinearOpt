function X_list = ch5_HooJee(f, X0, lb, ub, epsilon)
    % n元优化函数f(x1,...)，优化区间为[lb1,ub1]x...x[lbn,ubn]
    % fun = matlabFunction(f);  % 优化算法初始输入的符号函数
    syms s                  % 优化变量s
    [n,~] = size(lb);       % n variables  % n=2
    X_list = [];
    X_list = [X_list;X0];   % initialize
    while true
        Xplus = X_list(:,end);    % 最新一次迭代变量
        for j = 1:n
            % 针对xj求最小值
            ej = zeros([n,1]); ej(j) = 1;   % 单位向量
            fxj = matlabFunction(subs(f,symvar(f),[Xplus+s.*ej].'));    % 待最小化函数
            
            % s_list = ch2_BisectionMethod(fxj, lb(j)-Xplus(j), ub(j)-Xplus(j), epsilon);   % 可以很快给出解
            % s_list = ch2_SecantMethod(fxj, lb(j)-Xplus(j),ub(j)-Xplus(j), epsilon);      % 有点反常，但最终会到最优解
            s_list = ch2_ModifiedNewtonMethod(fxj, lb(j)-Xplus(j),ub(j)-Xplus(j), epsilon); % 可以很快到最优解
       
            Xplus(j) = Xplus(j) + s_list(end);       
        end
        % add a line search compare to Univariate
        f2xk = matlabFunction(subs(f,symvar(f),[Xplus+s.*(Xplus-X_list(:,end))].'));    % 待最小化函数 
        s_list = ch2_BisectionMethod(f2xk, 0, 1, epsilon);                         % s\in[0,1] 
        X_list = [X_list,Xplus+s_list(end).*(Xplus-X_list(:,end))];                % 选择最优的
        if norm(X_list(:,end)-X_list(:,end-1))<epsilon
            break;
        end
    end
end