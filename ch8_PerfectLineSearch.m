function sstar = ch8_PerfectLineSearch(f, xk, p)
    % 对某点xk沿着p进行完美线搜索
    syms s
    epsilon = 1e-6;
    phi = matlabFunction(subs(f,symvar(f),[xk+s.*p].'));   
    % 需要预先给定上下界
    lb0 = xk - ones(size(xk)); ub0 = xk + ones(size(xk));               
    lbANDub = [(lb0 - xk)./p,(ub0 - xk)./p];  % 计算步长s的容许上下界
    lb = max(min(lbANDub,[],2)); ub = min(max(lbANDub,[],2));
    s_list = ch2_BisectionMethod(phi, lb, ub, epsilon);
    sstar = s_list(end);    % 最优步长
end