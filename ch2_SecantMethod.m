function xkp2_list = ch2_SecantMethod(f0, a, b, epsilon)
    % 割线法针对的目标函数为一元函数
    syms X
    f = matlabFunction(diff(f0(X)));      % f'
    [xk, xkp1] = deal(a, b);
    xkp2_list = []; % 存储x_{k+2}
    xkp2 = xk - f(xk)/(f(xkp1)-f(xk))*(xkp1-xk);
    xkp2_list = [xkp2_list,xkp2];
    while abs(f(xkp2))>epsilon
        [xk, xkp1] = deal(xkp1, xkp2);
        xkp2 = xk - f(xk)/(f(xkp1)-f(xk))*(xkp1-xk);
        xkp2_list = [xkp2_list,xkp2];
    end
end