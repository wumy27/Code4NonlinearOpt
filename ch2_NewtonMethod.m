function xkp1_list = ch2_NewtonMethod(f, x0, epsilon)
    % 一元函数优化: f:=f(x) x0:=初值
    syms X
    fx = matlabFunction(diff(f(X)));      % f'
    fxx = matlabFunction(diff(fx(X)));    % f''
    xk = x0;
    xkp1 = xk - fx(xk)/fxx(xk);
    % disp([fxx(xk)/fxx(xk)]) % = 14.10141994717172,步长太大
    xkp1_list = [];
    xkp1_list = [xkp1_list,xkp1];
    while abs(fx(xkp1))>epsilon
        xk = xkp1;
        xkp1 = xk - fx(xk)/fxx(xk);
        xkp1_list = [xkp1_list,xkp1];
    end
end