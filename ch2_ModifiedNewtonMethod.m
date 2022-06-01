function xkp1_list = ch2_ModifiedNewtonMethod(f, a, b, epsilon)
    syms X
    fx = matlabFunction(diff(f(X)));      % f'
    fxx = matlabFunction(diff(fx(X)));    % f''
    % 一元函数优化: f:=f(x) fx:=f'(x) fxx:=f''(x)
    xk = b;
    xkp1_list = [];
    xkp1_list = [xkp1_list,xk];
    while true
        xk = xkp1_list(end);
        if fxx(xk)>0
            dx = -fx(xk)/fxx(xk);
        else
            dx = -fx(xk);
        end
        if dx<0
            alpha = min(1,(a-xk)/dx);
        end
        if dx>0
            alpha = min(1,(b-xk)/dx);
        end
        while f(xk+alpha*dx)>=f(xk)
            alpha = alpha/2;
        end
        xkp1 = xk + alpha*dx;
        xkp1_list = [xkp1_list,xkp1];
        if (abs(fx(xkp1))<epsilon)
            break;
        end
    end
 end