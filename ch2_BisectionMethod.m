function x_list = ch2_BisectionMethod(f, a, b, epsilon)
    % 二分法针对的目标函数为一元函数
    [xa, xb, xm] = deal(a, b, (a+b)/2);
    [Fa, Fb, Fm] = deal(f(xa), f(xb), f(xm));
    x_list = []; % 存储x_{min}
    if abs(xb-xa)<=epsilon
        x_list = [x_list,(xa+xb)/2];
    end
    while abs(xb-xa)>epsilon
        [xl, xr] = deal((xa+xm)/2, (xm+xb)/2);
        [Fl, Fr] = deal(f(xl), f(xr));
        x_candidate = [xa,xb,xm,xl,xr];
        [Fmin,xminInd] = min([Fa,Fb,Fm,Fl,Fr]);
        if (Fmin == Fa) || (Fmin == Fl)
            [xb, xm] = deal(xm, xl);
            [Fb, Fm] = deal(Fm, Fl);
        elseif Fmin == Fm
            [xa, xb] = deal(xl, xr);
            [Fa, Fb] = deal(Fl, Fr);
        else
            [xa, xm] = deal(xm, xr);
            [Fa, Fm] = deal(Fm, Fr);
        end
        x_list = [x_list,x_candidate(xminInd)];
    end
end