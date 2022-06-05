function sstar = ch8_ModifiedArmijoLS(f, xk, p)
    % 对某点xk沿着p进行线搜索
    syms s
    % 设计函数f'(x)
    df = jacobian(f,symvar(f));      % f的一阶导函数f'[vector]
    g = eval(subs(df,symvar(df),xk.')).';    % =f'(xk)要做一个转置
    
    % 输入的p需要满足(8.1.1)的搜索方向
    if ~exist('p','var')
        p = -g;   % =-f'(xk)
    end
    
    % 超参数
%     C = 1.1; c = 0.1; eta1 = 0.1; eta2 = 0.1;
    C = 1.5; c = 0.5; eta1 = 0.25; eta2 = 0.25;
    % 步长初始化
    sstar = 1; smin = 0;
    % 设定比率函数
    fxk = eval(subs(f,symvar(f),xk.'));
    D = matlabFunction(subs(f,symvar(f),[xk+s.*p].')-fxk)./(s.*p.'*g);
        % 分子 = f(xk+s.*p)-f(xk)
    while true    % 扩大步数，让s远离0
        toCompare = abs(0.5*sstar/(1-eval(subs(D,s,sstar))));
        if abs(1-eval(subs(D,s,sstar)))<eta2   % D(sstar)
            smin = sstar; 
            sstar = min([C*sstar,toCompare]);
        else
            break;
        end
    end
    
    while true    % 不断缩短步长，让s远离sbar
        toCompare = abs(0.5*sstar/(1-eval(subs(D,s,sstar))));
        sstar = max([smin + c*(sstar - smin),toCompare]);
        if eval(subs(D,s,sstar))>=eta1    % D(sstar)
            break;
        end
    end
end