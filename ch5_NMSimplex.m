function X_list = ch5_NMSimplex(f, lb, ub, epsilon, X0)
    % lb(nx1) ub(nx1) X0(nx{n+1})
    [n,~] = size(lb);                  % n variable -> n+1 Simplex
    if(~exist('X0','var'))
        Xs = lb+(ub-lb).*rand(n,n+1);  % 每一列即为一个simplex点  % if
    else
        Xs = X0;
    end
    X_list = [];
    while true
        % 计算当前Simplexs的值
        FX = [];
        for i = 1:n+1
            FX = [FX;eval(subs(f,symvar(f),Xs(:,i).'))];  
        end
        % 计算当前Simplexs种的最大值和最小值
        [Fxw,xwInd] = max(FX); xw = Xs(:,xwInd);
        [Fxb,xbInd] = min(FX); xb = Xs(:,xbInd);
        
        % 添加中间迭代点
        X_list = [X_list,Xs(:,xbInd)];
        
        % 计算去极值重心和重心
        xtilde = 1/(n+1-1).*(sum(Xs,2)-xw);
        xbar = 1/(n+1).*sum(Xs,2);
        % 计算反射点
        xplus = xw+2.*(xtilde-xw);
        % 更新
        if eval(subs(f,symvar(f),xplus.'))<Fxb
            xplusplus = xplus + (xplus-xw);
            if eval(subs(f,symvar(f),xplusplus.'))<eval(subs(f,symvar(f),xplus.'))
                xplus = xplusplus;
            end
            Xs(:,xwInd) = xplus;
        else
            if eval(subs(f,symvar(f),xplus.'))<Fxw
                Xs(:,xwInd) = xplus;
            else
                xplus = xw + 4/3.*(xtilde - xw);
                if eval(subs(f,symvar(f),xplus.'))<Fxw
                    Xs(:,xwInd) = xplus;
                else
                    for i = 1:n+1
                        if i == xbInd
                            continue;
                        end
                        Xs(:,i) = 1/2.*(Xs(:,i)+xb);
                    end
                end
            end
        end
        
        % 收敛条件
        if norm(Xs-xbar.*ones(size(Xs)))<epsilon
            break;
        end
    end

    % 添加最终迭代点
    FX = [];
    for i = 1:n+1
        FX = [FX,eval(subs(f,symvar(f),Xs(:,i).'))];  
    end
    [Fxb,xbInd] = min(FX); xb = Xs(:,xbInd);
    X_list = [X_list,Xs(:,xbInd)];
end