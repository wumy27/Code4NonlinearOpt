function [X_list,P_list,m,t] = ch5_MultiDIRECT(f, lb0, ub0, fg, epsilon)
    % X_list存储当前最小值的位置
    % P_list存储所有评估过点的位置
    % m=length(P_list)=3t+1
    % t为迭代次数

    m = 1;           % m用来count当前计算了几次函数值
    t = 0;           % iteration count
    T = 10;         % 最大迭代次数
    % 将初始矩形域加入候选集中
    S = [lb0,ub0];

    % 初始化最小值
    xmin  = (ub0+lb0)./2;
    fmin = eval(subs(f,symvar(f),xmin.'));   % f(xmin);
    X_list = xmin;
    P_list = xmin;

    while (~isempty(S)&&(t<T))
        lb = S(:,1); ub = S(:,2);
        S(:,1:2) = [];

        % 找到矩形域的重心
        c1 = (ub+lb)./2;     % 列向量

        % 获取矩形域的最大维度
        [len,Ind]=max(ub-lb);

        % 针对该维度对矩形域进行三分 
        eInd = zeros(size(lb)); eInd(Ind) = 1;
        c2 = c1+len/3.*eInd; c3 = c1-len/3.*eInd;

        % 分别计算这三个细分后的小矩形域重心处的函数值
        fc1 = eval(subs(f,symvar(f),c1.'));
        fc2 = eval(subs(f,symvar(f),c2.'));
        fc3 = eval(subs(f,symvar(f),c3.'));
        m = m+3;
        t = t+1;

        % 找出函数值最小的小矩形域加入候选集
        xmin_cand = [c3,c1,c2];
         [fmin_comp,fInd]= min([fc3,fc1,fc2]);
        if fmin_comp < fmin
            fmin = fmin_comp;
            xmin = xmin_cand(:,fInd);
        end
        X_list = [X_list,xmin];
        P_list = [P_list,c3,c2];

        lb(Ind) = (-1/2+(fInd-1)/3)*len + c1(Ind);
        ub(Ind) = (-1/2+fInd/3)*len + c1(Ind);
        S = [S,lb,ub];
        
        % fg = f_global
        if abs(fmin-fg)<epsilon
            break;
        end
        
    end
end