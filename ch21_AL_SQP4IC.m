function [X_list,xkp1] = ch21_AL_SQP4IC(f, ceq, cin, x0, epsilon)
    % f: 目标函数
    % ceq/cin: (不)等式约束函数[vector] 
    % min f(x) s.t. ceq_i(x)=0, i=1,...,l cin_j(x)>=0, j=l+1,...,m
    % x0: 初始点
    % epsilon: 误差限
    % xkp1 == x_star

    % 求解约束中存在未知变量的二次规划问题，无法用前面的方法求解

    % todo
    X_list = [];
    xkp1 = [];

end