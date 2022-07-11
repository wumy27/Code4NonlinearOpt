function [x_star,lambda_star] = ch17_RG4EQPprogv2(G, h, Aeq, beq)
    % Another implement version of Reduced gradient for Equality Constraint Quadratic Programming
    % 目标函数f的设计不需考虑常数项，因为其不影响求得的最优解，仅影响最优函数值
    % 输入的G必须为正定矩阵，才能保证输出的解是最小值点
    % l: 为所需要被n-l个变元代替的变元数量[hyper parameter]
    % 此版本不需要人工设置l的大小

    % Aeq*Q=[L;0]--Aeq=[L;0]*Q'
    % H的LQ分解可以借助QR分解完成:
    % [Q,R]= qr(H');H=(Q*R)'=R'*Q'=L*Q'; 
    [Q,R]= qr(Aeq.');L=R.'; % Q=Q;
    [l,~]=size(L);                   % l: size of L
    Y = Q(:,1:l); Z = Q(:,l+1:end);  % range/null space
    
    % Get (x_star,lambda_star)
    % 参数输入形式为Aeq*x=beq
    % 但是下面求解公式面对的是A*x+b=0,
    % 故需令beq=-b
    ybar = (Aeq*Y)\beq;   % 这里与原公式(17.3.5)差个负号
    zbar = -(Z.'*G*Z)\(Z.'*h+Z.'*G*Y*ybar);
    x_star = Y*ybar + Z*zbar;
    lambda_star = (Y.'*Aeq.')\(Y.'*h+Y.'*G*x_star);    
end