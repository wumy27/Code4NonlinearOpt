function [xkp1,k] = ch11_CG4LinearEqts( A, b, x0, epsilon)
    % for solving Ax + b = 0
    % Example: A = [1,2;3,4];b=[-5;-11];[x,k]=ch11_CG4LinearEqts( A, b, [0;0], 1e-6)
    % slow to converge

    % initialization
    xk = x0;
    gk = A*xk+b;
    pk = -gk;
    k = 0;             % k记录迭代次数
    while true  
        k = k + 1; 
        if k>=2   % k = 1不存在x/g/pkp1
            xk = xkp1;
            gk = gkp1;
            pk = pkp1;
        end
        
        % find s  
        sstar = (-pk.'*gk)/(pk.'*A*pk);
        
        % set xkp1 & gkp1
        xkp1 = xk + sstar.*pk;
        gkp1 = A*xkp1+b;
        % determine beta & pkp1
        beta = gkp1.'*gkp1/(gk.'*gk);        % Fletcher-Reeves公式 迭代次数 = 95
%         beta = gkp1.'*(gkp1-gk)/(gk.'*gk);   % Polak-Ribere公式 迭代次数 = 479
        pkp1 = -gkp1 + beta.*pk;

        if norm(gkp1)<epsilon   % 该收敛条件很难收敛
            break;
        end
    end
end