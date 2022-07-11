function [X_list,x,lambda]=ch21_IQP4pdH(G,h,Age,bge,x0)
    % Inequality QP algorithm for a positive-definite Hessian
    % min 1/2*(x.'*G*x) + h.'*x + c(can be ignored)   s.t. Age*x + bge >= 0
    % 可以拓展到 Aeq*x = beq
    % -Aeq*x-beq>=0 & Aeq*x+beq>=0  
    % 可这个idea不可以，因为Age会成为奇异矩阵
    
    % initialization 
    epsilon = 1e-6;           % 直接设置误差限（数值算法一定有误差）
    X_list = x0;
    lambda = zeros(size(bge));
    b = zeros(size(bge));

    while true
        x = X_list(:,end);

        % t设置为随active_ind变化而变化
        % 求有效约束的index 1~t              
        active_ind = [find(Age*x+bge<-epsilon).',find(abs(Age*x+bge)<epsilon & lambda>=0).'];
        % 删除lambda<0的约束
        active_ind = setdiff(active_ind,find(lambda<0).');
        t = length(active_ind);          % 有效约束的个数
        % inactive_ind即全集删去active_ind
        inactive_ind = [1:length(bge)]; inactive_ind(active_ind)=[];     
        
        % active_ind 对应 1~t
        g = G*x+h;
        b(active_ind) = Age(active_ind,:)*x+bge(active_ind);
        
        if t>0
            Aeq = Age(active_ind,:); beq = -b(active_ind);
            [p_star,mu_star] = ch17_EQPprog(G, g, Aeq, beq);
            lambda(active_ind)=mu_star; lambda(inactive_ind) = 0;
        else % 针对Aeq与beq为空集的情况，变为二次型函数
            % 一阶最优性条件: G*p+g=0  --> p = -G\g
            % 此时没有active_ind
            p_star = -G\g;
            lambda(inactive_ind) = 0; % equal to lambda(inactive_ind)=zeros(size(lambda(inactive_ind)));
        end
        
        % 找到一个可行的步长
        s = 1; 
        for i = 1:length(bge)-t
            ai_hat = Age(inactive_ind(i),:);
            bi_hat = bge(inactive_ind(i));
            if ai_hat*p_star<0
                s = min([s,-(ai_hat*x+bi_hat)/(ai_hat*p_star)]);
            end
        end
        x = x+s.*p_star;
        X_list = [X_list,x];

        % 收敛条件: optimity condition
        if ((Age*x+bge>=0)&(abs(G*x+h-Age.'*lambda)<epsilon)&(all(lambda>=0)))|(abs(x-X_list(:,end-1))<epsilon)
            break;
        end
    end
end


