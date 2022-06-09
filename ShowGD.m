function ShowGD(method_name,f,x_list,y_list)
    %% plot 2d function optimization on 3d space
    % input :
    %     method_name  : the name of iterative optimization method
    %     x_list       : created by iterative optimization method
    %     y_list       : created by iterative optimization method
    %
    % output:
    %     figure in one
    %
    % usage example :
    %     method_name = '某某法';x_list = [1:-0.2:0];y_list = [-1:0.2:0];ShowGD(method_name,x_list,y_list);
    %    
    % todo: 标注出起始点和终止点 completed

    %% 针对输入的(x,y,f(x,y))生成的梯度下降图
    [m,n] = size(x_list);      % default m = 1, n = length(x_list)
    z_list = f(x_list,y_list);
    D2_locz = floor(min(z_list))-1;     % 设定一个最低点，存放二维梯度平面

    % for quiver
    u = x_list(2:end)-x_list(1:end-1);
    v = y_list(2:end)-y_list(1:end-1);
    w = z_list(2:end)-z_list(1:end-1);
    
    % for scatter
    ori_tar_x = [x_list(1),x_list(end)];
    ori_tar_y = [y_list(1),y_list(end)];
    ori_tar_z = [z_list(1),z_list(end)];

    % 2d 梯度下降图
    scatter3(x_list,y_list,D2_locz.*ones(size(x_list)),20,'HandleVisibility','off');hold on;
    quiver3(x_list(1:end-1),y_list(1:end-1),D2_locz.*ones(m,n-1),u,v,zeros(m,n-1),...
        'LineWidth',2,'MaxHeadSize',1,'AutoScale','off','DisplayName',strcat(method_name,'//2d'));
    hold on;
    scatter3(ori_tar_x,ori_tar_y,D2_locz.*ones(size(ori_tar_x)),50,'p','filled','HandleVisibility','off');hold on;
    
    % 3d 梯度下降图
    scatter3(x_list,y_list,z_list,20,'HandleVisibility','off');hold on;
    quiver3(x_list(1:end-1),y_list(1:end-1),z_list(1:end-1),u,v,w,...
        'LineWidth',2,'MaxHeadSize',1,'AutoScale','off','DisplayName',strcat(method_name,'//3d'));
    hold on;
    scatter3(ori_tar_x,ori_tar_y,ori_tar_z,50,'p','filled','HandleVisibility','off');hold on;
    
    xlabel('x_1');ylabel('x_2');zlabel('f(x_1,x_2)');
    
end