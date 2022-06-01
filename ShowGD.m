function ShowGD(method_name,x_list,y_list)
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
    %     method_name = '改进牛顿法';
    %     x_list = [1:-0.2:0];
    %     y_list = [-1:0.2:0];
    %     ShowGD(method_name,x_list,y_list);
    % 
    % todo: 标注出起始点和终止点 completed
    
    %% f(X,Y)对应图形
    % 定义曲面图的范围
    f =@(X,Y) -2*cos(X) - cos(Y) - cos(X-Y);
    [X,Y] = meshgrid([-10:0.5:10]);    % origin: -10:10
    Z = f(X,Y);

    % 绘制曲面图和等高线图
    figure;
    sc = surfc(X,Y,Z,Z,'FaceAlpha',0.2,'EdgeColor','none');hold on;colorbar;
    sc(2).ZLocation = -5;   % sc(1)是surf, sc(2)是contour

    % 绘制负梯度系统
    [u,v] = gradient(Z);
    u = -u;v = -v;   % 负梯度
    l = streamslice(X,Y,u,v);   % 将梯度与坐标点对应
    % https://matplotlib.org/stable/gallery/images_contours_and_fields/plot_streamplot.html
    set(l,'Color','k');         % 设置梯度线为黑色
    for i=1:length(l) 
        % 绘制梯度图在曲面图上
        l(i).ZData = -2.*cos(l(i).XData) - cos(l(i).YData) - cos(l(i).XData-l(i).YData);
    end

    %% 针对输入的(x,y,f(x,y))生成的梯度下降图
    [m,n] = size(x_list);      % default m = 1, n = length(x_list)
    z_list = f(x_list,y_list);

    % for quiver
    u = x_list(2:end)-x_list(1:end-1);
    v = y_list(2:end)-y_list(1:end-1);
    w = z_list(2:end)-z_list(1:end-1);
    scale = 1;
    
    % for scatter
    ori_tar_x = [x_list(1),x_list(end)];
    ori_tar_y = [y_list(1),y_list(end)];
    ori_tar_z = [z_list(1),z_list(end)];

    % 2d 梯度下降图
    % plot3(x_list,y_list,-5.*ones(size(x_list)),'->','LineWidth',2);hold on;
    scatter3(x_list,y_list,-5.*ones(size(x_list)),20);hold on;
    quiver3(x_list(1:end-1),y_list(1:end-1),-5.*ones(m,n-1),u,v,zeros(m,n-1),'LineWidth',2,'MaxHeadSize',0.5,'AutoScale','off');hold on; % ,'AutoScaleFactor',scale,'AutoScale','on'
    scatter3(ori_tar_x,ori_tar_y,-5.*ones(size(ori_tar_x)),50,'p','filled');hold on;
    
    % 3d 梯度下降图
    % plot3(x_list,y_list,z_list,'->','LineWidth',2);hold on;
    scatter3(x_list,y_list,z_list,20);hold on;
    quiver3(x_list(1:end-1),y_list(1:end-1),z_list(1:end-1),u,v,w,'LineWidth',2,'MaxHeadSize',0.5,'AutoScale','off');hold on;   % ,'AutoScaleFactor',scale,'AutoScale','on'
    scatter3(ori_tar_x,ori_tar_y,ori_tar_z,50,'p','filled');hold on;
    
    xlabel('x_1');ylabel('x_2');zlabel('f(x_1,x_2)');title(method_name);  
    
end