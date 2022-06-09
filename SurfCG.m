function SurfCG(f)
    % plot surface contour gradient of f 

    % 绘制曲面图和等高线图
    [X,Y] = meshgrid([-10:0.5:10]);    % origin: -10:10
    Z = f(X,Y);

    sc = surfc(X,Y,Z,Z,'FaceAlpha',0.2,'EdgeColor','none');hold on;
    sc(1).HandleVisibility='off';sc(2).HandleVisibility='off';
    colorbar;
    D2_locz = min(Z(:))-1;            % 二维图像所处的z坐标值
    % sc(2).ZLocation = D2_locz;   % sc(1)是surf, sc(2)是contour

    % 绘制负梯度系统
    [u,v] = gradient(Z);
    u = -u;v = -v;   % 负梯度
    l = streamslice(X,Y,u,v);hold on;   % 将梯度与坐标点对应
    % https://matplotlib.org/stable/gallery/images_contours_and_fields/plot_streamplot.html
    set(l,'Color','k');         % 设置梯度线为黑色
    set(l,'HandleVisibility','off');
    for i=1:length(l) 
        % 绘制梯度图在曲面图上
        l(i).ZData = f(l(i).XData,l(i).YData);
        l(i).HandleVisibility='off';
    end
end