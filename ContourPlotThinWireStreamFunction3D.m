% 3D Plot of Stream functions generated by a thin-wire approach

function PTW3D = ContourPlotThinWireStreamFunction3D(CoilDefinition, ElementCurrents, n_contour)


for nP = 1:CoilDefinition(1).Partitions
    
    PlotCoord = (CoilDefinition(nP).thin_wire_nodes_start + CoilDefinition(nP).thin_wire_nodes_stop)/2;
    
    sx = reshape(PlotCoord(:,1),CoilDefinition(nP).num_elements);
    sy = reshape(PlotCoord(:,2),CoilDefinition(nP).num_elements);
    sz = reshape(PlotCoord(:,3),CoilDefinition(nP).num_elements);
    
    if CoilDefinition(1).StreamDirection == 1
        sx = (sx(1:end-1,:) + sx(2:end,:))/2;
        sy = (sy(1:end-1,:) + sy(2:end,:))/2;
        sz = (sz(1:end-1,:) + sz(2:end,:))/2;
    elseif CoilDefinition(1).StreamDirection == 2
        sx = (sx(:,1:end-1) + sx(:,2:end))/2;
        sy = (sy(:,1:end-1) + sy(:,2:end))/2;
        sz = (sz(:,1:end-1) + sz(:,2:end))/2;
    end
    
    
    if exist('n_contour') == 1
        n_cont = n_contour;
    else
        n_cont = 15;
    end
   
% Add two radial columns to not miss iso-contours at radial
% interconnections
    ElmtsPlot = reshape(ElementCurrents(nP).Stream,(CoilDefinition(nP).num_elements -[0 1]));
    ElmtsPlot = [ElmtsPlot(end,:); ElmtsPlot; ElmtsPlot(1,:)];
    
    
    figure; set(gcf,'Name','Stream Function','Position',[   1   1   500   500]);
    
    hold all
    imab(ElmtsPlot);
    cont_max_main = max(max(ElmtsPlot));
    [C,H] = contour(ElmtsPlot(:,:)',[-cont_max_main:(2*cont_max_main/n_cont):cont_max_main],'k','LineWidth', 2);
    % title('Stream Function');
    hold off
    
    ylabel('z-Axis [m]')
    xlabel('circumferential [rad]')
    
    
    %% 3D Plot of the contours
    
    S = contourdata(C);
    ccount = size(S);
    
    figure; set(gcf,'Name','3D coil','Position',[   1   1   1000   1000]);
    hold all
    for i = 1:ccount(2)
        sx=CoilDefinition(nP).Radius*cos(S(i).xdata/(CoilDefinition(nP).num_elements(1))*2*pi);
        sy=CoilDefinition(nP).Radius.*sin(S(i).xdata/(CoilDefinition(nP).num_elements(1))*2*pi);
        sz=S(i).ydata./length(ElmtsPlot(1,:)')*CoilDefinition(nP).Length - CoilDefinition(nP).Length/2;
        
        plot3(sx,sy,sz,'b','LineWidth', 1)
    end
    
    hold off
    view([-7 25]);
    
    axis tight equal
    font_size = 12;
    set(gca,'fontsize',font_size)
    
    
end

end

