% 3D Plot of Stream functions generated by a thin-wire approach
% 2019-11
% Sebastian Littin
% sebastian.littin@uniklinik-freiburg.de

function PTW3D = PlotThinWireStreamFunction3D(CoilDefinition, ElementCurrents)


for nP = 1:CoilDefinition(1).Partitions

PlotCoord = (CoilDefinition(nP).thin_wire_nodes_start);% + CoilDefinition(nP).thin_wire_nodes_stop)/2;

% sx = PlotCoord(:,1);
% sy = PlotCoord(:,2);
% sz = PlotCoord(:,3);

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

ElmtsPlot = reshape(ElementCurrents(nP).Stream,(CoilDefinition(nP).num_elements -[0 1]));

% tiny Hack to plot all radial elements
% sx_p = [sx; sx(1,:)];
% sy_p = [sy; sy(1,:)];
% sz_p = [sz; sz(1,:)];
% ElmtsPlot_p = [ElmtsPlot; ElmtsPlot(1,:)];

figure; set(gcf,'Name','3D coil','Position',[   1   1   1000   1000]);
hold all
surf(sx,sy,sz,ElmtsPlot,'EdgeColor','none');%,'FaceColor','interp' );
hold off
xlabel('x-Axis [m]');
ylabel('y-Axis [m]');
zlabel('z-Axis [m]');
view([-7 25]);

axis tight equal
font_size = 12;
set(gca,'fontsize',font_size)

end

end
