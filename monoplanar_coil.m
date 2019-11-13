% Find Wire Patterns for a given target fields using Streamlines of a
% simulated thin wire approximation.
% Works perfectly fine for simple geometries. However, no generalization 
% for arbitrary surfaces
% Example code for a shielded cylindrical geometry
% Units are SI: Meter, Ampere, Tesla, etc.
%
% 2019-11
% Sebastian Littin
% sebastian.littin@uniklinik-freiburg.de


clear all
close all


%% coil description: Cylindrical unshielded coil


%define coil-parameters
CoilDefinition(1).Partitions = 1;

length_z=1.2;
len_step_z = 0.025;
length_x=1.2;
len_step_x = 0.025; % 100mm => step-length-z 4.5cm
distance_y = -0.2; % distance to iso-center


[elm_x, elm_z] = ndgrid(-length_x/2:len_step_x:length_x/2, -length_z/2:len_step_z:length_z/2); 
num_elements=length(elm_x(:));



% Define Cylindrical Main Surface
CoilDefinition(1).thin_wire_nodes_start = [elm_x(:)-len_step_x/2,-ones(size(elm_x(:)))*distance_y,elm_z(:)];
CoilDefinition(1).thin_wire_nodes_stop = [elm_x(:)+len_step_x/2,-ones(size(elm_x(:)))*distance_y,elm_z(:)];
CoilDefinition(1).num_elements = size(elm_x);



%% Definition of target points in a 3D-volume

TargetDefinition.shape = 'sphere';
TargetDefinition.radius = 0.15;
TargetDefinition.resol_radial = 3;
TargetDefinition.resol_angular = 15;
TargetDefinition.strength = 5e-3;
TargetDefinition.direction = 'x';

target_main = Make_Target(TargetDefinition);

% % optionally plot main target
% figure; scatter3(target_main.points.x1(:), target_main.points.x2(:), target_main.points.x3(:), ones(size(target_main.points.x1(:)))*25, target_main.field(:))
% axis equal tight


x1 = [target_main.points.x1(:)];
x2 = [target_main.points.x2(:)];
x3 = [target_main.points.x3(:)];

Points=[x1(:),x2(:),x3(:)];
Target.Points=Points;
num_points=length(Points(:,1));
Target.num_points = num_points;

kn = length(x1)^2;
kp = length(x1);

num_points_main=length(target_main.points.x1);
% num_points_shield=length(target_shield.points.x1);


%% Calculate regularized Thin-wire solution
CoilDefinition(1).StreamDirection = 2;

Sensitivity = ThinWireSensitivity(CoilDefinition, Target);

%% Calculate the unregularized Solution

E_Mat = [Sensitivity(1).ElementFieldsStream];
btarget = [target_main.field(:)];

ElementCurrents_temp = pinv(E_Mat)*btarget;


%%

ElementCurrents(1).Stream = reshape(ElementCurrents_temp,size(elm_x)-[0 1]);

figure; set(gcf,'Name','3D coil','Position',[   1   1   1000   500]);
imab(ElementCurrents(1).Stream'); colorbar; title('a) coil currents surface 1');


% PlotThinWireStreamFunction3D(CoilDefinition, ElementCurrents)
%% Calculate the regularized Solution

ElementCurrents_temp=TikhonovReg(E_Mat, btarget, 0.0077); % regularisation automatically penelizes total power


%% Plot currents in 2D


ElementCurrents(1).Stream = reshape(ElementCurrents_temp,size(elm_x)-[0 1]);

figure; set(gcf,'Name','3D coil','Position',[   1   1   1000   500]);

n_cont = 9;
hold all
ElmtsPlot = [reshape(ElementCurrents(1).Stream,(CoilDefinition(1).num_elements -[0 1]))];
imab(ElementCurrents(1).Stream'); colorbar; title('coil currents Main');
contour(ElmtsPlot(:,:),n_cont,'k','LineWidth', 2);
hold off


%%
PlotThinWireStreamFunction3D(CoilDefinition, ElementCurrents)


% ContourPlotThinWireStreamFunction3D(CoilDefinition, ElementCurrents, 13)


%% Plot multi layer contours



nP = 1;

PlotCoord = (CoilDefinition(nP).thin_wire_nodes_start + CoilDefinition(nP).thin_wire_nodes_stop)/2;

n_cont = 15;

ElmtsPlot = [reshape(ElementCurrents(nP).Stream,(CoilDefinition(nP).num_elements -[0 1]))];
cont_max_main = max(max(ElmtsPlot));
[C1,H1] = contour(ElmtsPlot(:,:)',[-cont_max_main:(2*cont_max_main/n_cont):cont_max_main],'k','LineWidth', 2);

nP = 2;

PlotCoord = (CoilDefinition(nP).thin_wire_nodes_start + CoilDefinition(nP).thin_wire_nodes_stop)/2;

n_cont = 15;

ElmtsPlot = [reshape(ElementCurrents(nP).Stream,(CoilDefinition(nP).num_elements -[0 1]))];
cont_max_main = max(max(ElmtsPlot));
[C2,H2] = contour(ElmtsPlot(:,:)',[-cont_max_main:(2*cont_max_main/n_cont):cont_max_main],'k','LineWidth', 2);


%% 3D Plot of the contours

figure; set(gcf,'Name','3D coil','Position',[   1   1   1000   1000]);
hold all

S = contourdata(C1);
ccount = size(S);
nP =1;
for i = 1:ccount(2)
    sx=CoilDefinition(nP).Radius*cos(S(i).xdata/CoilDefinition(nP).num_elements(1)*2*pi);
    sy=CoilDefinition(nP).Radius.*sin(S(i).xdata/CoilDefinition(nP).num_elements(1)*2*pi);
    sz=S(i).ydata./length(ElmtsPlot(1,:)')*CoilDefinition(nP).Length - CoilDefinition(nP).Length/2;
    
    plot3(sx,sy,sz,'b','LineWidth', 1)
end

S = contourdata(C2);
ccount = size(S);
nP =2;
for i = 1:ccount(2)
    sx=CoilDefinition(nP).Radius*cos(S(i).xdata/CoilDefinition(nP).num_elements(1)*2*pi);
    sy=CoilDefinition(nP).Radius.*sin(S(i).xdata/CoilDefinition(nP).num_elements(1)*2*pi);
    sz=S(i).ydata./length(ElmtsPlot(1,:)')*CoilDefinition(nP).Length - CoilDefinition(nP).Length/2;
    
    plot3(sx,sy,sz,'r','LineWidth', 1)
end

hold off
view([-7 25]);

axis tight equal
font_size = 12;
set(gca,'fontsize',font_size)

