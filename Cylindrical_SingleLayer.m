% Find Wire Patterns for a given target fields using Streamlines of a
% simulated thin wire approximation.
% Works perfectly fine for simple geometries. However, no generalization 
% for arbitrary surfaces
% Example code for a single layer cylindrical geometry
% Units are SI: Meter, Ampere, Tesla, etc.
% 
% 2019-11
% Sebastian Littin

clear all
close all


%% coil description: Cylindrical unshielded coil

plot_all = 1; % optionally plot all intermediate step

% define coil-parameters of the matrix coil: segments_angular, half_length, len_step
CoilDefinition.Partitions = 1;
segments_angular=56;
half_length=0.75; % 500mm
len_step = 0.025; % 20mm
r_coil = 0.35;  % 700mm coil diameter


arc_angle = 360/(segments_angular);
[elm_angle, elm_z] = ndgrid((0:segments_angular-1)*arc_angle, (-half_length:len_step:half_length)); 
CoilDefinition(1).num_elements=size(elm_angle);
elm_angle_shift = elm_angle([2:end,1],:);


% Define Cylindrical Main Surface
CoilDefinition(1).thin_wire_nodes_start = [cosd(elm_angle(:))*r_coil,sind(elm_angle(:))*r_coil,elm_z(:)];
CoilDefinition(1).thin_wire_nodes_stop = [cosd(elm_angle_shift(:))*r_coil,sind(elm_angle_shift(:))*r_coil,elm_z(:)];

CoilDefinition(1).num_elements = size(elm_angle);




% possibility to plot thin wire elements
if plot_all == 1
figure;
hold all
for n = 1:length(CoilDefinition(1).thin_wire_nodes_start)
plot3([CoilDefinition(1).thin_wire_nodes_start(n,1) CoilDefinition(1).thin_wire_nodes_stop(n,1)], ...
    [CoilDefinition(1).thin_wire_nodes_start(n,2) CoilDefinition(1).thin_wire_nodes_stop(n,2)],...
    [CoilDefinition(1).thin_wire_nodes_start(n,3) CoilDefinition(1).thin_wire_nodes_stop(n,3)])
    
end
hold off
axis equal tight
title('Thin-wire current elements');
view([1 1 1])
end


% Some definitions for 3D contour plotting...
CoilDefinition(1).Radius = r_coil;

CoilDefinition(1).Length = half_length*2;


% Definition of main target points in a 3D-volume

TargetDefinition.shape = 'sphere';
TargetDefinition.radius = 0.15;
TargetDefinition.length = 0.3;
TargetDefinition.resol_radial = 3;
TargetDefinition.resol_angular = 32;
TargetDefinition.resol_length = 8;
TargetDefinition.strength = 5e-3;
TargetDefinition.direction = 'x';

target_points = Make_Target(TargetDefinition);

% plot target
if plot_all == 1
figure; scatter3(target_points.points.x1(:), target_points.points.x2(:), target_points.points.x3(:), ones(size(target_points.points.x1(:)))*25, target_points.field(:))
axis equal tight
title('Target Points and Field');
view([1 1 1])
end

x1 = target_points.points.x1(:);
x2 = target_points.points.x2(:);
x3 = target_points.points.x3(:);

Points=[x1(:),x2(:),x3(:)];
Target.Points=Points;
num_points=length(Points(:,1));
Target.num_points = num_points;

kn = length(x1)^2;
kp = length(x1);



%% Calculate regularized Thin-wire solution
CoilDefinition(1).StreamDirection = 2;

Sensitivity = ThinWireSensitivity(CoilDefinition, Target);

%% Calculate an unregularized Solution

btarget = target_points.field;

ElementCurrents_Unreg(1).Stream = pinv(Sensitivity(1).ElementFieldsStream(:,:))*btarget;
% ResultingField = ElementFields*ElementCurrents;

%% Plot the unregularized solution

figure; imab(reshape(ElementCurrents_Unreg(1).Stream,size(elm_angle)-[0 1])); colorbar; title('Unregularized Stream Function');

if plot_all == 1
PlotThinWireStreamFunction3D(CoilDefinition, ElementCurrents_Unreg)
ContourPlotThinWireStreamFunction3D(CoilDefinition, ElementCurrents_Unreg, 13)
end

%% Calculate the regularized Solution

lambda=1;

E_Mat = Sensitivity(1).ElementFieldsStream(:,:);

ElementCurrents_Reg=TikhonovReg(E_Mat, btarget, 0.0077); % regularisation automatically penelizes total power


%% Plot currents in 2D


ElementCurrentsReg(1).Stream = reshape(ElementCurrents_Reg,size(elm_angle)-[0 1]);

figure; imab(ElementCurrentsReg(1).Stream); colorbar; title('Regularized Stream Function');

if plot_all == 1
PlotThinWireStreamFunction3D(CoilDefinition, ElementCurrentsReg)
ContourPlotThinWireStreamFunction3D(CoilDefinition, ElementCurrentsReg, 19)
end

