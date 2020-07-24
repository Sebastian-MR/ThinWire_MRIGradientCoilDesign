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

plot_all = 1; % set to 1, to optionally plot intermediate steps

% define coil-parameters of the matrix coil: segments_angular, half_length, len_step
CoilDefinition.Partitions = 1;
segments_angular=56;
half_length=0.45; % 450mm
len_step = 0.015; % 15mm
r_coil = 0.4;  % 800mm coil diameter


arc_angle = 360/(segments_angular);
[elm_angle, elm_z] = ndgrid((0:segments_angular-1)*arc_angle, (-half_length:len_step:half_length)); 
CoilDefinition(1).num_elements=size(elm_angle);
elm_angle_shift = elm_angle([2:end,1],:);


% Define Cylindrical Surface
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
TargetDefinition.radius = 0.2;
TargetDefinition.resol_radial = 3;
TargetDefinition.resol_angular = 15;
TargetDefinition.resol_length = 12;
TargetDefinition.strength = 5e-3;
TargetDefinition.direction = 'y';

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



%% Calculate sensitivity matrix
CoilDefinition(1).StreamDirection = 2;

Sensitivity = ThinWireSensitivity(CoilDefinition, Target);

%% Calculate an unregularized Solution

btarget = target_points.field;

ElementCurrents_Unreg(1).Stream = pinv(Sensitivity(1).ElementFieldsStream(:,:))*btarget;

%% Plot the unregularized solution

if plot_all == 1
figure; imab(reshape(ElementCurrents_Unreg(1).Stream,size(elm_angle)-[0 1])); colorbar; title('Unregularized Stream Function');
PlotThinWireStreamFunction3D(CoilDefinition, ElementCurrents_Unreg)
% ContourPlotThinWireStreamFunction3D(CoilDefinition, ElementCurrents_Unreg, 13)
end

%% Calculate the (wrong) regularized Solution


btarget = target_points.field;
E_Mat = Sensitivity(1).ElementFieldsStream(:,:);

ElementCurrents_Reg=TikhonovReg(E_Mat, btarget, 0.2077); % regularisation automatically penelizes total power


figure; imab(reshape(ElementCurrents_Reg,size(elm_angle)-[0 1])); colorbar;

%% Use alternative regularization matrix for regularizing effective currents

E_Mat = Sensitivity(1).ElementFieldsStream(:,:);
W = eye(size(E_Mat,2));

w = W - circshift(W,segments_angular,2);
w(1:2*segments_angular,end-2*segments_angular+1:end) = zeros(2*segments_angular);


ElementCurrents_Reg=TikhonovReg_Weigh(E_Mat, btarget, 5e-1, w); 

figure; imab(reshape(ElementCurrents_Reg,size(elm_angle)-[0 1])); colorbar;



%% Add additional constraints to enforce peripheral elements to be 0


btarget = target_points.field;

lambda = 1e1;

E_Mat = Sensitivity(1).ElementFieldsStream(:,:);
W = eye(size(E_Mat,2));

w = W - circshift(W,segments_angular,2);
w(1:2*segments_angular,end-2*segments_angular+1:end) = zeros(2*segments_angular);


w_ext = [w; [lambda*eye(segments_angular) zeros(segments_angular,size(E_Mat,2)-segments_angular)];...
        [zeros(segments_angular,size(E_Mat,2)-segments_angular)  lambda*eye(segments_angular) ]];
    
   
btarget_ext = [target_points.field; zeros(CoilDefinition.Partitions*2,1)];

ElementCurrents_Reg = TikhonovReg_Weigh(E_Mat, btarget, 5e-1, w_ext); 

figure; imab(reshape(ElementCurrents_Reg,size(elm_angle)-[0 1])); colorbar;



%%
btarget = target_points.field;

red_count1 = 1:CoilDefinition(1).num_elements(1);
red_count2 = (CoilDefinition(1).num_elements(1)*(CoilDefinition(1).num_elements(2)-2)+1):(CoilDefinition(1).num_elements(1)*(CoilDefinition(1).num_elements(2)-1));
red_count = [red_count1 red_count2];
red_countd = [circshift(red_count1,1) circshift(red_count2,1)];

diff_constraint = Sensitivity(1).ElementFieldsStream(:,red_count)-Sensitivity(1).ElementFieldsStream(:,red_countd);

E_Mat = Sensitivity(1).ElementFieldsStream(:,:);

E_Mat_ext = [E_Mat 5e3*diff_constraint];
btarget_ext = [btarget; zeros(size(diff_constraint,2),1)];

n_stop = CoilDefinition(1).num_elements(1)*(CoilDefinition(1).num_elements(2)-1);


W = eye(size(E_Mat,2));

w = W - circshift(W,segments_angular,2);
w(1:2*segments_angular,end-2*segments_angular+1:end) = zeros(2*segments_angular);

w_ext = [w zeros(size(E_Mat,2),2*segments_angular);
        zeros(2*segments_angular,size(E_Mat,2)) eye(2*segments_angular)];

ElementCurrents_temp=TikhonovReg_Weigh(E_Mat_ext, btarget, 5e-1, w_ext); 

n_stop = CoilDefinition(1).num_elements(1)*(CoilDefinition(1).num_elements(2)-1);
ElementCurrents_Reg = ElementCurrents_temp(1:n_stop);

figure; imab(reshape(ElementCurrents_Reg,size(elm_angle)-[0 1])); colorbar;


%% Plot currents in 2D


ElementCurrentsReg(1).Stream = reshape(ElementCurrents_Reg,size(elm_angle)-[0 1]);


if plot_all == 1
figure; imab(ElementCurrentsReg(1).Stream); colorbar; title('Regularized Stream Function');
end

% PlotThinWireStreamFunction3D(CoilDefinition, ElementCurrentsReg)
ContourPlotThinWireStreamFunction3D(CoilDefinition, ElementCurrentsReg,19)


