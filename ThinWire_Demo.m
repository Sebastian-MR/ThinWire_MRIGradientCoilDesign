%  This code demonstrates the use of thin wires to approximate a current
% density. In gradient coil design for MRI usually only the the z-component
% of the magnetic field (Bz) is considered. Hence, only wires orthogonal to
% z are used in this simulation. A sensitivity matrix is used to calculate
% acurrent distribution. A regularization and an additional constraint is
% deployed to derive a ralizable coil design.
% This ,ethod may be used for simple geometries. However, no generalization
% for arbitrary surfaces.
%
% 2019-11
% Sebastian Littin
% sebastian.littin@uniklinik-freiburg.de

clear all
close all


%% coil description: Cylindrical unshielded coil

% define coil-parameters of the matrix coil: segments_angular, half_length, len_step

% segments_angular=56;
% half_length=0.6; % 500mm
% len_step = 0.02; % 20mm
% r_coil = 0.3;  % 500mm coil diameter
% circ_gap = 0.00; % 5mm (gap between the coils in the angular direction)
%
% arc_angle = 360/(segments_angular-1);
% [elm_angle, elm_z] = ndgrid((-0.5:segments_angular-1.5)*arc_angle, (-half_length:len_step:half_length));
% num_elements=length(elm_angle(:));
% elm_angle_shift = elm_angle([2:end,1],:);
%
% thin_wire_nodes_start = [cosd(elm_angle(:))*r_coil,sind(elm_angle(:))*r_coil,elm_z(:)];
% thin_wire_nodes_stop = [cosd(elm_angle_shift(:))*r_coil,sind(elm_angle_shift(:))*r_coil,elm_z(:)];
%

% define coil-parameters of the matrix coil: segments_angular, half_length, len_step
CoilDefinition.Partitions = 1;
segments_angular=56;
half_length=0.6; % 500mm
len_step = 0.02; % 20mm
r_coil = 0.3;  % 700mm coil diameter


arc_angle = 360/(segments_angular);
[elm_angle, elm_z] = ndgrid((-0:segments_angular-0)*arc_angle, (-half_length:len_step:half_length));
CoilDefinition(1).num_elements=size(elm_angle);
elm_angle_shift = elm_angle([2:end,1],:);


% Define Cylindrical Main Surface
CoilDefinition(1).thin_wire_nodes_start = [cosd(elm_angle(:))*r_coil,sind(elm_angle(:))*r_coil,elm_z(:)];
CoilDefinition(1).thin_wire_nodes_stop = [cosd(elm_angle_shift(:))*r_coil,sind(elm_angle_shift(:))*r_coil,elm_z(:)];

CoilDefinition(1).num_elements = size(elm_angle);




%% Definition of target points in a 3D-volume

TargetDefinition.shape = 'sphere';
TargetDefinition.radius = 0.15;
% TargetDefinition.length = 0.25;
TargetDefinition.resol_radial = 3;
TargetDefinition.resol_angular = 8;
% TargetDefinition.resolution = 12;
TargetDefinition.strength = 5e-3;
TargetDefinition.direction = 'y';

target_points = Make_Target(TargetDefinition);

% plot target
figure; scatter3(target_points.points.x1(:), target_points.points.x2(:), target_points.points.x3(:), ones(size(target_points.points.x1(:)))*25, target_points.field(:))
axis equal tight

x1 = target_points.points.x1(:);
x2 = target_points.points.x2(:);
x3 = target_points.points.x3(:);

Points=[x1(:),x2(:),x3(:)];
Target.Points=Points;
num_points=length(Points(:,1));
Target.num_points = num_points;

kn = length(x1)^2;
kp = length(x1);

%% plot the thin wire elements

figure;
hold all
for n = 1:length(CoilDefinition(1).thin_wire_nodes_start)
    plot3([CoilDefinition(1).thin_wire_nodes_start(n,1) CoilDefinition(1).thin_wire_nodes_stop(n,1)], ...
        [CoilDefinition(1).thin_wire_nodes_start(n,2) CoilDefinition(1).thin_wire_nodes_stop(n,2)],...
        [CoilDefinition(1).thin_wire_nodes_start(n,3) CoilDefinition(1).thin_wire_nodes_stop(n,3)])
    
end
hold off
axis equal tight

% Some definitions for 3D contour plotting...
CoilDefinition(1).Radius = r_coil;
CoilDefinition(1).Length = half_length*2;


%% Calculate regularized Thin-wire solution
CoilDefinition(1).StreamDirection = 2;

Sensitivity = ThinWireSensitivity(CoilDefinition, Target);


%% Calculate a unregularized solution using Moore-Penrose pseudoinverse (pinv)
btarget = target_points.field;
ElementCurrents = pinv(Sensitivity.ElementFields)*btarget;
ResultingField = Sensitivity.ElementFields*ElementCurrents;


figure; set(gcf,'Name','Currents unreg','Position',[   1   1   500   500]);
imab(reshape(ElementCurrents,size(elm_angle)));title('unregularized current distribution');
ylabel('z-Axis [m]')
xlabel('circumferential [rad]')


%% Calculate simple Tikhonov-regularized solution
% Unfortunately this results in a crrent distribution which is not 
% realizable with closed loops...

ElementCurrents=TikhonovReg(Sensitivity.ElementFields, btarget, 0.0077); % regularisation automatically penelizes total power
figure; imab(reshape(ElementCurrents,size(elm_angle))); colorbar; title('Regularized current distribution');


%% Add additional constraints to the regularisation
% add conditions for balancing along each angular column
% => This ensures that the sum of currents along z=0 -> realizable as
% closed loops. As many equations as angular segments are needed

eye_ang=eye(size(elm_angle,1));
ElementFields_Add3D=eye_ang(:,:,ones(1,size(elm_angle,2)));
ElementFields_Add=ElementFields_Add3D(:,:);
TargetFields_Add=zeros(size(elm_angle,1),1);

ElementFields_Balance = [Sensitivity.ElementFields; ElementFields_Add*5e-4]; % this coefficient controls the relative importance of the new equations
TargetField_Balance = [btarget; TargetFields_Add];

% Do regularization
ElementCurrents_Balance = TikhonovReg(ElementFields_Balance, TargetField_Balance, 0.00005); % << ~5% deviation, this time regularization is very different and depends on the above coefficient

% Calculate resulting field
ResultingField_Balance = Sensitivity.ElementFields*ElementCurrents_Balance;

figure; set(gcf,'Name','Currents reg','Position',[   1   1   500   500]);ElementFields_Balance = [Sensitivity.ElementFields; ElementFields_Add*5e-4]; % this coefficient controls the relative importance of the new equations
TargetField_Balance = [btarget; TargetFields_Add];
imab(reshape(ElementCurrents_Balance,size(elm_angle)));title('coil currents with balance');
ylabel('z-Axis [m]')
xlabel('circumferential [rad]')

ElementCurrents_Balance_reshape = reshape(ElementCurrents_Balance,size(elm_angle));
%%
% Calculate the Stream function from the current distribution, using the
% cumulative sum (=integral) along z. To balance unsymmetric solutions the
% interation is performed from both sides...

Stream_Reg=cumsum(ElementCurrents_Balance_reshape,2,'forward');
Stream_Reg_rev=cumsum(ElementCurrents_Balance_reshape,2,'reverse');

Stream = zeros(size(Stream_Reg)+[0 1]);
Stream(:,2:end) = Stream_Reg./2;
Stream(:,1:end-1) = Stream(:,1:end-1)-Stream_Reg_rev./2;

% Stream = -Stream_Reg_rev;

%% Plot the stream function and stream lines in 2D...

figure; set(gcf,'Name','Stream Function','Position',[   1   1   500   500]);

hold all
imab(Stream);
cont_max = 4500;
n_cont = 17;
[C,H] = contour(Stream(:,:)',[-cont_max:(2*cont_max/n_cont):cont_max],'k','LineWidth', 2);
title('Stream function from integration');
hold off

ylabel('z-Axis [m]')
xlabel('circumferential [rad]')

%% 2D Plot


figure; set(gcf,'Name','Stream Function','Position',[   1   1   1450   400]);

subplot(1,3,1)
imab(reshape(ElementCurrents,size(elm_angle))); colorbar; title('a) currents regularized');
ylabel('z-Axis [m]')
xlabel('circumferential [rad]')
subplot(1,3,2)
imab(reshape(ElementCurrents_Balance,size(elm_angle)));title('b) regularized w constraints');colorbar;
ylabel('z-Axis [m]')
xlabel('circumferential [rad]')

subplot(1,3,3)
hold all
imab(Stream);
cont_max = 4500;
n_cont = 13;
[C,H] = contour(Stream(:,:)',[-cont_max:(2*cont_max/n_cont):cont_max],'k','LineWidth', 2);
colorbar; title('c) Stream Function');
hold off
ylabel('z-Axis [m]')
xlabel('circumferential [rad]')

