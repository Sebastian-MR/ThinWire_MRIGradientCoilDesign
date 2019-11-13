% Generate target points for a gradient coil optimization
% Possible shape options are 'cube', 'cuboid', 'sphere', 'cylinder' and 
% 'ellipsoid'.
% 
% Following definitions are required:
% 'cube': skalars for length and resolution
% 'cuboid': 3 dimensional vector for length and resolution
% 'sphere': skalar for radius, resolution in angular and radial direction
% 'cylinder': skalar for radius and length along z; resolution in angular,
% radial and z- direction
% 
% A target field may be defined by specifying a direction ('x', 'y' or 
% 'z') and a gradient strength. If no target gradient is specified, 1mT/m 
% along x is given out.

function B_out = Make_Target(TargetDefinition)

% define different shapes

if isequal(TargetDefinition.shape, 'cube')
    r1 = TargetDefinition.length/2;
    d1 = TargetDefinition.length/(TargetDefinition.resolution-1);
    
    [x1,x2,x3] = ndgrid(-r1:d1:r1, -r1:d1:r1, -r1:d1:r1); 

elseif isequal(TargetDefinition.shape, 'cuboid')
    r1 = TargetDefinition.length(1)/2;
    r2 = TargetDefinition.length(2)/2;    
    r3 = TargetDefinition.length(3)/2;        
    d1 = TargetDefinition.length(1)/(TargetDefinition.resolution(1)-1);
    d2 = TargetDefinition.length(2)/(TargetDefinition.resolution(2)-1);    
    d3 = TargetDefinition.length(3)/(TargetDefinition.resolution(3)-1);
    
    [x1,x2,x3] = ndgrid(-r1:d1:r1, -r2:d2:r2, -r3:d:r3); 

elseif isequal(TargetDefinition.shape, 'sphere')
    r = TargetDefinition.radius;
    d1 = TargetDefinition.radius/(TargetDefinition.resol_radial-1);
    d2 = pi/(TargetDefinition.resol_angular-1);
    d3 = 2*pi/(TargetDefinition.resol_angular-1);
    [ra,theta,phi] = ndgrid(0:d1:r, 0:d2:pi, -pi:d3:pi);   
    
    x1 = ra.*sin(theta).*cos(phi);
    x2 = ra.*sin(theta).*sin(phi);
    x3 = ra.*cos(theta);    

elseif isequal(TargetDefinition.shape, 'cylinder')
    r = TargetDefinition.radius;
    d2 = 2*pi/(TargetDefinition.resol_angular-1);
    d3 = TargetDefinition.length/(TargetDefinition.resol_length-1);
    l1 = TargetDefinition.length/2;
    
        if TargetDefinition.resol_radial == 1;
                [ra,phi,z] = ndgrid(r, 0:d2:2*pi, -l1:d3:l1); 
        else
            d1 = TargetDefinition.radius/(TargetDefinition.resol_radial-1);
            [ra,phi,z] = ndgrid(0:d1:r, 0:d2:2*pi, -l1:d3:l1);    
        end
        
    x1 = ra.*cos(phi);
    x2 = ra.*sin(phi);
    x3 = z;     
end

B_out.points.x1 = x1(:);
B_out.points.x2 = x2(:);
B_out.points.x3 = x3(:);

B_field = zeros(size(x1));


% define field strength
if isequal(TargetDefinition.direction, 'x')
    B_field = x1(:)./max(x1(:)).*TargetDefinition.strength;
    
elseif isequal(TargetDefinition.direction, 'y')
    B_field = x2(:)./max(x2(:)).*TargetDefinition.strength;    

elseif isequal(TargetDefinition.direction, 'z')
    B_field = x3(:)./max(x3(:)).*TargetDefinition.strength;
    
else
    B_field = x1(:)./max(x1(:)).*1e-3;    

end

B_out.field = B_field;

end

