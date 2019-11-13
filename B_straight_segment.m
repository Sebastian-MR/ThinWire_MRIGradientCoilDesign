% 2019-11
% Sebastian Littin
% sebastian.littin@uniklinik-freiburg.de


function B = B_straight_segment(Pstart, Pend, Points)
% calculated field of a straight wire defined by
% its start and end at the locations given by points
% units: input in m, output in T/A

% simple scale: unit_scale=1/4/pi; 
% in SI shall be mu_0/(4*pi)
mu_0=4*pi*1e-7;
unit_scale=mu_0/4/pi; % SI units
epsilon=1e-12;

if vnorm(Pstart-Pend)<epsilon
    B=[0,0,0];
else
    a=Pend(ones(1,size(Points,1)),:)-Points;
    b=Pstart(ones(1,size(Points,1)),:)-Points;
    c=Pend(ones(1,size(Points,1)),:)-Pstart(ones(1,size(Points,1)),:);

    an=normalize(a);
    bn=normalize(b);
    cn=normalize(c);

    d=vnorm(vcross(a,cn));
    IndNon0=find(d>epsilon);
    inv_d=zeros(size(d));
    inv_d(IndNon0)=1./d(IndNon0);

    Babs=unit_scale*(vdot(an,cn)-vdot(bn,cn)).*inv_d;

    Dir_B = vcross(an,cn);
    LenDir_B=vnorm(Dir_B);

    % fix zero norm cases
    NormDir_B=zeros(size(Dir_B));
    NormDir_B(IndNon0,:)=Dir_B(IndNon0,:)./LenDir_B(IndNon0,ones(1,3)); 

    B=NormDir_B.*Babs(:,ones(1,3));
end

% "include" vectors.m
function vout=vnorm(vin) 
% norm of the vector (assuming the second dimension to be the vector coordinates)
vout=sqrt(dot(vin,vin,2));

function vout=vdot(v1,v2)
% dot product of two vectors (assuming the second dimension to be the vector coordinates)
vout=dot(v1,v2,2);

function vout=vcross(v1,v2)
% cross product of two vectors (assuming the second dimension to be the vector coordinates)
vout=cross(v1,v2,2);

function vout=normalize(vin) 
% normalized vector (assuming the second dimension to be the vector coordinates)
vn=sqrt(dot(vin,vin,2));
vout=vin./vn(:,ones(1,size(vin,2)));
