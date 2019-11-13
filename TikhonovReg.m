% calculates Currents approximately satisfying equation
% TargetField = ElementFields*Currents using simple Tikhonov regularization
% with the unity regularization matrix. Regularization parameter lambda is
% scaled to the norm of ElementFields
%
% 2019-11
% Sebastian Littin
% sebastian.littin@uniklinik-freiburg.de


function Currents = TikhonovReg(ElementFields, TargetField, lambda)

AtA=ElementFields'*ElementFields;
G=lambda*norm(ElementFields)*eye(size(ElementFields,2));
GtG=G'*G;
Currents=pinv(AtA+GtG)*(ElementFields')*TargetField;


