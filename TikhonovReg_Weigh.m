% calculates Currents approximately satisfying equation
% TargetField = ElementFields*Currents using simple Tikhonov regularization
% with the unity regularization matrix. Regularization parameter lambda is
% scaled to the norm of ElementFields
%
% 2019-11
% Sebastian Littin
% sebastian.littin@uniklinik-freiburg.de


function Currents = TikhonovReg_Weigh(ElementFields, TargetField, lambda, w)

AtA=ElementFields'*ElementFields;
G=lambda*norm(ElementFields)*w;
GtG=G'*G;
Currents=pinv(AtA+GtG)*(ElementFields')*TargetField;


