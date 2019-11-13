% A Sensitivity Matrix for partitions of thin wires is calculated and a 
% stream function is derived in this part
% 
% 2019-11
% Sebastian Littin
% sebastian.littin@uniklinik-freiburg.de

function TWS = ThinWireSensitivity(CoilDefinition, Target);

for nP = 1:CoilDefinition(1).Partitions; % repeat for n partitions
disp(['Partition ' num2str(nP)]);
    
ElementFields=zeros(Target.num_points, (CoilDefinition(nP).num_elements(1)*CoilDefinition(nP).num_elements(2)));
for e=1:(CoilDefinition(nP).num_elements(1)*CoilDefinition(nP).num_elements(2));
    Blead1=B_straight_segment(CoilDefinition(nP).thin_wire_nodes_start(e,:),CoilDefinition(nP).thin_wire_nodes_stop(e,:),Target.Points);

    ElementFields(:,e)=Blead1(:,3);
    
    if floor(e*10/(CoilDefinition(nP).num_elements(1)*CoilDefinition(nP).num_elements(2)))~=floor((e-1)*10/(CoilDefinition(nP).num_elements(1)*CoilDefinition(nP).num_elements(2)))
        disp([num2str(e) ' of ' num2str((CoilDefinition(nP).num_elements(1)*CoilDefinition(nP).num_elements(2))) ' done ']);
    end
end

ElementFieldsRe = reshape(ElementFields,[Target.num_points,CoilDefinition(nP).num_elements]);
TWS(nP).ElementFields = ElementFields;


% Define direction of performing the difference to get a stream function
if CoilDefinition(nP).StreamDirection == 1;
    ElementFieldsStreamRe = ElementFieldsRe(:,1:end-1,:) - ElementFieldsRe(:,2:end,:);
elseif CoilDefinition(nP).StreamDirection == 2;
    ElementFieldsStreamRe = ElementFieldsRe(:,:,1:end-1) - ElementFieldsRe(:,:,2:end);
% elseif CoilDefinition(nP).StreamDirection == 3;
%     ElementFieldsStreamRe = ElementFieldsRe(:,:,1:end-1) - ElementFieldsRe(:,:,2:end);
else
    ElementFieldsStreamRe = zeros(size(ElementFieldsRe));
end

TWS(nP).ElementFieldsStream = ElementFieldsStreamRe(:,:);

end % end partitions

end % end function