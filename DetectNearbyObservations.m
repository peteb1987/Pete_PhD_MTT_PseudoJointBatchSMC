function [ObsTargIndexes] = DetectNearbyObservations(t, Observs, Distn)
%DETECTNEARBYOBSERVATIONS Detect which observations are in the vicinity of
%each cluster. This acts as a crude gate and reduces accurate checking
%later.

%, Area

ObsTargIndexes = {(1:Observs(t).N)'};

% ObsTargIndexes = cell(Distn.N, 1);
% for c = 1:Distn.N
%     ObsTargIndexes{c} = (1:Observs(t).N)';
% end

% global Par;
% 
% ObsTargIndexes = cell(Distn.N, 1);
% % Area = cell(Distn.N, 1);
% 
% % Loop through targets and fetch state
% for c = 1:Distn.N
%     for j = 1:Distn.clusters{c}.N
%         Parts = zeros(Par.NumPart, 4);
%         States = zeros(Distn.clusters{c}.N, 2);
%         Spread = zeros(Distn.clusters{c}.N, 2);
%         PartsPol= zeros(Par.NumPart, 2);
%         StatesPol = zeros(Distn.clusters{c}.N, 2);
%         SpreadPol = zeros(Distn.clusters{c}.N, 2);
%         for ii = 1:Par.NumPart
%             Parts(ii, :) = Distn.clusters{c}.particles{ii}.tracks{j}.GetState(t)';
%             [PartsPol(ii, 1), PartsPol(ii, 2)] = cart2pol(Parts(ii, 1), Parts(ii, 2)); 
%         end
%         Parts(:, 3:4) = [];
%         
%         % Calculate the centre of the target particles and the spread
%         States(j, :) = mean(Parts);
%         Spread(j, :) = range(Parts);
%         
%         StatesPol(j, :) = mean(PartsPol);
%         SpreadPol(j, :) = range(PartsPol);
%         
%         % Find observations within these bounds and add their indexes to the list
%         for i = 1:Observs(t).N
% %             [obs_cart(1), obs_cart(2)] = pol2cart(Observs(t).r(i, 1), Observs(t).r(i, 2));
% %             if (abs(obs_cart(1)-States(j,1))<(Spread(j,1)+Par.Vlimit)) && (abs(obs_cart(2)-States(j,2))<(Spread(j,2)+Par.Vlimit))
% %                 ObsTargIndexes{c} = [ObsTargIndexes{c}; i];
% %             end
%             if StatesPol(j,2)>Par.Vlimit
%                 bng_span = SpreadPol(j,1)+Par.Vlimit/(StatesPol(j,2)-Par.Vlimit);
%             else
%                 bng_span = pi;
%             end
%             if (abs(Observs(t).r(i, 1)-StatesPol(j,1))<bng_span) && (abs(Observs(t).r(i, 2)-StatesPol(j,2))<(SpreadPol(j,2)+Par.Vlimit))
%                 if ~any(i==ObsTargIndexes{c})
%                     ObsTargIndexes{c} = [ObsTargIndexes{c}; i];
%                 end
%             end
%         end
%         
%     end
%     
% %     Area{c} = min(2*bng_span, 2*pi) * ((SpreadPol(j,2)+Par.Vlimit)+min(SpreadPol(j,2)+Par.Vlimit, StatesPol(j,2)));
%     
% end

end

