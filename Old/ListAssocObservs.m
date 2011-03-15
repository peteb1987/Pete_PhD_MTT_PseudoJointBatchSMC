function [ obs ] = ListAssocObservs( t, L, Target, Observs )
%LISTASSOCOBSERVS Draw up a list of observations associated with Target
%between t-L+1 and t.

obs = cell(L, 1);

% Loop through window
for k = 1:L
    
    % Get the association index
    ass = Target.GetAssoc(t-L+k);
    
    % Get observation if one is associated
    if ass ~= 0
        obs{k} = Observs(t-L+k).r(ass, :)';
    end
    
end

end

