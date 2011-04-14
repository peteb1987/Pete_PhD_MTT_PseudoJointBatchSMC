function PlotTracks( Distn, f )
%PLOTTRACKS Plot the output of the batch SMC multi-target tracker

global Par;

if nargin == 1
    % Create a window
    figure, hold on
    xlim([-Par.Xmax Par.Xmax]), ylim([-Par.Xmax Par.Xmax])

else
    figure(f)
end
    
% Loop through clusters
for c = 1:Distn.N

cellfun(@PlotParticle, Distn.clusters{c}.particles);

end

plot(0, 0, 'xk');

end



function PlotParticle(Part)

% Loop through targets
for j = 1:Part.N
    
    if Part.tracks{j}.num > 0
        
        % Choose a colour
        col = [rand, rand, 0];
        
        % Collate state
        state = cell2mat(Part.tracks{j}.state');
        x = state(1, :);
        y = state(2, :);
        
        % Plot track
        plot(x, y, '-', 'color', col);
        
    end
    
end

end

