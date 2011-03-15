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
for c = 1:Par.NumTgts

% Loop through particles
for ii = 1:Par.NumPart
    
    % Loop through targets
    for j = 1:Distn.clusters{c}.particles{ii}.N
        
        % Choose a colour
        col = [rand, rand, 0];
        
        % create an array
        num = Distn.clusters{c}.particles{ii}.tracks{j}.num;
        x = zeros(num,1);
        y = zeros(num,1);
        
        % Collate state
        for k = 1:num
            state = Distn.clusters{c}.particles{ii}.tracks{j}.state{k};
            x(k) = state(1);
            y(k) = state(2);
        end
            
        % Plot track
        plot(x, y, '-', 'color', col);
        
    end

end

end

plot(0, 0, 'xk');

end