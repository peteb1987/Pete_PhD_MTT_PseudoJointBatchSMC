function f = PlotTrueState( Tracks )
%PLOTTRUESTATE Plots the correct states

global Par;

% Make a figure of the right size
f = figure; hold on
xlim([-Par.Xmax Par.Xmax]), ylim([-Par.Xmax Par.Xmax])

% Loop through tracks
for j = 1:Par.NumTgts

    coords = zeros(Tracks{j}.num, 2);
    
    % Compile a list of positions
    coords(:, 1) = cellfun(@(x)(x(1)), Tracks{j}.state);
    coords(:, 2) = cellfun(@(x)(x(2)), Tracks{j}.state);
    
    % Plot it
    plot(coords(:,1), coords(:,2), '-x', 'color', [rand rand 0]);
    
end

pause(0.2);

end