function PlotObs( Observs, hits )
%PLOTOBS Plot observations

global Par;

T = length(Observs);

figure, hold on

if Par.FLAG_ObsMod == 0
    
    xlim([-Par.Xmax Par.Xmax]), ylim([-Par.Xmax Par.Xmax])
    
    for t = 1:T
        plot(Observs(t).r(:, 1), Observs(t).r(:, 2), 'xr');
    end
    for t = 1:T
        ind = hits(t, hits(t,:)~=0);
        plot(Observs(t).r(ind, 1), Observs(t).r(ind, 2), 'xb');
    end
    
elseif Par.FLAG_ObsMod == 1
    
    xlim([-Par.Xmax Par.Xmax]), ylim([-Par.Xmax Par.Xmax])
    
    for t = 1:T
        plot(Observs(t).r(:, 2).*cos(Observs(t).r(:, 1)), Observs(t).r(:, 2).*sin(Observs(t).r(:, 1)), 'xr');
        ind = hits(t, hits(t,:)~=0);
        plot(Observs(t).r(ind, 2).*cos(Observs(t).r(ind, 1)), Observs(t).r(ind, 2).*sin(Observs(t).r(ind, 1)), 'xb');
    end
    for t = 1:T
        ind = hits(t, hits(t,:)~=0);
        plot(Observs(t).r(ind, 2).*cos(Observs(t).r(ind, 1)), Observs(t).r(ind, 2).*sin(Observs(t).r(ind, 1)), 'xb');
    end
    
    figure, hold on
    xlim([-pi, pi]), ylim([0 Par.Xmax]);
    for t = 1:T
        plot(Observs(t).r(:, 1), Observs(t).r(:, 2), 'xr');
    end    
end

pause(0.1);

end