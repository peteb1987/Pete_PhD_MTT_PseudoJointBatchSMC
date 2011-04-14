function [ass, count, present] = AnalyseAss( correct, Distn, fr)
%ANALYSEASS Compare associations with correct values

global Par;

ass = cell(Distn.N, 1);
count = zeros(Distn.N, fr);
present = zeros(Distn.N, fr);

for c = 1:Distn.N
    
    for j = 1:Distn.clusters{c}.N
        
        jj = Distn.clusters{c}.members(j);
        ass{jj} = zeros(Par.NumPart, fr);
        
        for t = 1:fr
            for i=1:Par.NumPart
                if Distn.clusters{c}.particles{i}.tracks{j}.Present(t)
                    ass{jj}(i, t) = Distn.clusters{c}.particles{i}.tracks{j}.GetAssoc(t);
                    present(jj, t) = present(jj, t) + 1;
                end
            end
            
%             count(jj, t) = sum(ass{jj}(:, t)==correct(t, jj));
            
        end
        
    end
    
end

% figure, hold on
% for j = 1:Par.NumTgts
%     plot(count(j, :), 'color', [0, rand, rand])
% end

end