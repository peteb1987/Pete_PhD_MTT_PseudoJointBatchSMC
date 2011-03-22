function [ass, count] = AnalyseAss( correct, Distns, fr)
%ANALYSEASS Compare associations with correct values

global Par;

ass = cell(Par.NumTgts, 1);
count = zeros(Par.NumTgts, fr);

jj = 1;

for c = 1:Distns{fr}.N
    
    for j = 1:Distns{fr}.clusters{c}.N
        
        ass{jj} = zeros(Par.NumPart, fr);
        
        for t = 1:fr
            for i=1:Par.NumPart
                ass{jj}(i, t) = Distns{fr}.clusters{c}.particles{i}.tracks{j}.GetAssoc(t);
            end
            
            count(jj, t) = sum(ass{jj}(:, t)==correct(t, jj));
            
        end
        
        jj = jj + 1;
        
    end
    
end

figure, hold on
for j = 1:Par.NumTgts
    plot(count(j, :), 'color', [0, rand, rand])
end

end