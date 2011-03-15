function [ass, count] = AnalyseAss( correct, Distns)
%ANALYSEASS Compare associations with correct values

global Par;

ass = cell(Par.NumTgts, 1);
count = zeros(Par.NumTgts, Par.T);

for j = 1:Par.NumTgts
    
    ass{j} = zeros(Par.NumPart, Par.T);
    
    for t = 1:Par.T
        for i=1:Par.NumPart
            ass{j}(i, t) = Distns{t}.clusters{j}.particles{i}.tracks{1}.GetAssoc(t);
        end
        
        count(j, t) = sum(ass{j}(:, t)==correct(t, j));
        
    end
    
end

figure, hold on
for j = 1:Par.NumTgts
    plot(count(j, :), 'color', [0, rand, rand])
end

end