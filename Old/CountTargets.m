t = 20;
tt = 20;

Distn = Distns{t};

num_targs = zeros(Par.NumPart,1);

for ii = 1:Par.NumPart
    
    for j = 1:Distn.particles{ii}.N
        if Distn.particles{ii}.tracks{j}.Present(tt)
            num_targs(ii) = num_targs(ii) + 1;
        end
    end
end

mean(num_targs)