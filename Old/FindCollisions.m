for t = 1:Par.T
    
    used = cell(Par.NumTgts, 1);
    
    for j = 1:Par.NumTgts
        
        this_targ_ass = unique(ass{j}(:, t));
        this_targ_ass(this_targ_ass==0)=[];
        
        for i = 1:length(this_targ_ass)
            
            for jj = 1:j
                if any(this_targ_ass(i)==used{jj})
                    c1 = sum(ass{j}(:,t)==this_targ_ass(i));
                    c2 = sum(ass{jj}(:,t)==this_targ_ass(i));
                    if (c1>10)&&(c2>10)
                        disp(['Collision in frame ' num2str(t) ', target ' num2str(j) ' and target ' num2str(jj)]);
                        disp(['*** ' num2str(c1) ' particles and ' num2str(c2) ' particles respectively']);
                    end
                end
            end
            
        end
        
        used{j} = this_targ_ass;
        
    end
    
end