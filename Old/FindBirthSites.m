function BirthSites = FindBirthSites( t, Observs )
%FINDBIRTHSITES Locate points where a target could have been born in this
% frame.

global Par;

BirthSites = cell(0,1);

for i1 = 1:Observs(t).N
    coords = Pol2Cart(Observs(t).r(i1, 1), Observs(t).r(i1, 2));
    
    for i2 = 1:Observs(t-1).N
        prev_coords = Pol2Cart(Observs(t-1).r(i2, 1), Observs(t-1).r(i2, 2));
        dist1 = EuclidDist(coords, prev_coords);
        
        if (dist1 < (Par.Vlimit*Par.P))&&(Observs(t-1).r(i2, 2)>Par.BirthExclusionRadius)
            
            for i3 = 1:Observs(t-2).N
                prev_prev_coords = Pol2Cart(Observs(t-2).r(i3, 1), Observs(t-2).r(i3, 2));
                dist2 = EuclidDist(prev_coords, prev_prev_coords);
                
                if dist2 < (Par.Vlimit*Par.P)&&(Observs(t-2).r(i3, 2)>Par.BirthExclusionRadius)
                    
                    BirthSites = [BirthSites; {[i3, i2, i1]}];
                    
                end
                
            end
            
        end
        
    end
    
    
end

end


function dst = EuclidDist(x1, x2)
dst = sqrt( (x2(1)-x1(1))^2 + (x2(2)-x1(2))^2);
end