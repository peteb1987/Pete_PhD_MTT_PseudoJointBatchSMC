function norm = NormLogArr( arr )
%NORMLOGARR Normalise array of log probabilities

max_arr = max(arr);
max_arr = max_arr(1);
norm = arr - max_arr;
norm = exp(norm);
norm = norm/sum(norm);
norm = log(norm);

end

