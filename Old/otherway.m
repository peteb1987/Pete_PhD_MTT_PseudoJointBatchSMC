test = zeros(N,1);

N = Observs(t).N;

for i = 1:N
    
    obs_cart  = Pol2Cart(Observs(t).r(i, 1), Observs(t).r(i, 2));
    mean_cart = Pol2Cart(mean_obs(1), mean_obs(2));
    if Dist(obs_cart, mean_cart)<Par.Vmax
        test(i) = (Par.PDetect/Observs(t).N) * mvnpdf(Observs(t).r(i, :), mean_obs', 5*var_obs);
    else
        test(i) = 0;
    end

end

test(N+1) = Par.ClutDens * (1-Par.PDetect);

test = test/sum(test);