DefineParameters;

t = 20;
L = 20;
obs = cell(L, 1);
for k = 1:L
    obs{k} = Observs(t-L+k).r(1, :)';
end
init_track = Track(0, 1, {Par.TargInitState-[Par.TargInitState(3:4)' 0 0]'}, 0);
[m, v] = KalmanFilter(obs,init_track.GetState(0), Par.Q);
PlotTrueState(TrueState)
for k=1:20, plot(m{k}(1), m{k}(2), 'bx'), end
for k=1:20, plot_gaussian_ellipsoid(m{k}(1:2), v{k}(1:2,1:2), 2), end