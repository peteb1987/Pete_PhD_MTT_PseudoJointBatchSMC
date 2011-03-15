function [ AssocVector ] = AuctionAssoc( State, Obs )
%AUCTIONASSOC Use the auction algorithm to find the maximum likelihood
%assignment between targets and observations

% State is a cell array of target states
% Obs is an array of observations (row-wise)

global Par;

No = size(Obs, 1);
Ns = size(State, 1);

if (Ns == 0)||(No == 0)
    AssocVector = zeros(Ns, 1);
    return
end

% First create an association array with No extra elements corresponding
% to a clutter assignment for each observation
Assoc = zeros(No, No+Ns);

% Create an array of prices
Prices = zeros(1, No+Ns);

% Create an array of payoffs for an assignment, the log of the likelihood
% of each observation given an association with a target
temp = -inf*( xor(ones(No), eye(No)) );
temp(isnan(temp)) = 0;
Payoffs = [temp zeros(No, Ns)];
clear temp;

for j = 1:Ns
    
    x = State{j};
    [bng, rng] = Cart2Pol(x(1:2));
    range_squ = x(1)^2 + x(2)^2;
    range = sqrt(range_squ);
    jac = [-x(2)/range_squ, x(1)/range_squ, 0, 0; x(1)/range, x(2)/range, 0, 0];
    mean_obs = [bng; rng];
    var_obs = Par.R + jac * Par.Q * jac';
    
    for i = 1:No
%         Payoffs(i, No+j) = log(  mvnpdfFastSymm(Obs(i, :), State{j}(1:2), Par.AuctionVar) / Par.UnifPosDens  );
        if Dist(Pol2Cart(Obs(i, 1), Obs(i, 2)), Pol2Cart(mean_obs(1), mean_obs(2)))<Par.Vmax
            Payoffs(i, No+j) = log(  100000 * Par.PDetect*mvnpdf(Obs(i, :), mean_obs', var_obs) / (Par.ClutDens*(1-Par.PDetect)) );
        else
            Payoffs(i, No+j) = -inf;
        end
        
%         % If the payoff is less than that of a clutter assignment, disallow
%         % assignment by setting payoff to -inf.
%         if Payoffs(i, No+j) < 0
%             Payoffs(i, No+j) = -inf;
%         end
    end
end        

% Auction iterations
done = false;
while ~done
    
    % Iterate over observations ("persons")
    for i = 1:No
        
        if sum(Assoc(i, :))==0
            
            % Find the target or clutter ("object") that gives the biggest reward
            Rewards = Payoffs(i, :) - Prices;
            k = find(Rewards==max(Rewards));
            
            if length(k)>1
                disp('More than one best assignment!!!');
                k = k(1);
            end
            
            % Calculate the bidding increment
            BestReward = Rewards(k);
            Rewards(k) = -inf;
            
            SecondReward = Rewards(Rewards==max(Rewards));
            SecondReward = SecondReward(1);
            Bid = BestReward - SecondReward;
            
            if isnan(Bid)
                error('Bid is NaN');
            end
            
            % Increase the price and assign this target
            Prices(k) = Prices(k) + Bid;
            Assoc(:, k) = 0;
            Assoc(i, k) = 1;
            
        end
  
    end
    
    if all(sum(Assoc, 2))
        done = true;
    end
    
end

% Remove the associations corresponding to clutter
Assoc(:, 1:No) = [];

AssocVector = zeros(Ns, 1);
for j = 1:Ns
    if ~isempty(find(Assoc(:, j), 1))
        AssocVector(j) = find(Assoc(:, j));
    end
end

end %AuctionAssoc

function d = Dist(x1, x2)
d = sqrt((x1(1)-x2(1))^2+(x1(2)-x2(2))^2);
end