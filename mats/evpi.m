clear

% load costs and vaccine distributions
load('cost_tensor.mat')
load('jointdists.mat')

% number of time points
ntime = size(fs,3);

% number of vaccine points
nvacc = size(fs,4);

% seperate by strategies and assume w = 0.5
S1Costs = reshape(fs(51,1,:,:),ntime,nvacc);
S2Costs = reshape(fs(51,2,:,:),ntime,nvacc);
S3Costs = reshape(fs(51,3,:,:),ntime,nvacc);
S4Costs = reshape(fs(51,4,:,:),ntime,nvacc);

evpi_pos = zeros(1,4);
evpi_neg = zeros(1,4);

% pick vaccine distribution
for p = 1:size(P,1)
    VaccPD = reshape(P(p,:,:),ntime,nvacc);

    % average of optimal values conditional on each model (uncertainty)
    % optimum of average of values conditional on each model (uncertainty)
    p_Cik = [sum(VaccPD.*S1Costs,'all') sum(VaccPD.*S2Costs,'all') sum(VaccPD.*S3Costs,'all') sum(VaccPD.*S4Costs,'all')];
    evpi_neg(p) = max(-p_Cik);
    %evpi_neg(p) = min(p_Cik);

    for t = 1:ntime
        for v = 1:nvacc
            opt_Cik = max([-S1Costs(t,v) -S2Costs(t,v) -S3Costs(t,v) -S4Costs(t,v)]);
            %opt_Cik = min([S1Costs(t,v) S2Costs(t,v) S3Costs(t,v) S4Costs(t,v)]);
            evpi_pos(p) = evpi_pos(p) + VaccPD(t,v)*opt_Cik;
            
        end
    end
end

EVPI = evpi_pos - evpi_neg
