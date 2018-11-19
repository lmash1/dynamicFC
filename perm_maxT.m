%Adapted from Jeanette Mumford's permutation (maxT distribution) code
%available at https://www.dropbox.com/sh/kslgoqxwqda3hir/AABRbp1P8sI_bW2WPv2ZkTOfa?dl=0


%Steps of the max-statistic based permutation test
%Step 1:  Calculate test statistics for all ROI-ROI t-tests without use of
%permutations
%
%Step 2:  Randomly switch TD-ASD assignments.  Note, if the null were true in this case, then it shouldn't
%matter which group each observation is assigned.
%
%Step 3: Calculate the all ROI-ROI t-scores and store the value of the
%maximum of all t-scores
%
%Step 4: Repeat 2 & 3 5000 times
%Step 5:  Compare the original test statistics from step 1 to the cutoff
%corresponding to the 95th percentile of the max tstat distribution

function [sig_tVals, sig_pVals, sig_mes,p_perm, tmax] = perm_maxT(roi_comps_rho, roi_comps,dx)


%=================================
%Step 1: Calculate ROI-ROI t-tests

group_ttest = zeros(1,length(roi_comps));
group_pval = zeros(1,length(group_ttest));
group_mes = zeros(1,length(group_ttest));
disp('Calculating t-tests without permutation..')

[h,p,CI,stats]=ttest2(roi_comps_rho(dx==1,:),roi_comps_rho(dx==2,:),'tail','both');
group_ttest = stats.tstat;
group_pval = p;

%Measures of effect size- hedge's G
mesStats= mes(roi_comps_rho(dx==1,:), roi_comps_rho(dx==2,:), 'hedgesg');
group_mes = mesStats.hedgesg;

%=================================
%Steps 2 & 3 &4

nperm=5000;
tmax=zeros(1,nperm);
for i=1:nperm
    fprintf('Permutation Number %d \n',i);
    %randomly flip group assignments- sample group membership WITHOUT
    %REPLACEMENT- so we still have 62 ASD, 57 TD
    groupAssign_orig = dx;
    random_flip=randsample(groupAssign_orig, length(dx), false);
    
    %calculate all ROI-ROI t stats for permuted data
    
    [h,p,CI,stats]=ttest2(roi_comps_rho(random_flip==1,:),roi_comps_rho(random_flip==2,:),'tail','both');
    group_ttest_loop = stats.tstat;
    
    
    %Cacluate the max t stat for this iteration
    tmax(i) = max(abs(group_ttest_loop));
end

% Step 5:  Calculate the percentiles for each statistic based on
% permutation-based null dist

for i=1:length(roi_comps)
    p_perm(i)=sum(tmax>=abs(group_ttest(i)))/nperm;
end

sig_tVals = group_ttest;
sig_tVals(p_perm>.05)= 0;
sig_pVals = p_perm;
sig_pVals(p_perm>.05)= 0;
sig_mes = group_mes;
sig_mes(p_perm>.05)=0;

end