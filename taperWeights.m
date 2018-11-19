%ROIs_concat_mat is a cell length (# subjects) x (# windows)...
%...each cell in ROIs_mat is (#TRs x # ROIs)

% for each subject ?   
for s=1:length(subIDs)
    
    % taper each window
    for winNum = 1:NUM_WINDOWS_USED
        lengthWin = 1:size(ROIs_concat_mat{s,winNum},1);
        % gaussian curve with std 3, center at window median
        gaussWeights = gaussmf(lengthWin, [3, median(lengthWin)]); 
        % make sum of gaussWeights == 1
        new_weights = 1/sum(gaussWeights); 
        gaussWeights = gaussWeights * new_weights;
        
        %for each ROI timecourse, apply taper weights
        for i = 1:size(ROIs_concat_mat{s,winNum},2) 
            ROIs_concat_mat{s,winNum}(:,i) = ROIs_concat_mat{s,winNum}(:,i) .* gaussWeights'; %#ok<SAGROW>
        end
    end
end

