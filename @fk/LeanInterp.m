function Yi = LeanInterp(X, Y, Xi, splint)
%   Interpolate linearly between observations
%   X Y and Xi are 3d matrices, where i and j (i,j,k)
%   correspond to a Y(i,j)
       
    % get the slopes for each interval
    m = diff(Y,1,3)./diff(X,1,3);
    
    % we have 2 time scales: one of the resampled data and another with the
    % uniform time scale in the center array reference station time
    % [==================================]   % uniformly sampled
    % 200                                300
    %
    % [==================================]   % resampled data (no doppler)
    % 251.5                              441.5
    % 
    % To take scale 2 to agree with scale 1, we need to shift the data by
    % the integer number (200 - 251.5) mod sample rate.
    
    dif_sample = (X - Xi - mod(X - Xi, splint))/splint;
    
    % the absolute location on where to sample is in sample_at
    sample_at(1,1,:) = 1:size(Xi,3);
    sample_at = repmat(sample_at,[size(Xi,1) size(Xi,2) 1]) + dif_sample;
    % we can't sample outside the reference station time
    sample_at(sample_at <= 0) = 1;
    sample_at(sample_at > size(Xi,3)-1) = 1;
    
    sample_at = sample_at(:);
    row = (1:size(Xi,2))';
    row = repmat(row,[1 size(Xi,2) size(Xi,3)]);
    row = row(:);
    col = (1:size(Xi,2));
    col = repmat(col,[size(Xi,1) 1 size(Xi,3)]);
    col = col(:);
    
    % now, take it to index notation
    sample_at = sub2ind(size(Xi),row,col,sample_at);
    % the delta t var contains the difference between the uniform sampled
    % data point and the resampled data. This value will be multiplied by
    % the slope m to interpolate
    delta_t = mod(X - Xi, splint);
    
    % make the first sample 0. This is so that for the points that are
    % sample_at > size(Xi,3)-1 get a value of zero and don't distort the
    % borders
    Y(:,:,1) = zeros(size(Xi,1),size(Xi,2));
    
    %Yi = Y(sample_at); 
    Yi = m(sample_at).*delta_t(:) + Y(sample_at);
    Yi = reshape(Yi,size(Xi));
    %tt(:,:) = Yi(1,1,:);
    %plot(tt);

end