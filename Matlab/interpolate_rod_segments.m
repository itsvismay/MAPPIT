function v = interpolate_rod_segments(x, n_seg)
    if(n_seg < size(x,1) -1)
        sprintf("must have more segments")
    end
    n_seg = n_seg + 1;
    dx = x(2:end, :) - x(1:end-1,:);
    dxdx = dx.*dx;
    l_dx = sum(dxdx,2);
    
    [B, I] = sort(l_dx, 'descend');
    
    num_points_to_insert_between_I = zeros(size(I));
    
    if floor(n_seg/size(I,1))>0
        num_points_to_insert_between_I = num_points_to_insert_between_I + double(floor(n_seg/size(I,1)));
    end
    
    size(I)
    if mod(n_seg,size(I, 1))>0
        %if this fails:
        % - possibly the terrain mesh connectivity isn't enough? Subdivide the
        % terrain more.
        % - possibly increase the number of rod segments
        num_points_to_insert_between_I(1:mod(n_seg,size(I,1))) = num_points_to_insert_between_I(1:mod(n_seg,size(I,1))) + 1;
    end
    
    num_points_to_insert_between = zeros(size(I));
    num_points_to_insert_between(I) = num_points_to_insert_between_I;
    
    e = [(1:size(x,1) -1)' (2:size(x, 1))'];
    [v, Edges, J] = resample_rod(x, e, num_points_to_insert_between);
    
end


    