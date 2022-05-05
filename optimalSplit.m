x1 = 0.25; y1 = 0.25; v1 = [0.1, 0.5, 0.2];
x2 = 0.8; y2 = 0.1; v2 = [0.4, 0.1, 0.3];
x3 = 0.6; y3 = 0.8; v3 = [0.9, 0.2, 0.0];
x_data = [x1; x2; x3];
y_data = [y1; y2; y3];
vals = [v1; v2; v3];
R(1).s = [];
R(1).j = [];
R(1).left = [];
R(1).right = [];
R(1).I = 1:3;
R(1).x = [0;1]; % x-coordinates of the corners
R(1).y = [0;1]; % y-coordinates of the corners
R(1).val = 1/3*sum(vals); % average over RGB values of data in R(j)

leaves = [1];
size = 1;
splits = 2;
for split=1:splits
    k = leaves(1);
    leaves = leaves(2:end);
    [j_opt,s_opt,index_opt,s_vector] = OptimalSplitRegression(x_data,y_data,vals, R(k).I);
    valscurr = vals(R(k).I, :);
    xcurr = x_data(R(k).I);
    ycurr = y_data(R(k).I);
    %sort data in x-direction
    [sorted_x_data, SortedOrderX] = sort(xcurr);

    %sort data in y-direction
    [sorted_y_data, SortedOrderY] = sort(ycurr);

    begin_I = R(k).I(1);
    split_I = R(k).I(index_opt);
    last_I = R(k).I(end);    
    
    R(k).j = j_opt;
    R(k).s = s_opt;
    R(k).left = size + 1;
    R(k).right = size + 2;
    leaves = [leaves size+1 size+2];

    if j_opt == 1
        sorted_y_data = sorted_y_data(SortedOrderX);
        sorted_vals = valscurr(SortedOrderX,:);
        R(size+1).x = [R(k).x(1);s_opt];
        R(size+1).y = R(k).y;
        R(size+2).x = [s_opt;R(k).x(2)];
        R(size+2).y = R(k).y;
    else
        sorted_x_data = sorted_x_data(SortedOrderY);
        sorted_vals = vals(SortedOrderY,:);
        R(size+1).x = R(k).x;
        R(size+1).y = [R(k).y(1);s_opt];
        R(size+2).x = R(k).x;
        R(size+2).y = [s_opt;R(k).y(2)];
    end
    
    
    R(size+1).I = begin_I:split_I;
    R(size+1).val = mean(sorted_vals(1:index_opt,:),1);

    R(size+2).I = (split_I+1):last_I;
    R(size+2).val = mean(sorted_vals((index_opt+1):end,:),1);
    size = size + 2;
end
hold on
for idx_leave = leaves
    fill([R(idx_leave).x(1),R(idx_leave).x(2),R(idx_leave).x(2),R(idx_leave).x(1)],[R(idx_leave).y(1),R(idx_leave).y(1),R(idx_leave).y(2),R(idx_leave).y(2)],R(idx_leave).val);
end

function [j_opt,s_opt,index_opt, s_vector] = OptimalSplitRegression(xdata,ydata,valsdata, I_curr)
valscurr = valsdata(I_curr, :);
xcurr = xdata(I_curr);
ycurr = ydata(I_curr);
%sort data in x-direction
[sorted_x_data, SortedOrderX] = sort(xcurr);

%sort data in y-direction
[sorted_y_data, SortedOrderY] = sort(ycurr);

sorted_unique_x = sort(unique(xcurr));

sorted_unique_y = sort(unique(ycurr));

final_impurity = -1;
s_opt = 0;
j_opt = 0;
index_opt = 0;
s_vector = [];
%spliting in x-direction
for i=1:(size(sorted_unique_x)-1)
    s = mean([sorted_unique_x(i), sorted_unique_x(i+1)]);
    vals_sorted = valscurr(SortedOrderX,:);
    %vals = vals_sorted(I_curr, :);
    
    left_vals = vals_sorted(1:i,:);
    right_vals = vals_sorted((i+1):end,:);
    c1 = mean(left_vals,1);
    c2 = mean(right_vals,1); 
    
    %caculate purity
    half_1_impurity = sum((left_vals-ones(size(left_vals,1),1)*c1).^2, 'all');
    half_2_impurity = sum((right_vals-ones(size(right_vals,1),1)*c2).^2, 'all');
    impurity = half_1_impurity + half_2_impurity;
    s_vector = [s_vector s];
    if final_impurity == -1 | final_impurity > impurity
        final_impurity = impurity;
        j_opt = 1;
        s_opt = s;
        index_opt = i;
    end
end


%spliting in y-direction
for i=(1:size(sorted_unique_y)-1)
    s = mean([sorted_unique_y(i), sorted_unique_y(i+1)]);
    vals_sorted = valsdata(SortedOrderY,:);
    %vals = vals_sorted(I_curr, :);
    
    left_vals = vals_sorted(1:i,:);
    right_vals = vals_sorted((i+1):end,:);
    c1 = mean(left_vals,1);
    c2 = mean(right_vals,1); 
    
    %caculate purity
    half_1_impurity = sum((left_vals-ones(size(left_vals,1),1)*c1).^2, 'all');
    half_2_impurity = sum((right_vals-ones(size(right_vals,1),1)*c2).^2, 'all');
    impurity = half_1_impurity + half_2_impurity;
    s_vector = [s_vector s];
    if final_impurity == -1 | final_impurity > impurity
        final_impurity = impurity;
        j_opt = 2;
        s_opt = s;
        index_opt = i;
    end
end
end