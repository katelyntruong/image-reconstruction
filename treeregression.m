load MysteryImage.mat

m=1456;
n=2592;

x_data = cols/n;
y_data = (m/n)*(1 - rows/m);
n_data = 15000;

%structure of R(1)
R(1).s = [];
R(1).j = [];
R(1).left = [];
R(1).right = [];
R(1).I = (1:n_data);
R(1).x = [0;1]; % x-coordinates of the corners
R(1).y = [0;m/n]; % y-coordinates of the corners
R(1).val = 1/n_data*sum(vals); % average over RGB values of data in R(j)

leaves = [1];
size = 1;
splits = 5000;

for split=1:splits
    k = leaves(1);
    leaves = leaves(2:end);
    [j_opt,s_opt] = OptimalSplitRegression(x_data,y_data,vals, R(k).I);   
    R(k).j = j_opt;
    R(k).s = s_opt;
    if j_opt == 0
        continue
    end
    R(k).left = size + 1;
    R(k).right = size + 2;
    leaves = [leaves size+1 size+2];
    Ileft = [];
    Iright = [];
    
    if j_opt == 1
        for i=R(k).I
            if x_data(i) < s_opt
                Ileft = [Ileft i];
            else
                Iright = [Iright i];
            end
        end
        R(size+1).x = [R(k).x(1);s_opt];
        R(size+1).y = R(k).y;
        R(size+2).x = [s_opt;R(k).x(2)];
        R(size+2).y = R(k).y;
    else
        for i=R(k).I
            if y_data(i) < s_opt
                Ileft = [Ileft i];
            else
                Iright = [Iright i];
            end
        end
        R(size+1).x = R(k).x;
        R(size+1).y = [R(k).y(1);s_opt];
        R(size+2).x = R(k).x;
        R(size+2).y = [s_opt;R(k).y(2)];
    end
    
    
    R(size+1).I = Ileft;
    R(size+1).val = mean(vals(Ileft,:),1);

    R(size+2).I = Iright;
    R(size+2).val = mean(vals(Iright,:),1);
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