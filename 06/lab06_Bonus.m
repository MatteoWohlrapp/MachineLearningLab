load thyroid.mat

%Vector of outliers. 0 = regular, 1 = outlier.
outliers = y;

%Matrix of points, NX6
%Representation done with the 2 first values.
points = X;

%for ease of use, the values in y will be inverted.
%Thus making all representations "equal" 
for i = 1:length(y)
    if y(i) == 0
        outliers(i) = 1;
    else
        outliers(i) = 0;
    end
end

%We must calculate 3 things for each value we try:
%Precision, Recall and F-Score

epsilon = calculate_epsilon_plot_nn_graph(points, 4, 6);
%"Base" Clusters 
plot_cluster(points(:,1:2),outliers,epsilon,0)
experimental_results_plot(points, epsilon,outliers)

% function to plot the three different clusters for varying Min_Pts
% according to the results section in the assignment 
% input: Data vectors D
% Slight change on what we are printing, since we print 6 values per point
function experimental_results_plot(D, epsilon,true_classes)
    Min_Pts = [3,4,5];
    for i = 1:length(Min_Pts)
        Ps = DBSCAN.dbscan(D, epsilon, Min_Pts(i)); 
        %In here we have to add an additional thing, since all we want are
        %the outliers.
        %We must fuse all other "values" to become 1 in Ps.
        x = Only_outliers(Ps);
        %We then compare these values with the ones in the original for the
        %metrics:
        fprintf('Precision for K = %d\n', i+2)
        pr = precision(x,true_classes);
        disp(pr)
        
        fprintf('Recall for K = %d\n',i+2)
        r = recall(x,true_classes);
        disp(r)

        fprintf('F point for K = %d\n',i+2)
        f_p = f_point(pr,r);
        disp(f_p)

        plot_cluster(D(:,1:2), x, epsilon, Min_Pts(i));
    end
end

%Value = Outliers_correct/outliers_correct + outliers_false
function result = precision(test_classes,true_classes)

    correct_outliers = 0;
    false_outliers = 0;

    for i = 1:length(test_classes)
        
        if test_classes(i) == 0
            
            if test_classes(i) == true_classes(i)
                correct_outliers = correct_outliers + 1;
            else
                false_outliers = false_outliers + 1;
            end

        end

    end

    result = correct_outliers/(false_outliers + correct_outliers);
end


%Recall = TruePositives / (True Positives + False Negatives)
function result = recall(test_classes,true_classes)

    correct_positives = 0;
    false_negatives = 0;

    for i = 1:length(test_classes)
        
        if test_classes(i) == 0
            if test_classes(i) == true_classes(i)
                correct_positives = correct_positives +1;
            end
        end

        if test_classes(i) == 1
            if test_classes(i) ~= true_classes(i)
                false_negatives = false_negatives + 1;
            end
        end
    end

    result = correct_positives / (correct_positives + false_negatives);

end

%F-point = (2*precision * recall)/(precision + recall)
function result = f_point(precision,recall)
    result = (2*precision * recall)/(precision + recall);
end

function x = Only_outliers(Ps)
    x = zeros(length(Ps),1);

    for i = 1:length(Ps)
        if (Ps(i) ~= -1)
            x(i) = 1;
        end
    end    

end

% function to calculate the average epsilon and plot the k-NN graph for k
% from k_from to k_to according to the results section in the assignment
% input: Data vectors D, k-NN from, k-NN to
% output: average epsilon
function epsilon = calculate_epsilon_plot_nn_graph(D, k_from, k_to)
    epsilon = clusterDBSCAN.estimateEpsilon(D,k_from,k_to);
    clusterDBSCAN.estimateEpsilon(D,k_from,k_to)
end 


%function to plot the clusters
% 0 = normal | 1 = Outliers
function plot_cluster(D, Ps, eps, Min_Pts)
    figure_name = sprintf('DBSCAN_eps%d_minpts%d_bonus.pdf', eps, Min_Pts); 
    fig = figure('Name', figure_name); 
    gscatter(D(:,1), D(:,2), Ps,'rb','ox')
    xlabel('x');
    ylabel('y');
    grid
    save_plot(fig, figure_name);
end 

% function to save the plot
function save_plot(fig, name)
    set(fig, 'PaperPosition', [0 0 20 20]);
    set(fig, 'PaperSize', [20 20]);
    saveas(fig, name);
end 
