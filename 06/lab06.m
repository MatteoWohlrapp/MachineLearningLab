
D = readmatrix('data_clustering.csv');


epsilon = calculate_epsilon_plot_nn_graph(D, 4, 6);
experimental_results_plot(D, epsilon)
calculate_silhouette_score(D, epsilon)

% function to plot the three different clusters for varying Min_Pts
% according to the results section in the assignment 
% input: Data vectors D
function experimental_results_plot(D, epsilon)
    Min_Pts = [3,4,5];
    for i = 1:length(Min_Pts)
        Ps = DBSCAN.dbscan(D, epsilon, Min_Pts(i)); 
        plot_cluster(D, Ps, epsilon, Min_Pts(i));
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

% function to plot the three silhouette scores for varying Min_Pts
% according to the results section in the assignment 
% input: Data vectors D
function calculate_silhouette_score(D, epsilon)
    Min_Pts = [3,4,5];
   
    for i = 1:length(Min_Pts)
        Ps = DBSCAN.dbscan(D, epsilon, Min_Pts(i));
        s = silhouette(D, Ps); 
        disp(mean(s)); 
    end
end

%function to plot the clusters
function plot_cluster(D, Ps, eps, Min_Pts)
    figure_name = sprintf('DBSCAN_eps%d_minpts%d.pdf', eps, Min_Pts); 
    fig = figure('Name', figure_name); 
    gscatter(D(:,1), D(:,2), Ps)

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


  