% README: dendogram_cluster(..) creates all necessary plots
% if not all graphs need to be created, the different methods like
% plot_dendogram(..), plot_cluster(..), plot_silhouette(..) can be used
% wss and bss calculation also has seperate values
data = readmatrix('data_clustering.csv');

dendogram_cluster_silhouette(data);
plot_bss_wss(data);


% function creating all plots necessary for the exercise 
% input: data vectors
function dendogram_cluster_silhouette(data)
    methods = ["single", "complete", "average", "ward"];
    plot_original_data(data); 
    for i = 1:length(methods)
        link = plot_dendogram(methods(i), data); 
        for K = 2:4
            clust = plot_cluster(methods(i), link, data, K); 
            plot_silhouette(methods(i), clust, data, K);
       end
    end
end 

% function for plottig the dendogram
% input: data vectors, linkage method as string
% output: linkage
function link = plot_dendogram(method, data)
    
    link = linkage(data, method);
        
    fig = figure('Name', method + "_dendogram_");
    dendrogram(link)
    yline(median([link(end-1,3),link(end,3)]),'--b','K = 2')
    yline(median([link(end-2,3),link(end-1,3)]),'--b','K = 3')
    yline(median([link(end-3,3),link(end-2,3)]),'--b','K = 4')
    xlabel('Data vector');
    ylabel('Proximity');
    grid
    save_plot(fig, method + "_dendogram.pdf");
end    

% function to plot the clusters for given K
% input: linkage method as string, linkage, data vectors, K
function clust = plot_cluster(method, link, data, K)
    clust = cluster(link, 'Maxclust', K);
    
    figure_name = append(method, sprintf('_cluster_K%d', K));
    fig = figure('Name', figure_name);
    
    gscatter(data(:,1),data(:,2),clust);
    xlabel('x');
    ylabel('y');
    grid
    save_plot(fig, figure_name + ".pdf");
end 

% function to plot the clusters for given K and the silouhette
% input: linkage method as string, linkage, data vectors, K
function plot_silhouette(method, clust, data, K)  
    figure_name = append(method, sprintf('_silhouette_K%d', K));
    fig = figure('Name', figure_name);
    
    silhouette(data, clust);
    xlabel('Silhouette coefficient values');
    ylabel('Cluster label');
    grid
    save_plot(fig, figure_name + ".pdf");
end 

% function to plot the original data
% input: data vectors
function plot_original_data(data)
    fig = figure('Name', 'Original data');
    scatter(data(:,1),data(:,2));
    xlabel('x');
    ylabel('y');
    grid
    save_plot(fig, 'Original_data.pdf');
end 

% function to calculate and plot bss and wss 
% input: data vectors
function plot_bss_wss(data) 
   wss_bss_values = zeros(1, 2); 
   xs = zeros(1, 2);
   for k = 2:4
       link = linkage(data, 'ward'); 
       clust = cluster(link, 'MaxClust', k); 
       wss_bss_values(k-1, 1) = wss(data,clust,k); 
       wss_bss_values(k-1, 2) = bss(data,clust,k); 
       xs(k-1, 1) = k; 
       xs(k-1, 2) = k; 
   end
   fig = figure('Name', 'WSS-BSS');
    bar(xs, wss_bss_values, 'grouped'); 
    hold on
    h = zeros(1,1); 
    h(1) = plot(NaN,NaN, '.b');
    h(2) = plot(NaN,NaN, '.r');
    legend(h, 'WSS', 'BSS'); 
    xlabel('K');
    ylabel('Value');
    grid
    save_plot(fig, 'WSS_BSS.pdf');

end

% function to save the plot
function save_plot(fig, name)
    set(fig, 'PaperPosition', [0 0 20 20]);
    set(fig, 'PaperSize', [20 20]);
    saveas(fig, name);
end 

% function to calculate wss
% input: data vectors, clusters and K
function sum = wss(data, clust, K)
    sum = 0; 
    for k = 1:K 
       [m_x,m_y] = representative_point(data, clust, k);
       for i = 1:length(clust)
           if clust(i) == k 
            sum = sum + (data(i,1)- m_x).^2 + (data(i,2) - m_y).^2;
           end
       end
    end
end

% function to calculate bss
% input: data vectors, clusters and K
function sum = bss(data, clust, K)
   [m_x,m_y]  = mean_of_cluster_center(data, clust, K); 
    sum = 0; 
    for i = 1:K 
        size = cluster_size(clust, i); 
        [m_ix,m_iy] = representative_point(data, clust, i); 
        sum = sum + size * ((m_x- m_ix).^2 + (m_y - m_iy).^2);
    end 
end 

% function to calculate the mean of the center of all clusters 
% input: data vectors, clusters and K
function [x,y] = mean_of_cluster_center(data, clust, K) 
    x = 0; 
    y = 0; 
    for i = 1:K 
        [m_x,m_y] = representative_point(data, clust, i); 
        x = x + m_x; 
        y = y + m_y; 
    end 
    x = x / K; 
    y = y / K; 
end 

% function to calculate the size of cluster k 
% input: clusters and k
function size = cluster_size(clust, k) 
    size = 0; 
    for i=1:length(clust) 
        if clust(i) == k 
            size = size + 1; 
        end
    end
end

% function to find the representative point for cluster k 
% input: data vectors, clusters and K
function [x,y] = representative_point(data, clust, k)
    x = 0;
    y = 0;
    amount = 0; 
    for i = 1:length(clust) 
        if clust(i) == k 
            x = x + data(i,1);
            y = y + data(i, 2); 
            amount = amount + 1; 
        end
    end
    x = x / amount;
    y = y / amount;
end 
