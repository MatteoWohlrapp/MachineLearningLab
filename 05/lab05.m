
data = readmatrix('data_clustering.csv');

dendogram_cluster(data);


function dendogram_cluster(data)
    methods = ["single", "complete", "average", "ward"];
    plot_original_data(data); 
    for i = 1:length(methods)
        link = plot_dendogram(methods(i), data); 
        for K = 2:4
            plot_cluster_silhouette(methods(i), link, data, K); 
        end
    end
end 

%function for printing, and saving the clustering data as defined.
function link = plot_dendogram(method, data)
    
    link = linkage(data, method);
        
    fig = figure('Name', method + "_dendogram_");
    
    dendrogram(link)
    yline(median([link(end-1,3),link(end,3)]),'--b','K = 2')
    yline(median([link(end-2,3),link(end-1,3)]),'--b','K = 3')
    yline(median([link(end-3,3),link(end-2,3)]),'--b','K = 4')
    grid
    saveas(fig,method + "_dendogram.pdf");
end    

function plot_cluster_silhouette(method, link, data, K)
    clust = cluster(link, 'Maxclust', K);
    
    figure_name = append(method, sprintf('_cluster_K%d', K));
    fig = figure('Name', figure_name);
    
    gscatter(data(:,1),data(:,2),clust);
    grid
    saveas(fig, figure_name + ".pdf");
    
    figure_name = append(method, sprintf('_silhouette_K%d', K));
    fig = figure('Name', figure_name);
    
    silhouette(data, clust);
    grid
    saveas(fig, figure_name + ".pdf");
end 

function plot_original_data(data)
        fig = figure('Name', 'Original data');
        scatter(data(:,1),data(:,2)); 
        grid
        saveas(fig, 'Original_data.pdf');
end 
