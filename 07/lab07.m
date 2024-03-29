% README
% How to use: the function PCA is the algorithm retreived from the pseude
% code
% To visualize the data after reduction, call visualize_PCA(..), to see the
% eigen-values at the eigen-values indice, call visualize_eigenvalues(..),
% report_dimensionality(..) can be used to see which dimensionality is
% needed for given threshold values in variance

z = load('COIL20.mat'); 
% desired dimension for reduction
d = 30; 

visualize_eigenvalues(z.X);
report_dimensionality(z.X, d);
visualize_PCA(z.X, d, z.Y);

% function to implement the PCA algorithm according to the pseudeo code 
% input: data vectors X, dimensionality d
% output: d eigen-vectors U_d, eigen-values D, compressed data Z_d
function [U_d, D, Z_d] = PCA(X, d)
    % centralize
    X_z = centralize(X); 
    % create co-variance matrix
    C = co_variance_matrix(X_z); 
    %compute eigen-vectors and eigen-values
    [U, D] = eig(C);
    % pick only the necessary ones 
    U_d = pick_eigen_vectors(U, d); 
    % reduce the dimenstionality
    Z_d = reduce_dimensionality(U_d, X_z);
end 

% function to visualize the eigenvalues according to the results
% section in the assignments
% input: data vectors X
function visualize_eigenvalues(X)
    % centralize
    X_z = centralize(X); 
    % create co-variance matrix
    C = co_variance_matrix(X_z); 
    %compute eigen-vectors and eigen-values
    [~, D] = eig(C);
    eigenvalues = zeros(1, 1); 
    max_value = length(D(:,1));
    for i = 1:max_value
        % eigenvalues sorted ascending, we need highest ones first
        eigenvalues(i) = D(max_value - i + 1, max_value - i + 1); 
    end
    plot_eigenvalues(eigenvalues);
end 

% function to report the dimensionality for 3 given variances
% input: data vectors C, dimensions d
% output: to std-out: the variance value with the needed dimensionality 
function report_dimensionality(X, d)
    [U_d, D, ~] = PCA(X, d); 
    threshold_variances = [0.9,0.95,0.98];
    index = 1; 
    % retrieve all PCA_components
    PCA_components = zeros(1, 1); 
    max_value = length(U_d(:,1));
    for i = 1:max_value
        % PCA components sorted ascending, we need highest ones first
        PCA_components(i) = D(max_value - i + 1, max_value - i + 1); 
    end
    for i = 1:length(PCA_components)
        if index > length(threshold_variances)
            break
        end
        % calculate variance
        variance = sum(PCA_components(1:i)) / sum(PCA_components);
        % check if dimensionality is enough for the variance
        if variance > threshold_variances(index)
            fprintf('Variance: %d, dimensionality: %d\n', threshold_variances(index), i); 
            index = index + 1; 
        end
    end

end 

% function to visualize the data by using t-SNE
% input: data vectors X, dimensionality d, cluster labels 
function visualize_PCA(X, d, labels)
    [~, ~, Z_d] = PCA(X, d);
    %visualize using tsne 
    Y = tsne(transpose(Z_d));
    plot_t_SNE(Y, labels);
end

% function to centralize the data set 
% input: data vectors X 
% output: centralized data vectors X_z
function X_z = centralize(X)
    % calculate mean
    mean = find_mean(X); 
    X_z = X;
    for i = 1:length(X_z(:,1))
        % deduct mean from every vector
        X_z(i, :) = X_z(i, :) - mean; 
    end
end     

% function to find the mean of the data set 
% input: data vectors X
% output: mean of data
function mean = find_mean(X) 
    mean = 0; 

    for i = 1:length(X(:,1))
        mean = mean + X(i, :); 
    end 
    mean = mean / length(X(:,1)); 
end     

% function to calculate the co-variance matrix to data X_z 
% input: zentraliced data vectors X_z 
% output: co-variance matrix
function C = co_variance_matrix(X_z)
    C = transpose(X_z(1,:)) * X_z(1,:);
    
    for i = 2:length(X_z(:,1))
        C = C + transpose(X_z(i,:)) * X_z(i,:);
    end
    C = C / length(X_z(:,1));
end 

% function to pick the first d eigenvalues
% input: eigen-vectors U, dimensionality d
% output: d highest eigen-vectors
function U_d = pick_eigen_vectors(U, d)
    max_value = length(U(:,1));
    U_d = zeros(max_value, d); 
    % eigenvectors in U are sorted ascending, thats why max_value -1 +1 is
    % used
    for i = 1:d
        U_d(:,i) = U(:,max_value-i+1); 
    end
end

% function to reduce the dimensionality of X_z
% input: d eigen-vectors U_d, zentralized data vectors X_z
% output: dimensionality reduced vectors Z_d
function Z_d = reduce_dimensionality(U_d, X_z)
    % transpose matrix to get correct result
    Z_d = transpose(U_d) * transpose(X_z);
end 

% function to plot the compressed data using t-SNE
% input: zentralized, dimensionaly reduced data vectors Z_d, cluster labels
function plot_t_SNE(Y, labels)
    figure_name = 't-SNE_representation.pdf'; 
    fig = figure('Name', figure_name); 
    gscatter(Y(:,1), Y(:,2), labels)
    xlabel('x');
    ylabel('y');
    grid
    save_plot(fig, figure_name);
end

% function to plot the eigenvalues according to the results
% section in the assignments
% input: data vectors X
function plot_eigenvalues(eigenvalues)
    figure_name = 'eigenvalues.pdf'; 
    fig = figure('Name', figure_name); 
    plot(1:length(eigenvalues),eigenvalues, '-r.', 'LineWidth',1) 
    hold on 
    xlabel('eigen-values indices');
    ylabel('eigen-value');
    h = zeros(1,1); 
    h(1) = plot(NaN,NaN, '.r');
    legend(h, 'eigen-values sorted by value'); 
    grid
    save_plot(fig, figure_name);
end 

% function to save the plot
function save_plot(fig, name)
    set(fig, 'PaperPosition', [0 0 20 20]);
    set(fig, 'PaperSize', [20 20]);
    saveas(fig, name);
end 
