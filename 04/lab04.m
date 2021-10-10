edit lab04

% w6_1x is array in struct
% 1000x2 double array
vqdata = load('simplevqdata.mat'); 

K = 1; 

% try different values
n = 0.1; 
t_max = 800; 

%vector_quantization(vqdata.w6_1x, K, n, t_max);
%vector_quantization_K2_K4(vqdata.w6_1x, t_max);
%ellbow(vqdata.w6_1x, n, t_max);

%test_learning_rate(vqdata.w6_1x, t_max)
test_t_max(vqdata.w6_1x, n)


% function calculating vector quantization
function last_H_VQ = vector_quantization(points, K, n, t_max)
    prototypes = select_prototypes(points, K);
    H_VQ = do_training(points, prototypes, n, t_max); 
    %plot_learning_curve(H_VQ, n, K);
    last_H_VQ = H_VQ(length(H_VQ));
end 

function test_learning_rate(points, t_max)
    n = 0.0;
    learning_rate_H_VQ = zeros(1,2);
    index = 1; 
    while n < 1
        H_VQ = vector_quantization(points, 2, n, t_max);
        learning_rate_H_VQ(index,:) = [n, H_VQ];
        index = index + 1;
        n = n + 0.05;
    end
    
    plot_best_learning_rate(learning_rate_H_VQ);
end 

function test_t_max(points, n)
    t_max = 1;
    t_max_H_VQ = zeros(1,2);
    index = 1; 
    while t_max < 2000
        H_VQ = vector_quantization(points, 2, n, t_max);
        t_max_H_VQ(index,:) = [t_max, H_VQ];
        index = index + 1;
        t_max = t_max + 200;
    end
    
    plot_best_t_max(t_max_H_VQ);
end 

% function to plot and calculate for 3 learning rates for K = 2 and K = 4
function vector_quantization_K2_K4(points, t_max)
    ns = [0.01,0.1,0.9];
    
    for i = 1:length(ns)
        prototypes_K2 = select_prototypes(points, 2); 
        H_VQ_K2 = do_training(points, prototypes_K2, ns(i), t_max); 
        prototypes_K4 = select_prototypes(points, 4); 
        H_VQ_K4 = do_training(points, prototypes_K4, ns(i), t_max); 
        plot_learning_curve_K2_K4(H_VQ_K2, H_VQ_K4, ns(i));
    end 
end 

function ellbow(points, n, t_max)
    final_H_VQs = zeros(1,1);
    for K = 1:5
        prototypes = select_prototypes(points, K);
        H_VQ = do_training(points, prototypes, n, t_max); 
        final_H_VQs(K) = H_VQ(t_max);
    end
    plot_ellbow(final_H_VQs, n);
end 

% function to select a random point for the prototypes
function prototypes = select_prototypes(points, K)
    prototypes = zeros(K,2);
    for i = 1:K
        index = floor((length(points)-1).*rand(1,1) + 1);
        prototypes(i,:) = points(index,:);
    end   
end 

% function to train the prototypes, returns the error by means of the
% distance
function H_VQs = do_training(points, prototypes, n, t_max)
    
    H_VQs = zeros(t_max,1);
    
    
    for t = 1:t_max
        %randperm randomize the rows while keeping the columns intact. 
        %Source: https://www.mathworks.com/matlabcentral/answers/30345-swap-matrix-row-randomly
        permutated_points = points(randperm(size(points, 1)), :);
        
        for i = 1:length(permutated_points)
            point = permutated_points(i,:); 
            closest_prototype_index = find_closest_prototype(point, prototypes);
            closest_prototype = prototypes(closest_prototype_index, :);
            updated_prototype = update_prototype(point, closest_prototype, n); 
            prototypes(closest_prototype_index, :) = updated_prototype;
        end
        H_VQ = calculate_error(points, prototypes);
        H_VQs(t) = H_VQ;
    end
end 

% function to find the closest prototype according to the squared euclidiean distance
function closest_prototype_index = find_closest_prototype(point, prototypes)
        %starting with the first distance as the shortest
    min_dist = (point(1)- prototypes(1,1)).^2 + (point(2) - prototypes(1,2)).^2;
    closest_prototype_index = 1;
    %iterate over protoypes
    for i = 1:length(prototypes(:,1))
        %Calculate distance using the squared euclidean distance.
        new_distance = (point(1)- prototypes(i,1)).^2 + (point(2) - prototypes(i,2)).^2;
        %if the new distance is smaller than the previous one, change it
        if new_distance < min_dist
            min_dist = new_distance;
            closest_prototype_index = i;
        end
    end
end 

% function to update the prototype position
function updated_prototype = update_prototype(point, prototype, n)
    movement = point - prototype; 
    updated_prototype = prototype + n * movement;
end 

function H_VQ = calculate_error(points, prototypes)
    H_VQ = 0;
    for j = 1:length(prototypes(:,1))
        for u = 1:length(points)
            point = points(u, :);
            closest_prototype_index = find_closest_prototype(point, prototypes);
            if j == closest_prototype_index
                H_VQ = H_VQ + (point(1)- prototypes(j,1)).^2 + (point(2) - prototypes(j,2)).^2;
            end
        end 
    end    
end 

function plot_learning_curve(H_VQs, n, K) 
    figure_name = sprintf('LC_n%d_K%d.pdf', n, K);
    fig = figure('Name',figure_name);
    h = zeros(1,1);
    h(1) = plot(NaN,NaN, '.r');
    hold on
    
    % plot value for every epoch
    plot(1:length(H_VQs), H_VQs, '-r.', 'LineWidth',1);
    
    legend(h, 'H_{VQ}');
    xlabel('t');
    ylabel('Value of H_{VQ}');
    grid
    save_plot(figure_name, fig);
end 

function plot_learning_curve_K2_K4(H_VQs_K2, H_VQs_K4, n) 
    figure_name = sprintf('LC_n%d_K2_K4.pdf', n);
    fig = figure('Name',figure_name);
    h = zeros(2,1);
    h(1) = plot(NaN,NaN, '.r');
    hold on
    h(2) = plot(NaN,NaN, '.b');

    
    % plot value for every epoch
    plot(1:length(H_VQs_K2), H_VQs_K2, '-r.', 'LineWidth',1);
    plot(1:length(H_VQs_K4), H_VQs_K4, '-b.', 'LineWidth',1);
    
    legend(h, 'H_{VQ} for K = 2', 'H_{VQ} for K = 4');
    xlabel('t');
    ylabel('Value of H_{VQ}');
    grid
    save_plot(figure_name, fig);
end 

function plot_ellbow(final_H_VQ, n) 
    figure_name = sprintf('Ellbow_n%d.pdf', n);
    fig = figure('Name',figure_name);
    h = zeros(1,1);
    h(1) = plot(NaN,NaN, '.r');
    hold on
    
    % plot value for every epoch
    plot(1:length(final_H_VQ), final_H_VQ, '-r.', 'LineWidth',1);

    legend(h, 'Final H_{VQ}');
    xlabel('K');
    ylabel('Final value of H_{VQ}');
    grid
    save_plot(figure_name, fig);
end 

function plot_best_learning_rate(lr_HVQ)
    figure_name = sprintf('Best_learning_rate.pdf');
    fig = figure('Name',figure_name);
    h = zeros(1,1);
    h(1) = plot(NaN,NaN, '.r');
    hold on
    
    % plot value for every epoch
    plot(lr_HVQ(:,1), lr_HVQ(:,2), '-r.', 'LineWidth',1);

    legend(h, 'Final H_{VQ}');
    xlabel('learning rate');
    ylabel('Final value of H_{VQ}');
    grid
    save_plot(figure_name, fig);
end 

function plot_best_t_max(t_max_H_VQ)
    figure_name = sprintf('Best_t_max.pdf');
    fig = figure('Name',figure_name);
    h = zeros(1,1);
    h(1) = plot(NaN,NaN, '.r');
    hold on
    
    % plot value for every epoch
    plot(t_max_H_VQ(:,1), t_max_H_VQ(:,2), '-r.', 'LineWidth',1);

    legend(h, 'Final H_{VQ}');
    xlabel('t_{max}');
    ylabel('Final value of H_{VQ}');
    grid
    save_plot(figure_name, fig);
end


% function to save the plot
function save_plot(name, fig)
    set(fig, 'PaperPosition', [0 0 20 20]);
    set(fig, 'PaperSize', [20 20]);
    saveas(fig, name);
end 
