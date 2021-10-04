edit lab03

load lvqdata.mat
 

m = 5; 
cross_validation(lvqdata, m)
m=10;
cross_validation(lvqdata, m)


function cross_validation(lvqdata, m)
    % n: learning rate 
    n = 0.002; 
    % t_max : maximum number of epochs
    t_max = 100;
    % function to shuffle the data 
    randomized_data = randomize_data(classify_data(lvqdata));
    train_errors = zeros(1,1);
    test_errors = zeros(1,1);
    
    train_error_prototypes = zeros(5,100);
    
    
    for K = 1:5
        for i = 1:m 
            [train_data, test_data] = split_data(randomized_data, m, i);
            prototypes = choose_prototypes(train_data,K);
            [train_error, prototypes] = do_training(train_data, prototypes, n, t_max); 
            test_error = validation(prototypes, test_data);
            train_errors(K,i) = train_error(length(train_error)); 
            test_errors(K,i) = test_error;
            if i == 1
                train_error_prototypes(K, :) = train_error;
            end
        end
    end
    plot_cross_validation(train_errors, test_errors);
    plot_train_errors(train_error_prototypes);
end


function test_error = validation(prototypes, test_data)
    errors = 0; 
    
    for i = 1:length(test_data) 
        prototype_index = find_closest_prototype(test_data(i,:),prototypes);
        prototype_class = prototypes(prototype_index,3);
        
        if test_data(i,3) ~= prototype_class
            errors = errors +1; 
        end
    end 
    test_error = errors / length(test_data);

end

% function to split the data according to the current iteration
function [train_data, test_data] = split_data(randomized_data, m, k)
    interval_size = length(randomized_data)/m;
    train_data = zeros(length(randomized_data)-interval_size, 4);
    test_data = zeros(interval_size, 4);
    train_index = 1;
    test_index = 1; 

    for i = 1:length(randomized_data) 
        if ((k-1) * interval_size < i) && (i <= k * interval_size)
            test_data(test_index,:) = randomized_data(i,:);
            test_index = test_index + 1; 
           
        else 
            train_data(train_index,:) = randomized_data(i,:);
            train_index = train_index + 1; 
        end
    end

end

%Creates a randomized colum of data.
function new_datapoints = randomize_data(data_points)
  %randperm randomize the rows while keeping the columns intact. 
  %Source: https://www.mathworks.com/matlabcentral/answers/30345-swap-matrix-row-randomly
  new_datapoints = data_points(randperm(size(data_points, 1)), :);
end 

% Copying values for x,y,correctclass and the class we assume based on the 
%distance to the prototypes.
% We distinguish between two classes: 1 and 2 
function points_with_classes = classify_data(data_points)
    points_with_classes = zeros(100,4);
    points_with_classes(:,1:2) = data_points;
    points_with_classes(1:50,3) = 1;
    points_with_classes(51:100,3) = 2;
end

%Randomly choses a prototype amongst the points for each class
%Each returned prototype contains 3 things = x, y, class
function prototypes = choose_prototypes(points_with_classes, K)
    prototypes = zeros((2*K),3);
    
    for i = 1:(2*K)
        %if Class "1"
        if i <= K 
            class = 1; 
        %if Class "2"
        else 
            class = 2; 
        end
        %iterations = 0;
        correct = 0;
        %Search until we find a suitable value
        while correct == 0
            randomIndex = floor((length(points_with_classes)).*rand(1,1)+1);
            
            if points_with_classes(randomIndex,3) == class
                
                testsize_a = unique(prototypes,"rows");
                
                dummy = prototypes(i,:);
                
                prototypes(i,:) = points_with_classes(randomIndex,(1:3));
                
                testsize_b = unique(prototypes,"rows");
               
                if length(testsize_b(:,1)) >= length(testsize_a(:,1))
                      correct = 1;
                    
                %Case for last member of array      
                elseif length(testsize_b(:,1)) == length(testsize_a(:,1))
                    if(any(prototypes(:,3) == 0.0))
                    else
                        correct = 1;
                    end
                else
                    prototypes(i,:) = dummy;
                end
            end 
           % iterations = iterations + 1
        end

    end 
end

% Trains the prototypes based on squared euclidian distance. 
function [error, resulting_prototypes] = do_training(data_points, prototypes, n, t_max)
    
    %Set the randomized prototypes and the start of the errors.
    resulting_prototypes = prototypes;
    error = zeros(1,1);
    t = 1;
    
    %greater while loop, will only stop if t_max is reached.
    %For our case, t_max = 100
    while t <= t_max
    
    %randperm randomize the rows while keeping the columns intact. 
    %Source: https://www.mathworks.com/matlabcentral/answers/30345-swap-matrix-row-randomly
        random_data_points = data_points(randperm(size(data_points, 1)), :);
        %run for every data point
        for i = 1:length(data_points)
            %first find the closest prototype to the point:
            prototypeIndex = find_closest_prototype(random_data_points(i,:),resulting_prototypes);
            %now we find the direction in which we update the value
            direction = psi_function(random_data_points, i, resulting_prototypes);

            movement_dir = random_data_points(i,1:2) - resulting_prototypes(prototypeIndex,1:2);
            resulting_prototypes(prototypeIndex,1:2) = resulting_prototypes(prototypeIndex,1:2) + direction * movement_dir * n;   
        end 
        
        %calculating the last error percentage at t_max
        
        error(t) = 1 - find_correct_pctg(data_points, resulting_prototypes);
        t = t + 1;    

    end
      
end  

%finds the prototype closest to the dataPoint
%returns only the index of the closest prototype.
function closestPrototype = find_closest_prototype(dataPoint, prototypes)
    %starting with the first distance as the shortest
    minDist = (dataPoint(1,1)- prototypes(1,1)).^2 + (dataPoint(1,2) - prototypes(1,2)).^2;
    closestPrototype = 1;
    %iterate over protoypes
    for i = 1:length(prototypes(:,1))
        %Calculate distance using the squared euclidean distance.
        newDist = (dataPoint(1,1)- prototypes(i,1)).^2 + (dataPoint(1,2) - prototypes(i,2)).^2;
        %if the new distance is smaller than the previous one, change it
        if newDist < minDist
            minDist = newDist;
            closestPrototype = i;
        end
        
    end

end 

%caclulates the direction in which the prototype will move
%returns 1 if the class matches, -1 if not
function x = psi_function(dataPoints, index, prototypes)

    %first we find the closes of the members to the datapoint.
    winner = find_closest_prototype(dataPoints(index,:), prototypes);
    
    %Check if the prototype is the closest and has the same class
    if dataPoints(index, 3) == prototypes(winner,3)
        x = 1;
    else
        x = -1;
    end 
end 

%Finds the percentage of the dataPoints that are correct with the given prototypes
%Values from 0.0 to 1.0 
function correct_classification = find_correct_pctg(data_points, prototypes)
    size = length(data_points);
    curCorrect = 0;
    
    for i = 1:size
        %If your class is the right one, add to curError.
        if psi_function(data_points, i,prototypes) == 1
            curCorrect = curCorrect + 1; 
        end 
    end 
    
    correct_classification = curCorrect/size;
end 

function plot_cross_validation(train_errors, test_errors)
    % naming of figure
    fig = figure('Name', 'Error bar of training error');
    %plotting
    standard_deviation = std(train_errors);
    average = mean(train_errors);
    errorbar(1:length(train_errors),average, standard_deviation,'LineWidth',2); 
    hold on
    % legend
    h = zeros(1,1);
    h(1) = plot(NaN,NaN, 'b');
    legend(h, 'Train error');
    %labeling and adjusting of axis
    xlabel('K');
    ylabel('Value of error');
    max = length(train_errors) + 0.5;
    axis([0.5,max,0.0,0.5])
    grid
    % saving file
    saved_name = sprintf('Train_errors_P%d.pdf', length(train_errors(1, :)));
    save_plot(saved_name, fig)

     % naming of figure
    fig = figure('Name', 'Error bar of training error');
    %plotting
    standard_deviation = std(test_errors);
    average = mean(test_errors);
    errorbar(1:length(train_errors),average, standard_deviation, 'LineWidth',2);
    hold on
    % legend
    h = zeros(1,1);
    h(1) = plot(NaN,NaN, 'b');
    legend(h, 'Test error');
    %labeling and adjusting of axis
    xlabel('K');
    ylabel('Value of error');
    axis([0.5,max,0.0,0.5])
    grid
    % saving file
    saved_name = sprintf('Test_errors_P%d.pdf', length(train_errors(1, :)));
    save_plot(saved_name, fig)
end

% plots the learning curve for a given array of normalized errors
function plot_train_errors(train_error_prototypes)
    fig = figure('Name','Learning curve with differnt amounts of prototypes');
    h = zeros(2,1);
    colours = ['r', 'g', 'b', 'y', 'c'];
    
    % plot value for every prototype
    for i = 1:5 
        plot(1:100, train_error_prototypes(i, :), colours(i), 'LineWidth',1);
        hold on
        h(i) = plot(NaN,NaN, colours(i));
    end
    legend(h, '1 prototype', '2 prototypes', '3 prototypes', '4 prototypes', '5 prototypes');
    xlabel('t');
    ylabel('Value of error');
    axis([0.0,100.0,0.0,0.5])
    grid
    save_plot('Different_prototypes.pdf', fig);
end 

function save_plot(name, fig)
    set(fig, 'PaperPosition', [0 0 20 20]);
    set(fig, 'PaperSize', [20 20]);
    saveas(fig, name);
end 
