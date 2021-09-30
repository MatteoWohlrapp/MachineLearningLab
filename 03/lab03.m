edit lab03

load lvqdata.mat
 

m = 5; 
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
    for K = 1:5
        for i = 1:m 
            [train_data, test_data] = split_data(randomized_data, m, i);
            prototypes = choose_prototypes(train_data,K);
            [train_error, prototypes] = do_training(train_data, prototypes, n, t_max); 
            test_error = validation(prototypes, test_data);
            train_errors(K,i) = train_error; 
            test_errors(K,i) = test_error;
        end
    end
    %TODO
    plot_cross_validation(train_errors, test_errors);
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
    train_index = 0;
    test_index = 0; 
    
    for i = 1:length(randomized_data) 
        if ((k-1) * interval_size < i) && (i <= k * interval_size)
            train_data(train_index) = randomized_data(i);
            train_index = train_index + 1; 
        else 
            test_data(test_index) = randomized_data(i);
            test_index = test_index + 1; 
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
function prototypes = choose_prototypes(points_with_classes, number_of_prototypes)
    prototypes = zeros(number_of_prototypes,3);
    %divide by 2 since we have 2 classes
    numberOfPrototypesPerClass = number_of_prototypes / 2; 
    for i = 1:number_of_prototypes
        if i <= numberOfPrototypesPerClass 
            randomIndex = floor((49).*rand(1,1) + 1);
            prototypes(i,:) = points_with_classes(randomIndex,(1:3));
        else
            randomIndex = floor((49).*rand(1,1)+ 51);
            prototypes(i,:) = points_with_classes(randomIndex,(1:3));
        end 
    end   
end

% Trains the prototypes based on squared euclidian distance. 
function [error, resulting_prototypes] = do_training(data_points, prototypes, n, t_max)
    
    %Set the randomized prototypes and the start of the errors.
    resulting_prototypes = prototypes;
    error = 0;
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
        if t == t_max
            error = 1 - find_correct_pctg(data_points, resulting_prototypes, 100);
        end
      
        t = t + 1;    
        
    end  
end

%finds the prototype closest to the dataPoint
%returns only the index of the closest prototype.
function closestPrototype = find_closest_prototype(dataPoint, prototypes)
    %starting with the first distance as the shortest
    minDist = (dataPoint(1)- prototypes(1,1)).^2 + (dataPoint(2) - prototypes(1,2)).^2;
    closestPrototype = 1;
    %iterate over protoypes
    for i = 1:length(prototypes)
        %Calculate distance using the squared euclidean distance.
        newDist = (dataPoint(1)- prototypes(i,1)).^2 + (dataPoint(2) - prototypes(i,2)).^2;
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
    
    %Check if the prototype is the closest
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
        if psi_function(data_points, i,prototypes)
            curCorrect = curCorrect + 1; 
        end 
    end 
    
    correct_classification = curCorrect/size;
end 

% plots the classes and corresponding prototypes with
% different colours
function plotPrototypesAndPoints(pointsWithClasses, prototypes, numberOfPrototypes)
    %add the assumed class to the 4th column of the method
    for index = 1:length(pointsWithClasses) 
        closestPrototype = find_closest_prototype(pointsWithClasses(index,:), prototypes, numberOfPrototypes);
        pointsWithClasses(index, 4) = prototypes(closestPrototype, 3);
    end 
    %plotting class values 
    figureName = sprintf('Points with %d prototypes', numberOfPrototypes);
    figure('Name',figureName);
    for i = 1:100 
        if pointsWithClasses(i, 4) == 1
            plot(pointsWithClasses(i,1), pointsWithClasses(i,2), 'ro');
        elseif pointsWithClasses(i, 4) == 2
            plot(pointsWithClasses(i,1), pointsWithClasses(i,2), 'bo');
        end
        hold on
    end 
    
    %plotting the prototypes 
    for i = 1:numberOfPrototypes 
        if prototypes(i, 3) == 1 
            plot(prototypes(i, 1), prototypes(i, 2), 'rX');
        else 
            plot(prototypes(i, 1), prototypes(i, 2), 'bX');
        end 
        hold on
    end 
    
    h = zeros(4,1);
    h(1) = plot(NaN,NaN,'ro');
    h(2) = plot(NaN,NaN,'bo');
    h(3) = plot(NaN,NaN,'rX');
    h(4) = plot(NaN,NaN,'bX');
    legend(h, 'Class 1', 'Class 2', 'Prototype class 1', 'Prototype class 2');
    grid
end 

% plots the learning curve for a given array of normalized errors
function plotLearningCurve(errors, t_max, numberOfPrototypes) 
    figureName = sprintf('Learning curve with %d prototypes', numberOfPrototypes);
    figure('Name',figureName);
    modifiedErrors = zeros(t_max, 1);
    for i = 1:t_max 
       modifiedErrors(i) = errors(i); 
    end
    plot(1:t_max, modifiedErrors, 'r-', 1:t_max, modifiedErrors, 'r.')
    hold on
    h = zeros(2,1);
    h(1) = plot(NaN,NaN,'r.');
    h(2) = plot(NaN,NaN,'w.');
    legend(h, 'Error rate', 'X-axis: number of epochs');
    grid
end 

