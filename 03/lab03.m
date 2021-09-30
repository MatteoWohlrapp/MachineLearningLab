edit lab03

load lvqdata.mat

% N: dimensions of input vector: 100x2 
% P: number of examples: 100
% K: number of prototypes 
K_1 = 2;
K_2 = 4; 
K_3 = 8; 
% n: learning rate 
n = 0.002; 
% t_max : maximum number of epochs
t_max = 300;
% factor to consider change in epoch constant
constant = 0.01;

% variable containing dataPoints
pointsWithClasses = classify_data(lvqdata);

prototypesA = choosePrototype(pointsWithClasses, K_1);
prototypesB = choosePrototype(pointsWithClasses, K_2);
prototypesC = choosePrototype(pointsWithClasses, K_3);

do_training(pointsWithClasses, prototypesA, K_1, n, t_max);
do_training(pointsWithClasses, prototypesB, K_2, n, t_max);
do_training(pointsWithClasses, prototypesC, K_3, n, t_max);



function cross_validation(lvqdata, m)
    % function to shuffle the data 
    %TODO
    randomized_data = randomize_data(classify_data(lyqdata));
    train_errors = zeros(1,1);
    test_errors = zeros(1,1);
    for k = 1:5
        for i = 1:m 
            %TODO
            [train_data, test_data] = split_data(randomized_data, m, i);
            %TODO: adjust do_training
            prototypes = choosePrototype(train_data,(2*k));
            [train_error, prototypes] = do_training(); 
            %TODO
            test_error = validation(prototypes, test_data);
            train_errors(k,i) = train_error; 
            test_errors(k,i) = test_error;
        end
    end
    %TODO
    plot_cross_validation(train_errors, test_errors);
end

%Creates a randomized colum of data.
function new_datapoints = randomize_data(dataPoints)
  %randperm randomize the rows while keeping the columns intact. 
  %Source: https://www.mathworks.com/matlabcentral/answers/30345-swap-matrix-row-randomly
  new_datapoints = dataPoints(randperm(size(dataPoints, 1)), :);
end 

% Copying values for x,y,correctclass and the class we assume based on the 
%distance to the prototypes.
% We distinguish between two classes: 1 and 2 
function pointsWithClasses = classify_data(dataPoints)
    pointsWithClasses = zeros(100,4);
    pointsWithClasses(:,1:2) = dataPoints;
    pointsWithClasses(1:50,3) = 1;
    pointsWithClasses(51:100,3) = 2;
end

function prototypes = choosePrototype(pointsWithClasses, K)
    prototypes = zeros((2*K),3);
    
    for i = 1:(2*k)
        %if Class "1"
        if i <= K
            correct = 0;
            %Search until we find a suitable value
            while correct == 0
                randomIndex = floor((length(pointsWithClasses)).*rand(1,1)+1);
                if pointsWithClasses(randomIndex,3) == 1
                    testsize_a = unique(prototypes,"rows");
                    dummy = prototypes(i,:);
                    prototypes(i,:) = pointsWithClasses(randomIndex,(1:3));
                    testsize_b = unique(prototypes,"rows");
                    if length(testsize_b) > length(testsize_a)
                        correct = 1;
                    else
                        prototypes(i,:) = dummy;
                    end
                end 
            end
        %if Class "2"
        else
             correct = 0;
            %Search until we find a suitable value
            while correct == 0
                randomIndex = floor((length(pointsWithClasses)).*rand(1,1)+1);
                if pointsWithClasses(randomIndex,3) == 2
                    testsize_a = unique(prototypes,"rows");
                    dummy = prototypes(i,:);
                    prototypes(i,:) = pointsWithClasses(randomIndex,(1:3));
                    testsize_b = unique(prototypes,"rows");
                    if length(testsize_b) > length(testsize_a)
                        correct = 1;
                    else
                        prototypes(i,:) = dummy;
                    end
                end 
            end
        end

    end 
end

% Trains the prototypes based on squared euclidian distance. 
% Also plots the results in graphs at the end of the epochs
function [error, resultingPrototypes] = do_training(dataPoints, prototypes, K, learning_K, t_max)
    
    %Set the randomized prototypes and the start of the errors.
    resultingPrototypes = prototypes;
    errors = [];
    t = 1;
    
    %greater while loop, will only stop if t_max is reached.
    %For our case, t_max = 100
    while t <= t_max
    
    %randperm randomize the rows while keeping the columns intact. 
    %Source: https://www.mathworks.com/matlabcentral/answers/30345-swap-matrix-row-randomly
        randomDataPoints = dataPoints(randperm(size(dataPoints, 1)), :);
        %run for every data point
        for i = 1:length(dataPoints)
            %first find the closest prototype to the point:
            prototypeIndex = findClosestPrototype(randomDataPoints(i,:),resultingPrototypes, K);
            %now we find the direction in which we update the value
            direction = psiFunction(randomDataPoints, i, resultingPrototypes, K);
            
            %if the class is correct, we bring the prototype closer
            if direction == 1
                %Find the direction
                movement_dir = randomDataPoints(i,1:2) - resultingPrototypes(prototypeIndex,1:2);
                resultingPrototypes(prototypeIndex,1:2) = resultingPrototypes(prototypeIndex,1:2) + (movement_dir * learning_K);   
            %otherwise we increase the distance 
            else 
                %Find the direction
                movement_dir = randomDataPoints(i,1:2) - resultingPrototypes(prototypeIndex,1:2);
                resultingPrototypes(prototypeIndex,1:2) = resultingPrototypes(prototypeIndex,1:2) - (movement_dir * learning_K);    
            end 
            
        end 
        
      
        %with the "training" done, we calculate the error percentage
        errors(t) = 1 - findCorrectPCTG(dataPoints, resultingPrototypes, K, 100);
       
        t = t + 1;    
        
    end  
    error = errors(t_max);
    %plotLearningCurve(errors, t, K);
    %plotPrototypesAndPoints(randomDataPoints, resultingPrototypes, K);
end

%finds the prototype closest to the dataPoint
%returns only the index of the closest prototype.
function closestPrototype = findClosestPrototype(dataPoint, prototypes, K)
    %starting with the first distance as the shortest
    minDist = (dataPoint(1)- prototypes(1,1)).^2 + (dataPoint(2) - prototypes(1,2)).^2;
    closestPrototype = 1;
    %iterate over protoypes
    for i = 1:K
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
%returns 1 if the class matches, 0 if not
function x = psiFunction(dataPoints, index, prototypes, K)

    %first we find the closes of the members to the datapoint.
    winner = findClosestPrototype(dataPoints(index,:), prototypes,K);
    
    %Check if the prototype is the closest
    if dataPoints(index, 3) == prototypes(winner,3)
        x = 1;
    else
        x = 0;
    end 
end 

%Finds the percentage of the dataPoints that are correct with the given prototypes
%Values from 0.0 to 1.0 
function correctClassifications = findCorrectPCTG(dataPoints, prototypes, K, size)
    
    curCorrect = 0;
    
    for i = 1:size
        %If your class is the right one, add to curError.
        if psiFunction(dataPoints, i,prototypes,K)
            curCorrect = curCorrect + 1; 
        end 
    end 
    
    correctClassifications = curCorrect/size;

end 

% returns 1 if the number of errors is constant, 0 if not
% we declare the number of errors constant if the last 50 values are equal
function isConstant = isConstant(errors) 
    maxIndex = length(errors);
    isConstant = 1; 
    if maxIndex <= 50 
        isConstant = 0;
        return;
    end
    minIndex = maxIndex - 50 ; 
    
    targetValue = errors(minIndex); 
    
    for index = minIndex:maxIndex
        if errors(index)~= targetValue
            isConstant = 0; 
            return
        end
    end 
end 


% plots the classes and corresponding prototypes with
% different colours
function plotPrototypesAndPoints(pointsWithClasses, prototypes, numberOfPrototypes)
    %add the assumed class to the 4th column of the method
    for index = 1:length(pointsWithClasses) 
        closestPrototype = findClosestPrototype(pointsWithClasses(index,:), prototypes, numberOfPrototypes);
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


