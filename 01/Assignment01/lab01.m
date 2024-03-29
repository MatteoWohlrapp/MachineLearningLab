edit lab01

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
pointsWithClasses = classifyData(lvqdata);

prototypesA = choosePrototype(pointsWithClasses, K_1);
prototypesB = choosePrototype(pointsWithClasses, K_2);
prototypesC = choosePrototype(pointsWithClasses, K_3);

doTraining(pointsWithClasses, prototypesA, K_1, n, t_max);
doTraining(pointsWithClasses, prototypesB, K_2, n, t_max);
doTraining(pointsWithClasses, prototypesC, K_3, n, t_max);


% Copying values for x,y,correctclass and the class we assume based on the 
%distance to the prototypes.
% We distinguish between two classes: 1 and 2 
function pointsWithClasses = classifyData(dataPoints)
    pointsWithClasses = zeros(100,4);
    pointsWithClasses(:,1:2) = dataPoints;
    pointsWithClasses(1:50,3) = 1;
    pointsWithClasses(51:100,3) = 2;
end


%Randomly choses a prototype amongst the points for each class
%Each returned prototype contains 3 things = x, y, class
function prototypes = choosePrototype(pointsWithClasses, numberOfPrototypes)
    prototypes = zeros(numberOfPrototypes,3);
    %divide by 2 since we have 2 classes
    numberOfPrototypesPerClass = numberOfPrototypes / 2; 
    for i = 1:numberOfPrototypes
        if i <= numberOfPrototypesPerClass 
            randomIndex = floor((49).*rand(1,1) + 1);
            prototypes(i,:) = pointsWithClasses(randomIndex,(1:3));
        else
            randomIndex = floor((49).*rand(1,1)+ 51);
            prototypes(i,:) = pointsWithClasses(randomIndex,(1:3));
        end 
    end   
end

% Trains the prototypes based on squared euclidian distance. 
% Also plots the results in graphs at the end of the epochs
function doTraining(dataPoints, prototypes, K, learning_K, t_max)
    
    %Set the randomized prototypes and the start of the errors.
    resultingPrototypes = prototypes;
    errors = [];
    t = 1;
    
    %greater while loop, will only stop if t_max is reached.
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
        
        if isConstant(errors) == 1 || t == t_max
            break;  
        end 
        t = t + 1;    
        
    end  
    plotLearningCurve(errors, t, K);
    plotPrototypesAndPoints(randomDataPoints, resultingPrototypes, K);
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


