edit lab01

load lvqdata.mat



% N: dimensions of input vector: 100x2 
% P: number of examples: 100
% K: number of prototypes 
K_1 = 2;
K_2 = 4; 
% n: learning rate 
n = 0.002; 
% t_max : maximum number of epochs
t_max = 300;
% factor to consider change in epoch constant
constant = 0.01;


% This part was made by Lucas
% pointsWithClasses contains the points with the corresponding classes we will iterate over
% It has values x,y,correctclass and guessedclass.
% We distinguish between two classes: 1 and 2 
pointsWithClasses = createClasses(lvqdata);

prototypesA = choosePrototype(pointsWithClasses, K_1)
%prototypeB = choosePrototype(pointsWithClasses, K_2);

errorsA = zeros(1); 

modifiedPrototypesWithErrors = training(pointsWithClasses, prototypesA, K_1, constant, t_max);






function pointsWithClasses = createClasses(dataPoints)
    %create matrix with x,y values for points and the corresponding class
    pointsWithClasses = zeros(100,4);
    pointsWithClasses(:,1:2) = dataPoints;
    pointsWithClasses(1:50,3) = 1;
    pointsWithClasses(51:100,3) = 2;
end


%This functions randomly choses a prototype amongst the classs
%This one choses one per class
% numberOfPrototypes has to be even
%Each Prototype contains 3 things = x,y,class
%Therefore, information extracted from "pointsWithClasses" will 
%only be the first 3 arguments it contains.

%Changelog:
%Fixed the function so it behaves as intended, as in:
% --> Choses indexes within the desired range (the way it was previously
% implenented created floating numbers sometimes above the wanted level)
% --> Creates the proper prototype matrix. 
function prototypes = choosePrototype(pointsWithClasses, numberOfPrototypes)
    prototypes = zeros(numberOfPrototypes,3);
    
    split = numberOfPrototypes / 2; 
    for i = 1:numberOfPrototypes
        if i <= split 
            randomIndex = floor((49).*rand(1,1) + 1);
            %randomIndex = floor(mod(100 * rand, 50)) + 1;
            prototypes(i,:) = pointsWithClasses(randomIndex,[1:3]);
        else
            randomIndex = floor((49).*rand(1,1)+ 51);
            %randomIndex = floor(mod(100 * rand, 100)) + 1;
            prototypes(i,:) = pointsWithClasses(randomIndex,[1:3]);
        end 
    end 
    %a = 1;
    %b = 50;
    %c = 51;
    %d = 100;
    %k1 = floor((b-a).*rand(1,1) + a);
    %k2 = floor((d-c).*rand(1,1) + c);
    
    %prototypes(1,:) = mat_class(k1,1:3);
    %prototypes(2,:) = mat_class(k2,1:3);    
end


% pointsWithClasses is the data points with the corresponding classes 
% prototypes has the prototypes choosen randomly before 
% errors is an array containing all the error values in % from all epoches
% returns the modified prototypes in the first column and the array of
% errors in the second one
function prototypesAndErrors = training(dataPoints, prototypes, K_1, constant, t_max)
    pointsWithClasses = dataPoints;
    modifiedPrototypes = prototypes;
    numberOfPrototypes = K_1;
    numberOfErrors = 100;
    t = 1;
    errors = zeros(1); 
    
    while abs(numberOfErrors - errors(length(errors))) / 100 > constant || t > t_max
        t = t+1; 
        perm = randPerm(pointsWithClasses); 
        for i = 1:100 
            
            
        end
    end 

    
    plotPrototypesAndPoints(pointsWithClasses, modifiedPrototoypes, numberOfPrototyes)
    plotLearningCurve(errors, t-1)
    
end

%What this function does is the "iterable" part of the result.
%The actual training is done by another function that this one will call
%it returns the value of the errors (a vector of variable size) and the
%resulting prototypes. 
function [errors,resultingProtoypes] = doTraining(dataPoints, prototypes, K, learning_K, t_max)
    
    %Set the "first prototype" and the start of the errors.
    ResultingPrototypes = prototypes;
    errors = [];
    t = 1;
    
    %greater while loop, will only stop if t_max is reached.
    while t < t_max
    %Firstly we need to randomize our dataSet.
    
    %Found this line of code that is supossed to randomize the rows while
    %keeping the columns intact. 
    %LINK: https://www.mathworks.com/matlabcentral/answers/30345-swap-matrix-row-randomly
    %random_x = x(randperm(size(x, 1)), :)
    %Randomize data points
        randomDataPoints = dataPoints(randperm(size(dataPoints, 1)), :);
        %run for every data point
        for i = 1:100
            
            
        end 
        
    end  

end

%find, amongst the prototypes, which one is closest
%simple function that can be used for more than one purpose.
%Returns only the index of the closest prototype.
function winner = findClosest(dataPoint, prototypes, K)
    %some unimaginably far away distance.
    minDist = 100000;
    winner = 0;
    %iterate over protoypes
    for i = 1:K
        %Calculate distance
        %Here we are using the squared euclidean distance.
        newDist = (dataPoint(1)- prototypes(i,1)).^2 + (dataPoint(2) - prototypes(i,2)).^2;
        %if the new distance is smaller than the previous one, change it
        %and update the "winner"
        if newDist < minDist
            minDist = newDist;
            winner = i;
        end
        
    end

end 

%Finds the answer and if the answer is correct.
%Returns int x, 1 if the number corresponds, 0 if not
function x = correspondant(dataPoint, prototypes, K)

    %first we find the closes of the members to the datapoint.
    winner = findClosest(dataPoint, prototypes,K);
    
    %Does it match?
    if dataPoint(3) == prototypes(winner,3)
        x = 1;
    else
        x = 0;
    end 

end 

%Finds the percentage of the dataPoints that are correct with the given prototypes
%Values from 0.0 to 1.0 
% is the oposite of the returning 
function correct = findCorrectPCTG(dataPoints, prototypes, K, size)
    
    curCorrect = 0;
    
    for i = 1:size
        %If your class is the right one, add to curError.
        if correspondant(dataPoints(i,:),prototypes,K)
            curCorrect = curCorrect + 1; 
        end 
    end 
    
    correct = curCorrect/size;

end 




function plotPrototypesAndPoints(pointsWithClasses, prototypes, numberOfPrototypes)
    %plotting class values 
    %TODO: decide if we want to represent the actual class, or the one we
    %assigned depending on the vector
    figure
    for i = 1:100 
        if pointsWithClasses(i, 4) == 1
            plot(pointsWithClasses(i,1), pointsWithClasses(i,2), 'r.') 
        else 
            plot(pointsWithClasses(i,1), pointsWithClasses(i,2), 'bo') 
        end
        hold on
    end 
    
    %drawing the prototypes 
    for i = 1:numberOfPrototypes 
        if prototypes(i, 3) == 1 
            plot(prototypes(i, 1), prototypes(i, 2), 'rs')
        else 
            plot(prototypes(i, 1), prototypes(i, 2), 'bs')
        end 
        hold on
    end 
    hold off 
    grid
end 

function plotLearningCurve(errors, t_max) 
    figure
   for i = 1:t_max 
       plot(errors(i,1), errors(i,2), '-r.')
       hold on
   end
   hold off
   grid
end 


%firstClass = lvqdata(1:50,1:2);
%secondClass = lvqdata(51:100,1:2);

%prototypeOne = firstClass(1,:);
%prototypeTwo = secondClass(1,:);


%plot(pointsWithClasses(1:50,1), pointsWithClasses(1:50,2), 'r.') 
%hold on 
%plot(pointsWithClasses(51:100,1), pointsWithClasses(51:100,2), 'b.')
%hold on 


%plot(prototypeOne(1,1), prototypeOne(1,2), 'rs')
%hold on
%plot(prototypeTwo(1,1), prototypeTwo(1,2), 'bs')
%hold off


