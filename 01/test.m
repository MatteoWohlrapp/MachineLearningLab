clc

load lvqdata.mat

pointsWithClasses = createClasses(lvqdata)
random_x = pointsWithClasses(randperm(size(pointsWithClasses, 1)), :)



prototypes = [2 2 1;10 10 1; 0 0 2; 1 1 2];
dataPoint = [1 1 1];
dataPoints = [1 1 1; 2 2 1];

winner  = findClosest(dataPoint, prototypes,4);

is_it_right = correspondant(dataPoint,prototypes,4);

correctness = findCorrectPCTG(dataPoints,prototypes,4,2);

function pointsWithClasses = createClasses(dataPoints)
    %create matrix with x,y values for points and the corresponding class
    pointsWithClasses = zeros(100,4);
    pointsWithClasses(:,1:2) = dataPoints;
    pointsWithClasses(1:50,3) = 1;
    pointsWithClasses(51:100,3) = 2;
end

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
%Returns boolean x, True if the number corresponds, False if not
function x = correspondant(dataPoint, prototypes, K)

    %first we find the closes of the members to the datapoint.
    winner = findClosest(dataPoint, prototypes,K);
    
    %Does it match?
    if dataPoint(3) == prototypes(winner,3)
        x = true;
    else
        x = false;
    end 

end 


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


%Note = Functions do not change the value of previously assigned variables
%outside of themselves.s
function [a,b] = changeA(a)
    a = a + 10;
    b = "fart";
end 