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


