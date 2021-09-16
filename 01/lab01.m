edit lab01

load lvqdata.mat


firstClass = lvqdata(1:50,1:2);
secondClass = lvqdata(51:100,1:2);

prototypeOne = firstClass(1,:);
prototypeTwo = secondClass(1,:);

plot(firstClass(2:50,1), firstClass(2:50,2), 'r.') 
hold on 
plot(secondClass(2:50,1), secondClass(2:50,2), 'b.')
hold on 
plot(prototypeOne(1,1), prototypeOne(1,2), 'rs')
hold on
plot(prototypeTwo(1,1), prototypeTwo(1,2), 'bs')
hold off

%This part was made by Lucas
%classmat contains the class we will iterate over
%It has values x,y,correctclass and guessedclass.
classmat = createclass(lvqdata);
x = choose_prototype_1(classmat)

function x = createclass(mat)
    %create classes matrix
    x = zeros(100,4);
    x(:,1:2) = mat;
    x(1:50,3) = 1;
    x(51:100,3) = 2;

end

%N = dimension of input vector, 100x3 - x,y,class,"answer"
%number of examples, P = 100
%Number of protoypes = 2 or 4
%Training rate - 0.002

%t = 200 / 300

%This functions randomly choses a prototype amongst the classs
%This one choses one per class
function x = choose_prototype_1(mat_class)
    a = 1;
    b = 50;
    c = 51;
    d = 100;
    k1 = floor((b-a).*rand(1,1) + a);
    k2 = floor((d-c).*rand(1,1) + c);
    
    x = zeros(2,3);
    
    x(1,:) = mat_class(k1,1:3);
    x(2,:) = mat_class(k2,1:3);    
end


function x = epoch(classmatrix,prototypes)

end