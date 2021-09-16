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
