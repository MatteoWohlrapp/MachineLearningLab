a = 10

[a,b] = changeA(a);
[a,b] = changeA(a);
a
b

%Note = Functions do not change the value of previously assigned variables
%outside of themselves.s
function [a,b] = changeA(a)
    a = a + 10;
    b = "fart";
end 