edit lab02

linreg = load('linreg.mat');


calculateWeights(linreg, 40)


function w = calculateWeights(linreg, P) 
    % matrix of Px25 of input vectors
    X = linreg.xtrain(1:P,:); 
    % transposed input vectors of 25x500
    X_transposed = transpose(X);
    % matrix of 1x500 -> transposing it for right format of Px1
    Y = transpose(linreg.ytrain(1:P)); 
    % we are computing [X^T*X]^(-1)*X^T*Y -> vector of 25x1 of weights 
    X_X_transposed = X_transposed * X; 
    X_X_transposed_inverse = pinv(X_X_transposed);
    
    % calculating weights
    w = X_X_transposed_inverse * X_transposed * Y(1:P);
end 

function E_train = training_error(linreg, w, P)
    E_train = 0;
    

end 




