edit lab02

linreg = load('linreg.mat');

plot_errors(linreg);
plot_weights(linreg);

function w = calculate_weights(linreg, P) 
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

% function calculating the training error for given P
function E_train = training_error(linreg, w, P)
    Xs = linreg.xtrain;
    Ys = linreg.ytrain;
    E_train = 0;
    
    % sum formula
    for i = 1:P 
        % formula for each individual sum
        E_train = E_train + 1/2 * ((Xs(i,:) * w - Ys(i))^2);
    end
    % applying the average 
    E_train = 1/P * E_train;
end 

% function calculating the test error for P = 500
function E_test = test_error(linreg, w)
    Xs = linreg.xtest;
    Ys = linreg.ytest;
    E_test = 0;
    
    % sum formula
    for i = 1:500
        % formula for each individual sum
        E_test = E_test + ((Xs(i,:) * w - Ys(i))^2);
    end
    % applying the average 
    E_test = 1/1000 * E_test;
end 

% function plotting the training and testing error for different P 
function plot_errors(linreg)
    % initialising an array which holds all inspected P values
    p_values = zeros(1,1);
    % 50 values, from 10 to 500 with steps of 10
    for i = 1:50
        p_values(i) = i * 10;
    end 
    % initializing arrays for errors
    training_errors = zeros(1,1); 
    test_errors = zeros(1,1);
    
    % looping all p values
    for i = 1:length(p_values)
        %calculating the optimal weights
        w = calculate_weights(linreg, p_values(i));
        % calculating the training and test errors 
        tr_e = training_error(linreg, w, p_values(i));
        te_e = test_error(linreg, w);
        training_errors(i) = tr_e;
        test_errors(i) = te_e; 
    end 
    % naming of figure
    fig = figure('Name', 'Plot of training and testing errors');
    %plotting
    %plot(p_values, te_e, 'r.');
    plot(p_values,test_errors,'r. ');
    hold on 
    %plot(p_values, tr_e, 'b.');
    plot(p_values,training_errors, 'b. ');
    % legend
    h = zeros(2,1);
    h(1) = plot(NaN,NaN,'r.');
    h(2) = plot(NaN,NaN,'b.');
    legend(h, 'Test errors', 'Training errors');
    %labeling and adjusting of axis
    xlabel('P');
    ylabel('Value of error');
    axis([0,500,0.0,1.0])
    %axis([0,500,1.0,1.0])
    grid
    % saving file
    set(fig, 'PaperPosition', [0 0 25 25]);
    set(fig, 'PaperSize', [25 25]);
    saveas(fig, 'Testing_Errors.pdf');
end 

%function plotting the weights as bar graphs
function plot_weights(linreg)
    % values where we examine the resulting weights
    p_values = [30,40,50,75,100,500];
    % automated for every value in the array
    for i = 1:length(p_values) 
        %calculating optimal weight
        w = calculate_weights(linreg, p_values(i));
        
        % naming of figure
        figure_name = sprintf('Weights with P = %d', p_values(i));
        fig = figure('Name',figure_name);
        % plotting
        bar(1:25, w);
        hold on 
        % legend
        h = zeros(1,1);
        h(1) = plot(NaN,NaN,'b.');
        legend(h, 'Weight at the index w_i');
        % labeling of axis
        xlabel('Index w_i');
        ylabel('Value of weight');
        grid
        % saving file
        set(fig, 'PaperPosition', [0 0 25 25]);
        set(fig, 'PaperSize', [25 25]);
        saved_name = sprintf('W_P%d.pdf', p_values(i));
        saveas(fig, saved_name);
    end
end





