load data_clustering.csv;

linkage_and_show('single',data_clustering);
linkage_and_show('average',data_clustering);
linkage_and_show('complete',data_clustering);
linkage_and_show('ward',data_clustering);

%First, we need to create a matrix of classes.
%We do this using cell arrays.

%every point is a class at the start.
classes = [1:200]; %simple vector for the classes

%Class members vector. Starts with one value for each class (which is why
%it is N_classesxN_members_classX xy_values
class_members = zeros(200,1,2)

class_members(:,1,:) = data_clustering



%function for printing, and saving the clustering data as defined.
function x = linkage_and_show(method,data)
    
    x = linkage(data,method,'squaredeuclidean');
        
    fig = figure('Name',method);
    dendrogram(x)
    yline(median([x(end-1,3),x(end,3)]),'--b','K = 2')
    yline(median([x(end-2,3),x(end-1,3)]),'--b','K = 3')
    yline(median([x(end-3,3),x(end-2,3)]),'--b','K = 4')
    
    saveas(fig,method + ".pdf");

end    

%Distance functions
%Big Distance function:
%Get two (Nx2) matrix and get the resulting distance matrix (N1XN2).

%We can use this in single linkage, by choosing the smallest value found
%here. Same goes for complete linkage, this time, choosing the highest
%value found here. Average linkage can be done here too by simply adding
%the whole lower part of the matrix and averaging it by the 
%number of connections (|C1| * |C2|).
%Ward linkage probably cannot be done here.
function dist_mat = make_dist_matrix(dists_1, dists_2)

    %tecnically we only used the lower side of this matrix.
    dist_mat = zeros(length(dists_1(:,1)),length(dists_2(:,1)))

    for i = 1:length(dists_1(:,1))
        for j = 1:length(dists_2(:,1))
           
            dist_mat(i,j) = sum(dists_1(i,:)-dists_2(j,:).^2);
        end 
    end

end 


%Iteration
%Start with each cluster containing one data point - Class_and_members
%Compute the proximity matrix, for single,complete and average this means
%using make_dist_matrix, each finding a value.
%Create the whole matrix by doing this for each pair of classes.
%Find the smallest distance and fuse the members, that is, transfer from
%one class, to another and destroy the previous abandoned class

function final_result,final_members = iterator_single(class_members,K,method)

   
    final_members = class_members;
    %distance_matrix = zeros(length(classes),length(classes));

    while length(class_members(:,1,1)) > K
       
        %find distance for each member based on method.
        %Matrix is ClassesXClasses.
        %We can either create a new one each time, or modify it.
        %We will have to update the distance values anyway, since we joined
        %two different classes.
        distance_matrix = zeros(length(class_members(:,1,1)),length(class_members(:,1,1)));
        
        for i = 1:length(class_members(:,1,1))
            %We must remember to only use the lower triangle and ignore the
            %middle.
            for j = 1:length(class_members(:,1,1))
                %The method we`ll use can be passed into the function
                %This means we dont have to create a whole iteration cicle
                %for single, complete etc.
                distance_matrix(i,j) = method(final_members(i,:,:),final_members(j,:,:)); %Method used to find value we want. 
            end
        end 

        %Find the smallest, treat the middle values ixi (which will probably be
        %0) as to make them uneligible 
        m_i, m_j = findsmallest(distance_matrix);

        %With them found, we give the values in class j to class i.
        A = final_members(m_i,:,:);
        B = final_members(m_j,:,:);
        final_members(m_i,:,:) = [A;B]; %Values are added to class i
        final_members(m_j,:,:) = []; %we delete values in class j, deleting j.


    end 
end 


