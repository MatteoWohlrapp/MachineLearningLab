
load data_clustering.csv

A = zeros(200, 200,2);

A(:,1,:) = data_clustering;
A(:,2:end,:) = -1;

B = A(1,1,:)

A(1,2,:) = B;


%result_ = make_dist_matrix(firstmat,secondmat)

function dist_mat = make_dist_matrix(dists_1, dists_2)

    dist_mat = zeros(length(dists_1(:,1)),length(dists_2(:,1)));

    for i = 1:length(dists_1(:,1))
        for j = 1:length(dists_2(:,1))
        
            dist_mat(i,j) = sum((dists_1(i,:)-dists_2(j,:)).^2)
        end 
    end

end 