
load data_clustering.csv

cellBased = {};

for i = 1:200

    cellBased{i} = data_clustering(i,:); 

end

length(cellBased)

%testing: Add infor to cell 1

manipulator = cellBased{1};

manipulator2 = cellBased{2};

manipulator = [manipulator;manipulator2]

cellBased{1} = manipulator;

length(cellBased)

cellBased{1}

cellBased(2) = [];


manipulator = cellBased{2};

manipulator2 = cellBased{3};

manipulator = [manipulator;manipulator2]

cellBased{2} = manipulator;

length(cellBased)

cellBased{2}

cellBased(3) = [];


length(cellBased)

cell1 = cellBased{1}
cell2 = cellBased{2}

dist_result = make_dist_matrix(cellBased{1},cellBased{2})
dist_result = make_dist_matrix(cellBased{1},cellBased{1})
dist_result = make_dist_matrix(cellBased{2},cellBased{1})

[I,J] = find_smallest(dist_result)

dist_result(I,J)

%Minimal value from each Column
%[M,I] = min(dist_result)

%minimun value overall
%[M,I] = min(M)

%dist_result(1,1) = 0

%[M,I] = min(dist_result)

%[M,I] = min(M)


%Ok this above works.
%indexing = class

A = zeros(200, 200,2);

A(:,1,:) = data_clustering;
A(:,2:end,:) = -1;

B = A(1,1,:);

A(1,2,:) = B;


%result_ = make_dist_matrix(firstmat,secondmat)

function dist_mat = make_dist_matrix(dists_1, dists_2)

    dist_mat = zeros(length(dists_1(:,1)),length(dists_2(:,1)));

    for i = 1:length(dists_1(:,1))
        
        for j = 1:length(dists_2(:,1))
            
            dist_mat(i,j) = sum((dists_1(i,:)-dists_2(j,:)).^2);
        end 
    end

end 


function [m_i,m_j] = find_smallest(matrix)
    
    smallest = 100000;
    m_i = 0;
    m_j = 0;

    for i = 1:length(matrix)
        i
        for j = 1:(i-1)
            j
            if matrix(i,j) < smallest
                m_i = i;
                m_j = j;
                smallest = matrix(i,j);
            end 
        end
    end
end