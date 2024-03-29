
% class holding the methods to perform the DBSCAN
classdef DBSCAN 
    methods(Static)
        % function to perform the DBSCAN
        function Ps = dbscan(D, eps, Min_Pts) 
            C = 1; 
            Ps = zeros(length(D(:,1)),1); 
        
            index = exists_unvisited_point(Ps); 
            while index ~= -1 
                Ps(index) = C;
                Neighbor_Pts = region_query(D, D(index, :), eps); 
                if length(Neighbor_Pts) < Min_Pts
                    Ps(index) = -1;
                else 
                    new_Ps = expand_cluster(D, Ps, Neighbor_Pts, C, eps, Min_Pts); 
                    Ps = new_Ps; 
                    C = C+1;
                end
                index = exists_unvisited_point(Ps);
            end
        end 
    end 
end 

% function finding the first point which is not visited yet, if there is no
% such point, the function returns -1
function index = exists_unvisited_point(P)
    index = -1; 
    for i = 1:length(P)
        if P(i) == 0
            index = i; 
            break; 
        end
    end
end 

% function returning all points in the eps radius of P
% input: Data vectors D, point P and eps 
% output: indices of the neighbor vectors
function Neighbor_Pts = region_query(D, P, eps)
    Neighbor_Pts = []; 
    index = 1; 
    for i = 1:length(D(:,1))
        dist = squared_euclidian_distance(D(i, :), P);
        if dist < eps
            Neighbor_Pts(index) = i; 
            index = index + 1; 
        end
    end
end 

% function to expand the cluster 
% input: Data vectors D, clusters of points Ps, Neighbor points, Cluster
% index C, eps and MinPts
% output: modified cluster information for points
function new_Ps = expand_cluster(D, Ps, Neighbor_Pts, C, eps, Min_Pts)
    new_Ps = Ps;
    index = 1; 
    while index <= length(Neighbor_Pts)
        if new_Ps(Neighbor_Pts(index)) == 0
            new_Ps(Neighbor_Pts(index)) = C; 
            new_Neighbor_Pts = region_query(D, D(Neighbor_Pts(index), :), eps); 
            if length(new_Neighbor_Pts) >= Min_Pts
                Neighbor_Pts = [Neighbor_Pts, new_Neighbor_Pts];
            end
        end 
        index = index + 1; 
    end
end 

% function returning the squared euclidian distance for 2 points 
% input: 2 vectors in a 2D-plane, 1st value x, 2nd value y
function distance = squared_euclidian_distance(point_1, point_2) 
    distance = sqrt(sum((point_1 - point_2).^2)); 
end 