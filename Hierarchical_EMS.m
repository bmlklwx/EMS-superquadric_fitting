function [x, point_seg, point_outlier] = Hierarchical_SuperquadricsGaussian(point, para)

x = cell(1, para.maximum_layer);
point_seg = cell(1, para.maximum_layer);
point_outlier = cell(1, para.maximum_layer);
point_seg{1}{1} = point;

minDistance = para.minDistance;
minPoints = para.minPoints;

for h = 1 : para.maximum_layer
    
    point_seg{h + 1} = cell(0);
    point_outlier{h} = cell(0);
    
    for c = 1 : size(point_seg{h}, 2)
        
        [x_raw, p_raw] = SuperquadricsGaussian(point_seg{h}{c}, para);
        
        point_previous = point_seg{h}{c};
        x{h}{c} = x_raw;
        outlier = point_seg{h}{c}(:, p_raw < 0.1); % 0.1 outlier
        point_seg{h}{c} = point_seg{h}{c}(:, p_raw > 0.1); % 0.1 inlier update
        
        if sum(p_raw) < 0.8 * size(point_previous, 2) %0.9 / 0.8
                                                                 
            pc = pointCloud(outlier');
            [labels, numClusters] = pcsegdist(pc, minDistance, 'NumClusterPoints', minPoints);
            
            if numClusters >= 1                            
                for i = 1 : numClusters
                    point_seg{h + 1}{end + 1} = pc.Location(labels == i, :)';                    
                end
            end
            point_outlier{h}{end + 1} = pc.Location(labels == 0, :)'; 

        else
            point_outlier{h}{end + 1} = outlier;
        end
        
    end
end
end