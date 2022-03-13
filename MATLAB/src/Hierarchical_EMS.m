function [x, point_seg, point_outlier] = Hierarchical_EMS(point, varargin)

[para] = parseInputArgs(point, varargin{:});

x = cell(1, para.MaxLayer);
point_seg = cell(1, para.MaxLayer);
point_outlier = cell(1, para.MaxLayer);
point_seg{1}{1} = point;

minDistance = para.MinDistance;
minPoints = para.MinPoints;

for h = 1 : para.MaxLayer
    
    point_seg{h + 1} = cell(0);
    point_outlier{h} = cell(0);
    
    for c = 1 : size(point_seg{h}, 2)
        
        [x_raw, p_raw] = EMS(point_seg{h}{c}, ...
            'OutlierRatio', para.OutlierRatio, ...
            'MaxIterationEM', para.MaxIterationEM, ...
            'ToleranceEM', para.ToleranceEM, ...
            'RelativeToleranceEM', para.RelativeToleranceEM, ...
            'MaxOptiIterations', para.MaxOptiIterations, ...
            'Sigma', para.Sigma, ...
            'MaxSwitch', para.MaxSwitch, ...
            'AdaptiveUpperBound', para.AdaptiveUpperBound, ...
            'Rescale', para.Rescale);
        
        point_previous = point_seg{h}{c};
        x{h}{c} = x_raw;
        outlier = point_seg{h}{c}(:, p_raw < 0.1); 
        point_seg{h}{c} = point_seg{h}{c}(:, p_raw > 0.1); 
        
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

% ------------------parsing  input arguments-------------------------------
    function [para] = parseInputArgs(point, varargin)
        
        [n_row, n_col] = size(point);
        % check for dimension of input
        if n_row ~= 3
            error('Input point cloud shoud be an array of 3 x N.')
        end
        % check for minumum points
        if n_col < 11
            error('Number of points less than 11.')
        end
        
        % set input parser
        defaults = struct('OutlierRatio', 0.9, ...
            'MaxIterationEM', 20, ...
            'ToleranceEM', 1e-3, ...
            'RelativeToleranceEM', 2e-1, ...
            'MaxOptiIterations', 2, ...
            'Sigma', 0.3, ...
            'MaxSwitch', 2, ...
            'AdaptiveUpperBound', true, ...
            'Rescale', false, ...
            'MaxLayer', 3, ...
            'MinDistance', 1, ...
            'MinPoints', 60);
        
        parser = inputParser;
        parser.CaseSensitive = false;
        
        parser.addParameter('OutlierRatio', defaults.OutlierRatio, ...
            @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'nonnegative' , '>=', 0, '<', 1}));
        parser.addParameter('MaxIterationEM', defaults.MaxIterationEM, ...
            @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'integer', 'positive'}));
        parser.addParameter('ToleranceEM', defaults.ToleranceEM, ...
            @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'nonnegative'}));
        parser.addParameter('RelativeToleranceEM', defaults.RelativeToleranceEM, ...
            @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'nonnegative'}));
        parser.addParameter('MaxOptiIterations', defaults.MaxOptiIterations, ...
            @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'integer', 'positive'}));
        parser.addParameter('Sigma', defaults.Sigma, ...
            @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'nonnegative'}));
        parser.addParameter('MaxSwitch', defaults.MaxSwitch, ...
            @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'integer', 'nonnegative'}));
        parser.addParameter('AdaptiveUpperBound', defaults.AdaptiveUpperBound, @islogical);
        parser.addParameter('Rescale', defaults.Rescale, @islogical);
        parser.addParameter('MaxLayer', defaults.MaxLayer, ...
            @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'integer', 'positive'}));
        parser.addParameter('MinDistance', defaults.MinDistance, ...
            @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'integer', 'positive'}));
        parser.addParameter('MinPoints', defaults.MinPoints, ...
            @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'integer', 'positive'}));
        
        parser.parse(varargin{:});
        para = parser.Results;
    end
end