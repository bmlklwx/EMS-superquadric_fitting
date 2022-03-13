function [] = showPoints(point, varargin)

color = 'r';
MarkerSize = 10; % 20
ViewAxis = [0 0];
CamRoll = 0;
ShowAxis = 0;

for i = 1 : size(varargin, 2)
    if strcmp(varargin{i}, 'Color')
        color = varargin{i + 1};
    end
    if strcmp(varargin{i}, 'MarkerSize')
        MarkerSize = varargin{i + 1};
    end
    if strcmp(varargin{i}, 'ViewAxis')
        ViewAxis = varargin{i + 1};
    end
    if strcmp(varargin{i}, 'CamRoll')
        CamRoll = varargin{i + 1};
    end
    if strcmp(varargin{i}, 'ShowAxis')
        ShowAxis= varargin{i + 1};
    end
end

plot3(point(1, :), point(2, :), point(3, :), '.', 'Color', color,  'MarkerSize', MarkerSize)
view(ViewAxis)
camroll(CamRoll)
axis equal

if ShowAxis == 0
    axis off
end

hold off

end