function [t, point_fit] = showSuperquadrics(x, varargin)
color = 'r';
MarkerSize = 20;
ViewAxis = [0 0];
CamRoll = 0;
ShowAxis = 0;
Render = 1;
arclength = 0.1;
FaceAlpha = 1;
AmbientStrength = 0.5;
lighting = 0;

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
    if strcmp(varargin{i}, 'Render')
        Render= varargin{i + 1};
    end
    if strcmp(varargin{i}, 'Arclength')
        arclength= varargin{i + 1};
    end
    if strcmp(varargin{i}, 'FaceAlpha')
        FaceAlpha= varargin{i + 1};
    end
    if strcmp(varargin{i}, 'AmbientStrength')
        AmbientStrength= varargin{i + 1};
    end
    if strcmp(varargin{i}, 'Light')
        lighting= varargin{i + 1};
    end
end

[point_fit] = sphericalProduct_sampling(x, arclength);
if Render == 0
    showPoints(point_fit, 'Color', color, 'MarkerSize', MarkerSize)
else
    point_fit = unique(point_fit', 'rows');
    [t,~] = MyRobustCrust(point_fit);
    trisurf(t,point_fit(:,1),point_fit(:,2),point_fit(:,3), ...
           'facecolor', color, 'LineStyle', 'none', ...
           'FaceAlpha', FaceAlpha, 'AmbientStrength', AmbientStrength)
    if lighting == 1
        light
        material dull
    end
end

axis equal
view(ViewAxis)
camroll(CamRoll)

if ShowAxis == 0
    axis off
end

hold off
