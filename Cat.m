classdef Cat
    properties
        id
        weight = ones(2,1);
    end
    
    methods
        function obj = Cat(id)
            obj.id = id;
        end
    end
end
% % Create a new Cat object with id 1
% cat1 = Cat(1);

% % Create a new Cat object with id 2
% cat2 = Cat(2);
% cat1 = Cat(1);
