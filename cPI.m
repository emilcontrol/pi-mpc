classdef cPI
    %PI class

    properties
        I {mustBeNumeric}
    end

    methods
        function obj = cPI()
            % Constructor
            obj.I=0;
        end

        function [u,obj] = PI(obj,r,y,Kp,Ki)
            e=r-y;
            u = Kp*e + Ki*obj.I;
            obj.I = obj.I + e;
        end
    end
end