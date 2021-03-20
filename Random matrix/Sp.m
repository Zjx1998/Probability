function [ g ] = Sp( v )
% This function returns the spacing between eigenvalue which is calculated
% parallelly


    g     = v(2 : length(v)) - v(1 : length(v) - 1);
    
    
end