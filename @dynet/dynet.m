% MATLAB class definition of dynamic networks (Dynamical Structure Functions)
%
% Copyright (c) 2014-2017, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last update on 05 Dec 2016


classdef dynet < handle
    
% Private class memebers
properties (SetAccess = private)
    %   dynetObj: dynet object, structure of {Q, P, H, Ts}
    %        - Q : pxp tf
    %        - P : pxm tf
    %        - H : px1 tf, diagnoal
    %        - Ts: double
    %        - model: y(t) = Q(q)y(t) + P(q)u(t) + H(q)e(t)
    Q
    P
    H
    Ts
    Variable = 'q^-1'
end

% Public class members
properties (SetAccess = public)
    Report      % meta info
end

% Class member functions
methods
    % CONSTRUCTORS
    function obj = dynet(dsfobj)
        obj.Q = dsfobj.Q;
        obj.P = dsfobj.P;
        obj.Ts = dsfobj.Ts;
        obj.Variable = dsfobj.Variable;
    end        

    % MAIN OPERATIONS
    %   M-1: model comparation
    [fit, y, x0] = compare(net, data, plotFlag)
    
end    
    
end
