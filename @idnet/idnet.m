% MATLAB class definition of dynamic networks for identification
%
% Copyright (c) 2014-2017, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last update on 05 Aug 2016


classdef idnet < handle

    % === List of Private Class Members ===
    properties (SetAccess = private)
    %   A : px1 cell, A{i} = [1, a_{i1}, ..., a_{i n_i^a}], diagonal
    %   By: pxp cell, By{i,j} = [b^y_{ij1}, ..., b^y_{ij n_ij^by}]
    %   Bu: pxp cell, Bu{i,j} = [b^u_{ij1}, ..., b^u_{ij n_ij^bu}]
    %   C : px1 cell, C{i} = [1, c_{i1}, ..., c_{i n_i^c}], diagonal
    %   - model: A(q)y(t) = By(q)y(t) + Bu(q)u(t) + C(q)e(t)
        A
        By
        Bu
        C
        BoolNet = struct([])
        BoolNetProb = 'NA'
        NoiseVariance = 'eye'
        PredError
        Report
    end


    % === List of Public Class Members ===
    properties (SetAccess = public)
        Ts
        Units = 'sec'
        Variable = 'q^-1'
        GroundTruth = []
        ThresholdZero = 1e-6
    end


    % ==== List of Public Functions ====
    methods

        % ==== List of Fundamental Methods ====

        % Constructor
        function net = idnet(p,m)
            net.A = cell(1,p);
            net.By = cell(p,p);
            net.Bu = cell(p,m);
            net.C = cell(1,p);
        end

        % Propterty Access Methods
        function net = set.A(net, val)
            net.A = val;
        end
        function net = set.By(net, val)
            net.By = val;
        end
        function net = set.Bu(net, val)
            net.Bu = val;
        end
        function net = set.C(net, val)
            net.C = val;
        end


        % ==== List of SysId Functions ====

        % Inference via armax, arx
        %   - sys : idnet object
        %   - data: iddata object
        sys = armax(net, data, norder, algorithm)
        sys = arx(net, data, norder, algorithm)


        % ==== List of Model Validation Functions ====

        % boolnets = multi_runs(net, data, ModelOrders, solver_config)
        boolnets = multi_runs(net, data, varargin)


        % ==== List of Network Operation Functions ====

        % Generate Boolbean network
        boolnet = bool(net, dynet_type, zero)

        % Convert ident Object to dsf Object
        dynetObj = dsf(net)


        % ==== List of Data Visualization Functions ====

        % Reloading "plot" methods
        function pg = plot(net)
            if isempty(net.BoolNet)
                net.BoolNet = net.bool(net.ThresholdZero);
            end
            % plotting diagraph
            Gq = digraph(net.BoolNet.Q');
            pg = plot(Gq, 'Layout', 'force');
            pg.NodeColor = 'red';
        end

        % Draw performance curves, e.g. P-R, ROC, after multi_run()
        fig_hl = perfcurve(net, boolNets, curve_type, varParaVector)
    end


    % ==== List of Private Functions for Solvers ====
    methods (Access = private)
        % to be finished

    end

end