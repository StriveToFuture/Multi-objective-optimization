classdef DLF_ARVEA < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained/none>
% DualArchive - A Dual-Archive Co-evolutionary Algorithm with Dynamic Leader-Follower Selection.
% NV --- 100 --- Number of ReferenceVectors 


%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlATEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            NV = Algorithm.ParameterSet(100);
            fr = 0.1;

            %% Generate reference vectors for Leader Archive
            [V0, ~] = UniformPoint(NV, Problem.M);
            V = V0;
            
            %% Generate random population and initialize archives
            Population = Problem.Initialization();
            [Population, FrontNo, CrowdDis] = DA_EnvironmentalSelection(Population, Problem.N, V);

            %% Optimization

            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring = OperatorGA(Problem, Population(MatingPool));
                [Population, FrontNo, CrowdDis] = DA_EnvironmentalSelection([Population, Offspring], Problem.N, V);
                if ~mod(ceil(Problem.FE/Problem.N),ceil(fr*Problem.maxFE/Problem.N))
                    V = AdaptReferenceVectors(Population, V0);
                end
            end
        end
    end
end