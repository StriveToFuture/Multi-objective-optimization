function [Population, FrontNo, CrowdDis] = DA_EnvironmentalSelection(Population, N, V)

    [~,M] = size(Population.objs);

    %% Non-dominated sorting
    [FrontNo, MaxFNo] = NDSort(Population.objs, Population.cons, N);
    Next = FrontNo < MaxFNo;
    
    Selected = Population(Next);
    
    % Identify the last front to be considered
    LastFrontIdx = find(FrontNo == MaxFNo);
    LastFront = Population(LastFrontIdx);
    
    % Number of individuals to choose from the last front
    Needed = N - sum(Next);
    
    if Needed > 0
        %% --- Calculate Leader and Follower Scores for the Last Front ---
        
        % 1. Leader Score (based on APD)
        Zmin_last = min(LastFront.objs, [], 1);
        TranslatedObj_last = LastFront.objs - Zmin_last;
        NormP_last = vecnorm(TranslatedObj_last, 2, 2);
        NormP_last(NormP_last==0) = 1e-6;
        % Compute cosine similarity
        Cosine_last = 1 - pdist2(TranslatedObj_last ./ NormP_last, V, 'cosine');
        [~, association_last] = max(Cosine_last, [], 2);

        Angle_last = acos(Cosine_last);
        LeaderScore = min(NormP_last .* (1 + M * Angle_last), [], 2);

        % 2. Follower Score (based on angle to already selected elites)
        FollowerScore = zeros(length(LastFront), 1);
        if ~isempty(Selected)
            Zmin_sel = min(Selected.objs, [], 1);
            Norm_Sel = vecnorm(Selected.objs - Zmin_sel, 2, 2);
            Norm_Sel(Norm_Sel==0) = 1e-6;
            SelObjs_norm = (Selected.objs - Zmin_sel) ./ Norm_Sel;
            
            Norm_Last = vecnorm(LastFront.objs - Zmin_sel, 2, 2);
            Norm_Last(Norm_Last==0) = 1e-6;
            LastObjs_norm = (LastFront.objs - Zmin_sel) ./ Norm_Last;
            
            Cosine_to_Selected = 1 - pdist2(LastObjs_norm, SelObjs_norm, 'cosine');
            FollowerScore = acos(max(Cosine_to_Selected, [], 2));
        else
            FollowerScore = CrowdingDistance(LastFront.objs, ones(1,length(LastFront)));
        end
        
        %% --- Iteratively select the remaining individuals ---
        ChooseIdx = []; 
        LastFrontAvailable = 1:length(LastFront);
        
        represented_V = false(1, size(V,1));
        if ~isempty(Selected)
             Zmin_s = min(Selected.objs,[],1);
             TranslatedObj_s = Selected.objs - Zmin_s;
             NormP_s = vecnorm(TranslatedObj_s,2,2);
             NormP_s(NormP_s==0) = 1e-6;
             Cosine_s = 1 - pdist2(TranslatedObj_s./NormP_s, V, 'cosine');
             [~,association_s] = max(Cosine_s,[],2);
             represented_V(unique(association_s)) = true;
        end
        
        for i = 1:Needed
            unrepresented_V_idx = find(~represented_V);
            
            if ~isempty(unrepresented_V_idx)
                potential_leaders = find(ismember(association_last, unrepresented_V_idx));
                potential_leaders = intersect(potential_leaders, LastFrontAvailable);

                if ~isempty(potential_leaders)
                    [~, best_leader_idx_in_subset] = min(LeaderScore(potential_leaders));
                    chosen_local_idx = potential_leaders(best_leader_idx_in_subset);
                    represented_V(association_last(chosen_local_idx)) = true;
                else
                    [~, chosen_local_idx_in_subset] = max(FollowerScore(LastFrontAvailable));
                    chosen_local_idx = LastFrontAvailable(chosen_local_idx_in_subset);
                end
            else
                [~, chosen_local_idx_in_subset] = max(FollowerScore(LastFrontAvailable));
                chosen_local_idx = LastFrontAvailable(chosen_local_idx_in_subset);
            end
            
            ChooseIdx = [ChooseIdx, chosen_local_idx];
            LastFrontAvailable(LastFrontAvailable == chosen_local_idx) = [];
        end
        Next(LastFrontIdx(ChooseIdx)) = true;
    end
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdingDistance(Population.objs, FrontNo);
end