function [Trials] = Trials_gen(N_Repeats,N_Symbols,Color)
%this function creats a cell array with the locations of all symbols for
%all repeations. also included color randomization and equal number of
%trials with and without target.

Trials = cell(4, N_Repeats);

%locations without overlap
for i = 1:2
    for j = 1: N_Repeats
    Trials{i,j} = [0.001 * randperm(1000, N_Symbols)];
    end
end

%color randomly
for j = 1:N_Repeats
    Trials{3,j} = Color(randi(length(Color)));
end

%with or without target equaly.
Trials(4,:) = {0};
replace_position = randperm (N_Repeats,N_Repeats/2);
Trials(4,replace_position) = {1};

end

