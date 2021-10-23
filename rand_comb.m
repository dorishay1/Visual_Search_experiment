function [RandBlocks] = rand_comb(Set_Size,Exp_Type,Elements)
%This function preduce all the possible blocks without reapets in a random
%sort.
%exp_type and Target_Elment are converted to numbers according to their
%index
%Set size - vector with set size as numbers, each number represents the number
%of stimulus in the current experiment.
%Exp_Type - the input can be numbers or strings, each number represents a diffrent type of
%experiment
%Elements - the input can be numbers or strings, each number represents a diffrent target
%elment.

Exp_Type = 1:length(Exp_Type);

All_Blocks = combvec(Set_Size,Exp_Type);
cols = size(All_Blocks,2);
P = randperm(cols);
RandBlocks = All_Blocks(:,P);

Elements = 1:length(Elements);

TRG = [];                                                   %TRG - Target Element Vector
for i = 1:floor(length(RandBlocks)/length(Elements))
    TRG = [TRG Elements];
end

TRG = TRG(randperm(length(TRG)));

RandBlocks(3,:) = TRG;
end
