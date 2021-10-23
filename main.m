clear; close all; clc;

% Ex - 3

% Dor Ishay         312328339
% Liav Sommerfeld   312481468

%% settings
% task settings
Set_Size = [4 8 12 16];
Color = {'r' 'b'};
Elements =  { 'X' 'O' };
Exp_Type = ['f' 'c'];                                   % 1- feature search, 2- conjunction search
keys = {'A' 'L'};                                       % A - target present, L - target absent
N_Repeats = 30;

%Results settings
T_Max_filter = 3.75;                                    %max response time in seconds
Min_Trials_per_Block = 20;



%% Dry/Wet run (Bonus)
run_mode = 0;                                          %0- human, 1-robot
time_frame = 0.3 + 2;                                  %set average time response



%% Figure settings
h = figure;

h.Units = 'normalized';
h.Position = [0 0 1 1];
h.MenuBar = 'none';
h.Color = 'w';
h.Name = 'Visual Search Experiment';
axis off

%% Welcome Screen
welcome_screen = imread('welcome.bmp');
imshow(welcome_screen,'InitialMagnification','fit')

Press2Continue(h,' ',run_mode);


%% The Task

Rand_Blocks = rand_comb(Set_Size,Exp_Type,Elements);          %creates random blocks.
data = cell(1,length(Rand_Blocks));                           %data will be stored here - cell per Block.

%The main loop - the experiment goes according to Rand_Blocks. First we save the
%block data, than the splits into feature or conjunction search and than rather
%it is a target or a non target trail.

for N_Block = 1 :length(Rand_Blocks)
    data{N_Block} = zeros(N_Repeats,5);                   %sets the matrix in the cell for block data.
                                                          %5 = RT, ACC,Setsize,Target present,Exp_Type.
    
    %Preparing Block
    N_Symbols = Rand_Blocks(1,N_Block);                                 %set number of symbols.
    Target = Elements(Rand_Blocks(3,N_Block));                          %choose target symbol.
    
    %chosing distractor (choose one symbol from elments not included
    %target symbol.
    Block_Elements = Elements;
    Block_Elements((Rand_Blocks(3,N_Block))) = [];
    Distractor = Block_Elements(randi(length(Block_Elements)));
    
    %storing data
    data{N_Block}(:,5) = Rand_Blocks(2,N_Block);
    data{N_Block}(:,3) = N_Symbols;
    
    Trials = Trials_gen(N_Repeats,N_Symbols,Color);                     %generate trials for the block. 
    
    %Prepare screen
    text(0.5,0.5,['Your target is',(Target),'Please press spacebar to start'],...
        'HorizontalAlignment','center')
    
    Press2Continue(h,' ',run_mode);
    
    %(1)Feature Search
    if Rand_Blocks(2,N_Block) == 1
        
        for t = 1:N_Repeats
            
            Set_Color = cell2mat(Trials{3,t});
            
            %No Target
            
            if Trials{4,t} == 0                                                 
                text (Trials{1,t},Trials{2,t},Distractor,'Color',Set_Color);
                tic;
                
                if run_mode == 0                                      %Human mode
                    pause(); key = get(h,'CurrentCharacter');

                    
                elseif run_mode == 1                                  %Robot mode
                    time_sti = time_frame*rand();                     %random time according to settings.
                    key = keys(randi(length(keys)));                  %random key.
                    pause(time_sti);
                end
                
                %store the data
                data{N_Block}(t,1) = toc;
                data{N_Block}(t,2) = strcmpi(key,keys{2});
                data{N_Block}(t,4) = Trials{4,t}; 

                
            %with Target
            elseif Trials{4,t} == 1
                
                %the (n-1) first symbols location will be the distractors,
                %and the (n) symbol will be the target
                text (Trials{1,t}(1:(N_Symbols-1)),Trials{2,t}(1:(N_Symbols-1)),...
                Distractor,'Color',Set_Color);
                text (Trials{1,t}(N_Symbols),Trials{2,t}(N_Symbols),Target,...
                'Color',Set_Color);
                
                tic;
                if run_mode == 0                                      %Human mode
                    pause(); key = get(h,'CurrentCharacter');
                    
                elseif run_mode == 1                                  %Robot mode
                    time_sti = time_frame*rand();                     %random time according to settings.
                    key = keys(randi(length(keys)));                  %random key.
                    pause(time_sti);
                end
                
                
                %store data
                data{N_Block}(t,1) = toc;
                data{N_Block}(t,2) = strcmpi(key,keys{1});
                data{N_Block}(t,4) = Trials{4,t};
                
            end
            clf;
            axis off
        end
        
        
        %(2) Conjunction Search
    elseif Rand_Blocks(2,N_Block) == 2
        
        for t = 1:N_Repeats
            
            %setting colors
            Set_Color = cell2mat(Trials{3,t});
        
            dis_Color = cell2mat(Color(randi(length(Color))));         %choose random color.
            
            %if the color that have been chosen is the same as the first one the code will try to choose 
            %color once again until a different one is chosen. 
            while strcmp(Set_Color,dis_Color) == 1                     
                dis_Color = cell2mat(Color(randi(length(Color))));          
            end
            
            %No Target
            if Trials{4,t} == 0
                
                %Distractors (for example: in a set_size=4 without target, we
                %want 2 symbols (half) with the same shape and color(dis_num1 & color
                %1) and the other 2 symbols with 2nd color & shape (dis_num2
                %&color 2)
                Dis_num1 = [1:(N_Symbols/2)];
                Dis_num2 = [(N_Symbols/2)+1:(N_Symbols)];
                
                %color 1
                text (Trials{1,t}(Dis_num1),Trials{2,t}(Dis_num1),Distractor,'Color',Set_Color);
                
                %color 2
                text (Trials{1,t}(Dis_num2),Trials{2,t}(Dis_num2),Target,'Color',dis_Color);

                tic;
                
               if run_mode == 0                                       %Human mode
                    pause(); key = get(h,'CurrentCharacter');
               
               elseif run_mode == 1                                   %Robot mode
                   time_sti = time_frame*rand();                      %random time according to settings.
                   key = keys(randi(length(keys)));                   %random key.
                   pause(time_sti);
               end
                
                %store data      
                data{N_Block}(t,1) = toc;
                data{N_Block}(t,2) = strcmpi(key,keys{2});
                data{N_Block}(t,4) = Trials{4,t};
                

                %with target
            elseif Trials{4,t} == 1
                Dis_num1 = [1:(N_Symbols/2)];
                Dis_num2 = [(N_Symbols/2)+1:(N_Symbols-1)];
                target_num3 = [N_Symbols];
                
                %1st color distractors
                text (Trials{1,t}(Dis_num1),Trials{2,t}(Dis_num1),Distractor,'Color',Set_Color);
                
                %2nd color distractors
                text (Trials{1,t}(Dis_num2),Trials{2,t}(Dis_num2),Target,'Color',dis_Color);
                
                %target
                text (Trials{1,t}(target_num3),Trials{2,t}(target_num3),Target,'Color',Set_Color);
                
                tic;
                if run_mode == 0                                      %Human mode
                    pause(); key = get(h,'CurrentCharacter');
                    
                elseif run_mode == 1                                  %Robot mode
                    time_sti = time_frame*rand();                     %random time according to settings.
                    key = keys(randi(length(keys)));                  %random key.
                    pause(time_sti);
                end
                
                %store data
                data{N_Block}(t,1) = toc;
                data{N_Block}(t,2) = strcmpi(key,keys{1});
                data{N_Block}(t,4) = Trials{4,t};
            end
            clf;
            axis off
        end
    end
end

% Finish Screen
text(0.5,0.5,['Thank you!','Please press spacebar to close'],'HorizontalAlignment','center')
Press2Continue(h,' ',run_mode);


close

save('data.mat','data')


%% Filter The Data
filter_data = data;

%removes mistakes

for i = 1:N_Block
    idx = find(~(filter_data{i}(:,2)));
    filter_data{i}(idx,:) = [];
end

Valid_trials = zeros(1,N_Block);               %Number of valid trails per block will be interested here.

%time removal
for i = 1:N_Block
    idx = find((filter_data{i}(:,1))>T_Max_filter);
    filter_data{i}(idx,:) = [];
    
    Valid_trials(i) = length(filter_data{i});
    
    if Valid_trials(i)<Min_Trials_per_Block
        error('not enough trails for analysis');
    else
        continue
    end
end

%% Plots & Analysis

%creates memory space for the 2 plots to come. 1-mean RT, 2-STD, 3-set_size, 4-
%feature/conjunction ('1'/'2')
plot1 = zeros(4,(N_Block));                                         %no target
plot2 = zeros(4,(N_Block));                                         %with target

    
for i = 1:N_Block
    current_Block = filter_data{i};
    idx = find(~(current_Block(:,4)));                              %choosing zeros - no target
    plot1(1,i) = mean(current_Block(idx,1));                        
    plot1(2,i)= std((current_Block(idx,1)));
    plot1(3,i)= current_Block(1,3);
    plot1(4,i)= current_Block(1,5);
    
    idx = find((current_Block(:,4)));                               %choosing ones - with target
    plot2(1,i) = mean(current_Block(idx,1));
    plot2(2,i)= std((current_Block(idx,1)));
    plot2(3,i)= current_Block(1,3);
    plot2(4,i)= current_Block(1,5);
end

%splitting plot1 & plot2 by feature or conjunction into 4 matrix. each one represent 1 graph.
%F - feature,C - conjunction, A - absence, P - present

plot_F_A = (plot1(:,plot1(4,:) ==1))';                  
plot_F_A = sortrows(plot_F_A,3);                            %sorting the rows by setsize for the graph

plot_C_A = (plot1(:,plot1(4,:) ==2))';
plot_C_A = sortrows(plot_C_A,3);

plot_F_P = (plot2(:,plot2(4,:) ==1))';
plot_F_P = sortrows(plot_F_P,3);

plot_C_P = (plot2(:,plot2(4,:) ==2))';
plot_C_P = sortrows(plot_C_P,3);


%creating plots
k = figure;
k.Units = 'normalized';
k.Position = [0 0.2 1 0.7];

%creating memory. each row represents one graph(FA,CA,FP,CP), the 1st column is
%the slope and the 2nd is the coefficient.
Slopes = zeros(4,2);


%Absence of target plot
x=Set_Size;                                        %same x axes for all plots

%feature
subplot(1,2,1)
hold on


RT_F_A=plot_F_A(:,1);                              %mean reaction time for feature with absence of target
std_RT_F_A = plot_F_A(:,2);                        %standard deviation for the previous means.
poly_F_A = polyfit(x',RT_F_A,1);                   %best fit for 1st degree polynom.
Slopes(1,:) = poly_F_A;                            %saves the fit on slopes
poly_F_A_fit = polyval(poly_F_A,x);                %computes the values for the fitted function.

graph_F_A = plot(x,RT_F_A, 'color', 'b');          %graph for feature absence
errorbar(x,RT_F_A,std_RT_F_A, 'color', 'b')        %error bar for previous means.
graph_F_A_fits = plot(x,poly_F_A_fit,'b--');       %graph for feature absence fit


%conjunction

RT_C_A=plot_C_A(:,1);
std_RT_C_A = plot_C_A(:,2);
poly_C_A = polyfit(x',RT_C_A,1);
Slopes(2,:) = poly_C_A;
poly_C_A_fit = polyval(poly_C_A,x);

graph_C_A = plot(x, RT_C_A,'color', 'r');
errorbar(x,RT_C_A,std_RT_C_A,'color', 'r');
graph_C_A_fits = plot(x,poly_C_A_fit,'r--');

%plot properties
title('Correctly Identified The absence Of Target')

xlim([3.7 16.3]);
ylim([0 T_Max_filter]);

xlabel('Set Size');
ylabel('Mean Reaction Time (In Sec)');

legend([graph_F_A graph_F_A_fits graph_C_A graph_C_A_fits],...
    {'Feature Search','Feature Search fitted curve','Conjunction Search',...
    'Conjunction Search fitted curve'},'FontSize', 12, 'location', 'southeast');



%present of target plot
subplot(1,2,2)
hold on

%feature
RT_F_P=plot_F_P(:,1);
std_RT_F_P = plot_F_P(:,2);
poly_F_P = polyfit(x',RT_F_P,1);
Slopes(3,:) = poly_F_P;
poly_F_P_fit = polyval(poly_F_P,x);

graph_F_P = plot(x,RT_F_P, 'color', 'b');
errorbar(x,RT_F_P,std_RT_F_P, 'color', 'b');
graph_F_P_fit = plot(x,poly_F_P_fit,'b--');



%conjunction
RT_C_P=plot_C_P(:,1);
std_RT_C_P = plot_C_P(:,2);
poly_C_P = polyfit(x',RT_C_P,1);
Slopes(4,:) = poly_C_P;
poly_C_P_fit = polyval(poly_C_P,x);

graph_C_P = plot(x, RT_C_P,'color', 'r');
errorbar(x,RT_C_P,std_RT_C_P,'color', 'r');
graph_C_P_fit = plot(x,poly_C_P_fit,'r--');


%plot properties
title('Correctly Identified The Presence Of Target');

xlim([3.7 16.3]);
ylim([0 T_Max_filter]);

xlabel('Set Size');
ylabel('Mean Reaction Time (In Sec)');


sgtitle('Reaction Times In a Visual Search Experiment');
hold off


%% correlations checks

%creating memory for the table
correlation_sig = zeros(3,4);

%compute the r and the p-value for each condition
[correlation_sig(1,:), correlation_sig(2,:)] = corr(x',([RT_F_A,RT_C_A,RT_F_P,RT_C_P]));

%converts the correlation matrix into a nice table with: R, P-value and
%significance level for each condition of the experiment.
correlation_sig = num2cell(correlation_sig);
sig_p_value = [0.05; 0.01; 0.001];
for i = 1:length(correlation_sig)
    if correlation_sig{2,i}<sig_p_value(1) && correlation_sig{2,i}>sig_p_value(2)
        correlation_sig{3,i} = '*';
    elseif correlation_sig{2,i}<sig_p_value(2) && correlation_sig{2,i}>sig_p_value(3)
        correlation_sig{3,i} = '**';
    elseif correlation_sig{2,i}<sig_p_value(3)
        correlation_sig{3,i} = '***';
    elseif correlation_sig{2,i}>sig_p_value(1)
        correlation_sig{3,i} = 'n.s';
    end
end

sig_table = array2table(correlation_sig...
    ,'VariableNames',{'Feature Absence' 'Conjunction Absence' 'Feature Present' 'Conjunction Present'});
sig_table.Properties.RowNames =  {'R','P-Value', 'significance level'};

Slopes_table = array2table(Slopes'...
    ,'VariableNames',{'Feature Absence' 'Conjunction Absence' 'Feature Present' 'Conjunction Present'});
Slopes_table.Properties.RowNames =  {'Slope','Coefficient'};
