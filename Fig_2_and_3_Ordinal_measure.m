T = readtable('temp_data_full_batches.csv', 'PreserveVariableNames',true);
data = table2array(T(:, [24:28 32:36 40:44 48:52 56:60 64:68])).';
data(isnan(data)) = 0;
opt = 1:30;
opt_prefs = repmat(opt', 1, 484);
reshaped = Inf(30, 484);

batch_incorp = 1; % Set this equal to the max no. batches to incorporate in analysis

B = 5 * batch_incorp;

%% Data Reshaping

% The first for loop reassigns values for courses based on where in the
% preference list they are selected. For example, if a student's preference
% list is [2, 10, 15, 20, 25] in the original configuration, then in the
% reshaped data the vector element 2 is given a value of 1 to show that it 
% is chosen as the first preference, vector element 10 is given a value of
% 2 since it is the second preference and so on. All unselected preferences
% are defaulted to infinity (Inf) to signify that they are unchosen, or
% beyond the scope of the preference list. As such, each student's
% preference list is a 30-by-1 vector containing up to 5 real numbers for
% the baseline mechanism and up to 30 real numbers for the batched
% mechanism.

for p = 1:484
    for i = 1:B
        if data(i, p) > 0
            reshaped((data(i, p)), p) = i;
        end
    end
end

% The second for loop sets the optimal preferences depending on the
% number of batches a student undergoes in the batched mechanism
% treatments. The optimal vector into which we want to change the student's
% preference list is simply an increasing list of 1, 2, 3.... Thus we do not
% necessarily need to use the following loop to compare how different the
% preference lists are from the optimal.

data_batch_specific = [data([1:B], :);  zeros(30 - B, 484)];

for p = 1:484
    for i = 1:30
        if data_batch_specific(i, p) == 0
            opt_prefs(i, p) = Inf;
        end
    end
end

original_mat = reshaped; % Preserves the original reshaped data for comparison.

%% Kendall tau distance (counting discordant pairs à la https://en.wikipedia.org/wiki/Kendall_tau_distance)

% This loop counts the number of pairs of vector elements in the student's 
% preferences whose order is opposite to the order of pairs in the optimal 
% preferences. For example, if vector element (or course) 2 is ranked as 5
% but vector element (or course) 3 is ranked as 4, then this is a
% discordant pair because vector element 2 should be ranked lower than
% vector element 3 as per the optimal preferences. The loop then aggregates
% the amount of discordant pairs for each vector element, not just adjacent
% pairs. This loop compares element 1 against element 2, 3, 4... 30. This is
% equivalent to counting the number of interchanges of discordant adjacent pairs.

discordant_pairs = zeros(30, 484);

for p = 1:484
    for i = 1:29
        for k = 1:30-i
            if reshaped(i, p) > reshaped(i + k, p)
                discordant_pairs(i, p) = 1 + discordant_pairs(i, p);                  
            end
        end
    end
end

%% Scale factor

% The scale factor is the number of preferences that a student chose.
% Technically speaking, it is the number of non-Inf elements.

sf = sum(reshaped < Inf);

dp = sum(discordant_pairs);
denom = (0.5*30*(30-1))-(0.5*(30-sf).*(30-sf-1));

%% Kendall tau distances from optimal preferences (scaled)

Dist_opt = (dp./denom);

T1 = Dist_opt(1:121); % Treatment 1: Baseline
T2 = Dist_opt(122:242); % Treatment 2: Batches
T3 = Dist_opt(243:363); % Treatment 3: Baseline with advice
T4 = Dist_opt(364:484); % Treatment 4: Batches with advice

combine = [Dist_opt ; data]';
order = sortrows(combine , 1);

% %% Full sample figures
% 
% figure
% sgtitle("Full Sample ECDFs of Kendall Distances (incl. Batches 1-6)")
% 
% subplot(2,2,1)
%     ecdf(T1);
%     hold on
%     ecdf(T2);
%     hold off
%     legend(["Baseline", "Batched"], 'Location', 'southeast')
%     title('T1 vs. T2')
%     
% subplot(2,2,2)
%     ecdf(T1);
%     hold on
%     ecdf(T3);
%     hold off
%     legend(["Baseline", "Baseline (Advice)"], 'Location', 'southeast')
%     title('T1 vs. T3')
% 
%     
% subplot(2,2,3)
%     ecdf(T2);
%     hold on
%     ecdf(T4);
%     hold off
%     legend(["Batched", "Batched (Advice)"], 'Location', 'southeast')
%     title('T2 vs. T4')
%     
% subplot(2,2,4)
%     ecdf(T3);
%     hold on
%     ecdf(T4);
%     hold off
%     legend(["Baseline (Advice)", "Batched (Advice)"], 'Location', 'southeast')
%     title('T3 vs. T4')
%             
% figure
%     ecdf(T1);
%     hold on
%     ecdf(T2);
%     ecdf(T3);
% 	ecdf(T4);
%     hold off 
%     legend(["Baseline", "Batched", "Baseline (Advice)", "Batched (Advice)"], 'Location', 'southeast')
%     title('Full Sample ECDFs of Kendall Distances (incl. Batches 1-6)')

%% GEB and Non-GEB sub-sample splitting

geb_base = sum(data([1:4], [1:121 243:363]) == 26); % 
geb_batc = sum(data([1:B-1], [122:242 364:484]) == 26);
geb_combine = [geb_base(1:121) geb_batc(1:121) geb_base(122:242) geb_batc(122:242)];

att_ch = table2array(T(:, [19:22])).';
inatt_1 = att_ch(1, :) ~= 1;
inatt_2 = att_ch(2, :) ~= 2;
inatt_3 = att_ch(3, :) ~= 3;
inatt_4 = att_ch(4, :) > 2;
inatt_comb = [inatt_1; inatt_2; inatt_3; inatt_4]; % Acceptable attention response is either [1, 2, 3, 1] or [1, 2, 3, 2]
inatt_vec = sum(inatt_comb);

time_taken = table2array(T(:, 7)).';
timeframe = time_taken < 1500 & time_taken > 120; % Between 2 and 25 mins
timeframe = 1 - timeframe;

adv_fol = table2array(T(:, 5)).';
adv_fol(243:484) = 1 - adv_fol(243:484);

controls = [geb_combine; inatt_vec; timeframe; adv_fol];

non_geb = geb_combine;
rat_con = sum(controls);
bias_transpose = controls';

bias_prop = (sum(bias_transpose)/484); % NB The Baseline treatments do not need to follow advice! Last element should be x/242

geb_t = [sum(geb_combine(1:121))/121; sum(geb_combine(122:242))/121; sum(geb_combine(243:363))/121; sum(geb_combine(364:484))/121];
geb_treatments = geb_t';

disp(['Proportion of GEB/Inattention/Timeframe/Advice Followed (T1-T4): ', num2str(bias_prop)]);
disp(['Proportion of GEB per T1-T4: ', num2str(geb_treatments)]);


cont = find(rat_con == 0);

cont_T1 = Dist_opt(cont(cont >= 1 & cont <= 121));
cont_T2 = Dist_opt(cont(cont >= 122 & cont <= 242));
cont_T3 = Dist_opt(cont(cont >= 243 & cont <= 363));
cont_T4 = Dist_opt(cont(cont >= 364 & cont <= 484));

ngeb = find(non_geb == 0);

ngeb_T1 = Dist_opt(ngeb(ngeb >= 1 & ngeb <= 121));
ngeb_T2 = Dist_opt(ngeb(ngeb >= 122 & ngeb <= 242));
ngeb_T3 = Dist_opt(ngeb(ngeb >= 243 & ngeb <= 363));
ngeb_T4 = Dist_opt(ngeb(ngeb >= 364 & ngeb <= 484));


% figure
% sgtitle("Controls Sub-Sample ECDF of Kendall Distances (incl. Batches 1-6)")
% 
% subplot(2,2,1)
%     ecdf(cont_T1);
%     hold on
%     ecdf(cont_T2);
%     hold off
%     legend(["Baseline", "Batched"], 'Location', 'southeast')
%     title('T1 vs. T2')
%     
% subplot(2,2,2)
%     ecdf(cont_T1);
%     hold on
%     ecdf(cont_T3);
%     hold off
%     legend(["Baseline", "Baseline (Advice)"], 'Location', 'southeast')
%     title('T1 vs. T3')
%     
% subplot(2,2,3)
%     ecdf(cont_T2);
%     hold on
%     ecdf(cont_T4);
%     hold off
%     legend(["Batched", "Batched (Advice)"], 'Location', 'southeast')
%     title('T2 vs. T4')
%     
% subplot(2,2,4)
%     ecdf(cont_T3);
%     hold on
%     ecdf(cont_T4);
%     hold off
%     legend(["Baseline (Advice)", "Batched (Advice)"], 'Location', 'southeast')
%     title('T3 vs. T4')
% 
% figure
%     ecdf(cont_T1);
%     hold on
%     ecdf(cont_T2);
%     ecdf(cont_T3);
% 	ecdf(cont_T4);
%     hold off 
%     legend(["Baseline", "Batched", "Baseline (Advice)", "Batched (Advice)"], 'Location', 'southeast')
%     title('Controls Sub-Sample ECDFs of Kendall Distances (incl. Batches 1-6)')

    %% FSD Tests - FULL SAMPLE

dom=0:0.1:1;

T1=sort(T1);
T2=sort(T2);
T3=sort(T3);
T4=sort(T4);

FSD_T1=zeros(1,length(dom));
FSD_T2=zeros(1,length(dom));
FSD_T3=zeros(1,length(dom));
FSD_T4=zeros(1,length(dom));

for i=1:length(dom)
   FSD_T1(i)= mean(T1<=dom(i));
   FSD_T2(i)= mean(T2<=dom(i));
   FSD_T3(i)= mean(T3<=dom(i));
   FSD_T4(i)= mean(T4<=dom(i));

end

figure
sgtitle("Full Sample (incl. Batches 1-6)")

subplot(1,4,1)
    plot(dom,FSD_T1)
    hold on
    plot(dom,FSD_T2)
    hold off
    legend(["Baseline", "Batched"], 'Location', 'southeast')
    title('BSL vs. BAT')
    
subplot(1,4,2)
    plot(dom,FSD_T1)
    hold on
    plot(dom,FSD_T3)
    hold off
    legend(["Baseline", "Baseline (Advice)"], 'Location', 'southeast')
    title('BSL vs. BSL-ADV')
    
subplot(1,4,3)
    plot(dom,FSD_T2)
    hold on
    plot(dom,FSD_T4)
    hold off
    legend(["Batched", "Batched (Advice)"], 'Location', 'southeast')
    title('BAT vs. BAT-ADV')
    
subplot(1,4,4)
    plot(dom,FSD_T3);
    hold on
    plot(dom,FSD_T4)
    hold off
    legend(["Baseline (Advice)", "Batched (Advice)"], 'Location', 'southeast')
    title('BSL-ADV vs. BAT-ADV')

%% FSD Tests - CONTROLS

dom=0:0.1:1;

cont_T1=sort(cont_T1);
cont_T2=sort(cont_T2);
cont_T3=sort(cont_T3);
cont_T4=sort(cont_T4);

FSD_T1=zeros(1,length(dom));
FSD_T2=zeros(1,length(dom));
FSD_T3=zeros(1,length(dom));
FSD_T4=zeros(1,length(dom));

for i=1:length(dom)
   FSD_T1(i)= mean(cont_T1<=dom(i));
   FSD_T2(i)= mean(cont_T2<=dom(i));
   FSD_T3(i)= mean(cont_T3<=dom(i));
   FSD_T4(i)= mean(cont_T4<=dom(i));

end

figure
sgtitle("Controlled NGEB Sample (incl. Batches 1-6)")

subplot(1,4,1)
    plot(dom,FSD_T1)
    hold on
    plot(dom,FSD_T2)
    hold off
    legend(["Baseline", "Batched"], 'Location', 'southeast')
    title('BSL vs. BAT')
    
subplot(1,4,2)
    plot(dom,FSD_T1)
    hold on
    plot(dom,FSD_T3)
    hold off
    legend(["Baseline", "Baseline (Advice)"], 'Location', 'southeast')
    title('BSL vs. BSL-ADV')
    
subplot(1,4,3)
    plot(dom,FSD_T2)
    hold on
    plot(dom,FSD_T4)
    hold off
    legend(["Batched", "Batched (Advice)"], 'Location', 'southeast')
    title('BAT vs. BAT-ADV')
    
subplot(1,4,4)
    plot(dom,FSD_T3);
    hold on
    plot(dom,FSD_T4)
    hold off
    legend(["Baseline (Advice)", "Batched (Advice)"], 'Location', 'southeast')
    title('BSL-ADV vs. BAT-ADV')

%% FSD Tests - NON-GEB

dom=0:0.1:1;

ngeb_T1=sort(ngeb_T1);
ngeb_T2=sort(ngeb_T2);
ngeb_T3=sort(ngeb_T3);
ngeb_T4=sort(ngeb_T4);

FSD_T1=zeros(1,length(dom));
FSD_T2=zeros(1,length(dom));
FSD_T3=zeros(1,length(dom));
FSD_T4=zeros(1,length(dom));

for i=1:length(dom)
   FSD_T1(i)= mean(ngeb_T1<=dom(i));
   FSD_T2(i)= mean(ngeb_T2<=dom(i));
   FSD_T3(i)= mean(ngeb_T3<=dom(i));
   FSD_T4(i)= mean(ngeb_T4<=dom(i));

end

figure
sgtitle("NGEB Sample (incl. Batches 1-6)")

subplot(1,4,1)
    plot(dom,FSD_T1)
    hold on
    plot(dom,FSD_T2)
    hold off
    legend(["Baseline", "Batched"], 'Location', 'southeast')
    title('BSL vs. BAT')
    
subplot(1,4,2)
    plot(dom,FSD_T1)
    hold on
    plot(dom,FSD_T3)
    hold off
    legend(["Baseline", "Baseline (Advice)"], 'Location', 'southeast')
    title('BSL vs. BSL-ADV')
    
subplot(1,4,3)
    plot(dom,FSD_T2)
    hold on
    plot(dom,FSD_T4)
    hold off
    legend(["Batched", "Batched (Advice)"], 'Location', 'southeast')
    title('BAT vs. BAT-ADV')
    
subplot(1,4,4)
    plot(dom,FSD_T3);
    hold on
    plot(dom,FSD_T4)
    hold off
    legend(["Baseline (Advice)", "Batched (Advice)"], 'Location', 'southeast')
    title('BSL-ADV vs. BAT-ADV')