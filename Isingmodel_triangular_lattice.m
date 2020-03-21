function final_ising(L,conc_Mn,repeat)
clc
close all
disp('....MATLAB CODE FOR THE CALCULATION OF MAGNETIC PROPERTIES (MAGNETIZATION AND CRITICAL TEMPERATURE) OF DILUTE MAGNETIC SEMICONDUCTOR.......')
disp('Interaction terms are taken from: Ramasubramaniam, Ashwin, and Doron Naveh. Physical Review B 87.19 (2013): 195201.')
disp('.... Using 2D Ising Lattice Gas Model for the Triangular lattice ....')
disp('********** by Akash Jain ***********************')
%%
% f1_shell() will find Mn atoms in the 1st shell
    function  [count,first_shell_nn,nn_Mn] = f1_shell(site,S,L)
        % This function will find the 1st n.n. shell for any given site
        %
        i = site(1,1);
        j = site(1,2);
        count = 0; % number of Mn in the 1st nearest neighbor shell
        % four n.n are just like square lattice
        % (i+1,j), (i-1,j), (i,j+1), (i,j-1)
        % other 2 n.n. depends on the row
        if mod(i,2)==0 % for even numbered row
            z = 1; % n.n. in upper Right row and down Right row
        else % for ODD numbered row
            z = -1; % n.n. in upper left row and down left row
        end
        % Six nearest neighbors in hexagonal lattice
        nn = [i+1,j;i-1,j;i,j+1;i,j-1;i+z,j-1;i+z,j+1];
        first_shell_nn  = zeros(6,2);
        %
        for r = 1:1:6 % only 6 n.n. in 1st shell
            %
            m = nn(r,1);
            if m>L % PBCs
                m = m - L;
            elseif m<1
                m = m + L;
            end
            %
            k = nn(r,2);
            if k>L  % PBCs
                k = k - L;
            elseif k<=0
                k = k + L  ;
            end
            %
            first_shell_nn(r,:) = [m,k];
            if S(m,k)~=0  % Check for Mn in first shell
                % if there is a Mn in first shell, then make a list
                if count == 0
                    nn_Mn = [m,k];
                    count = count + 1;
                else
                    nn_Mn = [nn_Mn;[m,k]];
                    count = count + 1;
                end
            end
        end
        %
        if count == 0    % no Mn found
            nn_Mn=[];
        end
    end
%%
% f2_shell() will find Mn atoms in 2nd n.n shell
    function  [count,sec_next_shell,nn_2nd_Mn] = f2_shell(site,first_next_shell_nn,S,L)
        % count = total number of Mn in i+1_th shell
        % sec_next_shell = coordinates of i+1_th shell elements
        % nn_2nd_Mn, coordinates of all Mn atoms in the i+1_th shell
        n_nn = length(first_next_shell_nn);
        count = 0;
        flag = 0;
        % make list of  n.n. of the atoms in the 1st shell
        for i = 1:1:n_nn
            %
            [numberofMn,first_shell,nn_Mn_cord] = f1_shell(first_next_shell_nn(i,:),S,L);
            % make complete list to access next shell
            if flag == 0
                sec_next_shell = first_shell;
                flag = 1;
            else
                sec_next_shell = [sec_next_shell;first_shell];
            end
        end
        % remove duplicate atoms in the list
        sec_next_shell = unique(sec_next_shell,'rows') ;
        site2 = [site;first_next_shell_nn];
        % remove atoms present in the inner shells
        Lia2 = ismember(sec_next_shell,site2,'rows');
        sec_next_shell = sec_next_shell(~Lia2,:);
        nn2ndshell = length(sec_next_shell);
        % Find Mn in the 2nd Next shell
        for i = 1:1:nn2ndshell
            m = sec_next_shell(i,1);
            k = sec_next_shell(i,2);
            if S(m,k)~=0  % if there is a Mn in the shell, then make a list
                if count == 0
                    nn_2nd_Mn = [m,k];
                    count = count + 1;
                else
                    nn_2nd_Mn = [nn_2nd_Mn;[m,k]];
                    count = count + 1;
                end
            end
        end
        %
        if count == 0
            % no Mn found
            nn_2nd_Mn=[];
        end
    end
%%
% final_ising() to calculate Tc
sprintf(' %d x %d triangular lattice',L,L)
sprintf('Sample %d random distributions of Mn in Mo top layer', repeat)
%
N = 2*10^4 ; % total number of sweeps
%
sprintf('Total sweeps = %d', N)
% Unit conversions
Kb =  8.6173*10^-5; % boltzmann constant in eV/K
J0 = 0.17; % J1st = 1; J2 =0.17; J3rd = 0.11;J4th = 0.055 = 0.06; PBE
% J0 = 0.22 eV; J1st = 1; J2 =0.14; J3rd = 0.09;J4th = 0.055 = 0.05; HSE
T0 = J0/Kb;% for non-dimensional units
%
% conc_Mn = concentration of Mn in top layer
%
sprintf('%% of Mn in top layer of MoS2 = %0.1f ',100*conc_Mn)
%
m = 10; % bins
%
Temp = [0.05:0.01:0.70];
%
n_Mn = conc_Mn*L^2;
%
sprintf('%d Mn atoms and %d Total Mn+Mo atoms in top layer of MoS2',n_Mn,L^2)
%
Nspin = L*L;
%
% start with new random distribution of Mn atoms
average_T = 0;
max_Xt = zeros(repeat,1);
Tc =  zeros(repeat,1);
%%
% q different initial distributions of Mn to  get  average and std for Tc
for q = 1:1:repeat
    %
    % start with new random distribution of Mn atoms
    S = zeros(1,L^2); % initially, let's take all Mo atoms in top layer
    %
    for i = randperm(L^2,n_Mn)
        S(i) = 1; % 1-2*round(rand);
    end
    %
    S = reshape(S,[L,L]); % convert to square matrix
    Mn_coord = zeros(n_Mn,2);
    % figure
    % imshow(S,'InitialMagnification',1000);
    %     h1 = gcf;
    %     refresh(h1);
    %     pause(0.1);
    countMn = 0;
    %
    % store Mn coordinates
    %
    for i = 1:1:L
        for j = 1:1:L
            if S(i,j) ~= 0% Mo has spin = 0
                countMn = countMn + 1;
                Mn_coord(countMn,:) = [i,j];
            end
        end
    end
    %
    % store n.n. list for each Mn atom
    % 
    max_site_nn = 60; % maximum numbers of n.n. within 4 shells
    % store (i,j)
    all_nn = zeros(max_site_nn,2,n_Mn);
    %
    % to access nn list of say i_th Mn atom-->  all_nn(:,:,i)
    %
    number_Mn = zeros(n_Mn,3);
    total_Mn =  zeros(n_Mn,1);
    for ii = 1:1:n_Mn
        %
        site = Mn_coord(ii,:);
        %  First Shell
        [count1st,first_shell,nn_1st_Mn] = f1_shell(site,S,L);
        %  Second Shell
        [count2nd,sec_shell_nn,nn_2nd_Mn] = f2_shell(site,first_shell,S,L);
        %  Third Shell
        [count3rd,third_shell_nn,nn_3rd_Mn] = f2_shell([site;first_shell],sec_shell_nn,S,L);
        % Fourth Shell
        [count4th,four_shell_nn,nn_4th_Mn] = f2_shell([site;first_shell;sec_shell_nn],third_shell_nn,S,L);
        % Fifth Shell
        % [count5th,five_shell_nn,nn_5th_Mn] = f2_shell([site;first_shell;sec_shell_nn;third_shell_nn],four_shell_nn,S,L)
        %
        number_Mn(ii,:) = [count1st,count1st+count2nd,count1st+count2nd+count3rd];
        total_Mn(ii) =  count1st + count2nd + count3rd + count4th;
        if total_Mn(ii) > 0
            all_nn(1:total_Mn(ii),:,ii) = [nn_1st_Mn;nn_2nd_Mn;nn_3rd_Mn;nn_4th_Mn];
        end
    end
    %
    % Starting metropolis monte-carlo method
    %
    dE = 0;
    % assume initial energy = 0;
    energy = 0;
    %
    MCegy = zeros(1,length(Temp));
    MCegy2 = MCegy;
    MC_Mag = MCegy;
    MC_Mag2 = MCegy;
    MC_sp_heat = MCegy;
    MC_sp_heat2 = MCegy;
    MC_susept = MCegy;
    MC_susept2 = MCegy;
    %
    % Start Monte-carlo
    for T = 1:1:length(Temp)
        all_avg_egy = 0;
        all_avg_egy2 = 0;
        all_avg_mag = 0;
        all_avg_mag2 = 0;
        all_sp_heat = 0;
        all_sp_heat2 = 0;
        all_susept = 0;
        all_susept2 = 0;
        %
        for bin = 1:1:m % sample Tc from in 'm' bins
            %
            sweeps = 1;
            count  = 0;
            bin_avg_egy = 0;
            bin_avg_egy2 = 0;
            bin_avg_mag = 0;
            bin_avg_mag2 = 0;
            bin_sp_heat = 0;
            bin_susept = 0;
            avg_egy = 0;
            avg_egy2 = 0;
            Magnetization = 0;
            Magnetization2 = 0;
            accept_less = 0;
            accept_p = 0;
            accept_r = 0;
            %
            while(sweeps<=N) % N sweeps
                %
                % select any Mn atom randomly
                ii = randi(n_Mn);
                site = Mn_coord(ii,:);
                nnofMn = total_Mn(ii) ;
                %
                trial_dE = 0;
                for nn = 1:1:nnofMn
                    %
                    nm  = all_nn(nn,1,ii);
                    nk  = all_nn(nn,2,ii);
                    %
                    % find trial_dE to decide whether to flip spin or not
                    if nm ~=0 && nk ~= 0 % don't calculate trial_dE, if there is no Mn in the 1-4 n.n. shells
                        %
                        if nn <= number_Mn(ii,1) % first shell
                            %
                            trial_dE = trial_dE + 1*2*S(site(1,1),site(1,2))*S(nm,nk); % PBE,HSE
                            %
                        elseif nn > number_Mn(ii,1) &&  nn <= number_Mn(ii,2) % 2nd shell
                            %
                            trial_dE = trial_dE + 0.18*2*S(site(1,1),site(1,2))*S(nm,nk); % PBE
                            %                 trial_dE = trial_dE +2*0.14*S(site(1,1),site(1,2))*S(nm,nk); % HSE
                            %
                        elseif nn > number_Mn(ii,2) &&  nn <= number_Mn(ii,3) % 3rd shell
                            %
                            trial_dE = trial_dE + 2*0.12*S(site(1,1),site(1,2))*S(nm,nk); % PBE
                            %                 trial_dE = trial_dE + 2*0.09*S(site(1,1),site(1,2))*S(nm,nk); %HSE
                            %
                        else
                            % trial_dE is non-dimensional throughout
                            trial_dE = trial_dE + 2*0.06*S(site(1,1),site(1,2))*S(nm,nk);% PBE
                            %                 trial_dE = trial_dE + 2*0.05*S(site(1,1),site(1,2))*S(nm,nk);% HSE
                            %
                        end
                        %
                    end
                    %
                end
                %
                % Now, accept or reject flip/move
                dE =0;
                %
                if trial_dE < 0 % Accept low energy
                    %
                    S(site(1,1),site(1,2)) = -S(site(1,1),site(1,2));
                    dE = trial_dE; % flip spin
                    sweeps = sweeps + 1;
                    accept_less = accept_less+1;
                elseif trial_dE > 0
                    r = rand;
                    if exp(-trial_dE/(Temp(T)))> r
                        S(site(1,1),site(1,2)) = -S(site(1,1),site(1,2));
                        dE = trial_dE; % flip spin
                        sweeps = sweeps + 1;
                        accept_p = accept_p +1;
                    else
                        sweeps = sweeps + 1;% don't  flip
                        dE = 0;
                        accept_r = accept_r +1;
                    end
                end
                %
                temp = energy + dE;
                energy = temp;
                energy2 = energy^2;
                %
                % specific magnetization
                %
                mag_b = abs(sum(sum(S)))/n_Mn;
                mag_b2 = mag_b*mag_b;
                %
                % Accumulate data for a bin, n_Mn
                %
                if (sweeps>100) && (mod(sweeps,n_Mn^2)==1)
                    Magnetization =  Magnetization + mag_b;
                    Magnetization2 = Magnetization2 + mag_b2 ;
                    avg_egy = avg_egy + energy;
                    avg_egy2 = avg_egy2 + energy*energy;
                    count=count+1 ;
                end
            end
            %     sprintf('T = %0.02f, bin = %d, moves accepted with dE < 0 = %0.2f',Temp(T),bin,100*accept_less/N)
            %     sprintf('T = %0.02f, bin = %d, moves accepted with dE > 0 = %0.2f',Temp(T),bin,100*accept_p/N)
            %     sprintf('T = %0.02f, bin = %d, moves rejected with dE > 0 = %0.2f',Temp(T),bin,100*accept_r/N)
            %
            % save E,E^2,M,M^2 after energy minimization
            bin_avg_egy = avg_egy/count;
            bin_avg_egy2 = avg_egy2/count;
            bin_avg_mag = Magnetization/count;
            bin_avg_mag2 = Magnetization2/count;
            bin_sp_heat = (bin_avg_egy2-bin_avg_egy^2)/(n_Mn*(Temp(T)^2));
            bin_susept = n_Mn*n_Mn*(bin_avg_mag2-bin_avg_mag^2)/(Temp(T));
            %
            %  sum averages for the overall averaging
            %
            all_avg_egy = all_avg_egy +  bin_avg_egy ;
            all_avg_egy2 = all_avg_egy2 + bin_avg_egy^2;
            all_avg_mag = all_avg_mag + bin_avg_mag;
            all_avg_mag2 = all_avg_mag2 + bin_avg_mag^2;
            all_sp_heat = all_sp_heat + bin_sp_heat;
            all_sp_heat2 = all_sp_heat2 + bin_sp_heat^2;
            all_susept = all_susept + bin_susept;
            all_susept2 =  all_susept2 + bin_susept^2;
        end % for bin loop end here
        %
        % average over 'm' bins 
        %
        MCegy(1,T) = all_avg_egy/m;
        MCegy2(1,T) = all_avg_egy2/m;
        MC_Mag(1,T) = all_avg_mag/m;
        MC_Mag2(1,T) = all_avg_mag2/m;
        MC_sp_heat(1,T) = all_sp_heat/m;
        MC_sp_heat2(1,T) = all_sp_heat2/m;
        MC_susept(1,T) = all_susept/m;
        MC_susept2(1,T) = all_susept2/m;
        %
    end % for 'T' end here 
    uncertainity_egy = sqrt(MCegy2-MCegy.^2);
    un_magnetization = sqrt(MC_Mag2-MC_Mag.^2);
    un_specific_heat = sqrt(MC_sp_heat2-MC_sp_heat.^2);
    un_susceptibility = sqrt(MC_susept2-MC_susept.^2);
    %
    [max_Xt(q),index] =  max(MC_susept);
    Tc(q) = Temp(index);
end % end of 'q' repetition
average_T = mean(Tc);
uncertainity_T = std(Tc);
sprintf('Curie temperature, T_c  = %0.0f K, error = %0.2f K',average_T*T0,uncertainity_T*T0)
% end of overall function
end