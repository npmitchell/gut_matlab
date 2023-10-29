% findHCRcombinations.m 
%
% To find configurations of probes and amplifiers that avoid adjacent 
% probes having adjacent wavelengths in the amplifiers, we map the 
% problem into a graph theory problem.
% The probes graph shows which probes are physically adjacent or 
% overlapping.
% The amplifier graph shows which wavelengths are troublesome to have 
% nearby. 
% We then find all possible mappings from graph B nodes to Graph A nodes 
% such that connected nodes in Graph A are not connected in Graph B. 
% 
% Noah Mitchell 2023


%% Clear the workspace
clc; close all; clearvars

%% Define options 
filterHoxBasedOnProbeAvailability = true ;
ensureUniqueProbes = true ;

%% Create adjacency matrices for Probes and Amplifiers

% Adjacency/overlap matrix for probes. Connected nodes in the hox gene 
% graph are considered adjacent/overlapping. 
%            
%                   Lab     AbdA
%                  .-o-------o
%        Scr  Antp/   \\_  /  \
%          o-----o     | \.__  \ 
%                 \     \/   \__\ 
%                  .-----o-------o
%                       Ubx      AbdB
%

% % Default config ==> no solutions without adding amplifiers
% %          Scr Antp Lab Ubx AbdA AbdB             
probes =  [  0   1   0   0   0   0;       % Scr
             1   0   1   1   0   0;       % Antp
             0   1   0   1   1   1;       % Lab
             0   1   1   0   1   1;       % Ubx
             0   0   1   1   0   1;       % ABDA
             0   0   1   1   1   0];      % ABDB
pnames = {'Scr', 'Antp', 'Lab', 'Ubx', 'AbdA', 'AbdB'} ;

%% TESTING OTHER ADJACENCY CONSIDERATIONS IN AN ATTEMPT TO FIND VALID SOLUTIONS WITHOUT ADDITIONAL AMPLIFIERS
% % Sever Scr/Antp ==> no solutions without adding amplifiers
% %          Scr Antp Lab Ubx AbdA AbdB      
% probes =  [  0   0   0   0   0   0;       % Scr
%              0   0   1   1   0   0;       % Antp
%              0   1   0   1   1   1;       % Lab
%              0   1   1   0   1   1;       % Ubx
%              0   0   1   1   0   1;       % ABDA
%              0   0   1   1   1   0];      % ABDB
% 
% % Sever Lab/Antp ==> 2 solutions without adding amplifiers: 
% % 'Scr:647 Antp:568 Lab:546 Ubx:488 AbdA:405 AbdB:660 // Scr:3 Antp:7 Lab:1 Ubx:2 AbdA:6 AbdB:8'
% % 'Scr:647 Antp:568 Lab:546 Ubx:488 AbdA:594 AbdB:660 // Scr:3 Antp:7 Lab:1 Ubx:2 AbdA:5 AbdB:8'
% %          Scr Antp Lab Ubx AbdA AbdB      
% probes =  [  0   1   0   0   0   0;       % Scr
%              1   0   0   1   0   0;       % Antp
%              0   0   0   1   1   1;       % Lab
%              0   1   1   0   1   1;       % Ubx
%              0   0   1   1   0   1;       % ABDA
%              0   0   1   1   1   0];      % ABDB
% 
% % Sever Scr/Antp & Lab/Antp ==> 3 solutions without adding amplifiers: 
% % 'Scr:594 Antp:568 Lab:546 Ubx:488 AbdA:405 AbdB:660 // Scr:3 Antp:7 Lab:1 Ubx:2 AbdA:6 AbdB:8'
% % 'Scr:647 Antp:568 Lab:546 Ubx:488 AbdA:405 AbdB:660 // Scr:3 Antp:7 Lab:1 Ubx:2 AbdA:6 AbdB:8'
% % 'Scr:647 Antp:568 Lab:546 Ubx:488 AbdA:594 AbdB:660 // Scr:3 Antp:7 Lab:1 Ubx:2 AbdA:5 AbdB:8'
% %          Scr Antp Lab Ubx AbdA AbdB      
probes =  [  0   0   0   0   0   0;       % Scr
             0   0   0   1   0   0;       % Antp
             0   0   0   1   1   1;       % Lab
             0   1   1   0   1   1;       % Ubx
             0   0   1   1   0   1;       % ABDA
             0   0   1   1   1   0];      % ABDB
% 
% % Sever Lab/AbdB ==> no solutions without adding amplifiers
% %          Scr Antp Lab Ubx AbdA AbdB      
% probes =  [  0   1   0   0   0   0;       % Scr
%              1   0   1   1   0   0;       % Antp
%              0   1   0   1   1   0;       % Lab
%              0   1   1   0   1   1;       % Ubx
%              0   0   1   1   0   1;       % ABDA
%              0   0   0   1   1   0];      % ABDB
% pnames = {'Scr', 'Antp', 'Lab', 'Ubx', 'AbdA', 'AbdB'} ;


%% Adjacency matrix for Alexa Fluor amplifiers
%
%           o   o---o---o---o---o   o---o
%
%          405 488 514 546 568 594 647 660
ampfrs = [  0   0   0   0   0   0   0   0;       % 405 --> far enough away from 488 that no overlap
            0   0   1   0   0   0   0   0;       % 488
            0   1   0   1   0   0   0   0;       % 514
            0   0   1   0   1   0   0   0;       % 546
            0   0   0   1   0   1   0   0;       % 568
            0   0   0   0   1   0   0   0;       % 594 --> far enough away from 647 that no overlap in penultimate column
            0   0   0   0   0   0   0   1;       % 647 --> strong overlap with 660, could consider these the same channel
            0   0   0   0   0   0   1   0];      % 660 --> strong overlap with 647, could consider these the same channel
anames = {'405', '488', '514', '546', '568', '594', '647', '660'} ;

% Define colors to display these amplifiers (use excitation wavelength color)
colors = [ 0.00, 0.00, 1.00;  ... % 405
           0.00, 0.40, 0.75;  ... % 488
           0.00, 0.60, 0.50;  ... % 514
           0.00, 0.80, 0.25;  ... % 546
           0.25, 1.00, 0.00;  ... % 568
           0.50, 0.50, 0.00;  ... % 594
           0.75, 0.25, 0.00;  ... % 647
           1.00, 0.00, 0.00];     % 660
colors = colors ./ vecnorm(colors, 2, 2) ; % normalize the colors to unity

% Now input the constraints -- we have only certain probes in stock and validated
%            B1 B2 B3 B4 B5 B6 B7 B8              
HoxProbes = [ 0  0  1  0  0  0  0  0;   % Scr
              1  0  0  0  0  0  1  0;   % Antp
              1  0  0  0  0  0  0  0;   % Lab
              0  1  0  1  0  0  0  0;   % Ubx
              1  0  0  0  1  1  0  0;   % ABDA
              0  0  1  0  0  0  0  1];  % ABDB

%            405 488 514 546 568 594 647 660
ProbeAmps = [ 0   1   0   1   0   0   1   0  ;  % B1
              0   1   0   1   0   0   0   0  ;  % B2
              0   0   0   0   0   1   1   0  ;  % B3
              0   0   1   1   0   0   1   0  ;  % B4
              0   1   1   0   0   1   0   0  ;  % B5
              1   0   0   0   0   0   0   0  ;  % B6
              0   0   0   0   1   0   0   0  ;  % B7 --> buying AF647
              0   0   0   0   0   0   0   1  ]; % B8

nprobes = size(probes, 1) ;
nampfrs = size(ampfrs, 1) ;

% Find all permutations of the amplifiers 
permutations = perms(1:nampfrs);

% Filter unique permutations to have length nprobes
combinations = unique(permutations(:, 1:nprobes), 'rows');

% Each combination i contains combinations(i, j): the amplifier index for
% each probe j. So rows are combination numbers, columns are probe indices,
% and the elements are amplifier indices.

%% Check if each combination breaks adjacency rules: no spatially 
% adjacent/overlapping probes can have adjacent amplifiers.
goodOutput = {} ;
goodCombos = false(size(combinations, 1), 1) ;
dmyk = 1 ;
for ii = 1:size(combinations, 1)
    % combo has the indices of amplifiers that 
    combo = combinations(ii, :) ;
    disp(['----------------------new amplifier sequence: ' num2str(combo)])

    output = cell2mat(cellfun(@(p, a) [p, ':', a, ' '], pnames, anames(combo), 'UniformOutput', false));
    output = output(1:end-1);  % Remove the trailing character
    disp(output)

    isok = true ;

    % consider each probe
    probesWithSpecifiedAmp = cell(1, length(combo)) ;
    probesWithSpecifiedAmp(:) = {[]};
    for jj = 1:length(combo)
        % Stop checking if the combination is bad
        if isok
            % This probe has the amplifier given by the index probej_amp
            probej_amp = combo(jj) ;

            % Filter to see whether this amplifier is available for this hox
            % probe. This currently does not distinguish whether or not the
            % probe for the RNA is unique (ie B1 not used by another stain)
            if filterHoxBasedOnProbeAvailability || ensureUniqueProbes
                bbs = find(HoxProbes(jj, :)) ;
                amp_is_avail = false ;
                for bb = bbs
                    amps_avail = find(ProbeAmps(bb, :)) ;
                    if ismember(probej_amp, amps_avail)  
                        amp_is_avail = true ;
                        disp([pnames{jj} ':' anames{probej_amp} ' is available' ])
                        probesWithSpecifiedAmp{jj} = cat(2, probesWithSpecifiedAmp{jj}, bb) ;
                    end
                end
                isok = ~filterHoxBasedOnProbeAvailability || (isok && amp_is_avail) ;
            end
            
            % Stop checking if the combination is bad
            if isok
                % For this probe, look at adjacent probes
                adj_probes = find(probes(jj, :)) ;

                % Check if each adjacent probe is assigned an adjacent amplifier
                for adj_probe = adj_probes
                    % Stop checking if the combination is bad
                    if isok
                        % What amplifier does this adjacent probe have?
                        adj_amp = combo(adj_probe) ;

                        disp([pnames{jj} ':' anames{probej_amp} ' checking adjacent probe #' num2str(adj_probe) '--> ' pnames{adj_probe} ':' anames{adj_amp}])
                        % Does this adjacent probe have an adjacent amplifier?
                        bad_amps = find(ampfrs(adj_amp, :)) ;

                        % Consider this adjacent probe's amplifier. Is current 
                        % probe's amplifier a member of possible_adj_amps?
                        isok = isok && ~any(ismember(probej_amp, bad_amps)) ;

                        if isok
                            disp(['Ok: Probe ' pnames{jj} ':' anames{probej_amp} ' is not adjacent to ' pnames{adj_probe} ':' anames{combo(adj_probe)}])
                        else
                            disp(['Bad Combo: Probe ' pnames{jj} ':' anames{probej_amp} ' has probe adjacent to ' pnames{adj_probe} ':' anames{combo(adj_probe)}])
                            disp('...trying again with new combo...')
                        end

                    end
                end
            else
                disp(['Bad amplifier choice for ' pnames{jj} ': ' anames{probej_amp} ' not available for this probe.'])
                disp('...trying again with new combo...')
            end
        end
    end

    % Now check if the probes B1,B2, etc are uniquely assigned
    if ensureUniqueProbes
        result = findUniqueCombination(probesWithSpecifiedAmp, []) ;
        if isempty(result)
            disp('Invalid sequence: NO UNIQUE COMBINATION OF PROBES for target amplifier sequence.')
            isok = false ;
        else
            disp('Found unique combination of probes for target amplifier sequence')
        end
    end

    % Was this combo ok?
    if isok
        disp(['Found good combination: ' output])
        goodCombos(ii) = true ;
        goodOutput{dmyk} = output ;
        goodProbeSetWithSpecifiedAmps{dmyk} = probesWithSpecifiedAmp ; 
        dmyk = dmyk + 1 ;
    end
end

% These are all the valid combinations with adjacent colors spatially
% separated.
goodCombos = combinations(goodCombos, :) ;
nGoodCombos = size(goodCombos, 1) ;
disp(['Found ' num2str(nGoodCombos) ' valid combinations.'])

%% Plot valid combinations 
fig = figure('units', 'centimeters', 'position', [10,10,40,4]) ;
rr = colors(:, 1) ;
gg = colors(:, 2) ;
bb = colors(:, 3) ;
RR = rr(goodCombos) ;
GG = gg(goodCombos) ;
BB = bb(goodCombos) ;
imagesc(cat(3, RR, GG, BB))
% Make label on left: amplifiers
xticklabels(pnames)
yticks(1:nGoodCombos) ;
yticklabels(goodOutput)
set(gcf, 'color', 'w')
% Also label probes
ylabels = cell(1, nGoodCombos) ;
for ii = 1:dmyk
    hold on;
    plot([0, nprobes+1], [ii-0.5, ii-0.5], 'k-', 'linewidth', 2)
    if ii <= nGoodCombos
        formatNameValuePairs(goodProbeSetWithSpecifiedAmps{ii}, pnames)
        ylabels{ii} = [goodOutput{ii}, ' // ', formatNameValuePairs(goodProbeSetWithSpecifiedAmps{ii}, pnames)]  ;
    end
end
yticklabels(ylabels) ;

for ii = 1:nGoodCombos
    ylabels{ii}
end

saveas(gcf, 'validHCRcombinations.png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = findUniqueCombination(subsets, chosen_elements)
    % Define a recursive function to find a valid combination of elements
    % from the set of subsets, such that one element is chosen from each
    % subset with no elements repeated in multiple subsets. 
    % Note that this function returns the first valid combination rather
    % than looking for all valid combinations. 
    % 
    % Parameters
    % ----------
    % subsets : cell of N arrays 
    %   sets of elements for which to find unique element for each set
    % chosen_elements : 
    %   starting point for finding unique sets
    % 
    % Returns
    % -------
    % result : N x 1 array
    %   a valid unique combination of sets, with result(i) being the chosen
    %   element from the ith subset
    %
    % Example usage
    % -------------
    % % Define the set of subsets A
    % A = {[1,2], [3], [5,6], [8], [1,7]};
    % % Initialize a variable to keep track of the chosen elements
    % chosen_elements = [];
    % result = findUniqueCombination(A, chosen_elements)
    if isempty(subsets)
        result = chosen_elements;
        return;
    end

    for i = 1:length(subsets{1})
        element = subsets{1}(i);

        if ~any(element == chosen_elements)
            new_chosen_elements = [chosen_elements, element];
            new_subsets = subsets(2:end);
            result = findUniqueCombination(new_subsets, new_chosen_elements);
            if ~isempty(result)
                return;
            end
        end
    end

    result = [];
end

function formatted_string = formatNameValuePairs(cell_array, names)
    % Initialize an empty cell array to store formatted name-value pairs
    formatted_pairs = cell(1, numel(cell_array));

    % Loop through the cell_array and names cell arrays
    for i = 1:numel(cell_array)
        element = cell_array{i};
        name = names{i};

        if iscell(element)
            % Convert cell array to a string
            element_str = strjoin(cellfun(@num2str, element, 'UniformOutput', false), ',');
        else
            % Convert other data types to a string
            element_str = num2str(element);
            element_str = strrep(element_str, '  ', '/') ;
        end

        % Create the formatted name-value pair
        formatted_pairs{i} = [name, ':', element_str];
    end

    % Combine the formatted pairs into a single string
    formatted_string = strjoin(formatted_pairs, ' ');
end

