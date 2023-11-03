% Simulate contracting muscles 

clear all 
close all 

% Simulation parameters
numRows = 2;    % Number of rows in the grid
numCols = 20;    % Number of columns in the grid
numSteps = 1000;  % Number of simulation steps
k = 0.1;         % Spring constant
damping = 0.05;  % Damping factor
hh = 0.05 ;
stdDev = 0.1;    % Standard deviation for stochastic contractions
prob = 0.05 ;   % probability of contracting
plasticity = 0.001 / hh ;
refractoryPeriod = 20 * hh ;
pulseTime = 15 * hh ;
sigma = 5 ;

% Initialize grid positions and spring connections
[x, y] = meshgrid(0:(numCols-1), 0:(numRows-1));
xy = cat(2, x(:), y(:));
xy = xy - mean(xy) ;
xy0 = xy ; 
TRI = delaunay(xy) ;
BL = TRI2BL(TRI) ;
ba = bondAngles(BL, xy) ;
BL = BL(mod(abs(ba), pi/2) < 100*eps, :) ;
NL = BL2NL(BL) ;
KL = k * (NL > 0) ;
kL = k * ones(size(BL(:, 1))) ;
bL0 = vecnorm(xy(BL(:, 2), :) - xy(BL(:, 1), :), 2,2) ;
bL0matrix = bL2bLmatrix(NL,BL,bL0) ;
% check it
assert(all(all(bL0matrix == (KL > 0))))

% make boolean for only muscle springs
ba = bondAngles(BL, xy) ;
ismuscle = abs(ba - pi/2) < eps ;

% % Check it
% for row = 1:length(BL) 
%     pair = BL(row, :) ;
%     plot(xy(pair, 1), xy(pair, 2), '.-')
%     hold on;
% end
% close all

% Initialize figure
figure;
axis tight manual;
ax = gca;
ax.NextPlot = 'replaceChildren';

% Simulation loop
xyt = zeros(numSteps, size(xy, 1), size(xy, 2)) ;
contractionHistory = zeros(numSteps, size(BL, 1)) ;
refractoryHistory = zeros(numSteps, size(BL, 1)) ;

% timestamps
tt = (1:numSteps)*hh ;

for step = 1:numSteps
    disp(['t= ' num2str(tt(step))])
    
    % update rest bond length
    bLdiff = vecnorm(xy(BL(:, 2), :) - xy(BL(:, 1), :), 2,2) - bL0 ;
    bL0 = bL0 + bLdiff * plasticity ;
    
    % Apply stochastic contractions to springs
    pos = 0.5 * (xy(BL(:, 1), 1) + xy(BL(:, 2), 1)) ;
    contract_new = rand(length(bL0), 1)> (1-prob .* exp(-pos.^2/(2*sigma^2))) ;
    % contract_new = true(size(bL0)) ;
    
    contract_new(~ismuscle) = false ;
    
    if step == 1
        contract = contract_new ;
        refractory = zeros(size(contract)) ;
        contract_timer = contract*hh ;
        recovery_timer = refractoryPeriod * ones(size(contract)) ;
    else
        
        % refractory period
        refractory = recovery_timer < refractoryPeriod ;
        
        % Is this a new pulse?
        newPulse = (contract_new & ~refractory) ;
        
        
        % contract if new and after refractory period or if still pulsing
        contract = contract | newPulse & ~refractory ;
        
        
        % Update timers now
        % Note: recovery timer depends on contract timer state, but
        % contract timer does not depend on recovery timer. So update the
        % recovery timer first!        
        
        % count how long we've been contracting -- increment
        contract_timer_new = contract_timer + (contract>0)*hh;
        contract_timer = contract_timer_new ;
        % if we are done pulsing, reset contraction timer and turn off
        % contraction
        contract(contract_timer > pulseTime) = false ;
        contract_timer(contract_timer_new > pulseTime) = 0 ;
        
        % if we are done pulsing, set refractory timer
        recovery_timer(contract_timer_new > pulseTime) = 0 ;
        % if we are recovering, increment recovery timer
        recovery_timer(recovery_timer < refractoryPeriod) = recovery_timer(recovery_timer < refractoryPeriod) + hh ;
        
    end
    
    
    
    % Now define the rest lengths, but with active contractions as zero
    % rest length
    bL0now = bL0 ;
    bL0now(contract > 0) = 0 ;
    bL0matrix = bL2bLmatrix(NL,BL,bL0now) ;
    
    % Update positions and velocities
    [F, dF] = wFree_Energy4_dF(xy,NL,BL, bL0matrix, KL, bL0now,kL) ;
    
    % Euler
    % xy = xy + dF * damping ;
    
    % Runge-Kutta 4th Order
    % k1: tt(i)
    k1 = dF ;
    
    % k2: tt(i)+0.5*hh
    xy2 = xy + 0.5 * hh * k1 ;
    [F2, k2] = wFree_Energy4_dF(xy2, NL, BL, bL0matrix, KL,bL0now,kL);
    
    % k3: tt(i)+0.5*hh
    xy3 = xy + 0.5 * hh * k2 ;
    [F3, k3] = wFree_Energy4_dF(xy3, NL, BL, bL0matrix, KL,bL0now,kL);
   
    % k4: tt(i)+hh
    xy4 = xy + k3 * hh ;
    [F4, k4] = wFree_Energy4_dF(xy4, NL, BL, bL0matrix, KL,bL0now,kL); 
    
    % main equation
    xynew = xy - (1/6)*(k1 + 2*k2 + 2*k3 + k4) * hh;  
    xyt(step, :, :) = xynew ;    
    xy = xynew ;
    
    % Plot grid
    % clf;
    scatter(reshape(xyt(step, :, 1), [], 1), reshape(xyt(step, :, 2), [], 1), 'filled');
    ylim([-0.5, 0.5])
    xlim([min(xy0(:, 1)), max(xy0(:, 1))])
    drawnow;
    
    
    % Save all
    % Debugging
    contractionHistory(step, :) = contract ;
    refractoryHistory(step, :) = refractory ;
    
end
  
%%
figure; sgtitle('positions')
subplot(2, 2, 1) ;
plot(tt, xyt(:, 1, 1)-xy0(1, 1), '.-') ; hold on;
plot(tt, xyt(:, 1, 2)-xy0(1, 2), '.-') ; xlabel('time')
subplot(2, 2, 2) ;
plot(tt, xyt(:, 2, 1)-xy0(2, 1), '.-') ; hold on;
plot(tt, xyt(:, 2, 2)-xy0(2, 2), '.-') ; xlabel('time')
subplot(2, 2, 3) ;
plot(tt, xyt(:, 3, 1)-xy0(3, 1), '.-') ; hold on;
plot(tt, xyt(:, 3, 2)-xy0(3, 2), '.-') ; xlabel('time')
subplot(2, 2, 4) ;
plot(tt, xyt(:, 4, 1)-xy0(4, 1), '.-') ; hold on;
plot(tt, xyt(:, 4, 2)-xy0(4, 2), '.-') ; xlabel('time')

figure; sgtitle('Contraction state')
subplot(2, 2, 1) ; plot(tt, contractionHistory(:, 1), '.-') ; xlabel('time')
subplot(2, 2, 2) ; plot(tt, contractionHistory(:, 2), '.-') ; xlabel('time')
subplot(2, 2, 3) ; plot(tt, contractionHistory(:, 3), '.-') ; xlabel('time')
subplot(2, 2, 4) ; plot(tt, contractionHistory(:, 4), '.-') ; xlabel('time')

figure; sgtitle('Recovery state')
subplot(2, 2, 1) ; plot(tt, refractoryHistory(:, 1), '.-') ; xlabel('time')
subplot(2, 2, 2) ; plot(tt, refractoryHistory(:, 2), '.-') ; xlabel('time')
subplot(2, 2, 3) ; plot(tt, refractoryHistory(:, 3), '.-') ; xlabel('time')
subplot(2, 2, 4) ; plot(tt, refractoryHistory(:, 4), '.-') ; xlabel('time')

%%
save(fullfile('./muscle_contractions', 'output.mat'), 'xyt', 'contractionHistory', 'refractoryHistory', 'xy0')
  
