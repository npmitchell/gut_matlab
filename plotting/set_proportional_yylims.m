function set_proportional_yylims(y0, adjust_right_or_left)
%SET_PROPORTIONAL_YYLIMS(y0, adjust_right_or_left)
% Set the ylims of right and left dual axes so that they are proportional
% to each other wrt zero or supplied offset
%
% Parameters
% ----------
% y0 : float
%   the y value to align in both right and left yyaxes
% adjust_right_or_left : str ('right' or 'left' or 'both')
%   adjust the right axis or adjust the left axis   
%
% NPMitchell 2020

% Set yvalue at which both axis are aligned
if nargin < 1
    y0 = 0;
end

yyaxis right
rlims = ylim ;

yyaxis left
llims = ylim ;

% Set right proportionality to align with left
ll = y0 - llims(1) ;
LL = llims(2) - y0 ;
rr = y0 - rlims(2) ;
RR = rlims(2) - y0 ;

if strcmpi(adjust_right_or_left, 'right')
    
    % Do lims enclose y0?
    if all([ll, LL, rr, RR] > 0)
        % Adjust RR or rr depending on whether R/(R+r) >/< L/(L+l)
        Rfrac = RR / (RR + rr) ;
        Lfrac = LL / (LL + ll) ;
        if Rfrac < Lfrac
            RR = LL * rr / ll ;
        elseif Rfrac > Lfrac
            rr = RR * ll / LL ;
        else
            error('handle case here')
        end
    elseif all([ll, LL, RR] > 0)
        %    .______. R
        %    |      |
        % y0 -      |
        %    |      |
        %    .______. 
        %            }-r
        %           - y0
        %
        rr = ll * RR / LL ;
    else
        error('handle case here')
    end
    % Set the limits
    yyaxis right
    ylim([y0 - rr, y0 + RR])
    
elseif strcmpi(adjust_right_or_left, 'left')
    
    % Do lims enclose y0?
    if all([ll, LL, rr, RR] > 0)        
        % Adjust LL or ll depending on whether R/(R+r) >/< L/(L+l)
        Rfrac = RR / (RR + rr) ;
        Lfrac = LL / (LL + ll) ;
        if Lfrac < Rfrac
            LL = RR * ll / rr ;
        elseif Lfrac > Lfrac
            ll = LL * rr / RR ;
        else
            error('handle case here')
        end
    elseif all([rr, RR, LL] > 0)
        rr = ll * RR / ll ;
    end
    % Set the limits
    yyaxis left
    ylim([y0 - ll, y0 + LL])
    
elseif strcmpi(adjust_right_or_left, 'both')
    
    error('handle case -- use both left and right for optimal adjustment')

end
    



