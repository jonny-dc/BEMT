function [alpha_req] = Req_AoA(Cl,J,x,desig,m,c)

% Find the required angle of attack corresponding to a given lift
% coefficient for a certain hydrofoil section with chord c and maximum
% camber m.

% INPUTS:
% - Cl    : double
%           Required lift coefficient.
% - x     : double
%           Radius fraction.
% - desig : string
%           Section designation, valid entries are "WB" (Wageningen
%           B-series), "NACA" (NACA sections), or "C" (custom section
%           shape).
% - m     : double
%           Maximum section camber (m).
% - c     : double
%           Section chord length (m).

% OUTPUTS:
% - alpha_req : double
%               Required angle of attack corresponding to Cl (rad).

if desig == "WB"
    % O'Brien (1962) gives lift-camber slopes for Wageningen B-series
    % sections, which can be used to find zero-lift angles.
    rr = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
    Lm = [8.1 7.4 7.8 8.4 9.0 9.7 10.1 10.6];
    alpha0 = (interp1(rr,Lm,x,"makima")/(2*pi)) * (m/c);                    % Zero-lift angle
    alpha_req = (Cl/(2*pi)) - alpha0;                                       % Required angle of attack
elseif desig == "NACA"
    % Required angle of attack for NACA sections (still to do). Missing
    % zero-lift angle data for this section type.
elseif desig == "C"
    % Required angle of attack for custom section shapes (still to do).
else
    error("You have specified an invalid section name! Valid entries are " + ...
          "'WB' (Wageningen B-series sections), 'NACA' (NACA sections), or " + ...
          "'C' (custom section shapes, not recommended).")
end

% The lift-coefficient curve is not a true (one-to-one) function after
% stall, so caution needs to be exercised if the angle of attack is more
% than 10 degrees away from the zero-lift angle
if alpha_req > (10*(pi/180))-alpha0 || alpha_req < (-10*(pi/180))-alpha0
    warning("Required angle of attack may be outside the linear lift-curve " + ...
            "region! Proceed with caution!"+newline+ ...
            "J = "+num2str(J)+"x = "+num2str(x)+", alpha = "+num2str(alpha_req*180/pi))
end

end
