function [Pos,ang] = PosnAnlys_4barmech(len1,len2,len3,len4,angcrank)
% This function gives positions of all joints with respect to a fixed joint
% and the angles of all links with respect to the fixed link for a fourbar
% mechanism.
% len1 must be smallest link length followed by other link lengths in the
% kinematic chain in sequence.
% The input link lengths must satisfy Grashoff condition.
% A fourbar chain can be assembled in two different ways.
% Both possible points and corresponding angles of the coupler
% are given as output.
% Output
% Pos = [O1xy Axy Bxy O2xy B'xy]
% ang = [angcrank angcoup angfollow angcoup' angfollow']
drivercrank = len1;
fixedlink = len4;
couplerlink = len2;
followerlink = len3;
if(isempty(angcrank)==1)
angcrank = (0:360/100:360)' + 45;
end
for ii = 1:length(angcrank)
   O1xy = [0 0];
   O2xy = [fixedlink 0];
   Axy = drivercrank*[cosd(angcrank(ii,1)) sind(angcrank(ii,1))];
   [xp,yp] = circcirc(Axy(1,1),Axy(1,2),couplerlink,O2xy(1,1),O2xy(1,2),followerlink);
   if(ii==1)
   Bxy = [xp(1,1) yp(1,1)];
   Bpxy = [xp(1,2) yp(1,2)];
   elseif(ii>1)
       pnt1 = [xp(1,1) yp(1,1)];
       pnt2 = [xp(1,2) yp(1,2)];
       dist1 = norm(pnt1-Bxyold);
       dist2 = norm(pnt2-Bxyold);
       if(dist1<dist2)
           Bxy = pnt1;
           Bpxy = pnt2;
       elseif(dist1>dist2)
           Bxy = pnt2;
           Bpxy = pnt1;
       end
   end
   Bxyold = Bxy;
   Pos(ii,:) = [O1xy Axy Bxy O2xy Bpxy];
   angcoup(ii,1) = rad2deg(find_ang(Bxy'-Axy'));
   angfollow(ii,1) = rad2deg(find_ang(Bxy'-O2xy'));
   angcoupp(ii,1) = rad2deg(find_ang(Bpxy'-Axy'));
   angfollowp(ii,1) = rad2deg(find_ang(Bpxy'-O2xy'));

end
% ang = [angcrank angcoup angfollow angcoup' angfollow']
ang = [angcrank angcoup angfollow angcoupp angfollowp];
end
