function [lidx,ridx] = visual_field_idx(delta,cp,hvec,cl,extl,cr,extr,X,Y)
% Calculates points belonging to each eye Visual Field. Diagrams of vectors
% and points in the notebook.

delta = deg2rad(delta);   % Angle from eye horizontal line that corresponds to start of visual field (163deg total VF; Easter et al., 1996)

% Decide where eye vectors start from maximum distance to body centroid

hvec2 = hvec.*100;
p = [cp(1)+hvec2(1),cp(2)+hvec2(2)];
ddl = abs(sqrt((extl(:,1)-p(1)).^2 + (extl(:,2)-p(2)).^2));
ddr = abs(sqrt((extr(:,1)-p(1)).^2 + (extr(:,2)-p(2)).^2));
[~, lix] = min(ddl);
[~, rix] = min(ddr);

lpt1 = extl(lix,:);  % Left top eye point
rpt1 = extr(rix,:);  % Right ''
ul = [lpt1(1)-cl(1), lpt1(2)-cl(2)];    % Left eye vector
ur = [rpt1(1)-cr(1), rpt1(2)-cr(2)];

%   Rotate vectors around eye centroids

ul1 = [ul(1)*cos(delta) + ul(2)*sin(delta), ...
    -ul(1)*sin(delta) + ul(2)*cos(delta)];
ul2 = [ul(1)*cos(pi-delta) + ul(2)*sin(pi-delta), ...
    -ul(1)*sin(pi-delta) + ul(2)*cos(pi-delta)];
ur2 = [ur(1)*cos(pi+delta) + ur(2)*sin(pi+delta), ...
    -ur(1)*sin(pi+delta) + ur(2)*cos(pi+delta)];
ur1 = [ur(1)*cos(-delta) + ur(2)*sin(-delta), ...
    -ur(1)*sin(-delta) + ur(2)*cos(-delta)];

% Define rotated points 3 and 4

lpt2 = [cl(1)+ul1(1), cl(2)+ul1(2)];
lpt3 = [cl(1)+ul2(1), cl(2)+ul2(2)];
rpt2 = [cr(1)+ur1(1), cr(2)+ur1(2)];
rpt3 = [cr(1)+ur2(1), cr(2)+ur2(2)];

% Define coefficients of line passing through eye centroids and points 2 or 3

warning('off','all')
lcoeff1 = polyfit([cl(1),lpt2(1)], [cl(2),lpt2(2)], 1);
lcoeff2 = polyfit([cl(1),lpt3(1)], [cl(2),lpt3(2)], 1);
rcoeff1 = polyfit([cr(1),rpt2(1)], [cr(2),rpt2(2)], 1);
rcoeff2 = polyfit([cr(1),rpt3(1)], [cr(2),rpt3(2)], 1);
warning('on','all')

% Define angles of left and right eye vectors in standard coordenates

ul1ang = atan2d(ul1(2),ul1(1));
ul2ang = atan2d(ul2(2),ul2(1));
ur1ang = atan2d(ur1(2),ur1(1));
ur2ang = atan2d(ur2(2),ur2(1));

%%

% [X,Y] = meshgrid(1:stats.width,1:stats.height);
if abs(ul1ang) < 90
    lidx1 = Y <= lcoeff1(1).*X + lcoeff1(2);
else
    lidx1 = Y >= lcoeff1(1).*X + lcoeff1(2);
end

if abs(ul2ang) < 90
    lidx2 = Y >= lcoeff2(1).*X + lcoeff2(2);
else
    lidx2 = Y <= lcoeff2(1).*X + lcoeff2(2);
end

if abs(ur1ang) < 90
    ridx1 = Y >= rcoeff1(1).*X + rcoeff1(2);
else
    ridx1 = Y <= rcoeff1(1).*X + rcoeff1(2);
end

if abs(ur2ang) < 90
    ridx2 = Y <= rcoeff2(1).*X + rcoeff2(2);
else
    ridx2 = Y >= rcoeff2(1).*X + rcoeff2(2);
end

lidx = lidx1 & lidx2;
ridx = ridx1 & ridx2;

%% Debug

% plot(pts(:,1),pts(:,2),'go')
% hold on
% plot(cl(1),cl(2),'bo');
% plot(cr(1),cr(2),'ro');
% line([extl(1,1),extl(2,1)],[extl(1,2),extl(2,2)],'Color','b')
% line([extr(1,1),extr(2,1)],[extr(1,2),extr(2,2)],'Color','r')
% axis([0,stats.width,0,stats.height])
% set(gca,'YDir','reverse')
% plot(1:stats.width,[1:stats.width].*lcoeff1(1)+lcoeff1(2),'b')
% plot(1:stats.width,[1:stats.width].*lcoeff2(1)+lcoeff2(2),'b')
% plot(1:stats.width,[1:stats.width].*rcoeff1(1)+rcoeff1(2),'r')
% plot(1:stats.width,[1:stats.width].*rcoeff2(1)+rcoeff2(2),'r')
% plot(lpt1(1),lpt1(2),'k*')
% plot(rpt1(1),rpt1(2),'k*')
% hold off
    
end