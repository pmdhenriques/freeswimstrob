% Extraction of features from processed movies (after running
% free_swim_run)
%
% Pedro Henriques, UCL, September 2016

folders = uipickfiles('FilterSpec','G:\freeswim_wt');

for f = 1:length(folders)
    disp(['Processing folder ' num2str(f) '/' num2str(length(folders))])
    matfiles = dir([char(folders(f)) filesep '*.mat']); % Check if mat file of the video is already there
    if ~isempty(matfiles)
        for matf = 1:length(matfiles)
            filepath = [char(folders(f)) filesep matfiles(matf).name];
            disp(['Loading file ' num2str(matf) '/' num2str(length(matfiles))])
            load(filepath)
            if exist('stats','var')
                %% Interpolate tracking variables onto aquisition frequency
                
                disp('Interpolating...')
                
                intix = cellfun(@isempty,{stats.Body.Area});    % Find missing data frames (can occur both because of darkfield frames or lack of body or eye data
                
                bcnts = vertcat(stats.Body(~intix).Centroid);   % body centroids
                lcnts = vertcat(stats.Eyes.L(~intix).Centroid); % left eye centroids
                rcnts = vertcat(stats.Eyes.R(~intix).Centroid); % right eye centroids
                tp = horzcat(stats.Body(~intix).Tail_pts);  % tail points
                
                % Interpolate
                
                bxint = interp1(find(~intix),bcnts(:,1),1:numel(intix),'linear','extrap');
                byint = interp1(find(~intix),bcnts(:,2),1:numel(intix),'linear','extrap');
                bxyint = [bxint;byint]';
                bomjaint = interp1(find(~intix),vertcat(stats.Body(~intix).MajorAxisLength),1:numel(intix),'linear','extrap');
                bomnaint = interp1(find(~intix),vertcat(stats.Body(~intix).MinorAxisLength),1:numel(intix),'linear','extrap');
                boriint = interp1(find(~intix),vertcat(stats.Body(~intix).Orientation),1:numel(intix),'linear','extrap');
                bangint = interp1(find(~intix),vertcat(stats.Body(~intix).Angle),1:numel(intix),'linear','extrap');
                tangint = interp1(find(~intix),vertcat(stats.Body(~intix).Tail_angle),1:numel(intix),'linear','extrap');
                
                lxint = interp1(find(~intix),lcnts(:,1),1:numel(intix),'linear','extrap');
                lyint = interp1(find(~intix),lcnts(:,2),1:numel(intix),'linear','extrap');
                lxyint = [lxint;lyint]';
                lorisint = interp1(find(~intix),vertcat(stats.Eyes.L(~intix).Orientation),1:numel(intix),'linear','extrap');
                langsint = interp1(find(~intix),vertcat(stats.Eyes.L(~intix).Angle),1:numel(intix),'linear','extrap');
                lmajlint = interp1(find(~intix),vertcat(stats.Eyes.L(~intix).MajorAxisLength),1:numel(intix),'linear','extrap');
                lminlint = interp1(find(~intix),vertcat(stats.Eyes.L(~intix).MinorAxisLength),1:numel(intix),'linear','extrap');
                
                rxint = interp1(find(~intix),rcnts(:,1),1:numel(intix),'linear','extrap');
                ryint = interp1(find(~intix),rcnts(:,2),1:numel(intix),'linear','extrap');
                rxyint = [rxint;ryint]';
                rorisint = interp1(find(~intix),vertcat(stats.Eyes.R(~intix).Orientation),1:numel(intix),'linear','extrap');
                rmajlint = interp1(find(~intix),vertcat(stats.Eyes.R(~intix).MajorAxisLength),1:numel(intix),'linear','extrap');
                rminlint = interp1(find(~intix),vertcat(stats.Eyes.R(~intix).MinorAxisLength),1:numel(intix),'linear','extrap');
                rangsint = interp1(find(~intix),vertcat(stats.Eyes.R(~intix).Angle),1:numel(intix),'linear','extrap');
                
                tpxint = interp1(find(~intix),tp(:,1:3:end)',1:numel(intix),'linear','extrap');
                tpyint = interp1(find(~intix),tp(:,2:3:end)',1:numel(intix),'linear','extrap');
                tpaint = interp1(find(~intix),tp(:,3:3:end)',1:numel(intix),'linear','extrap');
                
                % Replace missing frames with interpolated data
                
                intix = find(intix);
                h = waitbar2(0.5,'Initializing');
                for i = intix
                    waitbar2(i/intix(end),h,'Interpolating missing frames');
                    
                    stats.Body(i).Angle = bangint(i);
                    stats.Body(i).Tail_angle = tangint(i);
                    stats.Body(i).Tail_pts = [tpxint(i,:);tpyint(i,:);tpaint(i,:)]';
                    stats.Body(i).Orientation = boriint(i);
                    stats.Body(i).Centroid = bxyint(i,:);
                    stats.Body(i).MajorAxisLength = bomjaint(i);
                    stats.Body(i).MinorAxisLength = bomnaint(i);
                    stats.Eyes.L(i).Centroid = lxyint(i,:);
                    stats.Eyes.L(i).Orientation = lorisint(i);
                    stats.Eyes.L(i).Angle = langsint(i);
                    stats.Eyes.L(i).MajorAxisLength = lmajlint(i);
                    stats.Eyes.L(i).MinorAxisLength = lminlint(i);
                    stats.Eyes.R(i).Centroid = rxyint(i,:);
                    stats.Eyes.R(i).Orientation = rorisint(i);
                    stats.Eyes.R(i).Angle = rangsint(i);
                    stats.Eyes.R(i).MajorAxisLength = rmajlint(i);
                    stats.Eyes.R(i).MinorAxisLength = rminlint(i);
                end
                close(h)
                
                clearvars -except stats f folders matfiles matf filepath
                
                %% Find moving object (paramecia) tracks
                
                disp('Finding tracks...')
                
                ndfr = size(stats.Paramecia,2);
                paraCentroids = cell(ndfr,1);
                for i = 1:ndfr
                    paraCentroids{i} = vertcat(stats.Paramecia(i).Parameters.Centroid);
                end
                
                ix = find(cellfun(@isempty,paraCentroids));
                if ~isempty(ix)
                    for i = ix'
                        paraCentroids{i} = [nan nan];
                    end
                end
                
                maxlinkdist = 10;   % maximum distance to look for tracks to link (pixels)
                maxgapclose = 3;    % maximum number of gap frames to link tracks before start a new one
                trakerdebug = false;
                
                % Use simple tracker to find object tracks
                tracks = simpletracker(paraCentroids, ...
                    'MaxLinkingDistance',maxlinkdist, ...
                    'MaxGapClosing',maxgapclose, ...
                    'Debug',trakerdebug);
                
                % Do some filtering by length of tracks found
                
                tlen = cellfun(@sum,cellfun(@not,cellfun(@isnan,tracks,'UniformOutput',false),'UniformOutput',false));
                %                 [y,x] = ecdf(tlen);
                %                 tracklenthr = knee_pt(y,x);
                tracks = tracks(tlen >= 8);
                ntrks = size(tracks,1);
                
                % Build cell array of centroids for each track
                
                cTracks = cell(ntrks,1);
                for j = 1:size(tracks,1)
                    t1 = NaN(ndfr,2);
                    for i = 1:ndfr
                        if ~isnan(tracks{j}(i))
                            t1(i,:) = paraCentroids{i}(tracks{j}(i),:);
                        end
                    end
                    cTracks{j} = t1;
                end
                
                % Do some more filtering based on net velocity for each track (try to
                % exclude particles that keep the same place
                
                % dxy = NaN(size(cTracks,1),1);
                % for i = 1:length(dxy)
                %     dxy(i) = nanmean(abs(diff(sqrt((diff(cTracks{i}(:,1)).^2)+(diff(cTracks{i}(:,2)).^2)))));
                % end
                %
                % ix = dxy >= prctile(dxy,1);
                % cTracks = cTracks(ix);
                % ntrks = size(cTracks,1);
                
                % Interpolate onto acquisition frequency
                
                if isfield(stats,'PTracksInterp')
                    stats = rmfield(stats,'PTracksInterp');
                end
                
                warning('off','all');
                h = waitbar2(0.5,'Initializing');
                for i = 1:ntrks
                    waitbar2(i/ntrks,h,'Interpolating PTracks');
                    xb = interp1(find(stats.FramesIdx),cTracks{i}(:,1),1:stats.nfr,'linear');
                    x = single(interp1(find(stats.FramesIdx),cTracks{i}(:,1),1:stats.nfr,'spline'));
                    x(isnan(xb)) = nan;
                    y = single(interp1(find(stats.FramesIdx),cTracks{i}(:,2),1:stats.nfr,'spline'));
                    y(isnan(xb)) = nan;
                    stats.PTracksInterp{i} = [x;y]';
                end
                close(h)
                warning('on','all');
                
                stats.TracksIx = tracks;    % Save tracks ixs
                
                clearvars -except stats f folders matfiles matf filepath
                
                %% Find bouts
                
                disp('Finding bouts...')
                
                stats.Bouts.btthr = 12;  % Bout threshold
                stats.Bouts.min_intbtlen = 2;   % min inter-bout length
                stats.Bouts.min_btlen = 4;  % min bout length
                
                ang = [stats.Body.Tail_angle];  % Tail angles
                angvf = smooth(abs([0 diff(ang)]),11,'lowess')';  %smooth tail velocity profile
                bv = angvf >= stats.Bouts.btthr;  % bouts
                btst = find([0 diff(bv)] == 1); % Bout start
                btend= find([0 diff(bv)] == -1);    % Bout end
                if bv(1) == 1; btst = [1 btst]; end     % Correct for when starts in a bout
                if bv(end) == 1; btend = [btend numel(bv)]; end     % Correct for when end in a bout
                intbtlen = diff([btend;[btst(2:end) 0]]); intbtlen(end) = [];   % inter-bout lenghts
                ix = find(intbtlen <= stats.Bouts.min_intbtlen);    % find bouts that don't obey min intbtlen
                btend(ix) = []; btst(ix+1) = [];    % delete those bouts
                btlen = diff([btst;btend]);   % bout lengths
                ix = find(btlen <= stats.Bouts.min_btlen);  % find bouts that don't obey min_btlen
                btst(ix) = []; btend(ix) = [];  % delete bouts
                
                stats.Bouts.Start = btst;
                stats.Bouts.End = btend;
                stats.Bouts.nbts = numel(btst); % number of bouts
                
                clearvars -except stats f folders matfiles matf filepath
                
                %% Find and filter convergences
                
                disp('Finding convergences...')
                
                stats.min_intconvlen = 10; % minimum number of frames between convergences
                stats.min_convlen = 20;   % minumum convergence length (frames)
                
                vang = sgolayfilt([stats.Eyes.L.Angle]+[stats.Eyes.R.Angle],1,7);   % filtered vergence angle
                stats.convthr = multithresh(vang); % convergencec threshold. Around 45 is a good value
                hvix = vang > stats.convthr;  % convergences
                hvst = find([0 diff(hvix)] == 1);   % start of convergences
                hvend = find([0 diff(hvix)] == -1); % end of convergences
                
                if hvix(1) == 1
                    hvst = [1 hvst]; end
                if hvix(end) == 1
                    hvend = [hvend numel(hvix)]; end
                
                intconvlen = diff([hvend;[hvst(2:end) 0]]); intconvlen(end) = [];   % inter-convergence lengths
                ix = find(intconvlen <= stats.min_intconvlen);    % find convergences that don't obey min intconvlen
                hvend(ix) = []; hvst(ix+1) = [];    % delete those convergences
                convlen = diff([hvst;hvend]);   % convergences length
                ix = find(convlen <= stats.min_convlen);  % find convergences that don't obey min_convlen
                hvst(ix) = []; hvend(ix) = [];  % delete convergences
                nconv = numel(hvst);    % number of convergences
                
                for i = 1:nconv
                    stats.convergences(i).Start = hvst(i);
                    stats.convergences(i).End = hvend(i);
                end
                
                %
                
                for c = 1:nconv
                    
                    st = findnearest(stats.convergences(c).Start,stats.Bouts.Start,0);
                    fin = findnearest(stats.convergences(c).End,stats.Bouts.End,1);
                    if isempty(fin)
                        fin = findnearest(stats.convergences(c).End,stats.Bouts.End,-1);
                    end
                    stf = stats.Bouts.Start(st);
                    endf = stats.Bouts.End(fin);
                    
                    btix = st:fin;
                    
                    k = 1;
                    for b = btix
                        stats.convergences(c).Ethog(k).Bout_St = stats.Bouts.Start(b);
                        stats.convergences(c).Ethog(k).Bout_End = stats.Bouts.End(b);
                        k = k+1;
                    end
                end
                
                clearvars -except stats f folders matfiles matf filepath
                
                %%  Process convergence bout features
                
                disp('Processing convergence features...')
                
                stats.vfdelta = 8.5;
                stats.deltalim = 0.250*stats.fr;
                stats.ldR = 230; % local density circle radius (centered on paramecia)
                stats.ldtau = 1/30;
                stats.objheadingfrms = 16;
                
                nconv = size(stats.convergences,2);
                ntrks = size(stats.PTracksInterp,2);
                
                langf = smooth([stats.Eyes.L.Angle],11,'lowess');
                rangf = smooth([stats.Eyes.R.Angle],11,'lowess');
                tang = [stats.Body.Tail_angle];
                tvf = smooth(abs([0 diff(tang)]),11,'lowess');
                bang = rad2deg(unwrap(deg2rad([stats.Body.Angle])));
                
                bmjaxlen = smooth([stats.Body.MajorAxisLength],11,'lowess');
                
                h = waitbar2(0.5,'Initializing');
                for c = 1:nconv
                    waitbar2(c/nconv,h,'Processing convergence features...');
                    for b = 1:length(stats.convergences(c).Ethog)
                        
                        st = stats.convergences(c).Ethog(b).Bout_St;
                        fin = stats.convergences(c).Ethog(b).Bout_End;
                        
                        if fin+stats.deltalim <= stats.nfr && st-stats.deltalim >= 0
                            
                            lcnts = vertcat(stats.Eyes.L(st-stats.objheadingfrms:st).Centroid);
                            rcnts = vertcat(stats.Eyes.R(st-stats.objheadingfrms:st).Centroid);
                            
                            cps = [(lcnts(:,1)+rcnts(:,1))./2,(lcnts(:,2)+rcnts(:,2))./2];   % eye central points
                            cp = cps(end,:);  % central point
                            
                            hvec = [cp(1)-stats.Body(st).Tail_pts(1,1), cp(2)-stats.Body(st).Tail_pts(1,2)];  % Orientation vector
                            
                            cl = stats.Eyes.L(st).Centroid; cr = stats.Eyes.R(st).Centroid;   % left eye centroid
                            
                            extl1 = [cl(1)+(cos(deg2rad(stats.Eyes.L(st).Orientation)).*stats.Eyes.L(st).MajorAxisLength./2), ...
                                cl(2)-(sin(deg2rad(stats.Eyes.L(st).Orientation)).*stats.Eyes.L(st).MajorAxisLength./2)];
                            extl2 = [cl(1)-(cos(deg2rad(stats.Eyes.L(st).Orientation)).*stats.Eyes.L(st).MajorAxisLength./2), ...
                                cl(2)+(sin(deg2rad(stats.Eyes.L(st).Orientation)).*stats.Eyes.L(st).MajorAxisLength./2)];
                            extr1 = [cr(1)+(cos(deg2rad(stats.Eyes.R(st).Orientation)).*stats.Eyes.R(st).MajorAxisLength./2), ...
                                cr(2)-(sin(deg2rad(stats.Eyes.R(st).Orientation)).*stats.Eyes.R(st).MajorAxisLength./2)];
                            extr2 = [cr(1)-(cos(deg2rad(stats.Eyes.R(st).Orientation)).*stats.Eyes.R(st).MajorAxisLength./2), ...
                                cr(2)+(sin(deg2rad(stats.Eyes.R(st).Orientation)).*stats.Eyes.R(st).MajorAxisLength./2)];
                            extl = [extl1;extl2]; extr = [extr1;extr2];
                            
                            pts = NaN(ntrks,2,'single');
                            for t = 1:ntrks
                                pts(t,:) = stats.PTracksInterp{t}(st,:);
                            end
                            
                            [lidx,ridx] = visual_field_idx(stats.vfdelta,cp,hvec,cl,extl,cr,extr,pts(:,1),pts(:,2));
                            
                            % Track related metrics
                            
                            if isfield(stats.convergences(c).Ethog,'Track')
                                
                                track = stats.convergences(c).Ethog(b).Track;
                                if any(track)
                                    % interpolate nan values of track
                                    trackpts = stats.PTracksInterp{track};
                                    trackpts = interp1(find(~isnan(trackpts(:,1))),trackpts(~isnan(trackpts(:,1)),:),1:size(trackpts,1),'splice','extrap');
                                    cnt = trackpts(st,:);
                                    
                                    % Exp local density
                                    
                                    dd = sqrt(((bsxfun(@minus,pts(lidx | ridx,1),cnt(1))).^2) + ...
                                        ((bsxfun(@minus,pts(lidx | ridx,2),cnt(2))).^2));   % Distance of all observable paramecia to target
                                    dd(dd == 0) = [];
                                    
                                    stats.convergences(c).Ethog(b).trk_expld = nansum(exp(-dd.*stats.ldtau));   % local target density
                                    
                                    tracku = [trackpts(st,1) - trackpts(st-stats.objheadingfrms,1), ...
                                        trackpts(st,2) - trackpts(st-stats.objheadingfrms,2)]; % target heading vector
                                    ufp = [cnt(1)-cp(1) cnt(2)-cp(2)];
                                    
                                    stats.convergences(c).Ethog(b).f2tang = atan2d(hvec(1).*ufp(:,2)-hvec(2).*ufp(:,1), ...
                                        hvec(1).*ufp(:,1)+hvec(2).*ufp(:,2));   % fish 2 track angle
                                    
                                    stats.convergences(c).Ethog(b).f2thang = atan2d(hvec(1).*tracku(:,2)-hvec(2).*tracku(:,1), ...
                                        hvec(1).*tracku(:,1)+hvec(2).*tracku(:,2)); % fish 2 target angle difference
                                    
                                    f2tpD = sqrt(((bsxfun(@minus,pts(:,1),cp(1))).^2) + ...
                                        ((bsxfun(@minus,pts(:,2),cp(2))).^2));  % fish 2 all paramecia distances
                                    
                                    [~,rnkix] = sort(f2tpD,'ascend');    % rank distances
                                    rnkix = setdiff(rnkix,find(~(lidx | ridx)),'stable');   % remove non-visible objects
                                    trnkix = find(rnkix == track);    % target rank
                                    if ~isempty(trnkix)
                                        stats.convergences(c).Ethog(b).trnk = trnkix;
                                    else
                                        stats.convergences(c).Ethog(b).trnk = nan; end
                                    
                                    stats.convergences(c).Ethog(b).f2tD = f2tpD(track);   % fish 2 target distance
                                    stats.convergences(c).Ethog(b).f2trD = min(f2tpD(lidx | ridx))/f2tpD(track);    % fish 2 target relative distance
                                    
                                    if lidx(track) && ridx(track); stats.convergences(c).Ethog(b).trackvf = 3; % Both vfs
                                    elseif lidx(track); stats.convergences(c).Ethog(b).trackvf = 1;    % left vf
                                    elseif ridx(track); stats.convergences(c).Ethog(b).trackvf = 2;    % right vf
                                    else stats.convergences(c).Ethog(b).trackvf = 0;   % not visible
                                    end
                                    
                                    stats.convergences(c).Ethog(b).abs_t_sp = ...
                                        nanmean(diff(sqrt((diff(trackpts(st-stats.objheadingfrms:st,1)).^2) + ...
                                        (diff(trackpts(st-stats.objheadingfrms:st,2)).^2))));    % mean absolute speed
                                    stats.convergences(c).Ethog(b).t_loomsp = ...
                                        nanmean(diff(sqrt(((trackpts(st-stats.objheadingfrms:st,1)-cps(:,1)).^2) + ...
                                        (trackpts(st-stats.objheadingfrms:st,2)-cps(:,2)).^2)));    % mean loom speed
                                else
                                    stats.convergences(c).Ethog(b).trk_expld = nan;
                                    stats.convergences(c).Ethog(b).abs_t_sp = nan;
                                    stats.convergences(c).Ethog(b).t_loomsp = nan;
                                    stats.convergences(c).Ethog(b).f2tD = nan;
                                    stats.convergences(c).Ethog(b).f2trD = nan;
                                    stats.convergences(c).Ethog(b).trackvf = nan;
                                    stats.convergences(c).Ethog(b).f2tang = nan;
                                    stats.convergences(c).Ethog(b).f2thang = nan;
                                    stats.convergences(c).Ethog(b).trnk = nan;
                                end
                            end
                            
                            ix = find(sqrt((bsxfun(@minus,pts(:,1),cp(1)).^2) + ...
                                (bsxfun(@minus,pts(:,2),cp(2)).^2)) <= stats.ldR & ...
                                (lidx | ridx));
                            
                            allld = NaN(ntrks,1);
                            
                            k = 1;
                            for p = find(~isnan(pts(:,1)))'
                                ld = sqrt(((bsxfun(@minus,pts(:,1),pts(p,1))).^2) + ...
                                    ((bsxfun(@minus,pts(:,2),pts(p,2))).^2));
                                ld(ld == 0) = [];
                                allld(p) = nansum(exp(-ld.*stats.ldtau));
                            end
                            
                            stats.convergences(c).Ethog(b).mean_expld_fish = nanmean(allld(ix));    % fish local density
                            stats.convergences(c).Ethog(b).mean_expld_total = nanmean(allld);    % arena density
                            stats.convergences(c).Ethog(b).nobj_fish = numel(ix);   % num fish centered visible objects
                            stats.convergences(c).Ethog(b).nobj_total = sum(~isnan(pts(:,1)));  % Total number of tracks
                            
                            % Body
                            
                            stats.convergences(c).Ethog(b).BodyAng = stats.Body(st).Angle;
                            
                            bcnts = vertcat(stats.Body(st:fin).Centroid);
                            stats.convergences(c).Ethog(b).BodyCnt = bcnts(1,:);
                            dd = sqrt((diff(bcnts(:,1)).^2)+(diff(bcnts(:,2)).^2));
                            ddv = abs(diff(dd));
                            
                            stats.convergences(c).Ethog(b).Total_Displacemet = nansum(dd);
                            stats.convergences(c).Ethog(b).Mean_Bout_Vel = nanmean(ddv);
                            stats.convergences(c).Ethog(b).Max_Bout_Vel = max(ddv);
                            
                            % Eyes
                            
                            stats.convergences(c).Ethog(b).Vang = nanmean(langf(st:fin)+rangf(st:fin));
                            stats.convergences(c).Ethog(b).Lang = nanmean(langf(st:fin));
                            stats.convergences(c).Ethog(b).Rang = nanmean(rangf(st:fin));
                            stats.convergences(c).Ethog(b).dVang = nanmean(langf(st:st+20)+rangf(st:st+stats.deltalim)) - ...
                                nanmean(langf(st-stats.deltalim:st)+rangf(st-stats.deltalim:st));
                            stats.convergences(c).Ethog(b).dLang = nanmean(langf(st:st+stats.deltalim)) - ...
                                nanmean(langf(st-stats.deltalim:st));
                            stats.convergences(c).Ethog(b).dRang = nanmean(rangf(st:st+stats.deltalim)) - ...
                                nanmean(rangf(st-stats.deltalim:st));
                            
                            % Tail
                            
                            stats.convergences(c).Ethog(b).TailAng = nanmean(tang(st:fin));
                            stats.convergences(c).Ethog(b).TailVig = nanmean(tvf(st:fin));
                            stats.convergences(c).Ethog(b).TailAcc = nanmean(abs(diff(tvf(st:fin))));
                            
                            % Tail frequency
                            
                            fs = stats.fr;
                            x = abs(tang(st:fin));
                            m = length(x); n = pow2(nextpow2(m));
                            y = fft(x,n);           % DFT
                            freqr = (0:n-1)*(fs/n);     % Frequency range
                            power = y.*conj(y)/n;   % Power of the DFT
                            [pks,locs] = findpeaks(power(1:floor(n/2)));
                            if ~isempty(pks)
                                [~,ix] = max(pks);
                                stats.convergences(c).Ethog(b).TailFreq = freqr(locs(ix));
                            else
                                stats.convergences(c).Ethog(b).TailFreq = nan; end
                            
                            % Major axis length
                            
                            stats.convergences(c).Ethog(b).MjAxLen = nanmean(bmjaxlen(st:fin));
                            stats.convergences(c).Ethog(b).dMjAxLen = nanmean(bmjaxlen(st:st+stats.deltalim)) -  ...
                                nanmean(bmjaxlen(st-20:st));
                            
                            % Body angle delta
                            
                            stats.convergences(c).Ethog(b).dBodyAng = bang(fin)-bang(st);
                            
                        end
                    end
                end
                close(h)
                
                stats.ProcessDate = str2double(char(datetime('today','Format','yyMMdd')));
                
                disp('Saving...')
                save(filepath,'stats','-v7.3');
                
                clearvars -except f folders matfiles matf filepath
            end
        end
    end
end
