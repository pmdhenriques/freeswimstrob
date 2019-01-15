% Video tracking script for free-swimming zebrafish larvae using the
% darkfield/brightfield stroboscopic illimination (still unnammed...)
%
% Pedro Henriques, September 2016
%
% 160912_1745 - tail peak finder now uses only the two largest peaks

version = '160912_1745';

prompt = {'Save videos?','Process Brightfield?','Process Darkfield','Debug plots?','D/B interval'};
defaultans = {'0','1','1','0','8'};
answer = inputdlg(prompt,'Input',1,defaultans);

savevideo = str2num(answer{1});  % Save darkfield//brightfield and background subtracted videos
processB = str2num(answer{2});   % Process brightfield videos
processD = str2num(answer{3});   % Process darkfield (paramaecia segmentation) videos
debug = str2num(answer{4});  % Plot fish tracking live
dfi = str2num(answer{5});    % Darkfield/Brightfield interval (frames)

if debug; [X,Y] = meshgrid(1:stats.width,1:stats.height); end

folders = uipickfiles('FilterSpec','G:\freeswim_wt');

for f = 1:length(folders)
    avi_files = dir([char(folders(f)) filesep '*.avi']);
    if isempty(avi_files)
        disp('No .avi files in directory')
        break
    else
        for a = 1:length(avi_files)
            
            stats.version = version;
            stats.filepath = [char(folders(f)) filesep avi_files(a).name];
            stats.filedir = char(folders(f));
            filename = strsplit(avi_files(a).name,'.'); stats.filename  = filename{1};
            filenameparts = strsplit(stats.filename,'_');
            
            matfiles = dir([char(folders(f)) filesep '*.mat']); % Check if mat file of the video is already there
            matix = find(~cellfun(@isempty,strfind({matfiles.name},filenameparts{end})));
            if ~isempty(matix)
                disp(['.mat file for ' stats.filepath ' was already found!'])
                load([char(folders(f)) filesep matfiles(matix(1)).name],'stats')
                if processD
                    disp('Replacing Darkfield structure')
                    stats = rmfield(stats,'Paramecia'); end
                if processB
                    disp('Replacing Brightfield structure')
                    stats = rmfield(stats,{'Body','Eyes'}); end
            end
            vsrc = VideoReader(stats.filepath);   % Open video file
            
            stats.height = vsrc.Height;  % Height
            stats.width = vsrc.Width;    % Width
            stats.fr = vsrc.FrameRate;   % Framerate
            stats.nfr = vsrc.Duration*vsrc.FrameRate;   % Number of frames
            npx = stats.height*stats.width;     % Number of pixels in a frame
            stats.dfi = dfi;        % darkfield/brightfield frame interval
            
            %% set frame separating threshold
            
            nvthrf = 100;   % Number of frames used to build D/B separating threshold
            stats.vthrlim = 200;    % Pixel size ROI on the border of the video to use for threshold determination
            thrs = zeros(nvthrf,1);
            
            k = 1; vsrc.CurrentTime = 0;
            h = waitbar(0.5,'Initializing');
            while hasFrame(vsrc)
                while k <= length(thrs)
                    waitbar(k/nvthrf,h,'Processing Threshold...');
                    video = readFrame(vsrc); video = video(1:stats.vthrlim,1:stats.vthrlim,1);
                    thrs(k) = median(reshape(video,stats.vthrlim^2,1));
                    k = k+1;
                end
                vsrc.CurrentTime = vsrc.Duration;
            end
            close(h)
            stats.BDthr = multithresh(thrs);
            
            fig0 = figure; histogram(thrs); hold on;
            ylimt = get(gca,'YLim'); line([stats.BDthr stats.BDthr],[ylimt(1) ylimt(2)]);
            title(stats.filename);
            
            clear thrs
            
            %% Build initial background models
            
            Bbkgfr = 40;      % Number of frames (in current format == seconds) used to contruct the brightfield background model
            Dbkgfr = 40;     % Frames for the darkfield background
            
            BkgMB = zeros(stats.height,stats.width,Bbkgfr,'single');     % Bright background
            BkgMD = zeros(stats.height,stats.width,Dbkgfr,'single');     % Dark background
            
            vsrc.CurrentTime = 0;
            k = 1; l = 1; b = 1; d = 1;
            h = waitbar(0.5,'Initializing');
            while hasFrame(vsrc)
                while l <= Dbkgfr || k <= Bbkgfr
                    waitbar((l+k-1)/(Dbkgfr+Bbkgfr),h,['Processing initial background model...' ...
                        ' k: ' num2str(k) ' l: ' num2str(l)]);
                    video = readFrame(vsrc); video = video(:,:,1);
                    videoresh = reshape(video(1:stats.vthrlim,1:stats.vthrlim),stats.vthrlim^2,1);
                    if median(videoresh) <= stats.BDthr
                        if k <= Bbkgfr && (b == 1 || mod(b,stats.fr) == 0)
                            BkgMB(:,:,k) = video;
                            k = k+1;
                        end
                        b = b+1;
                    else
                        if l <= Dbkgfr && (d == 1 || mod(d,stats.fr/stats.dfi) == 0)
                            BkgMD(:,:,l) = video;
                            l = l+1;
                        end
                        d = d+1;
                    end
                end
                vsrc.CurrentTime = vsrc.Duration;
            end
            close(h);
            
            %%
            
            bkgB = median(BkgMB,3); % Brightfield initial background model
            bkgD = median(BkgMD,3); % Darkfield backgrounf model
            
            fig1 = figure; colormap(fig1,'gray');
            subplot(1,2,1); imagesc(bkgB); title('Bright Background');
            subplot(1,2,2), imagesc(bkgD); title('Dark Background');
            mtit(stats.filename)
            
            %% Find arena and create arrays
            
            [arenaC, arenaR] = imfindcircles(bkgB,[420 440],'Sensitivity',0.99); % Determine the arena position
            if isempty(arenaC)
                arenaC = [stats.height/2 stats.width/2];
                arenaR = 434; end
            [arenaX, arenaY] = meshgrid(-(arenaC(1)-1):(stats.width-arenaC(1)),-(arenaC(2)-1):(stats.height-arenaC(2)));
            arena_mask = (arenaX.^2)+(arenaY.^2) <= (arenaR+95)^2;
            
            fig2 = figure; colormap(fig2,'gray');
            bkgBmask = bkgB; bkgBmask(~arena_mask) = 0; imagesc(bkgBmask); xlabel('Arena mask'); title(stats.filename);
            
            %%  Run the tracking script
            
            properties = {'Area','Centroid','Eccentricity','Orientation', ...
                'BoundingBox','MajorAxisLength','MinorAxisLength','PixelIdxList'};
            stats.bodythr = 0.02;     % Body threshold level
            stats.parathr = 0.035;   % Paramecia threshold level
            npts = 15;  % Number of points to build tail profile
            stats.tailsegR = 5;  % Radius of circle centered at each point (distance between each point)
            plotfrq = 5000; % Frame interval in which to plot tracking figures
            stats.bckupdf = 320;  % Background update frequency (frames)
            
            if savevideo
                vsnkdir = [stats.filedir filesep 'processed_videos'];
                if ~exist(vsnkdir,'file'); mkdir(vsnkdir); end
                vD = VideoWriter([vsnkdir filesep stats.filename '_D'],'Motion JPEG AVI');
                vDbkg = VideoWriter([vsnkdir filesep stats.filename '_Dbkg'],'Motion JPEG AVI');
                vB = VideoWriter([vsnkdir filesep stats.filename '_B'],'Motion JPEG AVI');
                vBbkg = VideoWriter([vsnkdir filesep stats.filename '_Bbkg'],'Motion JPEG AVI');
                open(vB); open(vD); open(vDbkg); open(vBbkg);   % Open video files to write to
            end
            
            Dix = false(stats.nfr,1);
            
            vsrc.CurrentTime = 0;
            k = 1; l = 1; n = 1;
            h = waitbar(0.5,'Initializing');
            while hasFrame(vsrc)
                %%  Brightfield frame
                waitbar(vsrc.CurrentTime/vsrc.Duration,h,['Processing ' num2str(vsrc.CurrentTime) '/' num2str(vsrc.Duration)]);
                video = readFrame(vsrc); video = video(:,:,1); % Get video frame and get only first dimention
                videoresh = reshape(video(1:stats.vthrlim,1:stats.vthrlim),stats.vthrlim^2,1);
                if median(videoresh) <= stats.BDthr    % Is Bright frame?
                    if processB
                        % Update background model
                        if mod(k,stats.bckupdf) == 0
                            BkgMB = circshift(BkgMB,-1,3);
                            BkgMB(:,:,end) = video;
                            bkgB = median(BkgMB,3);
                        end
                        imB = uint8(abs(single(video) - bkgB));     % Background subtraction
                        imB(~arena_mask) = 0;    % Get rid of stuff outside the arena
                        %         imB_filt = imgaussfilt(imB,3);    % Gaussian smooth filter
                        bodyBW = bwmorph(im2bw(imB,stats.bodythr),'clean');  % Morphological operations
                        bodyDBW = bwdist(~bodyBW);  % Distance trasform to use for midline segment analysis
                        if savevideo
                            writeVideo(vB,video);   % Write video frame
                            writeVideo(vBbkg,imB);  % Write image to video file
                        end
                        if any(any(bodyBW))     % Is there a body found?
                            %             imB(~bodyBW) = 0;    % Set background to 0
                            eyelvl = double(multithresh(imB(bodyBW),2));     % Get the eyes by histogram thresholding
                            eyesBW = imfill(im2bw(imB,eyelvl(2)/255),'holes');  % Apply eye threshold to find the eyes
                            s1 = regionprops(bodyBW,properties);    % Body props
                            s2 = regionprops(eyesBW,properties);    % Eyes props
                            [~,ix] = max([s1.Area]);    % Largest binary object
                            s1 = s1(ix);
                            s1.MeanI = mean(imB(s1.PixelIdxList));
                            
                            % Get heading direction from eyes and sb positions
                            
                            neyes = size(s2,1);
                            [~,ix] = sort([s2.Area],'descend');
                            s2 = s2(ix);
                            if neyes >= 3
                                s2 = s2(1:3); % Process only the 3 largest objects
                                p = combnk(1:3,2);   % Combinations of objects to calculate distances
                                eye_c = vertcat(s2.Centroid);   % Object centroids
                                dd = NaN(3,1);   % Object distance vector
                                for j = 1:3
                                    dd(j) = abs(sqrt((diff([eye_c(p(j,1),1),eye_c(p(j,2),1)]).^2) + ...
                                        ((diff([eye_c(p(j,1),2),eye_c(p(j,2),2)]).^2))));    % Calculate distances between each object centroid
                                end
                                [~, ix] = min(dd);  % Minumum distance should correspond to distance between both eyes
                                eye_ix = p(ix,:)';   % Eyes ix
                                sb_ix = setdiff(p,eye_ix);  % SB ix is the other
                                
                                %%
                                eye_xy1 = eye_c(eye_ix(1),:);   % eye 1 centroid
                                eye_xy2 = eye_c(eye_ix(2),:);   % eye 2 centroid
                                sb_cxy = eye_c(sb_ix,:); % sb centroid
                                cp = [(eye_xy1(1)+eye_xy2(1))/2, (eye_xy1(2)+eye_xy2(2))/2];    % Coordenate of middle point between both eyes
                                hvec = [cp(1)-sb_cxy(1), cp(2)-sb_cxy(2), 0];   % Vector between sb and eye middle point
                                b_ori = atan2d(hvec(2),hvec(1));    % Orientation of sb-eye vector
                                if b_ori < 0    % Get values in standart 0-360deg
                                    b_ori = abs(-180 - b_ori)+180; end
                                s1.Angle = b_ori;
                                s1.SB_Centroid = sb_cxy;
                                
                                % Which eye is which?
                                
                                s2 = s2(eye_ix);    % Process only eyes from now on
                                
                                eyeu = [eye_xy1(1)-eye_xy2(1), eye_xy1(2)-eye_xy2(2)];  % eye-eye vector
                                eye_head_ang = atan2d(hvec(1)*eyeu(2)-hvec(2)*eyeu(1), ...
                                    hvec(1)*eyeu(1)+hvec(2)*eyeu(2));
                                
                                % eye identity can be calculated from the polarity
                                % of the vector that connects both eyes and that of the body,
                                % as this depends on which eye the vector starts from.
                                % If starts from the left eye, direction is negative,
                                % and positive if starting from the right eye.
                                
                                if eye_head_ang >= 0
                                    side = logical([1,0]);
                                else
                                    side = logical([0,1]);
                                end
                                
                                % Calculate eye-body angle from orientation, major axis
                                % length and centroids of the regions
                                
                                eyes_ori = cat(1,s2.Orientation);
                                eyes_malen = cat(1,s2.MajorAxisLength);
                                eyes_cent = cat(1,s2.Centroid);
                                
                                x1 = eyes_cent(:,1)+(cos(deg2rad(eyes_ori)).*eyes_malen./2);
                                x2 = eyes_cent(:,1)-(cos(deg2rad(eyes_ori)).*eyes_malen./2);
                                y1 = eyes_cent(:,2)-(sin(deg2rad(eyes_ori)).*eyes_malen./2);
                                y2 = eyes_cent(:,2)+(sin(deg2rad(eyes_ori)).*eyes_malen./2);
                                
                                u = [x1-x2, y1-y2, [0;0]]; % major axis vectors for both eyes
                                
                                eye_body_ang = NaN(1,2);
                                for i = 1:2
                                    angle = atan2d(norm(cross(u(i,:),hvec)),dot(u(i,:),hvec));
                                    if abs(angle) > 90  % Correct for incorrect vector directions, which give high angles
                                        angle = sign(angle)*180-angle; end
                                    eye_body_ang(i) = angle;
                                end
                                
                                % Write to data structure
                                
                                L = s2(~side);  % Left eye stats
                                R = s2(side);   % Right eye stats
                                L.Angle = eye_body_ang(~side);
                                R.Angle = eye_body_ang(side);
                                
                                % Get tail curvature
                                
                                pts = NaN(npts,2);
                                pts(1,:) = sb_cxy;  % First point is sb centroid
                                
                                for i = 1:npts-1
                                    if i == 1
                                        cx = stats.tailsegR*2*cos(0:0.1:2*pi)+pts(i,1);  % Draw a circle
                                        cy = stats.tailsegR*2*sin(0:0.1:2*pi)+pts(i,2);
                                    else
                                        cx = stats.tailsegR*cos(0:0.1:2*pi)+pts(i,1);
                                        cy = stats.tailsegR*sin(0:0.1:2*pi)+pts(i,2); end
                                    [cxix,cyix,cI] = improfile(bodyDBW,cx,cy,'bilinear');   % Build intencity profile of the circle perimeter
                                    cI2 = circshift(cI,round(length(cI)/4));    % shift intensity profile to correct for peaks that are in the extremities
                                    warning('off','all')
                                    [pks,locs] = findpeaks(cI,'MinPeakDistance',pi*stats.tailsegR,'MinPeakHeight',1);  % Find peaks
                                    [pks2,locs2] = findpeaks(cI2,'MinPeakDistance',pi*stats.tailsegR,'MinPeakHeight',1);   % Find peaks in the shifted vector
                                    warning('on','all')
                                    
                                    if length(pks) < 2 && length(pks2) >= 2     % Add a peak in the extremety
                                        pks = pks2; locs = [1; locs]; end
                                    
                                    if length(pks) > 2
                                        % Use only two largest peaks
                                        [pks, ix] = sort(pks,'descend');
                                        pks = pks(1:2);
                                        locs = locs(ix(1:2)); end
                                    
                                    if length(pks) == 2
                                        if i == 1
                                            dd = abs(sqrt((cxix(locs)-cp(1)).^2 + (cyix(locs)-cp(2)).^2));  % Second point is defined by the maximum distance to the eye-eye center point
                                        else
                                            dd = NaN(length(pks),1);    % Next iterations look for points that maximize the distance to the previous point
                                            for j = 1:length(pks)
                                                dd(j) = abs(sqrt((cxix(locs(j))-pts(i-1,1))^2 + (cyix(locs(j))-pts(i-1,2))^2));
                                            end
                                        end
                                        [~, ix] = max(dd);
                                        pts(i+1,:) = [cxix(locs(ix)),cyix(locs(ix))];
                                    else
                                        % Terminate for loop if cant find more than one
                                        % point
                                        break
                                    end
                                end
                                
                                % Calculate angles between body segments                                
                                for i = 1:size(pts,1)-1
                                    if i == 1
                                        u1 = hvec;  % first segment is body heading axis
                                    else
                                        u1 = [pts(i-1,1)-pts(i,1), pts(i-1,2)-pts(i,2),0];
                                    end
                                    u2 = [pts(i,1)-pts(i+1,1), pts(i,2)-pts(i+1,2),0];
                                    
                                    angle = atan2d(u1(1)*u2(2)-u1(2)*u2(1),u1(1)*u2(1)+u1(2)*u2(2));
                                    pts(i,3) = angle;
                                end
                                
                                s1.Tail_pts = pts;
                                s1.Tail_angle = nansum(pts(:,3));
                                
                                % Write to data structure
                                
                                stats.Body(k+l-1) = s1;
                                stats.Eyes.L(k+l-1) = L;
                                stats.Eyes.R(k+l-1) = R;
                                
                                if k == 1 ...
                                        %                                         || mod(k,plotfrq) == 0
                                    fig3 = figure;
                                    imagesc(imB); colormap(fig3,'gray');
                                    hold on;
                                    scatter(pts(:,1),pts(:,2),'.g');
                                    scatter(s2(side).Centroid(1),s2(side).Centroid(2),'ro');
                                    scatter(s2(~side).Centroid(1),s2(~side).Centroid(2),'bo');
                                    axis([cp(1)-(150+10),cp(1)+150+10,cp(2)-(150+10),cp(2)+150+10]);
                                    xlabel(['k = ' num2str(k)]);
                                    title(stats.filename);
                                    drawnow
                                    hold off
                                end
                                
                                % plot stuff for debugging purposes
                                
                                if debug
                                    if exist('f1','var') == 1; delete(f1); end
                                    if exist('f2','var') == 1; delete(f2); end
                                    [lidx,ridx] = visual_field_idx(8.5,cp,hvec,stats.Eyes.L(k+l-1).Centroid, ...
                                        [x1(~side),y1(~side);x2(~side),y2(~side)], ...
                                        stats.Eyes.R(k+l-1).Centroid, ...
                                        [x1(side),y1(side);x2(side),y2(side)], ...
                                        X,Y);
                                    BL = bwboundaries(lidx);
                                    BR = bwboundaries(ridx);
                                    imagesc(video); colormap gray
                                    hold on;
                                    scatter(pts(:,1),pts(:,2),'.g');
                                    scatter(s2(side).Centroid(1),s2(side).Centroid(2),'r.');
                                    scatter(s2(~side).Centroid(1),s2(~side).Centroid(2),'b.');
                                    f1 = patch(BL{:}(:,2),BL{:}(:,1),'b','FaceAlpha',0.3);
                                    f2 = patch(BR{:}(:,2),BR{:}(:,1),'r','FaceAlpha',0.3);
                                    axis([0,stats.width,0,stats.height]);
                                    drawnow
                                    hold off;
                                end
                                
                            else
                                % Get orientation through BBox method or other
                                % Needs to be written
                            end
                        end
                    end
                    k = k+1;
                else
                    %% Darkfield frame
                    if processD
                        Dix(k+l-1) = 1;
                        % Update background model
                        if mod(l,stats.bckupdf/stats.dfi) == 0
                            BkgMD = circshift(BkgMD,-1,3);
                            BkgMD(:,:,end) = video;
                            bkgD = median(BkgMD,3);
                        end
                        imD = uint8(abs(single(video) - bkgD));
                        imD(~arena_mask) = 0;    % Get rid of stuff outside the arena
                        paraBW = bwmorph(bwmorph(bwmorph(im2bw(imD,stats.parathr),'clean'),'close',Inf),'fill');
                        s3 = regionprops(paraBW,properties);
                        s3 = s3([s3.Area] <= 150);
                        
                        for m = 1:size(s3,1)
                            s3(m).Px_Intensity = imD(s3(m).PixelIdxList);
                        end
                        
                        stats.Paramecia(l).Frame = k+l+1;
                        stats.Paramecia(l).Parameters = s3;
                        
                        % plot some stuff
                        
                        if l == 1 ...
                                %                                 || mod(l,round(plotfrq/stats.dfi)) == 0
                            fig4 = figure; colormap(fig4,'parula');
                            imagesc(imD);
                            hold on
                            cnts = vertcat(s3.Centroid);
                            if ~isempty(cnts)
                                scatter(cnts(:,1),cnts(:,2),100,'ro')
                            end
                            xlabel(['l = ' num2str(l)])
                            title(stats.filename);
                            hold off
                        end
                        
                        if savevideo
                            writeVideo(vD,video);
                            writeVideo(vDbkg,imD);
                        end
                    end
                    l = l+1;
                end
            end
            stats.FramesIdx = Dix;
            close(h)
            if savevideo
                close(vD); close(vB); close(vBbkg); close(vDbkg); end
            
            %% save data structure
            
            save([stats.filedir filesep stats.filename '_' char(datetime('today','Format','yyMMdd'))],'stats','-v7.3');
        end
        clear stats
    end
end