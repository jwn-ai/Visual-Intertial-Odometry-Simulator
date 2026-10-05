function [LM,S,info] = buildTrajectoryAwareMap(Rwb, rW, Rbc, fl, Wpx, Hpx, opts)
% Rwb : 3x3xN body-to-world rotations, e.g. RR
% rW  : 3xN vehicle world positions, e.g. x(1:3,:)
% Rbc : same extrinsic used by h(..., Rwb*Rbc, ...)
% fl  : focal length in pixels
% Wpx, Hpx : sensor width/height in pixels
%
% LM : 3xM selected static world landmarks

    if ~isfield(opts,'Kmin'),                opts.Kmin = 6; end
    if ~isfield(opts,'maxLandmarks'),        opts.maxLandmarks = 150; end
    if ~isfield(opts,'anchorStride'),        opts.anchorStride = 25; end
    if ~isfield(opts,'coverageStride'),      opts.coverageStride = 5; end
    if ~isfield(opts,'candidatesPerAnchor'), opts.candidatesPerAnchor = 30; end
    if ~isfield(opts,'depthRange'),          opts.depthRange = [5 20]; end
    if ~isfield(opts,'depthMin'),            opts.depthMin = 1.0; end
    if ~isfield(opts,'marginPx'),            opts.marginPx = 100; end
    if ~isfield(opts,'fovFraction'),         opts.fovFraction = 0.70; end

    N = size(rW,2);

    % Generate candidates from safe, central image regions at spaced poses.
    anchorIdx = 1:opts.anchorStride:N;
    if anchorIdx(end) ~= N
        anchorIdx(end+1) = N;
    end

    numCandidates = numel(anchorIdx) * opts.candidatesPerAnchor;
    candidates = zeros(3, numCandidates);

    safeHalfW = opts.fovFraction * (Wpx/2 - opts.marginPx);
    safeHalfH = opts.fovFraction * (Hpx/2 - opts.marginPx);

    c = 0;
    for ia = 1:numel(anchorIdx)
        n = anchorIdx(ia);

        Rwc = Rwb(:,:,n) * Rbc;

        for j = 1:opts.candidatesPerAnchor
            c = c + 1;

            u = (2*rand - 1) * safeHalfW;
            v = (2*rand - 1) * safeHalfH;
            z = opts.depthRange(1) + ...
                rand * (opts.depthRange(2) - opts.depthRange(1));

            % Back-project pixel/depth into camera coordinates.
            pc = [u*z/fl; v*z/fl; z];

            % Convert to a fixed world landmark.
            candidates(:,c) = rW(:,n) + Rwc * pc;

        end
    end

    % Visibility of every candidate at every time step.
    V = false(numCandidates, N);

    for n = 1:N
        Rwc = Rwb(:,:,n) * Rbc;

        pc = Rwc.' * (candidates - rW(:,n));

        z = pc(3,:);
        u = fl * pc(1,:) ./ z;
        v = fl * pc(2,:) ./ z;

        V(:,n) = (z > opts.depthMin & ...
                  abs(u) <= Wpx/2 - opts.marginPx & ...
                  abs(v) <= Hpx/2 - opts.marginPx).';
    end

    % Greedy multi-cover selection on a representative frame subset.
    coverIdx = 1:opts.coverageStride:N;
    if coverIdx(end) ~= N
        coverIdx(end+1) = N;
    end

    Vcover = V(:,coverIdx);
    count = zeros(1, numel(coverIdx));
    selected = false(numCandidates,1);

    while any(count < opts.Kmin) && sum(selected) < opts.maxLandmarks

        undercovered = count < opts.Kmin;

        % Candidate score: number of currently undercovered frames it covers.
        score = sum(Vcover(:,undercovered), 2);
        score(selected) = -Inf;

        [bestScore, bestId] = max(score);

        if bestScore <= 0
            warning(['No remaining candidate covers an undercovered frame. ', ...
                     'Increase candidatesPerAnchor or decrease anchorStride.']);
            break
        end

        selected(bestId) = true;
        count = count + double(Vcover(bestId,:));
    end

    LM = candidates(:,selected);

    % Verify across every time step, not only greedy-selection frames.
    coverageAll = sum(V(selected,:), 1);

    if any(coverageAll < opts.Kmin)
        warning('%d time steps have fewer than Kmin visible landmarks.', ...
            sum(coverageAll < opts.Kmin));
    end

    info.candidates = candidates;
    info.selected = find(selected);
    info.coverage = coverageAll;
    info.visibility = V(selected,:);

    S = nan(2,N,size(LM,2));
    for n=1:N
        n
        rn = rW(:,n);
        R = Rwb(:,:,n)*Rbc;    
        visibleIds = find(info.visibility(:,n)).';

        for k=visibleIds
        
            rk = LM(:,k);
            p = utils.TF(rk,R,rn);
            s = utils.PIz(p,fl);
            S(:,n,k) = s;
        end
    end
end