classdef EddyData
% EDDYDATA  Loads and preprocesses the eddy lifetime data files.
%
%   Usage:
%       d = eddyplot.EddyData('../../../data/Lifetime-last14days.txt');
%       d = d.load();              % parse, sort, zero-base time
%       d = d.filterHalf('first'); % keep first 14 days (optional)
%
%   Key properties (read-only after load):
%       time_min    – event times in minutes, zero-based, sorted
%       lifetime    – eddy lifetime (minutes)
%       x, y        – spatial coordinates (grid units)
%       mag         – amplitude  (aa column)
%       radius      – radius     (bb column)
%       label       – filesystem-safe stem, e.g. 'Lifetime-last14days'
%       formalLabel – human-readable title,  e.g. 'Last 14 Days'
%
%   Data file format (whitespace-delimited, 8 columns):
%       time_raw  lifetime  x  y  mag  radius  c7  c8
%   where time_raw is "day.hour.minute" (e.g. 205.01.30).

    properties (SetAccess = private)
        filepath    (1,1) string
        time_min    (:,1) double
        lifetime    (:,1) double
        x           (:,1) double
        y           (:,1) double
        mag         (:,1) double
        radius      (:,1) double
        windowMin   (1,1) double = NaN
        label       (1,1) string = "unknown"      % filesystem-safe stem
        formalLabel (1,1) string = "Unknown Data" % human-readable title
    end

    methods
        %------------------------------------------------------------------
        function obj = EddyData(filepath)
            arguments
                filepath (1,1) string
            end
            obj.filepath = filepath;
        end

        %------------------------------------------------------------------
        function obj = load(obj)
            if ~isfile(obj.filepath)
                error('EddyData:fileNotFound','File not found: %s', obj.filepath);
            end

            fid = fopen(obj.filepath, 'r');
            raw = textscan(fid, '%s %d %d %d %f %f %d %d');
            fclose(fid);

            time_raw     = raw{1};
            obj.lifetime = double(raw{2});
            obj.x        = double(raw{3});
            obj.y        = double(raw{4});
            obj.mag      = raw{5};
            obj.radius   = raw{6};

            % Parse timestamps using the same convention as the estimators.
            numericTimes = zeros(numel(time_raw), 1);
            for ii = 1:numel(time_raw)
                numericTimes(ii) = eddyplot.EddyData.parseTimeStr(time_raw{ii});
            end
            obj.time_min = numericTimes;

            % Sort events chronologically.
            [obj.time_min, sortIdx] = sort(obj.time_min);
            obj.lifetime = obj.lifetime(sortIdx);
            obj.x        = obj.x(sortIdx);
            obj.y        = obj.y(sortIdx);
            obj.mag      = obj.mag(sortIdx);
            obj.radius   = obj.radius(sortIdx);
            fprintf('Data sorted by timestamp (%d entries)\n', numel(obj.time_min));

            % Label (filesystem) and formalLabel (human-readable)
            [~, fname]      = fileparts(obj.filepath);
            obj.label       = string(fname);
            obj.formalLabel = eddyplot.EddyData.makeFormalLabel(fname);

            % Use midnight at the start day as time zero, preserving empty
            % 15-minute observation blocks before, between, and after events.
            originMin = floor(obj.time_min(1) / (24 * 60)) * (24 * 60);
            obj.time_min = obj.time_min - originMin;
            if contains(lower(fname), '14days')
                obj.windowMin = 14 * 24 * 60;
            else
                obj.windowMin = (floor(max(obj.time_min) / 15) + 1) * 15;
            end
        end

        %------------------------------------------------------------------
        function obj = filterHalf(obj, half)
            arguments
                obj  eddyplot.EddyData
                half (1,1) string {mustBeMember(half, ["first","last"])}
            end
            half_dur = 14 * 24 * 60;
            if half == "first"
                mask            = obj.time_min < half_dur;
                obj.label       = obj.label + "-first14days";
                obj.formalLabel = obj.formalLabel + " — First Half";
            else
                mask            = obj.time_min >= half_dur;
                obj.label       = obj.label + "-last14days";
                obj.formalLabel = obj.formalLabel + " — Second Half";
            end
            obj = obj.applyMask(mask);
        end

        %------------------------------------------------------------------
        function n = numEvents(obj);  n = numel(obj.time_min); end
        function T = totalDuration(obj); T = obj.windowMin; end
    end

    %----------------------------------------------------------------------
    methods (Access = private)
        function obj = applyMask(obj, mask)
            obj.time_min = obj.time_min(mask);
            obj.lifetime = obj.lifetime(mask);
            obj.x        = obj.x(mask);
            obj.y        = obj.y(mask);
            obj.mag      = obj.mag(mask);
            obj.radius   = obj.radius(mask);
        end
    end

    %----------------------------------------------------------------------
    methods (Static, Access = private)
        function t = parseTimeStr(s)
            % Convert day.hour.minute timestamps to absolute minutes.
            dotPos = strfind(s, '.');
            if numel(dotPos) >= 2
                day    = str2double(s(1            : dotPos(1)-1));
                hour   = str2double(s(dotPos(1)+1 : dotPos(2)-1));
                minute = str2double(s(dotPos(2)+1 : end));
                t = day*24*60 + hour*60 + minute;
            else
                t = NaN;
            end
        end

        function fl = makeFormalLabel(stem)
            % Map known filename stems to human-readable paper titles.
            map = { ...
                'Lifetime-first14days',         'First 14 Days'; ...
                'Lifetime-last14days',          'Last 14 Days'; ...
                'Lifetime-lasthalf-first1day',  'First Day (Last Half)'; ...
                'Lifetime-lasthalf-first2days', 'First Two Days (Last Half)'; ...
                'Lifetime-synthetic',           'Synthetic Dataset'; ...
            };
            for k = 1:size(map,1)
                if strcmpi(stem, map{k,1})
                    fl = map{k,2};
                    return;
                end
            end
            % Fallback: prettify
            fl = strtrim(strrep(strrep(stem, 'Lifetime-', ''), '-', ' '));
            fl(1) = upper(fl(1));
        end
    end
end
