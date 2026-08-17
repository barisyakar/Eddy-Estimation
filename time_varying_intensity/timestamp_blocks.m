function [timeBlocks, numericMinutes, originMinutes] = timestamp_blocks(timestamps)
%TIMESTAMP_BLOCKS Convert day.hour.minute timestamps to elapsed 15-min blocks.

if isstring(timestamps)
    timestamps = cellstr(timestamps);
end
if ~iscell(timestamps) || isempty(timestamps)
    error('timestamp_blocks:invalidInput', ...
          'timestamps must be a nonempty cell array or string array.');
end

numericMinutes = zeros(numel(timestamps), 1);
for i = 1:numel(timestamps)
    parts = sscanf(char(timestamps{i}), '%d.%d.%d');
    if numel(parts) ~= 3 || parts(2) < 0 || parts(2) >= 24 || ...
            parts(3) < 0 || parts(3) >= 60
        error('timestamp_blocks:invalidTimestamp', ...
              'Invalid day.hour.minute timestamp: %s', char(timestamps{i}));
    end
    numericMinutes(i) = parts(1) * 24 * 60 + parts(2) * 60 + parts(3);
end

originMinutes = floor(min(numericMinutes) / (24 * 60)) * (24 * 60);
rawBlocks = (numericMinutes - originMinutes) / 15;
timeBlocks = round(rawBlocks);

if any(abs(rawBlocks - timeBlocks) > 1e-10)
    error('timestamp_blocks:offGrid', ...
          'Every event timestamp must lie on the 15-minute observation grid.');
end

timeBlocks = double(timeBlocks);
end
