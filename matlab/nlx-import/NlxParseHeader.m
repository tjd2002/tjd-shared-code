function hdrstruct = NlxParseHeader(hdr, param);
hdrstruct = struct();
[parstr valstr] = strtok(hdr);

% Remove the varying amounts of whitespace at ends
parstr = strtrim(parstr);
valstr = strtrim(valstr);

for par_val = [parstr valstr hdr]',
    p = par_val{1}; % param token
    v = par_val{2}; % rest of line
    h = par_val{3}; % original hdr line

    % ignore blank lines
    if  isempty(p)
        continue
    end

    % handle comment lines
    if p(1) == '#'
        % special cases: file name, open/close times
        try
            if strmatch('## File Name', h)
                p = 'FileName';
                v = strtrim(h(14:end));
            elseif strmatch('## Time Opened', h)
                p = 'TimeOpened';
                v = datestr( h([25:35 findstr(h,'At Time:')+9:end-1]));
            elseif strmatch('## Time Closed', h)
                p = 'TimeClosed';
                v = datestr(h([25:35 findstr(h,'At Time:')+9:end-1]));
            else
                % ignore all other comments
                continue
            end
        catch
            warning('failed to extract filename/date line');
            continue
        end

        % handle non-comment, non-empty lines
    else
        if ~isempty(v)
            % try to convert value to number
            vnum = str2double(v);
            if ~isnan(vnum) % numeric conversion worked
                v = vnum;
            else
                % try to convert True/False to logical
                if strcmpi(v,'True')
                    v = true;
                elseif strcmpi(v,'False')
                    v = false;
                end
            end
        end
        % if no conversion was done, v is handled as a string

        % clean up param names
        % remove leading '-' from param name
        if p(1) == '-', p(1) = ''; end

        % remove trailing ':'
        if ~isempty(p) && p(end) == ':', p(end) = ''; end
    end

    % ensure p is a valid field name
    p = genvarname(p, fieldnames(hdrstruct));

    % write our param/value pair to the output struct
    hdrstruct.(p) = v;

end