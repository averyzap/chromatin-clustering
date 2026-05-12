% Returns list of state configurations in order used
% right now only doing p=3
function configs = state_configurations(p)
    %assert(p==3);
    %configs = [1 2 3; 2 1 3; 3 2 1; 1 3 2]';
    configs = sc_sub([], p)';
end

function out = sc_sub(determined, p)
    out = [];
    % determined is a single list of who each current is bound to
    
    if length(determined) == p
        out = determined;
        return
    end
    % figure out all possibilities for next value
    
    % first, check if next column is already determined by previous
    next_col = length(determined)+1;
    f = find(determined==next_col);
    if ~isempty(f)
        out = sc_sub([determined f], p);
        return
    else
        % we can use anything that isn't yet taken (all previous columns
        % and taken numbers
        for next_val=next_col:p
            if ismember(next_val, determined)
                continue
            else
                out = [out; sc_sub([determined next_val], p)];
            end
        end
    end
end