function [rIC, lIC, rTO, lTO] = read_gait_events(gaitFile)

n_rows = size(gaitFile,1);


pos_rIC=1;
pos_lIC=1;
pos_rTO=1;
pos_lTO=1;


for i=3:n_rows
        side = char(gaitFile(i,2)) 
        if side(1) == 'L'
            event = char(gaitFile(i,3))
            if event(1) == 'F'
                lIC(pos_lIC) = cell2mat(gaitFile(i,4));
                pos_lIC= pos_lIC+1;
            elseif event(1) == 'E'
                lTO(pos_lTO) = cell2mat(gaitFile(i,4));
                pos_lTO= pos_lTO+1;
            end
        elseif side(1) == 'R'
            event = char(gaitFile(i,3))
            if event(1) == 'F'
                rIC(pos_rIC) = cell2mat(gaitFile(i,4));
                pos_rIC= pos_rIC+1;
            elseif event(1) == 'E'
                rTO(pos_rTO) = cell2mat(gaitFile(i,4));
                pos_rTO= pos_rTO+1;
            end
        end
end

end