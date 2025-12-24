function [guess_offset] = update_offset(old_offset)
%UNTITLED5 Summary of this function goes here
    if old_offset==1 || old_offset==2 || old_offset==3 
        guess_offset=old_offset+1;
    else
        guess_offset=1;
    end
end