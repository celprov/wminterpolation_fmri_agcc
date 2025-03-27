%1 to work with resting state data
%0 to work with working memory data
RestingState_mode = 1;  

if RestingState_mode 
    whole_brain_cap_one('RestingState')
    clear all; close all;

    whole_brain_cap_one('ControlsRS')
    clear all; close all;

%     whole_brain_cap_both(RestingState_mode);
else
    whole_brain_cap_one('WorkingMemory')
    clear all; close all;

    whole_brain_cap_one('ControlsWM')
    clear all; close all;

    whole_brain_cap_both(RestingState_mode);
end