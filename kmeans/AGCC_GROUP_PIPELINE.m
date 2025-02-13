%1 to work with resting state data
%0 to work with working memory data
RestingState_mode = 1;  

if RestingState_mode 
    extractCAP_one('RestingState',30)
    clear all; close all;

    extractCAP_one('ControlsRS',30)
    clear all; close all;

    extractCAP_both(RestingState_mode);
else
    extractCAP_one('WorkingMemory',57)
    clear all; close all;

    extractCAP_one('ControlsWM',57)
    clear all; close all;

    extractCAP_both(RestingState_mode);
end