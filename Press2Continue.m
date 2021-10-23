function [] = Press2Continue(Screen,Key2Continue,run_mode)
%this function moves to next screen when the correct key is pressed.
%Screen - the screen that we record the key on
%Key2Continue - the only key that we want to make us move to the next
%screen.

if run_mode == 0
    pause(); Pressed_key = get(Screen,'CurrentCharacter');
    
    while strcmpi(Pressed_key,Key2Continue) == 0
        pause(); Pressed_key = get(Screen,'CurrentCharacter');
    end
    
    clf;
    axis off
    
elseif run_mode == 1
    time_sti = 2*rand();                       %random time according to settings.
    pause(time_sti);
    
    clf;
    axis off
end
    
