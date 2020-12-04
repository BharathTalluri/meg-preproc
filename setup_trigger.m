function trig = setup_trigger()

trig.zero = 0;
trig.width = 0.005; 

trig.block_start = 01;  % start of block
trig.block_end = 02;    % end of block

trig.trial_start = 11;
trig.int1_start = 12;
trig.interm_resp_cue = 13;
trig.interm_delay = 14;
trig.int2_start = 15;
trig.estim_resp_cue = 16;
trig.trial_end = 17;
trig.rest_end = 18;

trig.sample_on = 21;   % onset of individual sample
trig.sample_off = 22;  % offset of individual sample

trig.attention_cue_left = 31;
trig.attention_cue_right = 32;

trig.interm_resp_click = 41;
trig.interm_resp_left = 42;    % 'left' response
trig.interm_resp_right = 43;   % 'right' response
trig.interm_resp_no = 44;     % no response (either double press, or task-irrelevant button)
trig.interm_resp_bad = 45;     % bad response (either double press, or task-irrelevant button)

trig.estim_resp = 51;
trig.estim_resp_no = 52;
trig.estim_resp_bad = 53;

end
