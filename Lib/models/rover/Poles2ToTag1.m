%% transform position wr2 poles to pos wr2 tag1
function Pt = Poles2ToTag1(Pp)
    
    % shifts
    shx = 1.16 - (-2);
    shy = 0;
    shz = 2.20 - 2.10;

    % transform
    d = [shx shy shz].';
    A = eye(3);

    % compute
    Pt = A*Pp +d;

end