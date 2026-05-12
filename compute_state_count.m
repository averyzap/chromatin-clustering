% computes how many binding states p beads have
function n = compute_state_count(p)
    nprev2 = 1;
    nprev = 2;
    for i = 3:p
        n = nprev + (i-1)*nprev2;
        nprev2 = nprev;
        nprev = n;
    end
    enda