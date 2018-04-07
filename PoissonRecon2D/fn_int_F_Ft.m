function v = fn_int_F_Ft(F1, F2, t)
%int_F_Ft int F1(x) * F2(x - t) dx
F_int = fnmult(F1,fntrans(F2,t));
if F_int.pieces == 0
    v = 0;
else
    v = fnval(fnint(F_int),F_int.breaks(end));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tpp = fntrans(pp1, t)
%fnmult Translate function
tpp = pp1;
tpp.breaks = tpp.breaks - t;
end
