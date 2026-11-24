function tip=blade_tip_response(m,delta)
% Physical tip response from tip coordinate changes (one sample per column).
if ~strcmp(m.direction,'torsion')
    tip=delta;
    return
end
% Swing-twist decomposition of conj(e_eq)*e about the local span axis.
eq=m.q(m.tip); e=eq+delta; G=Base.genG(eq);
w=eq.'*e; x=G(1,:)*e;
assert(all(hypot(w,x)>1e-10),'Twist is undefined for a 180-degree pure swing.');
tip=atan2d(2*w.*x,w.^2-x.^2);
end
