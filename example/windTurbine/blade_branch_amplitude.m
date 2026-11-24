function amplitude=blade_branch_amplitude(bd,m)
% Convert the sampled tip trajectory before taking its cycle maximum.
if strcmp(m.direction,'torsion')
    samples=coco_bd_col(bd,'HT_x',false);
    amplitude=cellfun(@(x) max(abs(blade_tip_response(m,m.scale(m.tip).*x.'))),samples);
else
    amplitude=coco_bd_col(bd,'Max_x')*m.scale(m.tip);
end
end
