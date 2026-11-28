% Run COCO startup and QuasiP install before starting this batch.
for speed = [9.68, 10.89, 13.31, 14.52]
    fprintf('Starting edge FRC at %.2f rpm\n', speed);
    run_blade_frc(5, 'edge', 20, speed);
    fprintf('Finished edge FRC at %.2f rpm\n', speed);
end
