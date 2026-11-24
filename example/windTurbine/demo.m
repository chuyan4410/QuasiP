direction = 'edge';    % 'flap', 'edge', or 'torsion'
amplitude = 20;        % distributed load amplitude
                       % unit: N/m or N*m/m 
speed = 12.1;          % unit: rpm; default: 12.1 rpm
frc = run_blade_frc(5, direction, amplitude, speed);