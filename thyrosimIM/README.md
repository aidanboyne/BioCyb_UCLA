# Remote Desktop Guide

### Starting thyrosimIM:
1. Open Julia REPL package manager
2. Run command `activate thyrosimIM`
3. Return to main eval loop and run `using thyrosimIM`

### Run Sample Optimization
1. Load sample data: `sample_data = load_sample_data()`
2. Run parameter fitting: `solution = execute_fitting(sample_data, lsq_loss)`
3. View output: `view_estimate!(solution.minimizer, sample_data)`