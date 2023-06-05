% driver to get model valuse with compute_vars
clear all;
dat = load('./MealSim/05-Jun-2023_driverMeal_insulin-1_Kin-0_notes-MealOnly_computevars.mat');
do_insulin = 1;
do_FF = 0;

vals = compute_vars(dat.t, dat.y, dat.params, ...
                            'do_insulin', do_insulin,...
                            'do_FF', do_FF);

save_vals = input('save vals? (0/1) ');
if save_vals
    %fname = './MealSim/OriginalModel_MealOnly_computevars.mat';
    save(fname, 'vals', 'dat')
end