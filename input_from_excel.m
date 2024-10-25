clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cells = 'B13';%%$$$$$$$$$$$$$$$$$$$$
filename = 'C:\Users\Tanjib Ahmed\Pictures\Cancer Detection\variation of layer\variation of layer\1.385_32.csv';%$$$$$$$$
column_range = 'A1:A105';
datax = xlsread(filename, column_range);
datax = datax*(pi/180);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
column_range = 'B1:B105';
data1 = xlsread(filename, column_range);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'C:\Users\Tanjib Ahmed\Pictures\Cancer Detection\variation of layer\variation of layer\1.399_32.csv';%$$$$$$$$
column_range = 'B1:B105';
data2 = xlsread(filename, column_range);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
array_index = 105;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comb = 32;%$$$$$$$$
n_sm = [1.368, 1.392];%$$$$$$$$
%n_sm = [1.387, 1.401];
n_sm = [1.385, 1.399];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = length(n_sm);
FWHM = zeros(1,s);
for i = 1:s
    theta_x = datax*(180/pi);
    if i == 1
        R_y = data1;
    else
        R_y = data2;
    end
    figure(1)
    plot(theta_x(1:1:array_index),R_y(1:1:array_index), 'LineWidth', 2);
    xlabel('Incidence Angle (in degree)');
    ylabel('Reflectance(a.u.)');
    hold on;
    % 
    R_min = min(R_y(1:1:array_index));
    R_max = max(R_y(1:1:array_index));
    
    %Calculation of the FWHM
    R_halfmax = (R_max+R_min)/2;
    index1 = find(R_y(1:1:array_index) <= R_halfmax, 1, 'first');
    index2 = find(R_y(1:1:array_index) <= R_halfmax, 1, 'last');
    fwhm = theta_x(index2) - theta_x(index1);
    FWHM(i) = fwhm;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i == 1
        rminx = R_min;
    else
        rmaxx = R_max;
    end
    
    if i == 2
        dynamic_value = rmaxx - rminx;
        %disp(dynamic_value);  
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %finding theta at the dip
    for m=1:1:array_index
        if R_y(m)== R_min
           theta_dip(i) = theta_x(m);
           %disp(theta_dip);
        end
    end
end

for i = 1:s-1
    % Sensitivity
    sens(i) = (theta_dip(i+1)-theta_dip(1))/(n_sm(i+1)-n_sm(1));
    % Detection Accuracy
    da(i) = (theta_dip(i+1)-theta_dip(1))/FWHM(i+1);
    del_theta(i) = (theta_dip(i+1)-theta_dip(1));
    % Quality Factor
    qf(i) = sens(i)/FWHM(i+1);
    % LOD
    LOD(i) = (1/sens(i))*.0573;
end
disp('sensitivity:');
disp(sens);
disp('Detection Accuracy:');
disp(da);
disp('Quality Factor:');
disp(qf);
disp('del theta:');
disp(del_theta);
disp('LOD:');
disp(LOD);
parameters = [n_sm(1), n_sm(2), comb, sens, da, qf];
disp(parameters);
loc = 'C:\Users\Tanjib Ahmed\Pictures\Cancer Detection\variation of layer\structure\structure.csv';
xlswrite(loc, parameters, 'Sheet1', cells);