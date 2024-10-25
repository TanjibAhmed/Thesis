clc;
clear all;
close all;

%Defining Const.
neg = -1;
imgnum = 1*i;
c = 3e08;

%Number of value as refractive index of sensing medium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%n_s = [ 1.3739,  1.378,  1.3783];
n_s = [ 1.357, 1.378];
%n_s = [1.368, 1.392, 1.376, 1.390, 1.381, 1.395, 1.36, 1.38, 1.387, 1.401, 1.385, 1.399];%all in one
%n_s = [1.368, 1.392];%cervical cancer
%n_sm = [1.376, 1.390];%blood cancer
%n_sm = [1.381, 1.395];%adrenal gland cancer-----graph issue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%s = length(n_sm);
s=2;
FWHM = zeros(1,s);
%Initial variables
start_angle = 60;
start_angle=deg2rad(start_angle);
end_angle = 90; 
end_angle=deg2rad(end_angle);
anglestep = .05; 
anglestep=deg2rad(anglestep);

% Wavelength
lambda = 633e-9;

k0 = (2 * pi) / (lambda);                  
w = k0 * c;

%NanoComposite layer
n_d = 2.2;
e_d = n_d^2;
n_ag = 0.05626+4.2776i; 
e_ag = n_ag^2;

V = 3*10^-2; %range: 1 to 5

e_NC = ((2*e_d*V*(e_ag-e_d)) + (e_d*(e_ag+2*e_d)))/((2*e_d)+e_ag+(V*(e_d-e_ag)));
n_NC = e_NC^(1/2);


%Refractive index of each layer
n_prism = 1.51; %BK7                                    %Structure 1:1/2
n_first = 0.056206 + 4.2776i; %Ag                         %Structure 2:1/2/3
n_second = 2.58; %GaP                                      %Structure 3:1/2/3/4
n_third = n_NC; %Ag-CeO2 nanocomposite                    %Structure 4:1/2/3/4/5
n_fourth = 2.67; %PbTiO3

n = [n_prism, n_first, n_second, n_third, n_fourth, 0];

%Number of layers including prism and sensing medium
L=length(n);
theta_dip = zeros(1,s);

ind = 0;
width = [];
sensitivity = [];
%q_factor = [];
%d_accuracy = [];
fwhm = [];
dip = [];
Rmin = [];
www = [];
Q =[];
D = [];
for o = 1:6
    if o==1
        n_sm = n_s(1:2);
        disp('cervical cancer:');
    end
    if o==2
        n_sm = n_s(3:4);
        disp('blood cancer:');
    end
    if o==3
        n_sm = n_s(5:6);
        disp('adrenal gland cancer:');
    end
    if o==4
        n_sm = n_s(7:8);
        disp('Skin Cancer:');
    end
    if o==5
        n_sm = n_s(9:10);
        disp('Breast cancer(MDA-7):');
    end
    if o==6
        n_sm = n_s(11:12);
        disp('Breast cancer(MDA-231):');
    end

    for d_first = 50:50
        for d_second = 3.4:.1:3.4
            for d_third  = .5:.1:.5                    
                for d_fourth = .1:.1:.1
                    ind = ind+1;
                    ww = [d_third, d_fourth];
                    www = [www; ww];
                    %disp(ind);
                    d = [d_first*1e-9, d_second*1e-9, d_third*1e-9, d_fourth*1e-9];
                    width = [width; d];
                    for i = 1:s
                        array_index = 0;
                        %n_sensinglayer = input ('Enter the refractive index of sensing medium: ');
                        n(L) = n_sm(i);

                        for theta = start_angle:anglestep:end_angle
                            array_index = array_index + 1;
                            b = zeros(1,L);
                            B = zeros(1,L-2);

                            % Calculation of the variables
                            for k = 1:L        
                                b(k) = sqrt(n(k)^2-(n(1)^2 *(sin(theta))^2)) / n(k)^2;
                                if (1<k) && (k<L)
                                    B(k-1) = d(k-1)* k0 * sqrt(n(k)^2-(n(1)^2 *(sin(theta))^2));
                                end
                            end

                            % Calculation of the Matrix values
                            for j = 1:L-2
                                A1 = cos(B(j));
                                A2 = (neg * sin(B(j)) * (imgnum)) / b(j+1);
                                A3 = ((neg * imgnum * sin(B(j))) * b(j+1));
                                A4 = cos(B(j));
                                if j == 1
                                    A = [A1, A2; A3, A4];
                                end
                                if j > 1
                                    A_m = [A1, A2; A3, A4];
                                    A = A*A_m;
                                end
                            end


                            % Calculate the reflection value
                            Rnum = ((A(1,1) + A(1,2) * b(L)) * b(1)) - (A(2,1) + A(2,2) * b(L));
                            Rdem = ((A(1,1) + A(1,2) * b(L)) * b(1)) + (A(2,1) + A(2,2) * b(L));
                            R = (Rnum/Rdem);
                            phase = R;
                            R = sqrt(real(R)^2 + imag(R)^2);
                            theta_x(array_index)=theta; 
                            R_y(array_index) = R^2;
                            E(array_index) = 1/R^2;
                            P(array_index) = angle(phase);
                        end
                        
                        figure(o)
                        plot(theta_x(1:1:array_index)*(180/pi),R_y(1:1:array_index), 'LineWidth', 2);
                        xlabel('Incidence Angle (in degree)');
                        ylabel('Reflectance(a.u.)');
                        hold on;

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        R_min = min(R_y(1:1:array_index));
                        Rmn(i) = R_min; 
                        R_max = max(R_y(1:1:array_index));
                        %Calculation of the FWHM

                        R_halfmax = (R_max+R_min)/2;
                        index1 = find(R_y(1:1:array_index) <= R_halfmax, 1, 'first');
                        index2 = find(R_y(1:1:array_index) <= R_halfmax, 1, 'last');
                        fw = theta_x(index2) - theta_x(index1);
                        FWHM(i) = fw*(180/pi);

                        %finding theta at the dip
                        for m=1:1:array_index
                            if R_y(m)== R_min
                               theta_dip(i) = theta_x(m)*(180/pi);
                               %disp(theta_dip);
                            end
                        end
                    end
                    sens = zeros(1,s-1);
                    da = zeros(1,s-1);
                    qf = zeros(1,s-1);

                    for i = 1:s-1
                        % Sensitivity
                        sens(i) = (theta_dip(i+1)-theta_dip(1))/(n_sm(i+1)-n_sm(1));
                        % Detection Accuracy
                        da(i) = (theta_dip(i+1)-theta_dip(1))/FWHM(i+1);
                        % Quality Factor
                        qf(i) = sens(i)/FWHM(i+1);
                    end

                    fwhm = [fwhm; FWHM];
                    dip = [dip; theta_dip];
                    Rmin = [Rmin; Rmn];
                    sensitivity = [sensitivity, sens];
                    Q = [Q, qf];
                    D = [D, da];
                end
                %xlswrite(loc, transpose(max(sensitivity)), 'Sheet1', 'A1');
            end
        end 
    end
    loc = 'C:\Users\Tanjib Ahmed\Pictures\BG\layer\reflectance.csv';
    xlswrite(loc, transpose(R_y), 'Sheet1', 'A2');
    %disp(sensitivity);
    disp(sensitivity);
    %disp(Q);
    %disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    %disp(D);
    disp(Rmin);
    Rmin = [];
    sensitivity = [];
    Q =[];
    D = [];
end
