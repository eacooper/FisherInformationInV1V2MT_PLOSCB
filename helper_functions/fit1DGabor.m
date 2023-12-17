function[P, S, X, E] = fit1DGabor( tag, dat, P, S, X, E, n, N )

% number of iterations for fitting
niters = 200;

% load symbolic equations for Gabor and first derivate
define_Gabor_tuning_curve;

% convert disparity and mean responses to row vectors
x_orig  = dat(:,1)';
resp_orig = dat(:,2)';

% interpolate disparities and responses to constrain fits
xi = zeros(1,2*length(x_orig)-1);
xi(1:2:end) = x_orig;
for i = 2 : 2 : length(xi)-1
    xi(i) = (xi(i-1)+xi(i+1))/2.0; % interpolate between each measured disparity
end

ri = interp1(x_orig,resp_orig,xi,'linear');

ri = abs(ri); % all responses should be positive

% fit 1D Gabor to 1D data [Gabor: r0 + a*exp(-(d-d0)^2 / (2*s^2)) * cos(2*pi*f*(d-d0)+phi)]
% parameters: [pedestal amplitude Gaussian_mean Gaussian_sigma frequency phase]

% initialize optimization
options = optimset('Display','off');
minp = NaN;
minerr = inf;

lb = [0   0      -1.75 0 0 -2*pi]; % lower bounds
ub = [500 500    1.75  5 4.5   2*pi]; % upper bounds

for k = 1 : niters

    % random seed
    seed = [uniform_sample_in_range(lb(1),ub(1)) ...
        uniform_sample_in_range(lb(2),ub(2)) ...
        uniform_sample_in_range(lb(3),ub(3)) ...
        uniform_sample_in_range(lb(4),ub(4)) ...
        uniform_sample_in_range(lb(5),ub(5)) ...
        uniform_sample_in_range(lb(6),ub(6))];

    p = fmincon( 'errfun1DGabor', seed, [], [], [], [], lb, ub, [], options, [xi' ri'] );
    err = errfun1DGabor( p, [xi' ri'] );

    % if error is better than best Gaussian, go with it
    if( err < minerr )

        fprintf( '%d, Gabor iteration %d: found better fit: %.2f %.2f %.2f %.2f %.2f %.2f (%f) --> %.2f %.2f %.2f %.2f %.2f %.2f (%f)\n', ...
            n, k, minp, minerr, p, err );

        minp = p;
        minerr = err;

    end

end

% calculate R2 of best Gabor parameters

% calculate fitted spike rates for each tested disparity
TF = gabor_tuning_curve(x_orig,minp);

% then do R2
SSR = sum((resp_orig - TF).^2);
TSS = sum((resp_orig - mean(resp_orig)).^2);
r2   = 1 - (SSR/TSS);

fprintf('Final R2 = %.2f \n', r2)

% evaluate best fit Gabor and it's first derivative over a finer lattice
xg1 = [-2 : 0.01 : 2]; %fixed spatial range
g1 = gabor_tuning_curve(xg1,minp);
dg1 = first_deriv_of_GTC(xg1,minp);

% visualize data, model fit, and FI

% data and fit
figure(1);
cla; plot(xg1,g1,'b'); hold on;  plot(xi,ri,'ro');  hold off;
xlabel('horizontal disparity');
ylabel('model response');
set( gca, 'Ylim', [0 1.2*max(ri)] );
grid on;
h = title( sprintf('neuron = %d | err = %.3f', n, minerr) );

drawnow;

% Store across all neurons
if(n == 1) % first neuron, initialize
    P  = zeros( N, 6 ); % model parameters
    S(length(ri)).s = []; % spike count used for fitting
    X(length(xi)).x = []; % horizontal disparity used for fittimh
    E  = zeros( N, 1 );
end

P(n,:)  = minp;
S(n).s  = ri;
X(n).x  = xi;
E(n)    = minerr;



