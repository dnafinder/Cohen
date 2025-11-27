function [k, ci, stats] = kappa(x,varargin)
%KAPPA Compute Cohen's kappa coefficient (unweighted or weighted).
%
% Cohen's kappa coefficient is a statistical measure of inter-rater
% reliability. It is generally thought to be a more robust measure than
% simple percent agreement calculation since kappa takes into account the
% agreement occurring by chance.
%
% Kappa provides a measure of the degree to which two judges, A and B,
% concur in their respective sortings of N items into k mutually exclusive
% categories. A 'judge' in this context can be an individual human being, a
% set of individuals who sort the N items collectively, or some non-human
% agency, such as a computer program or diagnostic test, that performs a
% sorting on the basis of specified criteria.
%
% The original and simplest version of kappa is the unweighted kappa
% coefficient introduced by J. Cohen in 1960. When the categories are
% merely nominal, Cohen's simple unweighted coefficient is the only form of
% kappa that can meaningfully be used. If the categories are ordinal and if
% category 2 represents more of something than category 1, category 3 more
% than category 2, and so on, then it is potentially meaningful to take
% this ordering into account, weighting each cell of the matrix in
% accordance with how near it is to the cell on the diagonal that includes
% the absolutely concordant items.
%
% This function can compute unweighted, linearly weighted or quadratically
% weighted kappa, or use a custom weight matrix.
%
% Syntax:
%   kappa(X)
%   kappa(X,W)
%   kappa(X,W,ALPHA)
%   kappa(X,W,ALPHA,DISPLAY)
%   [K,CI] = kappa(...)
%   [K,CI,STATS] = kappa(...)
%
% Inputs:
%   X        - square data matrix (M-by-M) of nonnegative integer
%              frequencies.
%   W        - weight specification:
%                * scalar: 0 = unweighted (default)
%                          1 = linear weights
%                          2 = quadratic weights
%                * matrix (M-by-M): custom weight matrix for M categories.
%   ALPHA    - significance level for confidence interval (default 0.05).
%   DISPLAY  - logical flag (true/false) to control textual output
%              (default true).
%
% Outputs:
%   K        - Cohen's kappa.
%   CI       - 1-by-2 vector with confidence interval for K.
%   STATS    - structure with additional results:
%              .po                     observed agreement
%              .pe                     expected agreement by chance
%              .trueAgreement          po - pe
%              .residualAgreement      1 - pe
%              .kappa                  Cohen's kappa
%              .se                     standard error of kappa
%              .ci                     confidence interval for kappa
%              .alpha                  significance level
%              .z                      z statistic
%              .p                      p-value
%              .kappaMax               maximum possible kappa
%              .kappaRatio             kappa / kappaMax
%              .weightsType            'unweighted' | 'linear' | 'quadratic' | 'custom'
%              .weightMatrix           weight matrix used
%              .n                      total number of observations
%              .m                      number of categories
%              .maxObservableAgreement maximum observable agreement (pom)
%              .landisKochClass        qualitative classification
%
% Example:
%
%   x = [88 14 18; 10 40 10; 2 6 12];
%   kappa(x);
%
%   x = [88 14 18; 10 40 10; 2 6 12];
%   [k, ci, stats] = kappa(x,1,0.05,false); % linear weights, no display
%
% Author:   Giuseppe Cardillo
% Email:    giuseppe.cardillo.75@gmail.com
% GitHub:   https://github.com/dnafinder/Cohen
%
% To cite this file, this would be an appropriate format:
%   Cardillo G. (2007) Cohen's kappa: compute the Cohen's kappa ratio
%   on a square matrix. Available from:
%   https://github.com/dnafinder/Cohen
%
% -------------------------------------------------------------------------

% Input parsing and validation
p = inputParser;
p.FunctionName = 'kappa';

addRequired(p,'x', @(X) validateattributes(X,{'numeric'},...
    {'square','nonempty','integer','real','finite','nonnan','nonnegative'}));

addOptional(p,'w',0, @(W) isnumeric(W) && isreal(W) && ...
    all(isfinite(W(:))) && ~isempty(W));

addOptional(p,'alpha',0.05, @(a) validateattributes(a,{'numeric'},...
    {'scalar','real','finite','nonnan','>',0,'<',1}));

addOptional(p,'displayopt',true, @(d) islogical(d) && isscalar(d));

parse(p,x,varargin{:});
x          = p.Results.x;
w          = p.Results.w;
alpha      = p.Results.alpha;
displayopt = p.Results.displayopt;
clear p

% Basic size info
m = size(x,1);
n = sum(x(:));

if n <= 0
    error('kappa:NoData','The input matrix X contains no positive counts.');
end

% Build weight matrix
if isscalar(w)
    if ~ismember(w,[0 1 2])
        error('kappa:InvalidWeightCode',...
            'When W is scalar, it must be 0 (unweighted), 1 (linear), or 2 (quadratic).');
    end
    switch w
        case 0
            f = eye(m); % unweighted
            weightsType = 'unweighted';
        case 1
            J = repmat(1:m,m,1);
            I = flipud(rot90(J));
            f = 1 - abs(I-J)./(m-1); % linear weights
            weightsType = 'linear';
        case 2
            J = repmat(1:m,m,1);
            I = flipud(rot90(J));
            f = 1 - ((I-J)./(m-1)).^2; % quadratic weights
            weightsType = 'quadratic';
    end
else
    % Custom weight matrix
    if ~ismatrix(w) || any(size(w) ~= [m m])
        error('kappa:InvalidWeightMatrix',...
            'When W is a matrix, it must be M-by-M where M = size(X,1).');
    end
    f = w;
    weightsType = 'custom';
end

% Proportions
N = x;              % counts
P = N ./ n;         % proportions

% Row/column marginals
r = sum(P,2);       % rows sum (M-by-1)
s = sum(P,1);       % columns sum (1-by-M)

if any(r==0) || any(s==0)
    warning('kappa:ZeroMarginal',...
        'Some categories have zero marginal frequencies; interpretation of kappa may be problematic.');
end

% Expected proportions under independence
Ex = r * s;         % M-by-M

% Observed and expected agreement (weighted)
po = sum(sum(P .* f));  % observed agreement
pe = sum(sum(Ex .* f)); % expected agreement

denom = 1 - pe;

% Initialize outputs
k         = NaN;
ci        = [NaN NaN];
sek       = NaN;
z         = NaN;
pval      = NaN;
kmax      = NaN;
ratio     = NaN;
pom       = NaN;
trueAgree = po - pe;
residAgree= 1 - pe;

if abs(denom) < eps
    warning('kappa:ZeroDenominator',...
        '1 - pe is numerically zero; kappa cannot be computed reliably.');
else
    % Cohen's kappa
    k = (po - pe) / denom;

    % Maximum observable agreement
    pom = sum(min([r'; s]));
    % Maximum possible kappa, given the observed marginal frequencies
    kmax = (pom - pe) / denom;
    % Observed as proportion of maximum possible
    ratio = k / kmax;

    % Standard error
    sek = sqrt((po * (1 - po)) / (n * denom^2));

    % Confidence interval using normal approximation
    zcrit = abs(-sqrt(2)*erfcinv(alpha)); % two-sided
    ci    = k + [-1 1] * (zcrit * sek);

    % z statistic and p-value
    z    = k / sek;
    pval = (1 - 0.5*erfc(-abs(z)/sqrt(2))) * 2;
end

% Landis & Koch classification
if isnan(k)
    landis = 'Not defined';
elseif k < 0
    landis = 'Poor agreement';
elseif k <= 0.20
    landis = 'Slight agreement';
elseif k <= 0.40
    landis = 'Fair agreement';
elseif k <= 0.60
    landis = 'Moderate agreement';
elseif k <= 0.80
    landis = 'Substantial agreement';
else
    landis = 'Perfect agreement';
end

% Build stats struct
stats = struct();
stats.po                     = po;
stats.pe                     = pe;
stats.trueAgreement          = trueAgree;
stats.residualAgreement      = residAgree;
stats.kappa                  = k;
stats.se                     = sek;
stats.ci                     = ci;
stats.alpha                  = alpha;
stats.z                      = z;
stats.p                      = pval;
stats.kappaMax               = kmax;
stats.kappaRatio             = ratio;
stats.weightsType            = weightsType;
stats.weightMatrix           = f;
stats.n                      = n;
stats.m                      = m;
stats.maxObservableAgreement = pom;
stats.landisKochClass        = landis;

% Textual output (if requested)
if displayopt
    tr = repmat('-',1,80);
    switch weightsType
        case 'unweighted'
            disp('UNWEIGHTED COHEN''S KAPPA');
        case 'linear'
            disp('LINEAR WEIGHTED COHEN''S KAPPA');
        case 'quadratic'
            disp('QUADRATIC WEIGHTED COHEN''S KAPPA');
        otherwise
            disp('CUSTOM WEIGHTED COHEN''S KAPPA');
    end
    disp(tr)

    fprintf('Observed agreement (po) = %0.4f\n',po);
    fprintf('Random agreement (pe) = %0.4f\n',pe);
    fprintf('Agreement due to true concordance (po-pe) = %0.4f\n',trueAgree);
    fprintf('Residual not random agreement (1-pe) = %0.4f\n',residAgree);
    fprintf('Cohen''s kappa = %0.4f\n',k);
    fprintf('kappa error = %0.4f\n',sek);
    fprintf('kappa C.I. (alpha = %0.4f) = %0.4f     %0.4f\n',alpha,ci);
    fprintf('Maximum possible kappa, given the observed marginal frequencies = %0.4f\n',kmax);
    fprintf('k observed as proportion of maximum possible = %0.4f\n',ratio);
    disp(landis);
    fprintf('z (k/kappa error) = %0.4f    p = %0.4f\n',z,pval);
    if isnan(pval)
        disp('Null hypothesis test not available: kappa cannot be computed reliably.');
    elseif pval < alpha
        disp('Reject null hypothesis: observed agreement is not accidental');
    else
        disp('Accept null hypothesis: observed agreement is accidental');
    end
end

end
