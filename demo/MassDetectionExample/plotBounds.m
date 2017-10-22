function plotBounds(prior, posterior, labels)
%   PLOTBOUNDS(PRIOR, POSTERIOR) plots errorbars for PRIOR and POSTERIOR
%   vectors of size n-by-2 of minimum and maximum values of n number of 
%   observations or variables.
%
%   PLOTBOUNDS(PRIOR, POSTERIOR, LABELS) plots errorbars for PRIOR and
%   POSTERIOR as well as LABELS along the x-axis. LABELS is a cell array
%   of size n-by-1 containing descriptive names for the n number of 
%   observation or variables plotted.
%

if nargin < 3
    labels = 0;
end

for i = 1:size(prior,1)
    midPrior = mean(prior(i,:));
    h = errorbar(i, midPrior, prior(i,2)-midPrior, 'b');
    h.LineWidth = 2;
    hold on
    midPosterior = mean(posterior(i,:));
    h2 = errorbar(i, midPosterior, posterior(i,2)-midPosterior, 'r');
    h2.LineWidth = 4;
end

ax = gca;
ax.XTick = 1:size(prior,1);
legend('Prior Bounds', 'Posterior Bounds', 'Location', 'best');
grid on
if  iscell(labels)
    ax.XTickLabel = labels;
end
ax.FontSize = 13;

end
