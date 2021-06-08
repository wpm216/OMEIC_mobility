function sampleTable = bootstrapRateTable(rateTable, rowZMaxs)
% Takes a random sample from one of our rate tables.
% Average rate tables are of size 4x5.

% rateTable = table to sample.
% rowZMaxs = last index in each z-slice that isn't NaN.
%           because of how we ran the calculations, it's only a function of
%           the row (concentration), not the column (distance)

sampleTable = zeros(4,5);
for i=1:4
    for j=1:5
        sampleTable(i, j) = mean(datasample(rateTable(i, j, 1:rowZMaxs(i)), ...
                                    rowZMaxs(i), 'Replace', true), 'omitnan');
    end
end

end