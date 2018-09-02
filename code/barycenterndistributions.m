function baryCenterPoints=barycenterndistributions(allpoints, sampleSizes)

distanceMatrix=computesquareddistances(allpoints, true, 0);

if isOneInfinite(distanceMatrix)
    return
end

baryCenterPoints= computebarycenterdistributions(sampleSizes, allpoints, distanceMatrix);

end