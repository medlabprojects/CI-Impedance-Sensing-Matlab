function [polygons, areas] = ZsensingFindRegions(RGB)

% Segment to find centers of electrode markers and outlines of modiolus/EA masks
[E, EA, mod] = ZsensingSegmentFrameGPU(RGB);

% Find closest points along modiolus to each electrode
for ii = 1:length(E)
    % find index of mod.poly_smooth point that is closest to electrode marker
    E(ii).p_mod_ind = projPointOnPolygon(E(ii).center, mod.poly_smooth); 

    % find point on boundary of EA between electrode marker and the modiolus
    E(ii).p_EA = intersectPolylines( [E(ii).center; polylinePoint(mod.poly_smooth,E(ii).p_mod_ind)], EA.poly_smooth );
    E(ii).p_EA_ind = projPointOnPolygon(E(ii).p_EA, EA.poly_smooth);
end

% Create polygon for the regions between each pair of electrodes and the modiolus
polygons = cell(length(E)-1,1);
areas = zeros(size(polygons));

for ii = 2:length(E)
    polygons{ii-1} = [polylineSubcurve(EA.poly_smooth,E(ii-1).p_EA_ind,E(ii).p_EA_ind); polylineSubcurve(mod.poly_smooth,E(ii).p_mod_ind,E(ii-1).p_mod_ind)];
    areas(ii-1) = polygonArea(polygons(ii-1)); % area in pixels
end

end