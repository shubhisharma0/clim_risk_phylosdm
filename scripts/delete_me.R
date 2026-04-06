v <- terra::vect("~/Downloads/shubhi/richness_rarity_birds/richness_rarity_birds.shp")
one <- v[1]
res_guess <- sqrt(expanse(one, unit = "m")) # approximate cell size (m)
r <- rast(ext(v), resolution = res_guess, crs = crs(v))

r2 <- terra::rasterize(v, r, field = "fe")
writeRaster(r2, "~/Downloads/function_endemism_birds.tif", overwrite = TRUE)



v <- terra::vect("~/Downloads/global_facets_mol_20181016/global_facets_mol_20181016.shp")
library(terra)

birds <- v[v$TAXA == "birds", ]
birds

# infer approximate grid resolution from polygon geometry
# (takes a small sample to avoid heavy compute)
g <- geom(birds[1:2000, ], df = TRUE) # vertices
dx <- diff(sort(unique(round(g$x, 6))))
dy <- diff(sort(unique(round(g$y, 6))))
rx <- min(dx[dx > 0], na.rm = TRUE)
ry <- min(dy[dy > 0], na.rm = TRUE)

rx
ry

r_template <- rast(ext(birds), resolution = c(rx, ry), crs = crs(birds))

sr_r <- rasterize(birds, r_template, field = "SR", touches = TRUE)
rsr_r <- rasterize(birds, r_template, field = "SUM_RSR", touches = TRUE)


library(terra)

birds <- v[v$TAXA == "birds", ]

# Sane template (global lon/lat grid)
r_template <- rast(ext(-180, 180, -90, 90), resolution = 0.5, crs = "EPSG:4326")

sr_r <- rasterize(birds, r_template, field = "SR", touches = TRUE)
rsr_r <- rasterize(birds, r_template, field = "SUM_RSR", touches = TRUE)

writeRaster(sr_r, "~/Downloads/birds_SR.tif", overwrite = TRUE)
writeRaster(rsr_r, "~/Downloads/birds_SUM_RSR.tif", overwrite = TRUE)

r <- rast("~/Downloads/birds_SUM_RSR.tif")
world <- terra::vect("~/phylo-sdms2/raw_data/world.shp")
r <- terra::mask(r, world)
