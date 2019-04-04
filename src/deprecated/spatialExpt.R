## Code with experimental spatial analyses

## Spatial eigenvectors modelling ========
library(spdep)
tdwg_shp_raw <- readOGR("~/Dropbox/Projects/2019/palms/data/TDWG/level3/level3.shp")
glob_list <- tdwg_final_glob$LEVEL_3_CO
tdwg_shp_glob <- subset(tdwg_shp_raw, LEVEL_3_CO %in% glob_list)

# Generate neighbourhoods
# `poly2nb` not used as some polygons are not adjacent to anything else (i.e., islands)
k1_mat <- knearneigh(as.matrix(tdwg_final_glob[c("LONG","LAT")]), k = 1, longlat = TRUE)
k1_nb <- knn2nb(k1_mat)
n.comp.nb(k1_nb) # number of disjoint connected subgraphs

k2_mat <- knearneigh(as.matrix(tdwg_final_glob[c("LONG","LAT")]), k = 2, longlat = TRUE)
k2_nb <- knn2nb(k2_mat)
n.comp.nb(k2_nb) # number of disjoint connected subgraphs

k3_mat <- knearneigh(as.matrix(tdwg_final_glob[c("LONG","LAT")]), k = 3, longlat = TRUE)
k3_nb <- knn2nb(k3_mat)
n.comp.nb(k3_nb)

k4_mat <- knearneigh(as.matrix(tdwg_final_glob[c("LONG","LAT")]), k = 4, longlat = TRUE)
k4_nb <- knn2nb(k4_mat)
n.comp.nb(k4_nb)

k5_mat <- knearneigh(as.matrix(tdwg_final_glob[c("LONG","LAT")]), k = 5, longlat = TRUE)
k5_nb <- knn2nb(k5_mat)
n.comp.nb(k5_nb) # 3 subgroups are the different realms

k6_mat <- knearneigh(as.matrix(tdwg_final_glob[c("LONG","LAT")]), k = 6, longlat = TRUE)
k6_nb <- knn2nb(k6_mat)
n.comp.nb(k6_nb) # Number of subgraphs

pdf(file.path(fig.dir, "knn.pdf"), width = 12, height = 6)
par(mfrow = c(2,3))
plot(tdwg_shp_glob, col = "grey", border = "white", main = "k-nearest neighbours = 1")
plot(k1_nb, tdwg_final_glob[c("LONG","LAT")], add = T, col = "red")
plot(tdwg_shp_glob, col = "grey", border = "white", main = "k-nearest neighbours = 2")
plot(k2_nb, tdwg_final_glob[c("LONG","LAT")], add = T, col = "blue")
plot(tdwg_shp_glob, col = "grey", border = "white", main = "k-nearest neighbours = 3")
plot(k3_nb, tdwg_final_glob[c("LONG","LAT")], add = T, col = "gold")
plot(tdwg_shp_glob, col = "grey", border = "white", main = "k-nearest neighbours = 4")
plot(k4_nb, tdwg_final_glob[c("LONG","LAT")], add = T, col = "purple")
plot(tdwg_shp_glob, col = "grey", border = "white", main = "k-nearest neighbours = 5")
plot(k5_nb, tdwg_final_glob[c("LONG","LAT")], add = T, col = "green")
plot(tdwg_shp_glob, col = "grey", border = "white", main = "k-nearest neighbours = 6")
plot(k6_nb, tdwg_final_glob[c("LONG","LAT")], add = T, col = "pink")
dev.off()

# Generate spatial weighting matrices based on great circle distances
k6_nb_dist <- nbdists(k6_nb, coords = as.matrix(tdwg_final_glob[c("LONG","LAT")]), longlat = TRUE) # calculate distances between neighbours
fdist <- lapply(k6_nb_dist, function(x) 1 - x/max(dist(as.matrix(tdwg_final_glob[c("LONG","LAT")])))) # convert into spatial weights
listwgab <- nb2listw(neighbours = k6_nb, glist = fdist, style = "B") # essentially global weighting
mem.gab <- mem(listwgab) # generates Moran's eigenvector maps where correspionding eigenvalues are linearly related to Moran's index of spatial autocorrelation
barplot(attr(mem.gab, "values"), main = "Eigenvalues of spatial weighting matrix")

# Perform Moran's I on each eigenvector
round(attr(mem.gab, "values"), 3)
moranItest <- moran.randtest(x = mem.gab, 
                             list = listwgab, nrepet = 999)
signi <- which(moranItest$pvalue < 0.05) # identify eigenvectors that are positively significant (i.e., positive spatial autocorrelation)
plot(mem.gab[,signi[35:36]], SpORcoords = as.matrix(tdwg_final_glob[c("LONG","LAT")]), nb = k6_nb )

# F
tdwg_final_glob_sp <- cbind(tdwg_final_glob, mem.gab)
glob_med_f <- formula(paste0("logMedFS_scl ~ curr_logMedBS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl", "+", paste0(names(mem.gab)[signi], collapse = "+")))

glob_curr_medBS_spmod <- lm(logMedFS_scl ~ curr_logMedBS_scl + globalPC1_scl + globalPC2_scl + globalPC3_scl + lgm_ens_Tano_scl + lgm_ens_Pano_scl + MEM1, data =tdwg_final_glob_sp, na.action = "na.fail")
#glob_pnat_medBS_spmod <- update(glob_curr_medBS_spmod, ~.-curr_logMedBS_scl + pnat_logMedBS_scl)
summary(glob_curr_medBS_spmod)

# 
sp.cor <- sp.correlogram(k6_nb, tdwg_final_glob$logMedFS_scl, order=5,
                         method="I", randomisation=T)
sp.cor <- sp.correlogram(k6_nb, residuals(glob_curr_medBS_mod), order=5,
                         method="I", randomisation=T)
sp.cor <- sp.correlogram(k6_nb, residuals(glob_curr_medBS_spmod), order=5,
                         method="I", randomisation=T)
par(mfrow = c(1,1))
plot(sp.cor, ylim = c(-1,1))
abline(h = 0.1, col = "red")
# https://cran.r-project.org/web/packages/adespatial/vignettes/tutorial.html#ref-Dray2012

## Defining neighbourhoods
glob_coords <- as.matrix(tdwg_final_glob[c("LONG","LAT")])

# Gabriel triangulation
nb_gbrl <- graph2nb(gabrielneigh(glob_coords), sym = TRUE)

# Delaunay (symmetrical)


# Relative graph
nb_rltv <- graph2nb(relativeneigh(glob_coords), sym = TRUE)

# Create a list of neighbourhood definitions
nb_list <- list(nb_dnear, nb_knear, nb_dlny, nb_gbrl, nb_rltv, nb_soi)

# Plot neighbourhood definitions
pdf(file.path(fig.dir, "neighbourhood.pdf"), width = 10, height = 12)
par(mfrow = c(3,2))
plot(tdwg_shp_glob, col = "grey", border = "white", main = "k-nearest neighbours = 1")
plot(knear1, glob_coords, add = T, col = "#003262")
plot(tdwg_shp_glob, col = "grey", border = "white", main = "Distance")
plot(dnear1, glob_coords, add = T, col = "#3B7EA1")
plot(tdwg_shp_glob, col = "grey", border = "white", main = "Delaunay triangulation")
plot(nb_dlny, glob_coords, add = T, col = "#FDB515")
plot(tdwg_shp_glob, col = "grey", border = "white", main = "Gabriel")
plot(nb_gabr, glob_coords, add = T, col = "#EE1F60")
plot(tdwg_shp_glob, col = "grey", border = "white", main = "Relative")
plot(nb_rltv, glob_coords, add = T, col = "#D9661F")
plot(tdwg_shp_glob, col = "grey", border = "white", main = "Sphere of influence")
plot(nb_soi, glob_coords, add = T, col = "#00A598")
dev.off()
