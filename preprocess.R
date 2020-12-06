## Preprocessing

# load libraries
library("sf")
library("tidyr")
library("dplyr")

#multilevel modeling (test)
dat2 <- readRDS("data_for_COVIDmodel.rds")
names(dat2)

dat2 <- st_transform(dat2, 26912)

dat2 <- dat2[!is.na(dat2$TAVG), ]
dat2 <- dat2[!is.na(dat2$ln_empden_000_qtmi), ]
dat2 <- dat2[!is.na(dat2$ln_hhsize_qtmi), ]
dat2 <- dat2[!is.na(dat2$intden_hami), ]

signal.zero <- which(table(dat2$SIGNAL) == 0)

dat2$SIGNAL2 <- as.numeric(as.factor(dat2$SIGNAL))

dat3 <- dat2 %>% distinct(SIGNAL)

coords <- st_coordinates(dat3)
k1dist <- knn2nb(knearneigh(coords, k=1), sym = TRUE)
all.linked <- max(unlist(nbdists(k1dist, coords)))

nb <- dnearneigh(coords, d1 = 0, d2 = 1 * all.linked)

nb <- graph2nb(relativeneigh(coords), sym = TRUE, row.names = dat3$SIGNAL)

plot(nb, coords)

save(dat2, dat3, nb, file = "newdata.RData")
