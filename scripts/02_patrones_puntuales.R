# ==============================================================================
# 02_patrones_puntuales.R
# ==============================================================================

# 1. CONFIGURACIÓN Y CARGA DE DATOS
# Cargamos el dataset limpio generado por el Script 01
dataset_final <- readRDS("final/data/processed/dataset_final.rds")

# Librerías
library(spatstat)
library(sf)
library(ggplot2)
library(tidyverse)
library(kableExtra)
library(knitr)
library(ggthemes)
library(scales)

set.seed(2026) 

# 2. PREPARACIÓN DEL OBJETO PPP (spatstat) 
# Aseguramos CRS métrico
if(st_crs(dataset_final)$epsg != 25831) {
  dataset_final <- st_transform(dataset_final, 25831)
}

# Ventana de observación: límite municipal de Barcelona
bcn_muni_utm <- mapSpain::esp_get_munic() %>%
  filter(cpro == "08", cmun == "019") %>%
  st_transform(25831)

W <- as.owin(bcn_muni_utm)

ppp_cultura <- dataset_final %>%
  filter(categoria %in% c("Cultura", "Atracción")) %>%  # Solo teatros, museos, centros de arte
  st_transform(25831)

# Convertir a patrón puntual (ppp)
coords <- st_coordinates(ppp_cultura)
ppp_obj <- as.ppp(coords, W = W)

# 3. PRUEBAS DE ALEATORIEDAD ESPACIAL (CSR) 

# 3.1 Quadrant Test (3x3 celdas para evitar conteos bajos)
quad_res <- quadrat.test(ppp_obj, nx = 3, ny = 3)


# Función G con envelope (99 simulaciones, nrank=2 -> ~96% IC)
env_g <- envelope(ppp_obj, Gest, nsim = 99, nrank = 2, verbose = FALSE)

# Funciones K y L de Ripley
K_res <- Kest(ppp_obj)
L_res <- Lest(ppp_obj)


stats_csr <- tibble(
  Metodo = c("Quadrant Test (X²)", "Quadrant Test (p-value)",
             "Max G(r) observado", "G teórica CSR (100m)",
             "K(r) a 1000m", "L(r) a 1000m"),
  Valor = c(
    unname(quad_res$statistic),
    quad_res$p.value,
    max(env_g$obs, na.rm = TRUE),
    approx(env_g$r, env_g$theo, xout = 100)$y,
    approx(K_res$r, K_res$iso, xout = 1000)$y,
    approx(L_res$r, L_res$iso, xout = 1000)$y
  )
)

write.csv(stats_csr, "final/outputs/tabla_03_stats_csr.csv", row.names = FALSE)

# OPTIMIZACIÓN DE INTENSIDAD KERNEL

# bw.diggle minimiza el error cuadrático medio
h_optimo <- bw.diggle(ppp_obj)
h_optimo
# Estimación de densidad
dens <- density.ppp(ppp_obj, sigma = h_optimo)

# 4. EXPORTACIÓN DE FIGURAS 

# 4.1 Quadrant Test 

# Asegurar que los cuadrantes tengan el CRS correcto al crearlos
quad_count <- quadratcount(ppp_obj, nx = 5, ny = 5)

quad_sf <- st_as_sf(st_as_sfc(as.tess(quad_count))) %>%
  st_set_crs(25831) %>%
  mutate(
    conteo = as.vector(quad_count),
    lon = st_coordinates(st_centroid(.))[,1],
    lat = st_coordinates(st_centroid(.))[,2]
  )

# Generar el gráfico con ggplot2
p_quad <- ggplot() +
  # Dibujar el contorno del municipio de Barcelona de fondo
  geom_sf(data = bcn_muni_utm, fill = "#fcfcfc", color = "grey80", linewidth = 0.5) +
  
  # Dibujar la rejilla de los cuadrantes
  geom_sf(data = quad_sf, fill = NA, color = "grey50", linetype = "dashed", linewidth = 0.4) +
  
  # Dibujar los puntos del patrón puntual
  geom_sf(data = st_as_sf(ppp_obj) %>% st_set_crs(25831), color = "red", size = 0.6, alpha = 0.6) +
  
  # Añadir los números azules con el conteo en el centro de cada celda
  geom_text(data = quad_sf, aes(x = lon, y = lat, label = conteo), 
            color = "blue", fontface = "bold", size = 4) +
  labs(
    title = "Conteo por cuadrantes (5x5)",
    subtitle = "Número de puntos observados por celda",
    caption = "Elaboración propia. Basado en el patrón de puntos corregido"
  ) +
  theme_fivethirtyeight() +
  theme(
    plot.background = element_rect(fill = "white", color = NA), 
    legend.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(
      size = 16,
      face = "bold",
      color = "#2c3e50",
      hjust = 0.5,
      margin = margin(b = 12)
    ),
    plot.subtitle = element_text(
      size = 11,
      color = "#7f8c8d",
      hjust = 0.5,
      margin = margin(b = 18)
    ),
    plot.caption = element_text(
      size = 9,               
      color = "#7f8c8d",     
      hjust = 0.5,            
      margin = margin(t = 15) 
    ),
    axis.text.x = element_text(
      size = 10,
      color = "#34495e",
      margin = margin(t = 12)
    ),
    axis.title.x = element_text(
      size = 11,
      color = "#34495e",
      margin = margin(t = 18),
    ),
    axis.title.y = element_text(
      size = 11,
      color = "#34495e",
      margin = margin(r = 18)
    ),
    axis.text.y = element_text(
      size = 10, 
      color = "#34495e", 
      margin = margin(r = 10)
    ),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10.5),
    plot.margin = margin(30, 20, 20, 20)
  )

ggsave("final\\outputs\\figura_02_quadrant_test.tiff", 
       plot = p_quad, 
       width = 10, 
       height = 7, 
       dpi = 300, 
       bg = "white",
       compression = "lzw")

# 4.2 Figura 3: Función G con entorno de confianza
# Convertir el objeto envelope a dataframe estándar
env_g_df <- as.data.frame(env_g)

# Construir el gráfico
p_funcion_g <- ggplot(env_g_df, aes(x = r)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = "Entorno de simulación"), alpha = 0.3) +
  geom_line(aes(y = hi, color = "G_hi(r)"), linetype = "dotdash", linewidth = 0.6) +
  geom_line(aes(y = lo, color = "G_lo(r)"), linetype = "dotdash", linewidth = 0.6) +
  geom_line(aes(y = theo, color = "G_theo(r)"), linetype = "dashed", linewidth = 0.8) +
  geom_line(aes(y = obs, color = "obs"), linewidth = 1) +
  scale_color_manual(
    name = "Referencias",
    values = c(
      "obs" = "black",
      "G_theo(r)" = "#e74c3c",
      "G_hi(r)" = "grey40",
      "G_lo(r)" = "grey40"
    ),
    labels = c(
      "G_hi(r)" = expression(hat(G)[hi](r)),
      "G_theo(r)" = expression(G[theo](r)),
      "G_lo(r)" = expression(hat(G)[lo](r)),
      "obs" = expression(obs)
    )
  ) +
  scale_fill_manual(
    name = "Área",
    values = c("Entorno de simulación" = "#7647F5"),
    labels = c("Entorno de simulación" = "Entorno de confianza (96%)")
  ) +
  scale_x_continuous(breaks = seq(0, 200, by = 50)) +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.4)) +
  coord_cartesian(xlim = c(0, 220), ylim = c(0, 0.9)) +
  labs(
    title = "Test CSR: Función G con entorno de confianza (96%)",
    subtitle = "Valores observados superiores a la línea teórica confirman un patrón agrupado",
    x = "Distancia r (metros)",
    y = "G(r)",
    caption = "Elaboración propia. Basado en 99 simulaciones de Monte Carlo"
  ) +
  theme_fivethirtyeight() +
  theme(
    plot.background = element_rect(fill = "white", color = NA), 
    legend.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(
      size = 16,
      face = "bold",
      color = "#2c3e50",
      hjust = 0.5,
      margin = margin(b = 12)
    ),
    plot.subtitle = element_text(
      size = 11,
      color = "#7f8c8d",
      hjust = 0.5,
      margin = margin(b = 18)
    ),
    plot.caption = element_text(
      size = 9,               
      color = "#7f8c8d",     
      hjust = 0.5,            
      margin = margin(t = 15) 
    ),
    axis.text.x = element_text(
      size = 10,
      color = "#34495e",
      margin = margin(t = 12)
    ),
    axis.title.x = element_text(
      size = 11,
      color = "#34495e",
      margin = margin(t = 18),
    ),
    axis.title.y = element_text(
      size = 11,
      color = "#34495e",
      margin = margin(r = 18)
    ),
    axis.text.y = element_text(
      size = 10, 
      color = "#34495e", 
      margin = margin(r = 10)
    ),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10.5),
    plot.margin = margin(30, 20, 20, 20)
  )

ggsave("final\\outputs\\figura_03_funcion_g.tiff", 
       plot = p_funcion_g, 
       width = 10, 
       height = 7, 
       dpi = 300, 
       bg = "white",
       compression = "lzw")


# 4.3 Figura 4: K y L de Ripley
# Preparar los datos de la Función K

K_df <- as.data.frame(K_res) %>%
  select(r, theo, iso, trans, border) %>% 
  pivot_longer(cols = c(theo, iso, trans, border), names_to = "Variable", values_to = "Valor") %>%
  mutate(
    Funcion = "Función K de Ripley",
    Leyenda = case_when(
      Variable == "iso"    ~ "hat(K)[iso](r)",
      Variable == "trans"  ~ "hat(K)[trans](r)",
      Variable == "border" ~ "hat(K)[bord](r)",
      Variable == "theo"   ~ "K[pois](r)"
    )
  )

# Preparar los datos de la Función L
L_df <- as.data.frame(L_res) %>%
  select(r, theo, iso, trans, border) %>%
  pivot_longer(cols = c(theo, iso, trans, border), names_to = "Variable", values_to = "Valor") %>%
  mutate(
    Funcion = "Función L de Ripley",
    Leyenda = case_when(
      Variable == "iso"    ~ "hat(L)[iso](r)",
      Variable == "trans"  ~ "hat(L)[trans](r)",
      Variable == "border" ~ "hat(L)[bord](r)",
      Variable == "theo"   ~ "L[pois](r)"
    )
  )

# Combinar ambos conjuntos de datos
ripley_df <- bind_rows(K_df, L_df)

# Crear el gráfico
p_ripley <- ggplot(ripley_df, aes(x = r, y = Valor, color = Leyenda, linetype = Leyenda)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~Funcion, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA)) + 
  
  scale_color_manual(
    name = "Referencias",
    breaks = c("hat(K)[iso](r)", "hat(K)[trans](r)", "hat(K)[bord](r)", "K[pois](r)",
               "hat(L)[iso](r)", "hat(L)[trans](r)", "hat(L)[bord](r)", "L[pois](r)"),
    values = c(
      "hat(K)[iso](r)"   = "black",   "hat(L)[iso](r)"   = "black",
      "hat(K)[trans](r)" = "#e74c3c", "hat(L)[trans](r)" = "#e74c3c",
      "hat(K)[bord](r)"  = "#2ecc71", "hat(L)[bord](r)"  = "#2ecc71",
      "K[pois](r)"       = "#3498db", "L[pois](r)"       = "#3498db"
    ),
    labels = c(
      "hat(K)[iso](r)"   = expression(hat(K)[iso](r)),
      "hat(K)[trans](r)" = expression(hat(K)[trans](r)),
      "hat(K)[bord](r)"  = expression(hat(K)[bord](r)),
      "K[pois](r)"       = expression(K[pois](r)),
      "hat(L)[iso](r)"   = expression(hat(L)[iso](r)),
      "hat(L)[trans](r)" = expression(hat(L)[trans](r)),
      "hat(L)[bord](r)"  = expression(hat(L)[bord](r)),
      "L[pois](r)"       = expression(L[pois](r))
    )
  ) +
  
  scale_linetype_manual(
    name = "Referencias",
    breaks = c("hat(K)[iso](r)", "hat(K)[trans](r)", "hat(K)[bord](r)", "K[pois](r)",
               "hat(L)[iso](r)", "hat(L)[trans](r)", "hat(L)[bord](r)", "L[pois](r)"),
    values = c(
      "hat(K)[iso](r)"   = "solid",   "hat(L)[iso](r)"   = "solid",
      "hat(K)[trans](r)" = "dashed",  "hat(L)[trans](r)" = "dashed",
      "hat(K)[bord](r)"  = "dotted",  "hat(L)[bord](r)"  = "dotted",
      "K[pois](r)"       = "dotdash", "L[pois](r)"       = "dotdash"
    ),
    labels = c(
      "hat(K)[iso](r)"   = expression(hat(K)[iso](r)),
      "hat(K)[trans](r)" = expression(hat(K)[trans](r)),
      "hat(K)[bord](r)"  = expression(hat(K)[bord](r)),
      "K[pois](r)"       = expression(K[pois](r)),
      "hat(L)[iso](r)"   = expression(hat(L)[iso](r)),
      "hat(L)[trans](r)" = expression(hat(L)[trans](r)),
      "hat(L)[bord](r)"  = expression(hat(L)[bord](r)),
      "L[pois](r)"       = expression(L[pois](r))
    )
  ) +
  
  labs(
    title = "Funciones K y L de Ripley",
    subtitle = "Análisis multiescalar con correcciones de borde (Isotrópica, Transmisión y Borde)",
    x = "Distancia r (m)",
    y = "Valor de la función",
    caption = "Elaboración propia. Las funciones observadas por encima de la teórica (pois) indican agrupamiento."
  ) +
  theme_fivethirtyeight() +
  theme(
    plot.background = element_rect(fill = "white", color = NA), 
    legend.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(
      size = 16,
      face = "bold",
      color = "#2c3e50",
      hjust = 0.5,
      margin = margin(b = 12)
    ),
    plot.subtitle = element_text(
      size = 11,
      color = "#7f8c8d",
      hjust = 0.5,
      margin = margin(b = 18)
    ),
    plot.caption = element_text(
      size = 9,               
      color = "#7f8c8d",     
      hjust = 0.5,            
      margin = margin(t = 15) 
    ),
    axis.text.x = element_text(
      size = 10,
      color = "#34495e",
      margin = margin(t = 12)
    ),
    axis.title.x = element_text(
      size = 11,
      color = "#34495e",
      margin = margin(t = 18),
    ),
    axis.title.y = element_text(
      size = 11,
      color = "#34495e",
      margin = margin(r = 18)
    ),
    axis.text.y = element_text(
      size = 10, 
      color = "#34495e", 
      margin = margin(r = 10)
    ),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10.5),
    plot.margin = margin(30, 20, 20, 20)
  )

ggsave("final\\outputs\\figura_04_ripley_kl.tiff", 
       plot = p_ripley, 
       width = 10, 
       height = 7, 
       dpi = 300, 
       bg = "white",
       compression = "lzw")


# 4.4 Figura 5: Densidad Kernel (recortada a límites municipales)
# Convertir densidad a sf para ggplot2
dens_df <- as.data.frame(dens)
names(dens_df) <- c("x", "y", "densidad")
dx <- unique(sort(dens_df$x))[2] - unique(sort(dens_df$x))[1]
dy <- unique(sort(dens_df$y))[2] - unique(sort(dens_df$y))[1]

crear_celda <- function(x, y, dx, dy) {
  st_polygon(list(matrix(c(x-dx/2, y-dy/2, x+dx/2, y-dy/2,
                           x+dx/2, y+dy/2, x-dx/2, y+dy/2,
                           x-dx/2, y-dy/2), ncol=2, byrow=TRUE)))
}

# Limpiar valores problemáticos antes de crear el sf
dens_df <- as.data.frame(dens)
names(dens_df) <- c("x", "y", "densidad")

# Eliminar valores NA y muy cercanos a 0
dens_df <- dens_df %>%
  filter(!is.na(densidad), densidad > 0)

dx <- unique(sort(dens_df$x))[2] - unique(sort(dens_df$x))[1]
dy <- unique(sort(dens_df$y))[2] - unique(sort(dens_df$y))[1]

crear_celda <- function(x, y, dx, dy) {
  st_polygon(list(matrix(c(x-dx/2, y-dy/2, x+dx/2, y-dy/2,
                           x+dx/2, y+dy/2, x-dx/2, y+dy/2,
                           x-dx/2, y-dy/2), ncol=2, byrow=TRUE)))
}

dens_sf <- st_sf(
  densidad = dens_df$densidad,
  geometry = st_sfc(purrr::map2(dens_df$x, dens_df$y, ~crear_celda(.x, .y, dx, dy)), crs=25831)
)

dens_cortada <- st_intersection(dens_sf, bcn_muni_utm)

# gráfico

p_dens <- ggplot() +
  geom_sf(data = dens_cortada, aes(fill = densidad), color = NA) +
  geom_sf(data = bcn_muni_utm, fill = NA, color = "black", linewidth = 0.8) +
  scale_fill_viridis_c(
    option = "plasma", 
    name = "Densidad relativa",
    trans = "log1p",
    na.value = "transparent"
  ) +
  labs(
    title = "Estimación de intensidad kernel",
    subtitle = paste("h óptimo =", round(h_optimo, 0), "m | Patrón multimodal centro-periferia"),
    x = "Longitud",     
    y = "Latitud",      
    caption = "Elaboración propia. Fuente: Dataset integrado OSM + CKAN"
  ) +
  theme_fivethirtyeight() +
  theme(
    plot.background = element_rect(fill = "white", color = NA), 
    legend.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(
      size = 16,
      face = "bold",
      color = "#2c3e50",
      hjust = 0.5,
      margin = margin(b = 12)
    ),
    plot.subtitle = element_text(
      size = 11,
      color = "#7f8c8d",
      hjust = 0.5,
      margin = margin(b = 18)
    ),
    plot.caption = element_text(
      size = 9,               
      color = "#7f8c8d",     
      hjust = 0.5,            
      margin = margin(t = 15) 
    ),
    axis.text.x = element_text(
      size = 10,
      color = "#34495e",
      margin = margin(t = 12)
    ),
    axis.title.x = element_text(
      size = 11,
      color = "#34495e",
      margin = margin(t = 18),
    ),
    axis.title.y = element_text(
      size = 11,
      color = "#34495e",
      margin = margin(r = 18)
    ),
    axis.text.y = element_text(
      size = 10, 
      color = "#34495e", 
      margin = margin(r = 10)
    ),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10.5),
    plot.margin = margin(30, 20, 20, 20)
  )

# Exportar


ggsave("final\\outputs\\figura_05_densidad_kernel.tiff", 
       plot = p_dens, 
       width = 10, 
       height = 7, 
       dpi = 300, 
       bg = "white",
       compression = "lzw")
