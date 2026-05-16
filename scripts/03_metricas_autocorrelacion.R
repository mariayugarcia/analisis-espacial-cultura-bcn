# ==============================================================================
# 03_metricas_autocorrelacion.R
# ==============================================================================


library(sf)
library(tidyverse)
library(spdep)
library(ggplot2)
library(mapSpain)
library(kableExtra)
library(ggthemes)
library(knitr)
library(jsonlite)
library(dplyr)

set.seed(2026)

# 1. CONFIGURACIÓN Y CARGA DE DATOS 
dataset_final <- readRDS("final/data/processed/dataset_final.rds")

# Obtener barrios
# Definir el endpoint de la API SQL y la consulta
resource_id <- "b21fa550-56ea-4f4c-9adc-b8009381896e"
sql_query   <- paste0("SELECT * FROM \"", resource_id, "\"")

# URL codificada para la API
url_api <- paste0(
  "https://opendata-ajuntament.barcelona.cat/data/api/action/datastore_search_sql?sql=",
  URLencode(sql_query)
)

# Realizar la petición y leer el JSON
respuesta <- fromJSON(url_api)

# 3. Extraer los datos
datos_tabla <- respuesta$result$records
# CONVERTIR A OBJETO ESPACIAL (sf)
barrios_api <- st_as_sf(datos_tabla, wkt = "geometria_wgs84", crs = 4326) %>%
  st_transform(25831) %>%
  mutate(
    district_name = nom_districte, 
    district_id   = as.numeric(codi_districte),
    neighborhood_name = nom_barri,
    neighborhood_id   = as.numeric(codi_barri)
  )

summary(barrios_api)


# 2. PREPARACIÓN DE CAPAS ESPACIALES 

resource_id_padron <- "eb82adf2-a7b0-40e6-9624-b4b9eff23018"

# Traer todos los registros (usamos un límite alto para asegurar todo el padrón)
url_padron <- paste0(
  "https://opendata-ajuntament.barcelona.cat/data/api/action/datastore_search?resource_id=", 
  resource_id_padron, 
  "&limit=100000"
)

respuesta_padron <- fromJSON(url_padron)
tabla_padron <- respuesta_padron$result$records

# Normalizar nombres de columnas a minúsculas para evitar conflictos de sintaxis
names(tabla_padron) <- tolower(names(tabla_padron))

# Procesar y agrupar la población por barrio
poblacion_barrios <- tabla_padron %>%
  mutate(
    neighborhood_name = nom_barri,
    cantidad = as.numeric(valor)
  ) %>%
  group_by(neighborhood_name) %>%
  summarise(poblacion = sum(cantidad, na.rm = TRUE), .groups = "drop")

barrios_api <- barrios_api %>%
  mutate(neighborhood_name_clean = str_trim(tolower(neighborhood_name)))

poblacion_barrios <- poblacion_barrios %>%
  mutate(neighborhood_name_clean = str_trim(tolower(neighborhood_name)))

barrios_bcn_completo <- barrios_api %>%
  left_join(
    poblacion_barrios %>%
      select(neighborhood_name_clean, poblacion),
    by = "neighborhood_name_clean"
  ) %>%
  select(-neighborhood_name_clean)

barrios_bcn_completo <- barrios_bcn_completo %>%
  rename(geometry = geometria_wgs84) %>%
  st_as_sf()

# 3. AGREGACIÓN DE PUNTOS POR DISTRITO

# Filtrar solo cultura y atracciones
equipamientos_cultura <- dataset_final %>%
  filter(categoria %in% c("Cultura", "Atracción")) %>%
  st_transform(25831)

summary(barrios_bcn_completo)

# Spatial join: asignar cada punto a su barrio
puntos_por_barrio <- equipamientos_cultura %>%
  st_join(
    barrios_bcn_completo %>%
      select(
        neighborhood_id,
        neighborhood_name,
        geometry
      ),
    join = st_intersects,
    left = FALSE
  ) %>%
  st_drop_geometry() %>%
  count(
    neighborhood_id,
    neighborhood_name,
    name = "n_equipamientos"
  )


# Unir con población y calcular indicadores
metrics_df <- barrios_bcn_completo %>%
  mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) %>%
  st_drop_geometry() %>%
  
  select(
    neighborhood_id,
    neighborhood_name,
    poblacion,
    area_km2
  ) %>%
  
  left_join(
    puntos_por_barrio,
    by = c("neighborhood_id", "neighborhood_name")
  ) %>%
  mutate(n_equipamientos = replace_na(n_equipamientos, 0)) %>%
  
  mutate(
    ratio_10k    = round((n_equipamientos / poblacion) * 10000, 2),
    ratio_km2    = round(n_equipamientos / area_km2, 2),
    densidad_pob = round(poblacion / area_km2, 0)
  )



# 4. COBERTURA ESPACIAL (BUFFER 1000m)

# Crear buffers de 1000m alrededor de cada equipamiento
buffers_1km <- st_buffer(equipamientos_cultura, 1000)

# Calcular área cubierta por barrio
cobertura_por_barrio <- st_intersection(buffers_1km, barrios_bcn_completo) %>%
  group_by(nom_barri) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  mutate(area_cubierta_km2 = as.numeric(st_area(.)) / 1e6) %>%
  select(nom_barri, area_cubierta_km2)

# Unir cobertura al dataframe principal
metrics_df <- metrics_df %>%
  left_join(
    cobertura_por_barrio, 
    by = c("neighborhood_name" = "nom_barri") 
  ) %>%
  mutate(
    area_cubierta_km2 = replace_na(area_cubierta_km2, 0),
    cobertura_pct = round((area_cubierta_km2 / area_km2) * 100, 1)
  )


# 5. AUTOCORRELACIÓN ESPACIAL (MORAN, LISA) 
# Preparar objeto sf completo para análisis espacial

metrics_sf <- barrios_bcn_completo %>%
  left_join(
    metrics_df %>% select(neighborhood_name, area_km2, n_equipamientos, ratio_10k, 
                          ratio_km2, densidad_pob, area_cubierta_km2, cobertura_pct), 
    by = "neighborhood_name"
  )
names()

metrics_sf <- metrics_sf %>%
  select(
    district_id,
    district_name,
    neighborhood_id,
    neighborhood_name,
    
    poblacion,
    area_km2,
    densidad_pob,
    n_equipamientos,
    ratio_10k,
    ratio_km2,
    area_cubierta_km2,
    cobertura_pct,
    
    geometry
  )

# Matriz de vecindad Queen (contigüidad por vértice o arista)
vecindad <- poly2nb(metrics_sf, queen = TRUE)

# Pesos espaciales estandarizados por fila
pesos <- nb2listw(vecindad, style = "W", zero.policy = TRUE)

# Moran's I global
moran_res <- moran.test(metrics_sf$ratio_10k, pesos, zero.policy = TRUE)

# LISA (Local Moran's I)
lisa_res <- localmoran(metrics_sf$ratio_10k, pesos, zero.policy = TRUE)

# Calcular spatial lag
lag_espacial <- lag.listw(pesos, metrics_sf$ratio_10k, zero.policy = TRUE)

# Clasificación de clusters LISA
metrics_sf <- metrics_sf %>%
  mutate(
    Ii              = lisa_res[, 1],
    p_value         = lisa_res[, 5],
    lag             = lag.listw(pesos, ratio_10k, zero.policy = TRUE),
    media_ratio     = mean(ratio_10k, na.rm = TRUE),
    media_lag       = mean(lag, na.rm = TRUE),
    
    # Clasificación de clústeres espaciales (LISA)
    lisa_cluster = case_when(
      p_value > 0.05 ~ "No Significativo",
      ratio_10k > media_ratio & lag > media_lag ~ "Alto-Alto (Hotspot)",
      ratio_10k < media_ratio & lag < media_lag ~ "Bajo-Bajo (Coldspot)",
      ratio_10k > media_ratio & lag < media_lag ~ "Alto-Bajo (Atípico)",
      ratio_10k < media_ratio & lag > media_lag ~ "Bajo-Alto (Atípico)",
      TRUE ~ "No Significativo"
    )
  )

# 6. CLASIFICACIÓN DE OPORTUNIDADES DE NEGOCIO 

metrics_sf <- metrics_sf %>%
  mutate(
    # Clasificación basada en ratio y cobertura
    cluster_negocio = case_when(
      ratio_10k >= quantile(ratio_10k, 0.75) ~ "Saturado",
      ratio_10k <= quantile(ratio_10k, 0.25) & cobertura_pct < 70 ~ "Prioritario",
      ratio_10k <= quantile(ratio_10k, 0.25) & cobertura_pct >= 70 ~ "Equilibrado",
      TRUE ~ "Equilibrado"
    ),
    
    # Prioridad de inversión (1 = más prioritario)
    prioridad_inversion = case_when(
      cluster_negocio == "Prioritario" ~ 1,
      cluster_negocio == "Equilibrado" ~ 2,
      cluster_negocio == "Saturado" ~ 3
    )
  )

# 7. EXPORTACIÓN 

# Tabla 1: Indicadores por distrito
tabla_indicadores <- metrics_sf %>%
  st_drop_geometry() %>%
  select(neighborhood_name, poblacion, area_km2, n_equipamientos, 
         ratio_10k, ratio_km2, cobertura_pct, cluster_negocio, prioridad_inversion) %>%
  arrange(prioridad_inversion)

tabla_indicadores <- tabla_indicadores %>%
  rename(
    "Barrio"                     = neighborhood_name,
    "Población (hab)"              = poblacion,
    "Superficie (km^2)"             = area_km2,
    "Total Equipamientos"          = n_equipamientos,
    "Ratio por 10.000 hab"         = ratio_10k,
    "Equipamientos por km^2"        = ratio_km2,
    "Cobertura Espacial (%)"       = cobertura_pct,
    "Clasificación LISA"           = cluster_negocio,
    "Prioridad de Inversión"       = prioridad_inversion
  )

write.csv(tabla_indicadores, "final/outputs/tabla_05_indicadores_distrito.csv", row.names = FALSE)

# Tabla HTML formateada

tabla_html <- knitr::kable(tabla_indicadores, 
                           caption = "Tabla 5. Indicadores territoriales y clasificación de oportunidades por distrito",
                           format = "pipe", align = "c", digits = 2) %>%
  kableExtra::kable_styling(latex_options = "striped")


writeLines(tabla_html, "final/outputs/tabla_05_indicadores_distrito.html")

# Tabla 2: Resultados de autocorrelación espacial
tabla_moran <- tibble(
  "Indicador" = c("Moran's I", "Esperanza", "Varianza", "z-score", "p-value"),
  "Valor" = c(
    round(as.numeric(moran_res$estimate[1]), 4),    # Índice de Moran observado
    round(as.numeric(moran_res$estimate[2]), 4),    # Esperanza matemática teórica
    round(as.numeric(moran_res$estimate[3]), 4),    # Varianza teórica del test
    round(as.numeric(moran_res$statistic), 4),      # El estadístico estandarizado Z
    format.pval(moran_res$p.value, digits = 3)      # El P-valor del análisis
  )
)

write.csv(tabla_moran, "final/outputs/tabla_06_moran_global.csv", row.names = FALSE)

# 8. MAPAS 

# Mapa 1: Cluster LISA
mapa_lisa <- ggplot() +
  geom_sf(data = metrics_sf, aes(fill = lisa_cluster), color = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = c(
      "No Significativo" = "#3a3599",
      "Alto-Alto (Hotspot)" = "#edf54e",
      "Bajo-Bajo (Coldspot)" = "#9432b3",
      "Alto-Bajo (Atípico)" = "#f4a464",
      "Bajo-Alto (Atípico)" = "#d1678e"
    ),
    name = "Cluster LISA"
  ) +
  labs(
    title = "Autocorrelación espacial: Clusters LISA",
    subtitle = paste("Moran's I =", round(moran_res$estimate[1], 3), 
                     "| p =", format.pval(moran_res$p.value, digits = 2)),
    x = "Longitud",
    y = "Latitud",
    caption = "Elaboración propia. Fuente: Dataset integrado + Idescat"
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

ggsave("final/outputs/figura_06_lisa_clusters.tiff", mapa_lisa,
       width = 10, height = 7, dpi = 300, bg = "white", compression = "lzw")


# Mapa 2: Oportunidades de negocio

# Crear la capa de etiquetas
etiquetas_distritos <- metrics_sf %>%
  group_by(district_name) %>%
  summarise(.groups = "drop") %>% 
  st_point_on_surface()

# Generar el mapa
mapa_negocio <- ggplot() +
  geom_sf(data = metrics_sf, aes(fill = cluster_negocio), color = "white", linewidth = 0.4) +
  scale_fill_manual(
    values = c(
      "Saturado" = "#edf54e",
      "Equilibrado" = "#f4a464",
      "Prioritario" = "#d1678e"
    ),
    name = "Zona de Negocio"
  ) +
  # Capa de texto 
  geom_sf_text(
    data = etiquetas_distritos, 
    aes(label = district_name), 
    size = 2.8,              # Reducido ligeramente para que quepan todos
    color = "black", 
    fontface = "bold",
    check_overlap = FALSE   
  ) +
  labs(
    title = "Oportunidades de negocio cultural por distrito y barrio",
    subtitle = "Clasificación basada en ratio equipamientos/10k hab. y cobertura 1 km",
    x = "Longitud",
    y = "Latitud",
    caption = "Elaboración propia. Fuente: Dataset integrado + Idescat 2024"
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


ggsave("final/outputs/figura_07_oportunidades_negocio.tiff", mapa_negocio,
       width = 10, height = 7, dpi = 300, bg = "white", compression = "lzw")

# Scatterplot: Moran bivariante (ratio vs lag espacial)
moran_scatter <- ggplot(metrics_sf %>% st_drop_geometry(), 
                        aes(x = lag, y = ratio_10k)) +
  geom_point(aes(color = lisa_cluster), size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.8) +
  geom_hline(yintercept = mean(metrics_sf$ratio_10k), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = mean(metrics_sf$lag), linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = c(
      "No Significativo" = "#3a3599",
      "Alto-Alto (Hotspot)" = "#edf54e",
      "Bajo-Bajo (Coldspot)" = "#9432b3",
      "Alto-Bajo (Atípico)" = "#f4a464",
      "Bajo-Alto (Atípico)" = "#d1678e"
    ),
    name = "Cluster"
  ) +
  labs(
    title = "Moran Scatterplot: Ratio cultural vs. Vecinos",
    subtitle = "Eje X: Media de distritos vecinos | Eje Y: Ratio del distrito",
    x = "Lag espacial (media de vecinos)",
    y = "Ratio equipamientos/10k hab."
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

ggsave("final/outputs/figura_08_moran_scatter.tiff", moran_scatter,
       width = 9, height = 6, dpi = 300, bg = "white", compression = "lzw")

# Tabla 8: Barrios significativos (p-value < 0.05)
barrios_significativos <- metrics_sf %>%
  st_drop_geometry() %>%
  filter(p_value < 0.05) %>%
  select(
    neighborhood_name,
    district_name,
    ratio_10k,
    lisa_cluster,
    p_value
  ) %>%
  mutate(
    p_value = round(p_value, 4),
    ratio_10k = round(ratio_10k, 2)
  ) %>%
  arrange(p_value) %>%
  rename(
    "Barrio" = neighborhood_name,
    "Distrito" = district_name,
    "Ratio/10k hab" = ratio_10k,
    "Cluster" = lisa_cluster,
    "p-value" = p_value
  )

tabla_barrios_html <- knitr::kable(
  barrios_significativos,
  caption = "Tabla 8. Barrios con autocorrelación espacial significativa (p < 0.05)",
  format = "pipe",
  align = c("c", "c", "c", "c", "c"),
  digits = 4
) %>%
  kableExtra::kable_styling(latex_options = "striped")

writeLines(tabla_barrios_html, "final/outputs/tabla_08_barrios_significativos.html")
write.csv(barrios_significativos, "final/outputs/tabla_08_barrios_significativos.csv", row.names = FALSE)
