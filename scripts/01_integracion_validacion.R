# ==============================================================================
# 01_integracion_validacion
# ==============================================================================

# 1. CONFIGURACIÓN Y PARÁMETROS ------------------------------------------------
PARAMS <- list(
  crs_utm      = 25831,   # ETRS89/UTM zona 31N 
  crs_wgs      = 4326,    # WGS84 (para compatibilidad web/Leaflet)
  buffer_m     = 50,      # Umbral de coincidencia espacial (Baddeley et al., 2015)
  sig_level    = 0.05,    # Nivel de significancia estadística
  timeout_api  = 1200     # Timeout en segundos para consultas OSM/CKAN
)

library(tidyverse)
library(sf)
library(osmdata)
library(ggspatial)
library(tmap)
library(ggplot2)
library(ggrepel)
library(mapSpain)
library(leaflet)
library(dplyr)
library(httr)
library(jsonlite)
library(ggthemes)
library(leaflet)
library(htmlwidgets)
library(purrr)
library(htmltools)
library(stringr)


dir.create("final\\data", showWarnings = FALSE, recursive = TRUE)
dir.create("final\\data\\processed", showWarnings = FALSE)
dir.create("final\\outputs", showWarnings = FALSE)

# Semilla para reproducibilidad en simulaciones futuras
set.seed(2026)


# 2. LÍMITE MUNICIPAL DE BARCELONA
barcelona_mun <- esp_get_munic() %>%
  filter(cpro == "08", cmun == "019")

bbox_bcn_mun <- st_bbox(barcelona_mun)


# 3. EXTRACCIÓN OSM (API Overpass) 

# Consulta cultura y turismo
q_mun_turismo <- bbox_bcn_mun %>%
  opq(timeout = 1000) %>%
  add_osm_feature(key = "amenity", 
                  value = c("theatre", "arts_centre", "events_venue", 
                            "music_venue", "museum"))

mun_turismo <- osmdata_sf(q_mun_turismo)

mun_turismo_ptos <- mun_turismo$osm_points

q_mun_hoteles <- bbox_bcn_mun %>%
  opq(timeout = 1000) %>%
  add_osm_feature(key = "tourism", 
                  value = c("hotel", "hostel", "museum", "attraction"))

mun_hoteles <- osmdata_sf(q_mun_hoteles)
mun_hoteles_ptos <- mun_hoteles$osm_points

osm_puntos <- bind_rows(
  mun_turismo_ptos %>% mutate(categoria = "turismo"),
  mun_hoteles_ptos %>% mutate(categoria = "hoteles")
) %>%
  # Eliminar duplicados (mismo ID OSM)
  filter(!duplicated(osm_id)) %>%
  # Clasificar tipo 
  mutate(
    tipo = case_when(
      amenity == "theatre" ~ "Teatro",
      amenity == "arts_centre" ~ "Centro de arte", 
      amenity == "events_venue" ~ "Recinto de eventos",
      amenity == "music_venue" ~ "Sala de música",
      amenity == "museum" ~ "Museo",
      tourism == "hotel" ~ "Hotel",
      tourism == "hostel" ~ "Hostel",
      tourism == "attraction" ~ "Atracción",
      TRUE ~ NA_character_
    ),
    # Nombre limpio (prioriza name, luego tipo)
    nombre_clean = str_to_upper(str_trim(coalesce(name, tipo))),
    fuente_original = "OSM"
  ) %>%
  # Solo los que tienen clasificación válida
  filter(!is.na(tipo)) %>%
  select(nombre_clean, tipo, fuente_original, geometry)

# Se reproyectan ambas capas a WGS84 y se intersectan para conservar únicamente 
# los puntos que caen dentro del término municipal.

barcelona_mun <- st_transform(barcelona_mun, 4326)
osm_puntos <- st_transform(osm_puntos, 4326)

osm_bcn <- osm_puntos %>%
  st_intersection(barcelona_mun)


# Unificar puntos y polígonos (centroides) y proyectar a UTM
osm_bcn <- osm_bcn %>%
  mutate(
    categoria = case_when(
      tipo %in% c("Teatro", "Centro de arte", "Recinto de eventos") ~ "Cultura",
      tipo == "Atracción" ~ "Atracción", 
      tipo %in% c("Hotel", "Hostel") ~ "Alojamiento",
      TRUE ~ "Otro"
    )
  ) %>%
  filter(!is.na(tipo))


# 4. EXTRACCIÓN CKAN (Open Data BCN) 
# Función para la página API CKAN
get_ckan_data <- function(resource_id) {
  all_records <- list()
  limit <- 1000
  offset <- 0
  
  repeat {
    url <- paste0(
      "https://opendata-ajuntament.barcelona.cat/data/api/action/datastore_search?",
      "resource_id=", resource_id,
      "&limit=", limit,
      "&offset=", offset
    )
    response <- httr::GET(url)
    data_list <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
    records <- data_list$result$records
    if(length(records) == 0) break
    all_records <- append(all_records, records)
    offset <- offset + limit
  }
  return(dplyr::as_tibble(all_records))
}

# Teatros y auditorios

ckan_teatres <- get_ckan_data("0f706441-b9d8-47c9-9e71-ced453810a72") %>%
  filter(grepl("Teatres|Auditoris", secondary_filters_name, ignore.case = TRUE)) %>%
  filter(!is.na(geo_epgs_4326_lon) & !is.na(geo_epgs_4326_lat)) %>%
  mutate(
    nombre_clean = str_to_upper(str_trim(name)),
    tipo = ifelse(grepl("Teatre", secondary_filters_name, ignore.case = TRUE), 
                  "Teatro", "Auditorio"),
    fuente_original = "CKAN"
  ) %>%
  st_as_sf(coords = c("geo_epgs_4326_lon", "geo_epgs_4326_lat"), crs = 4326) %>%
  select(nombre_clean, tipo, fuente_original, geometry)

# Hoteles

ckan_hoteles <- get_ckan_data("9bccce1b-0b9d-4cc6-94a7-459cb99450de") %>%
  filter(grepl("Hotel", secondary_filters_name, ignore.case = TRUE)) %>%
  filter(!is.na(geo_epgs_4326_lon) & !is.na(geo_epgs_4326_lat)) %>%
  mutate(
    nombre_clean = str_to_upper(str_trim(name)),
    tipo = "Hotel",
    fuente_original = "CKAN"
  ) %>%
  st_as_sf(coords = c("geo_epgs_4326_lon", "geo_epgs_4326_lat"), crs = 4326) %>%
  select(nombre_clean, tipo, fuente_original, geometry)

# Hosteles y otros alojamientos

ckan_altres_allotjaments <- get_ckan_data("a57a2110-3249-450a-8d57-4e1848817d56") %>%
  filter(!is.na(geo_epgs_4326_lon) & !is.na(geo_epgs_4326_lat)) %>%
  mutate(
    nombre_clean = str_to_upper(str_trim(name)),
    categoria = case_when(
      secondary_filters_name == "Albergs juvenils" ~ "Hostel",
      secondary_filters_name == "Apartaments" ~ "Apartament",
      secondary_filters_name == "Residències d'estudiants" ~ "Residència",
      TRUE ~ "Altres allotjaments"
    ),
    fuente_original = "CKAN"
  ) %>%
  st_as_sf(coords = c("geo_epgs_4326_lon", "geo_epgs_4326_lat"), crs = 4326) %>%
  select(nombre_clean, categoria, fuente_original, geometry)

# Unir datos CKAN

ckan_bcn <- bind_rows(ckan_teatres, ckan_hoteles, ckan_altres_allotjaments)

ckan_bcn <- ckan_bcn %>%
  mutate(
    categoria = case_when(
      tipo %in% c("Teatro", "Auditorio", "Recinto de eventos") ~ "Cultura",
      tipo %in% c("Hotel", "Hostel") ~ "Alojamiento",
      TRUE ~ "Otro"
    )
  ) %>%
  filter(!is.na(tipo))

# 5. VALIDACIÓN GEOMÉTRICA

validar_sf <- function(df) {
  df %>%
    mutate(geometry = st_make_valid(geometry)) %>%
    filter(!st_is_empty(geometry), !is.na(nombre_clean), nombre_clean != "") %>%
    distinct(geometry, .keep_all = TRUE)
}

osm_sf  <- validar_sf(osm_bcn)
ckan_sf <- validar_sf(ckan_bcn)


# 6. FUSIÓN ESPACIAL Y DEDUPLICACIÓN 

# Después de validar_sf(), transforma a CRS métrico:
osm_sf  <- st_transform(osm_sf, PARAMS$crs_utm)   # 25831
ckan_sf <- st_transform(ckan_sf, PARAMS$crs_utm)

# Se aplica el buffer
ckan_buffer <- st_buffer(ckan_sf, PARAMS$buffer_m)

# Matriz de intersección (OSM ∩ CKAN_buffer)
match_matrix <- st_intersects(osm_sf, ckan_buffer, sparse = FALSE)
osm_sf$match_ckan <- apply(match_matrix, 1, any)
ckan_sf$match_osm <- apply(t(match_matrix), 1, any)

# Aplicar prioridad: CKAN > OSM en coincidencias
dataset_final <- bind_rows(
  # OSM exclusivos
  osm_sf %>% filter(!match_ckan) %>% mutate(fuente_fusion = "Solo_OSM"),
  # CKAN exclusivos
  ckan_sf %>% filter(!match_osm) %>% mutate(fuente_fusion = "Solo_CKAN"),
  # Coincidentes (se queda CKAN por fiabilidad institucional)
  ckan_sf %>% filter(match_osm) %>% mutate(fuente_fusion = "Ambas")
) %>%
  mutate(
    id_unico = paste0("BCN_", row_number()),
    tipo     = fct_lump(tipo, n = 4, other_level = "Otros") # Simplificar categorías para el informe
  ) %>%
  select(id_unico, nombre_clean, tipo, categoria, fuente_original, fuente_fusion, geometry)






# 7. EXPORTACIÓN Y MÉTRICAS PARA EL INFORME

# Guardar dataset unificado
st_write(dataset_final, "final\\data\\processed\\dataset_final.gpkg", delete_layer = TRUE, quiet = TRUE)
saveRDS(dataset_final, "final\\data\\processed\\dataset_final.rds")

# Tabla resumen 
tabla_categorias <- dataset_final %>%
  st_drop_geometry() %>%
  count(fuente_fusion, categoria, name = "N_registros") %>%
  pivot_wider(names_from = categoria, values_from = N_registros, values_fill = 0) %>%
  mutate(Total = rowSums(across(where(is.numeric)))) %>%
  rename(Fuente = fuente_fusion)  %>%
  select(Fuente, order(colnames(select(., -Fuente)), decreasing = TRUE))


tabla_tipo <- dataset_final %>%
  st_drop_geometry() %>%
  count(fuente_fusion, tipo, name = "N_registros") %>%
  pivot_wider(names_from = tipo, values_from = N_registros, values_fill = 0) %>%
  mutate(Total = rowSums(across(where(is.numeric)))) %>%
  rename(Fuente = fuente_fusion)  %>%
  select(Fuente, order(colnames(select(., -Fuente)), decreasing = TRUE))


# Exportar tabla formateada a HTML

tabla_html_categorias <- knitr::kable(tabla_categorias, 
                           format = "pipe", align = "c") %>%
  kableExtra::kable_styling(latex_options = "striped")

tabla_html_tipo <- knitr::kable(tabla_tipo,
                                      format = "pipe", align = "c") %>%
  kableExtra::kable_styling(latex_options = "striped")

writeLines(tabla_html_categorias, "final\\outputs\\tabla_02_dffinal_cat.html")
writeLines(tabla_html_tipo, "final\\outputs\\tabla_03_dffinal_tipo.html")



# Tabla de cobertura por categoría
tabla_cobertura <- dataset_final %>%
  st_drop_geometry() %>%
  count(categoria, tipo, fuente_fusion, name = "n") %>%
  pivot_wider(names_from = fuente_fusion, values_from = n, values_fill = 0) %>%
  mutate(
    total = rowSums(across(c("Solo_OSM", "Solo_CKAN", "Ambas"), ~replace_na(.x, 0))),
    pct_osm = round(Solo_OSM / total * 100, 1),
    pct_ckan = round(Solo_CKAN / total * 100, 1),
    pct_ambas = round(Ambas / total * 100, 1)
  ) %>%
  select(categoria, tipo, total, Solo_OSM, Solo_CKAN, Ambas, pct_osm, pct_ckan, pct_ambas)

write.csv(tabla_cobertura, "final/outputs/tabla_01_cobertura_fuentes.csv", row.names = FALSE)

# Gráfico de barras: % por categoría
p_cobertura <- tabla_cobertura %>%
  mutate(
    tipo_txt = as.character(tipo),
    categoria_eje = factor(tipo_txt, levels = c("Atracción", "Hostel", "Hotel", "Teatro", "Otros")),
    
    pct_del_total = round(total / sum(total) * 100, 1)
  ) %>%
  ggplot(aes(x = categoria_eje, y = pct_del_total, fill = categoria_eje)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_text(
    aes(label = paste0(pct_del_total, "%")), 
    vjust = -0.5, 
    size = 3.5, 
    color = "#2c3e50",
    fontface = "bold"
  ) +
  scale_fill_viridis_d(option = "plasma", guide = "none") +
  scale_y_continuous(limits = c(0, 65)) + 
  labs(
    title = "Distribución porcentual de registros por categoría",
    subtitle = "Proporción de cada tipo sobre el total del dataset integrado",
    x = "Categoría", 
    y = "% sobre el total de registros",
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

ggsave("final\\outputs\\figura_01_cobertura_fuentes.tiff", 
       plot = p_cobertura, 
       width = 10, 
       height = 7, 
       dpi = 300, 
       bg = "white",
       compression = "lzw")




# Mapa estático
mapa_base <- ggplot() +
  geom_sf(data = barcelona_mun, fill = "#f4f4f4", color = "grey50", linewidth = 0.5) +
  geom_sf(data = dataset_final, aes(color = tipo), size = 1.2, alpha = 0.7) +
  scale_color_viridis_d(option = "plasma", name = "Tipo de equipamiento") +
  labs(
    title = "Distribución espacial de equipamientos culturales y turísticos en Barcelona",
    subtitle = "Dataset integrado OSM + CKAN | n = 867 registros",
    caption = "Elaboración propia. Fuente: OpenStreetMap & Open Data BCN"
  ) +
  coord_sf(expand = FALSE) +
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


ggsave("final\\outputs\\figura_02_mapa_integrado.tiff", 
       plot = mapa_base, 
       width = 10, 
       height = 7, 
       dpi = 300, 
       bg = "white",
       compression = "lzw")


# MAPA INTERACTIVO LEAFLET 

# Preparación de datos y transformación a WGS84
df_web <- dataset_final %>%
  st_transform(4326) %>%
  mutate(
    popup_html = paste0(
      "<b>", ifelse(is.na(nombre_clean), "Sin nombre", nombre_clean), "</b><br/>",
      "Tipo: <b>", tipo, "</b><br/>",
      "Categoría: ", categoria, "<br/>",
      "Fuente: ", fuente_fusion
    )
  )

# Asegurar que bcn_muni esté en el sistema correcto
bcn_muni_leaf <- barcelona_mun %>% st_transform(4326)

colores_personalizados <- c(
  "Atracción" = "#3a3599", 
  "Hostel"    = "#9432b3", 
  "Hotel"     = "#d1678e", 
  "Teatro"    = "#f4a464", 
  "Otros"     = "#edf54e"  
)

pal_tipo <- colorFactor(
  palette  = colores_personalizados, 
  domain   = df_web$tipo,
  na.color = "grey80"
)

lista_tipos <- unique(df_web$tipo)

# Inicialización del mapa base
mapa_interactivo <- leaflet() %>%
  addProviderTiles("OpenStreetMap.Mapnik", group = "OpenStreetMap") %>%
  addProviderTiles(providers$CartoDB.Positron, group = "Claro") %>%
  addProviderTiles(providers$CartoDB.DarkMatter, group = "Oscuro")

# Separar los puntos en el mapa por cada Tipo de establecimiento
for (t in lista_tipos) {
  datos_tipo <- filter(df_web, tipo == t)
  
  mapa_interactivo <- mapa_interactivo %>%
    addCircleMarkers(
      data         = datos_tipo,
      radius       = 6,
      stroke       = TRUE,
      color        = "white",
      weight       = 1.5,
      fillColor    = ~pal_tipo(tipo), 
      fillOpacity  = 0.75,
      popup        = ~popup_html,
      group        = t, 
      clusterOptions = markerClusterOptions(
        showCoverageOnHover = FALSE,
        zoomToBoundsOnClick = TRUE
      )
    )
}

# Creación del mapa
mapa_interactivo <- mapa_interactivo %>%
  # Línea del límite municipal de Barcelona
  addPolylines(
    data         = bcn_muni_leaf,
    color        = "black",
    weight       = 2,
    opacity      = 0.8,
    smoothFactor = 0.5,
    label        = "Límite municipal de Barcelona"
  ) %>%
  
  # Fijar vista centrada en Barcelona
  setView(lng = 2.17, lat = 41.39, zoom = 12) %>%
  
  # Leyenda basada en Tipo 
  addLegend(
    pal      = pal_tipo, 
    values   = df_web$tipo, 
    title    = "Tipo de Equipamiento", 
    position = "bottomright",
    opacity  = 0.8,
    na.label = "Sin datos"
  ) %>%
  
  # caja de título
  addControl(
    HTML(paste0(
      "<div style='padding: 15px; background: rgba(255,255,255,0.95); 
                  border-radius: 8px; border: 2px solid #dee2e6; 
                  box-shadow: 0 4px 8px rgba(0,0,0,0.1); 
                  font-family: Arial, sans-serif; font-size: 13px;
                  min-width: 260px; margin-bottom: 5px;'>
         <div style='font-size: 16px; font-weight: bold; color: #2c3e50; 
                     margin-bottom: 6px; text-align: center;'>
           Oferta Turística y Cultural BCN
         </div>
         <div style='font-size: 12px; color: #7f8c8d; 
                     margin-bottom: 10px; text-align: center; line-height: 1.4;'>
           Distribución espacial de establecimientos hoteleros, teatros y atracciones en la ciudad.
         </div>
         <div style='margin-top: 8px; font-size: 11px; color: #7f8c8d; 
                     border-top: 1px solid #eee; padding-top: 6px; text-align: center;'>
           Fuente: Fusion Data BCN | ", Sys.Date(), "
         </div>
       </div>"
    )),
    position = "topright"
  ) %>%
  
  # control de capas
  addLayersControl(
    baseGroups   = c("OpenStreetMap", "Claro", "Oscuro"),
    overlayGroups = lista_tipos,
    options      = layersControlOptions(collapsed = TRUE)
  )

# Exportar
saveWidget(
  widget       = mapa_interactivo, 
  file         = "final/outputs/mapa_interactivo.html", 
  selfcontained = TRUE, 
  title        = "Mapa Interactivo - Equipamientos BCN"
)
