# 📊 Análisis Espacial de Oferta Cultural y Turística en Barcelona

**Autor:** María Yu García Muñoz  
**Asignatura:** Datos Espaciales y Espaciotemporales (DEE) – Universitat de València  
**Fecha:** Mayo 2026

## 📖 Descripción

Este proyecto integra datos espaciales abiertos para analizar la distribución de equipamientos culturales y turísticos en la ciudad de Barcelona. El objetivo es responder preguntas de investigación sobre patrones de aglomeración, identificar "desiertos culturales" y proponer aplicaciones estratégicas para la planificación territorial y el sector privado.

El flujo de trabajo combina **ingeniería de datos** (APIs y fusión espacial), **estadística espacial** (procesos puntuales y autocorrelación regional) y **visualización** de resultados.

## 🏗️ Estructura del Repositorio

```text
📦 Proyecto_Final_DEE
 ┣  01_integracion_validacion.R   # Extracción, limpieza y fusión OSM + CKAN
 ┣  02_patrones_puntuales.R       # Análisis CSR, Ripley y Kernel Density
 ┣  03_metricas_autocorrelacion.R # Moran, LISA y métricas de negocio por barrio
 ┣ 📂 final/                        # Directorio de salida (resultados)
 ┃  📂 data/processed/             # Datasets intermedios (.rds)
 ┃ ┗ 📂 outputs/                    # Figuras (.tiff) y tablas (.csv/html)
 ┗ 📄 README.md
```

## 🔧 Metodología

El análisis se ha estructurado en tres etapas secuenciales:

### 1. Integración de Datos y Validación (`01`)
- **Fuentes:** OpenStreetMap (Overpass API) y Open Data Barcelona (CKAN).
- **Proceso:** Limpieza de coordenadas, reproyección a UTM (EPSG:25831), validación geométrica y fusión espacial (buffer de 50m).
- **Output:** `dataset_final.rds` (dataset unificado y validado).

### 2. Análisis de Patrones Puntuales (`02`)
- **Técnicas:** Test de Aleatoriedad Espacial Completa (CSR), función G de vecino más cercano, funciones K y L de Ripley y estimación de intensidad por Kernel Density.
- **Objetivo:** Determinar si la distribución de equipamientos es aleatoria o responde a clustering estructural.

### 3. Autocorrelación Espacial y Negocio (`03`)
- **Escala:** Análisis agregado a nivel de barrio (73 unidades).
- **Técnicas:** Matriz de vecindad Queen, índice de Moran global y Local Indicators of Spatial Association (LISA).
- **Métricas:** Ratio de equipamientos por 10k habitantes y cobertura territorial (buffer 1km).

## 📂 Fuentes de Datos

- **OpenStreetMap:** Consultadas vía API Overpass (Amenities y Tourism).
- **Open Data BCN:** Portal oficial del Ayuntamiento de Barcelona (CKAN API).
- **Límites Administrativos:** Obtenidos vía paquete `mapSpain`.

## 📜 Licencia

Este proyecto se distribuye bajo licencia MIT. Consulte el archivo `LICENSE` para más detalles.

---
> *Este proyecto forma parte del Grado en Inteligencia y Analítica de Negocios (BIA) de la Universitat de València.*
