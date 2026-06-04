# Aula 13 — Obtenção de covariáveis ambientais para ambientipagem

Data: 03 de junho de 2026  
Carga horária: 4h
Professor convidado: Gustavo Marcatti (UFSJ)

## Conteúdo abordado

- Obtenção de covariáveis ambientais para ambientipagem
- Delimitação de uma TPE a partir de pontos genotípicos
- Criação de polígono, raster e grade espacial de pontos
- Download de covariáveis edáficas via SoilGrids
- Download de covariáveis climáticas via NASA POWER
- Reamostragem espacial de camadas ambientais
- Extração de covariáveis ambientais para amostras e grade espacial
- Exportação dos arquivos finais para uso em análises posteriores

## Material

- `Biometria_da_Interação_GxA_GMP0164_aula13_Gustavo.pdf` — roteiro da aula

## Ferramentas

- Python
- WinPython 3.12.10.1
- QGIS 3.44.11-1

## Dependências principais

```bash
pip install pandas geopandas rioxarray xarray geocube soilgrids rasterio matplotlib requests
