# Empirical-Orthogonal-Function-EOF-

Demo interactiva de descomposición en **Funciones Ortogonales Empíricas (EOF)** sobre datos reales de temperatura y precipitación (ERA5, 1940-2024) en la zona de Bío-Bío/Araucanía.

🔗 **Demo en vivo:** _pendiente de desplegar — ver instrucciones abajo_

## Qué hace

Elegí la variable (temperatura o precipitación) y el modo EOF, y la app muestra al instante el patrón espacial de ese modo y su serie de tiempo (componente principal), junto con el % de varianza explicada.

El método (SVD sobre las anomalías mensuales, ponderadas por área) es una implementación en Python fiel al script [`EOF.m`](matlab_exercise/EOF.m) original de la tarea del curso *Análisis Estadísticos en Climatología* — ver [`eof/eof_model.py`](eof/eof_model.py) y sus verificaciones matemáticas en [`eof/tests/test_eof_model.py`](eof/tests/test_eof_model.py) (ortogonalidad de los EOFs, no correlación de los PCs, reconstrucción exacta).

## Correr en local

```bash
cd eof
pip install -r requirements.txt
streamlit run app.py
```

## Desplegar gratis (Streamlit Community Cloud)

1. Andá a [share.streamlit.io](https://share.streamlit.io) y conectate con tu cuenta de GitHub.
2. "New app" → elegí este repo, branch `main`, main file path `eof/app.py`.
3. Deploy. En un par de minutos tenés el link público.

## Ejercicio original (MATLAB)

El script `EOF.m` y las tareas por estación del curso original quedaron en [`matlab_exercise/`](matlab_exercise/).
