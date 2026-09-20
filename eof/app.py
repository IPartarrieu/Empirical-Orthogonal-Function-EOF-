import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

from data_loader import VARIABLES, load_variable
from eof_model import compute_eof, lat_weights, monthly_anomalies

st.set_page_config(page_title="EOF Bío-Bío/Araucanía", page_icon="🌦️", layout="centered")

st.title("🌦️ Modos de variabilidad climática (EOF)")
st.markdown(
    "Descomposición en **Funciones Ortogonales Empíricas** de temperatura y "
    "precipitación (ERA5, 1940-2024) sobre la zona de Bío-Bío/Araucanía. "
    "Cada modo es un patrón espacial fijo cuya intensidad varía en el tiempo "
    "según su serie de coeficientes (PC)."
)


@st.cache_data
def get_eof(var_key, n_modes=6):
    values, lat, lon, time, months, unit, cmap = load_variable(var_key)
    nt, nlat, nlon = values.shape
    data2d = values.reshape(nt, nlat * nlon)
    anom = monthly_anomalies(data2d, months)
    weights = np.repeat(lat_weights(lat), nlon)
    eofs, pcs, var_frac = compute_eof(anom, weights, n_modes=n_modes)
    eofs_grid = eofs.T.reshape(n_modes, nlat, nlon)
    return eofs_grid, pcs, var_frac, lat, lon, time, unit, cmap


with st.sidebar:
    st.header("Controles")
    var_key = st.selectbox("Variable", list(VARIABLES.keys()))
    eofs_grid, pcs, var_frac, lat, lon, time, unit, cmap = get_eof(var_key)
    modo = st.slider("Modo EOF", 1, len(var_frac), 1)

k = modo - 1
fig, (ax_map, ax_pc) = plt.subplots(2, 1, figsize=(7, 7), height_ratios=[1.3, 1])

vmax = np.abs(eofs_grid[k]).max()
im = ax_map.pcolormesh(
    lon, lat, eofs_grid[k], cmap=cmap, vmin=-vmax, vmax=vmax, shading="auto"
)
ax_map.set_title(f"EOF modo {modo} — {var_frac[k]*100:.1f}% de varianza explicada")
ax_map.set_xlabel("Longitud")
ax_map.set_ylabel("Latitud")
fig.colorbar(im, ax=ax_map, label=f"Carga del patrón [{unit}]", shrink=0.85)

ax_pc.plot(pd.to_datetime(time), pcs[:, k], color="#2c3e50", lw=0.9)
ax_pc.axhline(0, color="0.7", lw=1)
ax_pc.set_title(f"Componente principal (PC) del modo {modo}")
ax_pc.set_xlabel("Año")
ax_pc.set_ylabel("Amplitud")
ax_pc.grid(alpha=0.3)

fig.tight_layout()
st.pyplot(fig)

st.markdown("**Varianza explicada por modo:**")
st.bar_chart(pd.Series(var_frac * 100, index=[f"Modo {i+1}" for i in range(len(var_frac))]))

with st.expander("¿Qué es una EOF?"):
    st.markdown(
        """
Una EOF descompone un campo que varía en espacio y tiempo en una suma de
patrones espaciales fijos, cada uno con una serie de tiempo propia (PC)
que indica cuánto de ese patrón está presente en cada instante. Se
calcula vía SVD sobre las anomalías mensuales (con el ciclo estacional
ya removido), ponderando cada celda de la grilla por su área real
(∝ cos(latitud)).

Implementación fiel al script `EOF.m` original de la tarea del curso
*Análisis Estadísticos en Climatología* — ver
[`eof_model.py`](https://github.com/IPartarrieu/Empirical-Orthogonal-Function-EOF-/blob/main/eof/eof_model.py)
y sus verificaciones matemáticas en
[`tests/test_eof_model.py`](https://github.com/IPartarrieu/Empirical-Orthogonal-Function-EOF-/blob/main/eof/tests/test_eof_model.py).
"""
    )
