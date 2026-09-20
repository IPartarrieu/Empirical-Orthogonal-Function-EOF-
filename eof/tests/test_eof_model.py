import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from eof_model import compute_eof, lat_weights, monthly_anomalies  # noqa: E402

rng = np.random.default_rng(0)
n_time, n_space = 240, 30
data = rng.normal(size=(n_time, n_space)) * 2 + 10

months = np.tile(np.arange(1, 13), n_time // 12)
anom = monthly_anomalies(data, months)

# 1) La anomalia mensual promedio por mes debe ser ~0 en cada punto
for m in range(1, 13):
    idx = months == m
    assert np.allclose(anom[idx].mean(axis=0), 0, atol=1e-10)

weights = np.ones(n_space)
eofs, pcs, var_frac = compute_eof(anom, weights, n_modes=n_space)

# 2) Los EOFs son ortonormales (columnas unitarias y perpendiculares)
gram = eofs.T @ eofs
assert np.allclose(gram, np.eye(n_space), atol=1e-8)

# 3) Los PCs son mutuamente no correlacionados (matriz pcs^T @ pcs diagonal)
pc_gram = pcs.T @ pcs
off_diag = pc_gram - np.diag(np.diag(pc_gram))
assert np.allclose(off_diag, 0, atol=1e-6)

# 4) La fraccion de varianza es decreciente y suma 1 al usar todos los modos
assert np.all(np.diff(var_frac) <= 1e-12)
assert abs(var_frac.sum() - 1.0) < 1e-10

# 5) Reconstruccion exacta usando todos los modos: EOFs @ PCs^T == anomalias
recon = pcs @ eofs.T
assert np.allclose(recon, anom, atol=1e-8)

# 6) lat_weights: peso maximo en el ecuador, cero en el polo
assert np.isclose(lat_weights(np.array([0.0]))[0], 1.0)
assert np.isclose(lat_weights(np.array([90.0]))[0], 0.0, atol=1e-6)

print("Todas las verificaciones de EOF pasaron OK")
print("Varianza explicada, primeros 3 modos:", var_frac[:3])
