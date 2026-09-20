"""Funciones Ortogonales Empiricas (EOF), via SVD.

Sigue la misma formulacion clasica de EOF.m (Kaplan, 2001) usada en el
curso "Analisis Estadisticos en Climatologia": para una matriz de datos
U (tiempo x espacio), el SVD U = C @ diag(s) @ EOFs^T da:

  - EOFs (patrones espaciales) = columnas de V (filas de EOFs^T)
  - EC / PCs (coeficientes de expansion) = C @ diag(s)
  - Varianza explicada por el modo k = s[k]^2 / sum(s^2)

Aqui se agrega ponderacion por latitud (sqrt(cos(lat))), practica
estandar para que las celdas de la grilla pesen segun su area real
antes de la descomposicion.
"""

import numpy as np


def monthly_anomalies(data, months):
    """Resta el ciclo climatologico mensual en cada punto de grilla.

    Parameters
    ----------
    data : ndarray (n_time, n_space)
    months : ndarray (n_time,) con el mes (1-12) de cada paso de tiempo.
    """
    anom = np.empty_like(data)
    for m in range(1, 13):
        idx = months == m
        anom[idx] = data[idx] - data[idx].mean(axis=0)
    return anom


def lat_weights(lat):
    """Peso por celda segun su area real (proporcional a cos(lat))."""
    return np.sqrt(np.cos(np.deg2rad(lat)))


def compute_eof(anomalies, weights=None, n_modes=10):
    """Calcula EOFs, PCs y varianza explicada via SVD.

    Parameters
    ----------
    anomalies : ndarray (n_time, n_space)
        Anomalias (con el ciclo estacional ya removido).
    weights : ndarray (n_space,), opcional
        Peso por punto de grilla (p.ej. de `lat_weights`, repetido a lo
        largo de longitud). Se aplica antes del SVD y se remueve de los
        patrones espaciales al final, para que EOFs quede en las
        unidades originales de la variable.
    n_modes : int
        Numero de modos a retornar.

    Returns
    -------
    eofs : ndarray (n_space, n_modes)
    pcs : ndarray (n_time, n_modes)
    variance_fraction : ndarray (n_modes,)
    """
    if weights is None:
        weights = np.ones(anomalies.shape[1])

    weighted = anomalies * weights[np.newaxis, :]

    U, s, Vt = np.linalg.svd(weighted, full_matrices=False)

    pcs = U * s
    eofs = Vt.T
    # Remueve la ponderacion de los patrones espaciales para dejarlos en
    # las unidades originales de la variable.
    eofs = eofs / weights[:, np.newaxis]

    variance_fraction = (s**2) / np.sum(s**2)

    n_modes = min(n_modes, eofs.shape[1])
    return eofs[:, :n_modes], pcs[:, :n_modes], variance_fraction[:n_modes]
