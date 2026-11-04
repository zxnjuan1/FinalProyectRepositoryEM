
## Cómo ejecutarlo

```
pip install -r requirements.txt
python rattlegen_sim.py
```

Genera en `resultados/`: `resumen_grafico.png/.pdf`, `resultados.json` y `resultados.csv`. La semilla aleatoria es fija, así que los resultados se reproducen.

## Parámetros

| Parámetro | Valor | Origen |
|---|---|---|
| Pico sin carga | 35 V | Página de Rimstar |
| Capacitor | 1000 µF / 16 V | Página de Rimstar |
| Puente rectificador, LED | sí | Página de Rimstar |
| N, τ, R_c, L, V_f, V_LED, R_LED, apretones/s | ver el código | **Supuestos** de los autores |

El flujo máximo Φ₀ se calibra para reproducir el pico de 35 V de la página.
