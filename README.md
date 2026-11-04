# RattleGen: verificación computacional de la ley de Faraday

Proyecto final de Electromagnetismo (Universidad de La Salle, Bogotá). Tema 10: generador AC/DC casero y verificación de la ley de Faraday.

**Autores:** Juan Camilo Almanza y Juan David Pinzón. **Profesor:** Fernando Vergara.

## Qué hace

`rattlegen_sim.py` simula el generador RattleGen de LabdeSyn publicado en rimstar.org
(https://rimstar.org/science_electronics_projects/make_led_flashlight_rattlegen_simple_squeeze_one_ba.htm):
un imán que se mueve junto a la bobina de un relé induce una fem; un puente rectificador y un capacitor de 1000 µF la convierten en corriente continua para encender un LED.

1. **Faraday:** calcula ε = −N dΦ/dt con derivada numérica y la compara con la expresión analítica.
2. **Escalamientos:** ε_pico ∝ N y ε_pico ∝ 1/τ, con barras de error por Monte Carlo (N ±2 %, τ ±3 %, Φ ±5 %).
3. **Línea base AC:** bobina que gira, V_rms = V_pico/√2.
4. **Circuito:** bobina (R, L) → puente de 4 diodos → capacitor → LED; 20 corridas Monte Carlo con parámetros inciertos.

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

## Limitaciones

Es una verificación **computacional**: los resultados confirman que el modelo cumple la ley de Faraday, no son una medición física independiente. Los parámetros marcados como supuestos no vienen de la fuente y definen los tiempos de carga del capacitor.

## Referencia

LabdeSyn, «Simple squeeze LED flashlight – RattleGen», rimstar.org.
