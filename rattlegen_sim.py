"""
Simulación del RattleGen (generador de inducción de Faraday con puente rectificador).
Proyecto final de Electromagnetismo, Universidad de La Salle - Tema 10.
Referencia base: https://rimstar.org/science_electronics_projects/make_led_flashlight_rattlegen_simple_squeeze_one_ba.htm

Uso:  python rattlegen_sim.py        (genera la carpeta resultados/)
"""
import json, csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

es = lambda x: x.replace(".", ",")      # coma decimal en los textos de las figuras
rng = np.random.default_rng(7)          # semilla fija: resultados reproducibles
OUT = "resultados"; os.makedirs(OUT, exist_ok=True)

# ---------------- PARÁMETROS ----------------
# [PÁGINA] = dato de la página de Rimstar/LabdeSyn · [SUPUESTO] = valor elegido por los autores
P = dict(
    V_pico_vacio=35.0,   # V, [PÁGINA] pico sin carga del relé pequeño con imán de neodimio
    C=1000e-6,           # F, [PÁGINA] capacitor de 1000 µF / 16 V
    N=2000,              # espiras de la bobina del relé, [SUPUESTO]
    tau=4e-3,            # s, duración del movimiento de la palanca, [SUPUESTO]
    Rc=300.0,            # Ω, resistencia de la bobina, [SUPUESTO]
    L=0.3,               # H, inductancia de la bobina, [SUPUESTO]
    Vf=0.7,              # V, caída por diodo de silicio, [SUPUESTO]
    V_led=3.0,           # V, umbral del LED blanco, [SUPUESTO]
    R_led=1000.0,        # Ω, resistencia en serie con el LED, [SUPUESTO]
    f_press=4.0,         # apretones por segundo, [SUPUESTO]
    hold=0.05,           # s, tiempo con la palanca apretada, [SUPUESTO]
)
# Flujo máximo calibrado para reproducir el pico de 35 V de la página: eps_pico = N·Φ0·π/(2τ)
P["Phi0"] = 2 * P["tau"] * P["V_pico_vacio"] / (np.pi * P["N"])

# ---------------- MODELO DE FLUJO Y FEM (Faraday) ----------------
def flujo_y_fem(t, Phi0, N, tau, T, hold):
    """Φ(t) de un apretón periódico (sube en tau, se mantiene, baja en tau) y ε = -N dΦ/dt analítica."""
    p = t % T
    up, dn = p < tau, (p >= hold) & (p < hold + tau)
    Phi = np.where(up, Phi0 * (1 - np.cos(np.pi * p / tau)) / 2,
          np.where(p < hold, Phi0,
          np.where(dn, Phi0 * (1 + np.cos(np.pi * (p - hold) / tau)) / 2, 0.0)))
    dPhi = np.where(up, Phi0 * np.pi / (2 * tau) * np.sin(np.pi * p / tau),
           np.where(dn, -Phi0 * np.pi / (2 * tau) * np.sin(np.pi * (p - hold) / tau), 0.0))
    return Phi, -N * dPhi

# ---------------- CIRCUITO: bobina (R,L) -> puente -> C -> LED ----------------
def circuito(Phi0, Rc, C, f_press, Vf, tmax=20.0, dt=5e-5, N=None, tau=None, L=None):
    N = N or P["N"]; tau = tau or P["tau"]; L = L or P["L"]
    T = 1 / f_press
    t = np.arange(0, tmax, dt)
    _, eps = flujo_y_fem(t, Phi0, N, tau, T, P["hold"])
    Vc = np.zeros_like(t); i = 0.0; v = 0.0; icoil = np.zeros_like(t)
    for k in range(1, len(t)):
        e = eps[k]; vb = v + 2 * Vf
        if i > 0 or e > vb:
            i += (e - Rc * i - vb) / L * dt; i = max(i, 0.0)
        elif i < 0 or e < -vb:
            i += (e - Rc * i + vb) / L * dt; i = min(i, 0.0)
        carga = max(0.0, v - P["V_led"]) / P["R_led"]
        v = max(0.0, v + (abs(i) - carga) / C * dt)
        Vc[k] = v; icoil[k] = i
    return t, eps, Vc, icoil

# ---------------- 1) Verificación de Faraday: ε numérica vs analítica ----------------
T0 = 1 / P["f_press"]; tf = np.arange(0, T0, 1e-5)
Phi, eps_an = flujo_y_fem(tf, P["Phi0"], P["N"], P["tau"], T0, P["hold"])
eps_num = -P["N"] * np.gradient(Phi, tf)
err_rel = np.max(np.abs(eps_num - eps_an)) / np.max(np.abs(eps_an))

# ---------------- 2) Escalamientos con barras de error (Monte Carlo) ----------------
def pico(N, tau, Phi0):
    t = np.arange(0, 0.02, 1e-5)
    _, e = flujo_y_fem(t, Phi0, N, tau, 1.0, 0.1)
    return np.max(np.abs(e))
def mc_pico(N, tau, n=300):
    v = [pico(N * rng.normal(1, .02), tau * rng.normal(1, .03), P["Phi0"] * rng.normal(1, .05)) for _ in range(n)]
    return np.mean(v), np.std(v, ddof=1)
Ns = np.array([500, 1000, 1500, 2000, 2500, 3000])
pN = np.array([mc_pico(n, P["tau"]) for n in Ns])
invt = np.array([100, 150, 200, 250, 300, 400.0]); taus = 1 / invt
pT = np.array([mc_pico(P["N"], tt) for tt in taus])
def ajuste_origen(x, y, s):
    w = 1 / s**2; m = np.sum(w * x * y) / np.sum(w * x * x); sm = 1 / np.sqrt(np.sum(w * x * x))
    chi = np.sum(w * (y - m * x)**2) / (len(x) - 1); return m, sm * max(1, np.sqrt(chi)), chi
mN, smN, cN = ajuste_origen(Ns, pN[:, 0], pN[:, 1]); mT, smT, cT = ajuste_origen(invt, pT[:, 0], pT[:, 1])
teo_N = P["Phi0"] * np.pi / (2 * P["tau"]); teo_T = P["N"] * P["Phi0"] * np.pi / 2

# ---------------- 3) Línea base AC senoidal: V_rms = V_pico/√2 ----------------
def rms_seno(f, n=200):
    out = []
    for _ in range(n):
        Nb, B, A = 200 * rng.normal(1, .02), 0.05 * rng.normal(1, .03), 4e-4 * rng.normal(1, .02)
        t = np.linspace(0, 10 / f, 20001)           # 10 ciclos completos
        e = -Nb * np.gradient(B * A * np.cos(2 * np.pi * f * t), t)
        out.append((np.sqrt(np.mean(e**2)), np.max(np.abs(e))))
    o = np.array(out); return o[:, 0].mean(), o[:, 0].std(ddof=1), np.mean(o[:, 0] / o[:, 1])
fs = np.array([2, 4, 6, 8, 10, 12.0]); sin_res = np.array([rms_seno(f) for f in fs])
mS, smS, cS = ajuste_origen(fs, sin_res[:, 0], sin_res[:, 1]); razon = sin_res[:, 2].mean()

# ---------------- 4) RattleGen: carga del capacitor con barras de banda (MC) ----------------
t, eps, Vc0, i0 = circuito(P["Phi0"], P["Rc"], P["C"], P["f_press"], P["Vf"])
runs = []
for _ in range(20):
    _, _, Vc, _ = circuito(P["Phi0"] * rng.normal(1, .05), P["Rc"] * rng.normal(1, .05), P["C"] * rng.normal(1, .10),
                           P["f_press"] * rng.normal(1, .10), P["Vf"] + rng.normal(0, .05))
    runs.append(Vc)
runs = np.array(runs); media, desv = runs.mean(0), runs.std(0, ddof=1)
def t_on(V): idx = np.argmax(V >= P["V_led"] + 0.05); return t[idx] if V.max() >= P["V_led"] + 0.05 else np.nan
ton = np.array([t_on(r) for r in runs]); ton_m, ton_s = np.nanmean(ton), np.nanstd(ton, ddof=1)
apretones = ton_m * P["f_press"]; Vmax = Vc0.max(); Eenerg = 0.5 * P["C"] * Vmax**2

# ---------------- FIGURA RESUMEN ----------------
plt.rcParams.update({"font.family": "serif", "font.serif": ["Liberation Serif", "Times New Roman", "DejaVu Serif"],
                     "font.size": 9, "axes.linewidth": .7, "mathtext.fontset": "stix"})
G, R_, B_, K = "#14532D", "#C0392B", "#1F4E9C", "#1A1A1A"
fig, ax = plt.subplots(2, 3, figsize=(12, 6.8)); a = ax.ravel()
w0 = tf < 0.1
a[0].plot(tf[w0] * 1e3, Phi[w0] * 1e6, color=B_, label=r"$\Phi(t)$ (µWb)"); a[0].set_xlabel("t (ms)"); a[0].set_ylabel("Φ (µWb)", color=B_)
b0 = a[0].twinx(); b0.plot(tf[w0] * 1e3, eps_an[w0], color=G, lw=2.2, label="ε analítica")
b0.plot(tf[w0][::40] * 1e3, eps_num[w0][::40], "o", color=R_, ms=3, label="ε = −N dΦ/dt (numérica)"); b0.set_ylabel("ε (V)", color=G)
b0.legend(fontsize=7, loc="upper right"); a[0].set_title(es(f"(a) Ley de Faraday: error relativo máx. = {err_rel*100:.2f} %"), fontsize=9)
a[1].errorbar(Ns, pN[:, 0], pN[:, 1], fmt="o", color=G, capsize=3, label="simulación ± σ"); a[1].plot(Ns, mN * Ns, "-", color=K, lw=1, label=es(f"ajuste: {mN:.4f}·N"))
a[1].set_xlabel("N (espiras)"); a[1].set_ylabel("ε$_{pico}$ (V)"); a[1].legend(fontsize=7); a[1].set_title("(b) ε$_{pico}$ $\propto$ N", fontsize=9)
a[2].errorbar(invt, pT[:, 0], pT[:, 1], fmt="s", color=B_, capsize=3, label="simulación ± σ"); a[2].plot(invt, mT * invt, "-", color=K, lw=1, label=es(f"ajuste: {mT:.3f}/τ"))
a[2].set_xlabel("1/τ (s$^{-1}$)"); a[2].set_ylabel("ε$_{pico}$ (V)"); a[2].legend(fontsize=7); a[2].set_title("(c) ε$_{pico}$ $\propto$ 1/τ (rapidez del apretón)", fontsize=9)
a[3].errorbar(fs, sin_res[:, 0], sin_res[:, 1], fmt="o", color=R_, capsize=3, label="V$_{rms}$ simulado ± σ"); a[3].plot(fs, mS * fs, "-", color=K, lw=1, label=es(f"ajuste: {mS:.4f}·f"))
a[3].set_xlabel("f (Hz)"); a[3].set_ylabel("V$_{rms}$ (V)"); a[3].legend(fontsize=7); a[3].set_title(es(f"(d) AC senoidal: V$_{{rms}}$/V$_{{pico}}$ = {razon:.4f} (1/√2 = 0.7071)"), fontsize=9)
z = (t > 1.0) & (t < 1.0 + T0)
a[4].plot(t[z], np.abs(eps[z]), color=B_, lw=1, label="|ε| (bobina)"); a[4].plot(t[z], Vc0[z], color=G, lw=2, label="V$_C$ (tras el puente)")
a[4].set_xlabel("t (s)"); a[4].set_ylabel("V"); a[4].legend(fontsize=7); a[4].set_title("(e) AC → DC: puente rectificador + capacitor", fontsize=9)
a[5].fill_between(t, media - desv, media + desv, color=G, alpha=.25, label="±1σ (Monte Carlo, 20 corridas)"); a[5].plot(t, media, color=G, lw=1.6, label="V$_C$ medio")
a[5].axhline(P["V_led"], color=R_, ls="--", lw=1, label="umbral del LED"); a[5].set_xlabel("t (s)"); a[5].set_ylabel("V$_C$ (V)")
a[5].legend(fontsize=7, loc="lower right"); a[5].set_title(es(f"(f) El LED enciende a t = ({ton_m:.1f} ± {ton_s:.1f}) s ≈ {apretones:.0f} apretones"), fontsize=9)
for x in a: x.grid(alpha=.25, lw=.4)
fig.suptitle("RattleGen: verificación computacional de la ley de Faraday y rectificación AC → DC", fontsize=11, color=G, weight="bold")
fig.tight_layout(rect=(0, 0, 1, .96)); fig.savefig(f"{OUT}/resumen_grafico.png", dpi=300); fig.savefig(f"{OUT}/resumen_grafico.pdf")

# ---------------- SALIDAS ----------------
res = dict(parametros={k: float(v) for k, v in P.items()},
           faraday_error_relativo=float(err_rel),
           pendiente_pico_vs_N=dict(sim=float(mN), sigma=float(smN), teoria=float(teo_N), chi2_nu=float(cN)),
           pendiente_pico_vs_invtau=dict(sim=float(mT), sigma=float(smT), teoria=float(teo_T), chi2_nu=float(cT)),
           senoidal=dict(pendiente_rms_vs_f=float(mS), sigma=float(smS), razon_rms_pico=float(razon), teoria_razon=float(1/np.sqrt(2))),
           rattlegen=dict(V_max_nominal=float(Vmax), t_LED_s=float(ton_m), sigma_t_LED=float(ton_s), apretones_LED=float(apretones), energia_J=float(Eenerg)))
json.dump(res, open(f"{OUT}/resultados.json", "w"), indent=2, ensure_ascii=False)
with open(f"{OUT}/resultados.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["N", "eps_pico_V", "sigma_V"]); w.writerows(zip(Ns, pN[:, 0], pN[:, 1]))
    w.writerow([]); w.writerow(["1/tau_1/s", "eps_pico_V", "sigma_V"]); w.writerows(zip(invt, pT[:, 0], pT[:, 1]))
    w.writerow([]); w.writerow(["f_Hz", "Vrms_V", "sigma_V"]); w.writerows(zip(fs, sin_res[:, 0], sin_res[:, 1]))
print(json.dumps(res, indent=2, ensure_ascii=False))
