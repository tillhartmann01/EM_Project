"""EE4C05 track assignment - 50 Hz steel / 20 kHz ferrite transformer, E65/32/27.
All formulas follow the lecture slides (Magnetic circuit, Diffusive fields, Design of magnetic components).
Items marked ASSUMPTION are design choices or values that must be checked against the datasheets.
"""
import numpy as np
from scipy.special import jv

mu0 = 4e-7 * np.pi
rho_cu20, a_cu = 1.72e-8, 0.00393
def rho_cu(T): return rho_cu20 * (1 + a_cu * (T - 20))

# ---------------- geometry (assignment drawing + TDK datasheet) ----------------
w_c, d_c = 20.0e-3, 27.4e-3          # centre leg width, stack depth
A_geo = w_c * d_c                    # 548 mm2 (steel, 100 % fill)
Ae_fer, le, Ve = 535e-6, 147e-3, 78650e-9   # TDK ferrite values
b_c = 44e-3                          # winding height (coil former)
A_w = 426e-6                         # winding area, TDK B66388B (1 core set)
l_T = 2*(20.5e-3 + 27.6e-3) + np.pi * 9.7e-3   # mean turn, own geometry (bobbin hole + half build)
J = 5e6                              # A/m2
V1, V2 = 230.0, 15.0
print(f"Geometric centre-leg area {A_geo*1e6:.0f} mm2, mean turn {l_T*1e3:.0f} mm, window {A_w*1e6:.0f} mm2")

# =============================== 50 Hz =====================================
print("\n===== 50 Hz, steel =====")
f = 50.0; Bpk = 1.6; kst = 0.95   # ASSUMPTION stacking factor
for name, Ac in [("100 % fill (A=548 mm2)", A_geo), ("kst=0.95", kst*A_geo)]:
    N1 = V1/(4.44*f*Bpk*Ac)
    print(f"{name}: N1 = {N1:.0f}, N2 = {N1*V2/V1:.1f}")
Ac = kst*A_geo
N1 = 1242; N2 = 81
print(f"chosen N1={N1}, N2={N2}: B = {V1/(4.44*f*N1*Ac):.3f} T ; with N1=1180: B(kst) = {V1/(4.44*f*1180*Ac):.3f} T")

kw = 0.40                             # ASSUMPTION copper fill, solid round wire
Acu = kw*A_w/2
print(f"Cu area per winding {Acu*1e6:.1f} mm2")
S_ideal = V1 * J*Acu/N1
print(f"ideal rating S = {S_ideal:.0f} VA")
d1, d2 = 0.28e-3, 1.12e-3             # standard enamelled wire sizes
A1, A2 = np.pi*d1**2/4, np.pi*d2**2/4
I1 = J*A1; I2 = I1*N1/N2
print(f"d1={d1*1e3} mm (A={A1*1e6:.4f} mm2) -> I1={I1:.3f} A ; I2={I2:.2f} A in d2={d2*1e3} mm (J2={I2/A2/1e6:.2f} A/mm2)")
S50 = V1*I1; print(f"rating S = {S50:.1f} VA   (S2 = {V2*I2:.1f} VA)")
print(f"wire turns per layer (b_c=44 mm, d_out={d1*1e3+0.03:.2f} mm): {b_c/(d1+0.03e-3):.0f}  -> layers {N1/(b_c/(d1+0.03e-3)):.1f}")
T = 75; R1 = rho_cu(T)*N1*l_T/A1; R2 = rho_cu(T)*N2*l_T/A2
print(f"R1={R1:.1f} ohm R2={R2*1e3:.0f} mohm ; Pcu = {I1**2*R1:.2f} + {I2**2*R2:.2f} W")
# skin depth in copper
for fr in (50, 20e3):
    for TT in (20, 100):
        dl = np.sqrt(2*rho_cu(TT)/(2*np.pi*fr*mu0))
        print(f"  delta_Cu({fr:g} Hz, {TT} C) = {dl*1e3:.2f} mm")

# core material: 23RK075 assumptions
print("\n-- 23RK075 (datasheet missing: ASSUMED 0.23 mm, P1.7/50 ~ 0.75 W/kg) --")
h = 0.23e-3; rho_fe = 46e-8; sigma_fe = 1/rho_fe; dens = 7650.0
for mur in (2000, 20000):
    dl = np.sqrt(2/(2*np.pi*f*mu0*mur*sigma_fe)); x = h/dl
    Ppl = (np.pi**2/2)*(Bpk*f)**2*h*dl*sigma_fe*(np.sinh(x)-np.sin(x))/(np.cosh(x)-np.cos(x))
    print(f"  mur={mur}: delta={dl*1e3:.2f} mm, h/delta={x:.2f}, classical eddy loss = {Ppl:.0f} W/m3 = {Ppl/dens:.3f} W/kg")
Ppl_cl = (np.pi**2/6)*sigma_fe*h**2*f**2*Bpk**2
print(f"  thin-sheet limit pi^2 sigma h^2 f^2 B^2/6 = {Ppl_cl/dens:.3f} W/kg")
m_core = Ve*dens*kst*(A_geo/Ae_fer)**0   # approx, same envelope
print(f"  core mass approx {m_core:.2f} kg; Pfe at 0.55-0.75 W/kg (1.6 T) = {0.55*m_core:.2f}-{0.75*m_core:.2f} W")

# magnetising and leakage
H16 = 56.0   # A/m ASSUMPTION: read from 23RK075 B-H curve at 1.6 T (implies mu_r ~ 22 600)
mur_sec = 1.6/(mu0*H16)
for N1_, N2_, lab in [(1180, 77, "N1=1180/N2=77"), (1242, 81, "N1=1242/N2=81")]:
    Ac_ = A_geo if N1_ == 1180 else Ac
    L1 = mu0*mur_sec*N1_**2*Ac_/le
    L2 = L1*(N2_/N1_)**2; M = L1*N2_/N1_
    print(f"  {lab}: mur_sec={mur_sec:.0f}  L1={L1:.1f} H  L2={L2:.3f} H  M={M:.2f} H  (k=1)")
    Im = H16*le/N1_
    print(f"     Im_pk={Im*1e3:.0f} mA ({Im/np.sqrt(2)*1e3:.0f} mA rms) = {Im/np.sqrt(2)/I1*100:.0f} % of I1")
print(f"  L2 for N2=154 (2x77) with L1=146.5 H: {146.5*(154/1180)**2:.2f} H ; M = {np.sqrt(146.5*146.5*(154/1180)**2):.1f} H")

# leakage (slide 17/28 style, Ampere's law in the window)
def Lsig_concentric(N, hp, hs, g): return mu0*N**2*l_T/b_c*(hp/3 + hs/3 + g)
def Lsig_sandwich(N, hp, hs, g):    return mu0*N**2*l_T/b_c*0.5*(hs/3 + g + hp/6)
hp, hs, g = 2.8e-3, 2.0e-3, 0.5e-3   # ASSUMPTION build-up (9 layers of 0.31 mm, 81 turns of 1.12 mm in 2 layers)
print(f"  50 Hz leakage (primary-referred): concentric {Lsig_concentric(N1,hp,hs,g)*1e3:.1f} mH ; sandwich {Lsig_sandwich(N1,hp,hs,g)*1e3:.1f} mH")
print(f"  ratio Lsig/L1 = {Lsig_concentric(N1,hp,hs,g)/(mu0*mur_sec*N1**2*Ac/le):.1e}")

# ======================= Tata Hi-Lite non-oriented alternative =======================
print("\n===== Tata Hi-Lite NO grades: split loss with modified Steinmetz Pv = kh f B^2 + ke f^2 B^2 (alpha=1, beta=2 assumed) =====")
grades = [
 ("NO10-1270N",0.10,12.7,23.6,1.52,1.63),("NO12-1180N",0.12,11.8,24.2,1.54,1.64),
 ("NO15-1230N",0.15,12.3,25.6,1.54,1.64),("NO18-1280N",0.18,12.8,27.5,1.55,1.65),
 ("NO20-1350N",0.20,13.5,29.1,1.55,1.65),("NO18-1160H",0.18,11.6,25.5,1.51,1.61),
 ("NO20-1200H",0.20,12.0,27.5,1.52,1.62),("NO20-1300H",0.20,13.0,28.0,1.53,1.63),
 ("NO25-1350H",0.25,13.5,31.5,1.52,1.62),("NO27-1400H",0.27,14.0,32.4,1.55,1.64),
 ("NO27-1450H",0.27,14.5,34.4,1.52,1.62),("NO30-1500H",0.30,15.0,36.0,1.52,1.62),
 ("NO35-1750H",0.35,17.5,42.2,1.52,1.62)]
print(f"{'grade':12s} {'kh':>8s} {'ke':>9s} {'P50@1.0T':>9s} {'P50@1.5T':>9s} {'H(1.6T)<=':>10s}")
rows = []
for n, th, p4, p7, j25, j50 in grades:
    A = np.array([[400, 400**2],[700, 700**2]]); kh, ke = np.linalg.solve(A, [p4, p7])
    p50 = kh*50 + ke*50**2
    H16_ub = 2500 + (1.6-j25)/(j50-j25)*2500 if j50 >= 1.6 else np.nan
    rows.append((n, kh, ke, p50, p50*1.5**2, H16_ub))
    print(f"{n:12s} {kh:8.4f} {ke:9.2e} {p50:9.2f} {p50*2.25:9.2f} {H16_ub:10.0f}")
best = min(rows, key=lambda r: r[3]); print("lowest 50 Hz loss:", best[0])
nh = min([r for r in rows if not np.isnan(r[5])], key=lambda r: r[5]); print("lowest H at 1.6 T (upper bound):", nh[0], f"{nh[5]:.0f} A/m")
Hno = 3900.0
Im = Hno*le/N1
print(f"NO steel at 1.6 T, H~{Hno:.0f} A/m (linear-interp upper bound): Im_pk = {Im:.2f} A ({Im/np.sqrt(2):.2f} A rms) vs rated I1 = {I1:.2f} A")

# =============================== 20 kHz =====================================
print("\n===== 20 kHz, ferrite (TDK E65/32/27) =====")
f = 20e3
for N1, N2 in [(31, 2), (46, 3), (15, 1), (16, 1)]:
    B = V1/(4.44*f*N1*Ae_fer); Bsq = V1/(4*f*N1*Ae_fer)
    print(f"N1={N1}, N2={N2}: ratio {N1/N2:.2f} (target 15.33), V2={V1*N2/N1:.2f} V, B sine = {B*1e3:.0f} mT, B if 230 V square = {Bsq*1e3:.0f} mT")
N1, N2 = 31, 2
B = V1/(4.44*f*N1*Ae_fer)
print(f"chosen N1={N1}, N2={N2}, B = {B*1e3:.0f} mT")
kw = 0.30                             # ASSUMPTION copper fill with litz
Acu = kw*A_w/2
I1 = J*Acu/N1; I2 = J*Acu/N2
print(f"Cu area per winding {Acu*1e6:.1f} mm2; I1={I1:.2f} A, I2={I2:.0f} A, S = {V1*I1:.0f} VA (S2 = {V1*N2/N1*I2:.0f} VA)")
S20 = V1*I1
print(f"area-product check: S = 4.44 f B Ae J kw Aw / 2 = {4.44*f*B*Ae_fer*J*kw*A_w/2:.0f} VA")
# materials
print("Core loss references (datasheet, per set, volume 78.65 cm3):")
for m, P, fr, Bm in [("N27",14.6,25e3,0.2),("N87",6.7,100e3,0.1),("N97",6.3,100e3,0.1),("N95",6.5,100e3,0.1),("N92",7.5,100e3,0.1),("N88",8.2,100e3,0.1)]:
    print(f"  {m}: {P} W/set @ {Bm*1e3:.0f} mT {fr/1e3:.0f} kHz,100C = {P/Ve*1e-3:.0f} kW/m3")
for a, b_ in [(1.3, 2.6), (1.5, 2.5)]:
    Pn27 = 14.6*(20e3/25e3)**a*(B/0.2)**b_
    Pn87 = 6.7*(20e3/100e3)**a*(B/0.1)**b_
    print(f"  Steinmetz-scaled to 20 kHz/{B*1e3:.0f} mT (alpha={a}, beta={b_}): N27 ~{Pn27:.1f} W, N87 ~{Pn87:.1f} W per set (ROUGH, check TDK curves)")
# magnetising
for m, AL in [("N27",7200e-9),("N87",7900e-9)]:
    L1 = AL*N1**2; Im = V1/(2*np.pi*f*L1)
    print(f"  {m}: Lm = AL N1^2 = {L1*1e3:.2f} mH, Im = {Im:.3f} A rms ({Im/I1*100:.1f} % of I1), Hpk={np.sqrt(2)*Im*N1/le:.0f} A/m")
# skin depth and strand
dl100 = np.sqrt(2*rho_cu(100)/(2*np.pi*f*mu0)); dl20 = np.sqrt(2*rho_cu(20)/(2*np.pi*f*mu0))
print(f"delta_Cu(20 kHz) = {dl20*1e3:.2f} mm (20 C), {dl100*1e3:.2f} mm (100 C)")
def FR_round(d, fr, T=75):
    r0 = d/2; k = np.sqrt(2*np.pi*fr*mu0/rho_cu(T)); z = 1j**1.5*k*r0
    Z = 1j**1.5*k*jv(0, z)/(2*np.pi*r0/rho_cu(T)*jv(1, z))   # slide 16, per m
    Rdc = rho_cu(T)/(np.pi*r0**2)
    return Z.real/Rdc
ds = 0.1e-3
print(f"skin factor Fr(strand 0.1 mm, 20 kHz) = {FR_round(ds, f):.4f}; solid 1.62 mm wire: {FR_round(1.62e-3, f):.2f}")
# litz sizing
n1 = (I1/J)/(np.pi*ds**2/4); print(f"primary litz: Acu/turn = {I1/J*1e6:.2f} mm2 -> {n1:.0f} strands of 0.1 mm")
fill_litz = 0.55
d_out = np.sqrt(4*(I1/J)/(np.pi*fill_litz)); print(f"  bundle OD ~ {d_out*1e3:.2f} mm (assumed bundle copper fill {fill_litz})")
n2 = (I2/J)/(np.pi*ds**2/4); print(f"secondary: Acu/turn = {I2/J*1e6:.1f} mm2 -> {n2:.0f} strands -> use several parallel litz / flat litz")
# DC copper loss
Vcu = kw*A_w*l_T
Pdc = rho_cu(75)*J**2*Vcu
print(f"Cu volume {Vcu*1e6:.1f} cm3 ; DC copper loss at 75 C = {Pdc:.1f} W")
# proximity in strands (lecture 'eddy loss in applied flux', round wire d<<delta)
sig = 1/rho_cu(75); w = 2*np.pi*f
Hpk_max = np.sqrt(2)*(N1*I1/2)/b_c   # sandwich: peak MMF = N I / 2
Hpk_max_nosand = np.sqrt(2)*(N1*I1)/b_c
for lab, Hm in [("sandwich", Hpk_max), ("plain concentric", Hpk_max_nosand)]:
    Bm = mu0*Hm
    Ps = np.pi*sig*w**2*Bm**2*ds**4/128          # W/m per strand at the peak field
    Lstr = N1*n1*l_T
    P = Ps*Lstr/3                                 # mean of H^2 over a linear profile = 1/3 of peak^2
    print(f"  {lab}: Hpk,max={Hm:.0f} A/m (B={Bm*1e3:.1f} mT) -> proximity loss in primary ~ {P:.3f} W")
# solid wire proximity for comparison: Dowell-like via Fr already, report DC and factor
print(f"  solid 1.62 mm primary would have AC/DC ~ {FR_round(1.62e-3,f):.1f}x on skin alone -> {Pdc/2*FR_round(1.62e-3,f):.0f} W for the primary")
# leakage at 20 kHz
hp, hs, g = 4.4e-3, 4.3e-3, 0.5e-3
print(f"20 kHz leakage (primary-referred): concentric {Lsig_concentric(N1,hp,hs,g)*1e6:.1f} uH ; sandwich S/2-P-S/2 {Lsig_sandwich(N1,hp,hs,g)*1e6:.2f} uH")
Ls = Lsig_sandwich(N1,hp,hs,g); L1m = 7900e-9*N1**2
print(f"  k = 1 - Lsig/(Lm+Lsig) = {1-Ls/(L1m+Ls):.5f};  voltage drop w*Lsig*I1 = {w*Ls*I1:.1f} V ({w*Ls*I1/V1*100:.1f} %)")
# summary
print(f"\nPOWER RATINGS: 50 Hz ~ {S50:.0f} VA ; 20 kHz ~ {S20:.0f} VA ; ratio {S20/S50:.0f}x (f B ratio = {20e3*B/(50*1.6):.0f}, kw ratio {0.3/0.4:.2f})")
