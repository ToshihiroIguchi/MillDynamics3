PBF = r"C:\Users\toshi\python\MillDynamics3\crates\mill-core\src\pbf.rs"
PRB = r"C:\Users\toshi\python\MillDynamics3\crates\mill-core\examples\drag_probe.rs"

def load(p):
    return open(p, encoding="utf-8").read()

def save(p, s):
    open(p, "w", encoding="utf-8", newline="\n").write(s)

s = load(PBF)
a = "const BALL_FORM_DRAG_COEFFICIENT: f32 = 1.0;\n"
assert s.count(a) == 1
s = s.replace(a, a + "\n/// TEMPORARY (exchange-rule experiment): reference speed u0 (m/s, f32 bits) of the step 3.5 exchange; 0 = shipped rule.\npub static EXCHANGE_U0: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);\n")
a = "                    let impulse_mag = full_impulse_mag.min(reduced_mass * arrest_speed);\n"
assert s.count(a) == 1
s = s.replace(a, """                    let u0 = f32::from_bits(EXCHANGE_U0.load(std::sync::atomic::Ordering::Relaxed));
                    let shape = if u0 > 0.0 { (arrest_speed / u0).min(1.0) } else { 1.0 };
                    let impulse_mag = full_impulse_mag.min(reduced_mass * arrest_speed * shape);
""")
save(PBF, s)
s = load(PRB)
s = s.replace("use mill_core::pbf::FluidParticles;", "use mill_core::pbf::{FluidParticles, EXCHANGE_U0};\nuse std::sync::atomic::Ordering;", 1)
a = '    let (wet, st) = (opt("--wet"), opt("--st"));\n'
assert a in s
s = s.replace(a, a + '    EXCHANGE_U0.store(one("--u0", 0.0).to_bits(), Ordering::Relaxed);\n', 1)
save(PRB, s)
print("ok")
