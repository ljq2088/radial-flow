import csv

with open('comparison_results.csv', 'r') as f:
    reader = csv.DictReader(f)
    data = [r for r in reader if float(r['omega']) >= 0.3]

# Sort by omega, then a
data.sort(key=lambda x: (float(x['omega']), float(x['a']), int(x['l']), int(x['m'])))

print(f"{'a':<6} {'ω':<6} {'l':<3} {'m':<3} {'λ':<8} {'Ref ratio':<12} {'Calc ratio':<12} {'Rel Error':<10}")
print("="*80)

for r in data:
    a = float(r['a'])
    omega = float(r['omega'])
    l = int(r['l'])
    m = int(r['m'])
    lam = float(r['lambda'])
    calc = float(r['calc_ratio'])

    if r['ref_ratio']:
        ref = float(r['ref_ratio'])
        rel_err = abs(calc - ref) / ref * 100 if ref > 1e-10 else 0
        print(f"{a:<6.1f} {omega:<6.1f} {l:<3} {m:<3} {lam:<8.2f} {ref:<12.4e} {calc:<12.4e} {rel_err:<10.2f}%")
    else:
        print(f"{a:<6.1f} {omega:<6.1f} {l:<3} {m:<3} {lam:<8.2f} {'N/A':<12} {calc:<12.4e} {'N/A':<10}")

# Save to file
with open('high_omega_comparison.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['a', 'omega', 'l', 'm', 'lambda', 'ref_ratio', 'calc_ratio', 'rel_error_%'])
    for r in data:
        calc = float(r['calc_ratio'])
        if r['ref_ratio']:
            ref = float(r['ref_ratio'])
            rel_err = abs(calc - ref) / ref * 100
            writer.writerow([r['a'], r['omega'], r['l'], r['m'], r['lambda'], ref, calc, rel_err])
        else:
            writer.writerow([r['a'], r['omega'], r['l'], r['m'], r['lambda'], 'N/A', calc, 'N/A'])

print(f"\nSaved to high_omega_comparison.csv")
