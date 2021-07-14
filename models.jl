using TTauUtils

star = Star("hart94", 2, 0.5, 4000, 15)
savestar(star)

Ṁs = [1e-7, 3e-8, 1e-8, 3e-9, 1e-9, 3e-10]
Ṁ_strings = ["70", "75", "80", "85", "90", "95"]

T_maxs = [7000, 7500, 8000, 8500, 9000, 9500, 10500, 11500]
T_max_strings = ["7000", "7500", "8000", "8500", "9000", "9500", "10500", "11500"]
mags = [(2.0, 3.0), (4.0, 6.0), (8.0, 10.0), (8.0, 12.0)]
mag_strings = ["2-3", "4-6", "8-10", "12-10"]

for (m, Ṁ) in enumerate(Ṁs)
    Ṁ_str = Ṁ_strings[m]
    for (T_max, T_str) in zip(T_maxs[m:m+2], T_max_strings[m:m+2])
        for ((r_mi, r_mo), mag_str) in zip(mags, mag_strings)
            name = star.name*'_'*Ṁ_str*'_'*T_str*'_'*mag_str
            println(name)
        end
    end
end