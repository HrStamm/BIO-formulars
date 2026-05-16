def mutations_per_division(genome_bp, error_rate_per_bp):
    # genome [bp], error rate [1/bp] -> mutations per division
    return genome_bp * error_rate_per_bp


def fraction_mutated_intra_species(mut_per_gen, T_years, gen_time_years, genome_bp):
    # accumulated mutations / genome
    n_gen = T_years / gen_time_years
    return mut_per_gen * n_gen / genome_bp


def fraction_diverged_two_species(mut_per_gen, T_years, gen_time_years, genome_bp):
    # two lineages -> factor 2
    return 2 * mut_per_gen * (T_years / gen_time_years) / genome_bp


if __name__ == "__main__":
    print(f"mut/div (3e9 bp, 1e-9/bp) = {mutations_per_division(3e9, 1e-9):.2f}")
    f = fraction_mutated_intra_species(mut_per_gen=3, T_years=1e6,
                                       gen_time_years=20, genome_bp=3e9)
    print(f"intra-species fraction over 1 Myr = {f:.3e}")
    d = fraction_diverged_two_species(mut_per_gen=3, T_years=6e6,
                                      gen_time_years=20, genome_bp=3e9)
    print(f"human-chimp divergence (6 Myr) = {d:.3e}  ({d*100:.2f} %)")
