import pandas as pd
from pathlib import Path

raw_path = Path('results/cassiopeia_512cell_30site_10state_missing0_03_50repeats.csv')
summary_path = Path('results/cassiopeia_512cell_30site_10state_missing0_03_50repeats_summary.csv')
wide_path = Path('results/cassiopeia_512cell_30site_10state_missing0_03_50repeats_wide.csv')

df = pd.read_csv(raw_path)
summary = df.groupby(['missing_rate', 'algorithm']).agg(
    n=('repeat', 'count'),
    triplet_mean=('triplet_accuracy', 'mean'),
    triplet_median=('triplet_accuracy', 'median'),
    triplet_q25=('triplet_accuracy', lambda x: x.quantile(0.25)),
    triplet_q75=('triplet_accuracy', lambda x: x.quantile(0.75)),
    rf_mean=('rf_similarity', 'mean'),
    rf_median=('rf_similarity', 'median'),
    rf_q25=('rf_similarity', lambda x: x.quantile(0.25)),
    rf_q75=('rf_similarity', lambda x: x.quantile(0.75)),
    coverage_mean=('sampled_triplet_coverage', 'mean'),
    solve_seconds_mean=('solve_seconds', 'mean'),
).reset_index()
summary.to_csv(summary_path, index=False)

wide = summary.pivot_table(
    index='algorithm',
    columns='missing_rate',
    values=['triplet_mean', 'triplet_median', 'rf_mean', 'rf_median'],
)
wide.to_csv(wide_path)

print('rows', len(df))
print('repeat counts')
print(df.groupby(['missing_rate', 'algorithm']).size().to_string())
print('\nsummary')
print(summary.round(4).to_string(index=False))
print(f'wrote {summary_path}')
print(f'wrote {wide_path}')
