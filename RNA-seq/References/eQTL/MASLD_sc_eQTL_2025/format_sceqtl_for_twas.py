import pandas as pd
import gzip
import argparse
import os

def format_sceqtl(input_file, output_file):
    print(f"Reading {input_file}...")

    # Check if empty or missing
    if not os.path.exists(input_file):
        print(f"File {input_file} not found.")
        return

    # Using chunksize for memory efficiency since the file has ~4.3 million rows
    chunksize = 10 ** 5
    first_chunk = True

    # Columns needed: require at least all_cells beta/se/pval to be non-NA
    required_cols = ['snp_chr_pos_ref_alt', 'Beta_all_cells', 'pval_all_cells']

    with pd.read_csv(input_file, sep='\t', compression='gzip', chunksize=chunksize) as reader:
        for chunk in reader:
            # Drop NaN rows in crucial columns (all_cells is primary)
            chunk = chunk.dropna(subset=required_cols)

            # Parse snp_chr_pos_ref_alt (e.g., 10_100000943_G_A_hg38)
            split_snp = chunk['snp_chr_pos_ref_alt'].str.split('_', expand=True)

            if split_snp.shape[1] < 4:
                print("Warning: Malformed snp string encountered, skipping chunk rows.")
                continue

            chunk['chromosome'] = 'chr' + split_snp[0]
            chunk['position'] = split_snp[1]
            chunk['ref'] = split_snp[2]
            chunk['alt'] = split_snp[3]

            chunk['variant_id'] = chunk['chromosome'] + '_' + chunk['position'] + '_' + chunk['ref'] + '_' + chunk['alt'] + '_b38'

            # Output all three conditions: all_cells (primary), MASLD_cells, control_cells
            out_df = pd.DataFrame({
                'gene_id': chunk['gene'],
                'variant_id': chunk['variant_id'],
                'chromosome': chunk['chromosome'],
                'position': chunk['position'],
                'ref': chunk['ref'],
                'alt': chunk['alt'],
                'beta_all': chunk['Beta_all_cells'],
                'se_all': chunk['Std_Error_all_cells'],
                'pval_all': chunk['pval_all_cells'],
                'beta_masld': chunk['Beta_MASLD_cells'],
                'se_masld': chunk['Std_Error_MASLD_cells'],
                'pval_masld': chunk['pval_MASLD_cells'],
                'beta_ctrl': chunk['Beta_control_cells'],
                'se_ctrl': chunk['Std_Error_control_cells'],
                'pval_ctrl': chunk['pval_control_cells'],
            })

            # Write to output file
            if first_chunk:
                out_df.to_csv(output_file, sep='\t', index=False, mode='w', compression='gzip')
                first_chunk = False
            else:
                out_df.to_csv(output_file, sep='\t', index=False, mode='a', header=False, compression='gzip')

    print(f"Successfully formatted {input_file} -> {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Format MASLD sc-eQTL summary stats for TWAS/MR.")
    parser.add_argument("-i", "--input", required=True, help="Input gzipped summary stats file")
    parser.add_argument("-o", "--output", required=True, help="Output formatted gzipped file")
    args = parser.parse_args()

    format_sceqtl(args.input, args.output)
