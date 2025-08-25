import os
import pandas as pd
import argparse

def calculate_genomad_score(row):
    try:
        if float(row['virus_score']) > 0.7:
            score1 = 1
        else:
            score1 = 0

        if float(row['virus_score']) > 0.9:
            score2 = 1
        else:
            score2 = 0

        total_scores = [score1, score2]
    except KeyError:
        total_scores = [0, 0]

    return total_scores
    
    

def process_genomad_results(input_file, output_file):
    df_vibrant = pd.read_csv(input_file, sep='\t', names=['Sequence ID', 'chromosome_score', 'plasmid_score', 'virus_score'], skiprows=1)
    scores_list = df_vibrant.apply(calculate_genomad_score, axis=1)
    df_vibrant['Score1'], df_vibrant['Score2'] = zip(*scores_list)
    df_vibrant.to_csv(output_file, index=False, sep='\t', columns=['Sequence ID', 'Score1', 'Score2'])



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process geNomad results and calculate scores for multiple files. \n Usage: python3 VPAC_part2_genomad.py /path/to/directory')
    parser.add_argument('input_dir', type=str, help='Path to the directory containing geNomad results files')
    args = parser.parse_args()

    files = [f for f in os.listdir(args.input_dir) if f.endswith('contigs_aggregated_classification.tsv')]

    for file in files:
        input_file = os.path.join(args.input_dir, file)
        output_file = os.path.join(args.input_dir, file.replace('.tsv', '_genomad.tsv'))  
        process_genomad_results(input_file, output_file)
