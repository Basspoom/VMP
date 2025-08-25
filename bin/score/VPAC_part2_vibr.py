import os
import pandas as pd
import argparse

def calculate_vibrant_score(row):
    try:
        if row['Quality'] == 'complete circular':
            score1 = 1
        else:
            score1 = 0

        if row['Quality'] == 'high quality draft':
            score2 = 1
        else:
            score2 = 0

        if row['Quality'] == 'medium quality draft':
            score3 = 1
        else:
            score3 = 0

        if row['Quality'] == 'low quality draft' and row['type'] == 'lytic':
            score4 = 1
        else:
            score4 = 0

        total_scores = [score1, score2, score3, score4]
    except KeyError:
        total_scores = [0, 0, 0, 0]

    return total_scores
    
    

def process_vibrant_results(input_file, output_file):
    df_vibrant = pd.read_csv(input_file, sep='\t', names=['Sequence ID', 'type', 'Quality'], skiprows=1)
    scores_list = df_vibrant.apply(calculate_vibrant_score, axis=1)
    df_vibrant['Score1'], df_vibrant['Score2'], df_vibrant['Score3'], df_vibrant['Score4'] = zip(*scores_list)
    df_vibrant.to_csv(output_file, index=False, sep='\t', columns=['Sequence ID', 'Score1', 'Score2', 'Score3', 'Score4'])



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process VIBRANT results and calculate scores for multiple files. \n Usage: python3 VPAC_part4_vibrant.py /path/to/directory')
    parser.add_argument('input_dir', type=str, help='Path to the directory containing VIBRANT results files')
    args = parser.parse_args()

    files = [f for f in os.listdir(args.input_dir) if f.endswith('.contigs.tsv')]

    for file in files:
        input_file = os.path.join(args.input_dir, file)
        output_file = os.path.join(args.input_dir, file.replace('.tsv', '_vibrant.tsv'))
        process_vibrant_results(input_file, output_file)
