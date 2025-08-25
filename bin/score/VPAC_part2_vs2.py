import os
import pandas as pd
import argparse

def calculate_score(score, hall_mark_gene, vs2_score1, vs2_score2, vs2_score3, hm_gene1, hm_gene2, hm_gene3):
    total_scores = []
    try:
        score = float(score)
        hall_mark_gene = float(hall_mark_gene)
        score1 = 1 if score >= float(vs2_score1) and hall_mark_gene >= float(hm_gene1) else 0
        score2 = 1 if score >= float(vs2_score2) and hall_mark_gene >= float(hm_gene2) else 0
        score3 = 1 if score >= float(vs2_score3) and hall_mark_gene >= float(hm_gene3) else 0
        total_scores = [score1, score2, score3]
    except ValueError:
        total_scores = [0, 0, 0]
    return total_scores

def process_virsorter2_results(input_file, output_file, vs2_score1, vs2_score2, vs2_score3, hm_gene1, hm_gene2, hm_gene3):
    df_virsorter2 = pd.read_csv(input_file, sep='\t', names=['Sequence ID', 'Some Other Column1', 'Some Other Column2', 'Score', 'Some Other Column3', 'Some Other Column4', 'Hall Mark Gene', 'Some Other Column5', 'Some Other Column6'], skiprows=1)
    df_virsorter2['Sequence ID'] = df_virsorter2['Sequence ID'].apply(lambda x: x.split('||')[0])
    scores_list = df_virsorter2.apply(lambda row: calculate_score(row['Score'], row['Hall Mark Gene'], vs2_score1, vs2_score2, vs2_score3, hm_gene1, hm_gene2, hm_gene3), axis=1)
    df_virsorter2['Score1'], df_virsorter2['Score2'], df_virsorter2['Score3'] = zip(*scores_list)
    df_virsorter2.to_csv(output_file, index=False, sep='\t', columns=['Sequence ID', 'Score1', 'Score2', 'Score3'])



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process VirSorter2 results and calculate scores for multiple files. \n Usage: python3 VPAC_part3_vs2.py /path/to/directory --vs2_score1 0.5 --vs2_score2 0.7 --vs2_score3 0.95 --hm_gene1 1 --hm_gene2 1 --hm_gene3 0')
    parser.add_argument('input_dir', type=str, help='Path to the directory containing VirSorter2 results files')
    parser.add_argument('--vs2_score1', type=float, default=0.5, help='Threshold score 1')
    parser.add_argument('--vs2_score2', type=float, default=0.7, help='Threshold score 2')
    parser.add_argument('--vs2_score3', type=float, default=0.95, help='Threshold score 3')
    parser.add_argument('--hm_gene1', type=int, default=1, help='Threshold hall mark gene 1')
    parser.add_argument('--hm_gene2', type=int, default=1, help='Threshold hall mark gene 2')
    parser.add_argument('--hm_gene3', type=int, default=0, help='Threshold hall mark gene 3')
    args = parser.parse_args()

    files = [f for f in os.listdir(args.input_dir) if f.endswith('final-viral-score.tsv')]

    for file in files:
        input_file = os.path.join(args.input_dir, file)
        output_file = os.path.join(args.input_dir, os.path.splitext(file)[0] + '_vs2.tsv')  
        process_virsorter2_results(input_file, output_file, args.vs2_score1, args.vs2_score2, args.vs2_score3, args.hm_gene1, args.hm_gene2, args.hm_gene3)