import os
import pandas as pd
import argparse


def calculate_dvf(score, p_value, dvf_score1, dvf_score2, dvf_pvalue1, dvf_pvalue2):
    total_scores = []
    try:
        score = float(score)
        p_value = float(p_value)
        score1 = 1 if score >= dvf_score1 and p_value <= dvf_pvalue1 else 0
        score2 = 1 if score >= dvf_score2 and p_value <= dvf_pvalue2 else 0
        total_scores = [score1, score2]
    except ValueError:
        total_scores = [0, 0]
    return total_scores

    
    

def process_deepvirfinder_results(input_file, output_file, dvf_score1, dvf_score2, dvf_pvalue1, dvf_pvalue2):
    df_deepvirfinder = pd.read_csv(input_file, sep='\t', names=['Sequence ID', 'Some Other Column1', 'score', 'pvalue'], skiprows=1)
    scores_list = df_deepvirfinder.apply(lambda row: calculate_dvf(row['score'], row['pvalue'], dvf_score1, dvf_score2, dvf_pvalue1, dvf_pvalue2), axis=1)
    df_deepvirfinder['Score1'], df_deepvirfinder['Score2'] = zip(*scores_list)
    df_deepvirfinder.to_csv(output_file, index=False, sep='\t', columns=['Sequence ID', 'Score1', 'Score2'])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process DeepVirFinder results and calculate scores for multiple files. \n Usage: python3 VPAC_part2_dvf.py /path/to/directory --dvf_score1 0.7 --dvf_score2 0.9 --dvf_pvalue1 0.05 --dvf_pvalue2 0.05')
    parser.add_argument('input_dir', type=str, help='Path to the directory containing DeepVirFinder results files')
    parser.add_argument('--dvf_score1', type=float, default=0.7, help='Threshold score 1')
    parser.add_argument('--dvf_score2', type=float, default=0.9, help='Threshold score 2')
    parser.add_argument('--dvf_pvalue1', type=float, default=0.05, help='Threshold p value 1')
    parser.add_argument('--dvf_pvalue2', type=float, default=0.05, help='Threshold p value 2')
    args = parser.parse_args()

    files = [f for f in os.listdir(args.input_dir) if f.endswith('dvfpred.txt')]

    for file in files:
        input_file = os.path.join(args.input_dir, file)
        output_file = os.path.join(args.input_dir, os.path.splitext(file)[0] + '_dvf.tsv')  
        process_deepvirfinder_results(input_file, output_file, args.dvf_score1, args.dvf_score2, args.dvf_pvalue1, args.dvf_pvalue2)


