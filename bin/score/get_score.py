import os
import re
import pandas as pd

folders_to_create = [
    "score",
    "summary",
    "./score/part1",
    "./score/part2",
    "./score/part3"
]

for folder_path in folders_to_create:
    os.makedirs(folder_path, exist_ok=True)

os.system("cp ./results/dvf/*.contigs.fasta_gt1bp_dvfpred_dvf.tsv ./score/part1")
os.system("cp ./results/vs2/*-final-viral-score_vs2.tsv ./score/part1")
os.system("cp ./results/vibr/VIBRANT_genome_quality_*.contigs_vibrant.tsv ./score/part1")

os.system("cp ./results/kaiju/*.kaiju.names_kaiju_scores.tsv ./score/part2")
os.system("cp ./results/checkV/*.quality_summary_checkV_add.tsv ./score/part2")
os.system("cp ./results/vs2/*-final-viral-score_vs2_add.tsv ./score/part2")

os.system("cp ./results/checkV/*.quality_summary_checkV_rmv.tsv ./score/part3")
os.system("cp ./results/vs2/*-final-viral-score_vs2_rmv.tsv ./score/part3")




path=os.getcwd()
path=path+'/raw_contigs/'
fasta_list=[]
for i in os.listdir(path):
	if os.path.isdir(i):
		continue
	else:
		if 'contigs.fa' in i:
			fasta_list.append(i[0:i.index('.')])
		else:
			continue
contigs_folder = "./filtered_contigs/"
score_folders = ["./score/part1", "./score/part2", "./score/part3"]

score_files_info = {
    "part1": ["*.contigs.fasta_gt1bp_dvfpred_dvf.tsv", "*-final-viral-score_vs2.tsv", "VIBRANT_genome_quality_*.contigs_vibrant.tsv"],
    "part2": ["*.kaiju.names_kaiju_scores.tsv", "*.quality_summary_checkV_add.tsv", "*-final-viral-score_vs2_add.tsv"],
    "part3": ["*.quality_summary_checkV_rmv.tsv", "*-final-viral-score_vs2_rmv.tsv"]
}

def custom_sort(seq_id):
    match = re.match(r"k\d+_(\d+)", seq_id)
    if match:
        return int(match.group(1))
    return 0

for sample_name in fasta_list:
    contigs_file = os.path.join(contigs_folder, f"{sample_name}.contigs.fasta")
    sequence_ids = set()
    with open(contigs_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                sequence_id = line.strip().split()[0][1:]  
                sequence_ids.add(sequence_id)
        
    sequence_ids = sorted(sequence_ids, key=custom_sort)
    summary_df = pd.DataFrame(columns=["Sequence ID"])
    summary_df["Sequence ID"] = sequence_ids


    for step, score_folder in enumerate(score_folders, start=1):
        score_files_pattern = score_files_info.get(f"part{step}")
        if score_files_pattern is None:
            print(f"未找到与步骤 part{step} 相关的得分文件信息")
            continue
        
        for pattern in score_files_pattern:
            pattern_with_sample = pattern.replace("*", sample_name)
            score_files = [file_name for file_name in os.listdir(score_folder) if file_name.endswith(pattern_with_sample)]
            
            for score_file in score_files:
                score_file_path = os.path.join(score_folder, score_file)
                
                with open(score_file_path, 'r') as f:
                    header = f.readline().strip().split('\t')  
                    score_columns = header[1:]
                    for line in f:
                        parts = line.strip().split('\t')
                        sequence_id = parts[0].split()[0] 
                        scores = list(map(float, parts[1:]))  #
                        if sequence_id in sequence_ids:
                            for col_name, score in zip(score_columns, scores):
                                col_name = f"{step}_{score_file.split('.')[0]}_score{col_name}"  
                                if col_name not in summary_df.columns:
                                    summary_df[col_name] = 0
                                summary_df.loc[summary_df["Sequence ID"] == sequence_id, col_name] = score
                    
    output_file = f"./summary/{sample_name}_summary.tsv"
    summary_df.to_csv(output_file, index=False, sep='\t')
