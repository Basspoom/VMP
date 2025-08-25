import pandas as pd
import joblib 
import sys

def predict_with_model(input_file, model_file, output_file):
    model = joblib.load(model_file)
    
    data_to_predict = pd.read_csv(input_file, sep='\t')
    X_to_predict = data_to_predict.iloc[:, 1:] 
    predicted_scores = model.predict_proba(X_to_predict)[:, 1]
    
    predicted_labels = ['Virus' if score >= 0.5 else 'Not Virus' for score in predicted_scores]
    
    results_df = pd.DataFrame({'ID': data_to_predict['Sequence ID'], 'Score': predicted_scores, 'Prediction': predicted_labels})
    
    results_df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python predict.py input_file model_file output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    model_file = sys.argv[2]
    output_file = sys.argv[3]
    
    predict_with_model(input_file, model_file, output_file)