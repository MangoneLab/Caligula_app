import pandas as pd
import joblib
import importlib.util
import sys

# === Load model and features ===
model = joblib.load("model.pkl")
features_list = pd.read_csv("features_list.csv")["Feature"].tolist()

# === Load feature extractor ===
spec = importlib.util.spec_from_file_location("extractor", "extract_strand_features_fixed.py")
extractor = importlib.util.module_from_spec(spec)
spec.loader.exec_module(extractor)
print("\n")
print("🧬 MangoneLab - Caligula miRNA Strand Prediction - Batch Version - v1.0 \n")

# === Input file ===
if len(sys.argv) < 2:
    print("Usage: python caligula_batch.py example.csv\n")
    sys.exit(1)

input_file = sys.argv[1]
output_file = input_file.replace(".csv", "_predicted.csv")

# === Load input ===
df = pd.read_csv(input_file)
if 'hairpin_seq' not in df.columns or 'hairpin_structure' not in df.columns:
    print("❌ Input CSV must have columns: hairpin_seq,hairpin_structure")
    sys.exit(1)

# === Extract features and predict ===
try:
    feats = extractor.extract_features_batch_robust(df)
    for f in features_list:
        if f not in feats.columns:
            feats[f] = 0.0
    feats = feats[features_list].fillna(0.0)

    probs = model.predict_proba(feats)
    preds = model.predict(feats)

    df["predicted_strand"] = ["3p" if p == 1 else "5p" for p in preds]
    df["prob_5p"] = [round(p[0], 2) for p in probs]
    df["prob_3p"] = [round(p[1], 2) for p in probs]

    # Reorder columns if 'name' exists
    cols = ['name', 'hairpin_seq', 'hairpin_structure', 'predicted_strand', 'prob_5p', 'prob_3p'] \
        if 'name' in df.columns else \
        ['hairpin_seq', 'hairpin_structure', 'predicted_strand', 'prob_5p', 'prob_3p']

    df[cols].to_csv(output_file, index=False)
    print(f"✅ Predictions saved to: {output_file}")
except Exception as e:
    print(f"❌ Prediction failed: {e}")
