import pandas as pd
import joblib
import importlib.util

# Load model
model = joblib.load("model.pkl")
features_list = pd.read_csv("features_list.csv")["Feature"].tolist()

# Load feature extractor
spec = importlib.util.spec_from_file_location("extractor", "extract_strand_features_fixed.py")
extractor = importlib.util.module_from_spec(spec)
spec.loader.exec_module(extractor)
print("\n")
print("🧬 MangoneLab - Caligula miRNA Strand Prediction - Terminal Version - v1.0 \n")

# Input from user
# Hardcoded test inputs (let-7a)
seq = "seq: ugggaUGAGGUAGUAGGUUGUAUAGUUuuuCUAUACAAUCUACUGUCUUUCcua"
struct = "struct: (((((((.((((((((((((((((.....))))))))))))))).)))))))"
print("Example - Human let-7a")
print(seq)
print(struct)
#seq = input("Enter hairpin sequence (5p and 3p regions in uppercase):\n\n> ").strip()
#struct = input("\nEnter dot-bracket structure:\n> ").strip()

# Validate
if not seq or not struct:
    print("\n❌ Error: Both sequence and structure are required.")
    exit(1)

# Run prediction
try:
    df = pd.DataFrame([{"hairpin_seq": seq, "hairpin_structure": struct}])
    feats = extractor.extract_features_batch_robust(df)
    for f in features_list:
        if f not in feats.columns:
            feats[f] = 0.0
    feats = feats[features_list]

    probs = model.predict_proba(feats)[0]
    pred = model.predict(feats)[0]

    print("\n✅ Prediction Complete:")
    print(f"Predicted Strand: {'3p' if pred == 1 else '5p'}")
    print(f"Probability 5p: {probs[0]:.2f}")
    print(f"Probability 3p: {probs[1]:.2f}")
except Exception as e:
    print(f"\n❌ Prediction failed: {e}")
