
import streamlit as st
import joblib
import numpy as np
from extract_strand_features_fixed import extract_strand_features
from sklearn.exceptions import NotFittedError

# Load model
@st.cache_resource
def load_model():
    return joblib.load("model.pkl")

model = load_model()

st.title("Caligula: miRNA Strand Selection Predictor")

st.markdown("""
This tool predicts whether the 5p or 3p strand of a miRNA hairpin is more likely to be selected based on structural and sequence-derived features.  
**Input your hairpin sequence and dot-bracket structure below.**  
Make sure the 5p and 3p regions are in uppercase.
""")

# User input
sequence = st.text_area("Hairpin sequence (5p and 3p strands must be in uppercase)", height=150)
structure = st.text_area("Hairpin dot-bracket structure (same length as sequence)", height=100)

# Run prediction
if st.button("Predict Strand"):
    if not sequence or not structure:
        st.warning("Please enter both sequence and structure.")
    elif len(sequence.strip()) != len(structure.strip()):
        st.error("Sequence and structure must be the same length.")
    else:
        try:
            features = extract_strand_features(sequence.strip(), structure.strip())
            if features is None:
                st.error("Feature extraction failed. Make sure there are exactly two uppercase regions.")
            else:
                X = np.array([features])
                prediction = model.predict(X)[0]
                prob = model.predict_proba(X)[0]
                conf_score = max(prob)

                st.success(f"**Predicted Strand: {prediction}**")
                st.info(f"Confidence: {conf_score:.3f}")
        except NotFittedError:
            st.error("Model not loaded correctly.")
        except Exception as e:
            st.error(f"An error occurred: {e}")
